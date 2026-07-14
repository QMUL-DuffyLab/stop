module mc
  use mpi_f08
  use io
  use lattice
  implicit none
  private
  logical(kind=CB) :: possible(6) ! 6 types of moves considered
  logical(kind=CB) :: bad_proposal
  integer(kind=CI) :: n_accepted(6)
  integer(kind=CI), allocatable :: gens(:), anns(:)
  real(kind=CF), allocatable :: pulse(:)
  real(kind=CF) :: mu, pulse_tmax
  real(kind=CF), allocatable :: gen_rates(:, :),&
    hop_rates(:, :, :), intra_rates(:, :, :), ann_rates(:, :, :)
  integer(kind=CI) :: rep
  type :: move_type
    integer(kind=CI) :: mt  = 0_CI
    integer(kind=CI) :: isi = 0_CI
    integer(kind=CI) :: ist = 0_CI
    integer(kind=CI) :: fsi = 0_CI
    integer(kind=CI) :: fst = 0_CI
    integer(kind=CI) :: loss_index = 0_CI
    logical(kind=CB) :: emissive = .false._CB
    real(kind=CF) :: rate = 0.0_CF
  end type
  type(move_type) :: cm ! current move
  public :: construct_pulse, do_run, rep, n_accepted, mc_deallocations
  
  contains

    subroutine init_random(salt)
      ! function to seed the actual RNG : cf.
      ! https://gcc.gnu.org/onlinedocs/gcc-4.9.1/gfortran/RANDOM_005fSEED.html
      use iso_fortran_env, only: int64
      integer(kind=CI), allocatable :: seed(:)
      integer(kind=CI) :: i, n, un, istat, dt(8), salt
      integer(int64) :: t
          
      call random_seed(size = n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream", &
           form="unformatted", action="read", status="old", &
           iostat=istat)
      if (istat == 0) then
         read(un) seed
         close(un)
      else
         ! Fallback to XOR:ing the current time and a salt.
         ! The salt is taken from the process rank and run
         ! number because this is being done in parallel.
         ! the gnu.org page uses getpid() for this but that's a
         ! GNU extension and i don't want to rely on it
         call system_clock(t)
         if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                 + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                 + dt(3) * 24_int64 * 60 * 60 * 1000 &
                 + dt(5) * 60 * 60 * 1000 &
                 + dt(6) * 60 * 1000 + dt(7) * 1000 &
                 + dt(8)
         end if
         t = ieor(t, int(salt, kind(t)))
         do i = 1, n
            seed(i) = lcg(t)
         end do
      end if
      call random_seed(put=seed)
      deallocate(seed)
      contains
        function lcg(s)
          ! shitty RNG to seed good RNG
          integer :: lcg
          integer(int64) :: s
          if (s == 0) then
            s = 104729
          else
            s = mod(s, 4294967296_int64)
          end if
          s = mod(s * 279470273_int64, 4294967291_int64)
          lcg = int(mod(s, int(huge(0), int64)), kind(0))
        end function lcg
    end subroutine init_random

    function rand_nonzero_int(arr) result(i)
      ! same as rand_nonzero_real but for an integer-valued 1d array
      integer(kind=CI), intent(in) :: arr(:)
      real(kind=CF) :: r
      integer(kind=CI) :: i
      integer(kind=CI), allocatable :: ind(:)
      ind = pack([(i, i = 1_CI, size(arr))], arr.gt.0_CI)
      call random_number(r)
      i = ind(ceiling(r * size(ind)))
      if ((i.lt.1_CI).or.(i.gt.size(arr))) then
        write(*, *) "invalid i returned by rand_nonzero_int"
        write(*, *) "input array: ", arr
        write(*, *) "index array: ", ind
        write(*, *) "rand(): ", r
        stop
      end if
    end function rand_nonzero_int

    function rand_nonzero_real(arr) result(i)
      ! return the index of a random nonzero value in the
      ! real-valued 1d array `arr`
      real(kind=CF), intent(in) :: arr(:)
      real(kind=CF) :: r
      integer(kind=CI) :: i
      integer(kind=CI), allocatable :: ind(:)
      ind = pack([(i, i = 1_CI, size(arr))], arr.gt.0.0_CF)
      call random_number(r)
      i = ind(ceiling(r * size(ind)))
      if ((i.lt.1_CI).or.(i.gt.size(arr))) then
        write(*, *) "invalid i returned by rand_nonzero_real"
        write(*, *) "input array: ", arr
        write(*, *) "index array: ", ind
        write(*, *) "rand(): ", r
        stop
      end if
    end function rand_nonzero_real

    function rand_true(arr) result(i)
      ! same again for a boolean 1d array
      logical(kind=CB), intent(in) :: arr(:)
      real(kind=CF) :: r
      integer(kind=CI) :: i
      integer(kind=CI), allocatable :: ind(:)
      ind = pack([(i, i = 1_CI, size(arr))], arr)
      if (size(ind).eq.0) then
        write(*, *) "size(ind) = 0. arr = ", arr
        stop
      end if
      call random_number(r)
      i = ind(ceiling(r * size(ind)))
      if ((i.lt.1_CI).or.(i.gt.size(arr))) then
        write(*, *) "invalid i returned by rand_true"
        write(*, *) "input array: ", arr
        write(*, *) "index array: ", ind
        write(*, *) "rand(): ", r
        stop
      end if
    end function rand_true

    function rand_int(imax) result(i)
      ! pick an integer from 1 to imax inclusive
      integer(kind=CI), intent(in) :: imax
      real(kind=CF) :: r
      integer(kind=CI) :: i
      call random_number(r)
      i = ceiling(r * imax)
    end function rand_int

    subroutine abundances()
      ! at the start of each run, randomise which sites have each state
      ! based on the abundance fractions
      integer(kind=CI) :: site, state, n_present, rand_site, nn, site2, i
      logical(kind=CB) :: current(n_sites)
      is_present = .true.
      ! turn all the rates on, then zero them as needed
      write(*, *) shape(gen_rates)
      do site = 1, n_sites
        gen_rates(site, :)      = xsec
        ann_rates(site, :, :)   = ann
        intra_rates(site, :, :) = intra
        do site2 = 1, coord
          hop_rates(site, :, site2) = hop
        end do
      end do

      do state = 1, n_s
        ! use the rand_true function above to keep track of which
        ! sites we've picked already so we don't double count
        current = .true.
        ! for each state, check its abundance
        ! if it's 1, we don't need to bother with all this
        if (abundance(state).lt.1.0_CF) then
          ! we've set everything to true and now we're turning
          ! sites off as necessary, so we need 1 - abundance
          n_present = nint((1.0_CF - abundance(state)) * n_sites, kind=CI)
          do site = 1, n_present
            rand_site = rand_true(current)
            is_present(rand_site, state) = .false.
            ! setting rand_site to false excludes it from the next
            ! call to rand_true so we always set the right number of falses
            current(rand_site) = .false.
            ! zero out the relevant rates for this state and site
            gen_rates(rand_site, state) = 0.0_CF
            intra_rates(rand_site, state, :) = 0.0_CF
            intra_rates(rand_site, :, state) = 0.0_CF
            ann_rates(rand_site, state, :) = 0.0_CF
            ann_rates(rand_site, :, state) = 0.0_CF

            ! to zero out the hop rates correctly, first zero out
            ! all the rates for hopping from this state on this site
            do site2 = 1, coord
              nn = neighbours(rand_site, site2)
              hop_rates(rand_site, state, nn) = 0.0_CF
            end do
            ! then loop over all sites and zero out all the hopping
            ! rates *to* this state from the site's neighbours
            do site2 = 1, n_sites
              do i = 1, coord
                if (neighbours(site2, i).eq.rand_site) then
                  hop_rates(site2, state, i) = 0.0_CF
                end if
              end do
            end do

          end do
        end if
      end do

    end subroutine abundances

    subroutine zero_move(m)
      type(move_type) :: m
      m%mt = 0_CI
      m%isi = 0_CI
      m%ist = 0_CI
      m%fsi = 0_CI
      m%fst = 0_CI
      m%loss_index = 0_CI
      m%emissive = .false._CB
      m%rate = 0.0_CF
    end subroutine zero_move

    subroutine construct_pulse(fwhm, dt, fluence)
      real(kind=CF) :: fwhm, sigma, dt, fluence
      integer :: i, pulse_len
      mu = 2.0 * fwhm
      pulse_tmax = 2.0 * mu
      pulse_len = int(pulse_tmax / dt)
      sigma = (fwhm) / (2.0 * (sqrt(2.0 * log(2.0))))
      allocate(pulse(pulse_len), source=0.0_CF)
      do i = 1, pulse_len
        pulse(i) = (fluence) / (sigma * sqrt(2.0 * pi)) * &
          exp(-1.0 * ((i * dt) - (mu))**2 / (sqrt(2.0) * sigma)**2)
      end do
    end subroutine construct_pulse

    subroutine possible_moves(site, ft)
      ! site - integer index into n_sites
      real(kind=CF), intent(in) :: ft 
      integer, intent(in) :: site
      integer :: state
      possible = .false._CB
      ! check that no process has reduced population where it shouldn't
      ! - don't think this could ever happen but checking it
      if (any(n_i(site, :).lt.0)) then
        write(*, *) "negative population on site", site, n_i(site, :)
        write(*, *) "current move", cm
        stop
      end if
      
      do state = 1, n_s
        if ((.not.is_present(site, state)).and.(n_i(site, state).gt.0_CI)) then
          write(*, *) "population on forbidden state!"
          write(*, '(a, I0, a, I0)') "site = ", site, " state = ", state
          write(*, *) is_present(site, :)
          stop
        end if
      end do
      ! if the pulse is on (ft) and there's a positive cross section
      ! then excitations can be generated. in general the second part
      ! should always be true, i guess, otherwise what are we doing here
      if ((ft.gt.0.0_CF).and.(any(gen_rates(site, :).gt.0.0_CF))) then
        ! there is a state on this site that can absorb
        possible(1) = .true.
        if (dot_product(n_i(site, :), gen_rates(site, :)).gt.0.0_CF) then
          ! for SE to be possible there must be population on a
          ! site that can also generate excitations; the dot product
          ! of n_i and gen_rates can only be > 0 if this is true
          possible(2) = .true.
        end if
      end if

      ! only bother checking for hops if there's a nonzero rate
      if (any(hop.gt.0.0_CF)) then
        hoploop: do state = 1, n_s
          if (any(n_i(site, state) * hop_rates(site, state, :).gt.0.0_CF)) then
            possible(3) = .true._CB
            ! we've found one, so it's possible; stop checking
            exit hoploop
          end if
        end do hoploop
      end if

      decayloop: do state = 1, n_s
        if (n_i(site, state)*intra_rates(site, state, state).gt.0.0_CF) then
          possible(5) = .true._CB
          exit decayloop
        end if
      end do decayloop

      if ((n_s.gt.1).and.(any(n_i(site, :).gt.0))) then
        ! for transfer between states to be possible, there
        ! must be more than one state on the protein
        possible(4) = .true._CB
      end if

      if ((sum(n_i(site, :)).gt.1_CI).and.&
        (any(ann_rates(site, :, :).gt.0.0_CF))) then
        ! this is not strictly true in all cases but much quicker
        possible(6) = .true._CB
      end if

    end subroutine possible_moves

    subroutine propose(site, ft)
      integer, intent(in) :: site
      real(kind=CF), intent(in) :: ft 
      integer :: s, s2, nn, n_eff, p
      ! make sure there's no leftover stuff in current_move
      ! that might mess with this
      call zero_move(cm)
      call possible_moves(site, ft)
      ! indices where possible is true
      if (.not.any(possible)) then
        ! something has gone badly wrong somewhere
        bad_proposal = .true.
        return
      end if

      cm%mt = rand_true(possible)
      cm%isi = site
      cm%fsi = site

      select case (cm%mt)
      case (1)
        ! generation
        ! need to pick a non-zero cross section and make note of
        ! which index that cross section is at
        s = rand_nonzero_real(gen_rates(site, :)) ! which state is it
        p = which_p(s) ! which pigment is this state on
        cm%rate = ft * gen_rates(site, s) * (n_tot(p) - n_i(site, s)) / n_tot(p)
        cm%ist = s
        cm%fst = s
        cm%loss_index = s
      case (2)
        ! stimulated emission
        s = rand_nonzero_real(gen_rates(site, :))
        p = which_p(s)
        cm%rate = ft * gen_rates(site, s) * n_i(site, s) / n_tot(p)
        cm%ist = s
        cm%fst = s
        cm%loss_index = n_s + s
      case (3)
        ! hop
        s = rand_nonzero_int(n_i(site, :))
        s2 = rand_nonzero_int(neighbours(site, :))
        nn = neighbours(site, s2)
        p = which_p(s)
        ! this will be zero if the state s does not exist on nn
        cm%rate = n_i(site, s) * hop_rates(site, s, nn)
        if (n_i(site, s).lt.n_i(nn, s)) then
          ! entropic penalty for hopping to a higher-occupied site
          cm%rate = cm%rate * (n_i(site, s) * (n_thermal(p) - n_i(nn, s))) /&
          ((n_i(nn, s) + 1) * (n_thermal(p) - n_i(site, s) + 1))
        end if
        cm%ist = s
        cm%fsi = nn
        cm%fst = s
        cm%loss_index = (2 * n_s) + s
      case (4)
        ! transfer
        ! pick a state with a nonzero rate
        s = rand_nonzero_int(n_i(site, :))
        s2 = rand_int(n_s)
        do while (s2.eq.s)
          ! intra has decays on the diagonal and transfer rates off-diagonal
          ! so ensure we've picked an off-diagonal element
          s2 = rand_int(n_s)
        end do
        cm%rate = n_i(site, s) * intra_rates(site, s, s2)
        cm%ist = s
        cm%fst = s2
        cm%loss_index = (4 + (s - 1)) * n_s + s2
      case (5)
        ! decay
        s = rand_nonzero_int(n_i(site, :))
        cm%rate = n_i(site, s) * intra_rates(site, s, s)
        cm%ist = s
        cm%fst = s
        cm%emissive = emissive(s)
        cm%loss_index = (3 * n_s) + s
      case (6)
        ! annihilation
        s = rand_nonzero_int(n_i(site, :))
        s2 = rand_nonzero_int(n_i(site, :))
        ! sort so that s < s2. if we generate a pair of columns in the
        ! histogram for each annihilation process, this ensures that
        ! all annihilation events go in one of them to make counting
        ! easier. if we're only generating one column per annihilation
        ! pair then the loss index calculation will have to be modified
        ! here - come back to this if so
        if (s.gt.s2) then
          n_eff = s
          s = s2
          s2 = n_eff
        end if
        if (s.eq.s2) then
          cm%rate = n_i(site, s) * ((n_i(site, s) - 1.0_CF) / 2.0) *&
            ann_rates(site, s, s)
        else
          cm%rate = n_i(site, s) * n_i(site, s2) * ann_rates(site, s, s2)
        end if
        cm%loss_index = (4 + n_s + (s - 1)) * n_s + s2
        cm%ist = s
        cm%fst = s2
      case default
        write(*, *) "shouldn't be here! cm type wrong in propose()"
        stop
      end select
    end subroutine propose
    
    subroutine do_move()
      integer :: a
      select case(cm%mt)
      case (1) ! gen
        n_i(cm%isi, cm%ist) = n_i(cm%isi, cm%ist) + 1
      case (2) ! se
        n_i(cm%isi, cm%ist) = n_i(cm%isi, cm%ist) - 1
      case (3) ! hop
        n_i(cm%isi, cm%ist) = n_i(cm%isi, cm%ist) - 1
        n_i(cm%fsi, cm%fst) = n_i(cm%fsi, cm%fst) + 1
      case (4) ! transfer
        n_i(cm%isi, cm%ist) = n_i(cm%isi, cm%ist) - 1
        n_i(cm%fsi, cm%fst) = n_i(cm%fsi, cm%fst) + 1
      case (5) ! decay
        n_i(cm%isi, cm%ist) = n_i(cm%isi, cm%ist) - 1
      case (6) ! annihilation
        a = ann_remainder(cm%ist, cm%fst) ! what state remains?
        n_i(cm%isi, cm%ist) = n_i(cm%isi, cm%ist) - 1
        n_i(cm%isi, cm%fst) = n_i(cm%isi, cm%fst) - 1
        n_i(cm%isi, a) = n_i(cm%isi, a) + 1
      case default
        write(*, *) "something wrong in do_move!", cm%mt
        stop
      end select
      if (any(n_i(:, :).lt.0)) then
        write(*, *) "do_move: negative population."
        write(*, *) "current move: ", cm
        write(*, *) "site population: ", n_i(cm%isi, :)
        stop
      end if
    end subroutine do_move

    subroutine mc_step(t, dt, bin)
      real(kind=CF) :: ft, prob, r, t, dt
      logical(kind=CB) :: bin
      integer(kind=CI) :: pind, i, j, s, n_poss, ibin

      n_accepted = 0_CI
      ft = 0.0_CF

      if (t.le.2.0*mu) then
        pind = int(t / dt) + 1
        if (pind.le.size(pulse)) then
          ft = pulse(pind)
        end if
      end if

      do i = 1, n_sites
        bad_proposal = .false.

        s = rand_int(n_sites)

        ! first pick a random site, check if any moves are possible
        call possible_moves(s, ft)
        n_poss = count(possible)
        if (n_poss.eq.0) then
          ! if not just try again
          cycle
        end if

        ! some moves are possible. try each one once *on average*
        moveloop: do j = 1, n_poss

          call propose(s, ft)

          if (bad_proposal) then
            write(*, *) "bad proposal. why?"
            write(*, *) "site = ", s
            write(*, *) "ni = ", n_i(i, :)
            write(*, *) "time = ", t
            write(*, *) "dt = ", dt
            write(*, *) "possible = ", possible
            write(*, *) "any(possible)", any(possible)
            write(*, *) "not(any(possible))", (.not.any(possible))
            write(*, *) "n_poss", n_poss
            stop
          end if

          prob = cm%rate * dt * exp(-1.0_CF * cm%rate * dt)
          call random_number(r)

          if (r.lt.prob) then
            n_accepted(cm%mt) = n_accepted(cm%mt) + 1
            call do_move()

            if (bin.and.cm%loss_index.gt.0) then
              ibin = ceiling(t / binwidth)
              if (ibin.eq.0) then
                ibin = 1 ! if t = 0 ceiling returns 0; put in first bin
              end if
              counts(ibin, cm%loss_index) = counts(ibin, cm%loss_index) + 1
              site_moves(cm%isi, cm%loss_index) = &
                site_moves(cm%isi, cm%loss_index) + 1

              if (cm%mt.eq.1) then
                gens(cm%isi) = gens(cm%isi) + 1
              end if
              if (cm%mt.eq.6) then
                anns(cm%isi) = anns(cm%isi) + 1
              end if
            end if
          end if

          ! check if any more moves are possible now
          call possible_moves(s, ft)

          if (.not.any(possible)) then
            exit moveloop
          end if

        end do moveloop

      end do

    end subroutine mc_step

    subroutine do_run(salt, max_counts, out_file_path)
      integer(kind=CI), intent(in) :: salt, max_counts
      real(kind=CF) :: t, interval
      integer(kind=CI) :: tot_accepted(6)
      integer(kind=CI) :: i, j, curr_maxcount, nunit
      integer(kind=CI), allocatable :: ec(:)
      character(200) :: out_file_path, outfile, pop_file
      logical(kind=CB) :: skip, bin_pulse

      call init_random(salt)

      allocate(gens(n_sites), source=0_CI)
      allocate(anns(n_sites), source=0_CI)
      allocate(gen_rates(n_sites, n_s), source=0.0_CF)
      allocate(hop_rates(n_sites, n_s, coord), source=0.0_CF)
      allocate(intra_rates(n_sites, n_s, n_s), source=0.0_CF)
      allocate(ann_rates(n_sites, n_s, n_s), source=0.0_CF)

      ! abundance call should go here, at the start of each run
      ! and after random init but before we do anything
      call abundances()

      counts = 0_CI
      tot_accepted = 0_CI

      ! first column in emissive columns isn't a real one -
      ! it's for the output file. could maybe replace this
      ! with just emissive??
      ec = pack([(i, i = 1_CI, size(emissive_columns) - 1)],&
        emissive_columns(2:))

      rep = 1_CI
      curr_maxcount = 0_CI
      interval = 1.0 / rep_rate

      write(pop_file, '(a, a, I0, a, I0, a)') trim(adjustl(out_file_path)),&
        "_salt_", salt, "_population.csv"
      open(newunit=nunit, file=pop_file)

      reploop: do while (curr_maxcount.lt.max_counts)

        gens = 0_CI
        anns = 0_CI

        bin_pulse = .true.
        if (rep.lt.burn_reps) then
          bin_pulse = .false.
        end if

        t = 0.0_CF
        skip = .false.

        pulseloop: do while (t.lt.tmax)
          call mc_step(t, dt1, bin_pulse)
          tot_accepted = tot_accepted + n_accepted
          t = t + dt1
          if ((t.gt.(pulse_tmax)).and.(sum(n_i).eq.0_CI)) then
            skip = .true.
            exit pulseloop ! start next rep now, nothing to simulate
          end if
        end do pulseloop

        if (.not.skip) then
          darkloop: do while (t.lt.interval)
            call mc_step(t, dt2, .false._CB)
            t = t + dt2
            if (sum(n_i).eq.0_CI) then
              exit darkloop ! start next rep now, nothing to simulate
            end if
          end do darkloop
        end if

        if (mod(rep, 10).eq.0) then
          do i = 1, n_bins
            do j = 1, size(ec)
              if (counts(i, ec(j)).gt.curr_maxcount) then
                curr_maxcount = counts(i, ec(j))
              end if
            end do
          end do
          write (*, '(a, I0, a, I0, a, I0, a, a)') "salt = ", salt,&
            " rep = ", rep, " current max count = ", curr_maxcount,&
            " type = ", labels(ec + 1) ! labels 1 is "Time (s)"
        end if

        if (debug) then
          do i = 1, n_sites
            j = gens(i)
            if (j.eq.0) then
              j = 1
            else if ((j.gt.0).and.(j.lt.hist_max)) then
              j = j + 1
            else if (j.ge.hist_max) then
              j = hist_max
            end if
            site_gen_hist(i, j) = site_gen_hist(i, j) + 1
            j = anns(i)
            if (j.eq.0) then
              j = 1
            else if ((j.gt.0).and.(j.lt.hist_max)) then
              j = j + 1
            else if (j.ge.hist_max) then
              j = hist_max
            end if
            site_ann_hist(i, j) = site_ann_hist(i, j) + 1
          end do
        end if

        write(nunit, *) rep, sum(n_i, dim=1), t
        rep = rep + 1

      end do reploop

      close(nunit)

      ! max count reached
      write(*, '(a, i0, a, i0, a, i0)') "Process with salt ", salt,&
        " finished at rep ", rep, "with max count ", curr_maxcount
      write(*, '(a, i0, a, G0.6)') "Total generated: ",&
        sum(counts(:, 1:n_s)), ", generated per rep per protein = ",&
        real(sum(counts(:, 1:n_s))) / (rep * n_sites)
      ! NB: the first number in tot_accepted will NOT be the same as
      ! total generated above, because that only counts the binned
      ! generations, i.e. it ignores those in the burn reps
      write(*, '(a, 6(I0, 1X))') "Total accepted: ",&
        tot_accepted
      write(*, *) "Annihilation to decay ratio:",&
        float(tot_accepted(6)) / tot_accepted(5)
      write(outfile, '(a, a, I0, a, I0, a)') trim(adjustl(out_file_path)),&
        "_salt_", salt, "_rep_", rep, "_final.csv"
      call write_histogram(outfile)

      if (debug) then
        write(outfile, '(a, a, I0, a, I0, a)') trim(adjustl(out_file_path)),&
          "_salt_", salt, "_rep_", rep, "_site_moves.csv"
        call write_site_moves(outfile)
        write(outfile, '(a, a, I0, a, I0, a)') trim(adjustl(out_file_path)),&
          "_salt_", salt, "_rep_", rep, "_"
        call write_move_hists(outfile)
      end if

      deallocate(gens)
      deallocate(anns)
      deallocate(gen_rates)
      deallocate(hop_rates)
      deallocate(intra_rates)
      deallocate(ann_rates)

    end subroutine do_run

    subroutine mc_deallocations()
      deallocate(pulse)
    end subroutine mc_deallocations

end module mc
