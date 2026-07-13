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
  real(kind=CF), allocatable :: gen_rates(:), se_rates(:),&
    hop_rates(:, :), intra_rates(:, :), ann_rates(:, :)
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

    function rand_nonzero_real_2d(arr) result(i)
      ! return the index of a random nonzero value in the
      ! real-valued 2d array `arr`. note that pack compresses
      ! the array into 1d, so this will return the 1d index;
      ! you have to figure out the indexing yourself
      real(kind=CF), intent(in) :: arr(:, :)
      real(kind=CF) :: r
      integer(kind=CI) :: i, sa(2)
      integer(kind=CI), allocatable :: ind(:)
      real(kind=CF), allocatable :: arr_1d(:)
      sa = shape(arr)
      arr_1d = reshape(arr, (/product(sa)/))
      ind = pack([(i, i = 1_CI, size(arr_1d))], arr_1d.gt.0.0_CF)
      call random_number(r)
      i = ind(ceiling(r * size(ind)))
      if ((i.lt.1_CI).or.(i.gt.size(arr))) then
        write(*, *) "invalid i returned by rand_nonzero_real"
        write(*, *) "input array: ", arr
        write(*, *) "index array: ", ind
        write(*, *) "rand(): ", r
        stop
      end if
    end function rand_nonzero_real_2d

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
      integer(kind=CI) :: site, state, n_present, rand_site
      logical(kind=CB) :: current(n_sites)
      is_present = .true.
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
          end do
          write(*, '(a, I0, a, F4.2)') "abundance(", state,&
            ") = ", abundance(state)
          write(*, '(a, I0)') "number of trues = ", count(is_present(:, state))
          write(*, '(a, I0)') "number there should be = ",&
            nint(abundance(state) * n_sites)
          write(*, *) is_present(:, state)
          write(*, *)
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

    subroutine possible_moves_and_rates(site, ft)
      real(kind=CF), intent(in) :: ft 
      integer(kind=CI), intent(in) :: site
      integer(kind=CI) :: i, j, s1, s2, nn, p
      integer(kind=CI), allocatable :: present_states(:)
      ! make sure arrays are zeroed
      gen_rates   = 0.0_CF
      se_rates    = 0.0_CF
      hop_rates   = 0.0_CF
      intra_rates = 0.0_CF
      ann_rates   = 0.0_CF
      possible = .false._CB
      ! only insert rates that correspond to states present on this site
      ! (cf. abundances subroutine)
      present_states = pack([(i, i = 1_CI, n_s)], is_present(site, :))

      do i = 1, size(present_states)
        s1 = present_states(i)

        p = which_p(s1)
        gen_rates(s1) = ft * xsec(s1) * (n_tot(p) - n_i(site, s1)) / n_tot(p)
        if (gen_rates(s1).gt.0.0_CF) then
          possible(1) = .true._CB
        end if

        se_rates(s1)  = ft * xsec(s1) * n_i(site, s1) / n_tot(p)
        if (se_rates(s1).gt.0.0_CF) then
          possible(2) = .true._CB
        end if

        if (any(hop.gt.0.0_CF)) then
          ! only bother with this if there is hopping in the first place
          do j = 1, coord
            nn = neighbours(site, j)
            if (nn.gt.0_CI.and.is_present(nn, s1)) then
              possible(3) = .true._CB
              hop_rates(s1, j) = n_i(site, s1) * hop(s1)
              p = which_p(s1)
              if (n_i(site, s1).lt.n_i(nn, s1)) then
                ! entropic penalty for hopping to a higher-occupied site
                hop_rates(s1, j) = hop_rates(s1, j) * &
                  (n_i(site, s1) * (n_thermal(p) - n_i(nn, s1))) /&
                ((n_i(nn, s1) + 1) * (n_thermal(p) - n_i(site, s1) + 1))
              end if
            end if
          end do
        end if

        do j = 1, size(present_states)
          s2 = present_states(j)

          intra_rates(s1, s2) = n_i(site, s1) * intra(s1, s2)
          if (s1 == s2) then
            ann_rates(s1, s2) = n_i(site, s1) * ((n_i(site, s1) - 1) / 2.0) * &
              ann(s1, s2)
            ! diagonal intra rate -> decay
            if (intra_rates(s1, s2).gt.0.0_CF) then
              possible(5) = .true._CB
            end if
          else
            ann_rates(s1, s2) = n_i(site, s1) * n_i(site, s2) * ann(s1, s2)
            ! off-diagonal intra rate -> transfer between states
            if (intra_rates(s1, s2).gt.0.0_CF) then
              possible(4) = .true._CB
            end if
          end if

          if (ann_rates(s1, s2).gt.0.0_CF) then
            possible(6) = .true._CB
          end if

        end do
      end do

      do i = 1, n_s
        if ((.not.is_present(site, i))) then
          if (any(ann_rates(i, :).gt.0.0_CF)&
            .or.any(ann_rates(:, i).gt.0.0_CF)) then
            write(*, *) "is_present(", site, ", ", i, ") = ", is_present(site, i)
            write(*, *) "nonzero ann_rate: "
            do j = 1, n_s
              write(*, *) ann_rates(j, :)
            end do
            stop
          end if
        end if
      end do
    end subroutine possible_moves_and_rates

    subroutine propose(site, ft)
      integer, intent(in) :: site
      real(kind=CF), intent(in) :: ft 
      integer :: s, s2, nn, n_eff
      ! make sure there's no leftover stuff in current_move
      ! that might mess with this
      call zero_move(cm)
      call possible_moves_and_rates(site, ft)
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
        s = rand_nonzero_real(gen_rates) ! which state is it
        cm%rate = gen_rates(s)
        cm%ist = s
        cm%fst = s
        cm%loss_index = s
      case (2)
        ! stimulated emission
        s = rand_nonzero_real(gen_rates) ! which state is it
        cm%rate = se_rates(s)
        cm%ist = s
        cm%fst = s
        ! this assumes there's only one protein in the system
        ! site.e. it'll need changing if there are multiple proteins
        cm%loss_index = n_s + s
      case (3)
        ! hop
        s2 = rand_nonzero_real_2d(hop_rates)
        ! split the returned index (hop is a n_s x coord matrix)
        ! fortran is column-major!!!
        s = mod(s2 - 1, n_s) + 1
        nn = ((s2 - 1) / n_s) + 1
        cm%rate = hop_rates(s, nn)
        cm%ist = s
        cm%fsi = nn
        cm%fst = s
        cm%loss_index = (2 * n_s) + s
      case (4)
        ! transfer
        ! pick a state with a nonzero rate
        nn = rand_nonzero_real_2d(intra_rates)
        ! split the returned 1d index (intra is a square n_s x n_s matrix)
        s = mod(nn - 1, n_s) + 1
        s2 = ((nn - 1) / n_s) + 1
        do while (s2.eq.s)
          ! intra has decays on the diagonal and transfer rates off-diagonal
          ! so ensure we've picked an off-diagonal element
          nn = rand_nonzero_real_2d(intra_rates)
          s = mod(nn - 1, n_s) + 1
          s2 = ((nn - 1) / n_s) + 1
        end do
        if (intra_rates(s, s2).eq.0.0_CF) then
          write(*, *) "picked a zero rate from nonzero?"
          write(*, *) nn, s, s2
          do n_eff = 1, n_s
            write(*, *) intra_rates(n_eff, :)
          end do
        end if
        cm%rate = intra_rates(s, s2)
        cm%ist = s
        cm%fst = s2
        cm%loss_index = (4 + (s - 1)) * n_s + s2
      case (5)
        ! decay
        nn = rand_nonzero_real_2d(intra_rates)
        s = mod(nn - 1, n_s) + 1
        s2 = ((nn - 1) / n_s) + 1
        do while (s2.ne.s)
          ! intra has decays on the diagonal, so
          ! ensure we've picked a diagonal element
          nn = rand_nonzero_real_2d(intra_rates)
          s = mod(nn - 1, n_s) + 1
          s2 = ((nn - 1) / n_s) + 1
        end do
        cm%rate = intra_rates(s, s2)
        cm%ist = s
        cm%fst = s
        cm%emissive = emissive(s)
        cm%loss_index = (3 * n_s) + s
      case (6)
        ! annihilation
        ! pick a state with a nonzero annihilation rate
        nn = rand_nonzero_real_2d(ann_rates)
        ! split the returned 1d index (ann is a square n_s x n_s matrix)
        s = mod(nn - 1, n_s) + 1
        s2 = ((nn - 1) / n_s) + 1
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
        cm%rate = ann_rates(s, s2)
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
        call possible_moves_and_rates(s, ft)
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
          call possible_moves_and_rates(s, ft)

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

      ! abundance call should go here, at the start of each run
      ! and after random init but before we do anything
      call abundances()

      allocate(gens(n_sites), source=0_CI)
      allocate(anns(n_sites), source=0_CI)
      allocate(gen_rates(n_s), source=0.0_CF)
      allocate(se_rates(n_s), source=0.0_CF)
      allocate(hop_rates(n_s, coord), source=0.0_CF)
      allocate(intra_rates(n_s, n_s), source=0.0_CF)
      allocate(ann_rates(n_s, n_s), source=0.0_CF)

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
      deallocate(se_rates)
      deallocate(hop_rates)
      deallocate(intra_rates)
      deallocate(ann_rates)

    end subroutine do_run

    subroutine mc_deallocations()
      deallocate(pulse)
    end subroutine mc_deallocations

end module mc
