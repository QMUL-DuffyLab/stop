module mc
  use mpi_f08
  use io
  use lattice
  implicit none
  private
  logical(kind=CB) :: possible(6) ! 6 types of moves considered
  logical(kind=CB) :: bad_proposal
  integer(kind=CI), public :: n_accepted(6)
  real(kind=CF), allocatable, public :: pulse(:)
  real(kind=CF), public :: mu, pulse_tmax
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
  type(move_type), public :: cm ! current move
  public :: construct_pulse, do_run
  
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

    function rand_nonzero_int(arr) result(i)
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

    function rand_true(arr) result(i)
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
      possible = .false._CB
      ! if the pulse is on (ft) and there's a positive cross section
      ! then excitations can be generated. in general the second part
      ! should always be true, i guess, otherwise what are we doing here
      if ((ft.gt.0.0).and.(any(xsec.gt.0.0))) then
        possible(1) = .true._CB
      end if
      ! product of xsec and current populations must be > 0
      ! in order for stimulated emission to occur
      if ((ft.gt.0.0).and.(any(xsec * n_i(site, :).gt.0.0))) then
        possible(2) = .true._CB
      end if
      ! any non-zero population allows for hopping, and decay
      ! (3, 5)
      if (any(n_i(site, :).gt.0)) then
        possible(3) = .true._CB
        possible(5) = .true._CB
      end if
      if ((n_s.gt.1).and.(any(n_i(site, :).gt.0))) then
        ! for transfer between states to be possible, there
        ! must be more than one state on the protein
        possible(4) = .true._CB
      end if
      if ((sum(n_i(site, :)).gt.1).and.(any(ann.gt.0.0))) then
        ! annihilation possible (actually this isn't strictly true
        ! and it might need fixing)
        possible(6) = .true._CB
      end if
    end subroutine possible_moves

    subroutine propose(i, ft)
      integer, intent(in) :: i
      real(kind=CF), intent(in) :: ft 
      integer :: s, s2, p, nn, n_eff

      ! make sure there's no leftover stuff in current_move
      ! that might mess with this
      call zero_move(cm)
      call possible_moves(i, ft)
      ! indices where possible is true
      if (.not.any(possible)) then
        ! something has gone badly wrong somewhere
        bad_proposal = .true.
        return
      end if

      cm%mt = rand_true(possible)
      cm%isi = i
      cm%fsi = i

      select case (cm%mt)
      case (1)
        ! generation
        ! need to pick a non-zero cross section and make note of
        ! which index that cross section is at
        s = rand_nonzero_real(xsec) ! which state is it
        p = which_p(s) ! which pigment is this state on
        cm%rate = ft * xsec(s) * (n_tot(p) - n_i(i, s)) / n_tot(p)
        cm%ist = s
        cm%fst = s
        cm%loss_index = s
      case (2)
        ! stimulated emission
        s = rand_nonzero_real(xsec) ! which state is it
        p = which_p(s) ! which pigment is this state on
        cm%rate = ft * xsec(s) * n_i(i, s) / n_tot(p)
        cm%ist = s
        cm%fst = s
        ! this assumes there's only one protein in the system
        ! i.e. it'll need changing if there are multiple proteins
        cm%loss_index = n_s + s
      case (3)
        ! hop
        s = rand_nonzero_int(n_i(i, :))
        nn = rand_nonzero_int(neighbours(i, :))
        p = which_p(s)
        cm%rate = n_i(i, s) * hop(s)
        if (n_i(i, s).lt.n_i(nn, s)) then
          ! entropic penalty for hopping to a higher-occupied site
          cm%rate = cm%rate * (n_i(i, s) * (n_thermal(p) - n_i(nn, s))) /&
          ((n_i(nn, s) + 1) * (n_thermal(p) - n_i(i, s) + 1))
        end if
        cm%ist = s
        cm%fsi = nn
        cm%fst = s
        cm%loss_index = (3 * n_s) + s
      case (4)
        ! transfer
        s = rand_nonzero_int(n_i(i, :))
        s2 = rand_int(n_s)
        do while (s2.eq.s)
          s2 = rand_int(n_s)
        end do
        cm%rate = n_i(i, s) * intra(s, s2)
        cm%ist = s
        cm%fst = s2
      case (5)
        ! decay
        s = rand_nonzero_int(n_i(i, :))
        cm%rate = n_i(i, s) * intra(s, s)
        cm%ist = s
        cm%fst = s
        cm%emissive = emissive(s)
        cm%loss_index = (2 * n_s) + s
      case (6)
        ! annihilation
        s = rand_nonzero_int(n_i(i, :))
        ! pick a state with a nonzero annihilation rate with s
        s2 = rand_nonzero_real(ann(s, :))
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
        if (dist(s, s2)) then
          cm%rate = n_i(i, s) * n_i(i, s2) * ann(s, s2)
        else
          n_eff = n_i(i, s) + n_i(i, s2)
          cm%rate = ann(s, s2) * (n_eff * (n_eff - 1)) / 2.0
        end if
        cm%loss_index = (4 + (s - 1)) * n_s + s2
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
        call possible_moves(s, ft)

        n_poss = count(possible)
        if (n_poss.eq.0) then
          cycle
        end if

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
      integer(kind=CI) :: i, j, curr_maxcount, rep, nunit
      integer(kind=CI), allocatable :: ec(:)
      character(100) :: out_file_path, outfile, pop_file
      logical(kind=CB) :: skip, bin_pulse

      call init_random(salt)
      counts = 0_CI
      tot_accepted = 0_CI

      ! first column in emissive columns isn't a real one -
      ! it's for the output file. could maybe replace this
      ! with just emissive??
      ec = pack([(i, i = 1_CI, size(emissive_columns) - 1)], emissive_columns(2:))

      rep = 1_CI
      curr_maxcount = 0_CI
      interval = 1.0 / rep_rate

      write(pop_file, '(a, a, I0, a, I0, a)') trim(adjustl(out_file_path)),&
        "_salt_", salt, "_population.csv"
      open(newunit=nunit, file=pop_file)

      reploop: do while (curr_maxcount.lt.max_counts)

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
          write (*, *) salt, rep, curr_maxcount, ec
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
      write(outfile, '(a, a, I0, a, I0, a)') trim(adjustl(out_file_path)),&
        "_salt_", salt, "_rep_", rep, "_final.csv"
      call write_histogram(outfile)

    end subroutine do_run

end module mc
