module mc
  use mpi_f08
  use io
  use lattice
  implicit none
  private
  logical(kind=CB) :: possible(6) ! 6 types of moves considered
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

    function rand_nonzero_real(arr) result(i)
      real(kind=CF), intent(in) :: arr(:)
      real(kind=CF) :: r
      integer(kind=CI) :: i
      integer(kind=CI), allocatable :: ind(:)
      ind = pack([(i, i = 1_CI, size(arr))], arr.gt.0.0_CF)
      call random_number(r)
      i = ceiling(r * size(ind))
      if ((i.lt.1_CI).or.(i.gt.size(ind))) then
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
      i = ceiling(r * size(ind))
      if ((i.lt.1_CI).or.(i.gt.size(ind))) then
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
      call random_number(r)
      i = ceiling(r * size(ind))
      if ((i.lt.1_CI).or.(i.gt.size(ind))) then
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
      m%mt = 0
      m%isi = 0
      m%ist = 0
      m%fsi = 0
      m%fst = 0
      m%loss_index = 0
      m%emissive = .false.
      m%rate = 0.0
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
      possible = .false.
      ! if the pulse is on (ft) and there's a positive cross section
      ! then excitations can be generated. in general the second part
      ! should always be true, i guess, otherwise what are we doing here
      if ((ft.gt.0.0).and.(any(xsec.gt.0.0))) then
        possible(1) = .true.
      end if
      ! product of xsec and current populations must be > 0
      ! in order for stimulated emission to occur
      if ((ft.gt.0.0).and.(any(xsec * n_i(site, :).gt.0.0))) then
        possible(2) = .true.
      end if
      ! any non-zero population allows for hopping, transfer, decay
      ! (3, 4, 5)
      if (any(n_i(site, :).gt.0)) then
        possible(3) = .true.
        possible(4) = .true.
        possible(5) = .true.
      end if
      if ((sum(n_i(site, :)).gt.1).and.(any(ann.gt.0.0))) then
        ! annihilation possible (actually this isn't strictly true
        ! and it might need fixing)
        possible(6) = .true.
      end if
    end subroutine possible_moves

    subroutine propose(i, ft)
      integer, intent(in) :: i
      real(kind=CF), intent(in) :: ft 
      integer :: s, s2, p, nn, n_eff

      ! make sure there's no leftover stuff in current_move
      ! that might mess with this
      ! write(*, *) cm%mt
      call zero_move(cm)
      ! write(*, *) cm%mt
      call possible_moves(i, ft)
      ! indices where possible is true
      cm%mt = rand_true(possible)
      ! write(*, *) cm%mt
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
        ! write(*, *) ft, s, xsec(s), p, cm%mt, cm%rate
      case (2)
        ! stimulated emission
        s = rand_nonzero_real(xsec) ! which state is it
        p = which_p(s) ! which pigment is this state on
        cm%rate = ft * xsec(s) * n_i(i, s) / n_tot(p)
        cm%ist = s
        cm%fst = s
        ! this assumes there's only one protein in the system
        ! i.e. it'll need changing if there are multiple proteins
        cm%loss_index = s
        if (cm%loss_index.gt.24) then
          write(*, *) "cm%loss_index is too big (SE)"
          write(*, *) "s = ", s
          write(*, *) "cm%loss_index = ", cm%loss_index
        end if
      case (3)
        ! hop
        ! write(*, *) i, n_i(i, :), possible, cm%mt
        s = rand_nonzero_int(n_i(i, :))
        nn = rand_nonzero_int(neighbours(i, :))
        p = which_p(s)
        cm%rate = n_i(i, s) * hop(s)
        if (n_i(i, s).lt.n_i(nn, s)) then
          cm%rate = cm%rate * (n_i(i, s) * n_thermal(p) - n_i(nn, s)) /&
          ((n_i(nn, s) + 1) * (n_thermal(p) - n_i(i, s) + 1))
        end if
        cm%ist = s
        cm%fsi = nn
        cm%fst = s
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
        cm%loss_index = n_s + s
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
        ! here - come back to this if i end up doing that!
        if (s.gt.s2) then
          n_eff = s
          s = s2
          s2 = n_eff
        end if
        ! write(*, *) "dist", s, s2
        if (dist(s, s2)) then
          cm%rate = n_i(i, s) * n_i(i, s2) * ann(s, s2)
        else
          n_eff = n_i(i, s) + n_i(i, s2)
          cm%rate = ann(s, s2) * (n_eff * (n_eff - 1)) / 2.0
        end if
        cm%loss_index = (2 + (s - 1)) * n_s + s2
        if (cm%loss_index.gt.24) then
          write(*, *) "cm%loss_index is too big (annihilation)"
          write(*, *) "s = ", s
          write(*, *) "s2 = ", s2
          write(*, *) "cm%loss_index = ", cm%loss_index
        end if
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
        a = which_ann(cm%ist, cm%fst)
        n_i(cm%isi, a) = n_i(cm%isi, a) - 1
      case default
        write(*, *) "something wrong in do_cm!"
        stop
      end select
    end subroutine do_move

    subroutine mc_step(t, bin)
      real(kind=CF) :: ft, prob, r, t
      logical(kind=CB) :: bin
      integer(kind=CI) :: pind, i, j, s, n_poss

      n_accepted = 0_CI

      if (t.le.2.0*mu) then
        pind = int(t / dt1) + 1
        if (pind.le.size(pulse)) then
          ft = pulse(pind)
        end if
      end if

      do i = 1, n_sites
        s = rand_int(n_sites)
        call possible_moves(s, ft)

        n_poss = count(possible)
        if (n_poss.eq.0) then
          cycle
        end if

        do j = 1, n_poss
          call possible_moves(s, ft)
          if (.not.any(possible)) then
            cycle
          end if
          call propose(i, ft)
          prob = cm%rate * dt1 * exp(-1.0 * cm%rate * dt1)
          call random_number(r)
          if (r.lt.prob) then
            n_accepted(cm%mt) = n_accepted(cm%mt) + 1
            call do_move()
            if (bin.and.cm%loss_index.gt.0) then
              counts(ceiling(t / binwidth), cm%loss_index) = &
                counts(ceiling(t / binwidth), cm%loss_index) + 1

            end if
          end if
        end do

      end do
    end subroutine mc_step

    subroutine do_run(seed, max_counts, out_file_path)
      integer(kind=CI), intent(in) :: seed(8), max_counts
      real(kind=CF) :: t, interval
      integer(kind=CI) :: i, j, curr_maxcount, rep, tot_acc(6)
      integer(kind=CI), allocatable :: ec(:)
      character(100) :: out_file_path, outfile
      logical :: skip

      call random_seed(put=seed)
      counts = 0_CI

      ec = pack([(i, i = 1_CI, size(emissive_columns))], emissive_columns)

      rep = 1_CI
      curr_maxcount = 0_CI
      interval = 1.0 / rep_rate
      write(*, *) seed

      reploop: do while (curr_maxcount.lt.max_counts)

        t = 0.0_CF
        skip = .false.

        pulseloop: do while (t.lt.tmax)
          call mc_step(t, .true._CB)
          tot_acc = tot_acc + n_accepted
          t = t + dt1
          if ((t.gt.(pulse_tmax)).and.(sum(n_i).eq.0_CI)) then
            skip = .true.
            exit pulseloop ! start the next rep now, nothing to simulate
          end if
        end do pulseloop

        if (.not.skip) then
          darkloop: do while (t.lt.interval)
            call mc_step(t, .false._CB)
            tot_acc = tot_acc + n_accepted
            t = t + dt2
            if (sum(n_i).eq.0_CI) then
              exit darkloop ! start the next rep now, nothing to simulate
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
          ! curr_maxcount = maxval(counts(:, ec))
          write(*,*) seed(1), rep, curr_maxcount
        end if

        rep = rep + 1

      end do reploop

      ! max count reached
      write(*, '(i0, 1X, i0, 1X, i0)') seed(1), rep, curr_maxcount
      write(outfile, '(a, a, I0, a, I0, a)') trim(adjustl(out_file_path)),&
        "_seed_", seed(1), "_rep_", rep, "_final.csv"
      call write_histogram(outfile)

    end subroutine do_run

end module mc
