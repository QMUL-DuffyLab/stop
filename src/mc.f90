module mc
  use mpi_f08
  use io
  use lattice
  implicit none
  private
  logical :: possible(5) ! 5 possible kinds of moves considered
  real(kind=CF), allocatable, public :: pulse(:)
  type :: move_type
    integer :: mt = 0
    integer :: isi = 0
    integer :: ist = 0
    integer :: fsi = 0
    integer :: fst = 0
    integer :: loss_index = 0
    logical :: emissive = .false.
    real :: rate = 0.0
  end type
  type(move_type), public :: cm ! current move
  
  contains

    function rand_nonzero(arr) result(i)
      real, intent(in) :: arr(:)
      real :: r
      integer :: i, ind(:)
      ind = pack([(i, i = 1_CI, size(arr))], arr.gt.0.0)
      call random_number(r)
      i = ceiling(r * size(arr))
    end function rand_nonzero

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

    function rand_int(imax) result(i)
      ! pick an integer from 1 to imax inclusive
      integer, intent(in) :: imax
      real :: r
      integer :: i
      call random_number(r)
      i = ceiling(r * imax)
    end function rand_nonzero

    function rand_true(arr) result(i)
      logical, intent(in) :: arr(:)
      real :: r
      integer :: i, ind(:)
      ind = pack([(i, i = 1_CI, size(arr))], arr)
      call random_number(r)
      i = ceiling(r * size(arr))
    end function rand_true

    function construct_pulse(mu, fwhm, dt, n_per_t,&
        n_sites) result(pulse)
      implicit none
      real(kind=CF) :: mu, sigma, pulse_tmax
      integer :: i, pulse_len
      mu = 2.0 * fwhm
      pulse_tmax = 2.0 * mu
      pulse_len = int(pulse_tmax / dt1)
      sigma = fwhm / (2.0_dp * (sqrt(2.0_dp * log(2.0_dp))))
      allocate(pulse(pulse_len))
      do i = 1, pulse_len
        pulse(i) = (fluence) / (sigma * sqrt(2.0_dp * pi)) * &
          exp(-1.0_dp * ((i * dt1) - mu)**2 / (sqrt(2.0_dp) * sigma)**2)
      end do
    end function construct_pulse

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
      ! any non-zero population allows for hopping, decay, transfer
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
      integer :: s, s2, p, nn, n_eff

      ! make sure there's no leftover stuff in current_move
      ! that might mess with this
      call zero_move(cm)
      call possible_moves(i, ft)
      ! indices where possible is true
      cm%mt = rand_true(possible)
      cm%isi = i
      cm%fsi = i

      select case (cm%mt)
      case (1)
        ! generation
        ! need to pick a non-zero cross section and make note of
        ! which index that cross section is at
        s = rand_nonzero(xsec) ! which state is it
        p = which_p(s) ! which pigment is this state on
        cm%rate = ft * xsec(s) * (n_tot(p) - n_i(i, s)) / n_tot(p)
        cm%ist = s
        cm%fst = s
      case (2)
        ! stimulated emission
        s = rand_nonzero(xsec) ! which state is it
        p = which_p(s) ! which pigment is this state on
        cm%rate = ft * xsec(s) * n_i(i, s) / n_tot(p)
        cm%ist = s
        cm%fst = s
        ! this assumes there's only one protein in the system
        ! i.e. it'll need changing if there are multiple proteins
        cm%loss_index = s
      case (3)
        ! hop
        s = rand_nonzero(n_i(i, :))
        nn = rand_nonzero(neighbours(i, :))
        p = which_p(s)
        cm%rate = n_i(i, s) * hop(s)
        if (n_i(i, s).lt.n(nn, s)) then
          cm%rate = cm%rate * (n_i(i, s) * n_thermal(p) - n_i(nn, s)) /&
          ((n_i(nn, s) + 1) * (n_thermal(p) - n_i(i, s) + 1))
        end if
        cm%ist = s
        cm%fsi = nn
        cm%fst = s
      case (4)
        ! transfer
        s = rand_nonzero(n_i(i, :))
        s2 = rand_int(n_states)
        do while (s2.eq.s)
          s2 = rand_int(n_states)
        end do
        cm%rate = n_i(i, s) * intra(s, s2)
        cm%ist = s
        cm%fst = s2
      case (5)
        ! decay
        s = rand_nonzero(n_i(i, :))
        cm%rate = n_i(i, s) * intra(s, s)
        cm%ist = s
        cm%fst = s
        cm%emissive = emissive(s)
        cm%loss_index = n_states + s
        end if
      case (6)
        ! annihilation
        s = rand_nonzero(n_i(i, :))
        ! pick a state with a nonzero annihilation rate with s
        s2 = rand_nonzero(ann(s, :))
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

        if (dist(s, s2)) then
          cm%rate = n(i, s) * n(i, s2) * ann(s, s2)
        else
          n_eff = n(i, s) + n(i, s2)
          cm%rate = ann(s, s2) * (n_eff * (n_eff - 1)) / 2.0
        end if
        cm%loss_index = (2 + (s - 1)) * n_states + s2
        ist = s
        fst = s2
      case default
        write(*, *) "shouldn't be here! cm type wrong in propose()"
        stop
      end if
    end function propose
    
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

    subroutine mc_step(pulse, t)
      real :: ft, prob, r, t
      integer :: pind, i, j, s, n_poss

      if (t.le.2.0*mu) then
        pind = ceiling(t / dt1)
        if (pind.le.size(pulse)) then
          ft = pulse(pind)
        end if
      end if

      do i = 1, n_sites
        s = rand_int(n_sites)
        call possible_moves(s, ft)
        if .not.any(possible) then
          continue
        end if
        n_poss = count(possible)
        do j = 1, n_poss
          call possible_moves(s, ft)
          if .not.any(possible) then
            continue
          end if
          call propose(i, ft)
          prob = cm%rate * dt1 * exp(-1.0 * cm%rate * dt1)
          call random_number(r)
          if (r.lt.prob) then
            call do_move()
            if (bin.and.cm%loss_index.gt.0) then
              counts(loss_index, ceiling(t / binwidth)) = &
                counts(loss_index, ceiling(t / binwidth)) + 1
          end if
        end do

      end do
    end subroutine mc_step

end module mc
