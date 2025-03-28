module mc
  use mpi_f08
  use io
  use lattice
  implicit none
  private
  logical :: possible(5) ! 5 possible kinds of moves considered
  real(kind=CF), allocatable, public :: pulse(:)
  type :: move
    integer :: mt = 0
    integer :: isi = 0
    integer :: ist = 0
    integer :: fsi = 0
    integer :: fst = 0
    integer :: loss_index = 0
    logical :: emissive = .false.
    real :: rate = 0.0
  end type
  
  contains

    function rand_nonzero(arr) result(i)
      real, intent(in) :: arr(:)
      real :: r
      integer :: i, ind(:)
      ind = pack([(i, i = 1_CI, size(arr))], arr.gt.0.0)
      call random_number(r)
      i = ceiling(r * size(arr))
    end function rand_nonzero

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

    subroutine possible_moves(n_i, i, ft)
      integer(kind=CI), allocatable, intent(in) :: n_i(:, :)
      real(kind=CF) :: ft 
      possible = .false.
      ! if the pulse is on (ft) and there's a positive cross section
      ! then excitations can be generated. in general the second part
      ! should always be true, i guess, otherwise what are we doing here
      if ((ft.gt.0.0).and.(any(xsec.gt.0.0))) then
        possible(1) = .true.
      end if
      ! product of xsec and current populations must be > 0
      ! in order for stimulated emission to occur
      if ((ft.gt.0.0).and.(any(xsec * n_i(i, :).gt.0.0))) then
        possible(2) = .true.
      end if
      ! any non-zero population allows for hopping, decay, transfer
      ! (3, 4, 5)
      if (any(n_i(i, :).gt.0)) then
        possible(3) = .true.
        possible(4) = .true.
        possible(5) = .true.
      end if
      if ((sum(n_i(i, :)).gt.1).and.(any(ann.gt.0.0))) then
        ! annihilation possible (actually this isn't strictly true
        ! and it might need fixing)
        possible(6) = .true.
      end if
    end subroutine possible_moves

    function propose(n_i, i, ft) result(move)
      integer :: s, s2, p, nn, n_eff
      type(move) :: move

      call possible_moves(n_i, i, ft)
      ! indices where possible is true
      move%mt = rand_true(possible)
      move%isi = i
      move%fsi = i

      select case (move%mt)
      case (1)
        ! generation
        ! need to pick a non-zero cross section and make note of
        ! which index that cross section is at
        s = rand_nonzero(xsec) ! which state is it
        p = which_p(s) ! which pigment is this state on
        move%rate = ft * xsec(s) * (n_tot(p) - n_i(i, s)) / n_tot(p)
        move%ist = s
        move%fst = s
      case (2)
        ! stimulated emission
        s = rand_nonzero(xsec) ! which state is it
        p = which_p(s) ! which pigment is this state on
        move%rate = ft * xsec(s) * n_i(i, s) / n_tot(p)
        move%ist = s
        move%fst = s
        ! this assumes there's only one protein in the system
        ! i.e. it'll need changing if there are multiple proteins
        move%loss_index = s
      case (3)
        ! hop
        s = rand_nonzero(n_i(i, :))
        nn = rand_nonzero(neighbours(i, :))
        p = which_p(s)
        move%rate = n_i(i, s) * hop(s)
        if (n_i(i, s).lt.n(nn, s)) then
          move%rate = move%rate * (n_i(i, s) * n_thermal(p) - n_i(nn, s)) /&
          ((n_i(nn, s) + 1) * (n_thermal(p) - n_i(i, s) + 1))
        end if
        move%ist = s
        move%fsi = nn
        move%fst = s
      case (4)
        ! transfer
        s = rand_nonzero(n_i(i, :))
        s2 = rand_int(n_states)
        do while (s2.eq.s)
          s2 = rand_int(n_states)
        end do
        move%rate = n_i(i, s) * intra(s, s2)
        move%ist = s
        move%fst = s2
      case (5)
        ! decay
        s = rand_nonzero(n_i(i, :))
        move%rate = n_i(i, s) * intra(s, s)
        move%ist = s
        move%fst = s
        move%emissive = emissive(s)
        move%loss_index = n_states + s
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
          move%rate = n(i, s) * n(i, s2) * ann(s, s2)
        else
          n_eff = n(i, s) + n(i, s2)
          move%rate = ann(s, s2) * (n_eff * (n_eff - 1)) / 2.0
        end if
        move%loss_index = (2 + (s - 1)) * n_states + s2
        ist = s
        fst = s2
      case default
        write(*, *) "shouldn't be here! move type wrong in propose()"
        stop
      end if
    end function propose
    
    subroutine do_move(move)
      integer :: a
      select case(move%mt)
      case (1) ! gen
        n_i(move%isi, move%ist) = n_i(move%isi, move%ist) + 1
      case (2) ! se
        n_i(move%isi, move%ist) = n_i(move%isi, move%ist) - 1
      case (3) ! hop
        n_i(move%isi, move%ist) = n_i(move%isi, move%ist) - 1
        n_i(move%fsi, move%fst) = n_i(move%fsi, move%fst) + 1
      case (4) ! transfer
        n_i(move%isi, move%ist) = n_i(move%isi, move%ist) - 1
        n_i(move%fsi, move%fst) = n_i(move%fsi, move%fst) + 1
      case (5) ! decay
        n_i(move%isi, move%ist) = n_i(move%isi, move%ist) - 1
      case (6) ! annihilation
        a = which_ann(move%ist, move%fst)
        n_i(move%isi, a) = n_i(move%isi, a) - 1
      case default
        write(*, *) "something wrong in do_move!"
        stop
      end select
    end subroutine do_move

    subroutine mc_step(pulse, t)
      real :: ft
    end subroutine mc_step

end module mc
