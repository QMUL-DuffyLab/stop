module mc
  use mpi_f08
  use input
  use lattice
  implicit none
  private
  logical :: possible(5) ! 5 possible kinds of moves considered
  real(kind=CF), allocatable, public :: pulse(:)
  
  contains

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
      ! any non-zero population allows for hopping and decay (3 & 4)
      if (any(n_i(i, :).gt.0)) then
        possible(3) = .true.
        possible(4) = .true.
      end if
      if ((sum(n_i(i, :)).gt.1).and.(any(ann.gt.0.0))) then
        possible(5) = .true.
      end if
    end subroutine possible_moves

    subroutine propose(n_i, i, ft)
      integer :: mt
      real :: r
      integer, allocatable :: p_i(:)

      r = 0.0
      call possible_moves(n_i, i, ft)
      ! indices where possible is true
      p_i = pack([(i, i = 1_CI, size(possible))], possible)
      call random_number(r)
      mt = ceiling(r * size(p_i))
      select case (mt)
      case (1)
        ! generation
        ! need to pick a non-zero cross section and make note of
        ! which index that cross section is at
      case (2)
        ! stimulated emission
      case (3)
        ! hop
      case (4)
        ! decay
      case (5)
        ! decay
      case default
        stop
      end if

      
    end subroutine propose
end module mc
