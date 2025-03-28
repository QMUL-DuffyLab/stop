module lattice
  use iso_c_binding
  use mpi_f08
  use io
  implicit none
  private
  real(kind=CF), parameter, public :: pi = 3.1415926535
  real(kind=CF), parameter :: atol = 1.0e-6
  integer(kind=CI), public :: coord
  integer(kind=CI), allocatable, public :: n_i(:, :), neighbours(:, :)
  integer(kind=CI), allocatable, public :: coords(:, :)
  public :: generate_lattice
  
  contains

    subroutine generate_lattice(lattice_name, nmax)
      character(len=10), intent(in) :: lattice_name
      integer(kind=CI), intent(in) :: nmax
      integer(kind=CI), dimension(:, :), allocatable :: lv
      integer(kind=CI) :: n_current, i, j, k, m, coord, ri(2)
      logical(kind=CB) :: add

      select case (lattice_name)
      case ("line")
        coord = 2
        allocate(lv(coord, 2))
        lv(1, :) = (/1, 0/)
        lv(2, :) = (/-1, 0/)
      case ("square")
        coord = 4
        allocate(lv(coord, 2))
        lv(1, :) = (/1, 0/)
        lv(2, :) = (/-1, 0/)
        lv(3, :) = (/0, 1/)
        lv(4, :) = (/0, -1/)
      case ("hex")
        coord = 6
        allocate(lv(coord, 2))
        lv(1, :) = (/1, 0/)
        lv(2, :) = (/-1, 0/)
        lv(3, :) = (/0, 1/)
        lv(4, :) = (/0, -1/)
        lv(5, :) = (/1, -1/)
        lv(6, :) = (/-1, 1/)
      case default
        coord = 0
        write(*, *) "invalid lattice name"
        stop
      end select

      allocate(neighbours(nmax, coord))
      allocate(n_i(nmax, n_s))
      allocate(coords(nmax, 2))
      neighbours = 0_CI
      n_i = 0_CI
      n_current = 1_CI
      coords = 0_CI
      do while (n_current.lt.nmax)
        ! hmmmmmmm. looping to nmax every time might not be necessary
        do i = 1, nmax
          do j = 1, coord
             ri = coords(i, :) + lv(j, :)
             add = .true.
             do k = 1, n_current
               if (all(coords(k, :) - ri.eq.0)) then
                 ! already in the lattice
                 add = .false.
                 exit
               end if
             end do
               if ((add).and.(n_current.lt.nmax)) then
                 n_current = n_current + 1
                 coords(n_current, :) = ri
               end if
          end do
        end do
      end do

      ! must be a better way of doing this lol
      do i = 1, nmax
        m = 1
        do j = 1, nmax
          do k = 1, coord
            if (all((coords(i, :) - coords(j, :)).eq.lv(k, :))) then
              neighbours(i, m) = j
              m = m + 1
            end if
          end do
        end do
      end do

    end subroutine generate_lattice

end module lattice
