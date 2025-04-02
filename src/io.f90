module io
  use mpi_f08
  use iso_c_binding
  implicit none
  private
  integer, parameter, public :: CF = c_double
  integer, parameter, public :: CI = c_int
  integer, parameter, public :: CB = c_bool
  character(len=10), public :: protein_name, lattice_name
  character(len=30), allocatable :: labels(:)
  character(len=100), public :: outdir
  integer(kind=CI), public :: n_p, n_s, n_sites, n_bins, n_counts, n_repeats
  real(kind=CF), public :: fwhm, fluence, rep_rate, tmax, dt1, dt2, binwidth
  character(len=10), allocatable, public :: p_names(:), s_names(:)
  integer(kind=CI), allocatable, public :: n_tot(:), n_thermal(:),&
    which_p(:), which_ann(:, :), counts(:, :)
  real(kind=CF), allocatable, public :: abundance(:), hop(:), xsec(:),&
                               intra(:, :), ann(:, :), bins(:)
  logical(kind=CB), allocatable, public :: emissive(:), dist(:, :),&
    emissive_columns(:)
  public :: get_protein_params, get_simulation_params,&
    print_lattice, generate_histogram, write_histogram

  contains

    subroutine get_protein_params(filename)
      ! get the protein parameters from a file. i had originally
      ! set this up with json-fortran but MPI segfaults trying to
      ! run that and i couldn't be bothered to figure out why, so
      ! run some python to set the params in a JSON file and print
      ! them out into a separate file for the fortran at runtime.
      character(len=*), intent(in) :: filename
      integer(kind=CI) :: nunit
      logical(kind=CB), allocatable :: dist_temp(:)
      integer(kind=CI), allocatable :: which_ann_temp(:)
      real(kind=CF), allocatable :: intra_temp(:), ann_temp(:)

      open(newunit=nunit, file=filename)
      read(nunit, *) protein_name
      read(nunit, *) n_p
      read(nunit, *) n_s
      
      allocate(p_names(n_p))
      allocate(s_names(n_s))
      allocate(which_p(n_s))
      allocate(abundance(n_s))
      allocate(dist_temp(n_s * n_s))
      allocate(dist(n_s, n_s))
      allocate(intra_temp(n_s * n_s))
      allocate(intra(n_s, n_s))
      allocate(ann_temp(n_s * n_s))
      allocate(ann(n_s, n_s))
      allocate(which_ann_temp(n_s * n_s))
      allocate(which_ann(n_s, n_s))
      allocate(n_tot(n_p))
      allocate(n_thermal(n_p))
      allocate(hop(n_s))
      allocate(emissive(n_s))
      allocate(xsec(n_s))

      read(nunit, *) p_names
      read(nunit, *) s_names
      read(nunit, *) which_p
      read(nunit, *) abundance
      read(nunit, *) dist_temp
      read(nunit, *) n_tot
      read(nunit, *) n_thermal
      read(nunit, *) hop
      read(nunit, *) intra_temp
      read(nunit, *) ann_temp
      read(nunit, *) which_ann_temp
      read(nunit, *) xsec
      read(nunit, *) emissive

      dist      = reshape(dist_temp,      (/ n_s, n_s /))
      intra     = reshape(intra_temp,     (/ n_s, n_s /))
      ann       = reshape(ann_temp,       (/ n_s, n_s /))
      which_ann = reshape(which_ann_temp, (/ n_s, n_s /))

      hop = 1.0_CF / hop
      intra = 1.0_CF / intra
      ann = 1.0_CF / ann
    end subroutine get_protein_params

    subroutine get_simulation_params(filename)
     ! get the simulation parameters from a file
     character(len=*), intent(in) :: filename
     integer(kind=CI) :: nunit

     open(newunit=nunit, file=filename)
     read(nunit, *) fwhm
     read(nunit, *) fluence
     read(nunit, *) n_sites
     read(nunit, *) lattice_name
     read(nunit, *) rep_rate
     read(nunit, *) tmax
     read(nunit, *) dt1
     read(nunit, *) dt2
     read(nunit, *) binwidth
     read(nunit, *) n_counts
     read(nunit, *) n_repeats
     read(nunit, '(a)') outdir

    end subroutine get_simulation_params

    subroutine print_lattice(filename, coords, neighbours)
      character(len=*) :: filename
      integer, intent(in) :: coords(:, :)
      integer, intent(in) :: neighbours(:, :)
      integer :: nunit, i
      open(newunit=nunit, file=filename)
      do i = 1, n_sites
        write(nunit, *) i, coords(i, :), neighbours(i, :)
      end do
      close(nunit)
    end subroutine print_lattice

    subroutine generate_histogram()
      integer :: n_losses, i, j, loss_index
      n_losses = n_s * (2 + n_s)
      n_bins = ceiling(tmax / binwidth)
      allocate(counts(n_bins, n_losses), source=0_CI)
      allocate(labels(n_losses + 1))
      allocate(bins(n_bins), source=0.0_CF)
      allocate(emissive_columns(n_losses), source=.false._CB)

      write(labels(1), '(a)') "Time (s)"

      do i = 1, n_bins
        bins(i) = (i - 1) * binwidth
      end do

      loss_index = 2
      do i = 1, n_s
        write(labels(loss_index), '(a, a, a)') trim(adjustl(protein_name)),&
          "_se_",&
          trim(adjustl(s_names(i)))
        loss_index = loss_index + 1
      end do

      do i = 1, n_s
        write(labels(loss_index), '(a, a, a)') trim(adjustl(protein_name)),&
          "_decay_",&
          trim(adjustl(s_names(i)))
        if (emissive(i)) then
          emissive_columns(loss_index - 1) = .true.
        end if
        loss_index = loss_index + 1
      end do

      do i = 1, n_s
        do j = 1, n_s
          write(labels(loss_index), '(a, a, a, a, a)') trim(adjustl(protein_name)),&
            "_ann_",&
            trim(adjustl(s_names(i))), "_", trim(adjustl(s_names(j)))
          loss_index = loss_index + 1
        end do
      end do

    end subroutine generate_histogram

    subroutine write_histogram(filename)
      character(len=*) :: filename
      integer :: nunit, i
      character(len=30) :: str_fmt

      write(str_fmt, '(a, i0, a)') "(ES10.4, ", n_s * (n_s + 2), "(1X, I0))"
      open(newunit=nunit, file=filename)
      write(nunit, *) labels
      write(nunit, *) emissive_columns
      do i = 1, n_bins
        write(nunit, str_fmt) bins(i), counts(i, :)
      end do
      close(nunit)
    end subroutine write_histogram

end module io
