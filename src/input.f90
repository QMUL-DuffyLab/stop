module input
  use mpi_f08
  use iso_c_binding
  implicit none
  private
  integer, parameter, public :: CF = c_double
  integer, parameter, public :: CI = c_int
  integer, parameter, public :: CB = c_bool
  character(len=10), public :: protein_name, lattice_name
  character(len=100), public :: outdir
  integer(kind=CI), public :: n_p, n_s, n_sites, n_counts, n_repeats
  real(kind=CF), public :: fwhm, fluence, rep_rate, tmax, dt1, dt2, binwidth
  character(len=10), allocatable, public :: p_names(:), s_names(:)
  integer(kind=CI), allocatable, public :: n_tot(:), n_thermal(:),&
    which_p(:), which_ann(:, :)
  real(kind=CF), allocatable, public :: abundance(:), hop(:), xsec(:),&
                               intra(:, :), ann(:, :)
  logical(kind=CB), allocatable, public :: emissive(:), dist(:, :)
  public :: get_protein_params, get_simulation_params

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
      read(nunit, *) emissive
      read(nunit, *) xsec

      dist      = reshape(dist_temp,      (/ n_s, n_s /))
      intra     = reshape(intra_temp,     (/ n_s, n_s /))
      ann       = reshape(ann_temp,       (/ n_s, n_s /))
      which_ann = reshape(which_ann_temp, (/ n_s, n_s /))

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

end module
