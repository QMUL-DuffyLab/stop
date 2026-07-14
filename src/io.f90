module io
  use mpi_f08
  use iso_c_binding
  use ieee_arithmetic
  implicit none
  private
  integer, parameter, public :: CF = c_double
  integer, parameter, public :: CI = c_int
  integer, parameter, public :: CB = c_bool
  integer, parameter, public :: hist_max = 10
  character(len=100), public :: protein_name, lattice_name
  character(len=50), allocatable :: labels(:)
  character(len=200), public :: outdir
  logical, public :: debug
  integer(kind=CI), public :: n_p, n_s, n_sites, n_bins,&
    n_counts, n_repeats, burn_reps
  real(kind=CF), public :: fwhm, fluence, rep_rate, tmax, dt1, dt2, binwidth
  character(len=10), allocatable, public :: p_names(:), s_names(:)
  integer(kind=CI), allocatable, public :: n_tot(:), n_thermal(:),&
    which_p(:), ann_remainder(:, :), counts(:, :), site_moves(:, :),&
    site_gen_hist(:, :), site_ann_hist(:, :)
  real(kind=CF), allocatable, public :: abundance(:), hop(:), xsec(:),&
    intra(:, :), ann(:, :), bins(:)
  logical(kind=CB), allocatable, public :: emissive(:), dist(:, :),&
    emissive_columns(:)
  public :: get_protein_params, get_simulation_params,&
    print_lattice, generate_histogram, write_histogram, write_site_moves,&
    write_move_hists, io_deallocations, labels

  contains

    subroutine get_protein_params(filename)
      ! get the protein parameters from a file. i had originally
      ! set this up with json-fortran but MPI segfaults trying to
      ! run that and i couldn't be bothered to figure out why, so
      ! run some python to set the params in a JSON file and print
      ! them out into a separate file for the fortran at runtime.
      character(len=*), intent(in) :: filename
      integer(kind=CI) :: nunit, i, j
      logical(kind=CB), allocatable :: dist_temp(:)
      integer(kind=CI), allocatable :: ann_rem_temp(:)
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
      allocate(ann_rem_temp(n_s * n_s))
      allocate(ann_remainder(n_s, n_s))
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
      read(nunit, *) ann_rem_temp
      read(nunit, *) xsec
      read(nunit, *) emissive

      write(*, *) "reading in rates"
      write(*, *) "----------------"
      write(*, *)
      write(*, *) "hop = ", hop
      write(*, *) "intra_temp = ", intra_temp
      write(*, *) "ann_temp = ", ann_temp
      write(*, *) "ann_rem_temp = ", ann_rem_temp

      ! fortran is column-major
      dist          = transpose(reshape(dist_temp,    (/ n_s, n_s /)))
      intra         = transpose(reshape(intra_temp,   (/ n_s, n_s /)))
      ann           = transpose(reshape(ann_temp,     (/ n_s, n_s /)))
      ann_remainder = transpose(reshape(ann_rem_temp, (/ n_s, n_s /)))

      ! parameters given in protein.json as time constants;
      ! convert them to rates in s^{-1}. be careful not to do 1./0
      do i = 1, n_s
        if (hop(i).eq.0.0_CF) then
          hop(i) = 0.0_CF
        else
          hop(i) = 1.0_CF / hop(i)
        end if
        do j = 1, n_s
          if (intra(i, j).eq.0.0_CF) then
            intra(i, j) = 0.0_CF
          else
            intra(i, j) = 1.0_CF / intra(i, j)
          end if
          if (ann(i, j).eq.0.0_CF) then
            ann(i, j) = 0.0_CF
          else
            ann(i, j) = 1.0_CF / ann(i, j)
          end if
        end do
      end do

      write(*, *)
      write(*, *) "hop (rates)"
      write(*, *) hop
      write(*, *) "intra (rates)"
      do i = 1, n_s
        write(*, *) intra(i, :)
      end do
      write(*, *)
      write(*, *) "ann (rates)"
      do i = 1, n_s
        write(*, *) ann(i, :)
      end do
      write(*, *)
      write(*, *) "ann remainder"
      do i = 1, n_s
        write(*, *) ann_remainder(i, :)
      end do
      write(*, *) "xsec"
      write(*, *) xsec

      deallocate(dist_temp)
      deallocate(intra_temp)
      deallocate(ann_temp)
      deallocate(ann_rem_temp)
    end subroutine get_protein_params

    subroutine get_simulation_params(filename)
      ! get the simulation parameters from a file
      character(len=*), intent(in) :: filename
      integer(kind=CI) :: nunit
      character(len=10) :: dstr

      open(newunit=nunit, file=filename)
      read(nunit, *) fwhm
      read(nunit, *) fluence
      read(nunit, *) n_sites
      read(nunit, *) lattice_name
      read(nunit, *) rep_rate
      read(nunit, *) burn_reps
      read(nunit, *) tmax
      read(nunit, *) dt1
      read(nunit, *) dt2
      read(nunit, *) binwidth
      read(nunit, *) n_counts
      read(nunit, *) n_repeats
      read(nunit, '(a)') outdir
      read(nunit, '(a)') dstr 
      if (dstr.eq.'F'.or.dstr.eq.'False'.or.dstr.eq.'false') then
        debug = .false.
      else if (dstr.eq.'T'.or.dstr.eq.'True'.or.dstr.eq.'true') then
        debug = .true.
      else
        write(*, *) "debug string not parsed in simulation.json"
        stop
      end if

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
      n_losses = n_s * (4 + 2 * n_s)
      n_bins = ceiling(tmax / binwidth)
      allocate(counts(n_bins, n_losses), source=0_CI)
      allocate(site_moves(n_sites, n_losses), source=0_CI)
      allocate(site_gen_hist(n_sites, hist_max), source=0_CI)
      allocate(site_ann_hist(n_sites, hist_max), source=0_CI)
      allocate(labels(n_losses + 1))
      allocate(bins(n_bins), source=0.0_CF)
      allocate(emissive_columns(n_losses + 1), source=.false._CB)

      write(labels(1), '(a)') "Time(s)"

      do i = 1, n_bins
        bins(i) = (i - 1) * binwidth
      end do

      loss_index = 2
      do i = 1, n_s
        write(labels(loss_index), '(a, a)') "gen_", trim(adjustl(s_names(i)))
        loss_index = loss_index + 1
      end do

      do i = 1, n_s
        write(labels(loss_index), '(a, a)') "se_", trim(adjustl(s_names(i)))
        loss_index = loss_index + 1
      end do

      do i = 1, n_s
        write(labels(loss_index), '(a, a)') "hop_", trim(adjustl(s_names(i)))
        loss_index = loss_index + 1
      end do

      do i = 1, n_s
        write(labels(loss_index), '(a, a)') "decay_", trim(adjustl(s_names(i)))
        if (emissive(i)) then
          emissive_columns(loss_index) = .true.
        end if
        loss_index = loss_index + 1
      end do

      do i = 1, n_s
        do j = 1, n_s
          write(labels(loss_index), '(a, a, a, a)') "transfer_",&
            trim(adjustl(s_names(i))), "_", trim(adjustl(s_names(j)))
          loss_index = loss_index + 1
        end do
      end do

      do i = 1, n_s
        do j = 1, n_s
          write(labels(loss_index), '(a, a, a, a)') "ann_",&
            trim(adjustl(s_names(i))), "_", trim(adjustl(s_names(j)))
          loss_index = loss_index + 1
        end do
      end do

    end subroutine generate_histogram

    subroutine write_histogram(filename)
      character(len=*) :: filename
      integer :: nunit, i
      character(len=30) :: str_fmt

      open(newunit=nunit, file=filename)

      write(str_fmt, '(a, i0, a)') "(", size(labels), "(a, 1X))"
      write(nunit, str_fmt) (trim(adjustl(labels(i))), i=1,size(labels))

      write(str_fmt, '(a, i0, a)') "(", size(labels), "(L1, 1X))"
      write(nunit, str_fmt) (emissive_columns(i), i=1, size(emissive_columns))

      write(str_fmt, '(a, i0, a)') "(ES10.4, ", size(labels), "(1X, I0))"
      do i = 1, n_bins
        write(nunit, str_fmt) bins(i), counts(i, :)
      end do
      close(nunit)
    end subroutine write_histogram

    subroutine write_site_moves(filename)
      character(len=*) :: filename
      integer :: nunit, i
      character(len=30) :: str_fmt

      open(newunit=nunit, file=filename)

      write(str_fmt, '(a, i0, a)') "(", size(labels), "(a, 1X))"
      write(nunit, str_fmt) (trim(adjustl(labels(i))), i=1,size(labels))

      write(str_fmt, '(a, i0, a)') "(", size(labels), "(L1, 1X))"
      write(nunit, str_fmt) (emissive_columns(i), i=1, size(emissive_columns))

      write(str_fmt, '(a, i0, a)') "(i0, ", size(labels), "(1X, I0))"
      do i = 1, n_sites
        write(nunit, str_fmt) i, site_moves(i, :)
      end do
      close(nunit)
    end subroutine write_site_moves

    subroutine write_move_hists(basename)
      character(len=*) :: basename
      character(len=200) :: gen_filename, ann_filename
      integer :: genunit, annunit, i
      character(len=30) :: str_fmt

      write(gen_filename, '(a, a)') trim(adjustl(basename)), "gen_hist.csv"
      write(ann_filename, '(a, a)') trim(adjustl(basename)), "ann_hist.csv"

      open(newunit=genunit, file=trim(adjustl(gen_filename)))
      open(newunit=annunit, file=trim(adjustl(ann_filename)))

      write(str_fmt, '(a, i0, a)') "(i0, ", hist_max, "(1X, I0))"
      do i = 1, n_sites
        write(genunit, str_fmt) i, site_gen_hist(i, :)
        write(annunit, str_fmt) i, site_ann_hist(i, :)
      end do
      close(genunit)
      close(annunit)
    end subroutine write_move_hists

    subroutine io_deallocations()
      deallocate(p_names)
      deallocate(s_names)
      deallocate(which_p)
      deallocate(abundance)
      deallocate(dist)
      deallocate(intra)
      deallocate(ann)
      deallocate(ann_remainder)
      deallocate(n_tot)
      deallocate(n_thermal)
      deallocate(hop)
      deallocate(emissive)
      deallocate(xsec)
      deallocate(counts)
      deallocate(site_moves)
      deallocate(labels)
      deallocate(bins)
      deallocate(emissive_columns)
      deallocate(site_gen_hist)
      deallocate(site_ann_hist)
    end subroutine io_deallocations

end module io
