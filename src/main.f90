program main
  use mpi_f08
  use io
  use lattice
  use mc
  implicit none
  character(100) :: protein_file, simulation_file, latt_file, hist_file
  integer(kind=CI) :: num_procs, rank, mpierr, run
  real(kind=CF) :: t
  logical(kind=CB) :: bin

  call MPI_Init(mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank,      mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, num_procs, mpierr)

  call get_command_argument(1, protein_file)
  call get_command_argument(2, simulation_file)

  run = 1
  t = 0.0
  bin = .true.

  call get_protein_params(protein_file)
  call get_simulation_params(simulation_file)

  write(latt_file, '(a, a, a, i0, a)') trim(adjustl(outdir)),&
    trim(adjustl(lattice_name)), "_lattice_", n_sites, ".txt"
  write(hist_file, '(a, a, a, i0, a)') trim(adjustl(outdir)),&
    trim(adjustl(protein_name)), "_run_", run, ".csv"

  call generate_lattice(lattice_name, n_sites)
  call generate_histogram()

  call construct_pulse(fwhm, dt1, fluence)
  do while (t.lt.tmax)
    call mc_step(t, bin)
    t = t + dt1
  end do

  if (rank.eq.1) then
    call print_lattice(latt_file, coords, neighbours)
    call write_histogram(hist_file)
  end if

  call MPI_Finalize(mpierr)

end program main
