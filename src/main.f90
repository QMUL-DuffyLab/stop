program main
  use mpi_f08
  use input
  use lattice
  implicit none
  character(100) :: protein_file, simulation_file, lattice_output_file
  integer :: num_procs, rank, mpierr, i

  call MPI_Init(mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank,      mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, num_procs, mpierr)

  rank = 0
  call get_command_argument(1, protein_file)
  call get_command_argument(2, simulation_file)

  call get_protein_params(protein_file)
  write(*, *) rank, n_s, p_names
  call get_simulation_params(simulation_file)
  write(lattice_output_file, '(a, a, a, i0, a)') trim(adjustl(outdir)),&
    trim(adjustl(lattice_name)), "_lattice_", n_sites, ".txt"
  write(*, *) rank, lattice_name, fwhm, dt1, n_counts, lattice_output_file
  call generate_lattice(lattice_name, n_sites)
  if (rank.eq.0) then
    do i = 1, n_sites
      write(*, *) i, coords(i, :), neighbours(i, :)
    end do
    write(*, '(a)') lattice_output_file
  end if

  call MPI_Finalize(mpierr)

end program main
