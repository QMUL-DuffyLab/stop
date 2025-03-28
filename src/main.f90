program main
  use mpi_f08
  use io
  use lattice
  implicit none
  character(100) :: protein_file, simulation_file, lattice_output_file
  integer :: num_procs, rank, mpierr

  call MPI_Init(mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank,      mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, num_procs, mpierr)

  call get_command_argument(1, protein_file)
  call get_command_argument(2, simulation_file)

  call get_protein_params(protein_file)
  call get_simulation_params(simulation_file)
  write(lattice_output_file, '(a, a, a, i0, a)') trim(adjustl(outdir)),&
    trim(adjustl(lattice_name)), "_lattice_", n_sites, ".txt"
  call generate_lattice(lattice_name, n_sites)
  if (rank.eq.1) then
    call print_lattice(lattice_output_file, coords, neighbours)
  end if

  call MPI_Finalize(mpierr)

end program main
