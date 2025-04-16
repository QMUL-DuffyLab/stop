program main
  use mpi_f08
  use io
  use lattice
  use mc
  implicit none
  character(100) :: protein_file, simulation_file, latt_file, hist_file, path
  integer(kind=CI) :: num_procs, rank, mpierr, i, salt
  real(kind=CF) :: t_start, t_end

  call MPI_Init(mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank,      mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, num_procs, mpierr)

  call get_command_argument(1, protein_file)
  call get_command_argument(2, simulation_file)

  call get_protein_params(protein_file)
  call get_simulation_params(simulation_file)

  call generate_lattice(lattice_name, n_sites)
  call generate_histogram()
  call construct_pulse(fwhm, dt1, fluence)

  if (rank.eq.0) then
    call cpu_time(t_start)
  end if

  do i = 1, n_repeats

    write(*, *) "Run ", i

    salt = rank + num_procs * i

    write(path, '(a, a, a, i0, a, i0, a)') trim(adjustl(outdir)),&
      trim(adjustl(protein_name)), "_run_", i, "_proc_", rank

    call do_run(salt, int(n_counts / num_procs), path)

    call MPI_Barrier(MPI_COMM_WORLD, mpierr)

    if (rank.eq.0) then
      call MPI_Reduce(MPI_IN_PLACE, counts, size(counts), MPI_INT,&
        MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    else
      call MPI_Reduce(counts, counts, size(counts), MPI_INT,&
        MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    end if

    write(latt_file, '(a, a, a, i0, a)') trim(adjustl(outdir)),&
      trim(adjustl(lattice_name)), "_lattice_", n_sites, ".txt"
    write(hist_file, '(a, a, a, i0, a)') trim(adjustl(outdir)),&
      trim(adjustl(protein_name)), "_run_", i, ".csv"

    if (rank.eq.0) then
      call print_lattice(latt_file, coords, neighbours)
      call write_histogram(hist_file)
    end if

  end do

  if (rank.eq.0) then
    call cpu_time(t_end)
    write(*, '(a, F8.3, a)') "Time elapsed: ", t_end - t_start, "s"
  end if

  call MPI_Finalize(mpierr)

end program main
