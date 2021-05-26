program tgyro_main

  use mpi
  use tgyro_globals
  use gyro_globals, only : path, GYRO_COMM_WORLD
  use ompdata

  !-----------------------------------------------------------------
  implicit none
  !
  integer :: supported
  integer, external :: omp_get_max_threads, omp_get_thread_num
  integer :: i, is
  !-----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Query OpenMP for dimensions
  !
  i_omp = omp_get_thread_num()
  n_omp = omp_get_max_threads()
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator, including support for 
  ! funneled threading (needed if OpenMP is enabled).
  !
  if (n_omp > 1) then
     call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,supported,ierr)
     if (supported < MPI_THREAD_FUNNELED) then
        call tgyro_catch_error('ERROR: (TGYRO) Multi-threaded MPI not supported.')
     endif
  else 
     call MPI_INIT_THREAD(MPI_THREAD_SINGLE,supported,ierr)
  endif
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc_global,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc_global,ierr)
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Get input paths and broadcast them
  !
  call tgyro_read_input   
  !
  ! At this stage, paths/procs information is only on process 0.
  !
  call tgyro_comm_setup
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Set GYRO path to local path (per instance). 
  path = lpath
  !
  ! Set GYRO communicator to local communicator (per instance).
  GYRO_COMM_WORLD = gyro_comm
  !-----------------------------------------------------------------

  select case (tgyro_mode)

  case (1)

     ! Local transport model

     call tgyro_iteration_driver 

     if (i_proc_global == 0) then
        open(unit=1,file=trim(runfile),position='append')
        write(1,*) error_msg
        close(1)
     endif

  case (3)

     ! Multi-job utility

     call tgyro_multi_driver

  case (4)

     allocate(cgyro_n_species_vec(n_inst))
     allocate(cgyro_tave_min_vec(n_inst))
     allocate(cgyro_tave_max_vec(n_inst))
     allocate(cgyro_flux_tave_vec(n_inst,11,3))
     allocate(cgyro_flux_tave_out(11,3))
     cgyro_nflux=33
     
     if (i_proc_global == 0) then
        open(unit=1,file=trim(runfile_cgyro_eflux),status='replace')
        write(1,*) '# tmin tmax Q'
        write(1,*)
        close(1)
     endif
     
     ! Multi-job cgyro utility with iteration
     call tgyro_multi_driver
     
     call MPI_GATHER(cgyro_n_species_out,1,MPI_INTEGER,&
          cgyro_n_species_vec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call MPI_GATHER(cgyro_tave_min_out,1,MPI_DOUBLE_PRECISION,&
          cgyro_tave_min_vec,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_GATHER(cgyro_tave_max_out,1,MPI_DOUBLE_PRECISION,&
          cgyro_tave_max_vec,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_GATHER(cgyro_flux_tave_out(1,2),cgyro_nflux,MPI_DOUBLE_PRECISION,&
          cgyro_flux_tave_vec(:,1,2),cgyro_nflux,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     
     if (i_proc_global == 0) then
        open(unit=1,file=trim(runfile_cgyro_eflux),position='append')
        do i=1,n_inst
           write(1,'(a)',advance='no')        trim(paths(i))
           write(1,'(f7.3)',advance='no')     cgyro_tave_min_vec(i)
           write(1,'(f7.3)',advance='no')     cgyro_tave_max_vec(i)
           write(1,'(i2)',  advance='no')     cgyro_n_species_vec(i)
           !do is=1,cgyro_n_species_vec(i)
           !   write(1,'(e11.5)',advance='no') cgyro_flux_tave_vec(is,2,i)
           !enddo
           !write(1,'(e11.5)',advance='no') cgyro_flux_tave_vec(i,1,2)
           write(1,*)
        enddo
        close(1)
     endif

     if(allocated(cgyro_n_species_vec)) deallocate(cgyro_n_species_vec)
     if(allocated(cgyro_tave_min_vec))  deallocate(cgyro_tave_min_vec)
     if(allocated(cgyro_tave_max_vec))  deallocate(cgyro_tave_max_vec)
     if(allocated(cgyro_flux_tave_vec)) deallocate(cgyro_flux_tave_vec)
     if(allocated(cgyro_flux_tave_out)) deallocate(cgyro_flux_tave_out)
     
  case default

     call tgyro_catch_error('ERROR: (TGYRO) Bad value for tgyro_mode')

  end select

  call MPI_FINALIZE(ierr)

end program tgyro_main
