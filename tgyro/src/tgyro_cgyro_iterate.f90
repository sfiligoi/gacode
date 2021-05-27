! Multi-job cgyro utility with feedback iteration

subroutine tgyro_cgyro_iterate
  use mpi
  use tgyro_globals
  implicit none
  integer :: i, j, p, is
  
  allocate(cgyro_n_species_vec(n_inst))
  allocate(cgyro_tave_min_vec(n_inst))
  allocate(cgyro_tave_max_vec(n_inst))
  allocate(cgyro_flux_tave_vec(11,3,n_inst))
  allocate(cgyro_flux_tave_out(11,3))
  cgyro_nflux=33
  
  if (i_proc_global == 0) then
     open(unit=1,file=trim(runfile_cgyro_eflux),status='replace')
     write(1,*) '# dir tmin tmax ns Gamma(ns) Q(ns) Pi(ns)'
     write(1,*)
     close(1)
  endif

  ! Initialize CGYRO
  call cgyro_init(lpath,gyro_comm)
  
  do i=0,tgyro_cgyro_n_iterate
  
     ! Run multiple cgyro's
     cgyro_var_in = 2.0        ! max time
     call cgyro_run(gyrotest_flag,cgyro_var_in,cgyro_n_species_out, &
          cgyro_flux_tave_out,cgyro_tave_min_out,cgyro_tave_max_out)

     ! Collect flux data
     call MPI_GATHER(cgyro_n_species_out,1,MPI_INTEGER,&
          cgyro_n_species_vec,1,MPI_INTEGER,0,gyro_adj,ierr)
     call MPI_GATHER(cgyro_tave_min_out,1,MPI_DOUBLE_PRECISION,&
          cgyro_tave_min_vec,1,MPI_DOUBLE_PRECISION,0,gyro_adj,ierr)
     call MPI_GATHER(cgyro_tave_max_out,1,MPI_DOUBLE_PRECISION,&
          cgyro_tave_max_vec,1,MPI_DOUBLE_PRECISION,0,gyro_adj,ierr)
     call MPI_GATHER(cgyro_flux_tave_out,cgyro_nflux,MPI_DOUBLE_PRECISION,&
          cgyro_flux_tave_vec,cgyro_nflux,MPI_DOUBLE_PRECISION,0,gyro_adj,ierr)
     
     if (i_proc_global == 0) then
        open(unit=1,file=trim(runfile_cgyro_eflux),position='append')
        do j=1,n_inst
           write(1,'(a)',advance='no')        trim(paths(j))
           write(1,'(f7.3)',advance='no')     cgyro_tave_min_vec(j)
           write(1,'(f7.3)',advance='no')     cgyro_tave_max_vec(j)
           write(1,'(i2)',  advance='no')     cgyro_n_species_vec(j)
           do p=1,3
              do is=1,cgyro_n_species_vec(j)
                 write(1,'(e14.5)',advance='no') cgyro_flux_tave_vec(is,p,j)
              enddo
           enddo
           write(1,*)
        enddo
        close(1)
     endif

  enddo
     
  if(allocated(cgyro_n_species_vec)) deallocate(cgyro_n_species_vec)
  if(allocated(cgyro_tave_min_vec))  deallocate(cgyro_tave_min_vec)
  if(allocated(cgyro_tave_max_vec))  deallocate(cgyro_tave_max_vec)
  if(allocated(cgyro_flux_tave_vec)) deallocate(cgyro_flux_tave_vec)
  if(allocated(cgyro_flux_tave_out)) deallocate(cgyro_flux_tave_out)

  
end subroutine tgyro_cgyro_iterate
