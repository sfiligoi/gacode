!-----------------------------------------------------------------
! tgyro_comm_sync.f90
!
! PURPOSE:
!  Synchronization (gather, broadcast) of variables returned by 
!  calls to flux routines.
!-----------------------------------------------------------------

subroutine tgyro_comm_sync

  use mpi
  use tgyro_globals

  implicit none

  real, dimension(n_r-1) :: collect
  integer :: i_ion

  ! Only a single flux and root (gradient) was computed 
  ! on a given processor.  Do an allgather to make these 
  ! quantities global.  A subsequent broadcast is in fact
  ! required when the number of processors in all 
  ! instances of gyro_comm is not equal.

  ! sync eflux_i

  call MPI_ALLGATHER(eflux_i_tot(i_r),&
       1,&
       MPI_DOUBLE_PRECISION,&
       collect,&
       1,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  eflux_i_tot(2:n_r) = collect(:)

  call MPI_BCAST(eflux_i_tot(2:n_r),&
       n_r-1,&
       MPI_DOUBLE_PRECISION,&
       0,&
       gyro_comm,&
       ierr)

  ! sync eflux_e

  call MPI_ALLGATHER(eflux_e_tot(i_r),&
       1,&
       MPI_DOUBLE_PRECISION,&
       collect,&
       1,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  eflux_e_tot(2:n_r) = collect(:)

  call MPI_BCAST(eflux_e_tot(2:n_r),&
       n_r-1,&
       MPI_DOUBLE_PRECISION,&
       0,&
       gyro_comm,&
       ierr)

  ! sync pflux_i

  call MPI_ALLGATHER(pflux_i_tot(i_r),&
       1,&
       MPI_DOUBLE_PRECISION,&
       collect,&
       1,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  pflux_i_tot(2:n_r) = collect(:)

  call MPI_BCAST(pflux_i_tot(2:n_r),&
       n_r-1,&
       MPI_DOUBLE_PRECISION,&
       0,&
       gyro_comm,&
       ierr)

  ! sync pflux_e

  call MPI_ALLGATHER(pflux_e_tot(i_r),&
       1,&
       MPI_DOUBLE_PRECISION,&
       collect,&
       1,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  pflux_e_tot(2:n_r) = collect(:)

  call MPI_BCAST(pflux_e_tot(2:n_r),&
       n_r-1,&
       MPI_DOUBLE_PRECISION,&
       0,&
       gyro_comm,&
       ierr)

  ! sync mflux

  call MPI_ALLGATHER(mflux_tot(i_r),&
       1,&
       MPI_DOUBLE_PRECISION,&
       collect,&
       1,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  mflux_tot(2:n_r) = collect(:)

  call MPI_BCAST(mflux_tot(2:n_r),&
       n_r-1,&
       MPI_DOUBLE_PRECISION,&
       0,&
       gyro_comm,&
       ierr)

  ! sync-0 GYRO diagnostics

  do i_ion=1,loc_n_ion

     call MPI_ALLGATHER(eflux_i_neo(i_ion,i_r),&
          1,&
          MPI_DOUBLE_PRECISION,&
          collect,&
          1,&
          MPI_DOUBLE_PRECISION,&
          gyro_adj,&
          ierr)

     eflux_i_neo(i_ion,2:n_r) = collect(:)

     call MPI_ALLGATHER(eflux_i_tur(i_ion,i_r),&
          1,&
          MPI_DOUBLE_PRECISION,&
          collect,&
          1,&
          MPI_DOUBLE_PRECISION,&
          gyro_adj,&
          ierr)

     eflux_i_tur(i_ion,2:n_r) = collect(:)

     call MPI_ALLGATHER(pflux_i_neo(i_ion,i_r),&
          1,&
          MPI_DOUBLE_PRECISION,&
          collect,&
          1,&
          MPI_DOUBLE_PRECISION,&
          gyro_adj,&
          ierr)

     pflux_i_neo(i_ion,2:n_r) = collect(:)

     call MPI_ALLGATHER(pflux_i_tur(i_ion,i_r),&
          1,&
          MPI_DOUBLE_PRECISION,&
          collect,&
          1,&
          MPI_DOUBLE_PRECISION,&
          gyro_adj,&
          ierr)

     pflux_i_tur(i_ion,2:n_r) = collect(:)

     call MPI_ALLGATHER(mflux_i_neo(i_ion,i_r),&
          1,&
          MPI_DOUBLE_PRECISION,&
          collect,&
          1,&
          MPI_DOUBLE_PRECISION,&
          gyro_adj,&
          ierr)

     mflux_i_neo(i_ion,2:n_r) = collect(:)

     call MPI_ALLGATHER(mflux_i_tur(i_ion,i_r),&
          1,&
          MPI_DOUBLE_PRECISION,&
          collect,&
          1,&
          MPI_DOUBLE_PRECISION,&
          gyro_adj,&
          ierr)

     mflux_i_tur(i_ion,2:n_r) = collect(:)

     call MPI_ALLGATHER(expwd_i_tur(i_ion,i_r),&
          1,&
          MPI_DOUBLE_PRECISION,&
          collect,&
          1,&
          MPI_DOUBLE_PRECISION,&
          gyro_adj,&
          ierr)

     expwd_i_tur(i_ion,2:n_r) = collect(:)

  enddo

  call MPI_ALLGATHER(eflux_e_neo(i_r),&
       1,&
       MPI_DOUBLE_PRECISION,&
       collect,&
       1,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  eflux_e_neo(2:n_r) = collect(:)

  call MPI_ALLGATHER(eflux_e_tur(i_r),&
       1,&
       MPI_DOUBLE_PRECISION,&
       collect,&
       1,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  eflux_e_tur(2:n_r) = collect(:)

  call MPI_ALLGATHER(pflux_e_neo(i_r),&
       1,&
       MPI_DOUBLE_PRECISION,&
       collect,&
       1,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  pflux_e_neo(2:n_r) = collect(:)

  call MPI_ALLGATHER(pflux_e_tur(i_r),&
       1,&
       MPI_DOUBLE_PRECISION,&
       collect,&
       1,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  pflux_e_tur(2:n_r) = collect(:)

  call MPI_ALLGATHER(mflux_e_neo(i_r),&
       1,&
       MPI_DOUBLE_PRECISION,&
       collect,&
       1,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  mflux_e_neo(2:n_r) = collect(:)

  call MPI_ALLGATHER(mflux_e_tur(i_r),&
       1,&
       MPI_DOUBLE_PRECISION,&
       collect,&
       1,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  mflux_e_tur(2:n_r) = collect(:)

  call MPI_ALLGATHER(expwd_e_tur(i_r),&
       1,&
       MPI_DOUBLE_PRECISION,&
       collect,&
       1,&
       MPI_DOUBLE_PRECISION,&
       gyro_adj,&
       ierr)

  expwd_e_tur(2:n_r) = collect(:)

end subroutine tgyro_comm_sync
