!-----------------------------------------------------------
! tgyro_multi_driver.f90
!
! PURPOSE:
!  Main driver for multi-job utility.
!----------------------------------------------------------

subroutine tgyro_multi_driver

  use mpi
  use tgyro_globals

  implicit none

  allocate(cgyro_status_vec(n_inst))
  allocate(cgyro_n_species_vec(n_inst))
  allocate(cgyro_tave_min_vec(n_inst))
  allocate(cgyro_tave_max_vec(n_inst))
  allocate(cgyro_flux_tave_vec(11,3,n_inst))
  allocate(cgyro_flux_tave_out(11,3))

  ! Initialize CGYRO
  call cgyro_init(lpath,gyro_comm)
  call cgyro_init_kernel

  ! Run CGYRO
  cgyro_var_in = -1.0
  call cgyro_run(gyrotest_flag,cgyro_var_in,cgyro_n_species_out, &
       cgyro_flux_tave_out,cgyro_tave_min_out,cgyro_tave_max_out,&
       cgyro_status_out)
  call cgyro_final_kernel

  if(allocated(cgyro_status_vec))    deallocate(cgyro_status_vec)
  if(allocated(cgyro_n_species_vec)) deallocate(cgyro_n_species_vec)
  if(allocated(cgyro_tave_min_vec))  deallocate(cgyro_tave_min_vec)
  if(allocated(cgyro_tave_max_vec))  deallocate(cgyro_tave_max_vec)
  if(allocated(cgyro_flux_tave_vec)) deallocate(cgyro_flux_tave_vec)
  if(allocated(cgyro_flux_tave_out)) deallocate(cgyro_flux_tave_out)

end subroutine tgyro_multi_driver

