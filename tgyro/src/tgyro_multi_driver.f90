!-----------------------------------------------------------
! tgyro_multi_driver.f90
!
! PURPOSE:
!  Main driver for multi-job utility.
!----------------------------------------------------------

subroutine tgyro_multi_driver

  use mpi
  use tgyro_globals
  use gyro_interface

  implicit none


  if (lcode == 'gyro') then

     ! See gyro/src/gyro_globals.f90 for definition of transport_method
     transport_method = 1

     ! Initialize GYRO
     call gyro_init(lpath,gyro_comm)

     ! Run GYRO
     call gyro_run(gyrotest_flag,gyro_restart_method,transport_method)

     ! These error variables part of gyro_interface
     call tgyro_trap_component_error(gyro_error_status_out,gyro_error_message_out)

  else 

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
     
  endif


end subroutine tgyro_multi_driver

