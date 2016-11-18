!-----------------------------------------------------------------------
! tgyro_flux_vector.f90
!
! PURPOSE:
!  Manage perturbations of experimental data to compute Jacobian.
!  ONLY this routine (or tgyro_flux_vector_dense) calls tgyro_flux.
!-----------------------------------------------------------------------

subroutine tgyro_flux_vector(x_vec,f_vec,dx,index)

  use tgyro_globals

  implicit none

  integer, intent(in) :: index

  real, intent(in), dimension(p_max) :: x_vec
  real, intent(inout), dimension(p_max) :: f_vec
  real, intent(in) :: dx

  call tgyro_profile_set(x_vec,dx,index)
  call tgyro_flux
  call tgyro_comm_sync
  call tgyro_flux_set(f_vec)

end subroutine tgyro_flux_vector
