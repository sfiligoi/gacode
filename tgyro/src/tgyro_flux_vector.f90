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

  integer :: i
  integer :: p
  integer, intent(in) :: index

  real, intent(in), dimension(p_max) :: x_vec
  real, intent(inout), dimension(p_max) :: f_vec
  real, intent(in) :: dx

  call tgyro_profile_set(x_vec,dx,index)

  call tgyro_flux
  call tgyro_comm_sync

  p = 0
  do i=2,n_r

     if (loc_ti_feedback_flag == 1) then
        p = p+1
        f_vec(p) = eflux_i_tot(i)
     endif

     if (loc_te_feedback_flag == 1) then
        p = p+1
        f_vec(p) = eflux_e_tot(i)
     endif

     if (loc_ne_feedback_flag == 1) then
        p = p+1
        f_vec(p) = pflux_e_tot(i)
     endif

     if (loc_er_feedback_flag == 1) then
        p = p+1
        f_vec(p) = mflux_tot(i)
     endif

     if (loc_he_feedback_flag == 1) then
        p = p+1
        f_vec(p) = pflux_he_tot(i)
     endif

  enddo

end subroutine tgyro_flux_vector
