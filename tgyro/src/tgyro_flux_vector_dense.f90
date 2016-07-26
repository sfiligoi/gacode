subroutine tgyro_flux_vector_dense(x_vec,f_vec)

  use tgyro_globals

  implicit none

  integer :: i
  integer :: p

  real, intent(in), dimension(p_max) :: x_vec
  real, intent(inout), dimension(p_max) :: f_vec

  call tgyro_profile_set(x_vec,0.0,0)

  ! Make profile consistent with new gradients
  call tgyro_profile_functions
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

end subroutine tgyro_flux_vector_dense
