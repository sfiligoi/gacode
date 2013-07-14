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

  p = 0
  do i=2,n_r

     if (loc_ti_feedback_flag == 1) then
        p = p+1
        if (index == 1) then
           dlntidr(therm_vec(:),i) = x_vec(p)+dx
        else
           dlntidr(therm_vec(:),i) = x_vec(p)
        endif
     endif

     if (loc_te_feedback_flag == 1) then
        p = p+1
        if (index == 2) then
           dlntedr(i) = x_vec(p)+dx
        else
           dlntedr(i) = x_vec(p)
        endif
     endif

     if (loc_ne_feedback_flag == 1) then 
        p = p+1
        if (index == 3) then
           dlnnedr(i) = x_vec(p)+dx
        else
           dlnnedr(i) = x_vec(p)
        endif
        ! Set dlnnidr(1,i) according to quasineutrality
        call tgyro_quasigrad(ne(i),dlnnedr(i),ni(:,i),dlnnidr(:,i),zi_vec(:),loc_n_ion,dlnridr(:,i))
     endif

     if (loc_er_feedback_flag == 1) then
        p = p+1
        if (index == 4) then
           f_rot(i) = x_vec(p)+dx
        else
           f_rot(i) = x_vec(p)
        endif
     endif

  enddo

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

  enddo

end subroutine tgyro_flux_vector
