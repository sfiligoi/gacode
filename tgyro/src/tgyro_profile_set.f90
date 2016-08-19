!-----------------------------------------------------------------------
! tgyro_profile_set.f90
!
! PURPOSE:
!  Manage common operation of mapping some modified profiles to internal 
!  profiles.
!-----------------------------------------------------------------------

subroutine tgyro_profile_set(x_vec,dx,index)

  use tgyro_globals

  implicit none

  integer :: i
  integer :: p
  integer, intent(in) :: index

  real, intent(in), dimension(p_max) :: x_vec
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
        ! Set dlnnidr for all ions at given radius according to quasineutrality
        call tgyro_quasigrad(ne(i),dlnnedr(i),ni(:,i),dlnnidr(:,i),zi_vec(1:loc_n_ion),loc_n_ion)
     endif

     if (loc_er_feedback_flag == 1) then
        p = p+1
        if (index == 4) then
           f_rot(i) = x_vec(p)+dx
        else
           f_rot(i) = x_vec(p)
        endif
     endif

     if (loc_he_feedback_flag == 1) then
        p = p+1
        if (index == 5) then
           dlnnidr(i_ash,i) = x_vec(p)+dx
        else
           dlnnidr(i_ash,i) = x_vec(p)
        endif
     endif

  enddo

end subroutine tgyro_profile_set
