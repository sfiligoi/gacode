subroutine tgyro_profile_reintegrate

  use tgyro_globals
  use tgyro_ped
  use tgyro_iteration_variables

  implicit none

  integer :: i_ion
  integer :: i_star  

  if (tgyro_ped_model > 1) then

     ! Map data over r(n_r) < r < a

     call tgyro_pedestal_map(dlnnedr(n_r),zn_top,n_top(1),nn_vec(:,2),i_star,exp_ne)
     call tgyro_pedestal_map(dlntedr(n_r),zt_top,t_top(1),t_vec(:),i_star,exp_te)
     call tgyro_pedestal_map(dlntidr(1,n_r),zt_top,t_top(1),t_vec(:),i_star,exp_ti(1,:))

     ! Map ion densities
     ! NOTE: Assumption is that ion profiles are 
     !
     !    ni(j,r) = ne(r)*ratio(j) for r > r(n_r)
     !
     ! where ratio(j) is the pivot density ratio at t=0.
     !
     do i_ion=1,loc_n_ion
        exp_ni(i_ion,i_star:n_exp) = exp_ne(i_star:n_exp)*n_ratio(i_ion)
     enddo
     
     ! Set thermal ion temperatures
     do i_ion=2,loc_n_ion
        if (therm_flag(i_ion) == 1) exp_ti(i_ion,i_star:n_exp) = exp_ti(1,i_star:n_exp)
     enddo
  endif

  ! Map data inside r < r(n_r)

  call tgyro_expro_map(r,dlnnedr,n_r,ne(n_r),rmin_exp,exp_ne,n_exp,'log')
  call tgyro_expro_map(r,dlntedr,n_r,te(n_r),rmin_exp,exp_te,n_exp,'log')
  do i_ion=1,loc_n_ion
     if (therm_flag(i_ion) == 1) then
        call tgyro_expro_map(r,dlnnidr(i_ion,:),n_r,ni(i_ion,n_r),rmin_exp,exp_ni(i_ion,:),n_exp,'log')
        call tgyro_expro_map(r,dlntidr(i_ion,:),n_r,ti(i_ion,n_r),rmin_exp,exp_ti(i_ion,:),n_exp,'log')
     endif
  enddo
  call tgyro_expro_map(r,w0p,n_r,w0(n_r),rmin_exp,exp_w0,n_exp,'lin')

  ptot_exp = exp_ne*exp_te
  do i_ion=1,loc_n_ion
     ptot_exp = ptot_exp + exp_ni(i_ion,:)*exp_ti(i_ion,:)
  enddo

  ! Convert to Pa: n[1/cm^3]*(kT[ev])/10  
  ptot_exp = ptot_exp*k/10.0

end subroutine tgyro_profile_reintegrate
