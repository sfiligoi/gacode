subroutine tgyro_profile_reintegrate

  use tgyro_globals
  use tgyro_ped

  implicit none

  integer :: i_ion
  integer :: i_star  

  ! Map data inside r < r(n_r)

  if (tgyro_ped_model > 1) then
     ! Map data past r(n_r)
     call tgyro_pedestal_map(dlnnedr(n_r),zn_top,n_top(1),nn_vec(:,2),i_star,exp_ne)
     exp_ne(i_star:n_exp) = exp_ne(i_star:n_exp)*1e-13
     call tgyro_pedestal_map(dlntedr(n_r),zt_top,t_top(1),t_vec(:),i_star,exp_te)
     exp_te(i_star:n_exp) = exp_te(i_star:n_exp)*1e-3
     call tgyro_pedestal_map(dlntidr(1,n_r),zt_top,t_top(1),t_vec(:),i_star,exp_ti(1,:))
     exp_ti(1,i_star:n_exp) = exp_ti(1,i_star:n_exp)*1e-3

     ! Set ion densities
     exp_ni(1,i_star:n_exp) = exp_ne(i_star:n_exp)
     if (loc_n_ion > 1) then
        exp_ni(1,i_star:n_exp) = exp_ni(1,i_star:n_exp)-zi_vec(2)*exp_ni(2,i_star:n_exp)
     endif
     if (loc_n_ion > 2) then
        exp_ni(1,i_star:n_exp) = exp_ni(1,i_star:n_exp)-zi_vec(3)*exp_ni(3,i_star:n_exp)
     endif
     ! Set thermal ion temperatures
     do i_ion=2,loc_n_ion
        if (therm_flag(i_ion) == 1) exp_ti(i_ion,i_star:n_exp) = exp_ti(1,i_star:n_exp)
     enddo
  endif

  call tgyro_expro_map(r,dlnnedr,n_r,rmin_exp,exp_ne,n_exp)
  call tgyro_expro_map(r,dlntedr,n_r,rmin_exp,exp_te,n_exp)
  do i_ion=1,loc_n_ion
     if (therm_flag(i_ion) == 1) then
        call tgyro_expro_map(r,dlnnidr(i_ion,:),n_r,rmin_exp,exp_ni(i_ion,:),n_exp)
        call tgyro_expro_map(r,dlntidr(i_ion,:),n_r,rmin_exp,exp_ti(i_ion,:),n_exp)
     endif
  enddo
  ptot_exp = exp_ne*exp_te
  do i_ion=1,loc_n_ion
     ptot_exp = ptot_exp + exp_ni(i_ion,:)*exp_ti(i_ion,:)
  enddo
  ! Convert to Pa (10^13 n)*(kT * 1e3)/10  
  ptot_exp = 1e15*ptot_exp*k 

end subroutine tgyro_profile_reintegrate
