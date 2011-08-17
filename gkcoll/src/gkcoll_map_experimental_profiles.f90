!-----------------------------------------------------------
! gkcoll_map_experimental_profiles.f90 
!
! PURPOSE:
!  Generate radial profile functions on the simulation 
!  (slice) r-grid.
!-----------------------------------------------------------

subroutine gkcoll_map_experimental_profiles

  use gkcoll_globals
  use gkcoll_profile_exp

  !-------------------------------------------------
  implicit none
  !
  integer :: i, j, ir
  logical :: prof_check_flag = .true.
  !-------------------------------------------------

  !------------------------------------------------------------------
  ! Use local cubic spline interpolation to get simulation 
  ! profiles from experimental (_p) ones.
  ! 
  call cub_spline(r_p,rhoN_torflux_exp,n_grid_exp,r,rhoN_torflux,n_gr)
  call cub_spline(r_p,psiN_polflux_exp,n_grid_exp,r,psiN_polflux,n_gr)
  call cub_spline(r_p,q_exp,n_grid_exp,r,q,n_gr)
  call cub_spline(r_p,shat_p,n_grid_exp,r,shat,n_gr)
  call cub_spline(r_p,shift_p,n_grid_exp,r,shift,n_gr)
  call cub_spline(r_p,kappa_exp,n_grid_exp,r,kappa,n_gr)
  call cub_spline(r_p,s_kappa_p,n_grid_exp,r,s_kappa,n_gr)
  call cub_spline(r_p,delta_exp,n_grid_exp,r,delta,n_gr)
  call cub_spline(r_p,s_delta_p,n_grid_exp,r,s_delta,n_gr)
  call cub_spline(r_p,zeta_exp,n_grid_exp,r,zeta,n_gr)
  call cub_spline(r_p,s_zeta_p,n_grid_exp,r,s_zeta,n_gr)
  call cub_spline(r_p,zmag_p,n_grid_exp,r,zmag,n_gr)
  call cub_spline(r_p,s_zmag_p,n_grid_exp,r,s_zmag,n_gr)

  if(geo_numeq_flag == 1) then
     do i=1,8
        do j=0,geo_ny
           call cub_spline(r_p,geo_yin_exp(i,j,:),n_grid_exp,r, &
                geo_yin(i,j,:),n_gr)
        enddo
     enddo
  else
     geo_yin(:,:,:) = 0.0
  endif

  call cub_spline(r_p,te_ade_exp,n_grid_exp,r,te_ade,n_gr)
  call cub_spline(r_p,ne_ade_exp,n_grid_exp,r,ne_ade,n_gr)

  do i=1,n_species
     ! Note: maping is only done for n_species (not n_species_exp)
     call cub_spline(r_p,den_exp(i,:),n_grid_exp,r,dens(i,:),n_gr)
     call cub_spline(r_p,tem_exp(i,:),n_grid_exp,r,temp(i,:),n_gr)
     call cub_spline(r_p,dlntdr_p(i,:),n_grid_exp,r,dlntdr(i,:),n_gr)
     call cub_spline(r_p,dlnndr_p(i,:),n_grid_exp,r,dlnndr(i,:),n_gr)
  enddo

  call cub_spline(r_p,b_unit_p,n_grid_exp,r,b_unit,n_gr)
  call cub_spline(r_p,rmaj_p,n_grid_exp,r,rmaj,n_gr)
  !------------------------------------------------------------------

end subroutine gkcoll_map_experimental_profiles
