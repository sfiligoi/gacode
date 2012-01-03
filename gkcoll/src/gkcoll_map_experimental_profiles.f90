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
  integer :: i, j
  integer, parameter :: n_gr=1
  !-------------------------------------------------

  !------------------------------------------------------------------
  ! Use local cubic spline interpolation to get simulation 
  ! profiles from experimental (_p) ones.
  ! 
  call cub_spline(r_p,rhoN_torflux_exp,n_grid_exp,rmin,rhoN_torflux,n_gr)
  call cub_spline(r_p,psiN_polflux_exp,n_grid_exp,rmin,psiN_polflux,n_gr)
  call cub_spline(r_p,q_exp,n_grid_exp,rmin,q,n_gr)
  call cub_spline(r_p,shat_p,n_grid_exp,rmin,shat,n_gr)
  call cub_spline(r_p,shift_p,n_grid_exp,rmin,shift,n_gr)
  call cub_spline(r_p,kappa_exp,n_grid_exp,rmin,kappa,n_gr)
  call cub_spline(r_p,s_kappa_p,n_grid_exp,rmin,s_kappa,n_gr)
  call cub_spline(r_p,delta_exp,n_grid_exp,rmin,delta,n_gr)
  call cub_spline(r_p,s_delta_p,n_grid_exp,rmin,s_delta,n_gr)
  call cub_spline(r_p,zeta_exp,n_grid_exp,rmin,zeta,n_gr)
  call cub_spline(r_p,s_zeta_p,n_grid_exp,rmin,s_zeta,n_gr)
  call cub_spline(r_p,zmag_p,n_grid_exp,rmin,zmag,n_gr)
  call cub_spline(r_p,s_zmag_p,n_grid_exp,rmin,s_zmag,n_gr)

  if(geo_numeq_flag == 1) then
     do i=1,8
        do j=0,geo_ny
           call cub_spline(r_p,geo_yin_exp(i,j,:),n_grid_exp,rmin, &
                geo_yin(i,j),n_gr)
        enddo
     enddo
  else
     geo_yin(:,:) = 0.0
  endif

  call cub_spline(r_p,te_ade_exp,n_grid_exp,rmin,te_ade,n_gr)
  call cub_spline(r_p,ne_ade_exp,n_grid_exp,rmin,ne_ade,n_gr)

  do i=1,n_species
     ! Note: maping is only done for n_species (not n_species_exp)
     call cub_spline(r_p,den_exp(i,:),n_grid_exp,rmin,dens(i),n_gr)
     call cub_spline(r_p,tem_exp(i,:),n_grid_exp,rmin,temp(i),n_gr)
     call cub_spline(r_p,dlntdr_p(i,:),n_grid_exp,rmin,dlntdr(i),n_gr)
     call cub_spline(r_p,dlnndr_p(i,:),n_grid_exp,rmin,dlnndr(i),n_gr)
  enddo

  call cub_spline(r_p,b_unit_p,n_grid_exp,rmin,b_norm,n_gr)
  call cub_spline(r_p,rmaj_p,n_grid_exp,rmin,rmaj,n_gr)
  !------------------------------------------------------------------

end subroutine gkcoll_map_experimental_profiles
