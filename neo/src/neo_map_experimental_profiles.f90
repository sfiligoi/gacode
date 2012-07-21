!-----------------------------------------------------------
! neo_map_experimental_profiles.f90 
!
! PURPOSE:
!  Generate radial profile functions on the simulation 
!  (slice) r-grid.
!-----------------------------------------------------------

subroutine neo_map_experimental_profiles

  use neo_globals
  use neo_profile_exp

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
  call cub_spline(r_p,rhoN_torflux_exp,n_grid_exp,r,rhoN_torflux,n_radial)
  call cub_spline(r_p,psiN_polflux_exp,n_grid_exp,r,psiN_polflux,n_radial)
  call cub_spline(r_p,dphi0dr_p,n_grid_exp,r,dphi0dr,n_radial)
  call cub_spline(r_p,omega_rot_p,n_grid_exp,r,omega_rot,n_radial)
  call cub_spline(r_p,omega_rot_deriv_p,n_grid_exp,r,omega_rot_deriv,n_radial)
  call cub_spline(r_p,q_exp,n_grid_exp,r,q,n_radial)
  call cub_spline(r_p,shat_p,n_grid_exp,r,shat,n_radial)
  call cub_spline(r_p,shift_p,n_grid_exp,r,shift,n_radial)
  call cub_spline(r_p,kappa_exp,n_grid_exp,r,kappa,n_radial)
  call cub_spline(r_p,s_kappa_p,n_grid_exp,r,s_kappa,n_radial)
  call cub_spline(r_p,delta_exp,n_grid_exp,r,delta,n_radial)
  call cub_spline(r_p,s_delta_p,n_grid_exp,r,s_delta,n_radial)
  call cub_spline(r_p,zeta_exp,n_grid_exp,r,zeta,n_radial)
  call cub_spline(r_p,s_zeta_p,n_grid_exp,r,s_zeta,n_radial)
  call cub_spline(r_p,zmag_p,n_grid_exp,r,zmag,n_radial)
  call cub_spline(r_p,s_zmag_p,n_grid_exp,r,s_zmag,n_radial)

  if(geo_numeq_flag == 1) then
     do i=1,8
        do j=0,geo_ny
           call cub_spline(r_p,geo_yin_exp(i,j,:),n_grid_exp,r, &
                geo_yin(i,j,:),n_radial)
        enddo
     enddo
  else
     geo_yin(:,:,:) = 0.0
  endif

  call cub_spline(r_p,te_ade_exp,n_grid_exp,r,te_ade,n_radial)
  call cub_spline(r_p,ne_ade_exp,n_grid_exp,r,ne_ade,n_radial)

  do i=1,n_species
     ! Note: maping is only done for n_species (not n_species_exp)
     call cub_spline(r_p,den_exp(i,:),n_grid_exp,r,dens(i,:),n_radial)
     call cub_spline(r_p,tem_exp(i,:),n_grid_exp,r,temp(i,:),n_radial)
     call cub_spline(r_p,dlntdr_p(i,:),n_grid_exp,r,dlntdr(i,:),n_radial)
     call cub_spline(r_p,dlnndr_p(i,:),n_grid_exp,r,dlnndr(i,:),n_radial)
  enddo

  call cub_spline(r_p,b_unit_p,n_grid_exp,r,b_unit,n_radial)
  call cub_spline(r_p,rmaj_p,n_grid_exp,r,rmaj,n_radial)
  !------------------------------------------------------------------

  ! EAB diagnostic
  if(prof_check_flag) then
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=30,file=trim(path)//'out.neo.profcheckp',status='replace')
        do ir=1,n_grid_exp
           write (30,'(e16.8)',advance='no') r_p(ir)
           write (30,'(e16.8)',advance='no') b_unit_p(ir)
           write (30,'(e16.8)',advance='no') rmaj_p(ir)
           write (30,'(e16.8)',advance='no') q_exp(ir)
           write (30,'(e16.8)',advance='no') kappa_exp(ir)
           write (30,'(e16.8)',advance='no') s_kappa_p(ir)
           write (30,'(e16.8)',advance='no') delta_exp(ir)
           write (30,'(e16.8)',advance='no') s_delta_p(ir)
           write (30,'(e16.8)',advance='no') zeta_exp(ir)
           write (30,'(e16.8)',advance='no') s_zeta_p(ir)
           write (30,'(e16.8)',advance='no') shift_p(ir)
           write (30,'(e16.8)',advance='no') zmag_p(ir)
           write (30,'(e16.8)',advance='no') s_zmag_p(ir)
           do i=1,n_species_exp
              write (30,'(e16.8)',advance='no') den_exp(i,ir)
              write (30,'(e16.8)',advance='no') tem_exp(i,ir)
              write (30,'(e16.8)',advance='no') dlnndr_p(i,ir)
              write (30,'(e16.8)',advance='no') dlntdr_p(i,ir)
           enddo
           write (30,'(e16.8)',advance='no') dphi0dr_p(ir)
           write (30,*)
        enddo
        close(30)
        
        open(unit=30,file=trim(path)//'out.neo.profchecks',status='replace')
        do ir=1,n_radial
           write (30,'(e16.8)',advance='no') r(ir)
           write (30,'(e16.8)',advance='no') b_unit(ir)
           write (30,'(e16.8)',advance='no') rmaj(ir)
           write (30,'(e16.8)',advance='no') q(ir)
           write (30,'(e16.8)',advance='no') kappa(ir)
           write (30,'(e16.8)',advance='no') s_kappa(ir)
           write (30,'(e16.8)',advance='no') delta(ir)
           write (30,'(e16.8)',advance='no') s_delta(ir)
           write (30,'(e16.8)',advance='no') zeta(ir)
           write (30,'(e16.8)',advance='no') s_zeta(ir)
           write (30,'(e16.8)',advance='no') shift(ir)
           write (30,'(e16.8)',advance='no') zmag(ir)
           write (30,'(e16.8)',advance='no') s_zmag(ir)
           do i=1,n_species
              write (30,'(e16.8)',advance='no') dens(i,ir)
              write (30,'(e16.8)',advance='no') temp(i,ir)
              write (30,'(e16.8)',advance='no') dlnndr(i,ir)
              write (30,'(e16.8)',advance='no') dlntdr(i,ir)
           enddo
           write (30,'(e16.8)',advance='no') dphi0dr(ir)
           write (30,*)
        enddo
        close(30)
        
        open(unit=30,file=trim(path)//'out.neo.profchecks_rhoN',status='replace')
        do ir=1,n_radial
           write (30,'(e16.8)',advance='no') r(ir)
           write (30,'(e16.8)',advance='no') rhoN_torflux(ir)
           write (30,'(e16.8)',advance='no') psiN_polflux(ir)
           write (30,*)
        enddo
        close(30)
        
     end if

  endif

end subroutine neo_map_experimental_profiles
