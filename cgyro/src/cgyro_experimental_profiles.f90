!----------------------------------------------------------------
! cgyro_experimental_profiles.f90
!
! PURPOSE:
!  Read experimental profiles. and generate additional profiles.
!
! NOTES: 
!  See http://fusion.gat.com/theory/INPUT_profiles for a 
!  complete description of the INPUT_profiles file structure.
! 
!  The files read by this routine are *generated* from 
!  INPUT_profiles by scripts. 
!----------------------------------------------------------------

subroutine cgyro_experimental_profiles(path,comm,numeq_flag,udsymmetry_flag,n_species,z,btccw,ipccw)

  use cgyro_globals, only : ae_flag,geo_yin,geo_ny,rmin, &
       dlnndr,dlntdr,sdlnndr,sdlntdr,dens,temp,&
       dlnndre_ade,dlntdre_ade,ne_ade,te_ade,is_ele
   
  use cgyro_experimental_globals
  use cgyro_io
  use EXPRO_interface

  implicit none

  character(len=*), intent(in) :: path 
  integer, intent(in) :: comm
  integer, intent(in) :: numeq_flag
  integer, intent(in) :: udsymmetry_flag
  integer, intent(in) :: n_species
  real, dimension(n_species) :: z
  real, intent(inout) :: btccw,ipccw
  
  integer :: i,j,is,i_ion
  
  n_species_exp = n_species+ae_flag
  
  !--------------------------------------------------------------
  ! use EXPRO routines to read data:
  !
  call EXPRO_palloc(comm,path,1)
  EXPRO_ctrl_quasineutral_flag = 1  ! quasi-neutrality density flag
  EXPRO_ctrl_numeq_flag = numeq_flag

  ! Number and charge of ion species
  EXPRO_ctrl_z(:) = 0.0
  EXPRO_ctrl_n_ion = n_species
  do is=1,n_species
     EXPRO_ctrl_z(is) = z(is)
  enddo

  call EXPRO_pread

  n_grid_exp = EXPRO_n_exp
  call cgyro_experimental_alloc(1)
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! Transfer data from read vector to individual arrays:
  !

  btccw = -EXPRO_signb
  ipccw = -EXPRO_signq*EXPRO_signb

  rmin_exp(:)    = EXPRO_rmin(:)
  rmaj_exp(:)    = EXPRO_rmaj(:)
  q_exp(:)       = EXPRO_q(:)
  s_exp(:)       = EXPRO_s(:)
  shift_exp(:)   = EXPRO_drmaj(:)
  kappa_exp(:)   = EXPRO_kappa(:)
  s_kappa_exp(:) = EXPRO_skappa(:)
  delta_exp(:)   = EXPRO_delta(:) 
  s_delta_exp(:) = EXPRO_sdelta(:) 
  zeta_exp(:)    = EXPRO_zeta(:)
  s_zeta_exp(:)  = EXPRO_szeta(:)
  if (udsymmetry_flag == 1) then
     zmag_exp(:)    = 0.0   
     dzmag_exp(:)   = 0.0
  else
     zmag_exp(:)    = EXPRO_zmag(:)   
     dzmag_exp(:)   = EXPRO_dzmag(:)
  endif

  ! Minor radius, a, in meters:
  a_meters = rmin_exp(n_grid_exp)

  rmin_exp(:) = rmin_exp(:)/a_meters
  rmaj_exp(:) = rmaj_exp(:)/a_meters
  zmag_exp(:) = zmag_exp(:)/a_meters

  if (numeq_flag == 1) then
     if (EXPRO_nfourier <= 0) then
        call cgyro_error('Geometry coefficients missing')
        return
     endif
     deallocate(geo_yin_exp)
     deallocate(geo_yin)
     geo_ny = EXPRO_nfourier
     allocate(geo_yin_exp(8,0:geo_ny,n_grid_exp))
     allocate(geo_yin(8,0:geo_ny))
     geo_yin_exp(1:4,:,:) = EXPRO_geo(:,:,:)/a_meters
     geo_yin_exp(5:8,:,:) = EXPRO_dgeo(:,:,:)
  endif

  te_ade_exp(:)      = EXPRO_te(:)
  ne_ade_exp(:)      = EXPRO_ne(:)
  dlntdre_ade_exp(:) = EXPRO_dlntedr(:) * a_meters
  dlnndre_ade_exp(:) = EXPRO_dlnnedr(:) * a_meters

  temp_exp(is_ele,:)   = EXPRO_te(:)
  dlntdr_exp(is_ele,:) = EXPRO_dlntedr(:) * a_meters 
  sdlntdr_exp(is_ele,:) = EXPRO_sdlntedr(:) * a_meters**2 

  dens_exp(is_ele,:)   = EXPRO_ne(:)
  dlnndr_exp(is_ele,:) = EXPRO_dlnnedr(:) * a_meters 
  sdlnndr_exp(is_ele,:) = EXPRO_sdlnnedr(:) * a_meters**2 
  
  do i_ion=1,n_species_exp-1
     ! ion temps should be equal, but not enforced 
     temp_exp(i_ion,:)   = EXPRO_ti(i_ion,:)
     dlntdr_exp(i_ion,:) = EXPRO_dlntidr(i_ion,:) * a_meters 
     sdlntdr_exp(i_ion,:) = EXPRO_sdlntidr(i_ion,:) * a_meters**2 

     ! First species density is re-set by quasi-neutrality [always 1 ?]
     if (EXPRO_ctrl_quasineutral_flag == 1 .and. i_ion == 1) then
        dens_exp(i_ion,:)   = EXPRO_ni_new(:)
        dlnndr_exp(i_ion,:) = EXPRO_dlnnidr_new(:) * a_meters
        sdlnndr_exp(i_ion,:) = EXPRO_sdlnnidr_new(:) * a_meters**2
     else
        dens_exp(i_ion,:)   = EXPRO_ni(i_ion,:)
        dlnndr_exp(i_ion,:) = EXPRO_dlnnidr(i_ion,:) * a_meters
        sdlnndr_exp(i_ion,:) = EXPRO_sdlnnidr(i_ion,:) * a_meters**2
     endif
  enddo

  ! Sanity check for densities
  do i=1,n_species_exp
     if (minval(dens_exp(i,:)) <= 0.0) then
        call cgyro_error('Nonpositive in exp. density profile')
        return
     endif
  enddo

  ! Z_eff
  z_eff_exp(:) = EXPRO_z_eff(:)

  ! Bunit 

  b_unit_exp(:) = EXPRO_bunit(:)

  ! Rotation

  gamma_e_exp(:) = -EXPRO_w0p(:) * (a_meters * rmin_exp(:)) / q_exp(:)

  gamma_p_exp(:) = -EXPRO_w0p(:) * (a_meters * rmaj_exp(:))

  mach_exp(:)    = EXPRO_w0(:) * (a_meters * rmaj_exp(:))

  rhos_exp(:)    = EXPRO_rhos(:)

  call EXPRO_palloc(comm,path,0)

  !------------------------------------------------------------------
  ! Use local cubic spline interpolation to get simulation 
  ! profiles from experimental (_exp) ones.
  ! 
  call cub_spline(rmin_exp,rmaj_exp,n_grid_exp,rmin,rmaj_loc,1)
  call cub_spline(rmin_exp,q_exp,n_grid_exp,rmin,q_loc,1)
  call cub_spline(rmin_exp,s_exp,n_grid_exp,rmin,s_loc,1)
  call cub_spline(rmin_exp,shift_exp,n_grid_exp,rmin,shift_loc,1)
  call cub_spline(rmin_exp,kappa_exp,n_grid_exp,rmin,kappa_loc,1)
  call cub_spline(rmin_exp,s_kappa_exp,n_grid_exp,rmin,s_kappa_loc,1)
  call cub_spline(rmin_exp,delta_exp,n_grid_exp,rmin,delta_loc,1)
  call cub_spline(rmin_exp,s_delta_exp,n_grid_exp,rmin,s_delta_loc,1)
  call cub_spline(rmin_exp,zeta_exp,n_grid_exp,rmin,zeta_loc,1)
  call cub_spline(rmin_exp,s_zeta_exp,n_grid_exp,rmin,s_zeta_loc,1)
  call cub_spline(rmin_exp,zmag_exp,n_grid_exp,rmin,zmag_loc,1)
  call cub_spline(rmin_exp,dzmag_exp,n_grid_exp,rmin,dzmag_loc,1)
  call cub_spline(rmin_exp,gamma_e_exp,n_grid_exp,rmin,gamma_e_loc,1)
  call cub_spline(rmin_exp,gamma_p_exp,n_grid_exp,rmin,gamma_p_loc,1)
  call cub_spline(rmin_exp,mach_exp,n_grid_exp,rmin,mach_loc,1)
  call cub_spline(rmin_exp,rhos_exp,n_grid_exp,rmin,rhos_loc,1)
  call cub_spline(rmin_exp,z_eff_exp,n_grid_exp,rmin,z_eff_loc,1)
  call cub_spline(rmin_exp,b_unit_exp,n_grid_exp,rmin,b_unit,1)

  call cub_spline(rmin_exp,te_ade_exp,n_grid_exp,rmin,te_ade,1)
  call cub_spline(rmin_exp,ne_ade_exp,n_grid_exp,rmin,ne_ade,1)
  call cub_spline(rmin_exp,dlntdre_ade_exp,n_grid_exp,rmin,dlntdre_ade,1)
  call cub_spline(rmin_exp,dlnndre_ade_exp,n_grid_exp,rmin,dlnndre_ade,1)

  do i=1,n_species
     ! Note: mapping is only done for n_species (not n_species_exp)
     call cub_spline(rmin_exp,dens_exp(i,:),n_grid_exp,rmin,dens(i),1)
     call cub_spline(rmin_exp,temp_exp(i,:),n_grid_exp,rmin,temp(i),1)
     call cub_spline(rmin_exp,dlntdr_exp(i,:),n_grid_exp,rmin,dlntdr(i),1)
     call cub_spline(rmin_exp,dlnndr_exp(i,:),n_grid_exp,rmin,dlnndr(i),1)
     call cub_spline(rmin_exp,sdlntdr_exp(i,:),n_grid_exp,rmin,sdlntdr(i),1)
     call cub_spline(rmin_exp,sdlnndr_exp(i,:),n_grid_exp,rmin,sdlnndr(i),1)
  enddo

  if (numeq_flag == 1) then
     do i=1,8
        do j=0,geo_ny
           call cub_spline(rmin_exp,geo_yin_exp(i,j,:),n_grid_exp,rmin, &
                geo_yin(i,j),1)
        enddo
     enddo
  else
     geo_yin(:,:) = 0.0
  endif

  call cgyro_experimental_alloc(0)

end subroutine cgyro_experimental_profiles
