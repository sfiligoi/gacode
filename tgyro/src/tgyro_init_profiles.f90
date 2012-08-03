!-----------------------------------------------------------
! tgyro_init_profiles.f90
!
! PURPOSE:
!  Manage generation of profiles on transport grid.
!----------------------------------------------------------

subroutine tgyro_init_profiles

  use mpi
  use tgyro_globals
  use EXPRO_interface

  implicit none

  integer :: i
  integer :: i_ion
  integer :: n
  integer :: n_exp
  real :: xh
  real :: arho

  !------------------------------------------------------
  ! PHYSICAL CONSTANTS
  !
  pi = 4.0*atan(1.0)
  !
  e_alpha = 3.5e6*1.6022e-12 ! eV*(erg/eV)
  e       = 4.8032e-10 ! statcoul
  k       = 1.6022e-12 ! erg/eV
  me      = 9.1094e-28 ! g
  mp      = 1.6726e-24 ! g
  c       = 2.9979e10  ! cm/s
  !------------------------------------------------------

  !------------------------------------------------------
  ! Convert dimensionless mass to grams.
  mi(:) = mi_vec(:)*mp
  !------------------------------------------------------

  !------------------------------------------------------
  ! Initialize integrated power vectors.
  !
  p_alpha(:) = 0.0
  p_exch(:)  = 0.0
  p_brem(:)  = 0.0
  !------------------------------------------------------

  !----------------------------------------------
  ! Generate radial vector:
  !
  if (tgyro_rmin > 0.0) then

     r(1) = 0.0
     do i=2,n_r

        ! Normalized r (dimensionless)
        r(i) = tgyro_rmin+(i-2)/(n_r-2.0)*(tgyro_rmax-tgyro_rmin)

     enddo

  else

     do i=1,n_r

        ! Normalized r (dimensionless)
        r(i) = (i-1)/(n_r-1.0)*tgyro_rmax

     enddo
  endif
  !----------------------------------------------

  !----------------------------------------------
  ! Radius where profiles will be matched.
  !
  i_bc = n_r-loc_bc_offset
  !----------------------------------------------

  EXPRO_ctrl_density_method = loc_quasineutral_flag+1
  EXPRO_ctrl_z = 0.0
  EXPRO_ctrl_z(1:loc_n_ion) = zi_vec(1:loc_n_ion)
  EXPRO_ctrl_numeq_flag = loc_num_equil_flag
  EXPRO_ctrl_signq = tgyro_ipccw_in*tgyro_btccw_in
  EXPRO_ctrl_signb = -tgyro_btccw_in
  EXPRO_ctrl_rotation_method = 1

  call EXPRO_palloc(MPI_COMM_WORLD,'./',1) 
  call EXPRO_pread

  n_exp = EXPRO_n_exp

  ! r_min in m:
  r_min = EXPRO_rmin(n_exp)

  ! arho in cm
  arho  = 100*EXPRO_arho

  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_rho(:),n_exp,r,rho,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_q(:),n_exp,r,q,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_s(:),n_exp,r,s,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_kappa(:),n_exp,r,kappa,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_delta(:),n_exp,r,delta,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_skappa(:),n_exp,r,s_kappa,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_sdelta(:),n_exp,r,s_delta,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_drmaj(:),n_exp,r,shift,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_zmag(:),n_exp,r,zmag,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dzmag(:),n_exp,r,dzmag,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_zeta(:),n_exp,r,zeta,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_szeta(:),n_exp,r,s_zeta,n_r)

  ! Convert V and dV/dr from m^3 to cm^3
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_vol(:)*1e6,n_exp,r,vol,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_volp(:)*1e4,n_exp,r,volp,n_r)

  ! Convert B to Gauss (from T):
  call cub_spline(EXPRO_rmin(:)/r_min,1e4*EXPRO_bunit(:),n_exp,r,b_unit,n_r)

  ! Convert r_maj to cm (from m):
  call cub_spline(EXPRO_rmin(:)/r_min,100*EXPRO_rmaj(:),n_exp,r,r_maj,n_r)

  ! Convert T to eV (from keV) and length to cm (from m):
  call cub_spline(EXPRO_rmin(:)/r_min,1e3*EXPRO_te(:),n_exp,r,te,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,1e13*EXPRO_ne(:),n_exp,r,ne,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dlnnedr(:)/100.0,n_exp,r,dlnnedr,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dlntedr(:)/100.0,n_exp,r,dlntedr,n_r)
  do i_ion=1,loc_n_ion
     call cub_spline(EXPRO_rmin(:)/r_min,1e3*EXPRO_ti(i_ion,:),n_exp,r,ti(i_ion,:),n_r)
     call cub_spline(EXPRO_rmin(:)/r_min,1e13*EXPRO_ni(i_ion,:),n_exp,r,ni(i_ion,:),n_r)
     call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dlntidr(i_ion,:)/100.0,n_exp,r,dlntidr(i_ion,:),n_r)
     call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dlnnidr(i_ion,:)/100.0,n_exp,r,dlnnidr(i_ion,:),n_r)
  enddo

  ! QUASINEUTRALITY:
  !
  ! Overwrite main ion density and gradient with corrected density and gradient 
  ! (done in EXPRO):
  if (loc_quasineutral_flag == 1) then
     call cub_spline(EXPRO_rmin(:)/r_min,1e13*EXPRO_ni_new(:),n_exp,r,ni(1,:),n_r)
     call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dlnnidr_new(:)/100.0,n_exp,r,dlnnidr(1,:),n_r)
  endif

  ! Enforce quasineutrality (set ni based on ne,nz)
  !call tgyro_quasineutral(ni,ne,dlnnidr,dlnnedr,zi_vec,loc_n_ion,n_r)

  ! w0 (1/s)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_w0(:),n_exp,r,w0,n_r)
  ! w0p = d(w0)/dr (1/s/cm)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_w0p(:)/100.0,n_exp,r,w0p,n_r)

  ! Determine Z_eff
  if (loc_zeff_flag == 1) then
     call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_z_eff(:),n_exp,r,z_eff,n_r)
  else
     z_eff = 1.0
  endif

  ! Below, convert powers to erg/s (from MW):

  ! Exchange power
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_ei(:)*1e13,n_exp,r,p_exch_in,n_r)

  if (loc_scenario == 3) then

     ! Assume auxiliary power is contained in pow_e and pow_i 

     call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_e(:)*1e13,n_exp,r,p_e_aux_in,n_r)
     call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_i(:)*1e13,n_exp,r,p_i_aux_in,n_r)

  else

     call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_e(:)*1e13,n_exp,r,p_e_in,n_r)
     call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_pow_i(:)*1e13,n_exp,r,p_i_in,n_r)

  endif

  ! Convert flows to 1/s (from MW/keV)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_flow_beam(:)*1e22/1.6022,n_exp,r,f_b_in,n_r)
  call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_flow_wall(:)*1e22/1.6022,n_exp,r,f_w_in,n_r)

  ! Fourier coefficients for plasma shape
  if (EXPRO_nfourier > 0) then

     if (i_proc_global == 0) then
        open(unit=1,file=trim(runfile),position='append')
        write(1,*) 'INFO: (TGYRO) Passing input.profiles.geo information to components'
        write(1,*)
        close(1)
     endif

     n_fourier_geo = EXPRO_nfourier

     do n=0,n_fourier_geo
        do i=1,4  

           ! aR_n = EXPRO_geo(1,n,:)
           ! bR_n = EXPRO_geo(2,n,:)
           ! aZ_n = EXPRO_geo(3,n,:)
           ! bZ_n = EXPRO_geo(4,n,:)
           ! d(aR_n)/dr
           ! d(bR_n)/dr
           ! d(aZ_n)/dr
           ! d(bZ_n)/dr

           call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_geo(i,n,:)/r_min,&
                n_exp,r,a_fourier_geo(i,n,:),n_r)
           call cub_spline(EXPRO_rmin(:)/r_min,EXPRO_dgeo(i,n,:),&
                n_exp,r,a_fourier_geo(i+4,n,:),n_r)

        enddo
     enddo

  else

     ! Numerical equilibrium not available

     n_fourier_geo      = 0
     loc_num_equil_flag = 0

  endif

  call EXPRO_palloc(MPI_COMM_WORLD,'./',0)

  ! Convert r_min to cm (from m):
  r_min = r_min*100.0

  !----------------------------------------------------------
  ! Restore dimension to r (cm)
  r = r*r_min
  !----------------------------------------------------------

  !----------------------------------------------------------
  if (loc_restart_flag == 1) then
     call tgyro_restart
  else
     i_tran = 0
  endif
  !----------------------------------------------------------

  !----------------------------------------------------------
  ! Convert to circle if parameter set
  !
  if (loc_circ_flag == 1) then
     delta(:) = 0.0
     kappa(:) = 1.0
     shift(:) = 0.0
     s_delta(:) = 0.0
     s_kappa(:) = 0.0
     beta_unit(:) = 0.0
  endif
  !----------------------------------------------------------

  ! Now compute all profiles:
  call tgyro_profile_functions

  !----------------------------------------------------
  ! Set conditions at r=0
  !
  dlnnidr(:,1) = 0.0
  dlnnedr(1) = 0.0
  dlntidr(:,1) = 0.0
  dlntedr(1) = 0.0

  pflux_e_neo(1) = 0.0
  pflux_e_tur(1) = 0.0
  pflux_i_neo(:,1) = 0.0
  pflux_i_tur(:,1) = 0.0

  eflux_e_neo(1) = 0.0
  eflux_e_tur(1) = 0.0
  eflux_i_neo(:,1) = 0.0
  eflux_i_tur(:,1) = 0.0

  mflux_e_neo(1) = 0.0
  mflux_e_tur(1) = 0.0
  mflux_i_neo(:,1) = 0.0
  mflux_i_tur(:,1) = 0.0

  pflux_i_tot(1) = 0.0
  pflux_e_tot(1) = 0.0

  eflux_i_tot(1) = 0.0
  eflux_e_tot(1) = 0.0

  mflux_tot(1) = 0.0

  eflux_i_target(1) = 0.0
  eflux_e_target(1) = 0.0
  pflux_e_target(1) = 0.0
  mflux_target(1) = 0.0

  p_i_aux_in(1) = 0.0
  p_e_aux_in(1) = 0.0
  !----------------------------------------------------

end subroutine tgyro_init_profiles
