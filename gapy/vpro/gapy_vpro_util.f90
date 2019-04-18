subroutine vpro_compute_derived

  use vpro
  use util
  use geo

  implicit none

  integer :: n
  integer :: i
  integer :: is

  real, parameter :: k  = 1.6022e-12 ! erg/eV
  real, parameter :: mp = 1.6726e-24 ! g
  real, parameter :: me = 0.00027230 ! m_ele/m_deuterium (dimensionless)
  real, parameter :: e  = 4.8032e-10 ! statcoul
  real, parameter :: c  = 2.9979e10  ! cm/s
  real, parameter :: pi = 3.1415926535897932

  real, dimension(:), allocatable :: rho
  real, dimension(:), allocatable :: dummy
  real, dimension(:), allocatable :: cc
  real, dimension(:), allocatable :: loglam

  real :: r_min
  real :: fa,fb
  real :: theta(1)

  EXPRO_ctrl_n_ion = EXPRO_n_ion
  EXPRO_ctrl_numeq_flag = 0
  EXPRO_nfourier = -1
  EXPRO_ctrl_quasineutral_flag = 0
  
  !---------------------------------------------------------------------
  ! Infer orientation
  ! 
  EXPRO_signb = nint(EXPRO_b_ref/abs(EXPRO_b_ref))
  EXPRO_signq = nint(EXPRO_q(1)/abs(EXPRO_q(1)))
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Derived quantities:
  !
  allocate(rho(EXPRO_n_exp))
  allocate(dummy(EXPRO_n_exp))

  rho(:) = EXPRO_arho*EXPRO_rho(:)

  ! b_unit
  call bound_deriv(dummy,rho**2,EXPRO_rmin**2,EXPRO_n_exp)
  EXPRO_bunit(:) = EXPRO_b_ref*dummy(:)

  ! s
  call bound_deriv(dummy,EXPRO_q,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_s(:) = (EXPRO_rmin(:)/EXPRO_q(:))*dummy(:)

  !         d(rmaj)
  ! drmaj = -------
  !           dr
  call bound_deriv(EXPRO_drmaj,EXPRO_rmaj,EXPRO_rmin,EXPRO_n_exp)

  !         d(zmag)
  ! dzmag = -------
  !           dr
  call bound_deriv(EXPRO_dzmag,EXPRO_zmag,EXPRO_rmin,EXPRO_n_exp)

  !             r   d(kappa)
  ! s_kappa = ----- -------- 
  !           kappa    dr
  call bound_deriv(dummy,EXPRO_kappa,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_skappa(:) = (EXPRO_rmin(:)/EXPRO_kappa(:))*dummy(:)

  !             d(delta)
  ! s_delta = r -------- 
  !                dr
  call bound_deriv(dummy,EXPRO_delta,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_sdelta(:) = EXPRO_rmin(:)*dummy(:) 

  !            d(zeta)
  ! s_zeta = r -------- 
  !              dr
  call bound_deriv(dummy,EXPRO_zeta,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_szeta(:) = EXPRO_rmin(:)*dummy(:) 

  ! 1/L_ne = -dln(ne)/dr (1/m)
  call bound_deriv(EXPRO_dlnnedr,-log(EXPRO_ne),EXPRO_rmin,EXPRO_n_exp)

  ! 1/L_Te = -dln(Te)/dr (1/m)
  call bound_deriv(EXPRO_dlntedr,-log(EXPRO_te),EXPRO_rmin,EXPRO_n_exp)

  ! NOTE: EXPRO_sdln* will be renormalized after calculation of rhos later

  ! sne = -ne''/ne (1/m^2) [not fully normalized yet]
  call bound_deriv(EXPRO_sdlnnedr,EXPRO_ne*EXPRO_dlnnedr,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_sdlnnedr = EXPRO_sdlnnedr/EXPRO_ne

  ! sTe = -Te''/Te (1/m^2) [not fully normalized yet]
  call bound_deriv(EXPRO_sdlntedr,EXPRO_te*EXPRO_dlntedr,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_sdlntedr = EXPRO_sdlntedr/EXPRO_te

  EXPRO_dlnnidr = 0.0
  EXPRO_dlntidr = 0.0
  EXPRO_sdlnnidr = 0.0
  EXPRO_sdlntidr = 0.0
  
  do is=1,EXPRO_n_ion
     if (minval(EXPRO_ni(is,:)) > 0.0) then
        ! 1/L_ni = -dln(ni)/dr (1/m)
        call bound_deriv(EXPRO_dlnnidr(is,:),-log(EXPRO_ni(is,:)),EXPRO_rmin,EXPRO_n_exp)

        ! 1/L_Ti = -dln(Ti)/dr (1/m)
        call bound_deriv(EXPRO_dlntidr(is,:),-log(EXPRO_ti(is,:)),EXPRO_rmin,EXPRO_n_exp)

        ! sni = -ni''/ni (1/m^2) [not fully normalized yet]
        call bound_deriv(EXPRO_sdlnnidr(is,:),EXPRO_ni(is,:)*EXPRO_dlnnidr(is,:),EXPRO_rmin,EXPRO_n_exp)
        EXPRO_sdlnnidr(is,:) = EXPRO_sdlnnidr(is,:)/EXPRO_ni(is,:)

        ! sTi = -Ti''/Ti (1/m^2) [not fully normalized yet]
        call bound_deriv(EXPRO_sdlntidr(is,:),EXPRO_ti(is,:)*EXPRO_dlntidr(is,:),EXPRO_rmin,EXPRO_n_exp)
        EXPRO_sdlntidr(is,:) = EXPRO_sdlntidr(is,:)/EXPRO_ti(is,:)
     endif
  enddo

  ! 1/L_Ptot = -dln(Ptot)/dr (1/m)
  if (minval(EXPRO_ptot) > 0.0) then
     call bound_deriv(EXPRO_dlnptotdr,-log(EXPRO_ptot),EXPRO_rmin,EXPRO_n_exp)
  else
     EXPRO_dlnptotdr = 0.0
  endif

  ! dr/d(rho)
  call bound_deriv(EXPRO_drdrho,EXPRO_rmin,rho,EXPRO_n_exp)
  !--------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Fourier coefficients for geometry (if they exist)
  !
  if (EXPRO_nfourier > 0) then
     do n=0,EXPRO_nfourier
        do i=1,4  

           ! aR_n = EXPRO_geo(1,n,:)
           ! bR_n = EXPRO_geo(2,n,:)
           ! aZ_n = EXPRO_geo(3,n,:)
           ! bZ_n = EXPRO_geo(4,n,:)

           ! d(aR_n)/dr, d(bR_n)/dr, d(aZ_n)/dr, d(bZ_n)/dr

           call bound_deriv(EXPRO_dgeo(i,n,:),EXPRO_geo(i,n,:),&
                EXPRO_rmin,EXPRO_n_exp)
        enddo
     enddo
  endif

  !-------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Geometry factors: 
  !
  ! - w0, w0p, vol, volp
  !
  GEO_nfourier_in = EXPRO_nfourier
  GEO_signb_in    = EXPRO_signb

  r_min = EXPRO_rmin(EXPRO_n_exp)

  do i=2,EXPRO_n_exp

     ! Parameters to be passed to GEO library   
     !
     ! NOTE: dp/dr set to zero without loss of generality.
     ! 
     GEO_rmin_in      = EXPRO_rmin(i)/r_min
     GEO_rmaj_in      = EXPRO_rmaj(i)/r_min
     GEO_drmaj_in     = EXPRO_drmaj(i)
     GEO_zmag_in      = EXPRO_zmag(i)/r_min
     GEO_dzmag_in     = EXPRO_dzmag(i)
     GEO_q_in         = EXPRO_q(i)
     GEO_s_in         = EXPRO_s(i)
     GEO_kappa_in     = EXPRO_kappa(i)
     GEO_s_kappa_in   = EXPRO_skappa(i)
     GEO_delta_in     = EXPRO_delta(i)
     GEO_s_delta_in   = EXPRO_sdelta(i)
     GEO_zeta_in      = EXPRO_zeta(i)
     GEO_s_zeta_in    = EXPRO_szeta(i)
     GEO_beta_star_in = 0.0
     !
     theta(1) = 0.0
     if (EXPRO_ctrl_numeq_flag == 0) then
        ! Call GEO with model shape
        GEO_model_in = 0
        call GEO_interp(1,theta,.true.)
     else
        ! Call GEO with general (numerical) shape
        GEO_model_in = 1
        GEO_fourier_in(1:4,0:GEO_nfourier_in) = EXPRO_geo(:,:,i)/r_min
        GEO_fourier_in(5:8,0:GEO_nfourier_in) = EXPRO_dgeo(:,:,i)
        call GEO_interp(1,theta,.true.)
        if (minval(GEOV_jac_r) <= 0.0) then
           print '(a,i3,a)','WARNING: (EXPRO) Negative GEO Jacobian for i =',i,' in input.profiles'
        endif
     endif

     ! V and dV/dr
     EXPRO_volp(i) = GEO_volume_prime*r_min**2
     EXPRO_vol(i)  = GEO_volume*r_min**3

     ! |grad r| at theta=0
     EXPRO_grad_r0(i) = GEO_grad_r0

     ! <|grad r|> 
     EXPRO_ave_grad_r(i) = GEO_fluxsurfave_grad_r

     ! B_poloidal and B_toroidal [T] at theta=0
     EXPRO_bp0(i) = GEO_bp(1)*EXPRO_bunit(i)
     EXPRO_bt0(i) = GEO_bt(1)*EXPRO_bunit(i)

     EXPRO_thetascale(i) = GEO_thetascale

     ! Plasma current [A] I = (1/mu0) Int[Bp dl] 
     EXPRO_ip(i) = 7.958e5*(GEO_bl*r_min*EXPRO_bunit(i))

  enddo

  !--------------------------------------------------------------
  ! Extrapolate some quantities to axis:
  !
  call bound_extrap(fa,fb,EXPRO_grad_r0,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_grad_r0(1) = fa

  call bound_extrap(fa,fb,EXPRO_ave_grad_r,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_ave_grad_r(1) = fa

  call bound_extrap(fa,fb,EXPRO_bp0,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_bp0(1) = fa

  call bound_extrap(fa,fb,EXPRO_bt0,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_bt0(1) = fa
  !
  ! Both V and dV/dr are zero on axis.
  !
  EXPRO_vol(1)  = 0.0
  EXPRO_volp(1) = 0.0  
  EXPRO_ip(1)   = 0.0
  EXPRO_thetascale(1) = EXPRO_thetascale(2)

  !--------------------------------------------------------------

  !-----------------------------------------------------------------
  ! CGS calculation of deuterium sound speed (cm/s) and 
  ! deuterium gyroradius (cm)
  !
  EXPRO_cs(:)   = sqrt( k*(1e3*EXPRO_te(:))/(2.0*mp) )   
  EXPRO_rhos(:) = EXPRO_cs(:)/(e*(1e4*EXPRO_bunit(:))/(2.0*mp*c))
  ! 
  ! Convert to m/s and m:
  !
  EXPRO_cs(:)   = EXPRO_cs(:)/100.0
  EXPRO_rhos(:) = EXPRO_rhos(:)/100.0
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Renormalize shearing parameters
  !
  ! sn = -n''/n*rhos (1/m) 
  ! sT = -T''/T*rhos (1/m)
  EXPRO_sdlnnedr = EXPRO_sdlnnedr*EXPRO_rhos(:)
  EXPRO_sdlntedr = EXPRO_sdlntedr*EXPRO_rhos(:)
  do is=1,EXPRO_n_ion
     EXPRO_sdlnnidr(is,:) = EXPRO_sdlnnidr(is,:)*EXPRO_rhos(:)
     EXPRO_sdlntidr(is,:) = EXPRO_sdlntidr(is,:)*EXPRO_rhos(:)
  enddo
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Compute the electron-electron collision frequency (1/s)
  allocate(cc(EXPRO_n_exp))
  allocate(loglam(EXPRO_n_exp))
  cc(:) = sqrt(2.0) * pi * 1.6022**4 * 1.0 / (4.0 * pi * 8.8542)**2 &
       * 1.0 / (sqrt(3.3452) * 1602.2**1.5) * 1e9
  loglam(:) = 24.0 - log(sqrt(EXPRO_ne(:)*1e13)/(EXPRO_te(:)*1e3))
  EXPRO_nuee(:) = cc(:) * loglam(:) * EXPRO_ne(:) &
       / (sqrt(me) * EXPRO_te(:)**1.5)
  deallocate(cc)
  deallocate(loglam)
  !-----------------------------------------------------------------

  !--------------------------------------------------------------
  ! Compute w0p, gamma_e, gamma_p and mach:
  !
  call bound_deriv(EXPRO_w0p,EXPRO_w0,EXPRO_rmin,EXPRO_n_exp)
  !  
  EXPRO_gamma_e(:) = -EXPRO_rmin(:)/EXPRO_q(:)*EXPRO_w0p(:)
  EXPRO_gamma_p(:) = -EXPRO_rmaj(:)*EXPRO_w0p(:)
  EXPRO_mach(:)    = EXPRO_rmaj(:)*EXPRO_w0(:)/EXPRO_cs(:)
  !--------------------------------------------------------------

  ! Clean up
  deallocate(rho)

  ! Density profile control

  if (EXPRO_ctrl_quasineutral_flag == 1) then

     EXPRO_ni_new(:) = 0.0
     do is=2,EXPRO_ctrl_n_ion
        EXPRO_ni_new(:) = EXPRO_ni_new(:)+EXPRO_z(is)*EXPRO_ni(is,:)
     enddo
     EXPRO_ni_new(:) = (EXPRO_ne(:)-EXPRO_ni_new(:))/EXPRO_z(1)

     ! 1/L_ni = -dln(ni)/dr (1/m)
     call bound_deriv(EXPRO_dlnnidr_new(:),-log(EXPRO_ni_new(:)),&
          EXPRO_rmin,EXPRO_n_exp)

     ! sni = -ni''/ni (1/m^2)
     call bound_deriv(EXPRO_sdlnnidr_new(:),EXPRO_ni_new(:)*EXPRO_dlnnidr_new(:),EXPRO_rmin,EXPRO_n_exp)
     EXPRO_sdlnnidr_new(:) = EXPRO_sdlnnidr_new(:)/EXPRO_ni_new(:)*EXPRO_rhos(:)

     !if (minval(EXPRO_ni_new(:)) <= 0.0) then
     !   EXPRO_error = 1
     !endif

  else

     EXPRO_ni_new(:) = EXPRO_ni(1,:)
     EXPRO_dlnnidr_new(:) = EXPRO_dlnnidr(1,:)
     EXPRO_sdlnnidr_new(:) = EXPRO_sdlnnidr(1,:)

  endif

  deallocate(dummy)

  !do is=1,EXPRO_ctrl_n_ion
  !   if (minval(EXPRO_ni(is,:)) <= 0.0) then
  !      EXPRO_error=1
  !   endif
  !enddo

end subroutine vpro_compute_derived

!----------------------------------------------------------------
! vpro_locsim_profiles.f90
!
! PURPOSE:
!  Read experimental profiles and generate local profile 
!  parameters at rmin = r/a.
!
! INPUTS:
!  path              : path to data
!  comm              : MPI communicator
!  numeq_flag        : Fourier series equilibrium (0=no,1=yes)
!  udsymmetry_flag   : enforce up-down symmetry (0=no,1=yes)
!  quasineutral_flag : enforce quasineutrality (0=no,1=yes)
!  n_species_in      : total species (e+i)
!  z                 : vector of charges (length n_species_in-1)
!  rmin              : r/a
!
! OUTPUTS:
!  btccw             : (+1 or -1) 
!  ipccw             : (+1 or -1)
!  a_meters          : minor radius i metres
!
! OUTPUTS (interface):
!  real :: shift_loc
!  real :: q_loc
!  real :: s_loc
!  real :: kappa_loc
!  real :: delta_loc
!  real :: zeta_loc
!  real :: s_kappa_loc
!  real :: s_delta_loc
!  real :: s_zeta_loc
!  real :: zmag_loc
!  real :: dzmag_loc
!  real :: gamma_e_loc
!  real :: gamma_p_loc
!  real :: mach_loc
!  real :: rmaj_loc
!  real :: rhos_loc [m]
!  real :: z_eff_loc
!  real :: b_unit_loc
!  real :: rho_norm_loc
!  real :: psi_norm_loc
!  real :: psi_a_loc
!  real :: cs_loc
!  real :: beta_star_loc
!
!  real, dimension(9) :: dens_loc
!  real, dimension(9) :: temp_loc
!  real, dimension(9) :: dlnndr_loc
!  real, dimension(9) :: dlntdr_loc
!  real, dimension(9) :: sdlnndr_loc
!  real, dimension(9) :: sdlntdr_loc
!----------------------------------------------------------------

subroutine vpro_locsim_profiles(&
     path,&
     comm,&
     numeq_flag,&
     udsymmetry_flag,&
     quasineutral_flag,&
     n_species_in,&
     rmin,&
     btccw,&
     ipccw,&
     a_meters)
  
  use vpro

  implicit none

  character(len=*), intent(in) :: path 
  integer, intent(in) :: comm
  integer, intent(in) :: numeq_flag
  integer, intent(in) :: udsymmetry_flag
  integer, intent(in) :: quasineutral_flag
  integer, intent(in) :: n_species_in
  real, intent(in) :: rmin
  real, intent(inout) :: btccw,ipccw,a_meters
  real, parameter :: pi=3.14159265358979323846
  integer :: i,j,i_ion

  rmin_loc = rmin
  
  n_species_exp = n_species_in

  !--------------------------------------------------------------
  ! use EXPRO routines to read data:
  !
  EXPRO_ctrl_quasineutral_flag = 1  ! quasi-neutrality density flag
  EXPRO_ctrl_numeq_flag = numeq_flag
  EXPRO_ctrl_n_ion = n_species_exp-1

  if (comm == -1) then
     call EXPRO_alloc(path,1)
     call EXPRO_read
  else
     call EXPRO_palloc(comm,path,1)
     call EXPRO_pread
  endif
  call EXPRO_locsim_alloc(1)
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! Transfer data from read vector to individual arrays:
  !
  btccw = -EXPRO_signb
  ipccw = -EXPRO_signq*EXPRO_signb

  rmin_exp(:) = EXPRO_rmin(:)

  if (udsymmetry_flag == 1) then
     EXPRO_zmag(:) = 0.0   
     EXPRO_dzmag(:) = 0.0
  endif

  ! Minor radius, a, in meters:
  a_meters = rmin_exp(EXPRO_n_exp)

  rmin_exp(:) = rmin_exp(:)/a_meters

  ! Pack electrons into top of species vector.
  temp_exp(n_species_exp,:)    = EXPRO_te(:)
  dlntdr_exp(n_species_exp,:)  = EXPRO_dlntedr(:)*a_meters 
  sdlntdr_exp(n_species_exp,:) = EXPRO_sdlntedr(:)*a_meters
  dens_exp(n_species_exp,:)    = EXPRO_ne(:)
  dlnndr_exp(n_species_exp,:)  = EXPRO_dlnnedr(:)*a_meters 
  sdlnndr_exp(n_species_exp,:) = EXPRO_sdlnnedr(:)*a_meters

  mass_loc(n_species_exp) = 1.0/1837
  z_loc(n_species_exp) = -1.0
 
  ! Pack ions from the bottom
  do i_ion=1,n_species_exp-1
     ! ion temps should be equal, but not enforced 
     temp_exp(i_ion,:)    = EXPRO_ti(i_ion,:)
     dlntdr_exp(i_ion,:)  = EXPRO_dlntidr(i_ion,:)*a_meters 
     sdlntdr_exp(i_ion,:) = EXPRO_sdlntidr(i_ion,:)*a_meters 

     mass_loc(i_ion) = EXPRO_mass(i_ion)
     z_loc(i_ion) = EXPRO_z(i_ion)
   
     ! First species density is reset by quasi-neutrality
     if (quasineutral_flag == 1 .and. i_ion == 1) then
        dens_exp(i_ion,:)    = EXPRO_ni_new(:)
        dlnndr_exp(i_ion,:)  = EXPRO_dlnnidr_new(:)*a_meters
        sdlnndr_exp(i_ion,:) = EXPRO_sdlnnidr_new(:)*a_meters
     else
        dens_exp(i_ion,:)    = EXPRO_ni(i_ion,:)
        dlnndr_exp(i_ion,:)  = EXPRO_dlnnidr(i_ion,:)*a_meters
        sdlnndr_exp(i_ion,:) = EXPRO_sdlnnidr(i_ion,:)*a_meters
     endif
  enddo

  ! Sanity check for densities
  do i=1,n_species_exp
     if (minval(dens_exp(i,:)) <= 0.0) then
        !call cgyro_error('Nonpositive in exp. density profile')
        return
     endif
  enddo

  ! Rotation
  gamma_e_exp(:) = -EXPRO_w0p(:)*(a_meters*rmin_exp(:))/EXPRO_q(:)
  gamma_p_exp(:) = -EXPRO_w0p(:)*EXPRO_rmaj(:)
  mach_exp(:)    = EXPRO_w0(:)*EXPRO_rmaj(:)

  !------------------------------------------------------------------
  ! Use local cubic spline interpolation to get simulation 
  ! profiles from experimental (_exp) ones.
  ! 
  call cub_spline(rmin_exp,EXPRO_rmaj/a_meters,EXPRO_n_exp,rmin,rmaj_loc,1)
  call cub_spline(rmin_exp,EXPRO_q,EXPRO_n_exp,rmin,q_loc,1)
  call cub_spline(rmin_exp,EXPRO_s,EXPRO_n_exp,rmin,s_loc,1)
  call cub_spline(rmin_exp,EXPRO_drmaj,EXPRO_n_exp,rmin,shift_loc,1)
  call cub_spline(rmin_exp,EXPRO_kappa,EXPRO_n_exp,rmin,kappa_loc,1)
  call cub_spline(rmin_exp,EXPRO_skappa,EXPRO_n_exp,rmin,s_kappa_loc,1)
  call cub_spline(rmin_exp,EXPRO_delta,EXPRO_n_exp,rmin,delta_loc,1)
  call cub_spline(rmin_exp,EXPRO_sdelta,EXPRO_n_exp,rmin,s_delta_loc,1)
  call cub_spline(rmin_exp,EXPRO_zeta,EXPRO_n_exp,rmin,zeta_loc,1)
  call cub_spline(rmin_exp,EXPRO_szeta,EXPRO_n_exp,rmin,s_zeta_loc,1)
  call cub_spline(rmin_exp,EXPRO_zmag/a_meters,EXPRO_n_exp,rmin,zmag_loc,1)
  call cub_spline(rmin_exp,EXPRO_dzmag,EXPRO_n_exp,rmin,dzmag_loc,1)
  call cub_spline(rmin_exp,gamma_e_exp,EXPRO_n_exp,rmin,gamma_e_loc,1)
  call cub_spline(rmin_exp,gamma_p_exp,EXPRO_n_exp,rmin,gamma_p_loc,1)
  call cub_spline(rmin_exp,mach_exp,EXPRO_n_exp,rmin,mach_loc,1)
  call cub_spline(rmin_exp,EXPRO_rhos,EXPRO_n_exp,rmin,rhos_loc,1)
  call cub_spline(rmin_exp,EXPRO_cs,EXPRO_n_exp,rmin,cs_loc,1)
  call cub_spline(rmin_exp,EXPRO_z_eff,EXPRO_n_exp,rmin,z_eff_loc,1)
  call cub_spline(rmin_exp,EXPRO_bunit,EXPRO_n_exp,rmin,b_unit_loc,1)
  call cub_spline(rmin_exp,EXPRO_rho,EXPRO_n_exp,rmin,rho_norm_loc,1)
  call cub_spline(rmin_exp,EXPRO_polflux,EXPRO_n_exp,rmin,psi_norm_loc,1)
  psi_norm_loc = psi_norm_loc/EXPRO_polflux(EXPRO_n_exp)
  psi_a_loc = EXPRO_polflux(EXPRO_n_exp)


  beta_star_loc = 0.0  
  do i=1,n_species_exp
     ! Note: mapping is only done for n_species (not n_species_exp)
     call cub_spline(rmin_exp,dens_exp(i,:),EXPRO_n_exp,rmin,dens_loc(i),1)
     call cub_spline(rmin_exp,temp_exp(i,:),EXPRO_n_exp,rmin,temp_loc(i),1)
     call cub_spline(rmin_exp,dlntdr_exp(i,:),EXPRO_n_exp,rmin,dlntdr_loc(i),1)
     call cub_spline(rmin_exp,dlnndr_exp(i,:),EXPRO_n_exp,rmin,dlnndr_loc(i),1)
     call cub_spline(rmin_exp,sdlntdr_exp(i,:),EXPRO_n_exp,rmin,sdlntdr_loc(i),1)
     call cub_spline(rmin_exp,sdlnndr_exp(i,:),EXPRO_n_exp,rmin,sdlnndr_loc(i),1)
     beta_star_loc = beta_star_loc+dens_loc(i)*temp_loc(i)*(dlnndr_loc(i)+dlntdr_loc(i))
  enddo
  ! CGS beta calculation
  betae_loc = 4.027e-3*dens_loc(n_species_exp)*temp_loc(n_species_exp)/b_unit_loc**2

  beta_star_loc = beta_star_loc*betae_loc/(dens_loc(n_species_exp)*temp_loc(n_species_exp))

  if (numeq_flag == 1) then

     if (EXPRO_nfourier <= 0) then
        return
     endif

     geo_ny_loc = EXPRO_nfourier
     allocate(geo_yin_exp(8,0:geo_ny_loc,EXPRO_n_exp))
     if(allocated(geo_yin_loc)) deallocate(geo_yin_loc)
     allocate(geo_yin_loc(8,0:geo_ny_loc))
     geo_yin_exp(1:4,:,:) = EXPRO_geo(:,:,:)/a_meters
     geo_yin_exp(5:8,:,:) = EXPRO_dgeo(:,:,:)

     do i=1,8
        do j=0,geo_ny_loc
           call cub_spline(rmin_exp,geo_yin_exp(i,j,:),EXPRO_n_exp,rmin, &
                geo_yin_loc(i,j),1)
        enddo
     enddo

  endif

  if (comm == -1) then
     call EXPRO_alloc(path,0)
  else
     call EXPRO_palloc(comm,path,0)
  endif
  call EXPRO_locsim_alloc(0)

end subroutine vpro_locsim_profiles
