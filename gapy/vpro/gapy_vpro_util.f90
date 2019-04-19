subroutine vpro_compute_derived

  use vpro
  use util
  use geo

  implicit none

  integer :: n
  integer :: i
  integer :: is

  double precision, parameter :: k  = 1.6022e-12 ! erg/eV
  double precision, parameter :: mp = 1.6726e-24 ! g
  double precision, parameter :: me = 0.00027230 ! m_ele/m_deuterium (dimensionless)
  double precision, parameter :: e  = 4.8032e-10 ! statcoul
  double precision, parameter :: c  = 2.9979e10  ! cm/s
  double precision, parameter :: pi = 3.1415926535897932

  double precision, dimension(:), allocatable :: rho
  double precision, dimension(:), allocatable :: dummy
  double precision, dimension(:), allocatable :: cc
  double precision, dimension(:), allocatable :: loglam

  double precision :: r_min
  double precision :: fa,fb
  double precision :: theta(1)

  if (expro_ctrl_n_ion == -1) expro_ctrl_n_ion = expro_n_ion
   
  !---------------------------------------------------------------------
  ! Sanity checks
  !
  if (expro_ctrl_numeq_flag == 1 .and. expro_nfourier == -1) then
     print '(a)','ERROR: (expro) input.profiles.geo missing'
     stop
  endif
  if (expro_ctrl_quasineutral_flag == -1) then
     print '(a)','ERROR: (expro) expro_ctrl_quasineutral_flag not set.'
     stop
  endif
  if (expro_ctrl_numeq_flag == -1) then
     print '(a)','ERROR: (expro) expro_ctrl_numeq_flag not set.'
     stop
  endif
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Infer orientation
  ! 
  expro_signb = nint(expro_b_ref/abs(expro_b_ref))
  expro_signq = nint(expro_q(1)/abs(expro_q(1)))
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Derived quantities:
  !
  allocate(rho(expro_n_exp))
  allocate(dummy(expro_n_exp))

  rho(:) = expro_arho*expro_rho(:)

  ! b_unit
  call bound_deriv(dummy,rho**2,expro_rmin**2,expro_n_exp)
  expro_bunit(:) = expro_b_ref*dummy(:)

  ! s
  call bound_deriv(dummy,expro_q,expro_rmin,expro_n_exp)
  expro_s(:) = (expro_rmin(:)/expro_q(:))*dummy(:)

  !         d(rmaj)
  ! drmaj = -------
  !           dr
  call bound_deriv(expro_drmaj,expro_rmaj,expro_rmin,expro_n_exp)

  !         d(zmag)
  ! dzmag = -------
  !           dr
  call bound_deriv(expro_dzmag,expro_zmag,expro_rmin,expro_n_exp)

  !             r   d(kappa)
  ! s_kappa = ----- -------- 
  !           kappa    dr
  call bound_deriv(dummy,expro_kappa,expro_rmin,expro_n_exp)
  expro_skappa(:) = (expro_rmin(:)/expro_kappa(:))*dummy(:)

  !             d(delta)
  ! s_delta = r -------- 
  !                dr
  call bound_deriv(dummy,expro_delta,expro_rmin,expro_n_exp)
  expro_sdelta(:) = expro_rmin(:)*dummy(:) 

  !            d(zeta)
  ! s_zeta = r -------- 
  !              dr
  call bound_deriv(dummy,expro_zeta,expro_rmin,expro_n_exp)
  expro_szeta(:) = expro_rmin(:)*dummy(:) 

  ! 1/L_ne = -dln(ne)/dr (1/m)
  call bound_deriv(expro_dlnnedr,-log(expro_ne),expro_rmin,expro_n_exp)

  ! 1/L_Te = -dln(Te)/dr (1/m)
  call bound_deriv(expro_dlntedr,-log(expro_te),expro_rmin,expro_n_exp)

  ! NOTE: expro_sdln* will be renormalized after calculation of rhos later

  ! sne = -ne''/ne (1/m^2) [not fully normalized yet]
  call bound_deriv(expro_sdlnnedr,expro_ne*expro_dlnnedr,expro_rmin,expro_n_exp)
  expro_sdlnnedr = expro_sdlnnedr/expro_ne

  ! sTe = -Te''/Te (1/m^2) [not fully normalized yet]
  call bound_deriv(expro_sdlntedr,expro_te*expro_dlntedr,expro_rmin,expro_n_exp)
  expro_sdlntedr = expro_sdlntedr/expro_te

  expro_dlnnidr = 0.0
  expro_dlntidr = 0.0
  expro_sdlnnidr = 0.0
  expro_sdlntidr = 0.0
  
  do is=1,expro_n_ion
     if (minval(expro_ni(is,:)) > 0.0) then
        ! 1/L_ni = -dln(ni)/dr (1/m)
        call bound_deriv(expro_dlnnidr(is,:),-log(expro_ni(is,:)),expro_rmin,expro_n_exp)

        ! 1/L_Ti = -dln(Ti)/dr (1/m)
        call bound_deriv(expro_dlntidr(is,:),-log(expro_ti(is,:)),expro_rmin,expro_n_exp)

        ! sni = -ni''/ni (1/m^2) [not fully normalized yet]
        call bound_deriv(expro_sdlnnidr(is,:),expro_ni(is,:)*expro_dlnnidr(is,:),expro_rmin,expro_n_exp)
        expro_sdlnnidr(is,:) = expro_sdlnnidr(is,:)/expro_ni(is,:)

        ! sTi = -Ti''/Ti (1/m^2) [not fully normalized yet]
        call bound_deriv(expro_sdlntidr(is,:),expro_ti(is,:)*expro_dlntidr(is,:),expro_rmin,expro_n_exp)
        expro_sdlntidr(is,:) = expro_sdlntidr(is,:)/expro_ti(is,:)
     endif
  enddo

  ! 1/L_Ptot = -dln(Ptot)/dr (1/m)
  if (minval(expro_ptot) > 0.0) then
     call bound_deriv(expro_dlnptotdr,-log(expro_ptot),expro_rmin,expro_n_exp)
  else
     expro_dlnptotdr = 0.0
  endif

  ! dr/d(rho)
  call bound_deriv(expro_drdrho,expro_rmin,rho,expro_n_exp)
  !--------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Fourier coefficients for geometry (if they exist)
  !
  if (expro_nfourier > 0) then
     do n=0,expro_nfourier
        do i=1,4  

           ! aR_n = expro_geo(1,n,:)
           ! bR_n = expro_geo(2,n,:)
           ! aZ_n = expro_geo(3,n,:)
           ! bZ_n = expro_geo(4,n,:)

           ! d(aR_n)/dr, d(bR_n)/dr, d(aZ_n)/dr, d(bZ_n)/dr

           call bound_deriv(expro_dgeo(i,n,:),expro_geo(i,n,:),&
                expro_rmin,expro_n_exp)
        enddo
     enddo
  endif

  !-------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Geometry factors: 
  !
  ! - w0, w0p, vol, volp
  !
  GEO_nfourier_in = expro_nfourier
  GEO_signb_in    = expro_signb

  r_min = expro_rmin(expro_n_exp)

  do i=2,expro_n_exp

     ! Parameters to be passed to GEO library   
     !
     ! NOTE: dp/dr set to zero without loss of generality.
     ! 
     GEO_rmin_in      = expro_rmin(i)/r_min
     GEO_rmaj_in      = expro_rmaj(i)/r_min
     GEO_drmaj_in     = expro_drmaj(i)
     GEO_zmag_in      = expro_zmag(i)/r_min
     GEO_dzmag_in     = expro_dzmag(i)
     GEO_q_in         = expro_q(i)
     GEO_s_in         = expro_s(i)
     GEO_kappa_in     = expro_kappa(i)
     GEO_s_kappa_in   = expro_skappa(i)
     GEO_delta_in     = expro_delta(i)
     GEO_s_delta_in   = expro_sdelta(i)
     GEO_zeta_in      = expro_zeta(i)
     GEO_s_zeta_in    = expro_szeta(i)
     GEO_beta_star_in = 0.0
     !
     theta(1) = 0.0
     if (expro_ctrl_numeq_flag == 0) then
        ! Call GEO with model shape
        GEO_model_in = 0
        call GEO_interp(1,theta,.true.)
     else
        ! Call GEO with general (numerical) shape
        GEO_model_in = 1
        GEO_fourier_in(1:4,0:GEO_nfourier_in) = expro_geo(:,:,i)/r_min
        GEO_fourier_in(5:8,0:GEO_nfourier_in) = expro_dgeo(:,:,i)
        call GEO_interp(1,theta,.true.)
        if (minval(GEOV_jac_r) <= 0.0) then
           print '(a,i3,a)','WARNING: (expro) Negative GEO Jacobian for i =',i,' in input.profiles'
        endif
     endif

     ! V and dV/dr
     expro_volp(i) = GEO_volume_prime*r_min**2
     expro_vol(i)  = GEO_volume*r_min**3

     ! |grad r| at theta=0
     expro_grad_r0(i) = GEO_grad_r0

     ! <|grad r|> 
     expro_ave_grad_r(i) = GEO_fluxsurfave_grad_r

     ! B_poloidal and B_toroidal [T] at theta=0
     expro_bp0(i) = GEO_bp(1)*expro_bunit(i)
     expro_bt0(i) = GEO_bt(1)*expro_bunit(i)

     expro_thetascale(i) = GEO_thetascale

     ! Plasma current [A] I = (1/mu0) Int[Bp dl] 
     expro_ip(i) = 7.958e5*(GEO_bl*r_min*expro_bunit(i))

  enddo

  !--------------------------------------------------------------
  ! Extrapolate some quantities to axis:
  !
  call bound_extrap(fa,fb,expro_grad_r0,expro_rmin,expro_n_exp)
  expro_grad_r0(1) = fa

  call bound_extrap(fa,fb,expro_ave_grad_r,expro_rmin,expro_n_exp)
  expro_ave_grad_r(1) = fa

  call bound_extrap(fa,fb,expro_bp0,expro_rmin,expro_n_exp)
  expro_bp0(1) = fa

  call bound_extrap(fa,fb,expro_bt0,expro_rmin,expro_n_exp)
  expro_bt0(1) = fa
  !
  ! Both V and dV/dr are zero on axis.
  !
  expro_vol(1)  = 0.0
  expro_volp(1) = 0.0  
  expro_ip(1)   = 0.0
  expro_thetascale(1) = expro_thetascale(2)

  !--------------------------------------------------------------

  !-----------------------------------------------------------------
  ! CGS calculation of deuterium sound speed (cm/s) and 
  ! deuterium gyroradius (cm)
  !
  expro_cs(:)   = sqrt( k*(1e3*expro_te(:))/(2.0*mp) )   
  expro_rhos(:) = expro_cs(:)/(e*(1e4*expro_bunit(:))/(2.0*mp*c))
  ! 
  ! Convert to m/s and m:
  !
  expro_cs(:)   = expro_cs(:)/100.0
  expro_rhos(:) = expro_rhos(:)/100.0
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Renormalize shearing parameters
  !
  ! sn = -n''/n*rhos (1/m) 
  ! sT = -T''/T*rhos (1/m)
  expro_sdlnnedr = expro_sdlnnedr*expro_rhos(:)
  expro_sdlntedr = expro_sdlntedr*expro_rhos(:)
  do is=1,expro_n_ion
     expro_sdlnnidr(is,:) = expro_sdlnnidr(is,:)*expro_rhos(:)
     expro_sdlntidr(is,:) = expro_sdlntidr(is,:)*expro_rhos(:)
  enddo
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Compute the electron-electron collision frequency (1/s)
  allocate(cc(expro_n_exp))
  allocate(loglam(expro_n_exp))
  cc(:) = sqrt(2.0) * pi * 1.6022**4 * 1.0 / (4.0 * pi * 8.8542)**2 &
       * 1.0 / (sqrt(3.3452) * 1602.2**1.5) * 1e9
  loglam(:) = 24.0 - log(sqrt(expro_ne(:)*1e13)/(expro_te(:)*1e3))
  expro_nuee(:) = cc(:) * loglam(:) * expro_ne(:) &
       / (sqrt(me) * expro_te(:)**1.5)
  deallocate(cc)
  deallocate(loglam)
  !-----------------------------------------------------------------

  !--------------------------------------------------------------
  ! Compute w0p, gamma_e, gamma_p and mach:
  !
  call bound_deriv(expro_w0p,expro_w0,expro_rmin,expro_n_exp)
  !  
  expro_gamma_e(:) = -expro_rmin(:)/expro_q(:)*expro_w0p(:)
  expro_gamma_p(:) = -expro_rmaj(:)*expro_w0p(:)
  expro_mach(:)    = expro_rmaj(:)*expro_w0(:)/expro_cs(:)
  !--------------------------------------------------------------

  ! Clean up
  deallocate(rho)

  ! Density profile control

  if (expro_ctrl_quasineutral_flag == 1) then

     expro_ni_new(:) = 0.0
     do is=2,expro_ctrl_n_ion
        expro_ni_new(:) = expro_ni_new(:)+expro_z(is)*expro_ni(is,:)
     enddo
     expro_ni_new(:) = (expro_ne(:)-expro_ni_new(:))/expro_z(1)

     ! 1/L_ni = -dln(ni)/dr (1/m)
     call bound_deriv(expro_dlnnidr_new(:),-log(expro_ni_new(:)),&
          expro_rmin,expro_n_exp)

     ! sni = -ni''/ni (1/m^2)
     call bound_deriv(expro_sdlnnidr_new(:),expro_ni_new(:)*expro_dlnnidr_new(:),expro_rmin,expro_n_exp)
     expro_sdlnnidr_new(:) = expro_sdlnnidr_new(:)/expro_ni_new(:)*expro_rhos(:)

     !if (minval(expro_ni_new(:)) <= 0.0) then
     !   expro_error = 1
     !endif

  else

     expro_ni_new(:) = expro_ni(1,:)
     expro_dlnnidr_new(:) = expro_dlnnidr(1,:)
     expro_sdlnnidr_new(:) = expro_sdlnnidr(1,:)

  endif

  deallocate(dummy)

  !do is=1,expro_ctrl_n_ion
  !   if (minval(expro_ni(is,:)) <= 0.0) then
  !      expro_error=1
  !   endif
  !enddo

end subroutine vpro_compute_derived
