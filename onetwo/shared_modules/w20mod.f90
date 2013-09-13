module w20mod

!                 >>>>>   N O T E   <<<<<
! {{{ and }}} are Vim folding markers. DO NOT REMOVE THEM.
! Vim is an improved version of the Vi text editor.

!----- DESCRIPTION -------------------------------------------------{{{
!
! This Fortran 90 module contains the Weiland20 anomalous transport
! subroutine w20main and its supporting codes. 
!
! This module contains a single public subroutines:
!
!    * Subroutine w20main
!
! Detailed description can be found in each subroutine.
!
! Revision History
! ----------------
! Apr 7, 2011   Lixiang Luo
!               Improved diagnostic output handling
!
! Mar 1, 2011   Lixiang Luo
!               Repackaging as a Fortran 90 module
!               V7.1 Original Release
!
! 2008          Federico Halpern
!               Original codes of MMM08
!
!-------------------------------------------------------------------}}}


! Module specification {{{
USE nrtype,                                     ONLY : DP,I4B
implicit none

private 

!Precision constants
INTEGER, PARAMETER :: R8 = SELECTED_REAL_KIND(12,100)

!Order of the dispersion relation
integer, parameter :: neq = 9

!Flags used for searching ion or electron drift modes
integer :: searchmode
integer, parameter :: S_ION = 2
integer, parameter :: S_ELC = 3
integer, parameter :: S_ANY = 0
integer, parameter :: N_NONE = 0
real(DP), parameter :: grdmin = 0.001_DP

!Printing flag
integer :: lprint
integer :: hprint      ! Unit number for diagnostic output

!Species mass
real(DP) :: aimp    ! Atomic mass a_z (a.m.u.)
real(DP) :: ahyd    ! Atomic mass a_i (a.m.u.)

!Species charge
real(DP) :: zimp    ! Impurity charge Z_z
real(DP) :: zhyd    ! Hydrogenic charge Z_i
real(DP) :: zfion   ! Fast ion charge Z_f
real(DP) :: zeff    ! Effective charge

!Density and temperature gradients
!Defined as -(R/U)(dU/dr), derivative respect to minor radius r (not rho)
real(DP) :: gte     ! Electron temperature gradient 
real(DP) :: gti     ! Ion temperature gradient
real(DP) :: gtz     ! Impurity temperature gradient
real(DP) :: gne     ! Electron density gradient
real(DP) :: gni     ! Ion density gradient
real(DP) :: gnz     ! Impurity density gradient
real(DP) :: gvt     ! Toroidal velocity gradient
real(DP) :: gvp     ! Poloidal velocity gradient
real(DP) :: gvpar   ! Parallel velocity gradient

!Correlation lengths
real(DP) :: kyrho   ! k_y rho = 0.316, ( k_y rho )^2 = 0.1
real(DP) :: flh     ! Hydrogenic modes correlation length ( k_y rho )**2~0.10
real(DP) :: fle     ! Electron modes correlation length
real(DP) :: zflh    ! Atomic mass * Correlation length / hydrogenic charge squared
real(DP) :: zflz    ! Same, for impurities
real(DP) :: epsncorr ! epsilon_n correlation factor ( \epsilon_n = 2 L_n / L_B = 2 / gne )

!Electron temperature and temperature ratios
real(DP) :: te      ! Te
real(DP) :: th      ! Th (hydrogenic ion temperature)
real(DP) :: tz      ! Tz (impurity temperature)
real(DP) :: tau     ! Te/Ti
real(DP) :: tau_inv ! Ti/Te
real(DP) :: tauh    ! Ti/Te
real(DP) :: tauz    ! Tz/Te

real(DP) :: ztauh   ! Z_i Ti/Te
real(DP) :: ztauz   ! Z_z Tz/Te

!Electron density and density ratios
real(DP) :: ne      ! Electron density n_e
real(DP) :: nh      ! Main hydrogenic ion density n_h
real(DP) :: nz      ! Impurity ion density n_z
real(DP) :: ns      ! Fast particle density n_s
real(DP) :: znz_ne  ! Z_imp nz/ne
real(DP) :: zns_ne  ! Z_(fast ion) ns/ne
real(DP) :: ne_nh   ! ne/nh
real(DP) :: nh_ne   ! ni/ne (fraction of main ions)
real(DP) :: nz_ne   ! nz/ne (fraction of impurity ions)
real(DP) :: ns_ne   ! ns/ne (fraction of fast ions)
real(DP) :: fte     ! n_e (trapped) / n_e

!Fractions
real(DP) :: fft
real(DP) :: fzft
real(DP) :: ftrt
real(DP) :: gm
real(DP) :: bta
real(DP) :: xt

!Collision frequency related
real(DP) :: vei     ! Electron-ion collision rate
real(DP) :: vef     ! nu effective
real(DP) :: bt
real(DP) :: bt1

!Geometry and equilibrium quantities
real(DP) :: rhos    ! Gyroradius
real(DP) :: rmaj    ! Local major radius
real(DP) :: rmaj0   ! Major radius at magnetic axis
real(DP) :: rmin    ! Local minor radius
real(DP) :: amin    ! Plasma minor radius a
real(DP) :: eps0    ! Global aspect ratio a/R_0
real(DP) :: eps     ! Local aspect ratio r/R
real(DP) :: btor    ! Local toroidal magnetic field
real(DP) :: q       ! Magnetic safety factor
real(DP) :: shear   ! Magnetic shear (r/q)(dq/dr)
real(DP) :: shat    ! Effective magnetic shear modified for elongation
real(DP) :: kappa   ! Plasma elongation
real(DP) :: alpha_MHD   ! MHD alpha

!Plasma beta
real(DP) :: beta    ! Plasma beta
real(DP) :: betae   ! Electron beta

!ExB shear rates
real(DP) :: wexb    ! Shear rate w_ExB

!Saved internal variables used in dispersion relation

real(DP) :: str = 7.0_DP/3.0_DP
real(DP) :: ftr = 5.0_DP/3.0_DP   
real(DP) :: tvr = 2.0_DP/3.0_DP
real(DP) :: em2 = 1.0_DP
real(DP) :: em1 = 1.0_DP
real(DP) :: em  = 1.0_DP
real(DP) :: kpc = 1.0_DP

!Plasma velocity normalized by the local ion soundspeed Cs = sqrt ( Te/mi )
real(DP) :: vtor       ! Toroidal velocity
real(DP) :: vpol       ! Poloidal velocity
real(DP) :: vpar       ! Parallel velocity

!Some other useful quantities
real(DP) :: wde        ! 2._DP * kyrho**2 * csound / R
real(DP) :: Csound     ! speed of sound
real(DP) :: Csound0    ! speed of sound at magnetic axis
real(DP) :: diffBohm   ! Bohm diffusion coefficient

public :: w20main
public :: DPTOMSQZ_MMM
!}}}

contains

Subroutine w20main &!{{{
   ( i_print, h_print,                                     &
     z_te,   z_ne,    z_vtor,  z_vpol,  z_vpar,            &
     z_btor, z_amin,  z_rmin,  z_rmaj,  z_eps0,            &
     z_aimp, z_ahyd,  z_zimp,  z_gte,   z_gti,    z_gtz ,  &
     z_gne,  z_gni,   z_gnz,   z_gvt,   z_gvp,    z_gvpar, &
     z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne,  z_fte,   &
     z_q,    z_shear, z_kappa, z_wexb,  z_Csound0,         &
     diffs,  vconv,   gamma,   omega)
!For definition of input parameters look at file w20_README.txt
!or at the variable definitions in w20data.f90

  implicit none

  !Variable declaration
  integer, intent(in) :: i_print, h_print
  
  real(DP), intent(in) :: &
       z_te,   z_ne,    z_vtor,  z_vpol,  z_vpar, &
       z_btor, z_amin,  z_rmin,  z_rmaj,  z_eps0, &
       z_aimp, z_ahyd,  z_zimp,  &
       z_gte,  z_gti,   z_gtz ,  &
       z_gne,  z_gni,   z_gnz,   &
       z_gvt,  z_gvp,   z_gvpar, & 
       z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
       z_q,    z_shear, z_kappa, z_wexb,  z_Csound0

  !Resulting fluxes
  real(DP), intent(out), dimension(6) :: diffs, vconv

  ! Growth rates and frequency are given in this order:
  ! 1 - ii - the most unstable ion mode in the positive-frequency direction
  ! 2 - ie - the most unstable electron mode in the positive-frequency direction
  ! 3 - ei - the most unstable ion mode in the negative-frequency direction
  ! 4 - ee - the most unstable electron mode in the negative-frequency direction
  Real(DP), Intent(Out), Dimension(1:) :: &
     gamma, &! Growth rate
     omega   ! Frequency

  integer :: ierr, iiu, ieu, iimax, iemax, eiu, eeu, eimax, eemax, j

  !Individual fluxes from ion and electron call
  real(DP), dimension(6) :: idiffs, ivconv, ediffs, evconv

  !Eigenvectors resulting from solved system of equations
  real(DP), dimension(neq,neq) :: izvr, izvi, ezvr, ezvi

  !Complex frequency of the mode 
  complex(DP), dimension(neq) :: izz, ezz, zz

  !Geometric average factor multiplying W_de, used to compute momentum thermoelectric pinch
  real(DP) :: G_ave_i, G_ave_e

  !Parallel wavenumber of drift mode, used to compute momentum pinches
  real(DP) :: kps_i, kps_e

  !Logical variable decides if separate correlation length needed
  logical :: l_elc

  !Initialize all internal variables
  call w20init( i_print, h_print, &
     z_te,   z_ne,    z_vtor,  z_vpol,  z_vpar, &
     z_btor, z_amin,  z_rmin,  z_rmaj,  z_eps0, &
     z_aimp, z_ahyd,  z_zimp, &
     z_gte,  z_gti,   z_gtz,  &
     z_gne,  z_gni,   z_gnz,  &
     z_gvt,  z_gvp,   z_gvpar,&
     z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
     z_q,    z_shear, z_kappa, z_wexb,  z_Csound0 )


  ! Algorithm starts here

  !Initialize geometric factor to 1.0 (circular)
  G_ave_e = 1._DP
  G_ave_i = 1._DP

  !Search for most unstable ion mode
  searchmode = S_ION

  !Variable correlation length for ion modes

  call w20solv( izz, izvr, izvi, iiu, ieu, iimax, iemax, G_ave_i, kps_i, ierr )

  ! Record the most signifiant modes
  gamma = 0._DP
  omega = 0._DP
  if ( iimax > 0 ) then
     gamma(1) = dimag( izz(iimax) )
     omega(1) = dreal( izz(iimax) )
  end if
  if ( iemax > 0 ) then
     gamma(2) = dimag( izz(iemax) )
     omega(2) = dreal( izz(iemax) )
  end if

  !Search for most unstable electron mode using varying correlation length
  searchmode = S_ELC

  !Different correlation length for electrons
  zflh = fle
  zflz = aimp * fle / zimp**2
  kyrho = sqrt(fle)
  wde  = 2._DP * kyrho * csound / rmaj

  call w20solv( ezz, ezvr, ezvi, eiu, eeu, eimax, eemax, G_ave_e, kps_e, ierr )

  ! Record the most signifiant modes
  if ( eimax > 0 ) then
     gamma(3) = dimag( ezz(eimax) )
     omega(3) = dreal( ezz(eimax) )
  end if
  if ( eemax > 0 ) then
     gamma(4) = dimag( ezz(eemax) )
     omega(4) = dreal( ezz(eemax) )
  end if

  do j=1,4
     if ( omega(j).lt.0._DP ) then ! Ion mode
        gamma(j) = max( 0._DP, gamma(j) - abs(wexb) )
     else
        gamma(j) = max( 0._DP, gamma(j)**2-0.25_DP*wexb**2 )
        if ( gamma(j) > 0._DP ) gamma(j) = sqrt(gamma(j))
     end if
  end do


  !Compute transport coefficients

  diffs = 0._DP
  vconv = 0._DP

  idiffs = 0._DP
  ediffs = 0._DP

  ivconv = 0._DP
  evconv = 0._DP

  !No unstable modes, then no transport
  if ( iiu .eq. 0 .and. ieu .eq. 0 .and. eiu .eq. 0 .and. eeu .eq. 0 ) then
     return
  end if

  !Consider whether a separate correlation length for electrons is needed
  !Conditions for yes:
  ! Fastest growing mode in electron w20solv call gives more transport than
  ! fastest growing mode in ion w20solv call
  ! (note that unstable modes must exist for the comparison to occur)
  ! No electron modes in ion w20solv call
  !Conditions for no
  ! All others
  l_elc = .false.
  if ( eeu > 0 ) then

     if ( ieu > 0 ) then
        if ( dimag(izz(iemax))/sqrt(flh) .lt. dimag(ezz(eemax))/sqrt(fle) ) then
           l_elc=.true.
        else
           l_elc=.false.
        end if
     else
        l_elc=.true.
     end if

  else if ( eeu .eq. 0 ) then
     l_elc=.false.
  end if

  if ( l_elc ) then
    
     !Set up call to flux computation using only
     ! -Eigenvector associated with fastest growing electron mode in electron w20solv call
     ! -Eigenvalue associated with fastest growing electron mode in electron w20solv call
     zz = ezz
     ezz= (0._DP,0._DP)
     ezz(1) = zz(eemax)
     
     ezvr(1:neq,1) = ezvr(1:neq,eemax) 
     ezvi(1:neq,1) = ezvi(1:neq,eemax) 

     !Use different correlation length for electron modes
     zflh = fle
     zflz = aimp * fle / zimp**2
     kyrho = sqrt(fle)
     wde  = 2._DP * kyrho * csound / rmaj

     !Compute fluxes
     call w20diff( ezz, ezvr, ezvi, G_ave_i, kps_e, ediffs, evconv )

     !Set up call to flux computation with ion modes
     ! -Remove fastest growing electron mode in ion w20solv call
     if (iemax>0) izz(iemax)=(0._DP,0._DP)

  end if

  !Use correlation length for ion modes
  kyrho = (dsqrt(flh)-0.1_DP*fte)*epsncorr
  zflh = kyrho**2._DP / zhyd**2._DP
  zflz = aimp * kyrho**2._DP / zimp**2._DP

  wde  = 2._DP * kyrho * csound / rmaj

  !Compute fluxes
  call w20diff( izz, izvr, izvi, G_ave_i, kps_i, idiffs, ivconv )

  !Add transport from two calls
  diffs = idiffs + ediffs
  vconv = ivconv + evconv

End Subroutine w20main!}}}

    subroutine w20init &!{{{
      (  i_print, h_print, &
         z_te,   z_ne,    z_vtor, z_vpol,   z_vpar, &
         z_btor, z_amin,  z_rmin, z_rmaj,   z_eps0, &
         z_aimp, z_ahyd,  z_zimp, &
         z_gte,  z_gti,   z_gtz , &
         z_gne,  z_gni,   z_gnz, &
         z_gvt,  z_gvp,   z_gvpar, &
         z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
         z_q,    z_shear, z_kappa, z_wexb , z_Csound0 )

      implicit none

      real(DP), intent(in) :: &
           z_te,   z_ne,   z_vtor, z_vpol, z_vpar,  &
           z_btor, z_amin, z_rmin, z_rmaj, z_eps0,  &
           z_aimp, z_ahyd, z_zimp, &
           z_gte,  z_gti,  z_gtz , &
           z_gne,  z_gni,  z_gnz,  &
           z_gvt,  z_gvp,  z_gvpar, &
           z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
           z_q,    z_shear, z_kappa, z_wexb, z_Csound0

      integer, intent(in) :: i_print, h_print
      
      lprint = i_print
      hprint = h_print
      te     = z_te
      ne     = z_ne*1E-19_DP
      btor   = z_btor
      amin   = z_amin
      rmin   = z_rmin
      rmaj   = z_rmaj
      eps0   = z_eps0
      aimp   = z_aimp
      ahyd   = z_ahyd
      zimp   = z_zimp
      kyrho  = z_kyrho
      tauh   = z_tauh
      tauz   = z_tauz
      nz_ne  = z_nz_ne
      ns_ne  = z_ns_ne
      fte    = z_fte
      q      = z_q
      shear  = z_shear
      kappa  = z_kappa
      Csound0= z_Csound0 
      
      !Initialize some additional variables
      zhyd  = 1._DP
      zfion = 1._DP

      !Aspect ratio and major radius at axis

      rmaj0 = amin / eps0
      eps   = max( 0.01_DP, rmin / rmaj0 )

      !Gradients
      gte    = sign( max( abs(z_gte), grdmin ), z_gte ) * rmaj / rmaj0
      gti    = sign( max( abs(z_gti), grdmin ), z_gti ) * rmaj / rmaj0
      gtz    = sign( max( abs(z_gtz), grdmin ), z_gtz ) * rmaj / rmaj0
      gne    = sign( max( abs(z_gne), grdmin ), z_gne ) * rmaj / rmaj0
      gni    = sign( max( abs(z_gni), grdmin ), z_gni ) * rmaj / rmaj0
      gnz    = sign( max( abs(z_gnz), grdmin ), z_gnz ) * rmaj / rmaj0
      gvp    = sign( max( abs(z_gvp), grdmin ), z_gvp ) * rmaj / rmaj0
      gvt    = sign( max( abs(z_gvt), grdmin ), z_gvt ) * rmaj / rmaj0
      gvpar  = sign( max( abs(z_gvpar), grdmin ), z_gvpar ) * rmaj / rmaj0

      !Effective magnetic shear
      shat = sqrt(max(2._DP*shear-1._DP+(kappa*(shear-1._DP))**2,0._DP))

      !Temperatures and temperature ratios
      th = te*tauh
      tz = te*tauz

      tau      = te  / max( th, 0.001_DP )    ! Te/Ti
      tau_inv  = th  / max( te, 0.001_DP )    ! Te/Ti

      ztauh = tauh / zhyd               ! T_h / ( Z_h Te )
      ztauz = tauz / zimp               ! T_z / ( Z_z Te )

      !Correlation lengths (k_y rho )**2
      epsncorr = min( 2._DP / abs(gne), 4._DP )**0.17_DP
      flh   = (0.7_DP+2.4_DP/(7.14_DP*q*shat+0.1_DP)) &
              *kyrho**2 * 2._DP/(1._DP+tau_inv)
      kyrho = (sqrt(flh)-0.1_DP*fte)*epsncorr
      zflh  = kyrho**2 / zhyd**2
      zflz  = aimp * kyrho**2 / zimp**2    ! A_z * ( k_y rho ) / Z_z**2._DP

      znz_ne = z_zimp*nz_ne
      zns_ne = zfion*ns_ne

      !Densities and density ratios
      nh_ne = (1._DP - znz_ne - zns_ne )
      ne_nh = 1._DP / nh_ne
      nh    = ne * nh_ne
      nz    = ne * nz_ne
      ns    = ne * ns_ne

      !Effective charge Zeff
      zeff = (nh_ne*zhyd**2 + nz_ne*zimp**2 + ns_ne*zfion**2 ) 

      fft   = (1._DP - znz_ne - zns_ne ) / ( 1._DP - fte )
      fzft  = znz_ne / ( 1._DP - fte )
      ftrt  = ( 5._DP/3._DP ) * tauh

      !Collision related quantities
      bt    = 1.5_DP
      bt1   = bt - 2.5_DP

      !Plasma beta
      beta = 4.027E-22_DP * ( z_ne*te*(1.0_DP + nh_ne*tauh + nz_ne*tauz ) ) / (btor**2)

      !Electron beta
      betae = 4.027E-22_DP * z_ne*te / (btor**2)      

      !MHD alpha
      alpha_MHD = em*(q**2.0_DP)*betae*( gne + gte + tauh*(gni + gti) )

      !Speed of sound
      Csound = sqrt( te * 1000._DP * 1.602E-19_DP / ( 1.6725E-27_DP * ahyd ) )

      !Ion gyroradius
      rhos = Csound / ( 1.e8_DP*0.957_DP*btor )

      !Toroidal velocity normalized by speed of sound
      vtor = z_vtor / Csound0

      !Poloidal velocity normalized by speed of sound
      vpol = z_vpol / Csound0

      !Parallel velocity normalized by speed of sound
      vpar = z_vpar / Csound0

      !Magnetic drift frequency with fixed k_y rho
      wde = 2.0_DP * sqrt(flh) * Csound / rmaj

      !ExB shear rates
      wexb   = abs(z_wexb / wde)

      !Bohm diffusion coefficient
      diffBohm = rhos * Csound

      !Electron ion collision rate (main ions only)
      vei = 0.09_DP*1.0E4_DP * nh * ( 15.95_DP - log( sqrt(ne)/te) ) / te**1.5_DP

      !Effective collision rate
      vef = vei / ( eps * wde )

      !Electron correlation length
      fle = max(0.005_DP,0.024_DP-0.0064_DP*(vei/10000.0_DP-0.61_DP)+0.4_DP*betae)

      gm   = 1.0_DP / (1.0_DP-fte)
      bta  = ftrt*(gm+tauh)
      xt   = 1.0_DP / (1.0_DP+tauh)

      lprint = i_print

    end subroutine w20init!}}}

subroutine w20diff &!{{{
   ( zz, zvr, zvi, G_ave, kps, diffs, conv )

  implicit none

  integer, parameter :: ichn = 6

  complex(DP), intent(in), dimension(neq) :: zz
  real(DP), intent(in), dimension(neq,neq) :: zvr, zvi
  real(DP), intent(out), dimension(6) :: diffs, conv
  real(DP), intent(in) :: G_ave, kps

  !Internal variables
  integer :: i,j               !Looping indeces

  !Unstable mode
  complex(DP) :: wu         !Unstable mode complex frequency
  real(DP) :: wu_sq         !Magnitude of wu squared
  real(DP) :: re_wu, im_wu  !Real and imaginary part of unstable mode
  real(DP) :: im_zz         !Imaginary part of unstable mode without ExB shear stabilization
  real(DP) :: im_wu_inv     !Inverse of the growth rate
  real(DP) :: wu_res        !Residual part of omega unstable?
  integer :: iunst             !Number of unstable modes

  !Some physics quantities
  real(DP) :: d1

  !Scaling factor for geometric correction (1/gdro2 in Jetto code)
  real(DP) :: gk = 1._DP

  !Fractions
  real(DP) :: thrd2
  real(DP) :: thrd = 0.3333333333_DP

  !eta_i threshold
  real(DP) :: eith

  !Various quantities related to diffusivities
  real(DP), dimension(ichn) :: xi, xe, xd, xq, xdq, xtv, hpt, dm, dmt, &
       chi, che, d, chq, dhq
  real(DP) :: hp, shp, shpe, dmi, hpe

  !Local diffusivities and accumulators for diffusivities
  real(DP) :: chic, schef, def, dqeff, chqeff, schq, ceft, deft
  real(DP) :: schi, gci, sche, gce, sd, gd, sdm, gcq, sdq, gnq, gei
  real(DP) :: xih, xqh, xeh, xdh, xhh

  real(DP) :: pe, pq, pn, pnq,tiq, pi, piq

  !Possibly, non-diagonal elements and pinches from EM effects
  real(DP) :: dms1, dms2, dms21, dms22, dmi21, dmi22, smp, smpt, bv, dmip, adisp, bdisp, &
       dt, dmit, dms, dmst, dmd, dmdn, dmdt, dmdt2, dmi1, dmi2, dadk, dbdk

  !Variables for output of poloidal, toroidal, perpendicular, thermoelectric
  !components of the pinch velocities
  real(DP)::v_p_tor,v_p_pol,v_p_per,v_p_thr

  !Collision related parameters
  complex(DP) :: ga, gb, gm2
  real(DP) :: re_ga, im_ga, re_gb, im_gb
  real(DP) :: hr, xde, xdi, yda

  !Related to momentum transport
  real(DP) :: dmh, dmef, dmef1, dmef2, dmef21, dmef22, dmdiag, dmdiagt, dmeftest
  real(DP) :: dmeff, dmtef, dmteff
  real(DP) :: sp, smef1, smef21, smef22
  real(DP) :: tsour
  real(DP) :: kxi
  real(DP) :: ainf, imf, tc, dp1, gp1, gp2
  real(DP) :: elfs

  !Lengthscales related to momentum transport not yet purified
  real(DP) :: lvf, lvft

  !Variables related to the parallel wavenumber in momentum transoprt computation
  real(DP) :: Kkap, Hm, kap1, acn, RDexp, ks

  !Magnitude square, real and imaginary parts of various perturbations
  real(DP) :: nehat_sq, re_nehat, im_nehat, nzhat_sq, re_nzhat, im_nzhat
  real(DP) :: re_eln, im_eln
  real(DP) :: dni, dne, dnz, di, de, dz, dnin
  real(DP) :: dznt, dzn1, dzn2
  real(DP) :: dtimp
  real(DP) :: dqn, dqt, dtiq
  real(DP) :: phs

  !More unknown variables
  !Random pieces of equations which need variables
  real(DP) :: tt, et, ht, h, k, t, ts, t1, t2, a3
  real(DP) :: c1, c2, e1, e2, h1, h2, at, ct, a2, a1, a, b, c, e, dh, f
  real(DP) :: stf, eqh

  !Divergences
  real(DP) :: divgat, divat, devgat, devat, diva3, devga1, devga2, devga3, deva1, deva2, deva3
  real(DP) :: divga1, divga2, divga3, diva1, diva2, deva, devb
  real(DP) :: divga, divgb, diva, divb, devga, devgb

  !Sources
  real(DP) :: svat, sva1, sva2, sva3, sva, svb

  !Effective collisionality normalized by electron perturbation magnitude square
  real(DP) :: vefn

  !Complex variables from Weiland's code
  complex(DP) :: Fm,hf,Vg,WRES
  complex(DP) :: FIH,NEF,NRAT,ELN,TEF,FRAT,AV,NIF,MRT
  complex(DP) :: TII,TICF,CN,AVRAT,MRT1,MRT2,NII
  complex(DP) :: nehat,BT2
  complex(DP) :: IU,HC
  complex(DP) :: DMSP,KPF,kpf1,kpf2,ELMS

  wde = 2._DP * csound * kyrho / rmaj0

  !Lengthscales of velocity as defined by J. Weiland
  !Note: these are replaced by gradients (gvt, gvp) in most places

  lvf   = rmaj / ( amin * gvp )
  lvft  = rmaj / ( amin * gvt )

  !Initialize some variables
  diffs = 0._DP
  iunst = 0
  wu    = ( 0._DP, 0._DP ) 

  THRD2=THRD*THRD
  IU = (0._DP,1._DP)
  SHP = 0._DP
  CHIC = 0._DP
  SHPE = 0._DP
  SCHEF = 0._DP
  DEF = 0._DP
  EQH = (gtz/gnz)-STR+FTR*(2._DP/gnz)
  d1 = 6.616_DP*dsqrt(ahyd)/(rmaj*btor**2)

  !Two missing variables
  RDexp = 1._DP

  !Initialize diffusion matrix elements 
  chi(1:6)=0._DP ; che(1:6)=0._DP ; d(1:6)  =0._DP ; chq(1:6)=0._DP 
  dm(1:6) =0._DP ; dhq(1:6)=0._DP ; dmt(1:6)=0._DP ; hpt(1:6)=0._DP

  !Initialize effective diffusivity
  SCHI=0._DP ; SCHE=0._DP ; SD  =0._DP ; SDM =0._DP ; SDQ =0._DP ; XHH =0._DP

  !Main loop over unstable modes
  !(further down is second loop for momentum transport)
  !Thermal and particle diffusivities are computed here

  main_modes: do j=1,neq

     re_wu=DREAL(zz(j))
     im_wu=DIMAG(zz(j))

     !Cycle when growth rate is negative (stable mode)
     if (im_wu.le.0._DP) cycle main_modes

     !Notice different substraction of ExB shear for
     !ion modes and electron modes
     if (re_wu.lt.0._DP) then
        im_wu = im_wu - dabs(wexb)
     else if ( dreal(zz(j)).gt.0._DP) then
        wu_res = im_wu**2 - 0.25_DP *wexb**2
        if (wu_res .lt. 0._DP) cycle
        im_wu = sqrt(wu_res)
     end if

     wu=re_wu+(0._DP,1._DP)*im_wu

     if(im_wu .le. 1E-3_DP) then
        cycle main_modes
     end if

     !Define complex frequency of mode after ExB shear effect
     wu=CMPLX(re_wu,im_wu)
     WRES=wu+FTRT
     gm2=1._DP+0.5_DP*gte/(wu-1._DP+IU*vef)
     BT2=BT-2.5_DP*GM2
     HC=wu-FTR+TVR*BT1
     nehat=wu*wu-2._DP*FTR*wu+FTR+IU*VEF*HC
     wu_sq=re_wu*re_wu+im_wu*im_wu
     im_wu_inv=1._DP/im_wu
     !
     !   ******   COLLISION PARAMETERS  *******
     !
     GA=wu-FTR+TVR*BT2
     GB=(wu-FTR)/(wu-1._DP+IU*VEF)
     re_ga=DREAL(GA)
     im_ga=DIMAG(GA)
     re_gb=DREAL(GB)
     im_gb=DIMAG(GB)

     HR=re_wu-FTR+TVR*BT1
     XDE=wu_sq-FTR*re_wu
     XDI=wu_sq+FTR*tau_inv*re_wu
     YDA=re_wu*(1._DP-(2._DP/gne))+(gte/gne)-STR
     !   ***************************************
     !

     re_nehat=DREAL(nehat)
     im_nehat=DIMAG(nehat)
     nehat_sq=(re_nehat)**2+(im_nehat)**2
     !
     !      IF(VEF.EQ.0.) GO TO 25
     DIVGA=FTRT*((2._DP/gne)*(re_ga*im_nehat-im_ga*re_nehat) &
          -YDA*im_wu+im_wu*HR*(1._DP-(2._DP/gne)))

     DIVGB=FTRT*(re_gb*im_nehat-im_gb*re_nehat-im_wu)

     DIVA=XDI*((2._DP/gne)*(re_ga*re_nehat+im_ga*im_nehat)-im_wu*im_wu*(1._DP-(2._DP/gne))-YDA*HR)

     DIVB=XDI*(re_gb*re_nehat+im_gb*im_nehat-HR)

     DEVGA=BT1*((YDA-VEF*im_ga*(2._DP/gne))*im_nehat    &
          -(im_wu*(1._DP-(2._DP/gne))+VEF*re_ga*(2._DP/gne)) &
          *re_nehat) &
          +FTR*(im_wu*YDA+(2._DP/gne)*(im_ga*re_nehat-re_ga*im_nehat)-im_wu*HR*(1._DP-(2._DP/gne)))

     DEVGB=BT1*((1._DP-VEF*im_gb)*im_nehat-VEF*re_gb*re_nehat) &
          +FTR*(im_wu+im_gb*re_nehat-re_gb*im_nehat)

     DEVA=XDE*((2._DP/gne)*(re_ga*re_nehat+im_ga*im_nehat)-im_wu*im_wu*(1._DP-(2._DP/gne))-YDA*HR) &
          -BT1*(re_wu-FTR)*((YDA-VEF*im_ga*(2._DP/gne))*re_nehat+(im_wu*(1._DP-(2._DP/gne)) &
          +VEF*re_ga*(2._DP/gne))*im_nehat)

     DEVB=XDE*(re_gb*re_nehat+im_gb*im_nehat-HR)-BT1*(re_wu-FTR)*((1._DP- &
          VEF*im_gb)*re_nehat+VEF*re_gb*im_nehat)

     SVA=(2._DP/gne)*(re_ga*re_nehat+im_ga*im_nehat)-im_wu*im_wu*(1._DP-(2._DP/gne))-YDA*HR

     SVB=re_gb*re_nehat+im_gb*im_nehat-HR

     VEFN=VEF/nehat_sq

     DIVGA=DIVGA*VEFN
     DIVGB=DIVGB*VEFN
     DEVGA=DEVGA*VEFN
     DEVGB=DEVGB*VEFN
     DIVA=DIVA*VEFN
     DIVB=DIVB*VEFN
     DEVA=DEVA*VEFN
     DEVB=DEVB*VEFN
     SVA=SVA*VEFN
     SVB=SVB*VEFN
     !   25 CONTINUE
     !
     STF=STR-FTR*(2._DP/gne)
     A1=wu_sq*((2._DP/gne)-1._DP)+2._DP*re_wu*STF
     A=wu_sq*(A1+FTR*(STR*(2._DP/gne)-11._DP*THRD-tau_inv*(1._DP- &
          FTR*(2._DP/gne))))+FTR*FTR*tau_inv*(2._DP*re_wu*(1._DP-(2._DP/gne))-STF)
     A=A/nehat_sq
     B=(wu_sq*(2._DP*(re_wu-FTR)+FTR*tau_inv)- &
          FTR*FTR*tau_inv)/nehat_sq
     C=wu_sq*(A1+TVR*FTR*((2._DP/gne)-4._DP))+ &
          FTR*FTR*(2._DP*re_wu*((2._DP/gne)-1._DP)+STF)
     C=C/nehat_sq
     DH=(wu_sq*(2._DP*re_wu-5._DP)+FTR*FTR)/nehat_sq
     E=wu_sq*(1._DP-(2._DP/gne))-2._DP*re_wu*STF+FTR*(11._DP*THRD &
          -STR*(2._DP/gne))
     E=E/nehat_sq
     F=2._DP*(-re_wu+FTR)/nehat_sq
     !
     re_nzhat=re_wu**2-im_wu**2+2._DP*FTR*ztauz*re_wu+FTR*ztauz*ztauz
     im_nzhat=2._DP*im_wu*(re_wu+FTR*ztauz)
     nzhat_sq=re_nzhat**2+im_nzhat**2
     !
     !   ****   SPLITTING IN  EN  AND EE  ************
     !
     !
     !  Ion thermal conductivity
     !
     DIVGA1=FTRT*((STR-re_wu)*im_wu+im_wu*HR)*VEFN
     DIVGA2=FTRT*(re_ga*im_nehat-im_ga*re_nehat+re_wu*im_wu-im_wu*HR)*VEFN
     DIVGA3=-FTRT*im_wu*VEFN

     DIVA1=XDI*(-im_wu*im_wu+(STR-re_wu)*HR)*VEFN
     DIVA2=XDI*(re_ga*re_nehat+im_ga*im_nehat+im_wu*im_wu+re_wu*HR)*VEFN
     DIVA3=-XDI*HR*VEFN

     !
     !   Electron thermal conductivity
     !
     DEVGA1=(BT1*((re_wu-STR)*im_nehat-im_wu*re_nehat)+FTR*(im_wu*(re_wu-STR)-im_wu*HR))&
          *VEFN
     DEVGA2=(BT1*((-re_wu-VEF*im_ga)*im_nehat+(im_wu-VEF*re_ga)*re_nehat)&
          +FTR*(-im_wu*re_wu+im_ga*re_nehat-re_ga*im_nehat+im_wu*HR))*VEFN
     DEVGA3=(BT1*im_nehat+FTR*im_wu)*VEFN

     DEVA1=(XDE*(-im_wu*im_wu-(re_wu-STR)*HR)-BT1*(re_wu-FTR)*((re_wu-STR)*re_nehat&
          +im_wu*im_nehat))*VEFN
     DEVA2=(XDE*(re_ga*re_nehat+im_ga*im_nehat+im_wu*im_wu+re_wu*HR)-BT1*(re_wu-FTR)&
          *((-re_wu-VEF*im_ga)*re_nehat-(im_wu-VEF*re_ga)*im_nehat))*VEFN
     DEVA3=(-XDE*HR-BT1*(re_wu-FTR)*re_nehat)*VEFN

     SVA1=(-im_wu*im_wu-(re_wu-STR)*HR)*VEFN
     SVA2=(re_ga*re_nehat+im_ga*im_nehat+im_wu*im_wu+re_wu*HR)*VEFN
     SVA3=-HR*VEFN

     !------------------------------------------------------------------------ 
     ! Contributions to the ion conductivity 
     !------------------------------------------------------------------------ 

     A2=(wu_sq*wu_sq+FTR*((STR+FTR*tau_inv-2._DP*re_wu)*wu_sq+FTR*tau_inv*(FTR-&
          2._DP*re_wu)))/nehat_sq
     A3=(wu_sq*(-wu_sq+2._DP*STR*re_wu-FTR*(11._DP*THRD+tau_inv))+FTR*FTR&
          *tau_inv*(2._DP*re_wu-STR))/nehat_sq

     !------------------------------------------------------------------------ 
     ! Contributions to the electron conductivity 
     !------------------------------------------------------------------------ 

     C1=(nehat_sq-wu_sq*wu_sq+14._DP*THRD*wu_sq*re_wu-40._DP*THRD2*wu_sq&
          -50._DP*THRD2*re_wu+175._DP/27._DP)/nehat_sq
     C2=(wu_sq*wu_sq-10._DP*THRD*wu_sq*re_wu+10._DP*THRD2*wu_sq&
          +50._DP*THRD2*re_wu-125._DP/27._DP)/nehat_sq

     !------------------------------------------------------------------------ 
     !   Contributions to the electron diffusivity 
     !------------------------------------------------------------------------ 

     E1=(wu_sq+11._DP*THRD*FTR-2._DP*re_wu*STR)/nehat_sq
     E2=-(-2._DP*re_wu*FTR+(wu_sq+STR*FTR))/nehat_sq

     !------------------------------------------------------------------------ 
     ! Contributions to the main ion Conductivity 
     !------------------------------------------------------------------------ 

     H1=(wu_sq*(-wu_sq-2._DP*ztauz*re_wu*STR&
          +FTR*ztauz*ztauz*(-11._DP*THRD)+FTR*tau_inv*ztauz)&
          +FTR*FTR*tau_inv*ztauz*ztauz*(2._DP*re_wu+STR*ztauz))/nzhat_sq

     H2=(wu_sq*(wu_sq+2._DP*ZTAUZ*re_wu*FTR+FTR*ZTAUZ*ZTAUZ*STR-FTR*tau_inv*ZTAUZ*FTR)&
          -FTR*FTR*tau_inv*ZTAUZ*ZTAUZ*(2._DP*re_wu+FTR*ZTAUZ))/nzhat_sq

     !------------------------------------------------------------------------ 
     ! Contributions to the impurity conductivity 
     !------------------------------------------------------------------------ 

     T1=(wu_sq*(-wu_sq-2._DP*ZTAUZ*re_wu*STR-FTR*ZTAUZ*ZTAUZ*8._DP*THRD)&
          +FTR*FTR*ZTAUZ**3*(2._DP*re_wu+ZTAUZ*STR))/nzhat_sq
     T2=(wu_sq*(wu_sq+2._DP*ZTAUZ*re_wu*FTR+FTR*ZTAUZ*ZTAUZ*TVR)&
          -FTR*FTR*ZTAUZ**3*(2._DP*re_wu+ZTAUZ*FTR))/nzhat_sq

     !------------------------------------------------------------------------ 
     ! Contributions to the impurity diffusivity 
     !------------------------------------------------------------------------ 

     DZN1=(-re_nzhat+2._DP*(re_wu+ZTAUZ*STR)*(re_wu+FTR*ZTAUZ))/nzhat_sq
     DZN2=(re_nzhat-2._DP*(re_wu+ZTAUZ*FTR)*(re_wu+FTR*ZTAUZ))/nzhat_sq
     !
     DIVGAT=DIVGA1+(2._DP/gne)*DIVGA2+(gte/gne)*DIVGA3
     DIVAT=DIVA1+(2._DP/gne)*DIVA2+(gte/gne)*DIVA3
     DEVGAT=DEVGA1+(2._DP/gne)*DEVGA2+(gte/gne)*DEVGA3
     DEVAT=DEVA1+(2._DP/gne)*DEVA2+(gte/gne)*DEVA3
     SVAT=SVA1+(2._DP/gne)*SVA2+(gte/gne)*SVA3
     AT=(2._DP/gne)*A2+A3
     CT=C1-1.+(2._DP/gne)*C2
     ET=E1+(2._DP/gne)*E2
     HT=H1+(2._DP/gnz)*H2
     TT=T1+(2._DP/gnz)*T2
     DZNT=DZN1+(2._DP/gne)*DZN2
     !
     ! **** IMPURITIES *****
     !
     nzhat_sq=(re_wu*(re_wu+2._DP*FTR*ZTAUZ)-im_wu*im_wu+FTR*ZTAUZ*ZTAUZ)**2 &
          +4._DP*im_wu*im_wu*(re_wu+FTR*ZTAUZ)**2
     DTIMP=(wu_sq*(wu_sq*((2._DP/gnz)-1._DP)+2._DP*ZTAUZ*re_wu*EQH+FTR*ZTAUZ &
          *ZTAUZ*(2._DP*(gtz/gnz)-11./3.+STR*(2._DP/gnz))+FTR*tau_inv*ZTAUZ*(1. &
          +(gtz/gnz)-FTR*(2._DP/gnz)))+FTR*tau_inv*(2._DP*FTR*ZTAUZ*ZTAUZ*re_wu*(1._DP-(2._DP/gnz)) &
          -FTR*ZTAUZ**3*EQH))/nzhat_sq
     !  *************
     !
     DNI=(re_wu+FTR*tau_inv)**2+im_wu*im_wu
     DNE=(re_wu-FTR)**2+im_wu*im_wu
     DNZ=(re_wu+FTR*ZTAUZ)**2+im_wu*im_wu
     DI=(gti/gni)*kyrho*DNI
     DE=(gte/gne)*kyrho*DNE
     DZ=(gtz/gnz)*kyrho*DNZ
     !
     ! **** IMPURITIES *****
     !
     nzhat_sq=re_nzhat**2+im_nzhat**2
     !
     H=(wu_sq*(wu_sq*((2._DP/gnz)-1._DP)-2._DP*ZTAUZ*re_wu*(STR-FTR*(2._DP/gnz)) &
          +FTR*ZTAUZ*ZTAUZ*(-11.*THRD+STR*(2._DP/gnz))+FTR*tau_inv*ZTAUZ*(1.-FTR &
          *(2._DP/gnz)))+FTR*FTR*tau_inv*ZTAUZ*ZTAUZ*(2._DP*re_wu*(1._DP-(2._DP/gnz))+(STR-FTR &
          *(2._DP/gnz))*ZTAUZ))/nzhat_sq
     K=ZTAUZ*(FTR*FTR*tau_inv*ZTAUZ*ZTAUZ-wu_sq*(2._DP*re_wu &
          +FTR*(2._DP*ZTAUZ+tau_inv)))/nzhat_sq
     !
     !  *************
     !
     T=(wu_sq*(wu_sq*((2._DP/gnz)-1._DP)-2._DP*ZTAUZ*re_wu*(STR-FTR*(2._DP/gnz)) &
          +FTR*ZTAUZ*ZTAUZ*(-8.*THRD+TVR*(2._DP/gnz)))+FTR*FTR*(ZTAUZ)**3 &
          *(2._DP*re_wu*(1.-(2._DP/gnz))+ZTAUZ*(STR-FTR*(2._DP/gnz))))/nzhat_sq
     TS=ZTAUZ*(FTR*FTR*(ZTAUZ)**3-wu_sq*(2._DP*re_wu+5._DP*ZTAUZ))/nzhat_sq
     !
     T=(GNZ/GNE)*T
     TS=(GNZ/GNE)*TS
     DQN=((2._DP/gnz)-1._DP)*re_nzhat+2._DP*((1._DP-(2._DP/gnz))*re_wu+ZTAUZ*(STR-FTR*(2._DP/gnz)) &
          )*(re_wu+FTR*ZTAUZ)
     DQT=2._DP*ZTAUZ*(re_wu+FTR*ZTAUZ)
     DQN=DQN/nzhat_sq
     DQT=DQT/nzhat_sq
     !
     DTIQ=H-(gtz/gnz)*K
     !      DTQQ=T-EQ*TS
     !
     PHS=(re_wu-FTR)*DREAL(BT2)+im_wu*DIMAG(BT2)
     XDH=D1*(TE**1.5_DP)*im_wu**3/kyrho
     XIH=XDH/DNI
     XEH=XDH/DNE
     XQH=XDH/DNZ

     !----------------------
     ! Ion transport channel
     !----------------------
     ! Diffusive part
     !
     XI(1)=XIH
     XI(2)=tvr*fte*xih*ne_nh*(B-DIVGA3-DIVGB-im_wu_inv*(DIVA3+DIVB))
     XI(3)=-tvr*xih*(gni/gne+fte*ne_nh*(A3+DIVGA1+DIVA1*im_wu_inv))
     !      XI(4)=-XIH*TVR*nz_ne*Z*ne_nh*K*TAUZ*tau_inv
     XI(4)=0._DP
     XI(5)=tvr*xih*zimp*h1/nh
     XI(6)=0._DP
     !
     ! Convective part
     !
     PI=fte*TVR*XIH*ne_nh*(A2+DIVGA2+im_wu_inv*DIVA2)*amin/rmaj0
     PIQ=-TVR*XIH*zimp*H2*nz_ne*ne_nh*amin/rmaj0
     HP=XIH*ne_nh*tau_inv*TVR*FTR*(1._DP-fte)*amin/rmaj0
     !
     ! Effective diffusivity
     !
     GCI = XI(1) &
          + XI(2)*gte/gti &
          + XI(3)*gne/gti &
          + XI(5)*nz*gne/gti &
          - (HP+PI+PIQ)*2._DP*gk*rmaj/amin/gti

     !-----------------------------------
     ! Electron thermal transport channel
     !-----------------------------------
     ! Diffusive part
     !
     XE(1)=0._DP
     XE(2)=fte*XEH*(1._DP+TVR*(DH-DEVGB-DEVGA3-im_wu_inv*(DEVB+DEVA3)))
     XE(3)=-TVR*fte*XEH*(C1+DEVGA1+im_wu_inv*DEVA1)
     XE(4)=0._DP
     XE(5)=0._DP
     XE(6)=0._DP
     !
     ! Convective part
     !
     PE=fte*TVR*XEH*(C2+DEVGA2+im_wu_inv*(DEVA2+VEF*PHS))*amin/rmaj0
     !
     ! Effective diffusivity
     !
     GCE = XE(2) &
          + XE(3)*gne/gte &
          - PE*gk*2._DP*rmaj/amin/gte

     !------------------------------------
     ! Electron particle transport channel
     !------------------------------------
     ! Diffusive part
     !
     XD(1)=0._DP
     XD(2)=-fte*XDH*(ne/te)*(F+im_wu_inv*(SVB+SVA3))
     XD(3)=fte*XDH*(E1-im_wu_inv*SVA1)
     XD(4)=0._DP
     XD(5)=0._DP
     XD(6)=0._DP
     !
     ! Convective part
     !
     PN=-fte*XDH*(E2-im_wu_inv*SVA2)*amin/rmaj0
     !
     ! Effective diffusivity
     !
     GD = XD(2)*te/ne*gte/gne &
          + XD(3) &
          -  PN*gk*2._DP*rmaj/amin/gne

     !-----------------------------------
     ! Impurity thermal transport channel
     !-----------------------------------
     ! Diffusive part
     !
     XQ(1)=0._DP
     XQ(2)=0._DP
     XQ(3)=0._DP
     XQ(4)=1._DP
     XQ(5)=0._DP
     XQ(6)=0._DP
     !
     ! Convective part
     !
     PQ=TVR*XQH*T2*2/rmaj
     !
     ! Effective diffusivity
     !
     !      GCQ=XQ(4)+XQ(5)*NQ/(TZ*(gtz/gnz))-PQ*GK*2._DP*LTZ

     GCQ = -XDH*DQT*NZ/TZ !XDQ(4)
     !------------------------------------
     ! Impurity particle transport channel
     !------------------------------------
     ! Diffusive part
     !
     XDQ(1)=0._DP
     XDQ(2)=0._DP
     XDQ(3)=0._DP
     XDQ(4)=-XDH*DQT*NZ/TZ
     XDQ(5)=XDH*DZN1
     XDQ(6)=0._DP
     !
     ! Convective part
     !
     PNQ=-XDH*DZN2*2/rmaj
     !
     ! Effective diffusivity
     !
     GNQ = XDQ(4)*TZ*(gtz/gnz)/NZ &
          + XDQ(5) &
          - PNQ*GK*rmaj/gnz

     XTV(1)=0._DP
     XTV(2)=0._DP
     XTV(3)=0._DP
     XTV(4)=0._DP
     XTV(5)=0._DP
     XTV(6)=1._DP
     !
     TIQ=1._DP/(TAUZ*(GNI/GNZ))

     SHP=SHP+HP
     !
     HPT(1)=HPT(1)+HP+PI+PIQ
     HPT(2)=HPT(2)+PE
     HPT(3)=HPT(3)+PN
     !      HPT(4)=HPT(4)+PQ
     HPT(5)=HPT(5)+PNQ
     !
     SCHI=SCHI+GCI
     SCHE=SCHE+GCE
     SD=SD+GD
     SDM=SDM+GCQ
     SDQ=SDQ+GNQ
     !
     !
     !      CHQEFF=D1*TE**1.5_DP*WI**3*((gtz/gnz)-TVR-TVR*DTQQ)/DQ
     chqeff=0.0_DP
     schq  =0.0_DP
     DQEFF=XDH*(DQN-(gtz/gnz)*DQT)

     XHH=XHH+XIH

     Do i=1,ichn
        CHI(I)=CHI(I)+XI(I)
        CHE(I)=CHE(I)+XE(I)
        D(I)=D(I)+XD(I)
        !      CHQ(I)=CHQ(I)+XQ(I)
        DHQ(I)=DHQ(I)+XDQ(I)
     end do

  End Do main_modes
  ! 
  !  IF electromagnetic effects or collisions on free electrons are included
  !  the transport coefficients are corrected for this.
  !------------------------------------------------------------
  DT=D1*TE**1.5_DP
  SHPE=0._DP
  SCHEF=0._DP
  DEF=0._DP
  DMI=0._DP
  TSOUR=0._DP
  DMIT=0._DP
  DMS=0._DP
  DMST=0._DP
  DMSP=(0._DP,0._DP)
  DMD=0._DP
  DMDT=0._DP
  DMDT2=0._DP
  DMI1=0._DP
  DMI2=0._DP
  DMS1=0._DP
  DMS2=0._DP
  DMI21=0._DP
  DMI22=0._DP
  SMP=0._DP
  SMPT=0._DP
  BV=0._DP
  DMIP=0._DP    !!  Used for summing up curvature pinch terms
  v_p_pol=0.0_DP
  v_p_tor=0.0_DP
  v_p_per=0.0_DP
  v_p_thr=0.0_DP
  adisp=((2._DP/gne)-1._DP-flh*(FTRT*(2._DP/gne)-tau_inv*(1._DP+(gti/gni))))/(1._DP+flh)
  bdisp=tau_inv*(2._DP/gne)/(1._DP+fte)
  dadk=-2._DP*kyrho*rhos*(FTRT*(2._DP/gne)-tau_inv*(1._DP+(gti/gni))+((2._DP/gne)-1._DP)/(1._DP+flh))
  dadk=dadk/(1._DP+flh)
  dbdk=-2._DP*(2._DP/gne)*tau_inv*kyrho*rhos/(1._DP+flh)**2
  EITH=TVR+10._DP*(2._DP/gne)*tau_inv/9._DP
  KXI=2._DP*q*kyrho/(eps*rmaj0)
  !------------------------------------------- loop

  main_momentum: DO J=1,NEQ

     re_wu=DREAL(zz(j))
     im_wu=DIMAG(zz(j))

     if ( im_wu.lt.0.0_DP ) cycle

     if ( re_wu.lt.0.0_DP ) then !Ion mode
        im_wu = im_wu - dabs(wexb)
     else if ( re_wu .gt.0.0_DP ) then
        wu_res =im_wu**2-0.25_DP*wexb**2
        if (wu_res.lt.0.0_DP) cycle
        im_wu = dsqrt(wu_res)
     end if

     wu=re_wu+(0._DP,1._DP)*im_wu

     if(im_wu .le. 1E-2_DP) then
        cycle
     end if

     wu=CMPLX(re_wu,im_wu)

     im_zz = max( 0.01_DP, dimag(zz(j)) )

     ! --- contr. to chii from em free electr. ----
     !
     DNI=(re_wu+FTR*tau_inv)**2+im_wu*im_wu
     DNIN=(re_wu+FTR*tau_inv)**2+im_zz*im_zz
     wu_sq=re_wu*re_wu+im_wu*im_wu
     XDH=DT*im_wu**3/kyrho
     XIH=XDH/DNI
     FIH=DCMPLX(ZVR(1,J),ZVI(1,J))
     ! F. Halpern: place a lower bound on the absolute value of FIH 
     IF ( abs(fih) .lt. 0.001 ) then
        fih = (0.001,0._DP)
     end if

     !-- Here NEF for disp 9 is defined --
     AV=DCMPLX(ZVR(8,J),ZVI(8,J))
     NEF=FIH-(ZZ(J)-0.5_DP*gne)*AV/KPC
     !--------------------------------------------

     IF(ABS(FIH).LT.0.0001_DP) FIH=(0.0001_DP,0._DP)
     NRAT=NEF/FIH
     ELN=NRAT-1._DP
     re_eln=DREAL(ELN)
     im_eln=DIMAG(ELN)
     AINF=TVR*(FTR*tau_inv*re_eln+im_eln*(im_wu*im_wu+re_wu*(re_wu+FTR*tau_inv))/im_zz)
     HPE=XIH*ne_nh*(1._DP-fte)*eps0*AINF

     SHPE=SHPE+HPE
     ! ****  Free electron heat flux *********************
     !
     !-----------------------------------------------------
     IF(NEQ.EQ.11) GO TO 10099
     IF(NEQ.EQ.9) GO TO 10098
     !--- Here TEF for disp10 is defined ---
     TEF=DCMPLX(ZVR(6,J),ZVI(6,J))
     GO TO 10099
10098 CONTINUE
     !-- Here TEF for disp9 is defined ---
     !      TEF=(gte/gne)*0.5*gne*AV/KPC
     TEF = 0.5_DP*gte*AV/KPC
     !-------------------------------------------------------------
10099 CONTINUE
     FRAT=TEF/FIH
     IMF=-IMAG(FRAT)*DNIN/DNI
     CEFT=(1._DP-fte)*IMF*DT*(2._DP/gte)*im_wu**3/(kyrho*im_zz)
     SCHEF=SCHEF+CEFT
     !**********************************************************
     ! ----Free electron particle flux -----------
     GEI=-IMAG(NRAT)/im_wu
     DEFT=(1._DP-fte)*GEI*(2._DP/gne)*XDH
     DEF=DEF+DEFT
     ! ----Momentum transport ------
     TII=DCMPLX(ZVR(2,J),ZVI(2,J))
     NIF=DCMPLX(ZVR(3,J),ZVI(3,J))
     AVRAT=AV/FIH
     MRT=(TII+NIF)/FIH
     MRT1=0.5_DP*(gti*gne/gni-gne*tvr)/(wu+ftrt*G_ave)
     TICF=TII/FIH
     NII=NIF/FIH
     MRT2=NIF*(1._DP+TVR*wu/(wu+FTRT*G_ave))/FIH
     DMS21=im_wu*im_wu*tau_inv*DREAL(MRT1)
     DMS22=im_wu*im_wu*tau_inv*DREAL(MRT2)
     DMS1=im_wu*im_wu
     DMS2=im_wu*im_wu*tau_inv*DREAL(MRT)
     DMS=DMS1+DMS2
     DMI1=DMI1+DMS1
     DMI2=DMI2+DMS2
     DMI21=DMI21+DMS21    !! part of temp. pert 
     DMI22=DMI22+DMS22    !! density pert including part from temp. pert.
     DMI=DMI+DMS
     DMD=DMD+im_wu**3/wu_sq
     !Comparison between old and new diagonal element
     DMDN=(re_wu+2._DP*tau_inv*G_ave)**2+im_wu*im_wu
     DMDT=DMDT+im_wu**3/DMDN  !! Diagonal element for toroidal momentum transport
     TC=(re_wu+2._DP*tau_inv*G_ave)/DMDN  !! Correl. time for TEP and Termoel. mom pinch (Hahm)
     DP1=2._DP*DT*im_wu*im_wu/kyrho

! By F.Halpern: 18-Feb-2009: Introduce gp1 and gp2 as stated in Jan Weiland's MOMCODE from 25-Jan-2009

! By A Kritz: 31-Mar-2009: Moved line defining ELMS before it is used in computing gp2
     ELMS = (wu+0.5*tau_inv*(gne+gti*gne/gni))*AVRAT/KPC ! From Weiland's momcode

     gp1 = - im_wu * dp1 * vtor * G_ave / dmdn
     gp2 = dreal( IU * ( MRT - tau*em2*ELMS )/( wu + 2._DP*tau_inv*G_ave ) ) * dp1 * vtor * G_ave

!     GP1=0.5*TC*DP1*vtor  !!  TEP  momentum pinch
!     GP2=DREAL(TICF/(wu+2._DP*tau_inv*G_ave))*DP1*vtor  !!  Termoel. momentum pinch


     ELFS = (1._DP+0.5*tau_inv*(gne+gti*gne/gni)/wu)*AVRAT/KPC ! From Weiland's momcode

! By F. Halpern: Add factor of average curvature as per 25-Jan-2009 MOMCODE
     DMSP=im_wu*im_wu*(IU/(wu+2._DP*tau_inv*G_ave))*(1+tau_inv*MRT-em2*ELMS) !! due to v_par incl TC
!     DMSP=im_wu*im_wu*(IU/(wu+2._DP*tau_inv))*(1+tau_inv*MRT-em2*ELMS) !! due to v_par incl TC

     !----- Computation of parallel wavenumber k_\parallel

     !K in Weiland's 2006 EPS paper
     !F. Halpern: replaced toroidal velocity gradient with parallel gradient
     !      Kkap = rmaj0/rmaj*q*kyrho*vtor*gvt
     Kkap = rmaj0/rmaj*q*kyrho*vpar*gvpar

     !\kappa_1 in Weiland's 2006 EPS paper
     !F. Halpern: replaced toroidal velocity gradient with parallel gradient

     kap1 = vpol*gvp/(gne*sign(max(abs(shear),0.1_DP),shear)*kyrho)
     !      kap1 = vpol*gvpar/(gne*sign(max(abs(shear),0.1_DP),shear)*kyrho*dsqrt(kappa))
     kap1 = sign(min(abs(kap1),1.0_DP),kap1)

     !ks in Weiland's Jan-2009 MOMCODE
     ks   = -2._DP*vtor*kyrho*q/(tau_inv*tau_inv)

     Fm=wu*(1.+FTRT)+tau_inv*0.5*(gti*gne/gni-TVR*gne)
     Hm=FTRT*(1.+tau_inv)
     hf=4.0*flh*q**2*wu*0.5*gne

     !Flux average of k_\parallel normalized, plus separate expressions 
     !poloidal and toroidal driving terms
     KPF=-(0.5*(wu+ftr)*(Kkap+ks) + IU*hf*im_wu*KAP1)/((Fm+Hm)*Q*Rmaj0) !!This is dimk_par
     kpf1=-(0.5*(wu+ftr)*(Kkap+ks))/((Fm+Hm)*Q*Rmaj0)
     kpf2=-(IU*hf*im_wu*KAP1)/((Fm+Hm)*Q*Rmaj0)

!     KPF=-(0.5*(wu+ftr)*Kkap+IU*hf*im_wu*KAP1)/((Fm+Hm)*Q*Rmaj0) !!This is dimk_par
!     kpf1=-(0.5*(wu+ftr)*Kkap)/((Fm+Hm)*Q*Rmaj0)

     !Normalization is removed here
     KPF=KPF*CSound/Wde
     kpf1=kpf1*CSound/Wde
     kpf2=kpf2*CSound/Wde

     !Compute the fluxes here
     !Split the different pinches due to velocity shear
     ! v_p_per ---> perpendicular part
     ! v_p_pol ---> poloidal part of parallel velocity pinch
     ! v_p_tor ---> toroidal part of parallel velocity pinch

     v_p_per = v_p_per -2._DP*d1*te**1.5_DP*(eps/q)*dms
     v_p_pol = v_p_pol +2._DP*d1*te**1.5_DP*dreal(DMSP*kpf2)/kyrho
     v_p_tor = v_p_tor +2._DP*d1*te**1.5_DP*dreal(DMSP*kpf1)/kyrho

     !Toroidal momentum pinch term
     DMST=-(eps/q)*dms +dreal(DMSP*kpf)/kyrho
     DMIT=DMIT+DMST

     !---- Summation of curvature pinch term fluxes here
     DMIP=DMIP+GP1+GP2     !!  Summation of curvature pinch terms

     v_p_thr = dmip

     !---- Quantities related to poloidal momentum transport

     !---- Obtained from Jan Weiland's momcode dated from Feb 4 2008
     !---- Not debugged or tested

     Vg=-((wu+ftrt)*dadk+0.5_DP*(gti*gne/gni-EITH*gne)*dbdk) & 
          /(2._DP*(wu+ftrt)+0.5_DP*gne*adisp)
     RDexp=DREAL( 0.25_DP*(gne*gti-eith*gne**2)*Vg/ (wu+ftrt)**2 )

     TSOUR=TSOUR+4.*RDexp*im_wu*im_wu*diffBohm*KXI*KPS/rmaj0**2
     BV=BV+im_wu*im_wu*DREAL(IU*NII*CONJG(NII-TICF))

     CN=wu+FTRT
     ACN=ABS(CN)

     ! End of looping through unstable modes
     ! -------
  End Do main_momentum

  !---- Final computation of momentum fluxes
  !---- First, some unknown computations

  CHIC=-SHPE*GK*2._DP*(rmaj/gti)
  HPT(1)=HPT(1)+EM*SHPE
  CHE(2)=CHE(2)+EM*SCHEF
  D(3)=D(3)+EM1*DEF
  !
  if ( abs(vpol) < 1E-6_DP ) vpol=1E-6_DP
  DMH=-2._DP*D1*TE**1.5_DP*eps0*(Csound/CSOUND0)/vpol
  DMEF=DMH*DMI*LVF
  DMEF1=DMH*DMI1*LVF
  DMEF2=DMH*DMI2*LVF
  DMEF21=DMH*DMI21*LVF  !! splitting of DMEF2
  DMEF22=DMH*DMI22*LVF   !! splitting of DMEF2

  !----- Toroidal momentum transport

  !Diagonal element of momentum diffusivity
  DMDIAGT=D1*TE**1.5_DP*DMDT/kyrho
  DMT(6)=DMDIAGT

  !Momentum pinches expressed as an effective diffusivity
  DMTEF = rmaj / rmaj0 * &
       ( DMIP  &
       + 2._DP*d1*te**1.5_DP*DMIT ) &
       /sign( max(abs(gvt*vtor),0.01_DP) , gvt*vtor )

  !Total momentum pinch expressed as a convection velocity in m/s
  !Diagnostic only, not output in interface
  HPT(6) = ( dmip + 2._DP*d1*te**1.5_DP*dmit ) / sign( max( abs(vtor) ,0.01_DP), vtor )

  !Parts of momentum pinches expressed as a convective velocity
  !These are output through the convective velocity array

  !Sum of thermoelectric and turbulent equipartition pinches
  v_p_thr = DMIP / sign( max(abs(rmaj0*vtor), 0.01_DP), rmaj0*vtor )
!  v_p_thr = DMIP / sign( max(abs(vtor), 0.01), vtor )

  !Sum of parallel velocity and Reynolds stress pinches
  v_p_tor = 2._DP*d1*te**1.5_DP*DMIT / sign( max(abs(rmaj0*vtor), 0.01_DP), rmaj0*vtor )

  !Total flux expressed as an effective diffusivity
  !WARNING: it is better to use diagonal element + convection
  !         Momentum pinches can produce large negative diffusion
  DMTEFF = DMT(6) + DMTEF

  !----- Poloidal momentum diffusivity
  ! WARNING: These computations have not been tested
  !          (There is no poloidal momentum eq. in PTRANSP,
  !           where the code was implemented and tested)

  !Diagonal element of poloidal momentum diffusivity
  !This should be fine
  DMDIAG=D1*TE**1.5_DP*DMD/kyrho
  DM(4)=DMDIAG

  !Poloidal momentum velocity pinches
  !This looks like the Reynolds source terms
  HPT(4)=0.5_DP*DMH*(DMI1+DMI21+DMI22)

  !Poloidal momentum pinches expressed as a diffusivity
  DMEFTEST=-2._DP*HPT(4)*LVF

  !Total flux expressed as a diffusivity 
  DMEFF=DM(4) -HPT(4)*2._DP*LVF

  !Unintelligible stuff. Quite likely, the Reynolds stress term
  !is being expressed as a source term here
  SP=1._DP
  IF (vpol.LT.0.) SP=-1._DP
  SMEF1=SP*D1*TE**1.5_DP*DMI1/rmaj
  SMEF21=SP*D1*TE**1.5_DP*DMI21/rmaj   !! splitting of DMEF2
  SMEF22=SP*D1*TE**1.5_DP*DMI22/rmaj   !! splitting of DMEF2
!  Vconv(IK)=2._DP*(SMEF1+SMEF21+SMEF22)*(Csound/CSOUND0)
!  SRE=Vconv(IK)/(rmin*LVF)

  !End of momentum transport computations
  !-----------------------------------------------------------
  !Write effective diffusion and convection to returned arrays
  !

  diffs(1) = SCHI   !Ion heat diffusivity
  diffs(2) = SD     !Electron particle diffusivity?
  diffs(3) = SCHE   !Electron heat diffusivity
  diffs(4) = SDQ    !Impurity particle diffusivity
  diffs(5) = DMT(6) !Toroidal momentum diffusivity
  diffs(6) = DM(4)  !Poloidal momentum diffusivity
                    !Note that with the curvature effect omitted
                    !DM(4) is also the toroidal momentum diffusivity
  conv(1) = 0._DP 
  conv(2) = 0._DP 
  conv(3) = 0._DP 
  conv(4) = v_p_tor ! Sum of parellel velocity + Reynolds stress pinch
  conv(5) = v_p_thr ! Thermoelectric pinch
  conv(6) = hpt(4) ! Sum of poloidal momentum related convection

  return
end subroutine w20diff!}}}

subroutine w20disp &!{{{
   ( zamr, zami, zbmr, zbmi, wz, G_ave, kps )
  !use w20data
  implicit none

  !Equation matrices
  real(DP), dimension(neq,neq) :: zamr, zami, zbmr, zbmi

  !Input frequency
  complex(DP), intent(in) :: wz

  real(DP) RAV, G_ave
  real(DP) H1,XH,R,HQR,HQI

  complex(DP) ALPC,ALPHA,HQ,IU,H2
  complex(DP) A_lpk

  real(DP) H

  real(DP) :: alp, alf, kps
!  real(DP) alp,kps,alf

  real(DP) zdenom

  integer :: j1, j2

!
!    variables i=1,6: e Phi/Te, dTi/Ti, dni/ni, dTe/Te, dnq/nq, dTq/Tq
!    variables j=1,6 same as for i
!

  IU=(0._DP,1._DP)
!
!..print header
!
  if ( lprint .gt. 2 ) then
     write(hprint,*)
     write(hprint,*) &
          'Weiland-Nordman eigenvalue equations, subroutine etawn6'
     write(hprint,*) '(all frequencies normalized by omega_{De})'
     write(hprint,*) '(all diffusivities normalized by ' &
          ,'omega_{De} / k_y^2'
  endif


!
!..diagnostic output
!

  if ( lprint .gt. 6 ) then
     write(hprint,*)
     write(hprint,*) '--------------------------------------'
     write(hprint,*)
     write(hprint,*)
     write(hprint,*) gne,' = gne'
     write(hprint,*) gni,' = gni'
     write(hprint,*) gnz,' = gnz'
     write(hprint,*) gte,' = gte'
     write(hprint,*) gti,' = gti'
     write(hprint,*) gtz,' = gtz'
     write(hprint,*)
     write(hprint,*) vef,' = vef'
     write(hprint,*) betae, '=betae'
     write(hprint,*) zimp,' = zimp'
     write(hprint,*) aimp,' = aimp'
     write(hprint,*) 
     write(hprint,*) zflh,' = zflh'
     write(hprint,*) zflz,' = zflz'
     write(hprint,*)
     write(hprint,*) bt,' =bt'
     write(hprint,*) bt1,'  =bt1'
     write(hprint,*) kappa,'  =kappa'
     write(hprint,*) wz,'   =wz'
  endif

  !Reset system of equations to zero
  zamr(1:neq,1:neq) = 0._DP
  zami(1:neq,1:neq) = 0._DP
  zbmr(1:neq,1:neq) = 0._DP
  zbmi(1:neq,1:neq) = 0._DP

!
!..Nine  equations with impurities, trapped electrons, parallel ion
!  motion, collisions,  FLR , finite beta and parallel motion of impurities
!
!  Equations for e phi/Te, Ti, ni, Tet, n_z, T_z, F, Av, K
!  Here F = omega*gamma (collisions), K = omega*Av and Av is the parallel
!  magnetic vector potential.
!

  H1=4._DP*q*q*zflh

  H=0.5_DP*ABS(shat)/q
  H2=IU*H

  A_lpk=0.5_DP*shat*SQRT(H1*XT*zflh*(1.0_DP+TAUH*(gni+gti)/(2._DP*WZ)) )

  IF ( DREAL(A_lpk) .LT. 0) then
     A_lpk=-A_lpk
  END IF

  if ( abs(dreal(A_lpk)).lt.0.01_DP ) then
     A_lpk = A_lpk - dreal(A_lpk) + 0.01_DP
  end if

  ALPC=-IU*A_lpk
  ALPHA=-IU*ABS(shat)*q*zflh
  XH=ABS(ALPHA/ALPC)
  ALPC=XH*ALPC
  R=2._DP*ABS(DREAL(WZ*ALPC))
  IF (R.LT.0.001_DP) R=0.001_DP    !! NEW 01.03.8
  HQ=2._DP*ALPC/H1

  !Geometric average multiplying \omega_{De}
  G_ave=(1.0_DP+0.5_DP*shear/R)*EXP(-0.25_DP/R)-0.5_DP*alpha_MHD*(1._DP-EXP(-1._DP/R))
  G_ave=max(G_ave,1E-2_DP)

  alp=max(0.1_DP,0.5_DP*R)

  zdenom = (2._DP*zflh*q*q*betae*(1._DP - fte))
  alf = (alp+1E-6_DP) / (zdenom+1E-6_DP)
  kps = 0.5*SQRT(alp/zflh)/q
  kpc = 1._DP
  RAV=1._DP+0.25_DP*shat**2._DP/alp

!
!  *********
  HQR=DREAL(HQ)
  HQI=DIMAG(HQ)
!  *********


      if ( lprint > 6 ) then
         write(hprint,10002) alp,shat**2,RAV
10002    FORMAT(2X,'alp=',ES15.6,'  SH2=',ES15.6,' RAV=',ES15.6)
         write(hprint,10003) XH,G_ave,alf
10003    FORMAT(2X,'XH=',ES15.6,' G_ave=',ES15.6,' alf=',ES15.6)
         write(hprint,10005) alpha_MHD
10005    FORMAT(2X,'MHD Alpha=',ES15.6)
         write(hprint,10006) DREAL(A_lpk),R
10006    FORMAT(2X,' A_lpk=',ES15.6,' R='ES15.6)
         write(hprint,*)"HQR=",hqr,"   HQI=",hqi
      end if

!
!--- WE NOW DEFINE THE MATRIX ELEMENTS --------------
! hydrogen density
!
!
!...Equations rewritten in terms of the profile gradients
!
  zamr(1,1) = - G_ave+HQR + 0.5_DP * ( gni - zflh * ztauh * ( gni + gti ) )
  zami(1,1) = HQI
  zamr(1,2) = (HQR-G_ave)*ztauh
  zami(1,2) = ztauh*HQI
  zamr(1,3) = (HQR-G_ave)*ztauh
  zami(1,3) = ztauh*HQI
  zamr(1,8) = -em*ztauh*HQR*0.5_DP*(gni+gti)/kpc
  zami(1,8) = -em*ztauh*HQI*0.5_DP*(gni+gti)/kpc
  zamr(1,9) = -em*HQR/kpc
  zami(1,9) = -em*HQI/kpc
!
  zbmr(1,1) = zflh
  zbmr(1,3) = 1._DP
!
!  hydrogen energy
!
  zamr(2,1) = 0.5_DP*(gti - tvr * gni)
  zamr(2,2) = - ztauh * ftr
!
  zbmr(2,2) = 1._DP
  zbmr(2,3) = - tvr
!
!  total electron density expressed in ion density and imp density
!
  zamr(3,1) = -1.0_DP + 0.5_DP*fte*gne
  zami(3,1) = vef*(1._DP-fte)
  zamr(3,3) = 1._DP - znz_ne -zns_ne
  zami(3,3) = -vef*(1._DP - znz_ne - zns_ne)
  zamr(3,4) = fte
  zamr(3,5) = znz_ne
  zami(3,5) = -vef*znz_ne
  zami(3,7) = vef*fte
  zamr(3,8) = -em*0.5_DP*gne*(1._DP - fte)/kpc
  zami(3,8) =  em*0.5_DP*gne*(1._DP - fte)*vef/kpc
  zamr(3,9) =  em* (1._DP - fte)* (1._DP+0.5_DP*gne) / kpc
  zami(3,9) = -em* (1._DP - fte)* vef / kpc

  zbmr(3,1) = fte - 1._DP
  zbmr(3,3) = 1._DP - znz_ne - zns_ne
  zbmr(3,5) = znz_ne
  zbmr(3,9) = em*(1._DP - fte)/kpc
!
!  trapped electron energy
!
  zamr(4,1) = fte * 0.5_DP * ( gte - tvr*gne )
  zami(4,1) = vef*tvr*(bt-2.5_DP*(1._DP-fte))
  zami(4,3) = -vef*tvr*bt1*(1._DP-znz_ne -zns_ne)
  zamr(4,4) = fte*ftr
  zami(4,5) = -vef*tvr*bt1*znz_ne
  zami(4,7) = -ftr*vef*fte
!
  zbmr(4,1) = (1._DP - fte)*tvr
  zbmr(4,3) = -(1._DP - znz_ne -zns_ne)*tvr
  zbmr(4,4) = fte
  zbmr(4,5) = -znz_ne*tvr
!
!  impurity density
!
  zamr(5,1) = - G_ave + zimp*HQR/aimp        &
       +0.5_DP*( gnz - zflz*ztauz*(gnz+gtz) )
  zami(5,1) = zimp*HQI/aimp

  zamr(5,5) = (HQR*zimp/aimp-G_ave)*ztauz
  zami(5,5) = zimp*ztauz*HQI/aimp

  zamr(5,6) = (HQR*zimp/aimp-G_ave)*ztauz
  zami(5,6) = zimp*ztauz*HQI/aimp

  zamr(5,8) = -em*HQR*zimp*ztauz*0.5_DP*(gnz+gtz)/(kpc*aimp)
  zami(5,8) = -em*HQI*zimp*ztauz*0.5_DP*(gnz+gtz)/(kpc*aimp)

  zamr(5,9) = -em*HQR*zimp/(kpc*aimp)
  zami(5,9) = -em*HQI*zimp/(kpc*aimp)
 
  zbmr(5,1) = zflz
  zbmr(5,5) = 1._DP
!
!  impurity energy
!
  zamr(6,1) = 0.5_DP*(gtz - tvr * gnz)
  zamr(6,6) = -ztauz*ftr
!
  zbmr(6,5) = -tvr
  zbmr(6,6) = 1._DP
!
!  variable F
!
  zamr(7,1) = 0.5_DP*gte - 1.0_DP
  zami(7,1) = vef
  ZAMR(7,7) = 1._DP
  zami(7,7) = -vef
!
  zbmr(7,1) = -1._DP
  zbmr(7,7) = 1._DP
!
!
!  electromagnetic parallel vectorpotential Av = e A_par/Te
!
  zamr(8,1) = em1*kpc*(0.5_DP*gne+HQR*(fft+zimp*fzft/aimp))
  zami(8,1) = em1*HQI*(fft+zimp*fzft/aimp)*kpc

  zamr(8,2) = em1*HQR*ztauh*fft*kpc
  zami(8,2) = em1*HQI*ztauh*fft*kpc

  zamr(8,3) = em1*HQR*ztauh*fft*kpc
  zami(8,3) = em1*HQI*ztauh*fft*kpc

  zamr(8,5) = em1*HQR*zimp*ztauz*fzft*kpc/aimp
  zami(8,5) = em1*HQI*zimp*ztauz*fzft*kpc/aimp

  zamr(8,6) = em1*HQR*zimp*ztauz*fzft*kpc/aimp
  zami(8,6) = em1*HQI*zimp*ztauz*fzft*kpc/aimp

  zamr(8,8) = em1*(0.5_DP*(gne + gte) - alf*zflh*RAV ) &
       - em1*HQR*(fft*ztauh*0.5_DP*(gni+gti)           &
       + zimp*fzft*ztauz*0.5_DP*(gtz+gnz)/aimp)

  zami(8,8) = -em1*HQI*(fft*ztauh*0.5_DP*(gni + gti)   &
       +zimp*fzft*ztauz*0.5_DP*(gnz + gtz)/aimp)

  zamr(8,9) = -em1*(0.5_DP*gne+HQR*(fft+zimp*fzft/aimp))
  zami(8,9) = -em1*HQI*(fft+zimp*fzft/aimp)

  zbmr(8,1) = em1*kpc
  zbmr(8,8) = em1
  zbmr(8,9)= -em1
!
!     K = omega*Av
!
  zamr(9,9) = em1
!
  zbmr(9,8) = em1
!
      if ( lprint .gt. 6 ) then
        write(hprint,*)
        write(hprint,*) ' zamr(j1,j2)  j2 ->'
        do j1=1,9
          write(hprint,192) (zamr(j1,j2),j2=1,9)
        enddo

        write(hprint,*)

        write(hprint,*) ' zami(j1,j2)  j2 ->'
        do j1=1,9
          write(hprint,192) (zami(j1,j2),j2=1,9)
        enddo

       write(hprint,*)
        write(hprint,*) ' zbmr(j1,j2)  j2->'
        do j1=1,9
          write(hprint,192) (zbmr(j1,j2),j2=1,9)
        enddo

       write(hprint,*)
       write(hprint,*) ' zbmi(j1,j2)  j2->'
       do j1=1,9
         write(hprint,192) (zbmi(j1,j2),j2=1,9)
       enddo

 192  format (1p10e12.4)
      endif

End subroutine w20disp!}}}

Subroutine w20solv &!{{{
   ( zz, zvr, zvi, iu, eu, imax, emax, G_ave, kps, ierr )

  implicit none

  !Maximum number of iterations
  integer, parameter :: nitmax = 50

  !Eigenvectors resulting from solved system of equations
  real(DP), intent(inout), dimension(neq,neq) :: zvr, zvi

  !Complex frequency of the mode 
  complex(DP), intent(inout), dimension(neq) :: zz

  !Pointing indices (im-->unstable mode used, iu-->number of unstable modes)
  integer, intent(out) :: iu, eu, imax, emax

  !External error checking flag
  integer, intent(out) :: ierr

  !Geometric average needed to compute momentum pinches
  real(DP) :: G_ave, G_ave_tmp, kps, kps_tmp

  !Internal error checking flag for eigenvalue solver
  integer :: ifail, it_err, it_conv

  !Current number of iterations
  integer :: niter

  !Iteration indices
  integer :: j, nmod

  !Skip mode flag
  integer :: iskip, foundmode

  !Complex frequencies of fastest growing eigenmode
  complex(DP), dimension(nitmax+1) :: wz, xbest, wzfnd

  !Input frequency for dispersion relation
  complex(DP) :: wzin, dw1, dw

  !Temporary storage for growth rates and eigenmodes
  real(DP)                     :: wimax, wemax, wamax
  real(DP), dimension(neq,neq) :: zvrtmp, zvitmp
  complex(DP), dimension(neq) :: zztmp

  !Mask array for active modes
  real(DP),  dimension(neq) :: mask

  !Variables dealing with error control
  real(DP) ::  zepsilon, ztol
  real(DP) ::  delWZ, delZ

  !Pointer to mode being used
  integer :: nmod_ptr

  !Initial count of unstable modes
  integer :: unst

  !Final test for largest growing mode in the other direction
  real(DP) :: womax, woomax, wiimax
  integer :: omax, ou, au, amax

  !Matrices and vectors for eigenvalue system
  real(DP), dimension(neq) :: omega, gamma, zbeta, zalfi, zalfr
  real(DP), dimension(neq,neq) :: zamr, zami, zbmr, zbmi
  
  au = 0
  iu = 0
  eu = 0
  ou = 0
  unst = 0
  amax = 0
  imax = 0
  emax = 0
  omax = 0
  ierr = 0
  zepsilon = 1E-6_DP
  ztol = 1E-3_DP
  zz(1:neq)=(0._DP,0._DP)
  mask(1:neq)=0._DP
  wamax = 0._DP
  wimax = 0._DP
  wemax = 0._DP
  womax = 0._DP
  wiimax = 0._DP
  woomax = 0._DP
  it_err = 0
  it_conv = 0

  modes:  do nmod=1,neq

     !Initialization of some parameters
     nmod_ptr = 0
     niter = 0
     ierr  = 0
     wz(1:nitmax+1)=(0._DP,0._DP)
     wzin=(0._DP,0._DP)
     wzfnd(1:nitmax+1)=(0._DP,0._DP)
     zztmp(1:neq)=(0._DP,0._DP)
     xbest(1:nitmax+1)=(0._DP,0._DP)
     delZ = 0._DP
     delWZ = 0._DP
     dw = (0._DP,0._DP)
     dw1 = (0._DP,0._DP)

     call w20wguess( wzin )

        if ( searchmode .eq. S_ELC .and. dreal(wzin).lt.0._DP ) then
           wzin=wzin+2*abs(dreal(wzin))
        end if

     wz(1)    = wzin
     xbest(1) = wzin     

     iskip = 0

     !Algorithm to find converged frequency begins here
     freq: Do
        !Initialize complex frequencies

        niter = niter + 1

        omega(1:neq)=0._DP
        gamma(1:neq)=0._DP


        !Maximum number of iterations 
        If ( niter .gt. nitmax ) Then
           it_err = 1
           If ( 0.1_DP*abs(delWZ) .lt. ztol ) then
              exit freq
           !Hard error if tolerance not found within an order of magnitude
           Else If ( 0.1_DP*abs(delWZ) .ge. ztol ) then
              ierr = 2
              it_err = 2
              Cycle modes
           End If
        End If

        !Obtain dispersion relation
        call w20disp(zamr, zami, zbmr, zbmi, wzin, G_ave_tmp, kps_tmp )

        !Diagonalize equation system
        call DPtomsqz_mmm(neq,neq,zamr,zami,zbmr,zbmi,zalfr,zalfi,zbeta,zvrtmp,zvitmp,ifail)

        !Check for errors
        If ( ifail .gt. 0 ) Then
           write(hprint,*) "MMM7.1|w20solv|ERROR in eigenvalue solver"
           ierr = 1
           exit modes
        End If

        zztmp = w20omg ( zalfr, zalfi, zbeta, neq )
        omega = dreal(zztmp)
        gamma = dimag(zztmp)

        !Begin carry out a series out tests to improve convergence and reliability

        !Test 1: if all modes are unstable, then return from this call
        unst = w20nunst ( zztmp, neq, S_ANY )

        if ( unst .eq. 0 .and. niter .ge. 3) then 
           it_conv = 1
           cycle modes
        end if

        !Test 2: lock onto the mode being looked at
        !        It should be nmod, but sometimes the solver rotates eigenvalues
        !        Thus find the mode that is most like the mode last used
        !        In the first iteration, choose the mode to be nmod
        if ( niter .eq. 1 ) then
           nmod_ptr = nmod
        else
           !For niter > 1, lock onto the mode last used
           nmod_ptr = w20wmatch ( wzfnd(niter-1), zztmp, neq )
        end if

        !Save the matched mode for the next iteration
        wzfnd(niter)=zztmp(nmod_ptr)

        !Test 3: if the mode is stable, go to the next mode
        If ( gamma(nmod_ptr) .lt. 0._DP .and. niter .ge. 3 ) then
           it_conv = 1
           cycle modes
        end if

        foundmode = w20s_fnd ( zztmp(nmod_ptr) )

        !Test 4: if the mode is in the wrong direction, go to the next mode
        if ( foundmode .ne. searchmode .and. niter .gt. 3 ) then
           it_conv = 1
           cycle modes
        end if
        !Choose the next trial frequency as the average of the old and the new mode
        if ( niter .le. 1 ) then
           wz(niter+1)   = zztmp(nmod_ptr)
        else
           wz(niter+1)   = 0.5_DP*(wzin+(0._DP,1._DP)*gamma(nmod_ptr)+omega(nmod_ptr))
        end if


        wzin          = wz(niter+1)
        Xbest(niter+1)= wzin

        !Compute estimate of error
        delWZ = w20delwz( Xbest(niter+1), Xbest(niter) )

        !Exit inner loop (freq) when convergence is found
        if ( delWZ .lt. ztol .and. niter .ge. 2) then
           it_conv = 1
           Exit freq
        end if

     End Do freq

     !If the unstable mode is in the direction we are looking for
     ! then evaluate whether this mode is faster growing than the 
     ! previously stored mode. If so, carry out the replacement

     foundmode = w20s_fnd ( zztmp(nmod_ptr) )
     
     if ( foundmode .eq. searchmode ) then
        if ( w20gamma( zztmp, neq, nmod_ptr ) .gt. wamax ) then       
           wamax = w20gamma( zztmp, neq, nmod_ptr )
           amax  = nmod_ptr
           G_ave = G_ave_tmp
           kps   = kps_tmp

           do j=1,neq
              foundmode = w20s_fnd( zztmp(j) )
              !if ( j .ne. amax .and. foundmode .eq. searchmode .and. w20gamma( zztmp, neq, j ) .gt. wamax ) then
              !  zztmp(j)=zztmp(j)-(0._DP,1.0)*zztmp(j)
              !end if
           end do

           imax  = w20wmunst ( zztmp, neq, S_ION ) !nmod_ptr
           wimax = w20gamma ( zztmp, neq, imax ) !wamax
           emax  = w20wmunst ( zztmp, neq, S_ELC )
           wemax = w20gamma ( zztmp, neq, emax )

           eu  = w20nunst( zztmp, neq, S_ELC )
           iu  = w20nunst( zztmp, neq, S_ION )

           zz  = zztmp
           zvr = zvrtmp
           zvi = zvitmp
        end if
              
     end if

  End Do modes

  if ( it_conv .ne. 1 ) then
     write(hprint,*) "MMM7.1|w20solv|ERROR: Excessive iterations in Weiland model"    

     if ( it_err .gt. 1 ) then
        write(hprint,*) "mode=",nmod,"  ptr=",nmod_ptr
        write(hprint,*) "delWZ=",delWZ
        write(hprint,*) "GRDNE=",gne
        write(hprint,*) "GRDTE=",gte
        write(hprint,*) "GRDNI=",gni
        write(hprint,*) "GRDTI=",gti
        write(hprint,*) "GRDNZ=",gnz
        write(hprint,*) "GRDTZ=",gtz
        write(hprint,*) "NE=",ne
        write(hprint,*) "TE=",te
        write(hprint,*) "NH=",nh
        write(hprint,*) "TI=",th
        write(hprint,*) "NZ=",nz
        write(hprint,*) "NS=",ns_ne*ne
        write(hprint,*) "Q=",q
        write(hprint,*) "Shear=",shear
        write(hprint,*) "Elong=",kappa
        write(hprint,*) "Aimp=",aimp
        write(hprint,*) "Zimp=",zimp
        write(hprint,*) "Ahyd=",ahyd
        write(hprint,*) "Rmaj=",rmaj
        write(hprint,*) "Rmin=",rmin
        write(hprint,*) "btor=",btor
     end if
  end if

End Subroutine w20solv!}}}

function w20wmunst &!{{{
  ( zz, nmod, S_SRCH ) result ( nmod_ptr )
!This function returns the position of the mode with
!the largest growthrate
!
!INPUT: zz(nmod): eigenvalues of drift mode dispersion relation )
!       nmode   : dimension of the vector
!       S_SRCH  : direction S_ION or S_ELC or S_ANY
!
!OUTPUT: nmod_ptr : position of the largest growing mode in
!                   direction S_SRCH in array zz
!                   will return 0 if no 
implicit none

integer, intent(in) :: nmod
complex(DP), intent(in), dimension(nmod) :: zz
integer, intent(in) :: S_SRCH

integer :: nmod_ptr, i, S_FND 
real(DP) :: gamma_max

gamma_max = -10e3_DP
S_FND = 0
nmod_ptr = N_NONE

do i=1,nmod

   if ( S_SRCH .ne. S_ANY ) then
      S_FND = w20s_fnd( zz(i) )
   end if
   
   if ( aimag(zz(i)) > gamma_max .and. S_FND .eq. S_SRCH ) then
      gamma_max = aimag(zz(i))
      nmod_ptr = i
   end if
end do

end function w20wmunst!}}}

function w20s_fnd &!{{{
   ( zz ) result ( S_FND )
!This function takes in a complex frequency and
!evaluates whether the mode is an electron or ion mode
!
! INPUT: zz = complex frequency ( w_r, gamma )
! OUTPUT: S_FND = S_ELC or S_ION ( integer flag for electron or ion mode)
!
implicit none

complex(DP), intent(in) :: zz

integer :: S_FND

if ( real(zz,kind=DP) > 0._DP ) then
   S_FND = S_ELC
else
   S_FND = S_ION
end if

end function w20s_fnd!}}}

function w20nunst &!{{{
   ( zz, nmod, S_SRCH ) result ( iunst )
!This function takes in a set of complex frequencies
!and counts the number of unstable modes in any direction
!
! INPUT: zz(nmod) = complex frequency
!        nmod     = dimension of zz vector
!        S_SRCH   = direction of search S_ANY, S_ION, S_ELC
!
! OUTPUT: iunst   = number of unstable modes in S_SRCH direction

implicit none

integer, intent(in) :: nmod, S_SRCH
complex(DP), intent(in), dimension(nmod) :: zz

integer :: iunst, i, S_FND

iunst = N_NONE
S_FND = S_ANY

do i=1,nmod
   
   if ( S_SRCH .ne. S_ANY ) then
      S_FND = w20s_fnd( zz(i) )
   end if
   
   if ( aimag(zz(i)) > 0._DP .and. S_SRCH .eq. S_FND ) then
      iunst = iunst + 1
   end if
   
end do

end function w20nunst!}}}

function w20gamma &!{{{
  ( zz, nmod, nmod_ptr ) result ( gamma )
!This function returns the growthrate of the drift mode
!at the position nmod_ptr. If nmod_ptr is not between
!1 and nmod, then the result is zero. This is useful
!when no unstable modes are found in one of the functions
!above
!
!INPUT: zz(nmod) = complex eigenfrequencies
!       nmode    = dimension of zz vector
!       nmod_ptr = position of the growthrate
!
implicit none

integer, intent(in) :: nmod, nmod_ptr
complex(DP), intent(in), dimension(nmod) :: zz

real(DP) :: gamma

if ( nmod_ptr .eq. N_NONE ) then
   gamma = 0._DP
else if ( nmod_ptr .ge. 1 .and. nmod_ptr .le. nmod ) then
   gamma = aimag(zz(nmod_ptr))
else
   write (*,*) "MMM7.1|w20gamma|ERROR: Invalid mode number: ",nmod_ptr
end if

end function w20gamma!}}}

function w20omg &!{{{
   ( zomega, zgamma, zbeta, nmod ) result ( zz )
!This function takes in the real and imaginary parts of the
!eigenvalues and folds them into a complex frequency
!The frequencies are then normalized

implicit none

integer, intent(in) :: nmod
real(DP) , intent(in), dimension(nmod) :: zomega, zgamma, zbeta

complex(DP), dimension(nmod) :: zz
real(DP), dimension(nmod) :: zb

zz = (0._DP,0._DP)
zb = max(1E-4_DP,zbeta)

zz(1:nmod) = ( zomega(1:nmod)+(0._DP,1.0_DP)*zgamma(1:nmod) ) / zb(1:nmod)

end function w20omg!}}}

function w20delwz &!{{{
   ( wz1, wz2 ) result ( delwz )
!This function computes the relative difference between two modes

implicit none

complex(DP) :: wz1, wz2

real(DP) :: delwz,denom
real(DP) :: eps = 0.001_DP
real(DP) :: wr1, wr2, wi1, wi2

eps = 0.001_DP
wr1 = real(wz1,kind=DP)
wr1 = sign( max( abs(wr1), eps ), wr1)
wi1 = aimag(wz1)
wi1 = sign( max( abs(wi1), eps ), wi1)

wr2 = real(wz2,kind=DP)
wi2 = aimag(wz2)

denom = ( max( abs(wz1), eps ) )**2
denom = 0.5_DP * ( abs(wz1) + abs(wz2) )

delwz = sqrt( ( wr1-wr2 )**2 + ( wi1-wi2 )**2 ) / denom

end function w20delwz!}}}

function w20wmatch &!{{{
   ( wz, zz, nmod ) result ( nmod_ptr )
!This function returns the mode closest to the input mode
implicit none

integer, intent(in) :: nmod
complex(DP), intent(in) :: wz
complex(DP), intent(in), dimension(nmod):: zz

integer :: nmod_ptr,i
real(DP) :: delwz, delwzmin

nmod_ptr = N_NONE
delwzmin = 10e6_DP
delwz = 0._DP

do i=1,nmod
   delwz = w20delwz( wz, zz(i) )
   if ( delwz .lt. delwzmin ) then
      delwzmin = delwz
      nmod_ptr = i
   end if
end do

end function w20wmatch!}}}

subroutine w20wguess &!{{{
   ( wz )

implicit none

complex(DP), intent(inout) :: wz

complex(DP) :: wz1, wz2, H_q, E1, E2
complex(DP) :: iu
real(DP)     :: G_ave

G_ave = 1.0_DP
IU = (0._DP, 1.0_DP)

H_q = cmplx(0._DP,0.5_DP*shat/q)

E1 = ftrt*(1.0_DP+zflh)-0.5_DP*gne+0.5_DP*zflh*tauh*(gne+gti) &
  + (gm+ftrt)*G_ave+H_q*(1.0_DP+ftr*tauh)

E1 = E1 * 0.5_DP / ( 1.0_DP + zflh )


E2 =  (  ( 0.5_DP * tauh * gm * ( gti - tvr*gne ) + bta )*( G_ave + H_q ) &
    - 0.5_DP * tauh * ftr * ( gni - zflh * tauh * (gti + gni))  ) &
    / ( 1.0_DP + zflh )

wz1=-E1+SQRT(E1*E1-E2)
wz2=-E1-SQRT(E1*E1-E2)

wz=wz1

IF(DIMAG(wz2).GT.DIMAG(wz1)) wz=wz2

if ( dimag(wz) .le. 0.01_DP ) then
  wz = wz - iu*dimag(wz) + iu*0.01_DP
end if

end subroutine w20wguess!}}}


! Original DPTOMSQZ.F90 {{{
!***********************************************************************
!
! Changes:
! May 04 2011  Lixiang Luo
!    subroutine remained to DPTOMSQZ_MMM to avoid namespace conflicts
!
!***********************************************************************
!                                                                       
!dmc -- real(DP) version generated using `fgtok'.                         
!  names changed:  tomsqz -> DPtomsqz, cqzhes -> DPcqzhes, etc.         
!  all constants converted to "D" expontent form                        
!  all declarations remapped to real(DP) / complex(DP)                     
!  cpp for standardizing REAL/COMPLEX intrinsics:                       
!                                                                       
!#include "f77_dcomplx.h"                                               
!                                                                       
SUBROUTINE DPTOMSQZ_MMM(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI,IFAIL) 
!-----------------------------------------------------------------------
! TOMSQZ  written by P. Strand 27-apr-98,       elfps@elmagn.chalmers.se
!-----------------------------------------------------------------------
! CQZHES, CQZVEC, and CQZVAL: Fortran subroutines implementing the QZ   
! algorithm for solving the generalized eigenvalue problem for complex  
! matrices. (See B.S.C Garbow, ACM TOMS 4 (1978) pp. 404-410.).         
!-----------------------------------------------------------------------
!                                                                       
! ON INPUT THE GENERALIZED EIGENVALUE PROBLEM IS DEFINED THROUGH THE    
! COMPLEX MATRICES                                                      
!                                                                       
!       A = cmplx (AR, AI)  AND   B = cmplx (BR, BI)                    
!                                                                       
! WHERE LEADING DIMENSION N IS AS DEFINED IN THE CALLING ROUTINE AND    
! WHERE NA IS THE ROW  RANGE IN THE CURRENT PROBLEM. THE EIGENVALUE     
! PROBLEM IS THEN DEFINED THROUGH                                       
!                                                                       
!       A x = w B x                                                     
!                                                                       
! WHERE  THE COMPLEX EIGENVECTORS                                       
!                                                                       
!       x = cmplx (ZVR, ZVI)                                            
!                                                                       
! TOGETHER WITH THE COMPLEX EIGENVALUE                                  
!                                                                       
!        w = cmplx(alfr, alfi)/beta                                     
!                                                                       
! IS OUTPUT FROM THE ROUTINE                                            
!                                                                       
! IFAIL WILL BE NONZERO IF CONVERGENCE HAS NOT BEEN REACH WITHIN 50     
! ITERATIONS                                                            
!-----------------------------------------------------------------------
! DECLARATIONS FOR INPUT VARIABLES                                      
!-----------------------------------------------------------------------
                                                                        
      IMPLICIT NONE 
      INTEGER N, NA 
      real(DP) AR(N,NA),AI(N,NA),BR(N,NA),BI(N,NA) 
                                                                        
!-----------------------------------------------------------------------
! DECALRATIONS FOR OUTPUT VARIABLES                                     
!-----------------------------------------------------------------------
                                                                        
      real(DP) ALFR(N),ALFI(N),BETA(N) 
      real(DP) ZVR(N,NA), ZVI(N,NA) 
      INTEGER IFAIL 
                                                                        
!-----------------------------------------------------------------------
! LOCAL VARIABLES                                                       
!-----------------------------------------------------------------------
                                                                        
      LOGICAL WANTX 
      real(DP) EPS1 
                                                                        
!-----------------------------------------------------------------------
! START OF ACTUAL CODING                                                
!-----------------------------------------------------------------------
                                                                        
      WANTX = .TRUE. 
      EPS1  = -0.0D0 
                                                                        
      CALL DPCQZHES_MMM(N,NA,AR,AI,BR,BI,WANTX,ZVR,ZVI) 
      CALL DPCQZVAL_MMM(N,NA,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,WANTX,         &
     &            ZVR,ZVI,IFAIL)                                        
      CALL DPCQZVEC_MMM(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI) 
      RETURN 
END SUBROUTINE DPTOMSQZ_MMM
                                                                        
!     ------------------------------------------------------------------
!                                                                       
SUBROUTINE DPCQZHES_MMM(NM,N,AR,AI,BR,BI,MATZ,ZR,ZI) 
!                                                                       
      IMPLICIT NONE 
      INTEGER I,J,K,L,N,K1,LB,L1,NM,NK1,NM1 
      real(DP) AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ZR(NM,N),ZI(NM,N) 
      real(DP) R,S,T,TI,U1,U2,XI,XR,YI,YR,RHO,U1I 
      LOGICAL MATZ 
      real(DP) ZERO,ZONE 
!                                                                       
!     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FIRST STEP OF THE    
!     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,  
!     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.      
!                                                                       
!     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX GENERAL MATRICES AND    
!     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM WITH REAL (AND NON-  
!     NEGATIVE) SUBDIAGONAL ELEMENTS AND THE OTHER TO UPPER TRIANGULAR  
!     FORM USING UNITARY TRANSFORMATIONS.  IT IS USUALLY FOLLOWED BY    
!     CQZVAL  AND POSSIBLY  CQZVEC.                                     
!                                                                       
!     ON INPUT-                                                         
!                                                                       
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
!          DIMENSION STATEMENT,                                         
!                                                                       
!        N IS THE ORDER OF THE MATRICES,                                
!                                                                       
!        A=(AR,AI) CONTAINS A COMPLEX GENERAL MATRIX,                   
!                                                                       
!        B=(BR,BI) CONTAINS A COMPLEX GENERAL MATRIX,                   
!                                                                       
!        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS 
!          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING             
!          EIGENVECTORS, AND TO .FALSE. OTHERWISE.                      
!                                                                       
!     ON OUTPUT-                                                        
!                                                                       
!        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS     
!          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO, AND THE   
!          SUBDIAGONAL ELEMENTS HAVE BEEN MADE REAL (AND NON-NEGATIVE), 
!                                                                       
!        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS     
!          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,               
!                                                                       
!        Z=(ZR,ZI) CONTAINS THE PRODUCT OF THE RIGHT HAND               
!          TRANSFORMATIONS IF MATZ HAS BEEN SET TO .TRUE.               
!          OTHERWISE, Z IS NOT REFERENCED.                              
!                                                                       
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
!     ********** INITIALIZE Z **********                                
                                                                        
      ZERO = 0.0D0 
      ZONE = 1.0D0 
      IF (.NOT. MATZ) GO TO 10 
!                                                                       
      DO 3 I = 1, N 
!                                                                       
         DO 2 J = 1, N 
            ZR(I,J) = 0.0D0 
            ZI(I,J) = 0.0D0 
    2    CONTINUE 
!                                                                       
         ZR(I,I) = 1.0D0 
    3 END DO 
!     ********** REDUCE B TO UPPER TRIANGULAR FORM WITH                 
!                TEMPORARILY REAL DIAGONAL ELEMENTS **********          
   10 IF (N .LE. 1) GO TO 170 
      NM1 = N - 1 
!                                                                       
      DO 100 L = 1, NM1 
         L1 = L + 1 
         S = 0.0D0 
!                                                                       
         DO 20 I = L, N 
            S = S + ABS(BR(I,L)) + ABS(BI(I,L)) 
   20    CONTINUE 
!                                                                       
         IF (S /= ZERO) Then  ! otherwise go to 100
         RHO = 0.0D0 
!                                                                       
         DO 25 I = L, N 
            BR(I,L) = BR(I,L) / S 
            BI(I,L) = BI(I,L) / S 
            RHO = RHO + BR(I,L)**2 + BI(I,L)**2 
   25    CONTINUE 
!                                                                       
         R = SQRT(RHO) 
         XR = abs(CMPLX(BR(L,L),BI(L,L))) 
         IF (XR .EQ. ZERO) GO TO 27 
         RHO = RHO + XR * R 
         U1 = -BR(L,L) / XR 
         U1I = -BI(L,L) / XR 
         YR = R / XR + 1.0D0 
         BR(L,L) = YR * BR(L,L) 
         BI(L,L) = YR * BI(L,L) 
         GO TO 28 
!                                                                       
   27    BR(L,L) = R 
         U1 = -1.0D0 
         U1I = 0.0D0 
!                                                                       
   28    DO 50 J = L1, N 
            T = 0.0D0 
            TI = 0.0D0 
!                                                                       
            DO 30 I = L, N 
               T = T + BR(I,L) * BR(I,J) + BI(I,L) * BI(I,J) 
               TI = TI + BR(I,L) * BI(I,J) - BI(I,L) * BR(I,J) 
   30       CONTINUE 
!                                                                       
            T = T / RHO 
            TI = TI / RHO 
!                                                                       
            DO 40 I = L, N 
               BR(I,J) = BR(I,J) - T * BR(I,L) + TI * BI(I,L) 
               BI(I,J) = BI(I,J) - T * BI(I,L) - TI * BR(I,L) 
   40       CONTINUE 
!                                                                       
            XI = U1 * BI(L,J) - U1I * BR(L,J) 
            BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J) 
            BI(L,J) = XI 
   50    CONTINUE 
!                                                                       
         DO 80 J = 1, N 
            T = 0.0D0 
            TI = 0.0D0 
!                                                                       
            DO 60 I = L, N 
               T = T + BR(I,L) * AR(I,J) + BI(I,L) * AI(I,J) 
               TI = TI + BR(I,L) * AI(I,J) - BI(I,L) * AR(I,J) 
   60       CONTINUE 
!                                                                       
            T = T / RHO 
            TI = TI / RHO 
!                                                                       
            DO 70 I = L, N 
               AR(I,J) = AR(I,J) - T * BR(I,L) + TI * BI(I,L) 
               AI(I,J) = AI(I,J) - T * BI(I,L) - TI * BR(I,L) 
   70       CONTINUE 
!                                                                       
            XI = U1 * AI(L,J) - U1I * AR(L,J) 
            AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J) 
            AI(L,J) = XI 
   80    CONTINUE 
!                                                                       
         BR(L,L) = R * S 
         BI(L,L) = 0.0D0 
!                                                                       
         DO 90 I = L1, N 
            BR(I,L) = 0.0D0 
            BI(I,L) = 0.0D0 
   90    CONTINUE 

        End If
!                                                                       
  100 END DO 
!     ********** REDUCE A TO UPPER HESSENBERG FORM WITH REAL SUBDIAGONAL
!                ELEMENTS, WHILE KEEPING B TRIANGULAR **********        
      DO 160 K = 1, NM1 
         K1 = K + 1 
!     ********** SET BOTTOM ELEMENT IN K-TH COLUMN OF A REAL ********** 
         IF (AI(N,K) .EQ. ZERO) GO TO 105 
         R = abs(CMPLX(AR(N,K),AI(N,K))) 
         U1 = AR(N,K) / R 
         U1I = AI(N,K) / R 
         AR(N,K) = R 
         AI(N,K) = 0.0D0 
!                                                                       
         DO 103 J = K1, N 
            XI = U1 * AI(N,J) - U1I * AR(N,J) 
            AR(N,J) = U1 * AR(N,J) + U1I * AI(N,J) 
            AI(N,J) = XI 
  103    CONTINUE 
!                                                                       
         XI = U1 * BI(N,N) - U1I * BR(N,N) 
         BR(N,N) = U1 * BR(N,N) + U1I * BI(N,N) 
         BI(N,N) = XI 
  105    IF (K .EQ. NM1) GO TO 170 
         NK1 = NM1 - K 
!     ********** FOR L=N-1 STEP -1 UNTIL K+1 DO -- **********           
         DO 150 LB = 1, NK1 
            L = N - LB 
            L1 = L + 1 
!     ********** ZERO A(L+1,K) **********                               
            S = ABS(AR(L,K)) + ABS(AI(L,K)) + AR(L1,K) 
            IF (S .EQ. ZERO) GO TO 150 
            U1 = AR(L,K) / S 
            U1I = AI(L,K) / S 
            U2 = AR(L1,K) / S 
            R = SQRT(U1*U1+U1I*U1I+U2*U2) 
            U1 = U1 / R 
            U1I = U1I / R 
            U2 = U2 / R 
            AR(L,K) = R * S 
            AI(L,K) = 0.0D0 
            AR(L1,K) = 0.0D0 
!                                                                       
            DO 110 J = K1, N 
               XR = AR(L,J) 
               XI = AI(L,J) 
               YR = AR(L1,J) 
               YI = AI(L1,J) 
               AR(L,J) = U1 * XR + U1I * XI + U2 * YR 
               AI(L,J) = U1 * XI - U1I * XR + U2 * YI 
               AR(L1,J) = U1 * YR - U1I * YI - U2 * XR 
               AI(L1,J) = U1 * YI + U1I * YR - U2 * XI 
  110       CONTINUE 
!                                                                       
            XR = BR(L,L) 
            BR(L,L) = U1 * XR 
            BI(L,L) = -U1I * XR 
            BR(L1,L) = -U2 * XR 
!                                                                       
            DO 120 J = L1, N 
               XR = BR(L,J) 
               XI = BI(L,J) 
               YR = BR(L1,J) 
               YI = BI(L1,J) 
               BR(L,J) = U1 * XR + U1I * XI + U2 * YR 
               BI(L,J) = U1 * XI - U1I * XR + U2 * YI 
               BR(L1,J) = U1 * YR - U1I * YI - U2 * XR 
               BI(L1,J) = U1 * YI + U1I * YR - U2 * XI 
  120       CONTINUE 
!     ********** ZERO B(L+1,L) **********                               
            S = ABS(BR(L1,L1)) + ABS(BI(L1,L1)) + ABS(BR(L1,L)) 
            IF (S .EQ. ZERO) GO TO 150 
            U1 = BR(L1,L1) / S 
            U1I = BI(L1,L1) / S 
            U2 = BR(L1,L) / S 
            R = SQRT(U1*U1+U1I*U1I+U2*U2) 
            U1 = U1 / R 
            U1I = U1I / R 
            U2 = U2 / R 
            BR(L1,L1) = R * S 
            BI(L1,L1) = 0.0D0 
            BR(L1,L) = 0.0D0 
!                                                                       
            DO 130 I = 1, L 
               XR = BR(I,L1) 
               XI = BI(I,L1) 
               YR = BR(I,L) 
               YI = BI(I,L) 
               BR(I,L1) = U1 * XR + U1I * XI + U2 * YR 
               BI(I,L1) = U1 * XI - U1I * XR + U2 * YI 
               BR(I,L) = U1 * YR - U1I * YI - U2 * XR 
               BI(I,L) = U1 * YI + U1I * YR - U2 * XI 
  130       CONTINUE 
!                                                                       
            DO 140 I = 1, N 
               XR = AR(I,L1) 
               XI = AI(I,L1) 
               YR = AR(I,L) 
               YI = AI(I,L) 
               AR(I,L1) = U1 * XR + U1I * XI + U2 * YR 
               AI(I,L1) = U1 * XI - U1I * XR + U2 * YI 
               AR(I,L) = U1 * YR - U1I * YI - U2 * XR 
               AI(I,L) = U1 * YI + U1I * YR - U2 * XI 
  140       CONTINUE 
!                                                                       
            IF (.NOT. MATZ) GO TO 150 
!                                                                       
            DO 145 I = 1, N 
               XR = ZR(I,L1) 
               XI = ZI(I,L1) 
               YR = ZR(I,L) 
               YI = ZI(I,L) 
               ZR(I,L1) = U1 * XR + U1I * XI + U2 * YR 
               ZI(I,L1) = U1 * XI - U1I * XR + U2 * YI 
               ZR(I,L) = U1 * YR - U1I * YI - U2 * XR 
               ZI(I,L) = U1 * YI + U1I * YR - U2 * XI 
  145       CONTINUE 
!                                                                       
  150    CONTINUE 
!                                                                       
  160 END DO 
!                                                                       
  170 RETURN 
!     ********** LAST CARD OF CQZHES **********                         
END SUBROUTINE DPCQZHES_MMM
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
SUBROUTINE DPCQZVAL_MMM(NM,N,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,         &
     &                                       MATZ,ZR,ZI,IERR)           
!                                                                       
      IMPLICIT NONE 
      INTEGER I,J,K,L,N,EN,K1,K2,LL,L1,NA,NM,ITS,KM1,LM1,               &
     &        ENM2,IERR,LOR1,ENORN                                      
      real(DP) AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),ALFI(N),       &
     &       BETA(N),ZR(NM,N),ZI(NM,N)                                  
      real(DP) R,S,A1,A2,EP,SH,U1,U2,XI,XR,YI,YR,ANI,A1I,A33,A34,A43,A44, &
     &       BNI,B11,B33,B44,SHI,U1I,A33I,A34I,A43I,A44I,B33I,B44I,     &
     &       EPSA,EPSB,EPS1,ANORM,BNORM,B3344,B3344I                    
      INTEGER max 
      LOGICAL MATZ 
      complex(DP) Z3 
!                                                                       
      real(DP) ZERO,ZONE,ZTWO 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF STEPS 2 AND 3 OF THE     
!     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,  
!     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,      
!     AS MODIFIED IN TECHNICAL NOTE NASA TN E-7305(1973) BY WARD.       
!                                                                       
!     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES, ONE OF THEM   
!     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM,  
!     THE HESSENBERG MATRIX MUST FURTHER HAVE REAL SUBDIAGONAL ELEMENTS.
!     IT REDUCES THE HESSENBERG MATRIX TO TRIANGULAR FORM USING         
!     UNITARY TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM     
!     OF THE OTHER MATRIX AND FURTHER MAKING ITS DIAGONAL ELEMENTS      
!     REAL AND NON-NEGATIVE.  IT THEN RETURNS QUANTITIES WHOSE RATIOS   
!     GIVE THE GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY      
!     CQZHES  AND POSSIBLY FOLLOWED BY  CQZVEC.                         
!                                                                       
!     ON INPUT-                                                         
!                                                                       
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
!          DIMENSION STATEMENT,                                         
!                                                                       
!        N IS THE ORDER OF THE MATRICES,                                
!                                                                       
!        A=(AR,AI) CONTAINS A COMPLEX UPPER HESSENBERG MATRIX           
!          WITH REAL SUBDIAGONAL ELEMENTS,                              
!                                                                       
!        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,          
!                                                                       
!        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.     
!          EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN      
!          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF   
!          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS    
!          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE      
!          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A     
!          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,       
!          BUT LESS ACCURATE RESULTS,                                   
!                                                                       
!        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS 
!          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING             
!          EIGENVECTORS, AND TO .FALSE. OTHERWISE,                      
!                                                                       
!        Z=(ZR,ZI) CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE        
!          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION              
!          BY  CQZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.       
!          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.        
!                                                                       
!     ON OUTPUT-                                                        
!                                                                       
!        A HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS     
!          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,               
!                                                                       
!        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS     
!          HAVE BEEN ALTERED.  IN PARTICULAR, ITS DIAGONAL HAS BEEN SET 
!          REAL AND NON-NEGATIVE.  THE LOCATION BR(N,1) IS USED TO      
!          STORE EPS1 TIMES THE NORM OF B FOR LATER USE BY  CQZVEC,     
!                                                                       
!        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE      
!          DIAGONAL ELEMENTS OF THE TRIANGULARIZED A MATRIX,            
!                                                                       
!        BETA CONTAINS THE REAL NON-NEGATIVE DIAGONAL ELEMENTS OF THE   
!          CORRESPONDING B.  THE GENERALIZED EIGENVALUES ARE THEN       
!          THE RATIOS ((ALFR+I*ALFI)/BETA),                             
!                                                                       
!        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS       
!          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE.,             
!                                                                       
!        IERR IS SET TO                                                 
!          ZERO       FOR NORMAL RETURN,                                
!          J          IF AR(J,J-1) HAS NOT BECOME                       
!                     ZERO AFTER 50 ITERATIONS.                         
!                                                                       
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ZTWO = 2.0D0 
      ZONE = 1.0D0 
      ZERO = 0.0D0 
                                                                        
      IERR = 0 
!     ********** COMPUTE EPSA,EPSB **********                           
      ANORM = 0.0D0 
      BNORM = 0.0D0 
!                                                                       
      DO 30 I = 1, N 
         ANI = 0.0D0 
         IF (I .NE. 1) ANI = ABS(AR(I,I-1)) 
         BNI = 0.0D0 
!                                                                       
         DO 20 J = I, N 
            ANI = ANI + ABS(AR(I,J)) + ABS(AI(I,J)) 
            BNI = BNI + ABS(BR(I,J)) + ABS(BI(I,J)) 
   20    CONTINUE 
!                                                                       
         IF (ANI .GT. ANORM) ANORM = ANI 
         IF (BNI .GT. BNORM) BNORM = BNI 
   30 END DO 
!                                                                       
      IF (ANORM .EQ. ZERO) ANORM = 1.0D0 
      IF (BNORM .EQ. ZERO) BNORM = 1.0D0 
      EP = EPS1 
      IF (EP .GT. ZERO) GO TO 50 
!     ********** COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO **********      
      EP = 1.0D0 
   40 EP = EP / ZTWO 
      IF (ZONE + EP .GT. ZONE) GO TO 40 
   50 EPSA = EP * ANORM 
      EPSB = EP * BNORM 
!     ********** REDUCE A TO TRIANGULAR FORM, WHILE                     
!                KEEPING B TRIANGULAR **********                        
      LOR1 = 1 
      ENORN = N 
      EN = N 
!     ********** BEGIN QZ STEP **********                               
   60 IF (EN .EQ. 0) GO TO 1001 
      IF (.NOT. MATZ) ENORN = EN 
      ITS = 0 
      NA = EN - 1 
      ENM2 = NA - 1 
!     ********** CHECK FOR CONVERGENCE OR REDUCIBILITY.                 
!                FOR L=EN STEP -1 UNTIL 1 DO -- **********              
   70 DO 80 LL = 1, EN 
         LM1 = EN - LL 
         L = LM1 + 1 
         IF (L .EQ. 1) GO TO 95 
         IF (ABS(AR(L,LM1)) .LE. EPSA) GO TO 90 
   80 END DO 
!                                                                       
   90 AR(L,LM1) = 0.0D0 
!     ********** SET DIAGONAL ELEMENT AT TOP OF B REAL **********       
   95 B11 = abs(CMPLX(BR(L,L),BI(L,L))) 
      IF (B11     .EQ. ZERO) GO TO 98 
      U1 = BR(L,L) / B11 
      U1I = BI(L,L) / B11 
!                                                                       
      DO 97 J = L, ENORN 
         XI = U1 * AI(L,J) - U1I * AR(L,J) 
         AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J) 
         AI(L,J) = XI 
         XI = U1 * BI(L,J) - U1I * BR(L,J) 
         BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J) 
         BI(L,J) = XI 
   97 END DO 
!                                                                       
      BI(L,L) = 0.0D0 
   98 IF (L .NE. EN) GO TO 100 
!     ********** 1-BY-1 BLOCK ISOLATED **********                       
      ALFR(EN) = AR(EN,EN) 
      ALFI(EN) = AI(EN,EN) 
      BETA(EN) = B11 
      EN = NA 
      GO TO 60 
!     ********** CHECK FOR SMALL TOP OF B **********                    
  100 L1 = L + 1 
      IF (B11 .GT. EPSB) GO TO 120 
      BR(L,L) = 0.0D0 
      S = ABS(AR(L,L)) + ABS(AI(L,L)) + ABS(AR(L1,L)) 
      U1 = AR(L,L) / S 
      U1I = AI(L,L) / S 
      U2 = AR(L1,L) / S 
      R = SQRT(U1*U1+U1I*U1I+U2*U2) 
      U1 = U1 / R 
      U1I = U1I / R 
      U2 = U2 / R 
      AR(L,L) = R * S 
      AI(L,L) = 0.0D0 
!                                                                       
      DO 110 J = L1, ENORN 
         XR = AR(L,J) 
         XI = AI(L,J) 
         YR = AR(L1,J) 
         YI = AI(L1,J) 
         AR(L,J) = U1 * XR + U1I * XI + U2 * YR 
         AI(L,J) = U1 * XI - U1I * XR + U2 * YI 
         AR(L1,J) = U1 * YR - U1I * YI - U2 * XR 
         AI(L1,J) = U1 * YI + U1I * YR - U2 * XI 
         XR = BR(L,J) 
         XI = BI(L,J) 
         YR = BR(L1,J) 
         YI = BI(L1,J) 
         BR(L1,J) = U1 * YR - U1I * YI - U2 * XR 
         BR(L,J) = U1 * XR + U1I * XI + U2 * YR 
         BI(L,J) = U1 * XI - U1I * XR + U2 * YI 
         BI(L1,J) = U1 * YI + U1I * YR - U2 * XI 
  110 END DO 
!                                                                       
      LM1 = L 
      L = L1 
      GO TO 90 
!     ********** ITERATION STRATEGY **********                          
  120 IF (ITS .EQ. 50) GO TO 1000 
      IF (ITS .EQ. 10) GO TO 135 
!     ********** DETERMINE SHIFT **********                             
      B33 = BR(NA,NA) 
      B33I = BI(NA,NA) 
      IF (abs(CMPLX(B33,B33I)) .GE. EPSB) GO TO 122 
      B33 = EPSB 
      B33I = 0.0D0 
  122 B44 = BR(EN,EN) 
      B44I = BI(EN,EN) 
      IF (abs(CMPLX(B44,B44I)) .GE. EPSB) GO TO 124 
      B44 = EPSB 
      B44I = 0.0D0 
  124 B3344 = B33 * B44 - B33I * B44I 
      B3344I = B33 * B44I + B33I * B44 
      A33 = AR(NA,NA) * B44 - AI(NA,NA) * B44I 
      A33I = AR(NA,NA) * B44I + AI(NA,NA) * B44 
      A34 = AR(NA,EN) * B33 - AI(NA,EN) * B33I                          &
     &    - AR(NA,NA) * BR(NA,EN) + AI(NA,NA) * BI(NA,EN)               
      A34I = AR(NA,EN) * B33I + AI(NA,EN) * B33                         &
     &     - AR(NA,NA) * BI(NA,EN) - AI(NA,NA) * BR(NA,EN)              
      A43 = AR(EN,NA) * B44 
      A43I = AR(EN,NA) * B44I 
      A44 = AR(EN,EN) * B33 - AI(EN,EN) * B33I - AR(EN,NA) * BR(NA,EN) 
      A44I = AR(EN,EN) * B33I + AI(EN,EN) * B33 - AR(EN,NA) * BI(NA,EN) 
      SH = A44 
      SHI = A44I 
      XR = A34 * A43 - A34I * A43I 
      XI = A34 * A43I + A34I * A43 
      IF (XR .EQ. ZERO .AND. XI .EQ. ZERO) GO TO 140 
      YR = (A33 - SH) / 2.0D0 
      YI = (A33I - SHI) / 2.0D0 
      Z3 = sqrt(CMPLX(YR**2-YI**2+XR,2.0D0*YR*YI+XI)) 
      U1 = dble(Z3) 
      U1I = dimag(Z3) 
      IF (YR * U1 + YI * U1I .GE. ZERO) GO TO 125 
      U1 = -U1 
      U1I = -U1I 
  125 Z3 = (CMPLX(SH,SHI) - CMPLX(XR,XI) / CMPLX(YR+U1,YI+U1I))         &
     &   / CMPLX(B3344,B3344I)                                          
      SH = dble(Z3) 
      SHI = dimag(Z3) 
      GO TO 140 
!     ********** AD HOC SHIFT **********                                
  135 SH = AR(EN,NA) + AR(NA,ENM2) 
      SHI = 0.0D0 
!     ********** DETERMINE ZEROTH COLUMN OF A **********                
  140 A1 = AR(L,L) / B11 - SH 
      A1I = AI(L,L) / B11 - SHI 
      A2 = AR(L1,L) / B11 
      ITS = ITS + 1 
      IF (.NOT. MATZ) LOR1 = L 
!     ********** MAIN LOOP **********                                   
      DO 260 K = L, NA 
         K1 = K + 1 
         K2 = K + 2 
         KM1 = max(K-1,L) 
!     ********** ZERO A(K+1,K-1) **********                             
         IF (K .EQ. L) GO TO 170 
         A1 = AR(K,KM1) 
         A1I = AI(K,KM1) 
         A2 = AR(K1,KM1) 
  170    S = ABS(A1) + ABS(A1I) + ABS(A2) 
         U1 = A1 / S 
         U1I = A1I / S 
         U2 = A2 / S 
         R = SQRT(U1*U1+U1I*U1I+U2*U2) 
         U1 = U1 / R 
         U1I = U1I / R 
         U2 = U2 / R 
!                                                                       
         DO 180 J = KM1, ENORN 
            XR = AR(K,J) 
            XI = AI(K,J) 
            YR = AR(K1,J) 
            YI = AI(K1,J) 
            AR(K,J) = U1 * XR + U1I * XI + U2 * YR 
            AI(K,J) = U1 * XI - U1I * XR + U2 * YI 
            AR(K1,J) = U1 * YR - U1I * YI - U2 * XR 
            AI(K1,J) = U1 * YI + U1I * YR - U2 * XI 
            XR = BR(K,J) 
            XI = BI(K,J) 
            YR = BR(K1,J) 
            YI = BI(K1,J) 
            BR(K,J) = U1 * XR + U1I * XI + U2 * YR 
            BI(K,J) = U1 * XI - U1I * XR + U2 * YI 
            BR(K1,J) = U1 * YR - U1I * YI - U2 * XR 
            BI(K1,J) = U1 * YI + U1I * YR - U2 * XI 
  180    CONTINUE 
!                                                                       
         IF (K .EQ. L) GO TO 240 
         AI(K,KM1) = 0.0D0 
         AR(K1,KM1) = 0.0D0 
         AI(K1,KM1) = 0.0D0 
!     ********** ZERO B(K+1,K) **********                               
  240    S = ABS(BR(K1,K1)) + ABS(BI(K1,K1)) + ABS(BR(K1,K)) 
         U1 = BR(K1,K1) / S 
         U1I = BI(K1,K1) / S 
         U2 = BR(K1,K) / S 
         R = SQRT(U1*U1+U1I*U1I+U2*U2) 
         U1 = U1 / R 
         U1I = U1I / R 
         U2 = U2 / R 
         IF (K .EQ. NA) GO TO 245 
         XR = AR(K2,K1) 
         AR(K2,K1) = U1 * XR 
         AI(K2,K1) = -U1I * XR 
         AR(K2,K) = -U2 * XR 
!                                                                       
  245    DO 250 I = LOR1, K1 
            XR = AR(I,K1) 
            XI = AI(I,K1) 
            YR = AR(I,K) 
            YI = AI(I,K) 
            AR(I,K1) = U1 * XR + U1I * XI + U2 * YR 
            AI(I,K1) = U1 * XI - U1I * XR + U2 * YI 
            AR(I,K) = U1 * YR - U1I * YI - U2 * XR 
            AI(I,K) = U1 * YI + U1I * YR - U2 * XI 
            XR = BR(I,K1) 
            XI = BI(I,K1) 
            YR = BR(I,K) 
            YI = BI(I,K) 
            BR(I,K1) = U1 * XR + U1I * XI + U2 * YR 
            BI(I,K1) = U1 * XI - U1I * XR + U2 * YI 
            BR(I,K) = U1 * YR - U1I * YI - U2 * XR 
            BI(I,K) = U1 * YI + U1I * YR - U2 * XI 
  250    CONTINUE 
!                                                                       
         BI(K1,K1) = 0.0D0 
         BR(K1,K) = 0.0D0 
         BI(K1,K) = 0.0D0 
         IF ( MATZ ) Then
!                                                                       
         DO 255 I = 1, N 
            XR = ZR(I,K1) 
            XI = ZI(I,K1) 
            YR = ZR(I,K) 
            YI = ZI(I,K) 
            ZR(I,K1) = U1 * XR + U1I * XI + U2 * YR 
            ZI(I,K1) = U1 * XI - U1I * XR + U2 * YI 
            ZR(I,K) = U1 * YR - U1I * YI - U2 * XR 
            ZI(I,K) = U1 * YI + U1I * YR - U2 * XI 
  255    CONTINUE 

         End If
!                                                                       
  260 END DO 
!     ********** SET LAST A SUBDIAGONAL REAL AND END QZ STEP ********** 
      IF (AI(EN,NA) .EQ. ZERO) GO TO 70 
      R = abs(CMPLX(AR(EN,NA),AI(EN,NA))) 
      U1 = AR(EN,NA) / R 
      U1I = AI(EN,NA) / R 
      AR(EN,NA) = R 
      AI(EN,NA) = 0.0D0 
!                                                                       
      DO 270 J = EN, ENORN 
         XI = U1 * AI(EN,J) - U1I * AR(EN,J) 
         AR(EN,J) = U1 * AR(EN,J) + U1I * AI(EN,J) 
         AI(EN,J) = XI 
         XI = U1 * BI(EN,J) - U1I * BR(EN,J) 
         BR(EN,J) = U1 * BR(EN,J) + U1I * BI(EN,J) 
         BI(EN,J) = XI 
  270 END DO 
!                                                                       
      GO TO 70 
!     ********** SET ERROR -- BOTTOM SUBDIAGONAL ELEMENT HAS NOT        
!                BECOME NEGLIGIBLE AFTER 50 ITERATIONS **********       
 1000 IERR = EN 
!     ********** SAVE EPSB FOR USE BY CQZVEC **********                 
 1001 IF (N .GT. 1) BR(N,1) = EPSB 
      RETURN 
!     ********** LAST CARD OF CQZVAL **********                         
END SUBROUTINE DPCQZVAL_MMM
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
SUBROUTINE DPCQZVEC_MMM(NM,N,AR,AI,BR,BI,ALFR,ALFI,BETA,ZR,ZI) 
!                                                                       
      IMPLICIT NONE 
      INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN 
      real(DP) AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),ALFI(N),       &
     &       BETA(N),ZR(NM,N),ZI(NM,N)                                  
      real(DP) R,T,RI,TI,XI,ALMI,ALMR,BETM,EPSB 
      complex(DP) Z3 
      real(DP) ZERO 
!                                                                       
!     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FOURTH STEP OF THE   
!     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,  
!     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.      
!                                                                       
!     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES IN UPPER       
!     TRIANGULAR FORM, WHERE ONE OF THEM FURTHER MUST HAVE REAL DIAGONAL
!     ELEMENTS.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM 
!     AND TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
!     IT IS USUALLY PRECEDED BY  CQZHES  AND  CQZVAL.                   
!                                                                       
!     ON INPUT-                                                         
!                                                                       
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
!          DIMENSION STATEMENT,                                         
!                                                                       
!        N IS THE ORDER OF THE MATRICES,                                
!                                                                       
!        A=(AR,AI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,          
!                                                                       
!        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX WITH REAL 
!          DIAGONAL ELEMENTS.  IN ADDITION, LOCATION BR(N,1) CONTAINS   
!          THE TOLERANCE QUANTITY (EPSB) COMPUTED AND SAVED IN  CQZVAL, 
!                                                                       
!        ALFR, ALFI, AND BETA ARE VECTORS WITH COMPONENTS WHOSE         
!          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED              
!          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  CQZVAL,        
!                                                                       
!        Z=(ZR,ZI) CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE   
!          REDUCTIONS BY  CQZHES  AND  CQZVAL, IF PERFORMED.            
!          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE            
!          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.                 
!                                                                       
!     ON OUTPUT-                                                        
!                                                                       
!        A IS UNALTERED,                                                
!                                                                       
!        B HAS BEEN DESTROYED,                                          
!                                                                       
!        ALFR, ALFI, AND BETA ARE UNALTERED,                            
!                                                                       
!        Z CONTAINS THE EIGENVECTORS.  EACH EIGENVECTOR IS NORMALIZED   
!          SO THAT THE MODULUS OF ITS LARGEST COMPONENT IS 1.0 .        
!                                                                       
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ZERO = 0.0D0 
      IF (N .LE. 1) GO TO 1001 
      EPSB = BR(N,1) 
!     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********              
      DO 800 NN = 2, N 
         EN = N + 2 - NN 
         NA = EN - 1 
         ALMR = ALFR(EN) 
         ALMI = ALFI(EN) 
         BETM = BETA(EN) 
!     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********            
         DO 700 II = 1, NA 
            I = EN - II 
            R = 0.0D0 
            RI = 0.0D0 
            M = I + 1 
!                                                                       
            DO 610 J = M, EN 
               T = BETM * AR(I,J) - ALMR * BR(I,J) + ALMI * BI(I,J) 
               TI = BETM * AI(I,J) - ALMR * BI(I,J) - ALMI * BR(I,J) 
               IF (J .EQ. EN) GO TO 605 
               XI = T * BI(J,EN) + TI * BR(J,EN) 
               T = T * BR(J,EN) - TI * BI(J,EN) 
               TI = XI 
  605          R = R + T 
               RI = RI + TI 
  610       CONTINUE 
!                                                                       
            T = ALMR * BETA(I) - BETM * ALFR(I) 
            TI = ALMI * BETA(I) - BETM * ALFI(I) 
            IF (T .EQ. ZERO .AND. TI .EQ. ZERO) T = EPSB 
            Z3 = CMPLX(R,RI) / CMPLX(T,TI) 
            BR(I,EN) = dble(Z3) 
            BI(I,EN) = dimag(Z3) 
  700    CONTINUE 
!                                                                       
  800 END DO 
!     ********** END BACK SUBSTITUTION.                                 
!                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.               
!                FOR J=N STEP -1 UNTIL 2 DO -- **********               
      DO 880 JJ = 2, N 
         J = N + 2 - JJ 
         M = J - 1 
!                                                                       
         DO 880 I = 1, N 
!                                                                       
            DO 860 K = 1, M 
               ZR(I,J) = ZR(I,J) + ZR(I,K) * BR(K,J) - ZI(I,K) * BI(K,J) 
               ZI(I,J) = ZI(I,J) + ZR(I,K) * BI(K,J) + ZI(I,K) * BR(K,J) 
  860       CONTINUE 
!                                                                       
  880 CONTINUE 
!     ********** NORMALIZE SO THAT MODULUS OF LARGEST                   
!                COMPONENT OF EACH VECTOR IS 1 **********               
      DO 950 J = 1, N 
         T = 0.0D0 
!                                                                       
         DO 930 I = 1, N 
            R = abs(CMPLX(ZR(I,J),ZI(I,J))) 
            IF (R .GT. T) T = R 
  930    CONTINUE 
!                                                                       
         DO 940 I = 1, N 
            ZR(I,J) = ZR(I,J) / T 
            ZI(I,J) = ZI(I,J) / T 
  940    CONTINUE 
!                                                                       
  950 END DO 
!                                                                       
 1001 RETURN 
!     ********** LAST CARD OF CQZVEC **********                         
END SUBROUTINE DPCQZVEC_MMM
!}}}

end module w20mod
