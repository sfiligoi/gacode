! 02feb2001 jek fixed section w/ xparam(24) for power law ExB
! 08may2000 jek merged ITG and ETG loops over ky for better MPI optimization
! 23nov1999 jek added zgeev eigenvalue solver (eigen_gf=2)
! 29mar1999 fgtok -s cgg.table "dmc:  rename eispack routines"
! 29mar1999 fgtok -s rr.table "dmc:  rename intrinsic REAL -> REAL"
!
! dmc -- Cray/workstation portable real*8<-->complex*16 conversion routines
 
!gms#include "f77_dcomplx.h"
 
!
!glf2d.f 02-feb-01 Kinsey
!
! GACDODE version Sep. 25, 2014 G.M. Staebler
! f90 version Oct. 10,2014 G.M. Staebler
!
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!
 
          subroutine glf23_TM
 
! version 3.3
!
!***********************************************************************
! questions  should be addressed to
!  Ron Waltz 619-455-4584  or email: waltz@gav.gat.com
!***********************************************************************
 
 
! 2D GLF equations with massless isothermal passing electrons from
!  Waltz et al, Phys. of Plasmas 6(1995)2408
 
! In Eq. 24  p_u_par is replaced by n_u/(1-reps) +t_u and
! the isothermal conditions t_u= (betae/2)*w_s/k_par*(rlt)*a_par is
! used. Thus n_u = (1-reps)*ph (adiabatic passing electrons) for betae=0
!
! In Eq. 23 (p_u_par+p_u_per) is replaced by n_u+t_u
! using the isothermal condition the mHD beta limit is too low by
! 1/(1+ reps(-0.25+0.75rlte)/(taui*(rln+rlti)+(rln+rlte))
! It is possible to patch this up by replacing with x*n_u+y*t_u
! then solving for x and y to obtain
! universal MHD betae_crit=2*k_par**2/(w_d*(taui*(rln+rlti)+w_d*(rln+rlti)))
! beta_crit=(1+taui)*betae_crit=(1/2)(1/q**2)*(L_p/rmajor)
! 1/2 is replaced by s_hat in models with shear
 
! EVERYTHING else including normalizing units follows the paper.
 
!  unit of microlength is rho_s, unit of macrolength is a
!   a is typically rho value at separatrix
!  unit of time is a/c_s; unit of diffusion is (c_s/a)*rho_s**2
!  c_s=sqrt(Te/M_i), omega=eB/cM_i, rho_s=c_s/omega
 
! example balance equations to clarify meaning  of diffusivities
!
!       chie_hat is effective energy diffusivity
!
!  (3/2) d ne Te/ dt =
! -1/V(rho)_prime d /d rho V(rho)_prime
!             |grad rho|**2 ne chie_hat (c_s rho_s**2/a)(-d Te/ d rho)
!      -exch_hat (ne Te) c_s/a (rho_s/a)**2 +heating density
!
! and similarly for ion equation
! note that no convective part is added, ie "convection" is included
! inside chie_hat
! note no impurity energy flow is computed although it could be easily done
 
!        d_hat is effective plasma diffusivity for ions
!
!        d ni / dt =
! -1/V(rho)_prime d /d rho V(rho)_prime
!               |grad rho|**2 ne d_hat (c_s rho_s**2/a) (-d ni/ d rho)
!        + plasma source density
 
!        d_im_hat is the effective plasma diffusivity for impurities
 
!        eta_phi is the toroidal viscosity or momentum diffusivity
!
!       M ni d v_phi/ dt =
! error found 5/12/98  should be d (M ni v_phi)/ dt =
! -1/V(rho)_prime d /d rho V(rho)_prime
!      |grad rho|**2 ( M ni eta_phi_hat (c_s rho_s**2/a) (-d v_phi/ d rho)
!                    +    M v_phi ne d_hat (c_s rho_s**2/a) (-d ne/ d rho))
!        + toroidal momentum source density
!
! note that a convective part was added
!
!  eta_par_hat and eta_per_hat are diagnostic. See CAUTION on eta_phi_hat
!  at large gamma_p=  (-d v_phi/ d rho) /(c_s/a)
!
!  chie_e_gf is the eta_e mode electron transport which is te <-> ti
!  and mi <-> me isomorphic to eta_i (ITG) ion transport
!  with adiabatic electrons.
!  these mode obtain at high-k where the ions are adiabatic from
!  the gyro cut-off.
!  their wave numbers are sqrt(mi/me) larger than ITG modes and
!  since their frequencies are sqrt(mi/me) larger, they are not
!  rotationally shaer satbilized.
!  when xparam_gf(10).eq.0 xparam_gf(10)*chie_e_gf is added to
!  chie_gf and chie_e_gf is a diagnostic.
 
! input

!  eigen_gf = 0 use cgg eigenvalue solver (default)
!           = 1 use generalized tomsqz eigenvalue solver
!           = 2 use zgeev eigenvalue solver
!  nroot number of equations
!  iflagin(1:20) control flags
!   iflagin(1) 0 use ky=ky0; 1 use landau damping point
!   iflagin(2) 0. local w_d and k_par "2d"; 1 fit to trial function "3d"
!   iflagin(3) 0,1,and 2 fix up park low high beta and beta lim elong factor
!   iflagin(4) 0 trapped electron Waltz EoS 1 weiland EoS
!   iflagin(5) rms_theta 0:fixed; 1 inverse to q/2 ; 2 inverse to root q/2
!                        3: inverse to xparam(13)*(q/2-1)+1.
!  xparam(1:20) control parameters
!   xparam(1:2): idelta=xi*xparam(1)+xparam(2) nonadiabatic electron response
!   xparam(3) multiplier park_gf(high betae)/ park_gf(low betae) -1
!   xparam(6)+1. is enhancement of xnueff
!   xparam(7) coef of resistivity
!   xparam(8) cut off on rotational stabilization
!   xparam(9)+1. is shape (triangularity) enhancement to beta_crit
!   xparam(10) is high k electron mode enhancement
!   xparam(11:12) lamda parameters
!   xparam(13) rms_theta q-dependence
!   xparam(14)  adjustment to gamma_p avoiding negative viscosity
!   xparam(15)   (1+xparam(15)*reps trapped electron fraction
!   xparam(16) rms_theta shat dependence
!   xparam(17) ""
!   xparam(18) rms_theta betae dependence
!   xparam(19:20)  extra
!   xparam(21) 1 add impurity energy diffusivity to ion energy diffusivity
!   xparam(22) >0 keeps gamma_e from changeing spectrum
!   xparam(23) 1. kills kx**2 in k_m**2
!   xparam(24) exb damping model
!  ky0=k_theta*rho_s; k_theta= nq/r; normally 0.3
!  rms_theta width of phi**2 mode function for best fit near pi/3
!  rlti=a/L_Ti   a/L_f= sqrt(kappa) a d ln f / d rho
!  rlte=a/L_Te
!  rlne= a/L_ne
!  rlni= a/L_ni
!  rlnimp= a/L_nim
!  dil=1.-ni_0/ne_0  dilution
!  apwt = ni_0/ne_0
!  aiwt = nim_0/ne_0
!  taui=Ti/Te
!  rmin=r/a
!  rmaj=Rmaj/a
!  xnu=nu_ei/(c_s/a)
!  betae=neTe/(B**2/(8pi))  0 is electrostatic
!  shat= dlnr/drho used only for parallel dynamics part
!  alpha local shear parameter or MHD pressure grad (s-alpha diagram)
!  elong= local elongation or kappa
!  xwell amount of magnetic well xwell*min(alpha,1)
!  park=1  (0) is a control parameter to turn on (off) parallel motion
!       0.405 best at zero beta and 2.5x larger at high beta..see iflagin(3)
!  ghat=1  (0) is a control parameter to turn on (off) curvature drift
!  gchat=1 (0) is a control parameter to turn on (off) div EXB motion
!  adamp= radial mode damping exponent  1/4 < adamp < 3/4
!       0.25 from direct fit of simulations varying radial mode damping
!   but 0.75 is better fit to rlti dependence
!  alpha_star O(1-3)  gyyrobohm breaking coef for diamg. rot. shear
!  gamma_star ion diamagnetic rot shear rate in units of c_s/a
!  alpha_e O(1-3)   doppler rot shear coef
!  gamma_e    doppler rot shear rate in units of c_s/a
!  alpha_p 1.5  fit for parallel velocity shear effect at rmaj=3 and q=2
!  gamma_p    parallel velocity shear rate (-d v_phi/ drho) in units of c_s/a
!  kdamp model damping normally 0.
 
! output
 
!  yparam(20) output diagnostics
! kyf  value of ky used
! gamma   leading mode growth rate in c_s/a
! freq    leading mode freq rate in c_s/a
! ph_m    (e phi /T_e)/(rho_s/a)  saturation value
! d_hat    plasma diffusivity for ions
! d_im_hat    plasma diffusivity for impurities
! chii_hat ion energy diffusivity
! chie_hat electron energy diffusivity
! exch_hat anomalous e to i energy exchange
! eta_par_hat parallel component of toroidal momentum diffusivity
! eta_per_hat perpendicular    ""
! eta_phi_hat toroidal momentun diffusivity
 
! internal definitions
! nroot = number of equations,
!   nroot=12 full impurity dynamics
!   nroot=9 exb convective impurity dynamics
!   nroot=8 full pure plasma, nrout=6 (betae=0), nrout=5 (betae=0 and park=0)
! v(i)  12 solution vector
!   v(1)=n_i,v(2)=p_par,v(3)=p_per,v(4)=n_t,v(5)=p_t
!   v(6)=u_par, v(7)=n_u, v(8)=a_par
!   v(9)=n_im, v(10)=p_im_par,v(11)=p_im_per,v(12)=u_im_par
! -i*omega v(i)= sum_j amat(i,j) v(j) where omega=freq+xi*gamma
! quasineitrality is
!  (-idelta+(1-dil)*(1/taui)*(1-g1))*ph=f0*ph=(1-dil)*n_i-n_t-n_u
!  or (1-idelta-reps+(1-dil)*(1/taui)*(1-g1))*ph=f1*ph=(1-dil)*n_i-n_t
!  for betae=0
! k_m  inverse mixing length
 
! numerical controls and flags
!
!
!...Dimensions
!
! neq maximum number of equations
! ieq actual number used
 
!***********************************************************************
!***********************************************************************
!
      USE glf23_gf
      USE glf23_cgg
      implicit none
!      include 'mpif.h'
!gms      include '../inc/glf.m'
!
! Glf is common block, which must contain all the _gf inputs and outputs
!
!      character*1 jobvr, jobvl
      integer :: nroot,iflagin(30), ilhmax, ilh, ikymaxtot
      integer :: lprint, ieq, j1, j2, j, i, jmax
      integer :: ifail, jroot(4), itheta, iky, iky0, iroot
      integer,parameter :: neq = 12
      real,parameter :: epsilon = 1.E-34 
!
      real::  pi, xparam(30),yparam(2*nmode)
      real :: ky0,rms_theta,rlti,rlte,rlne,rlni,dil,taui
      real :: rmin,rmaj,q,rlnimp,amassimp,zimp,mimp
      real :: aikymax, aiky, apwt, aiwt
      real :: alpha_mode, gamma_mode, alpha_p, gamma_p
      real :: amassgas, chi_par_1, chi_per_1, x43
      real::  anorm, ave_g, ave_g0
      real :: ave_cos, ave_theta2, ave_kxdky2
      real :: dtheta, theta, phi2, ave_k_par, chk
      real :: alpha_n, yk, byk, fnn, fnp, fnf, fpn, fpp, fpf
      real :: xnu_col, amass_e, chkf, xadiabat
      real :: eta_par_hat, eta_per_hat, anorm_k, del_k
      real :: gamma_k_max, gamma_gross_net
      real :: xnu,betae,shat,alpha,elong,xwell
      real :: park,ghat,gchat,kdamp
      real :: adamp,alpha_star,gamma_star,alpha_e,gamma_e
      real :: kyf,gamma,freq,ph_m,d_hat,d_im_hat
      real :: chii_hat,chi_im_hat,chie_hat,exch_hat
      complex ::  xi, idelta
      complex ::  v(1:12), amat(1:12,1:12)
      complex :: n_i,p_par,p_per,n_t,p_t,u_par,n_u,a_par,ph,t_u,n_e
      complex ::  n_im,p_im_par,p_im_per
!     complex u_im_par
      real :: b0,g0,g1,g2,g3,g12,g23
      real ::   b0i,g0i,g1i,g2i,g3i,g12i,g23i
      complex :: f0,f1
      real :: k_par,ky,kx,k_per,k_m
      real :: w_s, w_d, w_d0, w_cd
      real :: reps,xnueff,betae0,k_par0
      complex :: xmu,lamda_d
      complex :: xnu_par_par,xnu_par_per,xnu_per_par,xnu_per_per
      real :: gam_par,gam_per,x_par,x_per,xt_mhd,yt_mhd
      real :: th,tc,fh,fc
      real :: phi_norm,gamma_r
      complex :: chknu,chknt,chknt2
      real :: phi_renorm,gamma_net
!
!...Declarations for eigenvaluesolver
!
      real :: zgamax
!
!... solver varaibles
!
      integer,parameter :: iar=neq
!     integer iai, ivr, ivi, intger(neq) ! if NAG solver f02ake used
!     parameter ( iai=neq, ivr=neq, ivi=neq )
      real :: ar(iar,neq), ai(iar,neq), rr(neq), ri(neq)
      real :: vr(iar,neq), vi(iar,neq)
!      real*8 br(iar,neq), bi(iar,neq), betat_tom(neq),ztemp1
 
      integer :: matz
      real :: fv1(neq),fv2(neq),fv3(neq)
!
! amat(i,j) = complex matrix A
! zevec(j) = complex eigenvector
!
!      integer lwork
!      parameter ( lwork=198 )
!      complex*16 mata(iar,neq),cvr(iar,neq),cvl(iar,neq),w(neq)
!      complex*16 work(lwork)
!      double precision rwork(2*neq)
      complex :: zevec(neq,neq), zomega(neq)
      real :: gammaroot(4),freqroot(4),phi_normroot(4)
!
!      real*8 anorm_glob,diff_glob,diff_im_glob,chii_glob
!     &       ,chie_glob, exch_glob,eta_par_glob,eta_per_glob
!     &       ,eta_phi_glob,chie_e_glob
!      real*8 chie_e_k_gf_glob(20), chie_k_gf_glob(20),
!     &       chie_tot_gf(20), xkyf_k_gf_glob(20)
!      real*8 m(2), m_global(2)
!
!     external f02ake
!
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!
!
!   return zero if no unstable modes
!
!gms      if(ipert_gf.eq.1.and.ngrow_k_gf(0).eq.0)go to 888
!
!...initialize variables
!
      diff_gf=0.0
      diff_im_gf=0.0
      chii_gf=0.0
      chi_im_gf=0.0
      chie_gf=0.0
      exch_gf=0.0
      eta_par_gf=0.0
      eta_per_gf=0.0
      eta_phi_gf=0.0
      chie_e_gf=0.0
      do j1=1,4
        gamma_gf(j1)=0.0
        freq_gf(j1)=0.0
        xky_gf(j1)=0.0
      enddo
!
      iflagin_gf(1)=0
      iflagin_gf(2)=1
      iflagin_gf(3)=1
      iflagin_gf(4)=0
      iflagin_gf(5)=3
!
      xparam_gf(1)=0.0
      xparam_gf(2)=0
      xparam_gf(3)=.70
      xparam_gf(4)=0.0
      xparam_gf(6)=0.0
!      xparam_gf(6)=-0.250
      xparam_gf(7)=1.0
      xparam_gf(8)=0.0
      xparam_gf(9)=1.0
      xparam_gf(10)=0.0
      xparam_gf(11)=0.0
      xparam_gf(12)=0.0
      xparam_gf(13)=0.20
      xparam_gf(14)=1.0
      xparam_gf(15)=-0.10
      xparam_gf(16)=0.150
      xparam_gf(17)=0.10
      xparam_gf(18)=.00
      xparam_gf(19)=0.0
      xparam_gf(20)=0.0
      xparam_gf(21)=0.0
      if(ns_gf.eq.3)xparam_gf(21)=1.0  ! turn of impurity energy flux
      xparam_gf(22)=0.0
      xparam_gf(23)=1.0
      xparam_gf(24)=0.0
      xparam_gf(25)=0.0
      xparam_gf(26)=0.0
      xparam_gf(27)=0.0
      xparam_gf(28)=0.0
      xparam_gf(29)=0.0
      xparam_gf(30)=0.0
!
!      rms_theta_gf=zpi/3.0
      park_gf  =0.70
      ghat_gf  =1.0
      gchat_gf =1.0
!
      adamp_gf=0.50
      alpha_star_gf=0.0
      alpha_mode_gf=0.0
!      gamma_e_gf =-.0000000000010
      xkdamp_gf=0.0
      alpha_p_gf=0.50
!
! betae and collisionality multipliers
!
!      cbetae=1.E-6
!     cbetae=1.0 !  full electromagetic
!      cxnu=1.0   !  full collisionality
!
      cnorm_gf=50.0   ! normalization factor
      cnorm_p_gf=50.0 ! normalization factor for Di
!
! non glf23 parameter
!
!       cmodel=1.0
!       xalpha=x_alpha
 
! for original model
       if(version_gf.eq.1) then
         cnorm_gf=100.0
         cnorm_p_gf=100.0
         iflagin_gf(5)=3
         park_gf=0.7
         adamp_gf=0.50
         alpha_p_gf=0.500
         xparam_gf(10)=1.0
         xparam_gf(13)=0.20
         xparam_gf(15)=-0.1
         xparam_gf(16)=0.0
         xparam_gf(17)=0.10
         xparam_gf(19)=0.0
       endif
!
!
!... parameters for revised GLF23 models
!
        if (version_gf.eq.2) then   ! retuned model v1.61
          cnorm_gf=50.0         ! ITG/TEM normalization
          xparam_gf(10)=12.0    ! ETG normalization (cnorm*xparam(10))
          xparam_gf(13)=0.15      ! rms_theta q-dependence
          xparam_gf(15)=-0.1     ! trapped particle fraction reduction
          xparam_gf(16)=0.15     ! rms_theta shat dependence
          xparam_gf(17)=0.25     ! rms_theta shat dependence
          xparam_gf(19)=1.0      ! rms_theta alpha-dependence
          iflagin_gf(5)=5        ! rms_theta w/ tau**0.25
          park_gf=0.8
          adamp_gf=0.7
          alpha_p_gf=0.35
          alpha_e_gf=1.35
!          bt_flag=1
        endif
!
        if (version_gf.eq.3) then      ! renorm + real geometry fit
          cnorm_gf=27.0         ! ITG/TEM normalization
          xparam_gf(10)=17.80   ! ETG normalization
          xparam_gf(15)=0.100   ! trapped ptcle fraction enhancement
          xparam_gf(6)=-0.250   ! xnueff enhancement
          xparam_gf(13)=0.950
          xparam_gf(16)=0.200
!          bt_flag=2              ! real geometry in diffusion
        endif
!        if (iglf.eq.99) then      ! renormalization only
!          xparam_gf(10)=0
!          xky0_gf=0.30
!          ikymax_gf=1
!        endif
!        cnorm_gf=27.0         ! ITG/TEM normalization
!        xparam_gf(10)=17.80    ! ETG normalization
!        xpparam_gf(10)=0.0   ! ETG normalization
!        xparam_gf(15)=0.150   ! trapped ptcle fraction enhancement
!
!    relaxation can be turned on for one call per step
!       relx=0.0   ! relaxation off
!

! inputs.........................................................

      do i=1,30
       iflagin(i)=iflagin_gf(i)
       xparam(i)=xparam_gf(i)
      enddo
!      ilhmax=1
!      ikymaxtot=ikymax_gf
!     if (xparam_gf(10).gt.0.) ilhmax=2
!
! If ETG modes included, then double ky spectrum
! Inside ky loop, ilh=1 low k ion modes and ilh=2 high k electron modes
! For ETG mode, use complete te <-> ti and mi <-> me isomorphism
! with adiabatic electron ITG mode then chii_hat will be electron
! transport in gyrobohm electron units with T_i.
! chie_e_gf converted back to c_s*rho_s**2/a units and added
! to chie_gf after ky loop
!
!  transfer values from inputs
!
! ky spectrum parameters
!
      if(use_transport_model_gf)then 
        ikymax_gf=10
        xkymin_gf=.020
        xkymax_gf=.50
      else !single ky linear stability
        ikymax_gf=1
        xkymin_gf=xky0_gf
        xkymax_gf=xky0_gf       
      endif
!
      if(use_adiabatic_electrons_gf)then
         nroot_gf=6
         ns_gf=2
         apwt_gf = 1.0
         aiwt_gf = 0.0
         betae_gf = 1.0E-12 ! turn off EM 
         xnu_gf = 0.0       !  turn off electron collisions
         rlnimp_gf = epsilon
         rmin_gf = epsilon  ! remove trapped electrons
         ilhmax=1
         ikymaxtot=ikymax_gf
      else  !kinetic electrons
         nroot_gf=8
         if(ns_gf.gt.2)nroot_gf=12 !for full impurity dynamics
         ilhmax = 2
         ikymaxtot = 2*ikymax_gf
         if (xparam_gf(10).eq.0.0) then
           ilhmax=1
           ikymaxtot=ikymax_gf
         endif
      endif
      nroot=nroot_gf
      ky0=xky0_gf
      rms_theta=rms_theta_gf
      rlti=rlti_gf
      rlte=rlte_gf
      rlne=rlne_gf
      rlni=rlni_gf
      rlnimp=rlnimp_gf
      dil=dil_gf
      apwt=apwt_gf
      aiwt=aiwt_gf
      taui=taui_gf
      rmin=rmin_gf
      rmaj=rmaj_gf
      q=q_gf
      xnu=xnu_gf
      betae=betae_gf
      shat=shat_gf
      alpha=alpha_gf
      elong=elong_gf
      xwell=xwell_gf
      park=park_gf
      ghat=ghat_gf
      gchat=gchat_gf
      adamp=adamp_gf
      alpha_star=alpha_star_gf
      gamma_star=gamma_star_gf
      alpha_e=alpha_e_gf*alpha_e_mult_gf
      gamma_e=gamma_e_gf
      alpha_mode=alpha_mode_gf
      gamma_mode=gamma_mode_gf
      alpha_p=alpha_p_gf*alpha_p_mult_gf
      gamma_p=ABS(gamma_p_gf)
      kdamp=xkdamp_gf
      lprint=lprint_gf
      amassgas=amassgas_gf
      amassimp=amassimp_gf
      zimp=zimp_gf
!gms      if(ipert_gf.eq.0)then
       do j=0,nmode
        ngrow_k_gf(j) = 0
       enddo
!gms      endif
!
      idelta=0.0
!     if(ilh.eq.1) idelta=xi*xparam(1)+xparam(2)
 
!.................................................................
!
      if (lprint.gt.0)then
         open(unit=1,file='./out.glf23.debug',status='replace')
      endif
      ieq  = nroot
!
      if (lprint.eq.99) then
      write(1,*) 'ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q: ', &
        ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q
      write(1,*)'xnu,beta,shat,alpha,elong,xwell: ', &
        xnu,betae,shat,alpha,elong,xwell
      write(1,*)'park, ghat, gchat: ', &
        park, ghat, gchat
      write(1,*)'adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp: ', &
        adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp
 
      endif
      if (lprint.eq.98) then
      write(2,*) 'ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q: ', &
        ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q
        write(2,*)'xnu,betae,shat,alpha,elong,xwell: ', &
        xnu,betae,shat,alpha,elong,xwell
        write(2,*)'park, ghat, gchat: ', &
        park, ghat, gchat
        write(2,*)'adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp: ', &
        adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp
 
      endif
!
      xi=(0.0,1.0)
      pi=atan2 ( 0.00, -1.00 )
 
! GLF model coefficients
 
      chi_par_1=2.0*sqrt(2.0/pi)
      chi_per_1=sqrt(2.0/pi)
 
      gam_par=3.0
      gam_per=1.0
      x_par=2.0
      x_per=3.0/2.0
 
      xmu=(0.800+.570*xi)
 
      xnu_par_par=(1.0+xi)
      xnu_par_per=0.0
      xnu_per_par=0.0
      xnu_per_per=(1.0+xi)
 
      lamda_d=(-0.70-0.800*xi)
      if(xparam(11).ne.0..or.xparam(12).ne.0.)then
        lamda_d=xparam(11)-xi*xparam(12)
      endif
      x43=1.0
      if (iflagin(4).eq.1) then
       lamda_d=5.0/3.0
       x43=4.0/3.0
      endif
!
! 3d trial wave function analysis
!
      if(iflagin(2).ge.1) then
       if(rms_theta.eq.0.) rms_theta=pi/3.0
       if(iflagin(5).eq.1) rms_theta=rms_theta_gf*(2.0/q_gf)
       if(iflagin(5).eq.2) rms_theta=rms_theta_gf*(2.0/q_gf)**0.50
       if(iflagin(5).eq.3) then
         rms_theta=rms_theta_gf/  &
           (xparam(13)*(q_gf/2.0-1.0)+1.0)  &
          /sqrt(1.0+xparam(16)*(shat_gf**2-1.0)+xparam(17)*  &
          (shat_gf-1.0)**2)   &
          /(1.0+xparam(18)*sqrt(betae/.0060))
       endif
       if(iflagin(5).eq.4)then
         rms_theta=rms_theta_gf/  &
           (xparam(13)*(q_gf/2.0-1.0)+1.0)  &
       /sqrt(1.0+xparam(16)*(shat_gf**2-1.00)+xparam(17)*  &
         (shat_gf-1.0)**2+xparam(19)*(alpha-0.50)**2.0)  &
       /(1.0+xparam(18)*sqrt(betae/.0060))
       endif
       if(iflagin(5).eq.5)then
          rms_theta=rms_theta_gf/  &
          (xparam(13)*((q_gf/2.0)-1.0)+1.0)  &
         /sqrt(1.0+xparam(16)*((shat_gf-  &
         xparam(19)*alpha)**2-0.50)+xparam(17)*  &
         (shat_gf-xparam(19)*alpha-0.50)**2)  &
         /(1.0+xparam(18)*sqrt(betae/.0060))
       endif
       if(iflagin(5).eq.6)then
         rms_theta=rms_theta_gf/  &
         (xparam(13)*((q_gf/2.0)-1.0)+1.0)  &
         /sqrt(1.0+xparam(16)*((shat_gf-  &
         xparam(19)*alpha)**2-0.50)+xparam(17)*  &
         (shat_gf-xparam(19)*alpha-0.50)**2)/  &
         taui_gf**0.250
       endif
       if(lprint.eq.3) then
         write(*,'(1x,a12,1p3e13.5)')'rms_theta = ', rms_theta
       endif
! along the field line physics with wave function
! phi=exp(-theta**2/(4.*rms_theta**2))=W_even
! ave_F= [int(0 to inf) F phi**2 d_theta]/[int(0 to inf)  phi**2 d_theta]
! ave_theta**2=(rms_theta)**2
!
! phi, densities and pressures are even functions of theta
! the odd functions like u_par can be represented by
! W_odd= W*i*theta/rms_theta*W_even
! then for W=-1, the k_par operator = i*k_par=1/(rmaj*q) d/dtheta
! becomes 1/(rmaj*q)/(2.*rms_theta)*park  (park=1) in every equation
! ie ave_k_par=1/(2.*rms_theta)
! park is tuned to best fit.
!
! parallel velocity shear gamma_p breaks parity so wave functions
! become mixed W_even->(1-xi*gamma_p*alpha_n*theta/rms_theta)*W_even
!
! gamma_p*alpha_n mustbe dimensionless and independent of norm a
! hence alpha_n=(rmaj/3)*alpha_p since gamma_p is in units
! of c_s/a  and rmaj=rmajor/a. Since in the slab limit where
! where we can have the parallel shear drive, rmaj enters only with
! parameter rmaj*q, we further assume
! alpha_n=(rmaj/3.)*(q/2)*alpha_p as the appropriate scaling
! q=2 and rmaj=3 are the norm points for alpha_p=1.5
! For the extreme toroidal limit q-> infinity where rmaj and
! q are not product associated, we will lose the instability.
!
! to first order in gamma_p*alpha_n
! this leads to a weighting  factor [gamma_p*alpha_n] in the
! xi*ky*ph1*gamma_p linear drive and in the
! eta_phi_hat=conjg(u_par)*(-xi*ky*ph)/gamma_p toroidal vocosity.
!
! the correct dependence gamma-gamma0 going like gamma_p**2 is found
! but QLT eta_phi_hat goes like gamma_p**2 also
! CAUTION: this is worked out only to first order in gamma_p*alpha_n
! small.  It seems likely that there is a higher order saturation factor
! something like 1/(1+(gamma_p*alpha_n)**2) in  eta_phi_hat
!
! doppler (EXB) rotational shear also breaks parity
! thus there should a term egamma*alpha_n_e added to gamma_p*alpha_n
! but it is unclear how to weight alpha_n_e compared to alpha_n
!
! see R.R.Dominguez and G.M Staebler Phys. Fluids B5 (1993) 3876
! for a discussion of QLT theory of anomalous momentum transport in
! slab geometry
!
! compute weight factors
! fix later so these are computed only once per j grid point
!
!      ave_theta2=rms_theta**2
!
      anorm=0.0
      ave_g=0.0
      ave_g0=0.0
      ave_cos=0.0
      ave_theta2=0.0
      ave_kxdky2=0.0
!
      dtheta=4.0*rms_theta/100.0
      theta=0.0
!
      do itheta=1,100
       theta=theta+dtheta
       phi2=exp(-theta**2/(2.0*rms_theta**2))
       anorm=anorm+phi2*dtheta
       ave_theta2=ave_theta2+ theta**2*phi2*dtheta
       ave_g=ave_g + (-xwell*min(1.0,alpha)+cos(theta)+  &
        (shat*theta-alpha*sin(theta))*sin(theta))*phi2*dtheta
       ave_g0=ave_g0 + phi2*dtheta
       ave_kxdky2=ave_kxdky2+  &
        (abs(shat*theta-alpha*sin(theta)))**2*phi2*dtheta
       ave_cos=ave_cos +  &
        cos(theta)*phi2*dtheta
      enddo
!
      ave_theta2=ave_theta2/anorm
      ave_g=ave_g/anorm
      ave_g0=ave_g0/anorm
      ave_kxdky2=ave_kxdky2/anorm
      ave_cos=ave_cos/anorm
!
      ave_k_par=1/(2.0*rms_theta)
 
      chk=abs(ave_theta2-rms_theta**2)/rms_theta**2
      if (chk.gt..02) write (6,*) 'chk:', chk
 
      alpha_n=(rmaj/3.0)*(q/2.0)*alpha_p
 
      if(lprint.eq.2) then
       write(6,*) 'rms_theta,chk :', rms_theta, chk
       write(6,*) 'ave_theta2,ave_g,ave_k_par,ave_cos:',  &
         ave_theta2,ave_g,ave_k_par,ave_cos
      endif
      endif
      if(iflagin(2).eq.0) then
       shat=1.0
       ave_theta2=1.0
       ave_g=1.0
       ave_g0=1.0
       ave_kxdky2=1.0
       ave_k_par=1.0
       ave_cos=1.0
      endif
!
! start ky loop
! first half ITG, second half high-k ETG ... each with ikymax_gf modes
! ilh=1 for low k ion modes, ilh=2 high k electron modes
!
      do iky0=1,ikymaxtot
!
!gms      iky=iky0
      iky = ikymax_gf+1-iky0
      ilh=1
!
! offset iky if in high-k range and set ilh=2
!
      if (iky0.gt.ikymax_gf) then
!gms         iky=iky0-ikymax_gf
         iky = ikymaxtot+1-iky0
         ilh=2
      endif
!
      if (ilh.eq.2) then
       nroot=6
       ieq=nroot
       xnu=0.0
       betae=1.E-6
       rlte=rlti_gf
       rlti=rlte_gf
       rlne=rlni_gf
       rlni=rlne_gf
       rlnimp=epsilon
       dil=1.0-1.0/(1.0-dil_gf)
       apwt=1.0
       aiwt=0.0
       taui=1.0/taui_gf
       rmin=epsilon
       xparam(7)=0.0
       xparam(6)=-1.0
       alpha_star=0.0
       alpha_e=0.0
       alpha_p=0.0
       alpha_n=0.0
       alpha_mode=0.0
!       alpha_gf=0.0
! check this for current driven mode
      endif
!
      idelta=0.0
      if (ilh.eq.1) idelta=xi*xparam(1)+xparam(2)
!
! logarithmic ky grid
!
      if(ikymax_gf.gt.1) then
       aikymax=ikymax_gf
       aiky=iky
       yk=aiky/aikymax
       byk=log(xkymax_gf/xkymin_gf)/(1.0-1.0/aikymax)
       ky=xkymax_gf*exp(byk*(yk-1.0))
      endif
      if(ikymax_gf.eq.1) then
!     ky=sqrt(2.*taui)/rlti/(rmaj*q)
!     from w_star_ti=v_i_th*k_par
!  possible physics basis of q (ie current) scaling ..to be determined
       if(iflagin(1).eq.0) ky=ky0
       if(iflagin(1).eq.1)then
         ky=ky0*sqrt(taui)*(3.0/rlti)*(3.0*2.0/rmaj/q)
       endif
       if(iflagin(1).eq.2) ky=ky0*(2.0/q)
       if(ky0.eq.0.) ky=0.30
      endif
!
!
      kyf=ky
!
      kx=ky*sqrt(ave_kxdky2)
      k_per=sqrt(ky**2+kx**2)
      k_m=sqrt(ky**2+(1.0-xparam(23))*kx**2) ! inverse mixing length model
!
       do iroot=1,4
        gammaroot(iroot)=0.0
        freqroot(iroot)=0.0
        phi_normroot(iroot)=0.0
       enddo
        d_hat=0.0
        d_im_hat=0.0
        chie_hat=0.0
        chii_hat=0.0
        chi_im_hat=0.0
        exch_hat=0.0
        eta_par_hat=0.0
        eta_per_hat=0.0
        jroot(1)=0
        jroot(2)=0
        jroot(3)=0
!
! skip this k for perturbation if last call was stable
!
!gms      if(ipert_gf.eq.1.and.ngrow_k_gf(iky0).eq.0)go to 777
! skip this k if the previous k was stable and 4 k's have been done
!gms      if(iky.lt.ikymax_gf-4.and.ngrow_k_gf(iky0-1).eq.0)go to 777
!
! primary ions
!
      b0=taui*k_per**2
      b0=(1.0+xparam_gf(20))*b0
      if(lprint.eq.2) write(*,'(1x,a5,1p3e13.5)') 'b0 = ',b0
!
!     Pade aproximates...may use gamma functions later
      g0=1.0
      g1=1.0/(1+b0)
      g2=1.0/(1+b0)*g1
      g3=1.0/(1+b0)*g2
      g12=(g1+g2)/2.0
      g23=(g2+g3)/2.0
!
! impurity ions
!
      b0i=taui*k_per**2*amassimp/amassgas/zimp**2
!
!     Pade aproximates...may use gamma functions later
      g0i=1.0
      g1i=1.0/(1+b0i)
      g2i=1.0/(1+b0i)*g1i
      g3i=1.0/(1+b0i)*g2i
      g12i=(g1i+g2i)/2.0
      g23i=(g2i+g3i)/2.0
!
      mimp=amassimp/amassgas
!
!... jek alpha-stab. fix for NCS (renorm only)
!
!     if(ilh.eq.1 .and. iky0.eq.2) write(*,*) 'ave_g = ',ave_g,ky
!
      w_s=ky
      w_d=(ghat*2.0/rmaj)*ky*ave_g
      w_d0=(ghat*2.0/rmaj)*ky*ave_g0
      w_cd=(gchat*2.0/rmaj)*ky*ave_g
      if(lprint.eq.3) write(*,*) 'ilh,w_d = ',ilh,w_d,w_d0
!
      k_par=park/(rmaj*q)*ave_k_par*sqrt((1.0+elong**2)/2.0)
!
!     sqrt((1.+elong**2)/2.) to get higher beta_crit prop to k_par**2
!     roughly same as betae-> betae/((1.+elong**2)/2.)
!     physically like shortening the connection length to good curv.
!
      if (iflagin(3).eq.2) then
       betae=betae_gf/(1.0+xparam(3))**2/(1.0+xparam(9))
      endif
!
      if (iflagin(3).eq.1) then
       k_par=park/(rmaj*q)*ave_k_par
       betae=betae_gf/(1.0+xparam(3))**2/(1.0+xparam(9))  &
            /((1.0+elong**2)/2.0)
      endif
!
!     we put the park enhancement directy into betae
!     ie park=.5 best at low beta and 2.5x.5=1.25 at high beta
!
!     option iglagin(3)=1 puts beta_crit elongation enhancement
!     directly into betae
!
!     option iflagin(3)=2 puts beta_crit elongation factor into
!     the connection length
!
!     an extra shape  factor 2 (triangularity) enhancement
!     is optained by (1.+xparam(9))=2.
!     if w_d negative flip sign of dissipative parts
! error 12/21       lamda_d=-conjg(lamda_d)
      if(w_d.lt.0.) then
       lamda_d=conjg(lamda_d)
       xmu=-conjg(xmu)
       xnu_par_par=-conjg(xnu_par_par)
       xnu_par_per=-conjg(xnu_par_per)
       xnu_per_par=-conjg(xnu_per_par)
       xnu_per_per=-conjg(xnu_par_per)
      endif
!
      reps=(1.0+xparam(15))*  &
         sqrt((rmin/rmaj)*(1.0+ave_cos)/(1.0+(rmin/rmaj)*ave_cos))
!      if(iky0.eq.1 .and. ilh.eq.1) write(*,*) reps
!
      if(nroot.le.3) reps=0.0
!
! fix trapped eletron MHD limit
! 3/4*reps*(1+rlte) + xt_mhd*(1-reps)+yt_mhd*rlte=1+rlte
! solve for xt_mhd and yt_mhd
! 3/4*reps+yt_mhd=1; 3/4*reps+xt_mhd*(1-reps)=1
!
      yt_mhd=(1-x43*(3.0/4.0)*reps)
      xt_mhd=(1.0-x43*(3.0/4.0)*reps)/(1.0-reps)
!
! collision detrapping retrapping model
!
      xnueff=(1.0+xparam(6))*xnu/(reps**2+1.E-6)
!
! very difficult get xnueff correct hince add enhancement factor
! and fit to Beer or GKS
!
      th=4.080
      tc=0.9180
      fh=0.1840
      fc=0.8160
!
      fnn=xnueff*((th/tc**(3.0/2.0))-(tc/th**(3.0/2.0)))/(th-tc)
      fnp=xnueff*(3.0/2.0)*  &
         ((1.0/th)**(3.0/2.0)-(1.0/tc)**(3.0/2.0))/(th-tc)
      fnf=xnueff*((fh/th**(3.0/2.0))+(fc/tc**(3.0/2.0)))
!
      fpn=xnueff*(2.0/3.0)*  &
         ((th/tc**(1.0/2.0))-(tc/th**(1.0/2.0)))/(th-tc)
      fpp=xnueff*((1.0/th)**(1.0/2.0)-(1.0/tc)**(1.0/2.0))/(th-tc)
      fpf=xnueff*(2.0/3.0)*((fh/th**(1.0/2.0))+(fc/tc**(1.0/2.0)))
!
!  collisional modes added with xnu_col
!  must fix for atomic mass dependence other than deuterium
      xnu_col=xparam(7)
      amass_e=2.7E-4*(2.0/amassgas)
!
! check adiabatic property that chkf should be 1.0  (finite xnu)
!
      chkf=(fnn*fpp-fpn*fnp)/((fnf*fpp-fpf*fnp)+epsilon)
      if (lprint.eq.2) write(6,*) 'chkf:', chkf
!
      if(neq.le.3) reps=0.0
!
      f0=-idelta+(1.0-dil)*apwt*(1/taui)*(1.0-g1)  &
                +zimp**2*aiwt*(1/taui)*(1.0-g1i)
      f1=1.0-reps + f0
!
      xadiabat=0.0
      if(nroot.le.6) then
        betae=0.0
        f0=f1
        xadiabat=1.0
      endif
      if(nroot.le.5) k_par=0.0
!
      betae0=betae+epsilon
      k_par0=k_par+epsilon
!
      if (lprint.eq.98) then
        write(2,*) 'ky,g1,g2,g3,g12,g23,w_s,w_d,w_cd: ',  &
          ky,g1,g2,g12,g23,w_s,w_d,w_cd
        write(2,*) 'f0,taui,k_par,reps: ',  &
          f0,taui,k_par,k_per,reps
        write(2,*) 'chi_par_1,chi_per_1,gam_par,gam_per:',  &
          chi_par_1,chi_per_1,gam_par,gam_per
        write(2,*) 'x_par,x_per,xmu:',  &
          x_par,x_per,xmu
        write(2,*) 'xnu_par_par,xnu_par_per,xnu_per_par,xnu_per_per:',  &
          xnu_par_par,xnu_par_per,xnu_per_par,xnu_per_per
        write(2,*) 'lamda_d,betae,xadiabat:',  &
          lamda_d,betae,xadiabat
        write(2,*) 'yt_mhd,xt_mhd:',  &
          yt_mhd,xt_mhd 
      endif
!
! matrix in order
! note ph=(n_i-n_t-n_u)/f0 results in (i,1)-(i,4)-(i,7) parts
!
! n_i equ #1
!
      amat(1,1)= (1.0-dil)*apwt*  &
        (-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
!
      amat(1,2)=  +xi*w_d*taui*0.50
!
      amat(1,3)=  +xi*w_d*taui*0.50
!
      amat(1,4)= -(-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
!
      amat(1,5)= 0.0
!
      amat(1,6)=  -xi*k_par
!
      amat(1,7)= -(-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
!
      amat(1,8)= 0.0
!
      amat(1,9)= aiwt*zimp*  &
        (-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
!
      amat(1,10)=0.0
!
      amat(1,11)=0.0
!
      amat(1,12)=0.0
!
! p_par equ #2
!
      amat(2,1)= (1.0-dil)*apwt*  &
        (-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0  &
        +k_par*chi_par_1  &
        -(xi*w_d*taui*3.0/2.0-w_d*taui*xnu_par_par)  &
        -(xi*w_d*taui*1.0/2.0-w_d*taui*xnu_par_per)
!
      amat(2,2)=  -k_par*chi_par_1  &
        +xi*w_d*taui*x_par +  &
        (xi*w_d*taui*3.0/2.0-w_d*taui*xnu_par_par)
!
      amat(2,3)=(xi*w_d*taui*1.0/2.0-w_d*taui*xnu_par_per)
!
      amat(2,4)= -(-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
!
      amat(2,5)=0.0
!
      amat(2,6)= -xi*gam_par*k_par
!
      amat(2,7)= -(-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
!
      amat(2,8)=0.0
!
      amat(2,9)= aiwt*zimp*  &
        (-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
!
      amat(2,10)=0.0
!
      amat(2,11)=0.0
!
      amat(2,12)=0.0
!
! p_per equ #3
!
      amat(3,1)= (1.0-dil)*apwt*  &
        (-xi*w_s*((rlni-rlti)*g2+2.0*rlti*g3)+xi*x_per*w_cd*g23)/f0  &
        +k_par*chi_per_1  &
        -(xi*w_d*taui-w_d*taui*xnu_per_per)  &
        -(xi*w_d*taui*1.0/2.0-w_d*taui*xnu_per_par)
!
      amat(3,2)= +(xi*w_d*taui*1/2-w_d*taui*xnu_per_par)
!
      amat(3,3)= -k_par*chi_per_1  &
        +xi*w_d*taui*x_per   +(xi*w_d*taui-w_d*taui*xnu_per_per)
!
      amat(3,4)=   &
       -(-xi*w_s*((rlni-rlti)*g2+2.0*rlti*g3)+xi*x_per*w_cd*g23)/f0
!
      amat(3,5)=0.0
!
      amat(3,6)= -xi*gam_per*k_par
!
      amat(3,7)=  &
        -(-xi*w_s*((rlni-rlti)*g2+2.0*rlti*g3)+xi*x_per*w_cd*g23)/f0
!
      amat(3,8)=0.0
!
      amat(3,9)= aiwt*zimp*  &
         (-xi*w_s*((rlni-rlti)*g2+2.0*rlti*g3)+xi*x_per*w_cd*g23)/f0
!
      amat(3,10)=0.0
!
      amat(3,11)=0.0
!
      amat(3,12)=0.0
!
! n_t equ #4
!
      amat(4,1)=(1.0-dil)*apwt*  &
        (-xi*w_s*rlne*reps*g0+xi*x43*3.0/4*w_cd*reps*g0)/f0  &
        -(1.0-dil)*apwt*(-reps*fnf*(1.0-reps)*g0/f0*xadiabat)
!
      amat(4,2)=0.0
!
      amat(4,3)=0.0
!
      amat(4,4)=-(-xi*w_s*rlne*reps*g0+xi*x43*3.0/4*w_cd*reps*g0)/f0 &
        -((1.0-reps)*fnn)  &
        -(-(-reps*fnf*(1.0-reps)*g0/f0*xadiabat))
!
      amat(4,5)= -xi*w_d*x43*3.0/4.0  &
       -((1.0-reps)*fnp)
!
      amat(4,6)=0.0
!
      amat(4,7)=-(-xi*w_s*rlne*reps*g0+xi*x43*3.0/4*w_cd*reps*g0)/f0  &
         -(-reps*fnf)
!
      amat(4,8)=0.0
!
      amat(4,9)=aiwt*zimp*  &
         (-xi*w_s*rlne*reps*g0+xi*x43*3.0/4*w_cd*reps*g0)/f0  &
         -aiwt*zimp*(-reps*fnf*(1.0-reps)*g0/f0*xadiabat)
!
      amat(4,10)=0.0
!
      amat(4,11)=0.0
!
      amat(4,12)=0.0
!
! p_t equ #5
!
      amat(5,1)= (1.0-dil)*apwt*  &
        (-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5.0/4*w_cd*reps*g0)/f0  &
        -(1.0-dil)*apwt*(-reps*fpf*(1.0-reps)*g0/f0*xadiabat)
!
      amat(5,2)=0.0
!
      amat(5,3)=0.0
!
      amat(5,4)=  &
        -(-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5.0/4*w_cd*reps*g0)/f0  &
                  +xi*w_d*lamda_d  &
        -((1.0-reps)*fpn)  &
        -(-(-reps*fpf*(1.0-reps)*g0/f0*xadiabat))
!
      amat(5,5)= -xi*w_d*x43*5.0/4.0-xi*w_d*lamda_d  &
        -((1.0-reps)*fpp)
!
      amat(5,6)=0.0
!
      amat(5,7)=  &
        -(-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5.0/4*w_cd*reps*g0)/f0  &
        -(-reps*fpf)
!
      amat(5,8)=0.0
!
      amat(5,9)= aiwt*zimp*  &
        (-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5.0/4*w_cd*reps*g0)/f0  &
        -aiwt*zimp*(-reps*fpf*(1.0-reps)*g0/f0*xadiabat)
!
      amat(5,10)=0.0
!
      amat(5,11)=0.0
!
      amat(5,12)=0.0
!
! u_par equ #6
!
      amat(6,1)=(1.0-dil)*apwt*  &
         (-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0  &
         -(betae/2.0)*(-xi*k_par*(2.0/betae0)*g0)/f0)
!
      amat(6,2)=-xi*k_par*taui
!
      amat(6,3)=0.0
!
      amat(6,4)=  &
        -(-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0)  &
        -(-(betae/2.0)*(-xi*k_par*(2.0/betae0)*g0)/f0)
!
      amat(6,5)=0.0
!
      amat(6,6)=  &
       +xi*w_d*(gam_par+gam_per)/2.0 -w_d*xmu
!
      amat(6,7)=  &
        -(-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0)  &
        -(-(betae/2.0)*(-xi*k_par*(2.0/betae0)*g0)/f0)  &
       -(betae/2.0)*xi*k_par*(2.0/betae0)/(1.0-reps)
!
      amat(6,8)=  &
        -(betae/2.0)*(-xi*w_s*(rlni*g1+rlti*g2))  &
        -(betae/2.0)*(-xi*w_s*rlne)  &
        +amass_e*xnu*xnu_col*k_per**2  &
        -(betae/2.0)*(-2.0/betae0*amass_e*xnu*xnu_col*k_per**2)
! note there is no double counting in last two terms
!
      amat(6,9)=aiwt*zimp*  &
        (-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0  &
        -(betae/2.0)*(-xi*k_par*(2.0/betae0)*g0)/f0)
!
      amat(6,10)=0.0
!
      amat(6,11)=0.0
!
      amat(6,12)=0.0
!
! n_u equ #7
!
      amat(7,1)=(1.0-dil)*apwt*  &
        (-xi*w_s*rlne*(1.0-reps)*g0+xi*w_cd*  &
        (1.0-x43*(3.0/4.0)*reps)*g0)/f0  &
        +(1.0-dil)*apwt*(-reps*fnf*(1.0-reps)*g0/f0*xadiabat)
!
      amat(7,2)=0.0
!
      amat(7,3)=0.0
!
      amat(7,4)=  &
        -(-xi*w_s*rlne*(1.0-reps)*g0+xi*w_cd*  &
        (1.0-x43*(3.0/4.0)*reps)*g0)/f0  &
        +((1.0-reps)*fnn)  &
        +(-(-reps*fnf*(1.0-reps)*g0/f0*xadiabat))
!
      amat(7,5)= +((1.0-reps)*fnp)
!
      amat(7,6)=-xi*k_par
!
      amat(7,7)=  &
        -(-xi*w_s*rlne*(1.0-reps)*g0+xi*w_cd*  &
        (1.0-x43*(3.0/4.0)*reps)*g0)/f0  &
        -xi*w_d*xt_mhd  &
        +(-reps*fnf)
!
      amat(7,8)= -xi*k_par*(-k_per**2)  &
        -xi*w_d*yt_mhd*(w_s*(betae/2.0)/k_par0*rlte)
!
      amat(7,9)=aiwt*zimp*  &
        (-xi*w_s*rlne*(1.0-reps)*g0+xi*w_cd*  &
        (1.0-x43*(3.0/4.0)*reps)*g0)/f0  &
        +aiwt*zimp*(-reps*fnf*(1.0-reps)*g0/f0*xadiabat)
!
      amat(7,10)=0.0
!
      amat(7,11)=0.0
!
      amat(7,12)=0.0
!
! a_par equ #8
!
      amat(8,1)=(1.0-dil)*apwt*(-xi*k_par*(2.0/betae0)*g0/f0)
!
      amat(8,2)=0.0
!
      amat(8,3)=0.0
!
      amat(8,4)=-(-xi*k_par*(2.0/betae0)*g0/f0)
!
      amat(8,5)=0.0
!
      amat(8,6)=0.0
!
      amat(8,7)=-(-xi*k_par*(2.0/betae0)*g0/f0)  &
        +xi*k_par*(2.0/betae0)/(1.0-reps)
!
      amat(8,8)=-xi*w_s*rlne  &
        -(2.0/betae0)*amass_e*xnu*xnu_col*(k_per**2)
!
      amat(8,9)=aiwt*zimp*(-xi*k_par*(2.0/betae0)*g0/f0)
!
      amat(8,10)=0.0
!
      amat(8,11)=0.0
!
      amat(8,12)=0.0
!
! n_im equ #9
!
      amat(9,1)= (1.0-dil)*apwt*  &
        (-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
!
      amat(9,2)=0.0
!
      amat(9,3)=0.0
!
      amat(9,4)= -(-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
!
      amat(9,5)= 0.0
!
      amat(9,6)= 0.0
!
      amat(9,7)= -(-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
!
      amat(9,8)= 0.0
!
      amat(9,9) = aiwt*zimp*  &
        (-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
!
      amat(9,10)= +xi*w_d*taui*0.50/zimp
!
      amat(9,11)= +xi*w_d*taui*0.50/zimp
!
      amat(9,12)= -xi*k_par
!
! pim_par equ #10
!
      amat(10,1)= (1.0-dil)*apwt*  &
        (-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0
!
      amat(10,2)=0.0
!
      amat(10,3)=0.0
!
      amat(10,4)= -(-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0
!
      amat(10,5)=0.0
!
      amat(10,6)=0.0
!
      amat(10,7)= -(-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0
!
      amat(10,8)=0.0
!
      amat(10,9)= aiwt*zimp* &
         (-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0  &
         +k_par*chi_par_1/sqrt(mimp)  &
         -(xi*w_d*taui*3.0/2.0-w_d*taui*xnu_par_par)/zimp  &
         -(xi*w_d*taui*1.0/2.0-w_d*taui*xnu_par_per)/zimp
!
      amat(10,10)= -k_par*chi_par_1/sqrt(mimp)  &
        +xi*w_d*taui*x_par/zimp +(xi*w_d*taui*3.0/2.0  &
        -w_d*taui*xnu_par_par)/zimp
!
      amat(10,11)= (xi*w_d*taui*1.0/2.0-w_d*taui*xnu_par_per)/zimp
!
      amat(10,12)= -xi*gam_par*k_par
!
! pim_per equ #11
!
      amat(11,1)= (1.0-dil)*apwt*  &
        (-xi*w_s*((rlnimp-rlti)*g2i+2.0*rlti*g3i)+  &
        xi*x_per*w_cd*g23i)/f0
!
      amat(11,2)= 0.0
!
      amat(11,3)= 0.0
!
      amat(11,4)= -(-xi*w_s*((rlnimp-rlti)*g2i+2.0*rlti*g3i)+  &
        xi*x_per*w_cd*g23i)/f0
!
      amat(11,5)=0.0
!
      amat(11,6)=0.0
!
      amat(11,7)= -(-xi*w_s*((rlnimp-rlti)*g2i+2.0*rlti*g3i)+  &
        xi*x_per*w_cd*g23i)/f0
!
      amat(11,8)=0.0
!
      amat(11,9)= aiwt*zimp*  &
        (-xi*w_s*((rlnimp-rlti)*g2i+2.0*rlti*g3i)+  &
        xi*x_per*w_cd*g23i)/f0  &
        +k_par*chi_per_1/sqrt(mimp)  &
        -(xi*w_d*taui-w_d*taui*xnu_per_per)/zimp  &
        -(xi*w_d*taui*1.0/2.0-w_d*taui*xnu_per_par)/zimp
!
      amat(11,10)= +(xi*w_d*taui*1/2-w_d*taui*xnu_per_par)/zimp
!
      amat(11,11)= -k_par*chi_per_1/sqrt(mimp)  &
        +xi*w_d*taui*x_per/zimp &
        +(xi*w_d*taui-w_d*taui*xnu_per_per)/zimp
!
      amat(11,12)= -xi*gam_per*k_par
!
! uim_par equ #12
!gms 5/21/99 added gamma_p to amat(12,1),amat(12,4),amat(12,7),amat(12,9)
!    added xnu_col term to amat(12,8)
!    fixed mimp factor in amat(12,12)
!
      amat(12,1)=(1.0/mimp)*(1.0-dil)*apwt*  &
        ((-xi*k_par*g1i/f0)*zimp  &
        -xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0  &
        -(betae/2.0)*zimp*(-xi*k_par*(2.0/betae0)*g0i)/f0)
!
      amat(12,2)=0.0
!
      amat(12,3)=0.0
!
      amat(12,4)=  &
       -(1.0/mimp)*(-xi*k_par*g1i/f0)*zimp  &
       -(1.0/mimp)*(-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0)  &
       -(1.0/mimp)*(-(betae/2.0)*(-xi*k_par  &
       *(2.0/betae0)*g0i)/f0)*zimp
!
      amat(12,5)=0.0
!
      amat(12,6)=0.0
!
       amat(12,7)=  &
       -(1.0/mimp)*(-xi*k_par*g1i/f0)*zimp  &
       -(1.0/mimp)*(-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0)  &
       -(1.0/mimp)*(-(betae/2.0)*  &
        (-xi*k_par*(2.0/betae0)*g0i)/f0)*zimp  &
       -(1.0/mimp)*(betae/2.0)*xi*  &
        k_par*(2.0/betae0)/(1.0-reps)*zimp
!
      amat(12,8)=  &
        -(1.0/mimp)*(betae/2.0)*(-xi*w_s*(rlnimp*g1i+rlti*g2i))  &
        -(1.0/mimp)*(betae/2.0)*(-xi*w_s*rlne)*zimp  &
        +(1.0/mimp)*zimp*amass_e*xnu*xnu_col*k_per**2  &
        -(1.0/mimp)*(betae/2.0)*  &
         (-2.0/betae0*amass_e*xnu*xnu_col*k_per**2)*zimp
!
      amat(12,9)=(1.0/mimp)*aiwt*zimp*  &
         ((-xi*k_par*g1i/f0)*zimp  &
        -xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0  &
        -(betae/2.0)*zimp*(-xi*k_par*(2.0/betae0)*g0i)/f0)
!
      amat(12,10)=-(1.0/mimp)*xi*k_par*taui
!
      amat(12,11)=0.0
!
      amat(12,12)= (1.0/mimp)*(xi*w_d*(gam_par+gam_per)  &
        /2.0/zimp -w_d*xmu/zimp)
!
! put in rot shear stabilization and possible source of gyrobohm breaking
! and model damping kdamp
!
!***********************************************************************
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
! solve 12x12 complex
! -xi*omega*v(i)=sum_j amat(i,j)*v(j)  omega=freq+xi*gamma
! upto nroot
! order with max gamma and find eigenvector v(i) with ant fixed norm.
!
!...Fill matricies for eigenvalue equation
!
      do j1=1,neq
        rr(j1) = 0.00
        ri(j1) = 0.00
        do j2=1,neq
          ar(j1,j2) = REAL(  amat(j1,j2) )
          ai(j1,j2) = AIMAG( amat(j1,j2) )
!...test tmp
!         ai(j1,j2) = 0.0
!         ar(j1,j2) = 0.0
!         if (j1.eq.j2) ar(j1,j2)=j1
!
          vr(j1,j2) = 0.00
          vi(j1,j2) = 0.00
        enddo
      enddo
!
!...diagnostic output
!
      if ( lprint .gt. 6 ) then
        write (1,*)
        write (1,*) ' ar(j1,j2)  j2 ->'
        do j1=1,neq
          write (1,192) (ar(j1,j2),j2=1,neq)
        enddo
!
        write (1,*)
        write (1,*) ' ai(j1,j2)  j2->'
        do j1=1,neq
          write (1,192) (ai(j1,j2),j2=1,neq)
        enddo
 192    format (1p8e10.2)
! 193    format (1p8e12.4)
      endif
!
!..find the eigenvalues and eigenvectors
!
!.. eigen_gf = 0 use cgg solver (default)
!..          = 1 use tomsqz solver
!..          = 2 use zgeev solver
!.. not longer used:
!        call f02ake( ar,iar,ai,iai,ieq,rr,ri,vr,ivr,vi,ivi,
!     >               intger, ifail )
!
        ifail = 0
!
!gms        if (eigen_gf .eq. 2 ) then
!gms
!gms        jobvl = 'N'
!gms        jobvr = 'V'
!gms        do j1=1,neq
!gms         do j2=1,ieq
!gms           mata(j1,j2) = dcmplx(ar(j1,j2),ai(j1,j2))
!gms         enddo
!gms        enddo
!gms
!gms        call zgeev(jobvl,jobvr,ieq,mata,neq,w,cvl,neq,cvr,
!gms     &             neq,work,lwork,rwork,ifail)
!gms        do j1=1,neq
!gms         rr(j1) = real(w(j1))
!gms         ri(j1) = imag(w(j1))
!gms         do j2=1,ieq
!gms           vr(j1,j2) = real(cvr(j1,j2))
!gms           vi(j1,j2) = imag(cvr(j1,j2))
!gms         enddo
!gms        enddo
!gms
!gms        elseif (eigen_gf .eq. 1 ) then
!gms
!gms        do j2=1,neq
!gms           do j1=1,neq
!gms              bi(j1,j2)=0.00
!gms              if(j1.eq.j2) then
!gms                 br(j1,j2)=1.00
!gms              else
!gms                 br(j1,j2)=0.00
!gms              endif
!gms           enddo
!gms        enddo
!gms
!gms        call r8tomsqz(neq,ieq,ar,ai,br,bi, rr,ri,beta_tom, vr,vi, ifail)
!gms
!gms        do j1=1,ieq
!gms           ztemp1 = beta_tom(j1)
!gms           if ( abs(beta_tom(j1)) .lt. epsilon ) ztemp1 = epsilon
!gms           rr(j1)=rr(j1) / ztemp1
!gms           ri(j1)=ri(j1) / ztemp1
!gms        enddo
!gms
!gms        else
!gms
        matz=1
!       write(*,*) 'neq = ',neq
!       write(*,*) 'ieq = ',ieq
!       write(*,*) 'matz = ',matz
!       write (*,*) ' ar(j1,j2)  j2 ->'
!       do j1=1,neq
!         write (*,193) (ar(j1,j2),j2=1,neq)
!       enddo
!       write (*,*) ' ai(j1,j2)  j2 ->'
!       do j1=1,neq
!         write (*,193) (ai(j1,j2),j2=1,neq)
!       enddo
!
        call cgg_glf23(neq,ieq,ar,ai,rr,ri,matz,vr,vi,fv1,fv2,fv3,ifail)
!
!       write (*,*) ' wr(j1) and wi(j1)'
!       do j1=1,neq
!         write (*,193) rr(j1), ri(j1)
!       enddo
!       write (*,*) ' zr(j1,j2)  j2 ->'
!       do j1=1,neq
!         write (*,193) (vr(j1,j2),j2=1,neq)
!       enddo
!       write (*,*) ' zi(j1,j2)  j2 ->'
!       do j1=1,neq
!         write (*,193) (vi(j1,j2),j2=1,neq)
!       enddo
!
!gms        endif
!
        if ( lprint .gt. 1 ) then
          write (1,*) ifail,' = ifail routine '
        endif
!
!..print eigenvalues
!
        if ( lprint .gt. 6 ) then
          write (1,121)
          do j=1,ieq
            write (1,122) rr(j), ri(j)
          enddo
 121      format (/' Solution of the eigenvalue equations'  &
           /t4,'real   ',t18,'imag   ')
 122      format (1p2e14.5)
        endif
!
!...Store the complex eigenvectors and eigenvalues
!...Note the routines here solve A.v = lambda v
!...but that the equation solved is A.v = -i omega v
!...The i-th column of v is the i-th eigenvector
!
        do j1=1,ieq
          zomega(j1) = xi*(rr(j1)+xi*ri(j1))
          do j2=1,ieq
            zevec(j2,j1) = vr(j2,j1) + xi*vi(j2,j1)
          enddo
        enddo
!
        if ( lprint .gt. 6 ) then
          write (6,123)
          do j=1,ieq
            write (6,122) real(zomega(j)), aimag(zomega(j))
          enddo
 123      format (/' Multiplied by i: '  &
           /t4,'zomegar',t18,'zomegai')
        endif
        do iroot=1,4
!
!
!..save growth rates and frequencies in real variables
!
        zgamax = 0.00
!temp
        zgamax = -1.E10
        jmax=0
        gamma=0.0
        do j=1,ieq
         if(j.ne.jroot(1).and.j.ne.jroot(2).and.j.ne.jroot(3)) then
          if (aimag(zomega(j)).gt. zgamax) then
            zgamax = aimag(zomega(j))
            jmax=j
          endif
         endif
        enddo
!
        freq = REAL( zomega(jmax) )
        gamma = aimag( zomega(jmax) )
        if(lprint.eq.3) write(*,'(1x,a17,i3,1p3e12.4)') &
          'ilh,gamma,freq = ',ilh,gamma,freq
!
! skip stable modes
!        if(gamma.lt.0.0)go to 775
!
         jroot(iroot)=jmax
!
!temp        if(zgamax.lt.zepsqrt) gamma=0.
!
        if (jmax.ne.0) then
         gammaroot(iroot)=gamma
         freqroot(iroot)=freq
!
        do j=1,12
         v(j)=0.0
        enddo
        do j=1,ieq
          v(j) = zevec(j,jmax)
        enddo
!
!***********************************************************************
!
      n_i=0.0
      p_par=0.0
      p_per=0.0
      n_t=0.0
      p_t=0.0
      u_par=0.0
      n_u=0.0
      a_par=0.0
      n_im=0.0
      p_im_par=0.0
      p_im_per=0.0
!     u_im_par=0.0
!
      t_u=0.0
!
      n_i=v(1)
      p_par=v(2)
      p_per=v(3)
      if(ieq.ge.5) n_t=v(4)
      if(ieq.ge.5) p_t=v(5)
      if(ieq.ge.6) u_par=v(6)
      if(ieq.ge.8) n_u=v(7)
      if(ieq.ge.8) a_par=v(8)
      if(ieq.ge.9) n_im=v(9)
      if(ieq.eq.12) p_im_par=v(10)
      if(ieq.eq.12) p_im_per=v(11)
!     if(ieq.eq.12) u_im_par=v(12)
 
      if (ieq.ge.8) ph=((1.0-dil)*apwt*n_i+aiwt*zimp*n_im-n_t-n_u)/f0
      if (ieq.lt.8) ph= ((1.0-dil)*apwt*n_i-n_t)/f0
      if (ieq.le.3) ph= ((1.0-dil)*apwt*n_i)/f0
      t_u=(betae/2.0)*w_s/k_par0*rlte*a_par
!
      n_e=(1.0-dil)*apwt*(n_i-(g0-g1)/taui*ph)  &
           +aiwt*zimp*(n_im-zimp*(g0i-g1i)/taui*ph)
!
! impurity trace convective limit
!
      if (aiwt.lt.-epsilon) then
       n_im=0.0
       do j=1,8
        n_im=n_im+amat(9,j)*v(j)/(-xi*freq+gamma)
       enddo
      endif
!
! idelta=xi*yparam(1)+yparam(2)   for trapped electrons
!
      yparam(1)=aimag(-(n_t-reps*ph)/ph)
      yparam(2)=REAL(-(n_t-reps*ph)/ph)
!
      chknu=n_u/(1.0-reps)/ph
      chknt=n_t/reps/ph
      chknt2=n_t*(f0+reps)/(reps*((1.0-dil)*apwt*n_i+aiwt*zimp*n_im))
      if (lprint.eq.2) write (6,*) 'chknu,chknt,chknt2:',  &
          chknu,chknt,chknt2
!
! non linear saturation rule
!
      gamma_r= 0.20*3.0/2.0*abs(w_d)*taui   !only scaling important
      if(version_gf.eq.2) gamma_r= 0.20*3.0/2.0*abs(w_d0)*taui
      if(lprint.eq.3) write(*,'(1x,a17,i3,1p3e12.4)')  &
         'ilh,gamma_r    = ',ilh,gamma_r
!
      gamma_net=gamma-abs(alpha_star*gamma_star  &
               +alpha_e*gamma_e+alpha_mode*gamma_mode)-kdamp
!
      gamma_net=max(gamma_net,xparam(8)*gamma)
      if( gamma_net.gt.0.0)then
       ph_m=gamma_net**(1.0-adamp)*gamma_r**adamp/(k_m*ky)
       if(xparam(26).gt.0.0) then
        ph_m=ph_m*((xparam(26)*gamma_net/gamma_r)**xparam(27))  &
        /(1.0+(xparam(26)*gamma_net/gamma_r)**xparam(27))
       endif
! set flag ngrow_k_gf: found at least one unstable mode for this k
!gms       if(ipert_gf.eq.0)ngrow_k_gf(iky0)=1
      endif
!
      if( gamma_net.le.0.) ph_m=0.0
!
      if(xparam(24).gt.0.0) then
        if(gamma.gt.0.0)then
          ph_m=abs(gamma)**(1.0-adamp)*gamma_r**adamp/(k_m*ky)/  &
               sqrt(1.0+(abs(alpha_star*gamma_star+  &
               alpha_e*gamma_e+alpha_mode*gamma_mode)/  &
               (abs(gamma)+.000010))**xparam(24))
!gms          if(ipert_gf.eq.0)ngrow_k_gf(iky0)=1
        else
           ph_m=0.0
        endif
      endif
!
! 7.17.96
      if(xparam(22).gt.0) then
        if(gamma.gt.0.) ph_m=gamma**(1.0-adamp-xparam(22))  &
                             *gamma_r**adamp/(k_m*ky)
        if(gamma.le.0.) ph_m=0.0
      endif
!
         phi_norm=0
         phi_normroot(iroot)=0.0
!
       if( ph_m.gt.0.) then
!
       phi_norm=ph_m*ph_m/ABS((conjg(ph)*ph))
!
       phi_normroot(iroot)=phi_norm
!
! note only real part survives in diffusivities
!    ...units are c_s*rho_s**2/a
! magnetic futter component is too small to worry about
!
      d_hat  = phi_norm*REAL(conjg(n_i)*(-xi*ky*ph))  &
               +d_hat
!
      d_im_hat = phi_norm*REAL(conjg(n_im)*(-xi*ky*ph))  &
                 +d_im_hat
!
      chii_hat = phi_norm*3.0/2.0*  &
                 REAL(conjg((1.0/3.0)*p_par+  &
                 (2.0/3.0)*p_per)*(-xi*ky*ph))  &
                 +chii_hat
      chi_im_hat= aiwt/apwt*xparam(21)*phi_norm*3.0/2.0*  &
               REAL(conjg((1.0/3.0)*p_im_par+(2.0/3.0)*p_im_per)*  &
               (-xi*ky*ph)) + chi_im_hat
!
      chie_hat = phi_norm*3.0/2.0*  &
                 REAL(conjg(p_t+n_u+t_u)*(-xi*ky*ph))  &
                 +chie_hat
!
      if(lprint.eq.3) then
        write(*,'(1x,a25,i3,1p3e12.4)') 'ilh,abs((conjg(ph)*ph) = ', &
             ilh,ABS((conjg(ph)*ph))
        write(*,'(1x,a20,i3,1p3e12.4)') 'ilh,phi_norm,ph_m = ',  &
             ilh,phi_norm,ph_m
        write(*,'(1x,a20,i3,1p3e12.4)') 'ilh,gamma_net     = ',  &
              ilh,gamma_net
        write(*,'(1x,a20,i3,1p3e12.4)') 'ilh,chii,e_hat    = ',  &
             ilh,chii_hat,chie_hat
        write(*,'(1x,a12,i3,1p3e12.4)') 'ilh,ky,km = ',ilh,ky,k_m
        write(*,'(1x,a12,i3,1p3e12.4)') 'ilh,w_d0 = ',ilh,w_d0
      endif
!
! electron to ion energy exchange in units n0*t0*c_s/a*(rho_a/a)**2
! note here we interpret QLT d/dt=-xi*freq dropping gamma part to
! avoid getting nonzero result for n_e->ph adiabatic limit
! ie <(d n_e/dt)*conjg(ph)>_time ave -> 0 adiabatic limit
! note  (-1) means exch_hat is electron to ion rate or
! ion heating rate, ie positive exch_hat cools electrons and heats ions
!
      exch_hat = phi_norm*(-1.0)*  &
                 REAL(conjg(-xi*freq*n_e)*ph)  &
                 +exch_hat
!
      eta_par_hat=(1.0-xparam(14))*phi_norm*  &
                  REAL(conjg(u_par)  &
                  *(-xi*ky*ph))/(gamma_p+epsilon)*(-gamma_p*alpha_n)  &
                  +xparam(14)*phi_norm*REAL(conjg(  &
                  -xi*ky*gamma_p*(-gamma_p*alpha_n)*g1*ph/  &
                  (-xi*freq+gamma))  &
                  *(-xi*ky*ph))*(-gamma_p*alpha_n)  &
                  +eta_par_hat
!
      eta_per_hat = phi_norm*  &
                    REAL(conjg(-ky*(ky*shat*rms_theta)*ph)*  &
                    (ph+taui*((1.0/3.0)*p_par+(2.0/3.0)*p_per)))*  &
                    (-gamma_p*alpha_n)  &
                    +eta_per_hat
!
       endif
       endif
      enddo
! 777  continue
!
      if(ilh.eq.1) then
!
       do i=1,nmode
        yparam_k_gf(i,iky)=yparam(i)
       enddo
       xkyf_k_gf(iky)=kyf
!
       do iroot=1,4
        gamma_k_gf(iroot,iky)=gammaroot(iroot)
        freq_k_gf(iroot,iky)=freqroot(iroot)
        phi_norm_k_gf(iroot,iky)=phi_normroot(iroot)
       enddo
!
       diff_k_gf(iky)=d_hat
       diff_im_k_gf(iky)=d_im_hat
       chii_k_gf(iky)=chii_hat
       chi_im_k_gf(iky)=chi_im_hat
       chie_k_gf(iky)=chie_hat
       exch_k_gf(iky)=exch_hat
!
!test       exch_gf=-freq_gf(1)/xkyf_gf*diff_gf*rlni
!
       eta_par_k_gf(iky)=eta_par_hat
       eta_per_k_gf(iky)=eta_per_hat
!
! b_pol/b_phi=rmin/(rmaj*q)
!
       eta_phi_k_gf(iky)=eta_par_k_gf(iky)+  &
                         rmin/(rmaj*q)*eta_per_k_gf(iky)
!
      endif
!
! computed high-k ETG electron transport at each ky
! Note: not added to ITG chi-e here ... done after ky loop
      chie_e_k_gf(iky)=0.0
      if(ilh.eq.2) then
        chie_e_k_gf(iky)=xparam(10)*chii_hat*  &
                         taui_gf**(3.0/2.0)/  &
                         (1836.0*amassgas_gf)**.50
      endif
      if(lprint.eq.3) write(*,'(1x,a27,i3,1p3e12.4)')  &
          'ilh,chie_e_k_gf,chii_k_gf =',  &
          ilh,chie_e_k_gf(iky),chii_k_gf(iky)
! end ky loop
      enddo
!
! end big loop on ilh ... no longer used
!     enddo
!
! check to see if any unstable modes were found
!
!gms      if(ipert_gf.eq.0)then
!gms        do j=1,nmode
!gms         if(ngrow_k_gf(j).ne.0)ngrow_k_gf(0)=1
!gms        enddo
!gms        if(ngrow_k_gf(0).eq.0)go to 888
!gms      endif
!
!
!***********************************************************************
! initializations for summations over ky
!
      anorm_k=0.0
      diff_gf=0.0
      diff_im_gf=0.0
      chii_gf=0.0
      chi_im_gf=0.0
      chie_gf=0.0
      exch_gf=0.0
      eta_par_gf=0.0
      eta_per_gf=0.0
      eta_phi_gf=0.0
      chie_e_gf=0.0
!
!      anorm_glob=0.0
!      diff_glob=0.0
!      diff_im_glob=0.0
!      chii_glob=0.0
!      chie_glob=0.0
!      exch_glob=0.0
!      eta_par_glob=0.0
!      eta_per_glob=0.0
!      eta_phi_glob=0.0
!      chie_e_glob=0.0
!
! Sum ITG and ETG transport
! over logarithmic ky grid (d ky=ky*d yk)
!
      do iky=1,ikymax_gf
       del_k=xkyf_k_gf(iky)
       anorm_k=anorm_k+del_k
       diff_gf=diff_gf+diff_k_gf(iky)*del_k
       diff_im_gf=diff_im_gf+diff_im_k_gf(iky)*del_k
       chii_gf=chii_gf+chii_k_gf(iky)*del_k
       chi_im_gf=chi_im_gf+chi_im_k_gf(iky)*del_k
       chie_gf=chie_gf+chie_k_gf(iky)*del_k
       exch_gf=exch_gf+exch_k_gf(iky)*del_k
       eta_par_gf=eta_par_gf+eta_par_k_gf(iky)*del_k
       eta_per_gf=eta_per_gf+eta_per_k_gf(iky)*del_k
       eta_phi_gf=eta_phi_gf+eta_phi_k_gf(iky)*del_k
       chie_e_gf=chie_e_gf+chie_e_k_gf(iky)*del_k
       if(lprint.eq.3) then
         write(*,'(2x,i3,1p8e11.3)') iky,   &
             xkyf_k_gf(iky), chii_gf, chie_gf,  &
             chie_e_gf
         write(*,'(2x,i3,1p4e11.3," gamma")') iky,  &
            gamma_k_gf(1,iky),gamma_k_gf(2,iky),  &
            gamma_k_gf(3,iky),gamma_k_gf(4,iky)
       endif
      enddo
      if(lprint.eq.3) write(*,*) 'anorm = ',anorm_k
!
! Add ITG and ETG electron transport
!
      chie_gf=chie_gf + 1.0*chie_e_gf
!
      diff_gf=diff_gf/anorm_k
      diff_im_gf=diff_im_gf/anorm_k
      chii_gf=chii_gf/anorm_k
      chi_im_gf=chi_im_gf/anorm_k
      chie_gf=chie_gf/anorm_k
      exch_gf=exch_gf/anorm_k
      eta_par_gf=eta_par_gf/anorm_k
      eta_per_gf=eta_per_gf/anorm_k
      eta_phi_gf=eta_phi_gf/anorm_k
      chie_e_gf=chie_e_gf/anorm_k
!
! pick off maximum gamma
!
      do iroot=1,4
       gamma_k_max=-1.E6
       do iky=1,ikymax_gf
        if(gamma_k_gf(iroot,iky).gt.gamma_k_max) then
         gamma_k_max=gamma_k_gf(iroot,iky)
         gamma_gf(iroot)=gamma_k_gf(iroot,iky)
         freq_gf(iroot)=freq_k_gf(iroot,iky)
         xky_gf(iroot)=xkyf_k_gf(iky)
        endif
       enddo
      enddo
!
! pick off 2nd maximum gamma
!
       gamma_k_max=-1.E6
       do iky=1,ikymax_gf
        if( (gamma_k_gf(1,iky).gt.gamma_k_max) .and.  &
            (gamma_k_gf(1,iky).lt.gamma_gf(1)) ) then
         gamma_k_max=gamma_k_gf(1,iky)
         gamma_gf(2)=gamma_k_gf(1,iky)
         freq_gf(2)=freq_k_gf(1,iky)
         xky_gf(2)=xkyf_k_gf(iky)
        endif
       enddo
!
       gamma_r_gf=gamma_r
!
!       write(*,*) gamma_gf(1), gamma_gf(2), xky_gf(1), xky_gf(2)
!
! print to file log
!      write(*,66)chii_gf,(gamma_gf(j),j=1,4)
! 66    format(f14.9,4f14.9)
! 67    format(2i2,f14.9)
!
      if(xparam(22).gt.0.) then
       phi_renorm=1.0
       gamma_gross_net=gamma_gf(1)-abs(alpha_star*gamma_star  &
               +alpha_e*gamma_e+alpha_mode*gamma_mode)-kdamp
       if(gamma_gross_net.gt.0.)  &
        phi_renorm=gamma_gross_net**(xparam(22)*2.0)
       if(gamma_gross_net.le.0.) phi_renorm=0.0
!
      diff_gf=diff_gf*phi_renorm
      diff_im_gf=diff_im_gf*phi_renorm
      chii_gf=chii_gf*phi_renorm
      chi_im_gf=chi_im_gf*phi_renorm
      chie_gf=chie_gf*phi_renorm
      exch_gf=exch_gf*phi_renorm
      eta_par_gf=eta_par_gf*phi_renorm
      eta_per_gf=eta_per_gf*phi_renorm
      eta_phi_gf=eta_phi_gf*phi_renorm
      chie_e_gf=chie_e_gf*1.00
!
      endif
!
! put in cnorm_gf 12/22/96
!
      diff_gf=cnorm_p_gf*diff_gf
      diff_im_gf=cnorm_gf*diff_im_gf
      chii_gf=cnorm_gf*chii_gf
      chi_im_gf=cnorm_gf*chi_im_gf
      chie_gf=cnorm_gf*chie_gf
      exch_gf=cnorm_gf*exch_gf
      eta_par_gf=cnorm_gf*eta_par_gf
      eta_per_gf=cnorm_gf*eta_per_gf
      eta_phi_gf=cnorm_gf*eta_phi_gf
      chie_e_gf=cnorm_gf*chie_e_gf
!
! compute contributions from individual k's
!
      do iky=1,ikymax_gf
       del_k=xkyf_k_gf(iky)
       diff_k_gf(iky)=diff_k_gf(iky)*del_k/anorm_k*cnorm_p_gf
       chie_k_gf(iky)=chie_k_gf(iky)*del_k/anorm_k*cnorm_gf
       chie_e_k_gf(iky)=chie_e_k_gf(iky)*del_k/anorm_k*cnorm_gf
       chii_k_gf(iky)=chii_k_gf(iky)*del_k/anorm_k*cnorm_gf
       chi_im_k_gf(iky)=chi_im_k_gf(iky)*del_k/anorm_k*cnorm_gf
       eta_phi_k_gf(iky)=eta_phi_k_gf(iky)*del_k/anorm_k*cnorm_gf
!      write(*,'(2x,i3,1p8e11.3)') iky, xkyf_k_gf(iky),
!    >      chii_k_gf(iky), del_k, anorm_k
      enddo
!
      if(lprint.gt.0) then
          write(1,*) 'gamma_gf=',  gamma_gf
          write(1,*) 'freq_gf=',  freq_gf
          write(1,*) 'ph_m=',  ph_m
          write(1,*) 'diff_gf=',   diff_gf
          write(1,*) 'diff_im_gf=',   diff_im_gf
          write(1,*) 'chii_gf=', chii_gf
          write(1,*) 'chi_im_gf=', chi_im_gf
          write(1,*) 'chie_gf=', chie_gf
          write(1,*) 'exch_gf=', exch_gf
      endif
!
      if (lprint.eq.98) then
        write(6,*) 'rlti,rlte,rlne,rlni,rlnimp: ',  &
          rlti,rlte,rlne,rlni,rlnimp
       write(6,*) 'chii,chie,diff,diff_im: ',  &
          chii_gf,chie_gf,diff_gf,diff_im_gf
       write(6,*) 'gamma_gf,freq_gf,ph_m: ',  &
          gamma_gf,freq_gf,ph_m
       write(6,*) 'jmax: ',  &
          jmax
        write(2,*) 'rlti,rlte,rlne,rlni,rlnimp: ',  &
          rlti,rlte,rlne,rlni,rlnimp
       write(2,*) 'chii,chie,diff,diff_im: ',  &
          chii_gf,chie_gf,diff_gf,diff_im_gf
       write(2,*) 'gamma_gf,freq_gf,ph_m: ',  &
          gamma_gf,freq_gf,ph_m
       write(2,*) 'jmax: ',  &
          jmax
 
        write (2,*) ' ar(j1,j2)  j2 ->'
        do j1=1,neq
          write (2,*) (ar(j1,j2),j2=1,neq)
        enddo
!
        write (2,*)
        write (2,*) ' ai(j1,j2)  j2->'
        do j1=1,neq
          write (2,*) (ai(j1,j2),j2=1,neq)
        enddo
!
        write (2,*) ' vr(j1,j2)  j2 ->'
        do j1=1,neq
          write (2,*) (vr(j1,j2),j2=1,neq)
        enddo
!
        write (2,*)
        write (2,*) ' vi(j1,j2)  j2->'
        do j1=1,neq
          write (2,*) (vi(j1,j2),j2=1,neq)
        enddo
!
      endif
!
! 999  continue
      if (lprint.gt.0) close(1)
!
      return
!
! return for case with  no unstable modes
!
! 888  continue
      diff_gf=0.0
      diff_im_gf=0.0
      chii_gf=0.0
      chi_im_gf=0.0
      chie_gf=0.0
      exch_gf=0.0
      eta_par_gf=0.0
      eta_per_gf=0.0
      eta_phi_gf=0.0
      chie_e_gf=0.0
      do j1=1,4
        gamma_gf(j1)=0.0
        freq_gf(j1)=0.0
        xky_gf(j1)=0.0
      enddo
      return
!
      end subroutine glf23_TM
!
