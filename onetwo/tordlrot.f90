!
  MODULE tordlrot
  USE param, only : ksplin,kbctim,kj,ke,kb
!
! --- input data parameters (explained in subroutine INIT)
!
      character*128     rgc_string
!      real *8 lqn_mgbr
      integer,parameter  :: n_mgbr = 30
      real*8                                                         &
                   angrotin(ksplin,kbctim),rangrot(ksplin,kbctim),   &
                   angrm2d(4),angrcple,                              &
                   xmtmdifs(ksplin),rmtmdifs(ksplin),                &
                   rgc(kj),times_rgc,timee_rgc,rgce_mult,            &
                   time_mgbr(n_mgbr),berrqn_mgbr(n_mgbr),            &
                   cb_mgbr,qn_mgbr,lqn_mgbr,toff_mgbr,berrqn_val,    &
                   rgci_mult,rgca,rgcb,rgcc,rgcd
!
      integer                                                        &
                   njinang,irgc,iangrot,itangrot,iwangrot,nbeamtcx,  &
                   momtm_file, nt_mgbr, knotsang(kbctim)
 
!
!
! --- neutral beam angular momentum source arrays
!
      real *8                                                        &
                   spbrsav(kj,ke,kb),spbr(kj,ke,kb),                 &
                   pprb(kj,ke,kb),pprbsav(kj,ke,kb),                 &
                   pprbav(kj,ke,kb),spbolr(kj,ke,kb)                 
!
! --- source arrays for ion energy and angular momentum equations
!
      real *8                                                        &
                       sprbeam(kj),ssprcxl(kj),sprcxre(kj),          &
                       sprcx(kj),spreimpt(kj),sprbeame(kj),          &
                       sprbeami(kj),spr2d(kj),omegale(kj),           &
                       sprcxree(kj),sprcxe(kj),spreimpe(kj),         &
                       spbolt(kj),smagtorque(kj),sprbeame_intg,       &
                       sprbeami_intg ,storqueb_intg                  
                       !
! --- physics and miscellaneous quantities
!
      real *8                                                        &
                     vionz(kj),storque(kj),dangrot(kj),angmtm(kj),   &
                     angmtot,tauang(kj),tauangt,storquet,            &
                     xkangrot(kj),angrot(kj),storqueb(kj),beamtorq,  &
                     dangmtot,angmtotn,amtinrta(kj),totinrta ,       &
                     avgangmt(kj),avgangdt(kj),xkangrob,xkangroc,    &
                     rotenergy(kj),avrotjou,omegapi(kj),qomegapi(kj),&
                     vischeat(kj),omegdgam(kj),pvscheat(kj),         &
                     pomdgam(kj),fluxangv(kj),fluxangc(kj),          &
                     flxangce(kj),qangc(kj),qangv(kj),qangce(kj),    &
                     pwdnidt(kj),aniwdwdt(kj),pniwdwdt(kj),          &
                     psprcxe(kj),psprcxee(kj),pspreimp(kj),          &
                     vionzgrd(kj),wdnidt(kj),pqangce(kj),            &
                     dkapomeg(kj),rgc_mult(kj),angmtoto
                    !
! --- NTV momentum source term variables    C.K.Pan  03/01/2010
      integer mp_tornum,mp_polnum
      integer include_ntv
      real*8 c_p
      real*8 delta_b_sqr,bmn_theta_a_sqr,delta_b_o_b_sqr
      real*8 xlam,vth2,drdpsi,ddeltadpsi,dtidpsi,dpidpsi,dnedpsi,   &
            omega_gb,ntvtorquet,omegae,esupc
         
      real*8 xnu_ntv(kj),mu_ntv_p(kj),omegazero(kj),omegas(kj),     &                              
            angrot_ntv(kj),delta_ntv(kj),                           &
            sntvtorque(kj),sntvtorque_int(kj),                      &
            rntv(kj),tierg(kj)

      data  esupc    / 2.99792458e+9 /    ! Transforming from SI to 
                                           ! gaussian for Charge
    END   MODULE tordlrot
