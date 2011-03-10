      SUBROUTINE FORCEBAL(igrid_ncl,jmaxm,xneo,device,idshot,time,
     &                    ci0,cim1,cim2,nr_xptor,rhot_xptor,
     &                    rhop_r,q_r,frb_r,rm1_r,rm2_r,vol_r,iflag)
!***********************************************************************f
!FORCEBAL reads MHD equilibrium and profile data then uses it as input 
!  to NCLASS to calculate neoclassical transport properties, poloidal
!  rotation velocities and the radial electric field using the parallel
!  and radial force balance equations
!  Equilibrium option presently limited to EFIT data files
!  Profile option presently limited to 4D data files
!  W.A. Houlberg 3/2000
! 6/28/00 jek multiple changes:
!         changed PROGRAM to SUBROUTINE
!         changed STOP to RETURN
!         added nr_r,rhot_r to argument list
!         added GET_PRO_XPTOR, denfst_r(mxnr_r)
!***********************************************************************
c      IMPLICIT NONE
      include 'mpif.h'
      include '../inc/glf.m'
      INCLUDE '../inc/pamx_ms.inc'
      INCLUDE '../inc/pamx_mz.inc'
!Dimensions for radial parameters
      INTEGER        jmaxm
      INTEGER        mxnr_r,                   mx_ni,
     &               mxnpro
      PARAMETER     (mxnr_r=300,               mx_ni=5,
     &               mxnpro=250)
!FORCEBAL output
      INCLUDE '../inc/comfbl.inc'
!Declaration of functions
      INTEGER        IARRAY_MAX
!Declaration of I/O options, units and names
      CHARACTER*15   device,                  idtransp,
     &               idmhd
      CHARACTER*3    idrun
      CHARACTER*25   cn1d,                    cnout,
     &               cnqr,                    cnsum,
     &               cnwout
      CHARACTER*6    uid(6)
      CHARACTER*4    dda(6)
      INTEGER        iseqp(6)
      INTEGER        idshot,                  iyear,
     &               k_edotb,                 k_potato,
     &               k_sqz,                   n1d,
     &               nout,                    npro,
     &               nr_r,                    nsum
!Declaration of radial variables for 1D data file
      CHARACTER*15   namepro(mxnpro),         unitpro(mxnpro)
      REAL           valpro(mxnr_r,mxnpro)
!Declaration of geometry and plasma variables
!  0-D
      INTEGER        igrid_ncl,               nmhd_r        
      INTEGER        izi0(mx_ni),             izim1,
     &               izim2                    
      REAL           amui0(mx_ni),            amuim1,
     &               amuim2
      REAL           a0,                      bt0,
     &               c_den,                   r0,
     &               dt,                      time
!  Radial geometric variables
      REAL           b2_r(mxnr_r),            bm2_r(mxnr_r),
     &               bpout_r(mxnr_r),         btout_r(mxnr_r),
     &               elong_r(mxnr_r),         f_r(mxnr_r),
     &               frb_r(mxnr_r),           rm1_r(mxnr_r),
     &               rm2_r(mxnr_r),
     &               fhat_r(mxnr_r),          fm_r(3,mxnr_r),
     &               ftrap_r(mxnr_r),         gph_r(mxnr_r),
     &               gr2bm2_r(mxnr_r),        grho1_r(mxnr_r),
     &               grho2_r(mxnr_r),         grth_r(mxnr_r),      
     &               gth_r(mxnr_r),           phit_r(mxnr_r),   
     &               psi_r(mxnr_r),           psip_r(mxnr_r),                     
     &               q_r(mxnr_r),             rhop_r(mxnr_r),
     &               rhot_r(mxnr_r),          rin_r(mxnr_r),
     &               rout_r(mxnr_r),          vol_r(mxnr_r),
     &               vp_r(mxnr_r),            rhor_r(mxnr_r)
      REAL           b2_r_new(mxnr_r),        bm2_r_new(mxnr_r),
     &               bpout_r_new(mxnr_r),     btout_r_new(mxnr_r),
     &               elong_r_new(mxnr_r),     f_r_new(mxnr_r),
     &               frb_r_new(mxnr_r),
     &               fhat_r_new(mxnr_r),      fm_r_new(3,mxnr_r),
     &               ftrap_r_new(mxnr_r),     gph_r_new(mxnr_r),
     &               gr2bm2_r_new(mxnr_r),    grho1_r_new(mxnr_r),
     &               grho2_r_new(mxnr_r),     grth_r_new(mxnr_r),      
     &               gth_r_new(mxnr_r),       phit_r_new(mxnr_r),   
     &               psi_r_new(mxnr_r),       psip_r_new(mxnr_r),                     
     &               q_r_new(mxnr_r),         rhop_r_new(mxnr_r),
     &               rhot_r_new(mxnr_r),      rin_r_new(mxnr_r),
     &               rout_r_new(mxnr_r),      vol_r_new(mxnr_r),
     &               vp_r_new(mxnr_r)    
!  Radial plasma variables  
      REAL           rho_r(mxnr_r)
      REAL           cur_r(mxnr_r),           cur_bs_r(mxnr_r),
     &               cur_ex_r(mxnr_r),        cur_nb_ex_r(mxnr_r),
     &               den_r(mxnr_r),           deni_r(mxnr_r,mx_ni),
     &               denim_r(mxnr_r),         e_par_ex_r(mxnr_r),
     &               omega_ex_r(mxnr_r),      te_r(mxnr_r),
     &               ti_r(mxnr_r),            vp_im_ex_r(mxnr_r),
     &               vt_im_ex_r(mxnr_r),      xdum_r(mxnr_r),
     &               xj_r(mxnr_r),
     &               xj_nb_ex_r(mxnr_r),      xk_im_ex_r(mxnr_r),
     &               zeff_ex_r(mxnr_r),       zeff_r(mxnr_r)
      REAL           denim2_r(mxnr_r,mx_mz),  denfst_r(mxnr_r)
!Declaration of other local variables
      LOGICAL        ki0on(mx_ni),            kim1on,
     &               kim2on(mx_mz)
      CHARACTER*2    cn(mx_mz)
      CHARACTER*3    ci0(mx_ni),              cim1,
     &               cim2
      CHARACTER*5    ctimems
      CHARACTER*9    cidshot
      CHARACTER*120  message
      INTEGER        i,                       icid,
     &               ici0(mx_ni),             icim1,
     &               icim2,                   iflag,
     &               itimems,                 j,
     &               j1,                      jimp,
     &               k,                       k1,
     &               npropr,                  nrpr
      CHARACTER*15   cdum(6)
      INTEGER        idum(2)
      REAL           rdum(5)
      REAL           xidshot
      REAL           z_mu0,                   z_pi
c
      integer xneo, nr_xptor
      real rhot_xptor(mxnr_r)
      save nr_r, rhot_r, rho_r
c
!Namelist
      NAMELIST/inforce/device,idshot,iyear,idrun,cnwout,cnqr,time,dt,
     &                 a0,r0,bt0,
     &                 ci0,cim1,cim2,
     &                 idtransp,idmhd,uid,dda,iseqp,
     &                 k_edotb,k_potato,k_sqz
!Initialization
!  Warning and error flag and message
      iflag=0
      message=' '
!  Physical and conversion constants
      z_mu0=1.2566e-06
      z_pi=ACOS(-1.0)
!Set default input
c     device='diiid'
c     idshot=10000
      iyear=1999
c     DO i=1,mx_ni
c       ci0(i)='   '
c     ENDDO
c     cim1='   '
c     cim2='   '
      idrun='a01'
      IF(device.eq.'cmod') THEN      
!       CMOD
        device='cmod'
        idshot=960116024
        idmhd='efit'
        idrun='a01'
        time=0.920
        idtransp='transp17'
        k_edotb=0
        k_potato=0
        k_sqz=1
      ELSEIF(device.eq.'diiid') THEN
!       DIII-D
        device='diiid'
        idmhd='efit'
        idshot=84736
        idrun='a01'
        time=1.300
        k_edotb=0
        k_potato=0
        k_sqz=1
      ELSEIF(device.eq.'jet') THEN
!       JET
        device='jet'
        idshot=34476
        idrun='a01'
        time=52.913
!       Te profile source
        uid(1)='JETPPF'
        dda(1)='LIDR'
        iseqp(1)=0
!       Ti profile source
        uid(2)='JETPPF'
        dda(2)='CXSM'
        iseqp(2)=0
!       ne profile source
        uid(3)='JETPPF'
        dda(3)='LIDR'
        iseqp(3)=0
!       Impurity (1) density - single charge state
        uid(4)='JETPPF'
        dda(4)='CXSM'
        iseqp(4)=0
!       Impurity (1) toroidal rotation - single charge state
        uid(5)='JETPPF'
        dda(5)='CXSM'
        iseqp(5)=0
!       Impurity (2) density - multiple charge states
        uid(6)='JETRGI'
        dda(6)='NI  '
        iseqp(6)=0
        k_edotb=0
        k_potato=0
        k_sqz=1
      ELSEIF(device.eq.'start') THEN
!       START
        device='start'
        idshot=10000
        idrun='a01'
        time=0.5
        k_edotb=0
        k_potato=0
        k_sqz=1
      ELSEIF(device.eq.'textor') THEN
!       TEXTOR
        device='textor'
        idmhd='efit'
        idshot=70287
        idrun='a01'
        time=1.5
        k_edotb=0
        k_potato=0
        k_sqz=1
      ELSEIF(device.eq.'toresupra') THEN
!       TORE SUPRA
        device='toresupra'
        idmhd='circ'
        idshot=11111
        idrun='a01'
        cnqr='qrprof.dat'
        time=1.5
        a0=0.75
        r0=2.3
        bt0=2.8
        k_edotb=0
        k_potato=0
        k_sqz=1
      ELSEIF(device.eq.'tftr') THEN
!       TFTR
        device='tftr'
        idshot=94497
        iyear=1996
        idrun='a01'
        idtransp='t12'
        idmhd='vmec'
        cnwout='wout.104977d01t201'
        time=2.6
        dt=0.05
        k_edotb=0
        k_potato=0
        k_sqz=1
      ELSEIF(device.eq.'vmec') THEN
!       VMEC test
        device='vmec'
        idshot=10000
        idrun='a01'
        idmhd='vmec'
        cnwout='wout.sympp'
        ci0(1)='H'
        cim1=' '
        cim2=' '
        time=0.0
        k_edotb=0
        k_potato=0
        k_sqz=0
      ELSEIF(device.eq.'xptor') THEN
        idmhd='efit'
c        ci0(1)='D'
c        ci0(2)='N'
c        cim1='C'
        k_edotb=0
        k_potato=0
        k_sqz=0
      ENDIF
c
c     write(*,*) 'xneo = ',xneo
c     write(*,*) 'device = ',device
c     write(*,*) 'idshot = ',idshot
c     write(*,*) 'time = ',time
c     write(*,*) 'ci0 = ',ci0(1),ci0(2)
c     write(*,*) 'cim1 = ',cim1
c     write(*,*) 'cim2 = ',cim2
c
!Open namelist input file and read 
      if (device.ne.'xptor') then
      OPEN(unit=1,file='forcebal.nml',status='old',form='formatted')
      READ(1,inforce)
      CLOSE(unit=1)
      endif
!Set file names and device dependent options
      itimems=1.0e3*(time+.00001)
      WRITE(ctimems,'(i5)') itimems
      DO i=1,5
        IF(ctimems(i:i).eq.' ') ctimems(i:i)='0'
      ENDDO
      WRITE(cidshot,'(i9)') idshot
      xidshot=idshot
      icid=9-INT(LOG10(xidshot))
      n1d=10
      cn1d='1d_'//cidshot(icid:9)//'_'//ctimems//'_'//idrun//'.fbl'
      nsum=11
      cnsum='sum_'//cidshot(icid:9)//'_'//ctimems//'_'//idrun//'.fbl'
      nout=12
      cnout='msg_'//cidshot(icid:9)//'_'//ctimems//'_'//idrun//'.fbl'
!Open message output files
      OPEN(unit=nout,
     &     file=cnout,
     &     status='unknown',
     &     form='formatted')
!Set cutoff density
      c_den=1.0e10
!Set species charges and masses for species in namelist
!  Fully ionized species
      DO i=1,mx_ni
        CALL SPECIES(ci0(i),ici0(i),izi0(i),amui0(i),iflag,message)
!       Error checks
        IF(iflag.ne.0) THEN
          j=LEN(message)
          message=message(1:j-12)//', name = '//ci0(i)
          CALL WRITE_LINE(nout,message,0,0)
          IF(iflag.eq.1) GOTO 1000
          iflag=0
          message=' '
        ENDIF
        IF(izi0(i).gt.mx_mz) THEN
          message='ERROR:max charge state exceeded, increase mx_mz'
          CALL WRITE_LINE(nout,message,0,0)
          message='  name = '//ci0(i)
          GOTO 1000
        ENDIF
      ENDDO
!  Diagnostic specie
      CALL SPECIES(cim1,icim1,izim1,amuim1,iflag,message)
!     Error checks
      IF(iflag.ne.0) THEN
        CALL WRITE_LINE(nout,message,0,0)
        IF(iflag.eq.1) THEN
          message='  name = '//cim1
          CALL WRITE_LINE(nout,message,0,0)
          GOTO 1000
        ENDIF
        iflag=0
        message=' '
      ENDIF
      IF(izim1.gt.mx_mz) THEN
        message='ERROR:max charge state exceeded, increase mx_mz'
        CALL WRITE_LINE(nout,message,0,0)
        message='  name = '//cim1
        GOTO 1000
      ENDIF
!  Multiple charge state specie
      CALL SPECIES(cim2,icim2,izim2,amuim2,iflag,message)
!     Error checks
      IF(iflag.ne.0) THEN
        CALL WRITE_LINE(nout,message,0,0)
        IF(iflag.eq.1) THEN
          message='  name = '//cim2
          CALL WRITE_LINE(nout,message,0,0)
          GOTO 1000
        ENDIF
        iflag=0
        message=' '
      ENDIF
      IF(izim2.gt.mx_mz) THEN
        message='ERROR:max charge state exceeded, increase mx_mz'
        CALL WRITE_LINE(nout,message,0,0)
        message='  name = '//cim2
        GOTO 1000
      ENDIF
!Get geometry information
      iflag=0
      if(xneo.eq.0) then  ! only read geometry once from readxp.f
        if(i_proc.eq.0) 
     &     write(*,*) 'Reading EFIT geometry for NCLASS ...'
c
      IF(idmhd.eq.'vmec') THEN
        nr_r=31
        DO i=1,nr_r
          rhot_r(i)=FLOAT(i-1)/FLOAT(nr_r-1)
        ENDDO
        CALL GEOM_VMEC(cnwout,nr_r,rhot_r,rhop_r,f_r,q_r,rin_r,rout_r,
     &                 elong_r,vol_r,vp_r,phit_r,fm_r,grth_r,gph_r,
     &                 gth_r,b2_r,bm2_r,grho1_r,grho2_r,gr2bm2_r,
     &                 ftrap_r,fhat_r,bpout_r,btout_r,psi_r,a0,r0,bt0,
     &                 iflag,message)
      ELSEIF(idmhd.eq.'circ') THEN
        nr_r=31
        DO i=1,nr_r
          rhot_r(i)=FLOAT(i-1)/FLOAT(nr_r-1)
        ENDDO
        CALL GEOM_CIRC(cnqr,nr_r,a0,r0,bt0,rhot_r,rhop_r,f_r,q_r,rin_r,
     &                 rout_r,elong_r,vol_r,vp_r,phit_r,fm_r,grth_r,
     &                 gph_r,gth_r,b2_r,bm2_r,grho1_r,grho2_r,gr2bm2_r,
     &                 ftrap_r,fhat_r,bpout_r,btout_r,psi_r)
      ELSEIF(idmhd.eq.'xptor') THEN
        CALL GEOM_EFIT('efit',idshot,time,nr_r,rhop_r,f_r,q_r,rhot_r,
     &                 rin_r,rout_r,elong_r,vol_r,vp_r,phit_r,fm_r,
     &                 grth_r,gph_r,gth_r,b2_r,bm2_r,grho1_r,grho2_r,
     &                 gr2bm2_r,ftrap_r,fhat_r,bpout_r,btout_r,psi_r,a0,
     &                 r0,bt0,rm1_r,rm2_r,rhor_r,iflag,message)
      ELSE
        CALL GEOM_EFIT(device,idshot,time,nr_r,rhop_r,f_r,q_r,rhot_r,
     &                 rin_r,rout_r,elong_r,vol_r,vp_r,phit_r,fm_r,
     &                 grth_r,gph_r,gth_r,b2_r,bm2_r,grho1_r,grho2_r,
     &                 gr2bm2_r,ftrap_r,fhat_r,bpout_r,btout_r,psi_r,a0,
     &                 r0,bt0,rm1_r,rm2_r,rhor_r,iflag,message)
      ENDIF
c
      write(*,*) 'device = ',device
      write(*,*) 'jmaxm-forcebal = ',jmaxm
      write(*,*) 'igrid_ncl = ',igrid_ncl
      nmhd_r=nr_r
      if(igrid_ncl.eq.1) nr_r=jmaxm+1
      nr_xptor=nr_r
      do i=1,nr_r
        rhot_xptor(i)=rhot_r(i)
        frb_r(i)=f_r(i)*z_mu0/(2.0*z_pi)
      enddo
      endif
c
c     write(*,*) 'iflag = ',iflag
c     do i=1,41
c       write(*,100) i, rhot_r(i),rhop_r(i),q_r(i)
c     enddo
c
      IF(iflag.eq.1) GOTO 1000
      IF(iflag.eq.-1) THEN
        CALL WRITE_LINE(nout,message,0,0)
        iflag=0
        message=' '
      ENDIF
!Get plasma profile information
      IF(device.eq.'jet') THEN
c--        CALL READ_PRO_JET(idshot,time,uid,dda,iseqp,izi0,izim1,izim2,
c--     &                    nr_r,mxnr_r,rin_r,rout_r,te_r,ti_r,den_r,
c--     &                    deni_r,denim_r,denim2_r,vt_im_ex_r,zeff_r,
c--     &                    iflag,message)
      ELSEIF(device.eq.'start') THEN
        CALL READ_PRO_START(idshot,time,izi0,izim1,izim2,nr_r,mxnr_r,
     &                      rout_r,te_r,ti_r,den_r,deni_r,denim_r,
     &                      denim2_r,vt_im_ex_r,zeff_r,iflag,message)
      ELSEIF(device.eq.'textor') THEN
        CALL READ_PRO_TEXTOR(idshot,time,izi0,izim1,izim2,nr_r,mxnr_r,
     &                       rout_r,te_r,ti_r,den_r,deni_r,denim_r,
     &                       denim2_r,vt_im_ex_r,zeff_r,iflag,message)
      ELSEIF(device.eq.'tftr') THEN
c--        CALL READ_PRO_TFTR(idshot,iyear,idtransp,time,dt,ci0,ici0,izi0,
c--     &                     amui0,mx_ni,cim1,icim1,izim1,amuim1,nr_r,
c--     &                     mxnr_r,rhot_r,rout_r,te_r,ti_r,den_r,deni_r,
c--     &                     denim_r,omega_ex_r,vt_im_ex_r,zeff_r,iflag,
c--     &                     message)
      ELSEIF(device.eq.'vmec') THEN
        CALL READ_PRO_VMEC(nr_r,mxnr_r,rhot_r,te_r,ti_r,den_r,deni_r,
     &                     denim_r,denim2_r,zeff_r,iflag,message)
      ELSEIF(device.eq.'xptor') THEN
        CALL GET_PRO_XPTOR(nout,idshot,igrid_ncl,
     #                   time,bt0,izi0,izim1,izim2,nr_r,mxnr_r,
     #                   rho_r,rhot_r,te_r,ti_r,den_r,deni_r,
     #                   zeff_ex_r,denim_r,denim2_r,denfst_r,omega_ex_r,
     #                   vp_im_ex_r,vt_im_ex_r,xk_im_ex_r,zeff_r,
     #                   e_par_ex_r,xj_ex_r,xj_nb_ex_r,iflag)
      ELSE
        CALL READ_PRO_DIIID(idshot,time,bt0,izi0,izim1,izim2,nr_r,
     &                      mxnr_r,rhot_r,te_r,ti_r,den_r,deni_r,
     &                      zeff_ex_r,denim_r,denim2_r,omega_ex_r,
     &                      vp_im_ex_r,vt_im_ex_r,xk_im_ex_r,zeff_r,
     &                      e_par_ex_r,xj_ex_r,xj_nb_ex_r,iflag,
     &                      message)
      ENDIF
c
      if(igrid_ncl.eq.1) then
      CALL RARRAY_ZERO(nr_r,rhop_r_new)
      CALL RARRAY_ZERO(nr_r,f_r_new)
      CALL RARRAY_ZERO(nr_r,q_r_new)
      CALL RARRAY_ZERO(nr_r,rin_r_new)
      CALL RARRAY_ZERO(nr_r,rout_r_new)
      CALL RARRAY_ZERO(nr_r,elong_r_new)
      CALL RARRAY_ZERO(nr_r,vol_r_new)
      CALL RARRAY_ZERO(nr_r,vp_r_new)
      CALL RARRAY_ZERO(nr_r,phit_r_new)
      CALL RARRAY_ZERO(nr_r,fm_r_new)
c     write(*,*) 'nr_r = ',nr_r
c     write(*,*) 'nmhd_r = ',nmhd_r
c     write(*,*) 'before interpolation'
c     do i=1,nr_r
c       write(*,100) i, rhot_r(i), rhop_r(i), q_r(i)
c     enddo
        call w_lin_interp(nmhd_r,rhot_r,rhop_r,nr_r,rho_r,rhop_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,f_r,nr_r,rho_r,f_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,q_r,nr_r,rho_r,q_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,rin_r,nr_r,rho_r,rin_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,rout_r,nr_r,rho_r,rout_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,elong_r,nr_r,rho_r,elong_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,vol_r,nr_r,rho_r,vol_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,vp_r,nr_r,rho_r,vp_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,phit_r,nr_r,rho_r,phit_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,fm_r(1,:),nr_r,rho_r,
     &     fm_r_new(1,:),iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,fm_r(2,:),nr_r,rho_r,
     &     fm_r_new(2,:),iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,fm_r(3,:),nr_r,rho_r,
     &     fm_r_new(3,:),iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,grth_r,nr_r,rho_r,grth_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,gph_r,nr_r,rho_r,gph_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,gth_r,nr_r,rho_r,gth_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,b2_r,nr_r,rho_r,b2_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,bm2_r,nr_r,rho_r,bm2_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,grho1_r,nr_r,rho_r,grho1_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,grho2_r,nr_r,rho_r,grho2_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,gr2bm2_r,nr_r,rho_r,
     &     gr2bm2_r_new,iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,ftrap_r,nr_r,rho_r,ftrap_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,bpout_r,nr_r,rho_r,bpout_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,btout_r,nr_r,rho_r,btout_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,psi_r,nr_r,rho_r,psi_r_new,
     &     iflag,msg)
        call w_lin_interp(nmhd_r,rhot_r,fhat_r,nr_r,rho_r,fhat_r_new,
     &     iflag,msg)
c
c     write(*,*) 'after interpolation'
      do i=1,nr_r
        rhot_r(i)=rho_r(i)
        rhot_xptor(i)=rho_r(i)
        rhop_r(i)=rhop_r_new(i)
        f_r(i)=f_r_new(i)
        q_r(i)=q_r_new(i)
        rin_r(i)=rin_r_new(i)
        rout_r(i)=rout_r_new(i)
        elong_r(i)=vol_r_new(i)
        vp_r(i)=vp_r_new(i)
        phit_r(i)=phit_r_new(i)
        fm_r(1,i)=fm_r_new(1,i)
        fm_r(2,i)=fm_r_new(2,i)
        fm_r(3,i)=fm_r_new(3,i)
        grth_r(i)=grth_r_new(i)
        gph_r(i)=gph_r_new(i)
        gth_r(i)=gth_r_new(i)
        b2_r(i)=b2_r_new(i)
        bm2_r(i)=bm2_r_new(i)
        grho1_r(i)=grho1_r_new(i)
        grho2_r(i)=grho2_r_new(i)
        gr2bm2_r(i)=gr2bm2_r_new(i)
        ftrap_r(i)=ftrap_r_new(i)
        bpout_r(i)=bpout_r_new(i)
        btout_r(i)=btout_r_new(i)
        psi_r(i)=psi_r_new(i)
        fhat_r(i)=fhat_r_new(i)
        frb_r(i)=f_r(i)*z_mu0/(2.0*z_pi)
c       write(*,100) i, rho_r(i), grho1_r(i)
      enddo
      endif
c
      IF(iflag.eq.1) GOTO 1000
      IF(iflag.eq.-1) THEN
        CALL WRITE_LINE(nout,message,0,0)
        iflag=0
        message=' '
      ENDIF
!Set Psi'
      psip_r(1)=0.0
      DO i=2,nr_r
        psip_r(i)=z_mu0*f_r(i)/fhat_r(i)
      ENDDO
c
c
!Set integrated currents
      cur_r(1)=0.0
      cur_ex_r(1)=0.0
      cur_nb_ex_r(1)=0.0
      DO i=2,nr_r
        cur_r(i)=gth_r(i)*gph_r(i)/q_r(i)*f_r(i)
c       write(*,100) i, rhot_r(i), gth_r(i),
c    &               gph_r(i),f_r(i),q_r(i),cur_r(i)
        cur_ex_r(i)=f_r(i)*(cur_ex_r(i-1)/f_r(i-1)
     &              +a0*bt0*(rhot_r(i)-rhot_r(i-1))/z_mu0
     &              *0.5*(vp_r(i)  *xj_ex_r(i)  /f_r(i)**2
     &                   +vp_r(i-1)*xj_ex_r(i-1)/f_r(i-1)**2))
        cur_nb_ex_r(i)=f_r(i)*(cur_nb_ex_r(i-1)/f_r(i-1)
     &                 +a0*bt0*(rhot_r(i)-rhot_r(i-1))/z_mu0
     &                 *0.5*(vp_r(i)  *xj_nb_ex_r(i)  /f_r(i)**2
     &                      +vp_r(i-1)*xj_nb_ex_r(i-1)/f_r(i-1)**2))
      ENDDO
!  Reset total current using <J.B> from second to last node
      cur_r(nr_r)=f_r(nr_r)*(cur_r(nr_r-1)/f_r(nr_r-1)
     &            +(f_r(nr_r-1)/f_r(nr_r))**2
     &            *vp_r(nr_r)/vp_r(nr_r-1)
     &            *(rhot_r(nr_r)-rhot_r(nr_r-1))
     &            /(rhot_r(nr_r-1)-rhot_r(nr_r-2))
     &            *(cur_r(nr_r-1)/f_r(nr_r-1)
     &            -cur_r(nr_r-2)/f_r(nr_r-2)))
!  Total current density from equilibrium
      DO i=2,nr_r-1
        xj_r(i)=f_r(i)**2*z_mu0/vp_r(i)/bt0
     &          *(cur_r(i+1)/f_r(i+1)-cur_r(i-1)/f_r(i-1))
     &          /a0/(rhot_r(i+1)-rhot_r(i-1))
      ENDDO
      xj_r(1)=xj_r(2)
      xj_r(nr_r)=xj_r(nr_r-1)
c
c set jgrid=0 in readxp.f for no interpolation 
c to transport grid. Also, set nr_r=41 in 
c read_efit_eqdsk.f for testing
c Note: mx_ni set to 5 in efitdat.m
      call store_efit(nr_r,k_edotb,k_potato,
     &                a0,r0,bt0,
     &                izim1,izim2,amuim1,amium2,
     &                izi0,amui0,rhop_r,rhot_r,
     &                f_r,q_r,xj_r,xj_nb_ex_r,
     &                rin_r,rout_r,elong_r,vol_r,vp_r,
     &                phit_r,fm_r,grth_r,gph_r,gth_r,
     &                b2_r,bm2_r,grho1_r,grho2_r,
     &                gr2bm2_r,ftrap_r,fhat_r,
     &                bpout_r,btout_r,psi_r,psip_r,
     &                rm1_r,rm2_r)
c
c...jek diagnostic printout
      do i=1,nr_r
c        write(*,*) i, rhot_r(i), xj_r(i),'forcebal'
c        write(*,*) i, rhot_r(i)*a0, 'forcebal'
c        write(*,*) i, deni_r(i,1),deni_r(i,2)
      enddo
c
! jek diagnostic printout of metrics
!
!     DO i=1,nr_r
!        write(*,100) i, rhot_r(i), q_r(i), vp_r(i), f_r(i), 
!    &                gph_r(i), gth_r(i), xj_r(i), cur_r(i)
!        write(*,100) i, rhot_r(i), b2_r(i), bm2_r(i), fhat_r(i),
!    &        grth_r(i), gr2bm2_r(i), fm_r(1,i), fm_r(2,i), fm_r(3,i)
!     ENDDO
 100  format(i2,2x,0p1f8.6,1p10e13.6)
!Get rotation and transport profiles
      CALL NCLASS_DR(k_edotb,k_potato,k_sqz,amui0,izi0,amuim1,izim1,
     &               amuim2,izim2,c_den,a0,bt0,q_r(1),elong_r(1),nr_r,
     &               b2_r,bm2_r,fhat_r,ftrap_r,fm_r,grth_r,grho2_r,
     &               gr2bm2_r,rhot_r,den_r,deni_r,denim_r,denim2_r,te_r,
     &               ti_r,rout_r,btout_r,bpout_r,vp_im_ex_r,vt_im_ex_r,
     &               psi_r,psip_r,e_par_ex_r,xj_r,xj_nb_ex_r,iflag,
     &               message)
      IF(iflag.eq.1) GOTO 1000
      IF(iflag.eq.-1) THEN
        CALL WRITE_LINE(nout,message,0,0)
        iflag=0
        message=' '
      ENDIF
!  Integrated bootstrap current
      cur_bs_r(1)=0.0
      DO i=2,nr_r
        cur_bs_r(i)=f_r(i)*(cur_bs_r(i-1)/f_r(i-1)
     &              +bt0/z_mu0*a0*(rhot_r(i)-rhot_r(i-1))*0.5
     &              *(vp_r(i)  *xj_bs_r(i)  /f_r(i)**2
     &               +vp_r(i-1)*xj_bs_r(i-1)/f_r(i-1)**2))
      ENDDO
!Identify whether species are present
      DO j=1,mx_ni
        IF(deni_r(1,j).gt.c_den) THEN
          ki0on(j)=.true.
        ELSE
          ki0on(j)=.false.
        ENDIF
      ENDDO
      IF(denim_r(1).gt.c_den) THEN
        kim1on=.true.
      ELSE
        kim1on=.false.
      ENDIF
      DO j=1,mx_mz
        kim2on(j)=.false.
      ENDDO
      IF(izim2.gt.0) THEN
        DO j=1,izim2
          IF(denim2_r(1,j).gt.c_den) THEN
            kim2on(j)=.true.
            WRITE(cn(j),'(i2)') j
          ENDIF
        ENDDO
      ENDIF
!Load arrays for profile variables
      npro=0
!  Radial grids
      npro=npro+1
      namepro(npro)='rho_p_r'
      unitpro(npro)='-'
      CALL RARRAY_COPY(nr_r,rhop_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='rho_t_r'
      unitpro(npro)='-'
      CALL RARRAY_COPY(nr_r,rhot_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='r_p_r'
      unitpro(npro)='m'
      CALL RARRAY_COPY(nr_r,rhop_r,1,valpro(1,npro),1)
      CALL RARRAY_SCALE(nr_r,a0,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='r_t_r'
      unitpro(npro)='m'
      CALL RARRAY_COPY(nr_r,rhot_r,1,valpro(1,npro),1)
      CALL RARRAY_SCALE(nr_r,a0,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='R_in_r'
      unitpro(npro)='m'
      CALL RARRAY_COPY(nr_r,rin_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='R_o_r'
      unitpro(npro)='m'
      CALL RARRAY_COPY(nr_r,rout_r,1,valpro(1,npro),1)
!  Magnetic fields and fluxes
      npro=npro+1
      namepro(npro)='B_p_o_r'
      unitpro(npro)='T'
      CALL RARRAY_COPY(nr_r,bpout_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='B_t_o_r'
      unitpro(npro)='T'
      CALL RARRAY_COPY(nr_r,btout_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='Phi_t_r'
      unitpro(npro)='Wb'
      CALL RARRAY_COPY(nr_r,phit_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='Psi_r'
      unitpro(npro)='Wb/rad'
      CALL RARRAY_COPY(nr_r,psi_r,1,valpro(1,npro),1)
      npro=npro+1
      DO i=1,nr_r
        xdum_r(i)=ABS(q_r(i))
      ENDDO
      namepro(npro)='q_r'
      unitpro(npro)='-'
      CALL RARRAY_COPY(nr_r,xdum_r,1,valpro(1,npro),1)
!  Flux functions and metrics
      npro=npro+1
      namepro(npro)='F_r'
      unitpro(npro)='A'
      CALL RARRAY_COPY(nr_r,f_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='f_trap_r'
      unitpro(npro)='-'
      CALL RARRAY_COPY(nr_r,ftrap_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='grad(r_t)**2_r'
      unitpro(npro)='-'
      CALL RARRAY_COPY(nr_r,grho2_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='dV/d(r_t)_r'
      unitpro(npro)='m**2'
      CALL RARRAY_COPY(nr_r,vp_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='Vol_r'
      unitpro(npro)='m**3'
      CALL RARRAY_COPY(nr_r,vol_r,1,valpro(1,npro),1)
!  Plasma profiles
      npro=npro+1
      namepro(npro)='den_e_r'
      unitpro(npro)='/m**3'
      CALL RARRAY_COPY(nr_r,den_r,1,valpro(1,npro),1)
      DO j=1,mx_ni
        IF(ki0on(j)) THEN
          npro=npro+1
          namepro(npro)='den_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='/m**3'
          CALL RARRAY_COPY(nr_r,deni_r(1,j),1,valpro(1,npro),1)
        ENDIF
      ENDDO
      IF(kim1on) THEN
        npro=npro+1
        namepro(npro)='den_'//cim1(:icim1)//'_r'
        unitpro(npro)='/m**3'
        CALL RARRAY_COPY(nr_r,denim_r,1,valpro(1,npro),1)
      ENDIF
      DO j=1,izim2
        IF(kim2on(j)) THEN
          npro=npro+1
          namepro(npro)='den_'//cim2(:icim2)//'('//cn(j)//')_r'
          unitpro(npro)='/m**3'
          CALL RARRAY_COPY(nr_r,denim2_r(1,j),1,valpro(1,npro),1)
        ENDIF
      ENDDO
      npro=npro+1
      namepro(npro)='T_e_r'
      unitpro(npro)='keV'
      CALL RARRAY_COPY(nr_r,te_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='T_i_r'
      unitpro(npro)='keV'
      CALL RARRAY_COPY(nr_r,ti_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='prs_e_r'
      unitpro(npro)='keV/m**3'
      DO i=1,nr_r
        valpro(i,npro)=den_r(i)*te_r(i)
      ENDDO
      DO j=1,mx_ni
        IF(ki0on(j)) THEN
          npro=npro+1
          namepro(npro)='prs_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='keV/m**3'
          DO i=1,nr_r
            valpro(i,npro)=deni_r(i,j)*ti_r(i)
          ENDDO
        ENDIF
      ENDDO
      IF(kim1on) THEN
        npro=npro+1
        namepro(npro)='prs_'//cim1(:icim1)//'_r'
        unitpro(npro)='keV/m**3'
        DO i=1,nr_r
          valpro(i,npro)=denim_r(i)*ti_r(i)
        ENDDO
      ENDIF
      DO j=1,izim2
        IF(kim2on(j)) THEN
          npro=npro+1
          namepro(npro)='prs_'//cim2(:icim2)//'('//cn(j)//')_r'
          unitpro(npro)='keV/m**3'
          DO i=1,nr_r
            valpro(i,npro)=denim2_r(i,j)*ti_r(i)
          ENDDO
        ENDIF
      ENDDO
      i=IARRAY_MAX(nr_r,zeff_ex_r,1)
      IF(zeff_ex_r(i).ge.1.0) THEN
        npro=npro+1
        namepro(npro)='Z_eff_ex_r'
        unitpro(npro)='-'
        CALL RARRAY_COPY(nr_r,zeff_ex_r,1,valpro(1,npro),1)
      ENDIF
      npro=npro+1
      namepro(npro)='Z_eff_r'
      unitpro(npro)='-'
      CALL RARRAY_COPY(nr_r,zeff_r,1,valpro(1,npro),1)
!  Particle fluxes, diffusivities and pinches
!  g factor - first index is gradient and second is flux
!     Electrons
      npro=npro+1
      namepro(npro)='Gam_e_r'
      unitpro(npro)='#/m**2/s'
      CALL RARRAY_COPY(nr_r,gam_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='D_eff_e_r'
      unitpro(npro)='m**2/s'
      CALL RARRAY_COPY(nr_r,d_eff_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='D_eff_g_e_r'
      unitpro(npro)='m**2/s'
      CALL RARRAY_COPY(nr_r,d_eff_g_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='D_n_e_r'
      unitpro(npro)='m**2/s'
      CALL RARRAY_COPY(nr_r,d_n_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='v_n_e_r'
      unitpro(npro)='m/s'
      CALL RARRAY_COPY(nr_r,v_nt_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='v_wp_e_r'
      unitpro(npro)='m/s'
      CALL RARRAY_COPY(nr_r,v_eb_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='D_n_g_e_r'
      unitpro(npro)='m**2/s'
      CALL RARRAY_COPY(nr_r,d_n_g_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='v_n_g_e_r'
      unitpro(npro)='m/s'
      CALL RARRAY_COPY(nr_r,v_nt_g_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='v_wp_g_e_r'
      unitpro(npro)='m/s'
      CALL RARRAY_COPY(nr_r,v_eb_g_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='g_Te_e_r'
      unitpro(npro)='-'
      CALL RARRAY_COPY(nr_r,g_te_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='g_Ti_e_r'
      unitpro(npro)='-'
      CALL RARRAY_COPY(nr_r,g_ti_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='g_ne_e_r'
      unitpro(npro)='-'
      CALL RARRAY_COPY(nr_r,g_n_r(1,1,1),1,valpro(1,npro),1)
      j1=1
      DO j=1,mx_ni
        IF(ki0on(j)) THEN
          j1=j1+1
          npro=npro+1
          namepro(npro)='g_n'//ci0(j)(:ici0(j))//'_e_r'
          unitpro(npro)='-'
          CALL RARRAY_COPY(nr_r,g_n_r(1,j1,1),1,valpro(1,npro),1)
        ENDIF
      ENDDO
      IF(kim1on) THEN
        j1=j1+1
        npro=npro+1
        namepro(npro)='g_n'//cim1(:icim1)//'_e_r'
        unitpro(npro)='-'
        CALL RARRAY_COPY(nr_r,g_n_r(1,j1,1),1,valpro(1,npro),1)
      ENDIF
!     Ions
      j1=1
      DO j=1,mx_ni
        IF(ki0on(j)) THEN
          j1=j1+1
          npro=npro+1
          namepro(npro)='Sqz_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='-'
          CALL RARRAY_COPY(nr_r,sqz_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='Gam_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='#/m**2/s'
          CALL RARRAY_COPY(nr_r,gam_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='D_eff_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='m**2/s'
          CALL RARRAY_COPY(nr_r,d_eff_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='D_eff_g_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='m**2/s'
          CALL RARRAY_COPY(nr_r,d_eff_g_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='D_n_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='m**2/s'
          CALL RARRAY_COPY(nr_r,d_n_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='v_n_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='m/s'
          CALL RARRAY_COPY(nr_r,v_nt_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='v_wp_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='m/s'
          CALL RARRAY_COPY(nr_r,v_eb_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='D_n_g_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='m**2/s'
          CALL RARRAY_COPY(nr_r,d_n_g_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='v_n_g_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='m/s'
          CALL RARRAY_COPY(nr_r,v_nt_g_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='v_wp_g_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='m/s'
          CALL RARRAY_COPY(nr_r,v_eb_g_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='g_Te_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='-'
          CALL RARRAY_COPY(nr_r,g_te_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='g_Ti_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='-'
          CALL RARRAY_COPY(nr_r,g_ti_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='g_ne_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='-'
          CALL RARRAY_COPY(nr_r,g_n_r(1,1,j1),1,valpro(1,npro),1)
          k1=1
          DO k=1,mx_ni
            IF(ki0on(k)) THEN
              k1=k1+1
              npro=npro+1
              namepro(npro)='g_n'//ci0(k)(:ici0(k))//'_'
     &                      //ci0(j)(:ici0(j))//'_r'
              unitpro(npro)='-'
              CALL RARRAY_COPY(nr_r,g_n_r(1,k1,j1),1,valpro(1,npro),1)
            ENDIF
          ENDDO
          IF(kim1on) THEN
            npro=npro+1
            namepro(npro)='g_n'//cim1(:icim1)//'_'
     &                    //ci0(j)(:ici0(j))//'_r'
            unitpro(npro)='-'
            CALL RARRAY_COPY(nr_r,g_n_r(1,k1+1,j1),1,valpro(1,npro),1)
          ENDIF
        ENDIF
      ENDDO
!     Single charge state impurity
      IF(kim1on) THEN
        j1=j1+1
        npro=npro+1
        namepro(npro)='Sqz_'//cim1(:icim1)//'_r'
        unitpro(npro)='-'
        CALL RARRAY_COPY(nr_r,sqz_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='Gam_'//cim1(:icim1)//'_r'
        unitpro(npro)='#/m**2/s'
        CALL RARRAY_COPY(nr_r,gam_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='D_eff_'//cim1(:icim1)//'_r'
        unitpro(npro)='m**2/s'
        CALL RARRAY_COPY(nr_r,d_eff_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='D_eff_g_'//cim1(:icim1)//'_r'
        unitpro(npro)='m**2/s'
        CALL RARRAY_COPY(nr_r,d_eff_g_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='D_n_'//cim1(:icim1)//'_r'
        unitpro(npro)='m**2/s'
        CALL RARRAY_COPY(nr_r,d_n_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='v_n_'//cim1(:icim1)//'_r'
        unitpro(npro)='m/s'
        CALL RARRAY_COPY(nr_r,v_nt_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='v_wp_'//cim1(:icim1)//'_r'
        unitpro(npro)='m/s'
        CALL RARRAY_COPY(nr_r,v_eb_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='D_n_g_'//cim1(:icim1)//'_r'
        unitpro(npro)='m**2/s'
        CALL RARRAY_COPY(nr_r,d_n_g_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='v_n_g_'//cim1(:icim1)//'_r'
        unitpro(npro)='m/s'
        CALL RARRAY_COPY(nr_r,v_nt_g_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='v_wp_g_'//cim1(:icim1)//'_r'
        unitpro(npro)='m/s'
        CALL RARRAY_COPY(nr_r,v_eb_g_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='g_Te_'//cim1(:icim1)//'_r'
        unitpro(npro)='-'
        CALL RARRAY_COPY(nr_r,g_te_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='g_Ti_'//cim1(:icim1)//'_r'
        unitpro(npro)='-'
        CALL RARRAY_COPY(nr_r,g_ti_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='g_ne_'//cim1(:icim1)//'_r'
        unitpro(npro)='-'
        CALL RARRAY_COPY(nr_r,g_n_r(1,1,j1),1,valpro(1,npro),1)
        k1=1
        DO k=1,mx_ni
          IF(ki0on(k)) THEN
            k1=k1+1
            npro=npro+1
            namepro(npro)='g_n'//ci0(k)(:ici0(k))//'_'
     &                      //cim1(:icim1)//'_r'
            unitpro(npro)='-'
            CALL RARRAY_COPY(nr_r,g_n_r(1,k1,j1),1,valpro(1,npro),1)
          ENDIF
        ENDDO
        npro=npro+1
        namepro(npro)='g_n'//cim1(:icim1)//'_'//cim1(:icim1)//'_r'
        unitpro(npro)='-'
        CALL RARRAY_COPY(nr_r,g_n_r(1,j1,j1),1,valpro(1,npro),1)
      ENDIF
!     Multiple charge state impurity
      DO j=1,izim2
        IF(kim2on(j)) THEN
          npro=npro+1
          namepro(npro)='Gam_'//cim2(:icim2)//'('//cn(j)//')_r'
          unitpro(npro)='#/m**2/s'
          CALL RARRAY_COPY(nr_r,gam_im2_r(1,j),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='D_n_g_'//cim2(:icim2)//'('//cn(j)//')_r'
          unitpro(npro)='m**2/s'
          CALL RARRAY_COPY(nr_r,d_n_g_im2_r(1,j),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='v_n_g_'//cim2(:icim2)//'('//cn(j)//')_r'
          unitpro(npro)='m/s'
          CALL RARRAY_COPY(nr_r,v_nt_g_im2_r(1,j),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='vwp_g_'//cim2(:icim2)//'('//cn(j)//')_r'
          unitpro(npro)='m/s'
          CALL RARRAY_COPY(nr_r,v_eb_g_im2_r(1,j),1,valpro(1,npro),1)
        ENDIF
      ENDDO
!  Conduction heat fluxes and conductivities
      npro=npro+1
      namepro(npro)='q_con_e_r'
      unitpro(npro)='w/m**2'
      CALL RARRAY_COPY(nr_r,q_con_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='q_con_i_r'
      unitpro(npro)='w/m**2'
      CALL RARRAY_COPY(nr_r,q_con_r(1,2),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='chi_eff_e_r'
      unitpro(npro)='m**2/s'
      CALL RARRAY_COPY(nr_r,chi_eff_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='chi_eff_i_r'
      unitpro(npro)='m**2/s'
      CALL RARRAY_COPY(nr_r,chi_eff_r(1,2),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='chi_eff_g_e_r'
      unitpro(npro)='m**2/s'
      CALL RARRAY_COPY(nr_r,chi_eff_g_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='chi_eff_g_i_r'
      unitpro(npro)='m**2/s'
      CALL RARRAY_COPY(nr_r,chi_eff_g_r(1,2),1,valpro(1,npro),1)
!  Electrical resistivity
      npro=npro+1
      namepro(npro)='eta_par_r'
      unitpro(npro)='Ohm*m'
      CALL RARRAY_COPY(nr_r,eta_par_r,1,valpro(1,npro),1)
!  Parallel currents
      npro=npro+1
      namepro(npro)='J_r'
      unitpro(npro)='A/m**2'
      CALL RARRAY_COPY(nr_r,xj_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='I_r'
      unitpro(npro)='A'
      CALL RARRAY_COPY(nr_r,cur_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='J_bs_r'
      unitpro(npro)='A/m**2'
      CALL RARRAY_COPY(nr_r,xj_bs_r,1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='I_bs_r'
      unitpro(npro)='A'
      CALL RARRAY_COPY(nr_r,cur_bs_r,1,valpro(1,npro),1)
      IF(ABS(xj_ex_r(1)).gt.0.0) THEN
        npro=npro+1
        namepro(npro)='J_ex_r'
        unitpro(npro)='A/m**2'
        CALL RARRAY_COPY(nr_r,xj_ex_r,1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='I_ex_r'
        unitpro(npro)='A'
        CALL RARRAY_COPY(nr_r,cur_ex_r,1,valpro(1,npro),1)
      ENDIF
      IF(ABS(xj_nb_ex_r(1)).gt.0.0) THEN
        npro=npro+1
        namepro(npro)='J_nb_ex_r'
        unitpro(npro)='A/m**2'
        CALL RARRAY_COPY(nr_r,xj_nb_ex_r,1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='I_nb_ex_r'
        unitpro(npro)='A'
        CALL RARRAY_COPY(nr_r,cur_nb_ex_r,1,valpro(1,npro),1)
      ENDIF
!  Parallel electric field
      npro=npro+1
      namepro(npro)='E_par_r'
      unitpro(npro)='V/m'
      CALL RARRAY_COPY(nr_r,e_par_r,1,valpro(1,npro),1)
      IF(ABS(e_par_ex_r(1)).gt.0.0) THEN
        npro=npro+1
        namepro(npro)='E_par_ex_r'
        unitpro(npro)='V/m'
        CALL RARRAY_COPY(nr_r,e_par_ex_r,1,valpro(1,npro),1)
      ENDIF
!  Radial electric field
      IF(kim1on) THEN
        npro=npro+1
        namepro(npro)='E_rad_tot_r'
        unitpro(npro)='V/m'
        CALL RARRAY_COPY(nr_r,e_rad_r(1,1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='E_rad_t_r'
        unitpro(npro)='V/m'
        CALL RARRAY_COPY(nr_r,e_rad_r(1,2),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='E_rad_p_r'
        unitpro(npro)='V/m'
        CALL RARRAY_COPY(nr_r,e_rad_r(1,3),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='E_rad_prs_r'
        unitpro(npro)='V/m'
        CALL RARRAY_COPY(nr_r,e_rad_r(1,4),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='E_rad_ex_r'
        unitpro(npro)='V/m'
        CALL RARRAY_COPY(nr_r,e_rad_ex_r,1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='E_rad_tot_o_r'
        unitpro(npro)='V/m'
        CALL RARRAY_COPY(nr_r,e_rad_o_r(1,1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='E_rad_t_o_r'
        unitpro(npro)='V/m'
        CALL RARRAY_COPY(nr_r,e_rad_o_r(1,2),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='E_rad_p_o_r'
        unitpro(npro)='V/m'
        CALL RARRAY_COPY(nr_r,e_rad_o_r(1,3),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='E_rad_prs_o_r'
        unitpro(npro)='V/m'
        CALL RARRAY_COPY(nr_r,e_rad_o_r(1,4),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='E_rad_ex_o_r'
        unitpro(npro)='V/m'
        CALL RARRAY_COPY(nr_r,e_rad_ex_o_r,1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='omega_exb_o_r'
        unitpro(npro)='/s'
        CALL RARRAY_COPY(nr_r,omexb_o_r,1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='omega_exb_ex_o_r'
        unitpro(npro)='/s'
        CALL RARRAY_COPY(nr_r,omexb_ex_o_r,1,valpro(1,npro),1)
      ENDIF
!  Particle toroidal and poloidal flow velocities on the outside
      npro=npro+1
      namepro(npro)='v_t_o_e_r'
      unitpro(npro)='m/s'
      CALL RARRAY_COPY(nr_r,v_t_o_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='M_t_o_e_r'
      unitpro(npro)='-'
      CALL RARRAY_COPY(nr_r,xm_t_o_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='v_p_o_e_r'
      unitpro(npro)='m/s'
      CALL RARRAY_COPY(nr_r,v_p_o_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='M_p_o_e_r'
      unitpro(npro)='-'
      CALL RARRAY_COPY(nr_r,xm_p_o_r(1,1),1,valpro(1,npro),1)
      j1=1
      DO j=1,mx_ni
        IF(ki0on(j)) THEN
          j1=j1+1
          npro=npro+1
          namepro(npro)='v_t_o_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='m/s'
          CALL RARRAY_COPY(nr_r,v_t_o_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='M_t_o_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='-'
          CALL RARRAY_COPY(nr_r,xm_t_o_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='v_p_o_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='m/s'
          CALL RARRAY_COPY(nr_r,v_p_o_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='M_p_o_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='-'
          CALL RARRAY_COPY(nr_r,xm_p_o_r(1,j1),1,valpro(1,npro),1)
        ENDIF
      ENDDO
      IF(kim1on) THEN
        j1=j1+1
        jimp=j1
        npro=npro+1
        namepro(npro)='v_t_o_'//cim1(:icim1)//'_r'
        unitpro(npro)='m/s'
        CALL RARRAY_COPY(nr_r,v_t_o_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='M_t_o_'//cim1(:icim1)//'_r'
        unitpro(npro)='-'
        CALL RARRAY_COPY(nr_r,xm_t_o_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='M_t_o_'//cim1(:icim1)//'_ex_r'
        unitpro(npro)='-'
        CALL RARRAY_COPY(nr_r,xm_t_im_ex_o_r,1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='v_p_o_'//cim1(:icim1)//'_r'
        unitpro(npro)='m/s'
        CALL RARRAY_COPY(nr_r,v_p_o_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='M_p_o_'//cim1(:icim1)//'_r'
        unitpro(npro)='-'
        CALL RARRAY_COPY(nr_r,xm_p_o_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='M_p_o_'//cim1(:icim1)//'_ex_r'
        unitpro(npro)='-'
        CALL RARRAY_COPY(nr_r,xm_p_im_ex_o_r,1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='v_t_o_'//cim1(:icim1)//'_ex_r'
        unitpro(npro)='m/s'
        CALL RARRAY_COPY(nr_r,vt_im_ex_r,1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='v_p_o_'//cim1(:icim1)//'_ex_r'
        unitpro(npro)='m/s'
        CALL RARRAY_COPY(nr_r,vp_im_ex_r,1,valpro(1,npro),1)
      ENDIF
!  Flux surface particle flow velocities
      npro=npro+1
      namepro(npro)='Omega_e_r'
      unitpro(npro)='rad/s'
      CALL RARRAY_COPY(nr_r,omega_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='u_par_e_r'
      unitpro(npro)='T*m/s'
      CALL RARRAY_COPY(nr_r,u_par_r(1,1),1,valpro(1,npro),1)
      npro=npro+1
      namepro(npro)='u_p_e_r'
      unitpro(npro)='m/s/T'
      CALL RARRAY_COPY(nr_r,u_p_r(1,1),1,valpro(1,npro),1)
      j1=1
      DO j=1,mx_ni
        IF(ki0on(j)) THEN
          j1=j1+1
          npro=npro+1
          namepro(npro)='Omega_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='rad/s'
          CALL RARRAY_COPY(nr_r,omega_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='u_par_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='T*m/s'
          CALL RARRAY_COPY(nr_r,u_par_r(1,j1),1,valpro(1,npro),1)
          npro=npro+1
          namepro(npro)='u_p_'//ci0(j)(:ici0(j))//'_r'
          unitpro(npro)='m/s/T'
          CALL RARRAY_COPY(nr_r,u_p_r(1,j1),1,valpro(1,npro),1)
        ENDIF
      ENDDO
      IF(kim1on) THEN
        j1=j1+1
        npro=npro+1
        namepro(npro)='Omega_'//cim1(:icim1)//'_r'
        unitpro(npro)='rad/s'
        CALL RARRAY_COPY(nr_r,omega_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='Omega_'//cim1(:icim1)//'_ex_r'
        unitpro(npro)='rad/s'
        CALL RARRAY_COPY(nr_r,omega_ex_r,1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='u_par_'//cim1(:icim1)//'_r'
        unitpro(npro)='T*m/s'
        CALL RARRAY_COPY(nr_r,u_par_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='u_p_'//cim1(:icim1)//'_r'
        unitpro(npro)='m/s/T'
        CALL RARRAY_COPY(nr_r,u_p_r(1,j1),1,valpro(1,npro),1)
        npro=npro+1
        namepro(npro)='u_p_'//cim1(:icim1)//'_ex_r'
        unitpro(npro)='m/s/T'
        CALL RARRAY_COPY(nr_r,xk_im_ex_r,1,valpro(1,npro),1)
      ENDIF
!  Pad output (using 5 values per line) to fill lines for IDL
      npropr=((npro+4)/5)*5
      nrpr=((nr_r+4)/5)*5
      DO i=npro+1,npropr
        namepro(i)='dummy'
        unitpro(i)='dummy'
      ENDDO
!Write variable names, units and profiles to 1D data file
      OPEN(unit=n1d,file=cn1d,status='unknown',form='formatted')
      idum(1)=npro
      idum(2)=nr_r
      rdum(1)=0.0
      CALL WRITE_IR(n1d,2,idum,0,rdum,15,0)
!  1-D variable names
      DO j=1,npropr/5
        CALL WRITE_C(n1d,5,namepro(5*j-4),15)
      ENDDO
!  1-D variable units
      DO j=1,npropr/5
        CALL WRITE_C(n1d,5,unitpro(5*j-4),15)
      ENDDO
      DO j=1,npropr
        DO i=1,nrpr/5
          CALL WRITE_IR(n1d,0,idum,5,valpro(1+5*(i-1),j),15,2)
        ENDDO
      ENDDO
      CLOSE(unit=n1d)
!Summary of calculation
!  0-D
      OPEN(unit=nsum,file=cnsum,status='unknown',form='formatted')
      WRITE(nsum,'(1x,a,1(1x,a9))')
     &'Device                                 =',device
      WRITE(nsum,'(1x,a,1(1x,i9))')
     &'Shot number                            =',idshot
      WRITE(nsum,'(1x,a,1(1x,f9.5))')
     &'Time                               (m) =',time
      IF(device.eq.'cmod') THEN
        WRITE(nsum,'(1x,a,1(1x,a9))')
     &  'TRANSP run id for profiles             =',idtransp
      ENDIF
      IF(device.eq.'jet') THEN
        WRITE(nsum,'(1x,a,2(1x,a9),1(1x,i9))')
     &  'UID, DDA and ISEQP - Te profile        =',uid(1),dda(1),
     &                                             iseqp(1)
        WRITE(nsum,'(1x,a,2(1x,a9),1(1x,i9))')
     &  'UID, DDA and ISEQP - Ti profile        =',uid(2),dda(2),
     &                                             iseqp(2)
        WRITE(nsum,'(1x,a,2(1x,a9),1(1x,i9))')
     &  'UID, DDA and ISEQP - ne profile        =',uid(3),dda(3),
     &                                             iseqp(3)
        WRITE(nsum,'(1x,a,2(1x,a9),1(1x,i9))')
     &  'UID, DDA and ISEQ - imp 1 den profile  =',uid(4),dda(4),
     &                                             iseqp(4)
        WRITE(nsum,'(1x,a,2(1x,a9),1(1x,i9))')
     &  'UID, DDA and ISEQP - imp 1 tor rot     =',uid(5),dda(5),
     &                                             iseqp(5)
        WRITE(nsum,'(1x,a,2(1x,a9),1(1x,i9))')
     &  'UID, DDA and ISEQP - imp 2 den profile =',uid(6),dda(6),
     &                                             iseqp(6)
      ENDIF
      WRITE(nsum,'(1x,a,1(1x,i9))')
     &'Number of radial points (nr_r)         =',nr_r 
      WRITE(nsum,'(1x,a,1(1x,f9.5))')
     &'Reference minor radius (a0)        (m) =',a0  
      WRITE(nsum,'(1x,a,1(1x,f9.5))')
     &'Reference major radius (r0)        (m) =',r0  
      WRITE(nsum,'(1x,a,1(1x,f9.5))')
     &'Toroidal field at r0               (T) =',bt0
      WRITE(nsum,'(1x,a,1(1x,f9.5))')
     &'Elongation at edge                     =',elong_r(nr_r)
!  1-D radial grid and variables     
      DO j=1,npropr/5
        cdum(1)='          i '
        DO i=1,5
          cdum(i+1)=namepro(5*(j-1)+i)
        ENDDO
        CALL WRITE_C(nsum,6,cdum,15)
        cdum(1)='            '
        DO i=1,5
          cdum(i+1)=unitpro(5*(j-1)+i)
        ENDDO
        CALL WRITE_C(nsum,6,cdum,15)
        DO i=1,nr_r
          idum(1)=i
          rdum(1)=valpro(i,5*j-4)
          rdum(2)=valpro(i,5*j-3)
          rdum(3)=valpro(i,5*j-2)
          rdum(4)=valpro(i,5*j-1)
          rdum(5)=valpro(i,5*j)
          CALL WRITE_IR(nsum,1,idum,5,rdum,15,2)
        ENDDO
      ENDDO
 1000 IF(iflag.ne.0) THEN
        CALL WRITE_LINE(nout,message,0,0)
      ENDIF
      CLOSE(unit=nsum)
      if(xneo.eq.0 .and. i_proc.eq.0) then
         WRITE(*,*) 'FORCEBAL Run Completed'
      endif
      RETURN
      END
      SUBROUTINE NCLASS_DR(k_edotb,k_potato,k_sqz,amui0,izi0,amuim1,
     &                     izim1,amuim2,izim2,c_den,a0,bt0,q0,el0,nr_r,
     &                     b2_r,bm2_r,fhat_r,ftrap_r,fm_r,grth_r,
     &                     grho2_r,gr2bm2_r,rhot_r,den_r,deni_r,
     &                     denim_r,denim2_r,te_r,ti_r,rout_r,btout_r,
     &                     bpout_r,vp_im_ex_r,vt_im_ex_r,psi_r,psip_r,
     &                     e_par_ex_r,xj_r,xj_nb_ex_r,iflag,message)
!***********************************************************************
!NCLASS_DR generates profiles of NCLASS quantities using the MHD
!  equilibrium along with plasma density, temperature and toroidal
!  rotation profiles
!  W.A.Houlberg 3/2000
!Input:
!  k_edotb-option to use experimental or calculated <E.B>
!        =1 use experimental
!        =else use calculated <E.B>=eta(<J.B>_ex-<J.B>_bs-<J.B>_nb_ex)
!         where eta and <J.B>_bs are calculated
!  k_potato-option to include potato orbit effects (=0 is off)
!  k_sqz-option to include orbit squeezing (=0 is off)
!  amui0(i)-atomic mass of fully stripped ions
!  izi0(i)-charge of fully stripped ions
!  amuim1-atomic mass of single charge state diagnostic impurity
!  izim1-charge of single charge state diagnostic impurity
!  amuim2-atomic mass of multiple charge state impurity
!  izim2-max charge of multiple charge state impurity
!  c_den-density cutoff below which species is ignored (/m**3)
!  a0-plasma minor radius (m)
!  bt0-toroidal field at r0 (T)
!  q0-safety factor at magnetic axis
!  el0-elongation at magnetic axis
!  nr_r-number of radial points
!  b2_r(i)-<B**2> (T**2)
!  bm2_r(i)-<1/B**2> (1/T**2)
!  fhat_r(i)-RB_t/dpsi_r/drho_t/a0
!  ftrap_r(i)-trapped particle fraction
!  fm_r(3,i)-geometric factor
!  grth_r(i)-<n.grad(theta)> (1/m)
!  grho2_r(i)-a0**2*<|grad(rhot_r)|**2>
!  gr2bm2_r(i)-a0**2*<|grad(rhot_r)|**2/B**2> (1/T**2)
!  rhot_r(i)-normalized toroidal flux grid proportional to (Phi)**0.5
!  den_r(i)-electron density (/m**3)
!  deni_r(i,j)-fully stripped ion density (/m**3)
!  denim_r(i)-single charge state impurity density (/m**3)
!  denim2_r(i,k)-multiple charge state impurity density (/m**3)
!  te_r(i)-electron temperature (keV)
!  ti_r(i)-ion temperature (keV)
!  rout_r(i)-major radius on outside of flux surface (m)
!  btout_r(i)-toroidal field at rout_r(i) (T)
!  bpout_r(i)-poloidal field at rout_r(i) (T)
!  vp_im_ex_r(i)-expt pol velocity of impurity on outside (m/s)
!  vt_im_ex_r(i)-expt tor velocity of impurity on outside (m/s)
!  psi_r(i)-poloidal flux/2*pi (Wb/rad)
!  psip_r(i)-poloidal flux gradient (Wb/m)
!  e_par_ex_r(i)-expt <E.B>/bt0 (V/m)
!  xj_r(i)-expt <J.B>/bt0 (A/m**2)
!  xj_nb_ex_r(i)-expt <J.B>_nb/bt0 (A/m**2)
!Output:
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
!Declaration of parameters and commons
      include 'mpif.h'
      INCLUDE '../inc/pamx_mi.inc'
      INCLUDE '../inc/pamx_ms.inc'
      INCLUDE '../inc/pamx_mz.inc'
      include '../inc/glf.m'
c
      INTEGER        mx_ni
      PARAMETER     (mx_ni=5)
      INTEGER        mxnr_r
      PARAMETER     (mxnr_r=300)
      INCLUDE '../inc/comfbl.inc'
!Declaration of input variables
      CHARACTER*(*)  message
      INTEGER        izi0(*),                 izim1,
     &               izim2,                   k_edotb,
     &               k_sqz,                   nr_r                    
      REAL           amui0(*),                amuim1,
     &               amuim2
      REAL           a0,                      bt0,
     &               c_den,                   el0,
     &               q0
      REAL           b2_r(*),                 bm2_r(*),
     &               bpout_r(*),              btout_r(*),
     &               fhat_r(*),               ftrap_r(*),
     &               fm_r(3,*),               grth_r(*),
     &               grho2_r(*),              gr2bm2_r(*),
     &               psi_r(*),                psip_r(*),
     &               rhot_r(*),               rout_r(*)
      REAL           den_r(*),                deni_r(mxnr_r,*),               
     &               denim_r(*),              e_par_ex_r(*),
     &               te_r(*),                 ti_r(*),
     &               vp_im_ex_r(*),           vt_im_ex_r(*),
     &               xj_r(*),                 xj_nb_ex_r(*)
      REAL           denim2_r(mxnr_r,*)
!Declaration of input to NCLASS
      INTEGER        k_order,                 k_potato
      INTEGER        m_i,                     m_z
      REAL           c_potb,                  c_potl
      REAL           p_b2,                    p_bm2,
     &               p_eb,                    p_fhat,
     &               p_fm(3),                 p_ft,
     &               p_grbm2,                 p_grphi,
     &               p_gr2phi,                p_ngrth
      REAL           amu_i(mx_mi),            grt_i(mx_mi),
     &               temp_i(mx_mi)
      REAL           den_iz(mx_mi,mx_mz),     fex_iz(3,mx_mi,mx_mz),
     &               grp_iz(mx_mi,mx_mz)
!Declaration of output from NCLASS
      INTEGER        iflag,                   m_s
      INTEGER        jm_s(mx_ms),             jz_s(mx_ms)
      REAL           p_bsjb,                  p_etap,
     &               p_exjb
      REAL           calm_i(3,3,mx_mi)
      REAL           caln_ii(3,3,mx_mi,mx_mi),capm_ii(3,3,mx_mi,mx_mi),
     &               capn_ii(3,3,mx_mi,mx_mi)
      REAL           bsjbp_s(mx_ms),          bsjbt_s(mx_ms),
     &               chi_s(mx_ms),            dn_s(mx_ms),
     &               gfl_s(5,mx_ms),          qfl_s(5,mx_ms),
     &               sqz_s(mx_ms),
     &               upar_s(3,3,mx_ms),       utheta_s(3,3,mx_ms),
     &               vneb_s(mx_ms),           vnex_s(mx_ms),
     &               vnnt_s(mx_ms),
     &               vqeb_s(mx_ms),           vqex_s(mx_ms),
     &               vqnt_s(mx_ms),
     &               xi_s(mx_ms),             ymu_s(3,3,mx_ms)
      REAL           chip_ss(mx_ms,mx_ms),    chit_ss(mx_ms,mx_ms),
     &               dp_ss(mx_ms,mx_ms),      dt_ss(mx_ms,mx_ms)
!Declaration of local variables
      LOGICAL        ki0on(mx_ni),            kim1on,
     &               kim2on
      INTEGER        i,                       im,
     &               ima,                     iz,
     &               iza,                     j,
     &               jiter,                   jitermx,
     &               k,                       k2,
     &               kdiag,                   kk,l
      REAL           dr,                      drout,
     &               dent,                    erref,
     &               grad,                    omer,
     &               ppr,                     toler,
     &               vtherm,                  xfsqz
      REAL           z_coulomb,               z_electronmass,
     &               z_j7kv,                  z_protonmass 
      REAL           erold1_r(mxnr_r),        erold2_r(mxnr_r),
     &               ertemp_r(mxnr_r)
      REAL           omegas_r(mxnr_r,mx_mi)
      REAL           temp_ri(mxnr_r,mx_mi),   grt_ri(mxnr_r,mx_mi)
      REAL           den_riz(mxnr_r,mx_mi,mx_mz),
     &               grp_riz(mxnr_r,mx_mi,mx_mz)
!Declaration of external functions
      INTEGER        IARRAY_MAX
      REAL           RARRAY_SUM
!MPI variables
      REAL  chi_eff_sum(mxnr_r,2), d_eff_sum(mxnr_r,mx_ms),
     &      d_effr_sum(mxnr_r,mx_ms),
     &      v_p_o_sum(mxnr_r,mx_ms), v_t_o_sum(mxnr_r,mx_ms),
     &      e_rad_mp(mxnr_r), e_rad_mpsum(mxnr_r),
     &      e_rad_sum(mxnr_r,4)
!Set values for call to NCLASS
!  Options and constants
      jitermx=50
!     toler=0.01
      toler=0.02
      xfsqz=0.85
      k_order=2
      c_potb=el0*bt0/2.0/q0**2
      c_potl=rout_r(1)*q0
!  Physical and conversion constants
      z_coulomb=1.6022e-19
      z_electronmass=9.1095e-31
      z_j7kv=1.6022e-16
      z_protonmass=1.6726e-27
!  Zero out external force
      CALL RARRAY_ZERO(3*mx_mi*mx_mz,fex_iz)
!  Zero out temporary arrays
      CALL RARRAY_ZERO(mxnr_r*mx_mi,temp_ri)
      CALL RARRAY_ZERO(mxnr_r*mx_mi,grt_ri)
      CALL RARRAY_ZERO(mxnr_r*mx_mi*mx_mz,den_riz)
      CALL RARRAY_ZERO(mxnr_r*mx_mi*mx_mz,grp_riz)
!Check species and fill in ion densities if below cutoff
      m_i=1
      m_z=1
      amu_i(1)=z_electronmass/z_protonmass
      DO i=1,nr_r
        temp_ri(i,1)=te_r(i)
        den_riz(i,1,1)=den_r(i)
        IF(i.eq.1) THEN
          grt_ri(i,1)=0.0
          grp_riz(i,1,1)=0.0
        ELSEIF(i.eq.nr_r) THEN
          grt_ri(i,1)=grt_ri(nr_r-1,1)
          grp_riz(i,1,1)=grp_riz(nr_r-1,1,1)
        ELSE
          dr=(rhot_r(i+1)-rhot_r(i-1))*a0
          grt_ri(i,1)=(te_r(i+1)-te_r(i-1))/dr
          grp_riz(i,1,1)=(den_r(i+1)*te_r(i+1)
     &                   -den_r(i-1)*te_r(i-1))/dr
        ENDIF
c       write(*,88) i, rhot_r(i),deni_r(i,1), deni_r(i,2),
c    &             deni_r(i,3),deni_r(i,4)
c    &     den_r(i)*(te_r(i+1)-te_r(i-1))/dr+te_r(i)*
c    &     (den_r(i+1)-den_r(i-1))/dr
      ENDDO
!  Single charge state ions
      DO k=1,mx_ni
        ki0on(k)=.false.
        k2=IARRAY_MAX(nr_r,deni_r(1,k),1)
        IF(deni_r(k2,k).gt.c_den) THEN
          ki0on(k)=.true.
          m_i=m_i+1
          amu_i(m_i)=amui0(k)
          IF(izi0(k).gt.m_z) m_z=izi0(k)
          DO i=1,nr_r
            temp_ri(i,m_i)=ti_r(i)
            den_riz(i,m_i,izi0(k))=deni_r(i,k)
            IF(den_riz(i,m_i,izi0(k)).lt.2.0*c_den)
     &        den_riz(i,m_i,izi0(k))=2.0*c_den
            IF(i.eq.1) THEN
              grt_ri(i,m_i)=0.0
              grp_riz(i,m_i,izi0(k))=0.0
            ELSEIF(i.eq.nr_r) THEN
              grt_ri(i,m_i)=grt_ri(nr_r-1,m_i)
              grp_riz(i,m_i,izi0(k))=grp_riz(nr_r-1,m_i,izi0(k))
            ELSE
              dr=(rhot_r(i+1)-rhot_r(i-1))*a0
              grt_ri(i,m_i)=(ti_r(i+1)-ti_r(i-1))/dr
              grp_riz(i,m_i,izi0(k))=(deni_r(i+1,k)*ti_r(i+1)
     &                               -deni_r(i-1,k)*ti_r(i-1))/dr
c           write(*,89) i,rhot_r(i),k,m_i,izi0(k),
c    &                  grp_riz(i,m_i,izi0(k)),'forcebal'
            ENDIF
          ENDDO
        ENDIF
      ENDDO
!  Single charge state diagnostic impurity
      kim1on=.false.
      k2=IARRAY_MAX(nr_r,denim_r,1)
      kdiag=0
      IF(denim_r(k2).gt.c_den) THEN
        kim1on=.true.
        m_i=m_i+1
        kdiag=m_i
        amu_i(m_i)=amuim1
        IF(izim1.gt.m_z) m_z=izim1
        DO i=1,nr_r
          temp_ri(i,m_i)=ti_r(i)
          den_riz(i,m_i,izim1)=denim_r(i)
          IF(den_riz(i,m_i,izim1).lt.2.0*c_den)
     &      den_riz(i,m_i,izim1)=2.0*c_den
          IF(i.eq.1) THEN
            grt_ri(i,m_i)=0.0
            grp_riz(i,m_i,izim1)=0.0
          ELSEIF(i.eq.nr_r) THEN
            grt_ri(i,m_i)=grt_ri(nr_r-1,m_i)
            grp_riz(i,m_i,izim1)=grp_riz(nr_r-1,m_i,izim1)
          ELSE
            dr=(rhot_r(i+1)-rhot_r(i-1))*a0
            grt_ri(i,m_i)=(ti_r(i+1)-ti_r(i-1))/dr
            grp_riz(i,m_i,izim1)=(denim_r(i+1)*ti_r(i+1)
     &                           -denim_r(i-1)*ti_r(i-1))/dr
          ENDIF
        ENDDO
      ENDIF
!  Multiple charge state impurity
      kim2on=.false.
      IF(izim2.gt.0) THEN
        DO k=1,izim2
          k2=IARRAY_MAX(nr_r,denim2_r(1,k),1)
          IF(denim2_r(k2,k).gt.c_den) THEN
            kim2on=.true.
            IF(k.gt.m_z) m_z=k
            DO i=1,nr_r
              temp_ri(i,m_i+1)=ti_r(i)
              den_riz(i,m_i+1,k)=denim2_r(i,k)
              IF(den_riz(i,m_i+1,k).lt.2.0*c_den)
     &          den_riz(i,m_i+1,k)=2.0*c_den
              IF(i.eq.1) THEN
                grt_ri(i,m_i+1)=0.0
                grp_riz(i,m_i+1,k)=0.0
              ELSEIF(i.eq.nr_r) THEN
                grt_ri(i,m_i+1)=grt_ri(nr_r-1,m_i+1)
                grp_riz(i,m_i+1,k)=grp_riz(nr_r-1,m_i+1,k)
              ELSE
                dr=(rhot_r(i+1)-rhot_r(i-1))*a0
                grt_ri(i,m_i+1)=(ti_r(i+1)-ti_r(i-1))/dr
                grp_riz(i,m_i+1,k)=(deni_r(i+1,k)*ti_r(i+1)
     &                             -deni_r(i-1,k)*ti_r(i-1))/dr
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!
      IF(kim2on) THEN
        m_i=m_i+1
        amu_i(m_i)=amuim2
      ENDIF
!Parameter for orbit squeezing when k_sqz=0
      p_gr2phi=0.0
!  Initialize Er
      DO i=2,nr_r-1
!       Use dr for flux surface averaged balances
        dr=(rhot_r(i+1)-rhot_r(i-1))*a0
!       Use drout for outside balances
        drout=rout_r(i+1)-rout_r(i-1)
!       Toroidal rotation contribution
        e_rad_r(i,2)=btout_r(i)/fhat_r(i)*vt_im_ex_r(i)
        e_rad_o_r(i,2)=e_rad_r(i,2)*(dr/drout)
!       Poloidal rotation contribution
        e_rad_r(i,3)=0.0
        e_rad_o_r(i,3)=0.0
!       Pressure gradient contribution
        IF(denim_r(i).gt.c_den) THEN
          e_rad_r(i,4)=(denim_r(i+1)*ti_r(i+1)-denim_r(i-1)*ti_r(i-1))
     &                 *z_j7kv/dr/(izim1*z_coulomb*denim_r(i))
          e_rad_o_r(i,4)=e_rad_r(i,4)*(dr/drout)
        ELSE
          e_rad_r(i,4)=0.0
          e_rad_o_r(i,4)=0.0
        ENDIF
!       Totals
        e_rad_r(i,1)=e_rad_r(i,2)+e_rad_r(i,3)+e_rad_r(i,4)
        e_rad_o_r(i,1)=e_rad_o_r(i,2)+e_rad_o_r(i,3)+e_rad_o_r(i,4)
        e_rad_ex_r(i)=e_rad_r(i,2)+e_rad_r(i,4)-vp_im_ex_r(i)
     &                *btout_r(i)**2/bpout_r(i)/fhat_r(i)
        e_rad_ex_o_r(i)=e_rad_o_r(i,2)+e_rad_o_r(i,4)-vp_im_ex_r(i)
     &                  *btout_r(i)**2/bpout_r(i)/fhat_r(i)*(dr/drout)
        erold1_r(i)=e_rad_r(i,1)
        erold2_r(i)=e_rad_r(i,1)
        ertemp_r(i)=e_rad_r(i,1)
      ENDDO
!     Edge values
      e_rad_ex_r(nr_r)=e_rad_ex_r(nr_r-1)
      e_rad_ex_o_r(nr_r)=e_rad_ex_o_r(nr_r-1)
      e_rad_r(nr_r,2)=e_rad_r(nr_r-1,2)
      e_rad_r(nr_r,3)=e_rad_r(nr_r-1,3)
      e_rad_r(nr_r,4)=e_rad_r(nr_r-1,4)
      e_rad_r(nr_r,1)=e_rad_r(nr_r,2)+e_rad_r(nr_r,3)+e_rad_r(nr_r,4)
      erold1_r(nr_r)=e_rad_r(nr_r,1)
      erold2_r(nr_r)=e_rad_r(nr_r,1)
      ertemp_r(nr_r)=e_rad_r(nr_r,1)
!Call NCLASS to get transport properties with iteration on Er
      DO j=1,jitermx
        CALL RARRAY_ZERO(mxnr_r,e_par_r)
        CALL RARRAY_ZERO(mxnr_r,eta_par_r)
        CALL RARRAY_ZERO(mxnr_r,omexb_o_r)
        CALL RARRAY_ZERO(mxnr_r,omexb_ex_o_r)
        CALL RARRAY_ZERO(mxnr_r,xm_p_im_ex_o_r)
        CALL RARRAY_ZERO(mxnr_r,xm_t_im_ex_o_r)
        CALL RARRAY_ZERO(mxnr_r,xj_bs_r)
        CALL RARRAY_ZERO(2*mxnr_r,chi_eff_r)
        CALL RARRAY_ZERO(2*mxnr_r,chi_eff_g_r)
        CALL RARRAY_ZERO(2*mxnr_r,q_con_r)
        CALL RARRAY_ZERO(mxnr_r*mx_mi,omegas_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,gam_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,d_eff_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,d_n_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,v_nt_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,v_eb_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,v_ex_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,d_eff_g_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,d_n_g_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms**2,g_n_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,g_te_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,g_ti_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,v_nt_g_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,v_eb_g_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,v_ex_g_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,u_par_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,u_p_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,v_t_o_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,xm_t_o_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,xm_p_o_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,v_p_o_r)
        CALL RARRAY_ZERO(mxnr_r*mx_ms,omega_r)
        CALL RARRAY_ZERO(mxnr_r*mx_mz,gam_im2_r)
        CALL RARRAY_ZERO(mxnr_r*mx_mz,d_n_g_im2_r)
        CALL RARRAY_ZERO(mxnr_r*mx_mz,v_nt_g_im2_r)
        CALL RARRAY_ZERO(mxnr_r*mx_mz,v_eb_g_im2_r)
        CALL RARRAY_ZERO(mxnr_r*mx_mz,v_ex_g_im2_r)
!       write(*,*) 'jiter = ',j
!
! ** begin loop over grid **
!
        DO i=2+i_proc,nr_r-1,n_proc
          CALL RARRAY_ZERO(mx_mi,temp_i)
          CALL RARRAY_ZERO(mx_mi,grt_i)
          CALL RARRAY_ZERO(mx_mi*mx_mz,den_iz)
          CALL RARRAY_ZERO(mx_mi*mx_mz,grp_iz)
!  Radial step for gradients
          dr=(rhot_r(i+1)-rhot_r(i-1))*a0
          drout=rout_r(i+1)-rout_r(i-1)
!  Radial potential gradient
          p_grphi=-e_rad_r(i,1)
!  Squeezing factor
          IF(k_sqz.ne.0) THEN
            p_gr2phi=(ertemp_r(i)*(psip_r(i+1)-psip_r(i-1))/psip_r(i)
     &               -(ertemp_r(i+1)-ertemp_r(i-1)))/dr
          ENDIF
!  Parallel electric field
          p_eb=1.0
!  Geometric quantities
          p_b2=b2_r(i)
          p_bm2=bm2_r(i)
          p_fhat=fhat_r(i)
          p_ft=ftrap_r(i)
          p_grbm2=gr2bm2_r(i)
          p_ngrth=grth_r(i)
          DO k=1,3
            p_fm(k)=fm_r(k,i)
          ENDDO
!  Temperatures, pressures and gradients
          DO im=1,m_i
            temp_i(im)=temp_ri(i,im)
            grt_i(im)=grt_ri(i,im)
            DO iz=1,m_z
              IF(den_riz(i,im,iz).gt.c_den) THEN
                grp_iz(im,iz)=grp_riz(i,im,iz)
                den_iz(im,iz)=den_riz(i,im,iz)
!               Only for isotopes with single charge state
                IF((im.lt.m_i).or.((im.eq.m_i).and.(.not.kim2on)))
     &          omegas_r(i,im)=-grp_iz(im,iz)*z_j7kv*fhat_r(i)
     &                         /z_coulomb/iz/den_riz(i,im,iz)
     &                         /rout_r(i)/btout_r(i)
                IF(im.eq.1) omegas_r(i,im)=-omegas_r(i,im)
              ENDIF
            ENDDO
          ENDDO
c      write(*,88) i, rhot_r(i), den_iz(1,1), den_iz(2,1)
c      write(*,88) i, rhot_r(i),gfl_s(1,1),gfl_s(2,1),
c    &             gfl_s(3,1),gfl_s(4,1),'forcebal'
 88    format(2x,i2,2x,0p1f10.6,1p4e13.5,a10)
 89    format(2x,i2,2x,0p1f10.6,3i4,2x,1pe13.5,a10)
 91    format(2x,i2,2x,0p1f10.6,2x,1pe13.5,a10)
 92    format(2x,i2,2x,0p1f10.6,2x,1p2e13.5,a10)
!Call NCLASS
          CALL NCLASS(k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2,
     &                p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,
     &                p_gr2phi,p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,
     &                grp_iz,m_s,jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,
     &                caln_ii,capm_ii,capn_ii,bsjbp_s,bsjbt_s,chi_s,
     &                dn_s,gfl_s,qfl_s,sqz_s,upar_s,utheta_s,vnnt_s,
     &                vneb_s,vnex_s,vqnt_s,vqeb_s,vqex_s,xi_s,ymu_s,
     &                chip_ss,chit_ss,dp_ss,dt_ss,iflag,message)
          IF(iflag.eq.1) GOTO 1000
!  Parallel electrical resistivity
          eta_par_r(i)=p_etap
!  Bootstrap current density
          xj_bs_r(i)=p_bsjb/bt0
c      write(*,91) i, rhot_r(i), eta_par_r(i),'forcebal'
c
!  Parallel electric field from calculated resistivity and bootstrap
          e_par_r(i)=p_etap*(xj_r(i)-xj_nb_ex_r(i)-p_bsjb/bt0)
          IF(k_edotb.eq.1) THEN
            edotb_r(i)=e_par_ex_r(i)*bt0
          ELSE
            edotb_r(i)=e_par_r(i)*bt0
          ENDIF
!  Reset <E.B> dependent results
!    Neoclassical pinch
          CALL RARRAY_SCALE(m_s,edotb_r(i),vneb_s,1)
          DO l=1,m_s
!    Radial fluxes
            gfl_s(4,l)=edotb_r(i)*gfl_s(4,l)
            qfl_s(4,l)=edotb_r(i)*qfl_s(4,l)
!    Flow velocities
            DO k=1,k_order
              upar_s(k,2,l)=edotb_r(i)*upar_s(k,2,l)
              utheta_s(k,2,l)=edotb_r(i)*utheta_s(k,2,l)
            ENDDO
          ENDDO
!Save profiles
!     Poloidal flow contribution to radial electric field
          IF(kdiag.ne.0) THEN
            e_rad_r(i,3)=-btout_r(i)**2/fhat_r(i)
     &                   *(utheta_s(1,1,kdiag)+utheta_s(1,2,kdiag)
     &                     +utheta_s(1,3,kdiag))
            e_rad_o_r(i,3)=e_rad_r(i,3)*(dr/drout)
          ENDIF
!         Update total radial electric field
          e_rad_r(i,1)=e_rad_r(i,2)+e_rad_r(i,3)+e_rad_r(i,4)
          e_rad_o_r(i,1)=e_rad_o_r(i,2)+e_rad_o_r(i,3)+e_rad_o_r(i,4)
          e_rad_mp(i)=e_rad_r(i,2)+e_rad_r(i,3)+e_rad_r(i,4)
!     Map species
          IF(kim2on) THEN
            k2=m_i-1
          ELSE
            k2=m_i
          ENDIF
          DO k=1,k2
            im=jm_s(k)
            iz=jz_s(k)
            iza=IABS(iz)
            gam_r(i,k)=RARRAY_SUM(5,gfl_s(1,k),1)
            grad=(den_riz(i+1,im,iza)-den_riz(i-1,im,iza))/dr
            IF(ABS(grad*a0/den_riz(i,im,iza)).gt.1.0e-5)
     &             d_eff_r(i,k)=-gam_r(i,k)/grad
            d_n_r(i,k)=dn_s(k)
            v_nt_r(i,k)=vnnt_s(k) 
            v_eb_r(i,k)=vneb_s(k)
            v_ex_r(i,k)=vnex_s(k)
            d_eff_g_r(i,k)=d_eff_r(i,k)/grho2_r(i)
            d_n_g_r(i,k)=d_n_r(i,k)/grho2_r(i)
            v_nt_g_r(i,k)=v_nt_r(i,k)/SQRT(grho2_r(i))
            v_eb_g_r(i,k)=v_eb_r(i,k)/SQRT(grho2_r(i))
            v_ex_g_r(i,k)=v_ex_r(i,k)/SQRT(grho2_r(i))
!     Orbit squeezing factor
            sqz_r(i,k)=sqz_s(k)
!     Particle flow velocities
!           Average parallel flow
            u_par_r(i,k)=upar_s(1,1,k)+upar_s(1,2,k)+upar_s(1,3,k)
            ppr=grp_iz(im,iza)*z_j7kv/(z_coulomb*iz*den_iz(im,iza))
!           Average poloidal flow
            u_p_r(i,k)=utheta_s(1,1,k)+utheta_s(1,2,k)+utheta_s(1,3,k)
!           Toroidal flow velocity on outside
            v_t_o_r(i,k)=u_p_r(i,k)*btout_r(i)-(ppr-e_rad_r(i,1))
     &                   *fhat_r(i)/btout_r(i)
!           Toroidal Mach number on outside
            vtherm=SQRT(2.0*temp_i(im)*z_j7kv/amu_i(im)/z_protonmass)
            xm_t_o_r(i,k)=v_t_o_r(i,k)/vtherm
!           Poloidal flow velocity on outside
            v_p_o_r(i,k)=u_p_r(i,k)*bpout_r(i)
!           Poloidal Mach number on outside
            xm_p_o_r(i,k)=xm_t_o_r(i,k)-e_rad_o_r(i,1)/bpout_r(i)/vtherm
!     Exponent factors
            g_te_r(i,k)=-(dp_ss(1,k)+dt_ss(1,k))/dp_ss(k,k)
            DO kk=1,k2
              g_ti_r(i,k)=g_ti_r(i,k)
     &                    -(dp_ss(kk,k)+dt_ss(kk,k))/dp_ss(k,k)
              g_n_r(i,kk,k)=-dp_ss(kk,k)/dp_ss(k,k)
            ENDDO
          ENDDO
          IF(kdiag.ne.0) THEN
            vtherm=SQRT(2.0*temp_i(kdiag)*z_j7kv/amu_i(kdiag)
     &                  /z_protonmass)
            xm_t_im_ex_o_r(i)=vt_im_ex_r(i)/vtherm 
            xm_p_im_ex_o_r(i)=xm_t_im_ex_o_r(i)
     &                        -e_rad_ex_r(i)/bpout_r(i)/vtherm
          ENDIF
! Multiple charge state impurity
          IF(kim2on) THEN
!     Orbit squeezing factor-use same as previous isotope
            sqz_r(i,k2+1)=sqz_r(i,k2)
            DO k=1,izim2
              IF(denim2_r(i,k).gt.c_den) THEN
                k2=k2+1
!     Particle fluxes and diffusivities
                gam_im2_r(i,k)=RARRAY_SUM(5,gfl_s(1,k2),1)
                d_n_g_im2_r(i,k)=dn_s(k2)/grho2_r(i)
                v_nt_g_im2_r(i,k)=vnnt_s(k2)/SQRT(grho2_r(i))
                v_eb_g_im2_r(i,k)=vneb_s(k2)/SQRT(grho2_r(i))
                v_ex_g_im2_r(i,k)=vnex_s(k2)/SQRT(grho2_r(i))
              ENDIF
            ENDDO
          ENDIF
!  Conduction fluxes and conductivities
!     Electrons
          q_con_r(i,1)=RARRAY_SUM(5,qfl_s(1,1),1)
c      write(*,91) i, rhot_r(i), qfl_s(i,1),'forcebal'
c    &    qfl_s(2,1),qfl_s(3,1),qfl_s(4,1),qfl_s(5,1)
          IF(ABS(grt_i(1)*a0/temp_i(1)).gt.1.0e-5)
     &      chi_eff_r(i,1)=-q_con_r(i,1)/den_r(i)/z_j7kv/grt_i(1)
          chi_eff_g_r(i,1)=chi_eff_r(i,1)/grho2_r(i)
!     Ions - sum over all species
          dent=0.0
          DO k=2,m_s
            ima=jm_s(k)
            iza=jz_s(k)
            q_con_r(i,2)=q_con_r(i,2)+RARRAY_SUM(5,qfl_s(1,k),1)
     	      dent=dent+den_iz(ima,iza)
          ENDDO
          IF(ABS(grt_i(2)*a0/temp_i(2)).gt.1.0e-5)
     &      chi_eff_r(i,2)=-q_con_r(i,2)/dent/z_j7kv/grt_i(2)
          chi_eff_g_r(i,2)=chi_eff_r(i,2)/grho2_r(i)
        ENDDO
!
!  MPI Reductions
!
        do i=1,2
          call MPI_REDUCE(chi_eff_g_r(:,i),chi_eff_sum(:,i),mxnr_r,
     &         MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,i_err)
          call MPI_BCAST(chi_eff_sum(:,i),mxnr_r,
     &         MPI_REAL,0,MPI_COMM_WORLD,i_err)
        enddo
          call MPI_REDUCE(e_rad_r(:,3),e_rad_sum(:,3),mxnr_r,
     &         MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,i_err)
          call MPI_BCAST(e_rad_sum(:,3),mxnr_r,
     &         MPI_REAL,0,MPI_COMM_WORLD,i_err)
          call MPI_REDUCE(e_rad_mp,e_rad_mpsum,mxnr_r,
     &         MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,i_err)
          call MPI_BCAST(e_rad_mpsum,mxnr_r,
     &         MPI_REAL,0,MPI_COMM_WORLD,i_err)
        do i=1,k2
          call MPI_REDUCE(d_eff_r(:,i),d_effr_sum(:,i),mxnr_r,
     &         MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,i_err)
          call MPI_REDUCE(d_eff_g_r(:,i),d_eff_sum(:,i),mxnr_r,
     &         MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,i_err)
          call MPI_REDUCE(v_p_o_r(:,i),v_p_o_sum(:,i),mxnr_r,
     &         MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,i_err)
          call MPI_REDUCE(v_t_o_r(:,i),v_t_o_sum(:,i),mxnr_r,
     &         MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,i_err)
          call MPI_BCAST(d_effr_sum(:,i),mxnr_r,
     &         MPI_REAL,0,MPI_COMM_WORLD,i_err)
          call MPI_BCAST(d_eff_sum(:,i),mxnr_r,
     &         MPI_REAL,0,MPI_COMM_WORLD,i_err)
          call MPI_BCAST(v_p_o_sum(:,i),mxnr_r,
     &         MPI_REAL,0,MPI_COMM_WORLD,i_err)
          call MPI_BCAST(v_t_o_sum(:,i),mxnr_r,
     &         MPI_REAL,0,MPI_COMM_WORLD,i_err)
        enddo
!
!  Recopy reduced arrays
!
        do i=1,mxnr_r
          do k=1,2
            chi_eff_g_r(i,k)=chi_eff_sum(i,k)
            d_eff_g_r(i,k)=d_eff_sum(i,k)
            d_eff_r(i,k)=d_effr_sum(i,k)
          enddo
          e_rad_r(i,1)=e_rad_mpsum(i)
          e_rad_r(i,3)=e_rad_sum(i,3)
          v_t_o_r(i,2)=v_t_o_sum(i,2)
          v_p_o_r(i,2)=v_p_o_sum(i,2)
        enddo
!
!       write(*,*) 'mx_ms = ',mx_ms,k2
!       do i=1,mxnr_r
!         write(*,50) i, i_proc, rhot_r(i), d_eff_g_r(i,1),
!    &                d_eff_g_r(i,2)
!       enddo
!       do i=2,4
!         write(*,50) i, i_proc, rhot_r(i), v_p_o_r(i,2),
!    &                v_p_o_sum(i,2)
!       enddo
!     stop
 50   format(1x,i2,2x,i2,0p1f10.6,1p2e15.7,' after loop')
!
! ** end loop over grid **
!
!Check convergence of poloidal contribution to Er
        jiter=0
        erref=1.0e3*te_r(1)/a0
        DO i=2,nr_r-1
          IF(ABS((e_rad_r(i,1)-erold1_r(i))/erref).gt.toler) jiter=1
          ertemp_r(i)=0.25*e_rad_r(i,1)+0.5*erold1_r(i)+0.25*erold2_r(i)
          erold2_r(i)=erold1_r(i)
          erold1_r(i)=e_rad_r(i,1)
        ENDDO
!Calculations for all but multiple charge state impurity
        IF(kim2on) THEN
          k2=m_i-1
        ELSE
          k2=m_i
        ENDIF
!  Toroidal rotation constant            
        DO i=2,nr_r-1
          omer=e_rad_r(i,1)*fhat_r(i)/rout_r(i)/btout_r(i)
          DO k=1,k2
            omega_r(i,k)=omegas_r(i,k)+omer
          ENDDO
        ENDDO
!  Axial and edge values
!    Fluxes, transport coefficients and flow velocities
        DO k=1,k2
!         Axis
          d_eff_r(1,k)=d_eff_r(2,k)
          d_n_r(1,k)=d_n_r(2,k)
          v_nt_r(1,k)=0.0
          v_eb_r(1,k)=0.0
          v_ex_r(1,k)=0.0
          d_eff_g_r(1,k)=d_eff_r(1,k)/grho2_r(1)
          d_n_g_r(1,k)=d_n_r(1,k)/grho2_r(1)
          v_nt_g_r(1,k)=0.0
          v_eb_g_r(1,k)=0.0
          v_ex_g_r(1,k)=0.0
          sqz_r(1,k)=sqz_r(2,k)
          u_par_r(1,k)=u_par_r(2,k)
          u_p_r(1,k)=u_p_r(2,k)
          v_t_o_r(1,k)=v_t_o_r(2,k)
          xm_t_o_r(1,k)=xm_t_o_r(2,k)
          v_p_o_r(1,k)=0.0
          xm_p_o_r(1,k)=xm_p_o_r(2,k)
          omega_r(1,k)=omega_r(2,k)
          g_te_r(1,k)=g_te_r(2,k)
          g_ti_r(1,k)=g_ti_r(2,k)
          DO kk=1,k2
            g_n_r(1,kk,k)=g_n_r(2,kk,k)
          ENDDO
!         Edge
          gam_r(nr_r,k)=gam_r(nr_r-1,k)
          d_eff_r(nr_r,k)=d_eff_r(nr_r-1,k)
          d_n_r(nr_r,k)=d_n_r(nr_r-1,k)
          v_nt_r(nr_r,k)=v_nt_r(nr_r-1,k)
          v_eb_r(nr_r,k)=v_eb_r(nr_r-1,k)
          v_ex_r(nr_r,k)=v_ex_r(nr_r-1,k)
          d_eff_g_r(nr_r,k)=d_eff_r(nr_r,k)/grho2_r(nr_r)
          d_n_g_r(nr_r,k)=d_n_r(nr_r,k)/grho2_r(nr_r)
          v_nt_g_r(nr_r,k)=v_nt_r(nr_r,k)/SQRT(grho2_r(nr_r))
          v_eb_g_r(nr_r,k)=v_eb_r(nr_r,k)/SQRT(grho2_r(nr_r))
          v_ex_g_r(nr_r,k)=v_ex_r(nr_r,k)/SQRT(grho2_r(nr_r))
          sqz_r(nr_r,k)=sqz_r(nr_r-1,k)
          u_par_r(nr_r,k)=u_par_r(nr_r-1,k)
          u_p_r(nr_r,k)=u_p_r(nr_r-1,k)
          v_t_o_r(nr_r,k)=v_t_o_r(nr_r-1,k)
          xm_t_o_r(nr_r,k)=xm_t_o_r(nr_r-1,k)
          v_p_o_r(nr_r,k)=v_p_o_r(nr_r-1,k)
          xm_p_o_r(nr_r,k)=xm_p_o_r(nr_r-1,k)
          omega_r(nr_r,k)=omega_r(nr_r-1,k)
          g_te_r(nr_r,k)=g_te_r(nr_r-1,k)
          g_ti_r(nr_r,k)=g_ti_r(nr_r-1,k)
          DO kk=1,k2
            g_n_r(nr_r,kk,k)=g_n_r(nr_r-1,kk,k)
          ENDDO
        ENDDO
!Multiple charge state impurity
!  Axial and edge values
!    Fluxes, transport coefficients and flow velocities
        DO k=1,izim2
          IF(denim2_r(i,k).gt.c_den) THEN
!     Particle fluxes and diffusivities
!           Axis
            d_n_g_im2_r(1,k)=d_n_g_im2_r(2,k)
!           Edge
            gam_im2_r(nr_r,k)=gam_im2_r(nr_r-1,k)
            d_n_g_im2_r(nr_r,k)=d_n_g_im2_r(nr_r-1,k)
            v_nt_g_im2_r(nr_r,k)=v_nt_g_im2_r(nr_r-1,k)
            v_eb_g_im2_r(nr_r,k)=v_eb_g_im2_r(nr_r-1,k)
            v_ex_g_im2_r(nr_r,k)=v_ex_g_im2_r(nr_r-1,k)
          ENDIF
        ENDDO
!Experimental flow velocities for diagnostic impurity
        xm_t_im_ex_o_r(1)=xm_t_im_ex_o_r(2)
        xm_t_im_ex_o_r(nr_r)=xm_t_im_ex_o_r(nr_r-1)
        xm_p_im_ex_o_r(1)=xm_p_im_ex_o_r(2)
        xm_p_im_ex_o_r(nr_r)=xm_p_im_ex_o_r(nr_r-1)
!Electron conduction
!       Axis
        q_con_r(1,1)=0.0
        chi_eff_r(1,1)=chi_eff_r(2,1)
        chi_eff_g_r(1,1)=chi_eff_r(1,1)/grho2_r(1)
!       Edge
        q_con_r(nr_r,1)=q_con_r(nr_r-1,1)
        chi_eff_r(nr_r,1)=chi_eff_r(nr_r-1,1)
        chi_eff_g_r(nr_r,1)=chi_eff_r(nr_r,1)/grho2_r(nr_r)
!Ion conduction
!       Axis
        q_con_r(1,1)=0.0
        chi_eff_r(1,2)=chi_eff_r(2,2)
        chi_eff_g_r(1,2)=chi_eff_r(1,2)/grho2_r(1)
!       Edge
        q_con_r(nr_r,2)=q_con_r(nr_r-1,2)
        chi_eff_r(nr_r,2)=chi_eff_r(nr_r-1,2)
        chi_eff_g_r(nr_r,2)=chi_eff_r(nr_r,2)/grho2_r(nr_r)
!Parallel electrical resistivity
        eta_par_r(1)=eta_par_r(2)
        eta_par_r(nr_r)=eta_par_r(nr_r-1)
     &                  +(eta_par_r(nr_r-1)-eta_par_r(nr_r-2))
     &                  *(rhot_r(nr_r)  -rhot_r(nr_r-1))
     &                  /(rhot_r(nr_r-1)-rhot_r(nr_r-2))
!Bootstrap current density
        xj_bs_r(nr_r)=xj_bs_r(nr_r-1)
     &                +(xj_bs_r(nr_r-1)-xj_bs_r(nr_r-2))
     &                *(rhot_r(nr_r)  -rhot_r(nr_r-1))
     &                /(rhot_r(nr_r-1)-rhot_r(nr_r-2))
!Parallel electric field
        e_par_r(1)=e_par_r(2)
        e_par_r(nr_r)=e_par_r(nr_r-1)
     &                +(e_par_r(nr_r-1)-e_par_r(nr_r-2))
     &                *(rhot_r(nr_r)  -rhot_r(nr_r-1))
     &                /(rhot_r(nr_r-1)-rhot_r(nr_r-2))
!Radial electric field
        e_rad_r(nr_r,1)=e_rad_r(nr_r-1,1)
        e_rad_r(nr_r,2)=e_rad_r(nr_r-1,2)
        e_rad_r(nr_r,3)=e_rad_r(nr_r-1,3)
        e_rad_r(nr_r,4)=e_rad_r(nr_r-1,4)
        e_rad_o_r(nr_r,1)=e_rad_o_r(nr_r-1,1)
        e_rad_o_r(nr_r,2)=e_rad_o_r(nr_r-1,2)
        e_rad_o_r(nr_r,3)=e_rad_o_r(nr_r-1,3)
        e_rad_o_r(nr_r,4)=e_rad_o_r(nr_r-1,4)
!Calculate ExB shear damping on outside
        DO i=3,nr_r-1
          omexb_o_r(i)=(rout_r(i)*bpout_r(i))**2
     &                 /SQRT(btout_r(i)**2+bpout_r(i)**2)
     &                 *(e_rad_o_r(i+1,1)/rout_r(i+1)/bpout_r(i+1)
     &                  -e_rad_o_r(i-1,1)/rout_r(i-1)/bpout_r(i-1))
     &                 /(psi_r(i+1)-psi_r(i-1))
          omexb_ex_o_r(i)=(rout_r(i)*bpout_r(i))**2
     &                 /SQRT(btout_r(i)**2+bpout_r(i)**2)
     &                 *(e_rad_ex_o_r(i+1)/rout_r(i+1)/bpout_r(i+1)
     &                  -e_rad_ex_o_r(i-1)/rout_r(i-1)/bpout_r(i-1))
     &                 /(psi_r(i+1)-psi_r(i-1))
        ENDDO
        omexb_o_r(2)=omexb_o_r(3)
        omexb_o_r(1)=0.0
        omexb_o_r(nr_r)=omexb_o_r(nr_r-1)
        omexb_ex_o_r(2)=omexb_ex_o_r(3)
        omexb_ex_o_r(1)=0.0
        omexb_ex_o_r(nr_r)=omexb_ex_o_r(nr_r-1)
        IF(jiter.eq.0) GOTO 1000
      ENDDO
 1000 continue
c     do i=1,20
c       write(*,150) i, rhot_r(i), q_con_r(i,1)
c     enddo
c       do i=2,4
c         write(*,60) i, i_proc, rhot_r(i), v_p_o_r(i,2)
c       enddo
c      stop
 60   format(2x,i2,2x,i2,2x,0p1f10.6,1p6e15.7)
 150  format(2x,i2,2x,0p1f10.6,1p6e13.5)
!1000 RETURN
      END
      SUBROUTINE READ_PRO_DIIID(idshot,time,bt0,izi0,izim1,izim2,nr_r,           
     &                          mxnr_r,rhot_r,te_r,ti_r,den_r,deni_r,
     &                          zeff_ex_r,denim_r,denim2_r,omega_ex_r,
     &                          vp_im_ex_r,vt_im_ex_r,xk_im_ex_r,
     &                          zeff_r,e_par_ex_r,xj_ex_r,xj_nb_ex_r,
     &                          iflag,message)
!***********************************************************************
!READ_PRO_DIIID reads 4D data files to get plasma profiles
!  W.A.Houlberg 3/2000
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        idshot,                  izi0(*),
     &               izim1,                   izim2,
     &               mxnr_r,                  nr_r
      REAL           bt0,                     rhot_r(*),
     &               time
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           den_r(*),                deni_r(mxnr_r,*),
     &               denim_r(*),              denim2_r(mxnr_r,*),
     &               e_par_ex_r(*),           omega_ex_r(*),
     &               te_r(*),                 ti_r(*),
     &               vp_im_ex_r(*),           vt_im_ex_r(*),
     &               zeff_ex_r(*),            zeff_r(*),
     &               xj_ex_r(*),              xj_nb_ex_r(*),
     &               xk_im_ex_r(*)
!Declaration of local variables
!  Profiles
      INTEGER        mxnr_4d,                 nr_4d
      PARAMETER     (mxnr_4d=101)
      REAL           r_4d(mxnr_4d),           y_4d(mxnr_4d)         
      INTEGER        mxn
      PARAMETER     (mxn=300)
      REAL           denb_ex_r(mxn)         
!  Other
      CHARACTER      cidshot*9,               cnin*25,
     &               ctimems*5
      INTEGER        i,                       ierr,
     &               itimems,                 jpro,
     &               k,                       nin
      REAL           fscale,                  yim,
     &               yim2
      nin=20
      CALL RARRAY_ZERO(nr_r,te_r)
      CALL RARRAY_ZERO(nr_r,ti_r)
      CALL RARRAY_ZERO(nr_r,den_r)
      CALL RARRAY_ZERO(nr_r,zeff_ex_r)
      CALL RARRAY_ZERO(nr_r,omega_ex_r)
      CALL RARRAY_ZERO(nr_r,denim_r)
      CALL RARRAY_ZERO(mxnr_r*izim2,denim2_r)
      CALL RARRAY_ZERO(nr_r,vt_im_ex_r)
      CALL RARRAY_ZERO(nr_r,vp_im_ex_r)
      CALL RARRAY_ZERO(nr_r,xk_im_ex_r)
      CALL RARRAY_ZERO(nr_r,e_par_ex_r)
      CALL RARRAY_ZERO(nr_r,xj_ex_r)      
      CALL RARRAY_ZERO(nr_r,xj_nb_ex_r)      
!Set up file naming convention
      itimems=1.0e3*(time+.00001)
      WRITE(ctimems,'(i5)') itimems
      DO i=1,5
        IF(ctimems(i:i).eq.' ') ctimems(i:i)='0'
      ENDDO
      WRITE(cidshot,'(i9)') idshot
      DO i=1,9
        IF(cidshot(i:i).eq.' ') cidshot(i:i)='0'
      ENDDO
!Read data
      DO jpro=1,14
        IF(jpro.eq.1) THEN
!         Electron temperature
          cnin='te'//cidshot(4:9)//'.'//ctimems
          ierr=1
          fscale=1.0
        ELSEIF(jpro.eq.2) THEN
!         Ion temperature
          cnin='ti'//cidshot(4:9)//'.'//ctimems
          ierr=1
          fscale=1.0
        ELSEIF(jpro.eq.3) THEN
!         Electron density
          cnin='ne'//cidshot(4:9)//'.'//ctimems
          ierr=1
          fscale=1.0e19
        ELSEIF(jpro.eq.4) THEN
!         Zeff
          cnin='zf'//cidshot(4:9)//'.'//ctimems
          ierr=0
          fscale=1.0
        ELSEIF(jpro.eq.5) THEN
!         Impurity flux surface average toroidal rotation frequency
          cnin='om'//cidshot(4:9)//'.'//ctimems
          ierr=0
          fscale=1.0
        ELSEIF(jpro.eq.6) THEN
!         Impurity density
          cnin='nc'//cidshot(4:9)//'.'//ctimems
          ierr=0
          fscale=1.0e19
        ELSEIF(jpro.eq.7) THEN
!         Impurity toroidal rotation velocity on outside
          cnin='vt'//cidshot(4:9)//'.'//ctimems
          ierr=0
          fscale=1.0
        ELSEIF(jpro.eq.8) THEN
!         Impurity poloidal rotation velocity on outside
          cnin='vp'//cidshot(4:9)//'.'//ctimems
          ierr=0
          fscale=1.0
        ELSEIF(jpro.eq.9) THEN
!         Impurity flux surface average toroidal rotation frequency
          cnin='kk'//cidshot(4:9)//'.'//ctimems
          ierr=0
          fscale=1.0
        ELSEIF(jpro.eq.10) THEN
!         Impurity 2 density - usually neon
          cnin='nn'//cidshot(4:9)//'.'//ctimems
          ierr=0
          fscale=1.0e19
        ELSEIF(jpro.eq.11) THEN
!         <E.B>
          cnin='eb'//cidshot(4:9)//'.'//ctimems
          ierr=0
          fscale=1.0/bt0
        ELSEIF(jpro.eq.12) THEN
!         <J.B>
          cnin='jb'//cidshot(4:9)//'.'//ctimems
          ierr=0
          fscale=1.0/bt0
        ELSEIF(jpro.eq.13) THEN
!         <J_NBI.B>
          cnin='jnb'//cidshot(4:9)//'.'//ctimems
          ierr=0
          fscale=1.0/bt0
        ELSEIF(jpro.eq.14) THEN
!         n_NBI
          cnin='nb'//cidshot(4:9)//'.'//ctimems
          ierr=0
          fscale=1.0e19
        ENDIF
        CALL READ_REC(nin,cnin,nr_4d)
        IF(nr_4d.eq.0.and.ierr.eq.1) THEN
!         Critical data missing
          iflag=1
          message='READ_PRO_DIIID(1)/ERROR:missing file '//cnin
          GOTO 1000
        ELSEIF(nr_4d.eq.0.and.ierr.eq.0) THEN
!         Non-critical data missing
          iflag=0
        ELSEIF(nr_4d.gt.mxnr_4d) THEN
          iflag=1
          message='READ_PRO_DIIID(2)/WARNING:grid dim in file'//cnin
          GOTO 1000
        ELSE
          OPEN(unit=nin,status='old',file=cnin,form='formatted')
          DO i=1,nr_4d
            READ(nin,*) r_4d(i),y_4d(i)
          ENDDO
          CLOSE(unit=nin)
          CALL RARRAY_SCALE(nr_4d,fscale,y_4d,1)
          IF(jpro.eq.1) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,te_r,iflag,
     &                        message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(te_r)/'//message(1:i-21)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.2) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,ti_r,iflag,
     &                        message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(ti_r)/'//message(1:i-21)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.3) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,den_r,iflag,
     &                        message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(den_r)/'//message(1:i-22)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.4) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,zeff_ex_r,
     &                        iflag,message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(zeff_ex_r)/'//message(1:i-26)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.5) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,omega_ex_r,
     &                        iflag,message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(omega_ex_r)/'//message(1:i-27)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.6) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,denim_r,iflag,
     &                        message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(denim_r)/'//message(1:i-24)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.7) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,vt_im_ex_r,
     &                        iflag,message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(vt_im_ex_r)/'//message(1:i-27)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.8) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,vp_im_ex_r,
     &                        iflag,message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(vp_im_ex_r)/'//message(1:i-27)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.9) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,xk_im_ex_r,
     &                        iflag,message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(xk_im_ex_r)/'//message(1:i-27)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.10) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,
     &                        denim2_r(1,izim2),iflag,message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(denim2_r/'//message(1:i-25)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.11) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,e_par_ex_r,
     &                        iflag,message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(e_par_ex_r)/'//message(1:i-27)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.12) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,xj_ex_r,iflag,
     &                        message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(xj_ex_r)/'//message(1:i-24)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.13) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,xj_nb_ex_r,
     &                        iflag,message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(xj_nb_ex_r)/'//message(1:i-27)
              GOTO 1000
            ENDIF
          ELSEIF(jpro.eq.14) THEN
            CALL W_LIN_INTERP(nr_4d,r_4d,y_4d,nr_r,rhot_r,denb_ex_r,
     &                        iflag,message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='READ_PRO_DIIID(denb_ex_r)/'//message(1:i-26)
              GOTO 1000
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!For other 'test' species, use a small percentage of electron density
      DO k=2,5
        IF(izi0(k).ne.0) THEN
          DO i=1,nr_r
            deni_r(i,k)=1.0e-5*den_r(i)
          ENDDO
        ENDIF
      ENDDO      
!Break electron density into ion and impurity components
      DO i=1,nr_r
!       Diagnostic impurity
        yim=izim1*denim_r(i)
        yim2=izim1**2*denim_r(i)
!       Other than main ion
        DO k=2,5
          yim=yim+izi0(k)*deni_r(i,k)
          yim2=yim2+izi0(k)**2*deni_r(i,k)
        ENDDO
!       Multiple charge state impurity
        DO k=1,izim2
          yim=yim+k*denim2_r(i,k)
          yim2=yim2+k**2*denim2_r(i,k)
        ENDDO
        deni_r(i,1)=(den_r(i)-yim-denb_ex_r(i))/izi0(1)
        zeff_r(i)=(izi0(1)*deni_r(i,1)+yim2+denb_ex_r(i))/den_r(i)
      ENDDO
 1000 RETURN
      END
      SUBROUTINE READ_PRO_VMEC(nr_r,mxnr_r,rhot_r,te_r,ti_r,den_r,
     &                         deni_r,denim_r,denim2_r,zeff_r,iflag,
     &                         message)
!***********************************************************************
!READ_PRO_VMEC reads a data file get plasma profiles for benchmarking
!  W.A.Houlberg 3/2000
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        mxnr_r,                  nr_r
      REAL           rhot_r(*)
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           den_r(*),                deni_r(mxnr_r,*),
     &               denim_r(*),              denim2_r(mxnr_r,*),
     &               te_r(*),                 ti_r(*),
     &               zeff_r(*)
!Declaration of local variables
!  Profiles
      CHARACTER*4    adum
      INTEGER        mxnr_vm,                 nr_vm
      PARAMETER     (mxnr_vm=32)
      REAL           y_vm(mxnr_vm,10)         
!  Other
      INTEGER        i,                       j,
     &               k,                       nin
      REAL           delphi
!Declaration of external functions
      REAL           W_POLY_INTERP
      nin=20
      nr_vm=mxnr_vm
      CALL RARRAY_ZERO(nr_r,te_r)
      CALL RARRAY_ZERO(nr_r,ti_r)
      CALL RARRAY_ZERO(nr_r,den_r)
      CALL RARRAY_ZERO(nr_r,zeff_r)
      CALL RARRAY_ZERO(nr_r,denim_r)
      CALL RARRAY_ZERO(nr_r,denim2_r)
!Open data file
      OPEN(unit=nin,
     &     status='old',
     &     file='sympp.dat',
     &     form='formatted',
     &     iostat=iflag)
      IF(iflag.ne.0) THEN
        iflag=1
        message='READ_PRO_VMEC/ERROR:file sympp.dat'
        GOTO 1000
      ENDIF
!Read data
!  Uniform mesh in toroidal flux
!  Convert to square root toroidal flux
!  Starts and ends a half mesh point from the origin and edge
      delphi=1.0/(nr_vm-2)
      READ(nin,'(a4)') adum
      DO i=2,nr_vm-1
        READ(nin,*) j,(y_vm(i,k),k=1,10)
      ENDDO
!Extrapolate values to axis and edge
      DO k=1,10
        y_vm(1,k)=1.5*y_vm(2,k)-0.5*y_vm(3,k)
        y_vm(nr_vm,k)=1.5*y_vm(nr_vm-1,k)-0.5*y_vm(nr_vm-2,k)
      ENDDO
      y_vm(1,6)=0.0
      y_vm(1,10)=0.0
!Convert grid to (Phi/Phi_a)**0.5
      DO i=2,nr_vm-1
        y_vm(i,7)=SQRT(y_vm(i,7))
      ENDDO
      y_vm(1,7)=0.0
      y_vm(nr_vm,7)=1.0
!Interpolate from VMEC grid to working grid
      DO i=1,nr_r
        te_r(i)=W_POLY_INTERP(y_vm(1,7),y_vm(1,1),nr_vm,3,rhot_r(i))
        ti_r(i)=W_POLY_INTERP(y_vm(1,7),y_vm(1,2),nr_vm,3,rhot_r(i))
        den_r(i)=W_POLY_INTERP(y_vm(1,7),y_vm(1,3),nr_vm,3,rhot_r(i))
      ENDDO
!Electron density is normalized to 10^20/m**3, rescale
      CALL RARRAY_SCALE(nr_r,1.0e20,den_r,1)
!Ion density
      CALL RARRAY_COPY(nr_r,den_r,1,deni_r,1)
!Zeff=1
      DO i=1,nr_r
        zeff_r(i)=1.0
      ENDDO
 1000 CLOSE(unit=nin)
      RETURN
      END
      SUBROUTINE READ_REC(nin,cnin,nrec)
!***********************************************************************
!READ_REC reads a data file to determine the number of records 
!  W.A. Houlberg 6/97
!Input:
!  nin-input file unit number
!  cnin-input file name
!Output:
!  nrec-number of records
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      CHARACTER*(*)  cnin
      INTEGER        nin
!Declaration of output variables
      INTEGER        nrec
!Declaration of local variables
      INTEGER        iflag
      REAL           x
      nrec=0
      iflag=0
      OPEN(unit=nin,status='old',file=cnin,form='formatted',
     &     iostat=iflag)
      IF(iflag.ne.0) GOTO 1000
      DO WHILE (iflag.eq.0)
        READ(nin,*,iostat=iflag) x
        nrec=nrec+1
      ENDDO
      nrec=nrec-1
 1000 CLOSE(unit=nin)
      RETURN
      END
      SUBROUTINE READ_PRO_START(idshot,time,izi0,izim1,izim2,nr_r,
     &                          mxnr_r,rout_r,te_r,ti_r,den_r,deni_r,
     &                          denim_r,denim2_r,vt_im_ex_r,zeff_r,
     &                          iflag,message)
!***********************************************************************
!READ_PRO_START reads START u-file profile data then interpolates onto
!  the working grid
!  W.A.Houlberg 3/2000
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        idshot,                  izi0(*),                  
     &               izim1,                   izim2,
     &               mxnr_r,                  nr_r
      REAL           time,                    rout_r(*)
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag,                   nin
      REAL           te_r(*),                 ti_r(*),
     &               den_r(*),                deni_r(mxnr_r,*),
     &               denim_r(*),              vt_im_ex_r(*),
     &               zeff_r(*),               denim2_r(mxnr_r,*)
!Declaration of local variables
!  Profiles
      INTEGER        mxnr_uf
      PARAMETER     (mxnr_uf=101)
      REAL           r_uf(mxnr_uf),           y_uf(mxnr_uf)
!  Other
      CHARACTER      cidshot*6,               cnin*25,
     &               ctimems*5
      INTEGER        i,                       icid,
     &               ierr,                    itimems,
     &               jpro,                    nr_uf
      REAL           fscale,                  xidshot
      nin=20
      CALL RARRAY_ZERO(nr_r,te_r)
      CALL RARRAY_ZERO(nr_r,ti_r)
      CALL RARRAY_ZERO(nr_r,den_r)
      CALL RARRAY_ZERO(nr_r,denim_r)
      CALL RARRAY_ZERO(nr_r,vt_im_ex_r)
      izim2=0
      CALL RARRAY_ZERO(mxnr_r,denim2_r)
!Set up file naming convention
      itimems=1.0e3*(time+.00001)
      WRITE(ctimems,'(i5)') itimems
      DO i=1,5
        IF(ctimems(i:i).eq.' ') ctimems(i:i)='0'
      ENDDO
      WRITE(cidshot,'(i6)') idshot
      xidshot=idshot
      icid=7-LOG10(xidshot)
!Read data
      DO jpro=1,4
        IF(jpro.eq.1) THEN
!         Electron temperature
          cnin='f'//cidshot(icid:6)//'.ter'
          ierr=1
          fscale=1.0e-3
        ELSEIF(jpro.eq.2) THEN
!         Ion temperature
          cnin='f'//cidshot(icid:6)//'.ti2'
          ierr=1
          fscale=1.0e-3
        ELSEIF(jpro.eq.3) THEN
!         Electron density
          cnin='f'//cidshot(icid:6)//'.ner'
          ierr=1
          fscale=1.0e6
        ELSEIF(jpro.eq.4) THEN
!         Impurity toroidal rotation velocity on outside
          cnin='f'//cidshot(icid:6)//'.vp2'
          ierr=0
          fscale=1.0e-2
        ENDIF
        CALL READ_UF_2D(nin,cnin,mxnr_uf,time,nr_uf,r_uf,y_uf,iflag,
     &                  message)
        IF(iflag.ne.0.and.ierr.eq.1) THEN
          iflag=1
          i=LEN(message)
          message='READ_PRO_START(1)/'//message(1:i-18)
          GOTO 1000
        ELSE
          CALL RARRAY_SCALE(nr_uf,fscale,y_uf,1)
          IF(jpro.eq.1) THEN
            CALL W_LIN_INTERP(nr_uf,r_uf,y_uf,nr_r,rout_r,te_r,iflag,
     &                        message)
          ELSEIF(jpro.eq.2) THEN
            CALL W_LIN_INTERP(nr_uf,r_uf,y_uf,nr_r,rout_r,ti_r,iflag,
     &                        message)
          ELSEIF(jpro.eq.3) THEN
            CALL W_LIN_INTERP(nr_uf,r_uf,y_uf,nr_r,rout_r,den_r,iflag,
     &                        message)
          ELSEIF(jpro.eq.4) THEN
            CALL W_LIN_INTERP(nr_uf,r_uf,y_uf,nr_r,rout_r,vt_im_ex_r,
     &                        iflag,message)
          ENDIF
          IF(iflag.ne.0) THEN
            i=LEN(message)
            message='READ_PRO_START(2)/'//message(1:i-18)
            IF(iflag.eq.1) GOTO 1000
          ENDIF
        ENDIF
      ENDDO
!Single charge state impurity density
!     Temporarily scaled with electron density
      CALL RARRAY_COPY(nr_r,den_r,1,denim_r,1)
      fscale=1.0/(izim1*(izim1-1))
      CALL RARRAY_SCALE(nr_r,fscale,denim_r,1)
!Set main ion density from electron and impurity densities
!     Impurity has single charge state izim1
      DO i=1,nr_r
        deni_r(i,1)=(den_r(i)-izim1*denim_r(i))/izi0(1)
        zeff_r(i)=(izi0(1)*deni_r(i,1)+izim1**2*denim_r(i))/den_r(i)
      ENDDO
 1000 RETURN
      END
      SUBROUTINE READ_PRO_TEXTOR(idshot,time,izi0,izim1,izim2,nr_r,
     &                           mxnr_r,rout_r,te_r,ti_r,den_r,deni_r,
     &                           denim_r,denim2_r,vt_im_ex_r,zeff_r,
     &                           iflag,message)
!***********************************************************************
!READ_PRO_TEXTOR reads TEXTOR u-file profile data then interpolates onto
!  the working grid
!  W.A. Houlberg 3/2000
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        idshot,                  izi0(*),                  
     &               izim1,                   izim2,
     &               mxnr_r,                  nr_r
      REAL           time,                    rout_r(*)
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag,                   nin
      REAL           te_r(*),                 ti_r(*),
     &               den_r(*),                deni_r(mxnr_r,*),
     &               denim_r(*),              vt_im_ex_r(*),
     &               zeff_r(*),               denim2_r(mxnr_r,*)
!Declaration of local variables
!  Profiles
      INTEGER        mxnr_uf
      PARAMETER     (mxnr_uf=101)
      REAL           r_uf(mxnr_uf),           y_uf(mxnr_uf)
!  Other
      CHARACTER      cidshot*6,               cnin*25,
     &               ctimems*5
      INTEGER        i,                       icid,
     &               ierr,                    itimems,
     &               jpro,                    nr_uf
      REAL           fscale,                  xidshot
      nin=20
      CALL RARRAY_ZERO(nr_r,te_r)
      CALL RARRAY_ZERO(nr_r,ti_r)
      CALL RARRAY_ZERO(nr_r,den_r)
      CALL RARRAY_ZERO(nr_r,denim_r)
      CALL RARRAY_ZERO(nr_r,vt_im_ex_r)
      izim2=0
      CALL RARRAY_ZERO(mxnr_r,denim2_r)
!Set up file naming convention
      itimems=1.0e3*(time+.00001)
      WRITE(ctimems,'(i5)') itimems
      DO i=1,5
        IF(ctimems(i:i).eq.' ') ctimems(i:i)='0'
      ENDDO
      WRITE(cidshot,'(i6)') idshot
      xidshot=idshot
      icid=7-LOG10(xidshot)
!Read data
      DO jpro=1,4
        IF(jpro.eq.1) THEN
!         Electron temperature
          cnin='f'//cidshot(icid:6)//'.ter'
          ierr=1
          fscale=1.0e-3
        ELSEIF(jpro.eq.2) THEN
!         Ion temperature
          cnin='f'//cidshot(icid:6)//'.ti2'
          ierr=1
          fscale=1.0e-3
        ELSEIF(jpro.eq.3) THEN
!         Electron density
          cnin='f'//cidshot(icid:6)//'.ner'
          ierr=1
          fscale=1.0e6
        ELSEIF(jpro.eq.4) THEN
!         Impurity toroidal rotation velocity on outside
          cnin='f'//cidshot(icid:6)//'.vp2'
          ierr=0
          fscale=1.0e-2
        ENDIF
        CALL READ_UF_2D(nin,cnin,mxnr_uf,time,nr_uf,r_uf,y_uf,iflag,
     &                  message)
        IF(iflag.ne.0.and.ierr.eq.1) THEN
          iflag=1
          i=LEN(message)
          message='READ_PRO_TEXTOR(1)/'//message(1:i-19)
          GOTO 1000
        ELSE
          CALL RARRAY_SCALE(nr_uf,fscale,y_uf,1)
          IF(jpro.eq.1) THEN
            CALL W_LIN_INTERP(nr_uf,r_uf,y_uf,nr_r,rout_r,te_r,iflag,
     &                        message)
          ELSEIF(jpro.eq.2) THEN
            CALL W_LIN_INTERP(nr_uf,r_uf,y_uf,nr_r,rout_r,ti_r,iflag,
     &                        message)
          ELSEIF(jpro.eq.3) THEN
            CALL W_LIN_INTERP(nr_uf,r_uf,y_uf,nr_r,rout_r,den_r,iflag,
     &                        message)
          ELSEIF(jpro.eq.4) THEN
            CALL W_LIN_INTERP(nr_uf,r_uf,y_uf,nr_r,rout_r,vt_im_ex_r,
     &                        iflag,message)
          ENDIF
          IF(iflag.ne.0) THEN
            i=LEN(message)
            message='READ_PRO_TEXTOR(2)/'//message(1:i-19)
            IF(iflag.eq.1) GOTO 1000
          ENDIF
        ENDIF
      ENDDO
!Single charge state impurity density
!     Temporarily scaled with electron density
      CALL RARRAY_COPY(nr_r,den_r,1,denim_r,1)
      fscale=1.0/(izim1*(izim1-1))
      CALL RARRAY_SCALE(nr_r,fscale,denim_r,1)
!Set main ion density from electron and impurity densities
!     Impurity has single charge state izim1
      DO i=1,nr_r
        deni_r(i,1)=(den_r(i)-izim1*denim_r(i))/izi0(1)
        zeff_r(i)=(izi0(1)*deni_r(i,1)+izim1**2*denim_r(i))/den_r(i)
      ENDDO
 1000 RETURN
      END
      SUBROUTINE READ_UF_2D(nin,cnin,nrmax,time,nr,r,v,iflag,message)
!***********************************************************************
!READ_UF_2D reads a 2D U-file (r,t) to get the first radial slice
!  following the requested time
!  W.A.Houlberg 3/2000
!Input:
!  nin-input file unit number
!  cnin-input file name
!  nrmax-maximum number of radial points
!  time-time of data requested (s)
!Output:
!  nr-number of radial points
!  r(i)-radial grid (m)
!  v(i)-radial values
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      CHARACTER*(*)  cnin
      INTEGER        nin,                     nrmax
      REAL           time
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag,                   nr
      REAL           r(*),                    v(*)
!Declaration of local variables
      INTEGER        mxnt,                    mxnr
      PARAMETER      (mxnt=2000,              mxnr=400)
      CHARACTER*72   line
      INTEGER        i,                       j,
     &               ifollow,                 nt
      REAL           t(mxnt),                 z(mxnr,mxnt)
      CALL RARRAY_ZERO(mxnr,r)
      CALL RARRAY_ZERO(mxnr,v)
      OPEN(unit=nin,
     &     file=cnin,
     &     status='old',
     &     iostat=iflag)
      IF(iflag.ne.0) THEN
        iflag=1
        message='READ_UF_2D(1)/ERROR:opening'//cnin
        GOTO 1000
      ENDIF
      ifollow=0
      DO WHILE (ifollow.eq.0)
        READ(nin,'(a72)',iostat=iflag) line
        IF(iflag.ne.0) THEN
          iflag=1
          message='READ_UF_2D(2)/ERROR:reading line in '//cnin
          GOTO 1000
        ENDIF
        IF(INDEX(line,'# OF X PTS').gt.0) THEN
          READ(line,'(1x,i10)',iostat=iflag) nt
          IF(iflag.ne.0) THEN
            iflag=1
            message='READ_UF_2D(3)/ERROR:reading nt in '//cnin
            GOTO 1000
          ELSEIF(nt.gt.mxnt) THEN
            iflag=1
            message='READ_UF_2D(4)/ERROR:too many nt in '//cnin
            GOTO 1000
          ENDIF
        ENDIF
        IF(INDEX(line,'# OF Y PTS').gt.0) THEN
          READ(line,'(1x,i10)',iostat=iflag) nr
          IF(iflag.ne.0) THEN
            iflag=1
            message='READ_UF_2D(5)/ERROR:reading nr in '//cnin
            GOTO 1000
          ELSEIF(nr.gt.mxnr.or.nr.gt.nrmax) THEN
            iflag=1
            message='READ_UF_2D(6)/ERROR:too many nr in '//cnin
            GOTO 1000
          ENDIF
        ENDIF
        ifollow=INDEX(line,'FOLLOW')
      ENDDO
!Time
      READ(nin,'(6(1x,1pe11.4))',iostat=iflag) (t(i),i=1,nt)
      IF(iflag.ne.0) THEN
        iflag=1
        message='READ_UF_2D(7)/ERROR:not enough nt in '//cnin
        GOTO 1000
      ENDIF
      CALL RARRAY_SCALE(nt,1.0e-3,t,1)
!Radius
      READ(nin,'(6(1x,1pe11.4))',iostat=iflag) (r(i),i=1,nr)
      IF(iflag.ne.0) THEN
        iflag=1
        message='READ_UF_2D(8)/ERROR:not enough nr in '//cnin
        GOTO 1000
      ENDIF
      CALL RARRAY_SCALE(nr,1.0e-2,r,1)
!Values
      READ(nin,'(6(1x,1pe11.4))',iostat=iflag) ((z(j,i),i=1,nt),j=1,nr)
      IF(iflag.ne.0) THEN
        iflag=1
        message='READ_UF_2D(9)/ERROR:data error in'//cnin
        GOTO 1000
      ENDIF
!Find index following requested time
      i=i
      DO WHILE (time.lt.t(i))
        i=i+1
        IF(i.gt.nt) THEN
          iflag=1
          message='READ_UF_2D(7)/ERROR:t out of range in '//cnin
          GOTO 1000
        ENDIF
      ENDDO
!Copy profile
      CALL RARRAY_COPY(nr,z(1,i-1),1,v,1)
 1000 CLOSE(nin)
      RETURN
      END
      SUBROUTINE SPECIES(cs,ics,izs,amus,iflag,message)
!***********************************************************************
!SPECIES sets the mass number and atomic charge of a species identified
!  by its name
!  W.A.Houlberg 3/2000
!Input:
!  cs-species name
!Output:
!  ics-number of non-null characters in cs (-)
!  izs-nuclear charge (-)
!  amus-atomic mass number (-)
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      CHARACTER*(*)  cs
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        ics,                     iflag,
     &               izs
      REAL           amus
!Declaration of local variables
      INTEGER        mxd
      PARAMETER      (mxd=98)
      CHARACTER*3    cd(mxd)
      INTEGER        icd(mxd),                izd(mxd)
      INTEGER        i,                       idone
      REAL           amud(mxd)
      SAVE           cd,                      icd,
     &               izd,                     amud
      DATA cd/           'H',      'D',      'T',    'He3',    'He4',
     &                  'Li',     'Be',      'B',      'C',      'N',
     &                   'O',      'F',     'Ne',     'Na',     'Mg',
     &                  'Al',     'Si',      'P',      'S',     'Cl',
     &                   'A',      'K',     'Ca',     'Sc',     'Ti',
     &                   'V',     'Cr',     'Mn',     'Fe',     'Co',
     &                  'Ni',     'Cu',     'Zn',     'Ga',     'Ge',
     &                  'As',     'Se',     'Br',     'Kr',     'Rb',
     &                  'Sr',      'Y',     'Zr',     'Nb',     'Mo',
     &                  'Tc',     'Ru',     'Rh',     'Pd',     'Ag',
     &                  'Cd',     'In',     'Sn',     'Sb',     'Te',
     &                   'I',     'Xe',     'Cs',     'Ba',     'La',
     &                  'Ce',     'Pr',     'Nd',     'Pm',     'Sm',
     &                  'Eu',     'Gd',     'Tb',     'Dy',     'Ho',
     &                  'Er',     'Tm',     'Yb',     'Lu',     'Hf',
     &                  'Ta',      'W',     'Re',     'Os',     'Ir',
     &                  'Pt',     'Au',     'Hg',     'Tl',     'Pb',
     &                  'Bi',     'Po',     'At',     'Rn',     'Fr',
     &                  'Ra',     'Ac',     'Th',     'Pa',      'U',
     &                  'Np',     'Pu',     'Ar'/
      DATA icd/            1,        1,        1,        3,        3,
     &                     2,        2,        1,        1,        1,
     &                     1,        1,        2,        2,        2,
     &                     2,        2,        1,        1,        2,
     &                     1,        1,        2,        2,        2,
     &                     1,        2,        2,        2,        2,
     &                     2,        2,        2,        2,        2,
     &                     2,        2,        2,        2,        2,
     &                     2,        1,        2,        2,        2,
     &                     2,        2,        2,        2,        2,
     &                     2,        2,        2,        2,        2,
     &                     1,        2,        2,        2,        2,
     &                     2,        2,        2,        2,        2,
     &                     2,        2,        2,        2,        2,
     &                     2,        2,        2,        2,        2,
     &                     2,        1,        2,        2,        2,
     &                     2,        2,        2,        2,        2,
     &                     2,        2,        2,        2,        2,
     &                     2,        2,        2,        2,        1,
     &                     2,        2,        2/
      DATA izd/            1,        1,        1,        2,        2,
     &                     3,        4,        5,        6,        7,
     &                     8,        9,       10,       11,       12,
     &                    13,       14,       15,       16,       17,
     &                    18,       19,       20,       21,       22,
     &                    23,       24,       25,       26,       27,
     &                    28,       29,       30,       31,       32,
     &                    33,       34,       35,       36,       37,
     &                    38,       39,       40,       41,       42,
     &                    43,       44,       45,       46,       47,
     &                    48,       49,       50,       51,       52,
     &                    53,       54,       55,       56,       57,
     &                    58,       59,       60,       61,       62,
     &                    63,       64,       65,       66,       67,
     &                    68,       69,       70,       71,       72,
     &                    73,       74,       75,       76,       77,
     &                    78,       79,       80,       81,       82,
     &                    83,       84,       85,       86,       87,
     &                    88,       89,       90,       91,       92,
     &                    93,       94,       18/
      DATA amud/       1.000,    2.000,    3.000,    3.000,    4.000,
     &                 6.939,    9.012,   10.811,   12.011,   14.007,
     &                15.999,   18.998,   20.183,   22.990,   24.312,
     &                26.982,   28.086,   30.974,   32.064,   35.453,
     &                39.948,   39.102,   40.080,   44.956,   47.900,
     &                50.942,   51.996,   54.938,   55.847,   58.933,
     &                58.710,   63.540,   65.370,   69.720,   72.590,
     &                74.922,   78.960,   79.909,   83.800,   85.470,
     &                87.620,   88.905,   91.220,   92.906,   95.940,
     &                99.000,  101.070,  102.905,  106.400,  107.870,
     &               112.400,  114.820,  118.690,  121.750,  127.600,
     &               126.904,  131.300,  132.905,  137.340,  138.910,
     &               140.120,  140.907,  144.240,  145.000,  150.350,
     &               151.960,  157.250,  158.924,  162.500,  164.930,
     &               167.260,  168.934,  173.040,  174.970,  178.490,
     &               180.948,  183.850,  186.200,  190.200,  192.200,
     &               195.090,  196.967,  200.590,  204.370,  207.190,
     &               208.980,  210.000,  210.000,  222.000,  223.000,
     &               226.000,  227.000,  232.038,  231.000,  238.030,
     &               237.000,  242.000,   39.948/
!Set defaults
      ics=0
      izs=0
      amus=0.0
!Find number of characters in cs by looking for first null
      i=1
      idone=0
      DO WHILE(idone.eq.0.and.i.le.3)
        IF(cs(i:i).ne.' ') THEN
          ics=i
        ELSE
          idone=1
        ENDIF
        i=i+1
      ENDDO
!Check for null input
      IF(ics.eq.0) THEN
        iflag=-1
        message='SPECIES/WARNING:null 1st character of species name'
        GOTO 1000
      ENDIF
      i=1
      idone=0
      DO WHILE(idone.eq.0.and.i.le.mxd)
        IF(ics.eq.icd(i)) THEN
          IF(cs(:ics).eq.cd(i)(:ics)) THEN
            izs=izd(i)
            amus=amud(i)
            idone=1
          ENDIF
        ENDIF
        i=i+1
      ENDDO
!Check for invalid name
      IF(idone.eq.0) THEN
        iflag=1
        message='SPECIES/ERROR:invalid species name'
      ENDIF
 1000 RETURN
      END
