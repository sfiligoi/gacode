      SUBROUTINE FLUXAV_EFIT(bmag,rmag,zmag,nx_xy,ny_xy,x_xy,y_xy,
     &                       psi_xy,n_lim,x_lim,y_lim,nr_r,psi_r,f_r,
     &                       ffp_r,q_r,rin_r,rout_r,rhot_r,elong_r,
     &                       gth_r,gph_r,vol_r,vp_r,r2_r,rm2_r,phit_r,
     &                       fm_r,grth_r,b2_r,bm2_r,grho1_r,grho2_r,
     &                       gr2bm2_r,ftrap_r,fhat_r,bpout_r,btout_r,
     &                       a0,rm1_r,rhor_r,iflag,message)
!***********************************************************************
!FLUXAV_EFIT calculates flux surface quantities from EFIT equilibria
!References:
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  Lin-Liu, Miller, Phys Plasmas 2 (1995) 1666
!  W.A.Houlberg 3/2000
!Input:
!  bmag-axial magnetic field (T)
!  rmag-major radius of magnetic axis (m)
!  zmag-vertical position of magnetic axis (m)
!  nx_xy-number of x points on psi(x,y) grid
!  ny_xy-number of y points on psi(x,y) grid
!  x_xy(i)-horizontal grid for 2-D poloidal flux (m)
!  y_xy(i)-vertical grid for 2-D poloidal flux (m)
!  psi_xy(i,j)-poloidal flux/(2*pi) on 2-D grid (Wb/rad)
!  n_lim-number of points on limiter
!  x_lim(i)-horizontal positions of limiter points (m)
!  y_lim(i)-vertical positions of limiter points (m)
!  nr_r-number of radial points
!  psi_r(i)-poloidal flux/2*pi (Wb/rad)
!  f_r(i)-R*B_t (m*T)
!  ffp_r(i)-f_r*df_r/dpsi_r (rad*T)
!  q_r(i)-safety factor
!Output:
!  rin_r(i)-major radius of intersection of surface on inside (m)
!  rout_r(i)-major radius of intersection of surface on outside (m)
!  rhot_r(i)-normalized toroidal flux grid proportional to (Phi)**0.5
!  rhor_r(i)-rho, sqrt(toroidal flux/pi*Bt)
!  elong_r(i)-plasma elongation
!  gth_r(i)-<gtt/sqrt(g)>-theta average (-)
!  gph_r(i)-<sqrt(g)/R**2>-theta average (-)
!  vol_r(i)-volume enclosed (m**3)
!  vp_r(i)-d vol_r/d rhot_r/a0 (m**2)
!  r2_r(i)-<R**2> (m**2)
!  rm2_r(i)-<1/R**2> (/m**2)
!  phit_r(i)-toroidal flux (Wb)
!  fm_r(3,i)-geometric factor
!  grth_r(i)-<n.grad(theta)> (1/m)
!  b2_r(i)-<B**2> (T**2)
!  bm2_r(i)-<1/B**2> (/T**2)
!  grho1_r(i)-a0*<|grad(rhot_r)|>
!  grho2_r(i)-a0**2*<|grad(rhot_r)|**2>
!  gr2bm2_r(i)-a0**2*<|grad(rhot_r)|**2/B**2> (1/T**2)
!  ftrap_r(i)-trapped particle fraction
!  fhat_r(i)-RB_t/dpsi_r/drho_t/a0
!  bpout_r(i)-poloidal field at rout_r(i) (T)
!  btout_r(i)-toroidal field at rout_r(i) (T)
!  a0-minor radius, half diameter of boundary flux surface (m)
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!Comments:
!  Modified from routines in ONETWO created by H. StJohn
!***********************************************************************
      IMPLICIT NONE
!Declaration of parameters
!  x,y grid dimensions for EFIT psi values
      INTEGER        mxnx_xy,                 mxny_xy
      PARAMETER     (mxnx_xy=130,             mxny_xy=130)
!Declaration of input variables
      INTEGER        nx_xy,                   ny_xy,
     &               n_lim,                   nr_r
      REAL           bmag,                    rmag,
     &               zmag
      REAL           x_xy(*),                 y_xy(*),
     &               psi_xy(mxnx_xy,*)
      REAL           x_lim(*),                y_lim(*) 
      REAL           psi_r(*),                f_r(*),
     &               ffp_r(*),                q_r(*)
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           a0
      REAL           b2_r(*),                 bm2_r(*),
     &               bpout_r(*),              btout_r(*),
     &               elong_r(*),              fhat_r(*),
     &               fm_r(3,*),               ftrap_r(*),
     &               gph_r(*),                gr2bm2_r(*),
     &               grho1_r(*),              grho2_r(*),
     &               grth_r(*),               gth_r(*),
     &               phit_r(*),               r2_r(*),
     &               rhot_r(*),               rin_r(*),
     &               rm2_r(*),                rout_r(*),
     &               vol_r(*),                vp_r(*),
     &               rm1_r(*),                rhor_r(*)
!Declaration of local variables
      LOGICAL        use_cnt1
      REAL           psifctr
!  Radial grid for EFIT profiles
      INTEGER        mxnr_r
      PARAMETER     (mxnr_r=300)
      REAL           bpl(mxnr_r),             psivlcpy(mxnr_r)
!  Spline coefficients
      REAL           cspln(2,mxnx_xy,2*mxny_xy),
     &               wnoperm(2*mxnx_xy*mxny_xy+2*mxny_xy),
     &               pds(6)
!  Number of contours
      INTEGER        mxncon
      PARAMETER     (mxncon=2000)
      REAL           arclen(mxncon),          bpinv(mxncon),
     &               bp(mxncon),              h(mxncon),
     &               hlin(mxncon),            ydum(mxncon),
     &               xp(mxncon),              yp(mxncon),
     &               xp0(mxncon),             yp0(mxncon),
     &               bp0(mxncon),             bt(mxncon),
     &               btot(mxncon),            ndgradb(mxncon),
     &               nthinv(mxncon),          ctheta(mxncon)
      INTEGER        i,                       iflip,
     &               ii,                      iold,
     &               j,                       k,
     &               kflarc,                  kflauto,
     &               kflauto_s,               kflconv,
     &               mp,                      mpmin
      REAL           a,                       arcl,
     &               avbp,                    bb,
     &               bmax,                    bperr,
     &               bpinteg,                 capbr,
     &               capbz,                   dang,
     &               dbdr,                    dbdz,
     &               dbrdr,                   dbrdz,
     &               dbtdr,                   dbtdz,
     &               dbzdr,                   dbzdz,
     &               dpsia,                   dpsidrho,
     &               fmc,                     fms,
     &               ftl,                     ftu,
     &               h2a,                     ha,
     &               hca,                     sum2,
     &               taxis,                   tlim,
     &               xmax,                    xmin,
     &               ydif,                    ymax,
     &               ymin,                    yref,
     &               zcmax,                   zcmin
      REAL           z_pi
!Initialization
!  Physical and conversion constants
      z_pi=ACOS(-1.0)
!Set options
!  Set first contour method as default (more robust)
      use_cnt1=.true.
!  Set outer boundary at surface
      psifctr=0.0
!  Set minimum number of acceptable points on a surface
      mpmin=50
!Get extremes of limiter positions
      xmin=x_lim(1)
      xmax=x_lim(1)
      ymin=y_lim(1)
      ymax=y_lim(1)
      DO i=2,n_lim
        IF(x_lim(i).gt.xmax) xmax=x_lim(i)
        IF(x_lim(i).lt.xmin) xmin=x_lim(i)
        IF(y_lim(i).gt.ymax) ymax=y_lim(i)
        IF(y_lim(i).lt.ymin) ymin=y_lim(i)
      ENDDO   
      x_lim(n_lim+1)=xmin
      x_lim(n_lim+2)=xmax
      y_lim(n_lim+1)=ymin
      y_lim(n_lim+2)=ymax
!Get bicubic spline coefficients for poloidal flux
      iflip=0
      IF(psi_r(2).gt.psi_r(1)) THEN
        iflag=-1
        message='FLUXAV_EFIT(1)/WARNING:flip sign of Psi'
        iflip=1
!       Flip signs on 2-D poloidal flux array
        DO i=1,nx_xy
          DO j=1,ny_xy
            psi_xy(i,j)=-psi_xy(i,j)
          ENDDO
        ENDDO
      ENDIF
!  Set up bicubic spline coefficient array
      CALL V_BICUBIC(psi_xy,x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,wnoperm,
     &               iflag,message)
      IF(iflag.ne.0) THEN
        i=LEN(message)
        message='FLUXAV_EFIT(2)/'//message(1:i-15)
        IF(iflag.eq.1) GOTO 1000
      ENDIF
      IF(iflip.eq.1) THEN
!       Flip signs back on 2-D poloidal flux array
        DO i=1,nx_xy
          DO j=1,ny_xy
            psi_xy(i,j)=-psi_xy(i,j)
          ENDDO
        ENDDO
      ENDIF
!All contours to be found must be inside the box defined by limiters
      xmin=x_lim(n_lim+1)
      xmax=x_lim(n_lim+2)
      ymin=y_lim(n_lim+1)
      ymax=y_lim(n_lim+2)
      arcl=0.02
      taxis=5.0
      tlim=30.0
!Make copy of flux levels
!The target value of the flux will be overwritten with the acutal value
   10 DO i=1,nr_r
        IF(iflip.eq.1) THEN
          psivlcpy(i)=-psi_r(i)
        ELSE
          psivlcpy(i)=psi_r(i)
        ENDIF
      ENDDO
!From the plasma edge to the contour next to the magnetic axis
!  - find the contour describing the surface (xp(i),yp(i)),i = 1...mp
!  - form the required flux surface averages and get other surface info
      DO j=nr_r,2,-1 
        kflarc=0
   20   IF(j.eq.nr_r) THEN
          kflauto=1
          kflconv=1
          psivlcpy(nr_r)=psifctr*(psivlcpy(nr_r-1)-psivlcpy(nr_r))
     &                   +psivlcpy(nr_r)
          a=(tlim-taxis)/(psivlcpy(nr_r)-psivlcpy(1))
          dpsia=psivlcpy(nr_r-1)-psivlcpy(nr_r)
        ELSE
          kflauto=0
          kflconv=0
        ENDIF
        dang=(a*(psivlcpy(j)-psivlcpy(1))+taxis)*(kflarc+1)
        bperr=0.01
        kflauto_s=kflauto
        IF(use_cnt1) THEN
!         Start with first contour method
          kflauto=kflauto_s
          CALL FLUXAV_EC1(kflconv,x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,
     &                    xmin,xmax,ymin,ymax,rmag,zmag,mxncon,arcl,
     &                    dpsia,kflauto,psivlcpy(j),bperr,xp,yp,mp,bp,
     &                    iflag,message)
          IF(iflag.eq.1) THEN
            i=LEN(message)
            message='FLUXAV_EFIT(3)/'//message(1:i-15)
            GOTO 1000
          ELSEIF(iflag.lt.0) THEN
!           First method failed, try to recover with second
            iflag=0
            kflauto=kflauto_s
            CALL FLUXAV_EC2(kflconv,x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,
     &                      xmin,xmax,ymin,ymax,rmag,zmag,mxncon,dang,
     &                      dpsia,kflauto,psivlcpy(j),bperr,xp,yp,mp,bp,
     &                      iflag,message)
            IF(iflag.eq.1) THEN
              i=LEN(message)
              message='FLUXAV_EFIT(4)/'//message(1:i-15)
              GOTO 1000
            ENDIF
          ENDIF
        ELSE
!         Start with second contour method
          kflauto=kflauto_s
          CALL FLUXAV_EC2(kflconv,x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,
     &                    xmin,xmax,ymin,ymax,rmag,zmag,mxncon,dang,
     &                    dpsia,kflauto,psivlcpy(j),bperr,xp,yp,mp,
     &                    bp,iflag,message)
          IF(iflag.eq.1) THEN
            i=LEN(message)
            message='FLUXAV_EFIT(5)/'//message(1:i-15)
            GOTO 1000
          ELSEIF(iflag.lt.0) THEN
!           Second method failed, try first
            iflag=0
            kflauto=kflauto_s
            CALL FLUXAV_EC1(kflconv,x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,
     &                      xmin,xmax,ymin,ymax,rmag,zmag,mxncon,arcl,
     &                      dpsia,kflauto,psivlcpy(j),bperr,xp,yp,mp,bp,
     &                      iflag,message)
          ENDIF
        ENDIF
        IF(iflag.eq.-1) THEN
!         Both contour methods fail
          iflag=1
          i=LEN(message)
          message='FLUXAV_EFIT(6)/'//message(1:i-15)
          GOTO 1000
        ELSEIF(iflag.eq.1) THEN
!         Nonrecoverable error in contour routines
          i=LEN(message)
          message='FLUXAV_EFIT(7)/'//message(1:i-15)
          GOTO 1000
        ELSEIF(iflag.lt.-1) THEN
!         Too many points; reduce accuracy requirement and try again
          iflag=0
          kflarc=kflarc+1
          arcl=arcl*1.3
          IF(kflarc.lt.5) GOTO 20
        ENDIF
        IF(j.eq.nr_r.and.mp.lt.mpmin) THEN
!         Too few points on boundary
          IF(psifctr.lt.1.0e-3) THEN
!           Pull in outer boundary a little and retry
            psifctr=0.1
            GOTO 20
          ELSE
            iflag=1
            i=LEN(message)
            message='FLUXAV_EFIT(8)/ERROR:boundary reduction failed'
            GOTO 1000
          ENDIF
        ENDIF
        IF(iflip.eq.0) THEN
!         Change sign of poloidal field
          CALL RARRAY_SCALE(mp,-1.0,bp,1)
        ENDIF
!Compute volume enclosed
        vol_r(j)=0.0
        DO i=2,mp
          vol_r(j)=vol_r(j)+(xp(i)-xp(i-1))
     &                     *(xp(i)*yp(i)+xp(i-1)*yp(i-1))/2.0
        ENDDO
        vol_r(j)=ABS(vol_r(j))*2.0*z_pi
!       See if volume of this surface is less than next outer surface
        IF(j.lt.nr_r.and.vol_r(j)-vol_r(j+1).ge.0.0) THEN
!         Volume is not correct, see if other contour method can be used
          IF(use_cnt1) THEN
!           Start all over with the other contour method as default
            use_cnt1=.false.
            GOTO 10
          ELSE
!           Both methods of finding the contours have failed
            IF(psifctr.lt.1.0e-3) THEN
!             Pull in outer boundary a little and start all over
              psifctr=0.1
              GOTO 10
            ELSE
!             Everything failed
              iflag=1
              message='FLUXAV_EFIT(9)/ERROR:total failure'
              GOTO 1000
            ENDIF
          ENDIF
        ENDIF
!Re-order contour arrays so that contour begins on the outside
!  Copy arrays
        DO i=1,mp
          xp0(i)=xp(i)
          yp0(i)=yp(i)
          bp0(i)=bp(i)
        ENDDO
!  Find index of point nearest the plane of the axis on outside
        ii=0
        yref=1.0e6
        DO i=1,mp
          ydif=ABS(zmag-yp0(i))
          IF(ydif.le.yref.and.xp0(i).ge.rmag) THEN
            ii=i
            yref=ydif
          ENDIF
        ENDDO
!  Reorder arrays from new starting point
        DO i=1,mp-1
          iold=i+ii-1
          IF(iold.gt.mp) iold=iold-mp+1
          xp(i)=xp0(iold)
          yp(i)=yp0(iold)
          bp(i)=bp0(iold)
        ENDDO
        xp(mp)=xp(1)
        yp(mp)=yp(1)
        bp(mp)=bp(1)
!Compute arclength around the contour
        arclen(1)=0.0
        DO i=1,mp-1
          arclen(i+1)=arclen(i)+SQRT((xp(i+1)-xp(i))**2
     &                              +(yp(i+1)-yp(i))**2)
        ENDDO
!Find index of point nearest the plane of the axis on inside
!Find maximum and minimum z on a contour
        ii=0
        yref=1.0e6
        zcmin=0.0
        zcmax=0.0
        DO i=1,mp
          ydif=ABS(yp(i)-zmag)
          IF(xp(i).le.rmag.and.ydif.le.yref) THEN
            ii=i
            yref=ydif
          ENDIF
          IF(yp(i).le.zcmin) zcmin=yp(i)
          IF(yp(i).ge.zcmax) zcmax=yp(i)
        ENDDO
!Inside and outside major radii 
        rout_r(j)=xp(1)
        rin_r(j)=xp(ii)
!Minor radius
        IF(j.eq.nr_r) a0=0.5*(rout_r(j)-rin_r(j))
!Elongation
        elong_r(j)=(zcmax-zcmin)/(rout_r(j)-rin_r(j))
!Outside poloidal and toroidal field
        bpout_r(j)=bp(1)
        btout_r(j)=f_r(j)/rout_r(j)
!Magnetic field
        DO i=1,mp
          bt(i)=f_r(j)/xp(i)
          bpinv(i)=1.0/bp(i)
!         Get psi and derivatives
          CALL W_BICUBIC(x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,xp(i),
     &                   yp(i),6,pds,iflag,message)
          IF(iflag.ne.0) THEN
            iflag=1
            k=LEN(message)
            message='FLUXAV_EFIT(10)/'//message(1:k-16)
            IF(iflag.eq.1) GOTO 1000
          ENDIF
          IF(iflip.eq.1) THEN
!           Real psi is opposite sign to what is returned
            pds(1)=-pds(1)
            pds(2)=-pds(2)
            pds(3)=-pds(3)
            pds(4)=-pds(4)
            pds(5)=-pds(5)
            pds(6)=-pds(6)
          ENDIF
!  B field components
          capbr=-pds(3)/xp(i)
          capbz=pds(2)/xp(i)
          btot(i)=SQRT(capbr*capbr+capbz*capbz+bt(i)*bt(i))
!  B field derivatives
          dbrdr=pds(3)/(xp(i)*xp(i))-pds(4)/xp(i)
          dbzdr=-pds(2)/(xp(i)*xp(i))+pds(5)/xp(i)
          dbtdr=-f_r(j)/(xp(i)*xp(i))+ffp_r(j)/f_r(j)*pds(2)/xp(i)
          dbdr=(capbr*dbrdr+capbz*dbzdr+bt(i)*dbtdr)/btot(i)
          dbrdz=-pds(6)/xp(i)
          dbzdz=pds(4)/xp(i)
          dbtdz=ffp_r(j)/f_r(j)*pds(3)/xp(i)
          dbdz=(capbr*dbrdz+capbz*dbzdz+bt(i)*dbtdz)/btot(i)
!  n dot grad(B)
          ndgradb(i)=-(capbr*dbdr+capbz*dbdz)/btot(i)
!  Use the identity Lp/2pi*B dot grad(theta) = bp
          nthinv(i)=ABS(btot(i))/(bp(i))
        ENDDO
!  Calculate n dot grad(CapTheta)=grth_r
        CALL W_LIN_INTEG(arclen,nthinv,mp,sum2)
        grth_r(j)=2.0*z_pi/sum2
!  Calculate CapTheta array
        ctheta(1)=0.0
        DO i=2,mp
          CALL W_LIN_INTEG(arclen,nthinv,i,sum2)
          ctheta(i)=sum2*grth_r(j)
        ENDDO
!Surface integrals
!  Integral of mag(bpoloidal)
        CALL W_LIN_INTEG(arclen,bp,mp,bpl(j))
!  Integral of 1.0/bpoloidal
        CALL W_LIN_INTEG(arclen,bpinv,mp,bpinteg)
!  dV/dpsi-convert with d(psi)/drho
        vp_r(j)=2.0*z_pi*bpinteg
!  <B**2>
        DO i=1,mp
          ydum(i)=bpinv(i)*btot(i)*btot(i)
        ENDDO
        CALL W_LIN_INTEG(arclen,ydum,mp,sum2)
        b2_r(j)=sum2/bpinteg
!  <1/B**2>
        DO i=1,mp
          ydum(i)=bpinv(i)/(btot(i)*btot(i))
        ENDDO
        CALL W_LIN_INTEG(arclen,ydum,mp,sum2)
        bm2_r(j)=sum2/bpinteg
!  <(grad psi/B)**2>
        DO i=1,mp
          ydum(i)=xp(i)*xp(i)*bp(i)/(btot(i)*btot(i))
        ENDDO
        CALL W_LIN_INTEG(arclen,ydum,mp,sum2)
        gr2bm2_r(j)=sum2/bpinteg
!  Poloidal expansion factors Fm(1...3)
!  <n dot grad(small theta)>
        DO i=1,mp
          ydum(i)= bpinv(i)*btot(i)/nthinv(i)
        ENDDO
        CALL W_LIN_INTEG(arclen,ydum,mp,sum2)
        avbp=2.0*z_pi*sum2/arclen(mp)/bpinteg
        DO k=1,3
          DO i=1,mp
            ydum(i)=bpinv(i)*SIN(k*ctheta(i))*ndgradb(i)
          ENDDO
          CALL W_LIN_INTEG(arclen,ydum,mp,sum2)
          fms=sum2/bpinteg
          DO i=1,mp
            ydum(i)=ydum(i)*grth_r(j)*ABS(btot(i))
          ENDDO
          CALL W_LIN_INTEG(arclen,ydum,mp,sum2)
          fms=fms*sum2/bpinteg
          DO i=1,mp
            ydum(i)=bpinv(i)*COS(k*ctheta(i))*ndgradb(i)
          ENDDO
          CALL W_LIN_INTEG(arclen,ydum,mp,sum2)
          fmc=sum2/bpinteg
          DO i=1,mp
            ydum(i)=ydum(i)*grth_r(j)*ABS(btot(i))
          ENDDO
          CALL W_LIN_INTEG(arclen,ydum,mp,sum2)
          fmc=fmc*sum2/bpinteg
          fm_r(k,j)=2.0/(b2_r(j)*avbp)*(fms+fmc)
        ENDDO 
!  rm2_r=<1/R**2>
        DO i=1,mp
          ydum(i)=bpinv(i)/(xp(i)**2)
        ENDDO
        CALL W_LIN_INTEG (arclen,ydum,mp,sum2)
        rm2_r(j)=sum2/bpinteg
!  r2_r=<R**2>:
        DO i=1,mp
          ydum(i)=bpinv(i)*xp(i)**2
        ENDDO
        CALL W_LIN_INTEG(arclen,ydum,mp,sum2)
        r2_r(j)=sum2/bpinteg
!  rm1_r=<1/R>
        DO i=1,mp
          ydum(i)=bpinv(i)/(xp(i))
        ENDDO
        CALL W_LIN_INTEG (arclen,ydum,mp,sum2)
        rm1_r(j)=sum2/bpinteg
!  <|(grad psi)|>-convert with d(psi)/drho
        CALL W_LIN_INTEG(arclen,xp,mp,sum2)
        grho1_r(j)=sum2/bpinteg
c        write(*,*) j, grho1_r(j), 'XX'        
!  <(grad psi)**2>-convert with d(psi)/drho
        DO i=1,mp
          ydum(i)=xp(i)*xp(i)*bp(i)
        ENDDO
        CALL W_LIN_INTEG(arclen,ydum,mp,sum2)
        grho2_r(j)=sum2/bpinteg
!  <(grad psi)**2/R**2>-convert with d(psi)/drho
        gph_r(j)=bpl(j)/bpinteg             
! Max total B field on contour
        bmax=-1.0e30
        DO i=1,mp
!         Convert to grad psi squared
          bp(i)=(xp(i)*bp(i))**2
          bb=SQRT(f_r(j)**2+bp(i))/xp(i)
          IF(bb.gt.bmax) bmax=bb
        ENDDO
        DO i=1,mp
          h(i)=ABS(bmag)/(SQRT(f_r(j)**2+bp(i))/xp(i))
          ydum(i)=h(i)**2*bpinv(i)
          hlin(i)=ABS(bmag/(bmax*h(i)))
          hlin(i)=MIN(hlin(i),1.0)   
        ENDDO
!  Trapped fraction, use approximation of Lin-Liu and Miller
        DO i=1,mp
          ydum(i)=(hlin(i)**2)*bpinv(i)
        ENDDO
!       <h**2>
        CALL W_LIN_INTEG(arclen,ydum,mp,h2a)
        h2a=h2a/bpinteg           
!       <h>   
        DO i=1,mp
          ydum(i)=hlin(i)*bpinv(i)
        ENDDO
        CALL W_LIN_INTEG(arclen,ydum,mp,ha)
        ha=ha/bpinteg             
        DO i=1,mp
          ydum(i)=(1.0-(SQRT(1.0-hlin(i))*(1.0+0.5*hlin(i))))
     &            *bpinv(i)/(hlin(i)**2)                     
        ENDDO
        CALL W_LIN_INTEG(arclen,ydum,mp,hca)
        hca=hca/bpinteg
        ftl=1.0-h2a*hca
        ftu=1.0-(1.0-(1.0+0.5*ha)*SQRT(1.0-ha))
     &      *h2a/(ha*ha)
        ftrap_r(j)=0.75*ftu+0.25*ftl
      ENDDO
!Set axial values
      elong_r(1)=elong_r(2)
      rin_r(1)=rmag
      rout_r(1)=rmag
      bpout_r(1)=0.0
      btout_r(1)=bmag
      rm2_r(1)=1.0/rmag**2
      rm1_r(1)=1.0/rmag
      r2_r(1)=1.0/rm2_r(1)
      vol_r(1)=0.0
      bpl(1)=0.0
      b2_r(1)=bmag*bmag
      bm2_r(1)=1.0/b2_r(1)
      gr2bm2_r(1)=0.0
      grth_r(1)=grth_r(2)
      ftrap_r(1)=0.0
      DO k=1,3
        fm_r(k,1)=0.0
      ENDDO
      vp_r(1)=2.0*vol_r(2)/(psivlcpy(1)-psivlcpy(2))-vp_r(2)
!Set edge values
      rm2_r(nr_r)=rm2_r(nr_r-1)
      rm1_r(nr_r)=rm1_r(nr_r-1)
      r2_r(nr_r)=r2_r(nr_r-1)
!Toroidal flux at interior nodes, integrate dV/dPhi=2*pi/F/<R^-2>
      phit_r(1)=0.0
      DO j=2,nr_r
        phit_r(j)=phit_r(j-1)+(f_r(j-1)*rm2_r(j-1)+f_r(j)*rm2_r(j))
     &                        /(4.0*z_pi)*(vol_r(j)-vol_r(j-1))
      ENDDO
!Toroidal flux grid
      DO j=1,nr_r
        rhot_r(j)=SQRT(phit_r(j)/phit_r(nr_r))
        rhor_r(j)=rhot_r(j)*SQRT(phit_r(nr_r)/(z_pi*bmag))
      ENDDO
!Convert from grad(Psi) to grad(rho)
      a0=rhor_r(nr_r)
      DO j=2,nr_r-1
! <ABS (grad rho)>
        dpsidrho=ABS(rhot_r(j)*phit_r(nr_r)/(z_pi*q_r(j)*a0))
        grho1_r(j)=grho1_r(j)/dpsidrho
! <(grad rho)**2>
        grho2_r(j)=grho2_r(j)/dpsidrho**2
!       write(*,*) j, rhot_r(j), rhor_r(j), grho1_r(j)
! dV/drho
        vp_r(j)=vp_r(j)*dpsidrho
! Gph
        gph_r(j)=gph_r(j)/dpsidrho**2*vp_r(j)/(2.0*z_pi)**2
! Gth=<1/R**2>V'/(2*pi)**2
        gth_r(j)=rm2_r(j)*vp_r(j)/(2.0*z_pi)**2
      ENDDO
      grho1_r(1)=grho1_r(2)             
      grho1_r(nr_r)=grho1_r(nr_r-1)             
      grho2_r(1)=grho2_r(2)
      grho2_r(nr_r)=grho2_r(nr_r-1)
      vp_r(1)=0.0
      vp_r(nr_r)=vp_r(nr_r-1)
      gph_r(1)=0.0
      gph_r(nr_r)=gph_r(nr_r-1)*rhot_r(nr_r)/rhot_r(nr_r-1)
      gth_r(1)=0.0
      gth_r(nr_r)=gth_r(nr_r-1)*rhot_r(nr_r)/rhot_r(nr_r-1)
!Calculate fhat
      DO j=2,nr_r-1
        fhat_r(j)=f_r(j)/(psi_r(j+1)-psi_r(j-1))*a0
     &            *(rhot_r(j+1)-rhot_r(j-1))
      ENDDO
      fhat_r(1)=fhat_r(2)
      fhat_r(nr_r)=fhat_r(nr_r-1)
 1000 RETURN
      END
      SUBROUTINE FLUXAV_EC1(kflconv,x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,
     &                      xmin,xmax,ymin,ymax,rmag,zmag,mxncon,arcl,
     &                      dpsia,kflauto,psivl,bperr,xp,yp,mp,bp,iflag,
     &                      message)
!***********************************************************************
!FLUXAV_EC1 generates a set of points lying on a flux surface (contour)
!  W.A. Houlberg, H. St. John 3/2000
!Input:
!  kflconv-convergence option when kflauto=1 (-)
!         =0 not used
!         =1 contour is converged to the outermost contour that never
!            leaves the search box given by xmin,xmax,ymin,ymax
!  x_xy(nx_xy)-array of x values for MHD grid (m)
!  nx_xy-number of x_xy grid points (-)
!  y_xy(ny_xy)-array of y values for MHD grid (m)
!  ny_xy-number of y_xy grid points (-)
!  cspln(2,mxnx_xy,*)-bicubic spline coefficient array
!  mxnx_xy-dimension of x_xy grid in spline coefficients (-)
!  xmin-minimum allowable x value for search (m)
!  xmax-maximum allowable x value for search (m)
!  ymin-minimum allowable y value for search (m)
!  ymax-maximum allowable y value for search (m)
!  rmag-x position of interior (axis) point (m)
!  zmag-y position of interior (axis) point (m)
!  mxncon-max number of points to be returned in vectors xp,yp (-)
!  arcl-arc length of step (m)
!  dpsia-psi grid spacing at edge (Wb/radian)
!       =psi(penultimate value)-psi(edge value)
!Input/Output:
!  kflauto-option for error recovery (-)
!         =0 no recovery
!         =1 automatic error recovery
!           (may be set to 0 in if a failure occurs)
!  psivl-psi value for which contour is desired (Wb/radian)
!  bperr-relative change in poloidal b field between (xp(i),yp(i)) and
!        (xp(i+1),yp(i+1))
!Output:
!  xp(mp)-array of x values on contour (m)
!  yp(mp)-array of y values on contour (m)
!  mp-number of points on surface (-)
!  bp(mp)-poloidal bfield at xp,yp (T)
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!Comments:
!  Not intended to be a general purpose contouring routine
!  Designed specifically for finding plasma flux surfaces assuming that
!    - The plasma surface is convex (bean shapes may be problematic)
!    - The function to be contoured is monotonically decreasing in the
!      region of interest
!  Uses a bicubic spline representation of psi to get psi values at
!    arbitrary points
!  The spline coefficients must exist (in cspln) before entry
!  The contour must fully encircle (rmag,zmag)
!  Given (rmag,zmag) and a psi value, psivl, generate a contour of
!    ordered points,(xp(i),yp(i)), i = 1,mp, which has (rmag,zmag) as
!    as an interior point
!  The search is limited to a rectangle bounded by (xmin,xmax,ymin,ymax)
!  Define the search box (xmin,xmax,ymin,ymax) by a rectangle, within
!    which the plasma must reside under all circumstances:
!    - If the search box is set to the entire MHD grid, this will
!      normally not work satisfactorily because there are local minima
!      and maxima around the f coils which confuse the contour tracing
!    - Use the rectangle defined by the extremes of the limiters
!  If the number of elements (xp,yp) exceeds the limit (mxncon) bperr is
!    relaxed (up to twice) by the increment dbperr before an error exit
!    - bperr=0.03 is suggested
!  There are three ways to run:
!    a) kflauto=0, kflconv not used
!       Should be ok for all interior plasma surfaces
!       The value of psi to be traced, psivl, is not altered
!       Either the routine returns with the desired contour or a failure
!         is indicated with iflag when the surface passes outside box
!    b) kflauto=1, kflconv=0
!       Can be used on the plasma boundary near a separatrix
!       The value of psi to be traced, psivl, may be altered
!       If the contour with the value of psivl could not be traced or
!         passes outside the box the value of psivl will be repeatedly
!         brought closer to the value on the magnetic axis, either by
!         some fraction of dpsia if it is set or some fraction of
!         psiaxis-psilim if it is not
!       The first contour succesfully traced will be returned
!    c) kflauto=1, kflconv=1
!       Same as b) except that a binary search is done to find the
!         plasma surface
!       For this purpose the plasma surface is defined as the largest
!         closed flux surface that can be found inside the search box
!       The actual search box used is not critical for diverted plasmas
!         but is very important for limited plasmas
!  If iflag<0:
!    - The error condition on bperr has been relaxed before exit
!    - To continue this routine could be called again with a larger arcl
!    - This has not to be a problem if mxncon is large enough to begin
!      with to give a reasonable representation of the contour
!  If the value of psi to be traced leads to an error and kflauto=1 then
!    dpsia is used to control the amount by which the value of psi is
!    adjusted, otherwise dpsia is not used
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        kflconv,                 mxncon,
     &               mxnx_xy,                 nx_xy,
     &               ny_xy
      REAL           arcl,                    dpsia,
     &               rmag,                    xmax,
     &               xmin,                    ymax,
     &               ymin,                    zmag
      REAL           cspln(2,mxnx_xy,*),      x_xy(*),
     &               y_xy(ny_xy)
!Declaration of input/output variables
      INTEGER        kflauto
      REAL           bperr,                   psivl
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag,                   mp
      REAL           bp(*),                   xp(*),
     &               yp(*)
!Declaration of local variables
      LOGICAL        reset_psi
      INTEGER        i,                       iflage,
     &               iflg,                    iqpt,
     &               isgn,                    isgnsave,
     &               itry,                    newti,
     &               nitr8,                   nsteps
      REAL           a,                       arclmin,
     &               b,                       bp1,
     &               bp2,                     bperrsave,
     &               bpmintol,                cmult,
     &               cost,                    dbperr,
     &               derrt,                   dth,
     &               dtharcl,                 dthmin,
     &               dx,                      dxmin,
     &               dxx,                     dy,
     &               dymin,                   dyy,
     &               psi1,                    psi2,
     &               psiin,                   psimag,
     &               psiout,                  psisave,
     &               serr,                    serrt,
     &               sint,                    step,
     &               theta,                   thend,
     &               thnew,                   thstart,
     &               x1,                      x2,
     &               xmult,                   xn,
     &               xns,                     xsave,
     &               y1,                      y2,
     &               ymult,                   yn,
     &               yns,                     ysave
      REAL           pds(6)
      DATA nitr8/    10/
      REAL           z_pi
!Initialization
!  Physical and conversion constants
      z_pi=ACOS(-1.0)
!  Number of tries to pull in flux surface to avoid limiter
      itry=0
      psiin=0.0
      bperrsave=bperr
      dbperr=bperr/2.0
      reset_psi=.false.
!  If min bp drops below bpmintol and kflauto=1, the contour is searched
!    for an x point
      bpmintol=0.05
!  Set minimum arc length for determining minimum angle
      arclmin=0.10*arcl
!  Get psi at (rmag,zmag)
      CALL W_BICUBIC(x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,rmag,zmag,1,
     &               pds,iflag,message)
      IF(iflag.eq.-1) THEN
        i=LEN(message)
        message='FLUXAV_EC1(1)/'//message(1:i-14)
      ELSEIF(iflag.eq.1) THEN
!       Magnetic axis is not in the domain of (x_xy,y,xy)
        message='FLUXAV_EC1(2)/ERROR:mag axis not in domain'
        GOTO 1000
      ENDIF
      psimag=pds(1)
      IF(psimag.lt.psivl) THEN
        iflag=1
        message='FLUXAV_EC1(3)/ERROR:psimag less that psivl'
        GOTO 1000
      ENDIF
   10 mp=0
      thstart=1.25*z_pi
      thend=thstart+2.0*z_pi
      theta=thstart
!  Minimum increment in theta between rays
!  May be increased later, after the first point is found
      dthmin=1.0e-5
      dth=0.0
!  Set step to twice grid spacing
      dxx=x_xy(3)-x_xy(1)
      dyy=y_xy(3)-y_xy(1)
!  Set minimum acceptable step
      dxmin=0.05*(x_xy(2)-x_xy(1))
      dymin=0.05*(y_xy(2)-y_xy(1))
      dx=dxx
      dy=dyy
!  Set absolute error convergence criteria for Newton's method
      serrt=3.5e-06
      derrt=0.5e-07
!Loop over theta from thstart to thend (=thstart+2*pi)
   20 theta=theta+dth
      IF(theta.gt.2.0*z_pi.and.thstart.ne.0.0) THEN
        theta=theta-2.0*z_pi
        thend=thend-2.0*z_pi
      ENDIF
      theta=MIN(theta,thend)
      IF(theta.lt.thend) THEN
        dxx=dx
        dyy=dy
        iqpt=0
!  Get equation of ray emanating from (rmag,zmag)
        CALL FLUXAV_ELN(rmag,zmag,theta,a,b,iflg,isgn)
!  Now have y=a*x+b (iflg=0)
!  or       x=a*y+b (iflg=1)
!  Start search from axis point
        x1=rmag
        y1=zmag
        psi1=psimag
        cost=COS(theta)
        sint=SIN(theta)
!  Sliding interval search, max width of interval ~1.41*(dx or dy)
        nsteps=0
        cmult=1.0
        xsave=-1.0e20
   30   IF(iflg.eq.0) THEN
!         Search in x
          xmult=1.0
          y2=ymin-dymin
          DO WHILE(y2.lt.ymin.or.y2.gt.ymax)
!           Increment x
            x2=x1+isgn*dxx*xmult*cmult
!           Limit the search
            x2=MAX(x2,xmin)
            x2=MIN(x2,xmax)
            IF(x2.eq.x1) THEN
!             If x values are the same, psivl is not within reach
              reset_psi=.true.
              GOTO 1000
            ENDIF
!           Get corresponding value of y
            y2=a*x2+b
            xmult=0.5*xmult
          ENDDO
        ELSE
!         Search in y
          ymult=1.0
          x2=xmin-dxmin
          DO WHILE(x2.lt.xmin.or.x2.gt.xmax)
!           Increment y
            y2=y1+isgn*dyy*ymult*cmult
!           Limit the search
            y2=MAX(y2,ymin)
            y2=MIN(y2,ymax)
!           If y values are the same, psivl is not within reach
            IF(y2.eq.y1) THEN
              reset_psi=.true.
              GOTO 1000
            ENDIF
!           Get corresponding value of x
            x2=a*y2+b
            ymult=ymult*0.5
          ENDDO
        ENDIF
        CALL W_BICUBIC(x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,x2,y2,1,pds,
     &                 iflag,message)
        IF(iflag.eq.-1) THEN
          i=LEN(message)
          message='FLUXAV_EC1(4)/'//message(1:i-14)
        ELSEIF(iflag.eq.1) THEN
!         (x2,y2) not in domain
          message='FLUXAV_EC1(5)/ERROR:(x2,y2) is not in domain'
          GOTO 1000
        ENDIF
        psi2=pds(1)
        IF((psivl-psi1)*(psivl-psi2).gt.0.0) THEN
!         Since psi is increasing, we are in the vicinity of an x-point
!         To be sure that we stay on the contour which envelopes the
!         plasma, we first search for the minimum in psi along the ray
!         We guarantee that the point we find is on the contour by not
!         allowing the ray to extend past this minimum
!         Use Newton's method to search for the minimum -- the point
!         along the ray where the directional derivative of psi along
!         the ray is zero
          IF(cmult.lt.0.99) THEN
            cmult=cmult*0.5
            GOTO 30
          ENDIF
          IF(psi2-psi1.gt.0.0) THEN    ! psi is increasing
            xsave=x1   ! save current guess for possible later restore
            ysave=y1
            psisave=psi1
            isgnsave=isgn
            xn=x1
            yn=y1
            newti=0
            step=serrt+2
            DO WHILE(ABS(step).gt.serrt.and.newti.le.nitr8)
              CALL W_BICUBIC(x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,xn,yn,
     &                       6,pds,iflag,message)
              IF(iflag.eq.-1) THEN
                i=LEN(message)
                message='FLUXAV_EC1(6)/'//message(1:i-14)
              ELSEIF(iflag.eq.1) THEN
!               (xn,yn) not in domain        
                message='FLUXAV_EC1(7)/ERROR:(xn,yn) not in domain'
                GOTO 1000
              ENDIF
              step=-(pds(2)*cost+pds(3)*sint)/(cost*cost*pds(5)
     &             +2.0*cost*sint*pds(4)+sint*sint*pds(6))
              xn=xn+step*cost
              yn=yn+step*sint
              newti=newti+1
            ENDDO
            IF(newti.le.nitr8) THEN
!             Newton's method converged
              isgn=-isgn
              IF(iflg.eq.1) dyy=0.5*dyy
              IF(iflg.eq.0) dxx=0.5*dxx
              x2=xn
              y2=yn
              psi2=pds(1)
!             The minimum in psi along the ray must be less than or
!             equal to the value for which we are searching, psivl
!             If this is not the case then we can't find psivl on the
!             ray so abandon the search
              IF(psi2.gt.psivl) THEN
                reset_psi=.true.
                GOTO 1000
              ENDIF
            ELSE
!             Set cmult<1
              cmult=0.5
              x1=xsave
              y1=ysave
              psi1=psisave
!             Restore search direction
              isgn=isgnsave
              GOTO 30
            ENDIF
          ENDIF
          nsteps=nsteps+1
          x1=x2
          y1=y2
          psi1=psi2
          GOTO 30
        ENDIF
!Now have psivl between psi1 and psi2, converge using Newton-Raphson
        newti=0
        IF(iflg.eq.0) THEN
          xn=x1+isgn*dxx*0.5
          yn=a*xn+b
        ELSE
          yn=y1+isgn*dyy*0.5
          xn=a*yn+b
        ENDIF
   40   iflag=0
        CALL W_BICUBIC(x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,xn,yn,3,pds,
     &                 iflag,message)
        IF(iflag.gt.0) THEN
!         (xn,yn) is not in the domain of (x_xy,y,xy)
!         Set number of Newton iterations to maximum to restart
          newti=nitr8
          iflag=0
        ELSE
          serr=-(pds(1)-psivl)/(pds(2)*cost+pds(3)*sint)
          IF(ABS(serr).lt.serrt) GOTO 50
          IF(psivl.eq.0.0) THEN
            IF(ABS(pds(1)-psivl).lt.derrt) GOTO 50
          ELSE
            IF(ABS((pds(1)-psivl)/psivl).lt.derrt) GOTO 50
          ENDIF
          xn=xn+serr*cost
          yn=yn+serr*sint
          newti=newti+1
        ENDIF
        IF(newti.ge.nitr8) THEN
          IF(dxx.le.dxmin) GOTO 50
          IF(dyy.le.dymin) GOTO 50
          IF(iflg.eq.0) dxx=0.5*dxx
          IF(iflg.eq.1) dyy=0.5*dyy
          IF(iqpt.eq.0) GOTO 30
          theta=theta-dth
          GOTO 20
        ENDIF
        GOTO 40
!End of Newton iteration
!Check for sufficient accuracy in point spacing as determined by theta
!Accuracy test is based on a relative error in poloidal b field of bperr
   50   bp2=SQRT(pds(2)**2+pds(3)**2)/xn
        IF(theta.ne.thstart) THEN
!         Not the first point of the contour
          IF(ABS(bp2-bp1)/MAX(bp2,bp1).ge.bperr) THEN
!           Relative change in bp too large
            IF(dth.ne.dthmin) THEN
!             Spacing too large for grad psi, decrease theta and retry
              dth=dth*0.5
!             dth lower limit must be observed
              dth=MAX(dth,dthmin)
            ENDIF
          ENDIF
        ENDIF
!Acceptable point found, collect it and set up for next point
        mp=mp+1
        IF(mp.eq.mxncon) THEN
!         Ran out of storage for the contour points
!         If mxncon is reasonable then bperr may be too small
          IF(bperr-bperrsave.ge.2.0*dbperr) THEN
!           Increasing bperr has been tried and failed
            iflag=-2
            message='FLUXAV_EC1(8)/WARNING:bperr increased and failed'
            GOTO 1000
          ELSE
!           Increase bperr by dbperr and retry to generate the contour
            bperr=bperr+dbperr
            GOTO 10
          ENDIF
        ENDIF
!Found new point
!  Save it
        xp(mp)=xn
        yp(mp)=yn
        bp(mp)=bp2
!  Set up for next ray
        IF(mp.eq.1) THEN
!         Estimate the radius and relate the theta increment to arc length
          dthmin=arclmin/SQRT((xn-rmag)**2+(yn-zmag)**2)
          dtharcl=arcl/arclmin*dthmin
        ENDIF
        IF(ABS(bp2-bp1)/MAX(bp2,bp1).le.0.25*bperr) dth=2.0*dth
!  Set bp1 for error test on next ray
        bp1=bp2
!  Don't let dth get out of range
        dth=MIN(dth,dtharcl)
        dth=MAX(dth,dthmin)
!  For contours sufficiently far removed from rmag,zmag, use the
!    following approach:
!    - An approximation to the new point,(xns,yns) is found by moving
!      along the tangent line (after 2nd step), for a distance arcl
!    - If near an x point, as signaled by bp < bpmintol, then switch
!      to the ray method with small dx,dy increments to avoid crossing
!      the x point
        IF(nsteps.gt.2.and.bp(mp).gt.bpmintol) THEN
!         Use tangent line method
          iflage=0
          CALL FLUXAV_ENP(pds,arcl,xn,yn,theta,xmin,xmax,ymin,ymax,zmag,
     &                    rmag,sint,cost,xns,yns,thnew,iflage)
!         If this failed, go back and use the ray method
          IF(iflage.gt.0) GOTO 20
          dth=thnew-theta
!         Do not use thnew here
          theta=theta+dth
          IF(theta.gt.2.0*z_pi.and.thstart.ne.0.0) THEN
            theta=theta-2.0*z_pi
            thend=thend-2.0*z_pi
          ENDIF
          theta=MIN(theta,thend)
          IF(theta.lt.thend) THEN
            sint=SIN(theta)
            cost=COS(theta)
            xn=xns
            yn=yns
            newti=0
            iqpt=1
!           Skip directly to Newton's method
            GOTO 40
          ENDIF
        ELSE
!         Use ray method
          GOTO 20
        ENDIF
      ENDIF
!Normal exit
      mp=mp+1
      xp(mp)=xp(1)
      yp(mp)=yp(1)
      bp(mp)=bp(1)
      bperr=bperrsave
      IF(itry.ne.0.and.kflconv.eq.1) THEN
        psiin=psivl
        IF(ABS(0.5*(psiout+psiin)-psivl).lt.1.0e-6) GOTO 1000
        psivl=0.5*(psiout+psiin)
        GOTO 10
      ENDIF
!Check for open contour with possible recovery
 1000 IF(reset_psi) THEN
        IF(kflauto.eq.0) THEN
!         No recovery
          iflag=1
          message='FLUXAV_EC1(9)/ERROR:no recovery'
        ELSE
          psiout=psivl
!         Try to correct improper value of psivl
          IF(dpsia.gt.0.0) THEN
            psivl=psivl+dpsia/FLOAT(99)
          ELSE
            psivl=psivl+0.0005*(psimag-psivl)
          ENDIF
          IF((itry.gt.0).and.(psiin.ne.0.0).and.(kflconv.eq.1)) THEN
            psivl=0.5*(psiin+psiout)
          ENDIF
          itry=itry+1
          IF(itry.le.100) THEN
!           Go back and try again
            reset_psi=.false.
            GOTO 10
          ELSE
!           Allow error exit to be taken
            iflag=1
          message='FLUXAV_EC1(10)/ERROR:no recovery'
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE FLUXAV_EC2(kflconv,x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,
     &                      xmin,xmax,ymin,ymax,rmag,zmag,mxncon,dang,
     &                      dpsia,kflauto,psivl,bperr,xp,yp,mp,bp,iflag,
     &                      message)
!***********************************************************************
!FLUXAV_EC2 generates a set of points lying on a flux surface (contour)
!  W.A. Houlberg, H. St. John 3/2000
!Input:
!  kflconv-convergence option when kflauto=1 (-)
!         =0 not used
!         =1 contour is converged to the outermost contour that never
!            leaves the search box given by xmin,xmax,ymin,ymax
!  x_xy(nx_xy)-array of x values for MHD grid (m)
!  nx_xy-number of x_xy grid points (-)
!  y_xy(ny_xy)-array of y values for MHD grid (m)
!  ny_xy-number of y_xy grid points (-)
!  cspln(2,mxnx_xy,*)-bicubic spline coefficient array
!  mxnx_xy-dimension of x_xy grid in spline coefficients (-)
!  xmin-minimum allowable x value for search (m)
!  xmax-maximum allowable x value for search (m)
!  ymin-minimum allowable y value for search (m)
!  ymax-maximum allowable y value for search (m)
!  rmag-x position of interior (axis) point (m)
!  zmag-y position of interior (axis) point (m)
!  mxncon-max number of points to be returned in vectors xp,yp (-)
!  dang-angular step for the rays emmanating from (rmag,zmag) (degrees)
!  dpsia-psi grid spacing at edge (Wb/radian)
!       =psi(penultimate value)-psi(edge value)
!Input/Output:
!  kflauto-option for error recovery (-)
!         =0 no recovery
!         =1 automatic error recovery
!           (may be set to 0 in if a failure occurs)
!  psivl-psi value for which contour is desired (Wb/radian)
!  bperr-relative change in poloidal b field between (xp(i),yp(i)) and
!        (xp(i+1),yp(i+1))
!Output:
!  xp(mp)-array of x values on contour (m)
!  yp(mp)-array of y values on contour (m)
!  mp-number of points on surface (-)
!  bp(mp)-poloidal bfield at xp,yp (T)
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!Comments:
!  Not intended to be a general purpose contouring routine
!  Designed specifically for finding plasma flux surfaces assuming that
!    - The plasma surface is convex (bean shapes may be problematic)
!    - The function to be contoured is monotonically decreasing in the
!      region of interest
!  Uses a bicubic spline representation of psi to get psi values at
!    arbitrary points
!  The spline coefficients must exist (in cspln) before entry
!  The contour must fully encircle (rmag,zmag)
!  Given (rmag,zmag) and a psi value, psivl, generate a contour of
!    ordered points,(xp(i),yp(i)), i = 1,mp, which has (rmag,zmag) as
!    as an interior point
!  The search is limited to a rectangle bounded by (xmin,xmax,ymin,ymax)
!  Define the search box (xmin,xmax,ymin,ymax) by a rectangle, within
!    which the plasma must reside under all circumstances:
!    - If the search box is set to the entire MHD grid, this will
!      normally not work satisfactorily because there are local minima
!      and maxima around the f coils which confuse the contour tracing
!    - Use the rectangle defined by the extremes of the limiters
!  If the number of elements (xp,yp) exceeds the limit (mxncon) bperr is
!    relaxed (up to twice) by the increment dbperr before an error exit
!    - bperr=0.03 is suggested
!  There are three ways to run:
!    a) kflauto=0, kflconv not used
!       Should be ok for all interior plasma surfaces
!       The value of psi to be traced, psivl, is not altered
!       Either the routine returns with the desired contour or a failure
!         is indicated with iflag when the surface passes outside box
!    b) kflauto=1, kflconv=0
!       Can be used on the plasma boundary near a separatrix
!       The value of psi to be traced, psivl, may be altered
!       If the contour with the value of psivl could not be traced or
!         passes outside the box the value of psivl will be repeatedly
!         brought closer to the value on the magnetic axis, either by
!         some fraction of dpsia if it is set or some fraction of
!         psiaxis-psilim if it is not
!       The first contour succesfully traced will be returned
!    c) kflauto=1, kflconv=1
!       Same as b) except that a binary search is done to find the
!         plasma surface
!       For this purpose the plasma surface is defined as the largest
!         closed flux surface that can be found inside the search box
!       The actual search box used is not critical for diverted plasmas
!         but is very important for limited plasmas
!  If iflag<0:
!    - The error condition on bperr has been relaxed before exit
!  If the value of psi to be traced leads to an error and kflauto=1 then
!    dpsia is used to control the amount by which the value of psi is
!    adjusted, otherwise dpsia is not used
!  dang should be set according to shape of equilibrium
!    - If dang is too small, many unnecessary points will be generated
!    - If dang is too large for bperr, the routine wastes cpu time
!      reducing the magnitude of dang
!    - dang=10deg near psilim and dang=30deg near psimax is suggested
!    - For highly elongated plasmas control with dang is difficult to
!      set for all contours.
!    - If dang yields bperr(computed)<bperr(input), then dang is used,
!      otherwise the angle is successively reduced by 0.5 until this
!      condition is met
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        kflconv,                 mxncon,
     &               mxnx_xy,                 nx_xy,
     &               ny_xy
      REAL           dang,                    dpsia,
     &               rmag,                    xmax,
     &               xmin,                    ymax,
     &               ymin,                    zmag
      REAL           cspln(2,mxnx_xy,*),      x_xy(*),
     &               y_xy(ny_xy)
!Declaration of input/output variables
      INTEGER        kflauto
      REAL           bperr,                   psivl
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag,                   mp
      REAL           bp(*),                   xp(*),
     &               yp(*)
!Declaration of local variables
      LOGICAL        reset_psi
      INTEGER        i,                       iflg,
     &               isgn,                    itry,
     &               newti,                   nitr8
      REAL           a,                       b,
     &               bp1,                     bp2,
     &               bperrsave,               cost,
     &               dbperr,                  derrt,
     &               dth,                     dth0,
     &               dthmin,                  dx,
     &               dxmin,                   dxx,
     &               dy,                      dymin,
     &               dyy,                     psi1,
     &               psi2,                    psiin,
     &               psimag,                  psiout,
     &               serr,                    serrt,
     &               sint,                    theta,
     &               thend,                   thstart,
     &               x1,                      x2,
     &               xmult,                   xn,
     &               y1,                      y2,
     &               ymult,                   yn
      REAL           pds(6)
      DATA nitr8/    10/
      REAL           z_pi
!Initialization
!  Physical and conversion constants
      z_pi=ACOS(-1.0)
!  Number of tries to pull in flux surface to avoid limiter
      itry=0
      psiin=0.0
      bperrsave=bperr
      dbperr=bperr/2.0
      reset_psi=.false.
!  Get psi at (rmag,zmag)
      CALL W_BICUBIC(x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,rmag,zmag,1,
     &               pds,iflag,message)
      IF(iflag.eq.-1) THEN
        i=LEN(message)
        message='FLUXAV_EC2(1)/'//message(1:i-14)
      ELSEIF(iflag.eq.1) THEN
!       Magnetic axis is not in the domain of (x_xy,y,xy)
        message='FLUXAV_EC2(2)/ERROR:mag axis not in domain'
        GOTO 1000
      ENDIF
      psimag=pds(1)
      IF(psimag.lt.psivl) THEN
        iflag=1
        message='FLUXAV_EC2(3)/ERROR:psimag less that psivl'
        GOTO 1000
      ENDIF
   10 mp=0
      thstart=1.25*z_pi
      thend=thstart+2.0*z_pi
      theta=thstart
!  Minimum increment in theta between rays
      dth0=2.0*z_pi*dang/360.0
      dthmin=1.0e-5
      IF(dthmin.gt.dth0) dthmin=0.5*dth0
      dth=0.0
!  Set step to twice grid spacing
      dxx=x_xy(3)-x_xy(1)
      dyy=y_xy(3)-y_xy(1)
!  Set minimum acceptable step
      dxmin=0.05*(x_xy(2)-x_xy(1))
      dymin=0.05*(y_xy(2)-y_xy(1))
      dx=dxx
      dy=dyy
!  Set absolute error convergence criteria for Newton's method
      serrt=3.5e-06
      derrt=0.5e-07
!Loop over theta from thstart to thend (=thstart+2*pi)
   20 theta=theta+dth
      IF(theta.gt.2.0*z_pi.and.thstart.ne.0.0) THEN
        theta=theta-2.0*z_pi
        thend=thend-2.0*z_pi
      ENDIF
      theta=MIN(theta,thend)
      IF(theta.lt.thend) THEN
        dxx=dx
        dyy=dy
!  Get equation of ray emanating from (rmag,zmag)
        CALL FLUXAV_ELN(rmag,zmag,theta,a,b,iflg,isgn)
!  Now have y=a*x+b (iflg=0)
!  or       x=a*y+b (iflg=1)
!  Start search from axis point
        x1=rmag
        y1=zmag
        psi1=psimag
        cost=COS(theta)
        sint=SIN(theta)
!  Sliding interval search, max width of interval ~1.41*(dx or dy)
   30   IF(iflg.eq.0) THEN
!         Search in x
          xmult=1.0
          y2=ymin-dymin
          DO WHILE(y2.lt.ymin.or.y2.gt.ymax)
!           Increment x
            x2=x1+isgn*dxx*xmult
!           Limit the search
            x2=MAX(x2,xmin)
            x2=MIN(x2,xmax)
            IF(x2.eq.x1) THEN
!             If x values are the same, psivl is not within reach
              reset_psi=.true.
              GOTO 1000
            ENDIF
!           Get corresponding value of y
            y2=a*x2+b
            xmult=0.5*xmult
          ENDDO
        ELSE
!         Search in y
          ymult=1.0
          x2=xmin-dxmin
          DO WHILE(x2.lt.xmin.or.x2.gt.xmax)
!           Increment y
            y2=y1+isgn*dyy*ymult
!           Limit the search
            y2=MAX(y2,ymin)
            y2=MIN(y2,ymax)
!           If y values are the same, psivl is not within reach
            IF(y2.eq.y1) THEN
              reset_psi=.true.
              GOTO 1000
            ENDIF
!           Get corresponding value of x
            x2=a*y2+b
            ymult=ymult*0.5
          ENDDO
        ENDIF
        CALL W_BICUBIC(x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,x2,y2,1,pds,
     &                 iflag,message)
        IF(iflag.eq.-1) THEN
          i=LEN(message)
          message='FLUXAV_EC2(4)/'//message(1:i-14)
        ELSEIF(iflag.eq.1) THEN
!         (x2,y2) not in domain
          message='FLUXAV_EC2(5)/ERROR:(x2,y2) is not in domain'
          GOTO 1000
        ENDIF
        psi2=pds(1)
        IF((psivl-psi1)*(psivl-psi2).gt.0.0) THEN
          x1=x2
          y1=y2
          psi1=psi2
          GOTO 30
        ENDIF
!Now have psivl between psi1 and psi2, converge using Newton-Raphson
        newti=0
        IF(iflg.eq.0) THEN
          xn=x1+isgn*dxx*0.5
          yn=a*xn+b
        ELSE
          yn=y1+isgn*dyy*0.5
          xn=a*yn+b
        ENDIF
   40   iflag=0
        CALL W_BICUBIC(x_xy,nx_xy,y_xy,ny_xy,cspln,mxnx_xy,xn,yn,3,pds,
     &                 iflag,message)
        IF(iflag.gt.0) THEN
!         (xn,yn) is not in the domain of (x_xy,y,xy)
!         Set number of Newton iterations to maximum to restart
          newti=nitr8
          iflag=0
        ELSE
          serr=-(pds(1)-psivl)/(pds(2)*cost+pds(3)*sint)
          IF(ABS(serr).lt.serrt) GOTO 50
          IF(psivl.eq.0.0) THEN
            IF(ABS(pds(1)-psivl).lt.derrt) GOTO 50
          ELSE
            IF(ABS((pds(1)-psivl)/psivl).lt.derrt) GOTO 50
          ENDIF
          xn=xn+serr*cost
          yn=yn+serr*sint
          newti=newti+1
        ENDIF
        IF(newti.ge.nitr8) THEN
          IF(dxx.le.dxmin) THEN
            iflag=1
            message='FLUXAV_EC2(6)/ERROR:dxx lesss than dxmin'
            GOTO 1000
          ENDIF
          dxx=0.5*dxx
          dyy=0.5*dyy
          GOTO 30
        ENDIF
        GOTO 40
!End of Newton iteration
!Check for sufficient accuracy in point spacing as determined by theta
!Accuracy test is based on a relative error in poloidal b field of bperr
   50   bp2=SQRT(pds(2)**2+pds(3)**2)/xn
        IF(theta.ne.thstart) THEN
!         Not the first point of the contour
          IF(ABS(bp2-bp1)/MAX(bp2,bp1).ge.bperr) THEN
!           Relative change in bp too large
            IF(dth.ne.dthmin) THEN
!             Spacing too large for grad psi, decrease theta and retry
              theta=theta-dth
              dth=dth*0.5
!             dth lower limit must be observed
              dth=MAX(dth,dthmin)
              GOTO 20
            ENDIF
          ENDIF
        ENDIF
        bp1=bp2
        mp=mp+1
        IF(mp.eq.mxncon) THEN
!         Ran out of storage for the contour points
!         If mxncon is reasonable then bperr may be too small
          IF(bperr-bperrsave.ge.2.0*dbperr) THEN
!           Increasing bperr has been tried and failed
            iflag=-2
            message='FLUXAV_EC2(7)/WARNING:bperr increased and failed'
            GOTO 1000
          ELSE
!           Increase bperr by dbperr and retry to generate the contour
            bperr=bperr+dbperr
            GOTO 10
          ENDIF
        ENDIF
!Found new point
!  Save it
        xp(mp)=xn
        yp(mp)=yn
        bp(mp)=bp2
        IF(ABS(bp2-bp1)/MAX(bp2,bp1).le.0.25*bperr) dth=2.0*dth
!  Don't let dth get out of range
        dth=MIN(dth,dth0)
        dth=MAX(dth,dthmin)
!       Start next ray
        GOTO 20
      ENDIF 
!Normal exit
      mp=mp+1
      xp(mp)=xp(1)
      yp(mp)=yp(1)
      bp(mp)=bp(1)
      bperr=bperrsave
      IF(itry.ne.0.and.kflconv.eq.1) THEN
        psiin=psivl
        IF(ABS(0.5*(psiout+psiin)-psivl).lt.1.0e-6) GOTO 1000
        psivl=0.5*(psiout+psiin)
        GOTO 10
      ENDIF
!Check for open contour with possible recovery
 1000 IF(reset_psi) THEN
        IF(kflauto.eq.0) THEN
!         No recovery
          iflag=1
          message='FLUXAV_EC2(8)/ERROR:no recovery'
        ELSE
          psiout=psivl
!         Try to correct improper value of psivl
          IF(dpsia.gt.0.0) THEN
            psivl=psivl+dpsia/FLOAT(99)
          ELSE
            psivl=psivl+0.0005*(psimag-psivl)
          ENDIF
          IF((itry.gt.0).and.(psiin.ne.0.0).and.(kflconv.eq.1)) THEN
            psivl=0.5*(psiin+psiout)
          ENDIF
          itry=itry+1
          IF(itry.le.100) THEN
!           Go back and try again
            reset_psi=.false.
            GOTO 10
          ELSE
!           Allow error exit to be taken
            iflag=1
            message='FLUXAV_EC2(9)/ERROR:no recovery'
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE FLUXAV_ELN(xaxis,yaxis,theta,a,b,keqxy,keqdir)
!***********************************************************************
!FLUXAV_ELN finds the equation of a straight line passing through
!  (xaxis,yaxis) and with slope whose angle is theta
!  W.A. Houlberg, H. St. John 3/98
!Input:
!  xaxis-starting x coordinate, nominally the magnetic axis (m)
!  yaxis-starting y coordinate, nominally the magnetic axis (m)
!  theta-poloidal angle defining line direction (radians)
!Output:
!  a-slope of line (-)
!  b-intercept of line (m)
!  keqxy-switch indicating form of equation for line
!       =0 y=a*x+b
!       =1 x=a*y+b
!  keqdir-switch indicating direction to increment abscissa to move
!         along the line from (xaxis,yaxis) toward the boundary
!        =+1 x or y increases from axis to surface
!        =-1 x or y decreases from axis to surface
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      REAL           theta,                   xaxis,
     &               yaxis
!Declaration of output variables
      INTEGER        keqdir,                  keqxy
      REAL           a,                       b
!Declaration of local variables
      REAL           theta1
      REAL           z_pi
!Initialization
!  Physical and conversion constants
      z_pi=ACOS(-1.0)
!Get equation of ray emanating from (xaxis,yaxis)
      IF(((0.25*z_pi.le.theta).and.(theta.le.0.75*z_pi)).or.
     &   ((1.25*z_pi.le.theta).and.(theta.le.1.75*z_pi))) THEN
!       Express x as a function of y, x=a*y+b
        keqxy=1
        IF(theta.gt.z_pi) THEN
!         Bottom of plasma
          keqdir=-1
          theta1=1.5*z_pi-theta
          IF(theta.gt.1.5*z_pi) theta1=z_pi-ABS(theta1)
        ELSE
!         Top of plasma
          keqdir=1
          theta1=0.5*z_pi-theta
          IF(theta.gt.0.5*z_pi) theta1=2.0*z_pi-ABS(theta1)
        ENDIF
        a=TAN(theta1)
        b=xaxis-a*yaxis
      ELSE
!       Express y as a function of x, y=a*x+b
        keqxy=0
        IF((theta.lt.0.25*z_pi).or.(theta.gt.1.75*z_pi)) THEN
!         Outside of plasma
          keqdir=1
        ELSE
!         Inside of plasma
          keqdir=-1
        ENDIF
        a=TAN(theta)
        b=yaxis-a*xaxis
      ENDIF
      RETURN
      END
      SUBROUTINE FLUXAV_ENP(pds,arcl,xn,yn,theta,xmin,xmax,ymin,ymax,
     &                      yaxd,xaxd,sint,cost,xns,yns,thnew,iflag)
!***********************************************************************
!FLUXAV_ENP makes an approximation to get the next point on a flux
!  surface contour given the current point and the gradient in psi
!  W.A. Houlberg, H. St. John 3/98
!Input:
!pds(6)-bicubic spline coefficients for polidal flux
!arcl-radius of circle at xn,yn (m)
!xn-x value of reference point (m)
!yn-y value of reference point (m)
!theta-reference poloidal angle (radians)
!xmin-minimum x value for new point (m)
!xmax-maximum x value for new point (m)
!ymin-minimum y value for new point (m)
!ymax-maximum y value for new point (m)
!yaxd-y value for axis (m)
!xaxd-x value for axis (m)
!sint-sin(theta) (-)
!cost-cos(theta) (-)
!Output:
!xns-x value of new point (m)
!yns-y value of new point (m)
!thnew-vlaue of theta at new point (radians)
!iflag-error flag (-)
!     =0 no errors
!     =1 error in finding new point
!Comments:
!  The tangent line at the point (xn,yn) is df = 0 = df/dx*dx+df/dy*dy
!  A circle of radius arcl, centered at (xn,yn) is arcl**2 = dr**2+dz**2
!  This gives us two equations in two unkowns: dr and dz .
!  The new point is xn +/- dx and yn +/- dy where the signs have to
!    be picked so that theta increases.
!  Note the special treatment required as thnew crosses zero (outside
!    this routine)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      REAL           pds(*)
      REAL           arcl,                    cost,
     &               sint,                    theta,
     &               xaxd,                    xmax,
     &               xmin,                    xn,
     &               yaxd,                    ymax,
     &               ymin,                    yn
!Declaration of output variables
      INTEGER        iflag
      REAL           thnew,                xns,
     &               yns
!Declaration of local variables
      REAL           alpha,                   dpx,
     &               dpy,                     dx,
     &               dy
      REAL           z_pi
!Initialization
!  Error flag
      iflag=0
!  Physical and conversion constants
      z_pi=ACOS(-1.0)
      dpx=pds(2)
      dpy=pds(3)
      IF(ABS(dpx).gt.ABS(dpy)) THEN
        IF(dpx.eq.0.0) THEN
          iflag=1
          GOTO 1000
        ENDIF
        alpha=dpy/dpx
        dy=arcl/SQRT(1.0+alpha*alpha)
        dx=-alpha*dy
      ELSE
        IF(dpy.eq.0.0) THEN
          iflag=1
          GOTO 1000
        ENDIF
        alpha=dpx/dpy
        dx=arcl/SQRT(1.0+alpha*alpha)
        dy=-alpha*dx
      ENDIF
!The sign on dx,dy must be taken so that theta increases
!A unit vector in the direction of increasing theta (i.e., theta
!  counterclockwise) is (-SIN (theta),COS (theta))
!The displacement vector is (dx,dy)
!Its projection on the above vector must be positive and equals
!  -dx*SIN (theta)+dy*COS (theta)
      IF(-dx*sint+dy*cost.lt.0.0) THEN
        dx=-dx
        dy=-dy
      ENDIF
      xns=xn+dx
      yns=yn+dy
      IF(xns.lt.xmin.or.xns.gt.xmax.or.yns.lt.ymin.or.yns.gt.ymax) THEN
        iflag=1
        GOTO 1000
      ENDIF
      thnew=ATAN2(yns-yaxd,xns-xaxd)
      IF(thnew.lt.0.0) THEN
        thnew=thnew+2.0*z_pi
      ELSEIF(theta.gt.7.0*z_pi/4.0) THEN
        thnew=thnew+2.0*z_pi
      ENDIF
 1000 RETURN
      END
