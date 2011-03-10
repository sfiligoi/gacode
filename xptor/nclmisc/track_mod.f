      SUBROUTINE TRACK_BCON(rho,theta,zeta,bar_iota,phiprm,sqrtg,
     &                      b_theta,b_zeta,iflag,message)
!***********************************************************************
!TRACK_BCON evaluates the contravariant components of B and the
!  derivatives of the stream function for specified flux coordinates
!References:
!  W.A.Houlberg, P.I. Strand 5/2000
!Input:
!  rho-radial flux coordinate (-)
!  theta-poloidal angle (radians)
!  zeta-toroidal angle (radians)
!  bar_iota-iota-bar (-)
!  phiprm-d(Phi)/d(rho) (Wb)
!  sqrtg-Jacobian of transformation from cyl to flux coord (m**3)
!Output:
!  b_theta-contravariant poloidal component of B (T/m)
!  b_zeta-contravzraint toroidal component of B (T/m)
!  iflag-error and warning flag (-)
!       =-1 warning
!       =0 none
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
      INCLUDE '../inc/pamxkm.inc'
      INCLUDE '../inc/pamxnx.inc'
      INCLUDE '../inc/comm30.inc'
!Declaration of input variables
      REAL           bar_iota,                phiprm,
     &               rho,                     sqrtg,
     &               theta,                   zeta                    
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           b_theta,                 b_zeta
!Declaration of local variables
      INTEGER        i,                       k
      INTEGER        k_vopt(3)
      REAL           value(3)
      REAL           cosk,                    el,
     &               el_theta,                el_zeta,
     &               x
      REAL           z_pi,                    z_precision
      DATA i/        1/
      SAVE           i
!Initialization
!  Physical and conversion constants
      z_pi=ACOS(-1.0)
      z_precision=2.0e-7
!  Other
      k_vopt(1)=1
      k_vopt(2)=0
      k_vopt(3)=0
      IF(rho.gt.1.0+10.0*z_precision) THEN
!       This point is outside separatrix
        iflag=-1
        message='TRACK_BCON/WARNING:outside MHD solution domain'
        GOTO 1000
      ENDIF      
!This point is inside the MHD solution domain
      el_theta=0.0
      el_zeta=0.0
      x=AMAX1(rho,1.0e2*z_precision)
      x=AMIN1(x,1.0-10.0*z_precision)
      DO k=1,nk3d
        cosk=COS(m3d(k)*theta-n3d(k)*zeta)
        CALL W_SPLINE(k_vopt,nx3d,x,xm3d,elmna(1,1,k),i,value)
        el=value(1)*x**ABS(m3d(k))
        el_theta=el_theta+el*m3d(k)*cosk
        el_zeta=el_zeta-el*n3d(k)*cosk
      ENDDO
      b_theta=phiprm/(2.0*z_pi*sqrtg)*(bar_iota-el_zeta)
      b_zeta=phiprm/(2.0*z_pi*sqrtg)*(1.0+el_theta)  
 1000 RETURN
      END
c      SUBROUTINE TRACK_BCYL(r,r_theta,r_zeta,z_theta,z_zeta,b_theta,
c     &                      b_zeta,b_r,b_z,b_phi,iflag,message)
!***********************************************************************
!TRACK_BCYL evaluates the cylindrical components of B
!References:
!  W.A.Houlberg, P.I. Strand 5/2000
!Input:
!  rho-radial flux coordinate (-)
!  theta-poloidal angle (radians)
!  phi-geometric toroidal angle (radians)
!  xiota-iota-bar (-)
!  phiprm-d(Phi)/d(rho) (Wb)
!  sqrtg-Jacobian of transformation from cyl to flux coord (m**3)
!Output:
!  el_phi-d(lambda)/d(phi)
!  el_theta-d(lambda)/d(theta)
!  b_theta-contravariant poloidal component of B (T/m)
!  b_phi-contravzraint toroidal component of B (T/m)
!  iflag-error and warning flag (-)
!       =-1 warning
!       =0 none
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      SUBROUTINE TRABFL(kbinv,r,z,phi,tolinv,br,bz,btor,rho,theta,
     &                  rx,rt,rp,zx,zt,zp,xiota,elp,elt,iflag,message)
!***********************************************************************
!TRABFL computes the cylindrical components of the magnetic field at a
!  specified position in either a) cylindrical coordinates and also
!  returns the corresponding flux coordinates, or b) flux coordinates
!References:
!  Hirshman, ORNL Theory Section Memo 85/19 (1985)
!  W.A.Houlberg, P.I. Strand 3/2000
!Input:
!  kbinv-option to do inversion (-)
!       =0 no inversion
!       =else do inversion
!  r-major radius (m)
!  z-vertical coordinate (m)
!  phi-geometric toroidal angle (radians)
!  tolinv-tolerance for inverse transformation (m)
!        =1.0e-5*a0 is suggested
!Output:
!  br-b dot grad(r) (T)
!  bz-b dot grad(z) (T)
!  btor-toroidal field (T)
!  rho-radial flux coordinate (-)
!  theta-poloidal angle (radians)
!  rx-dr/dx (m)
!  rt-dr/dtheta (m/radian)
!  rp-dr/dphi (m/radian)
!  zx-dz/dx (m)
!  zt-dz/dtheta (m/radian)
!  zp-dz/dphi (m/radian)
!  xiota-iota-bar (-)
!  elp-d(lambda)/d(phi)
!  elt-d(lambda)/d(theta)
!  iflag-error and warning flag (-)
!       =-1 warning
!       =0 none
!       =1 error
!  message-warning or error message (character)
!Comments:
!  The toroidal field is given by r*b dot grad(phi)
!***********************************************************************
      IMPLICIT NONE
      INCLUDE '../inc/pamxkm.inc'
      INCLUDE '../inc/pamxnx.inc'
      INCLUDE '../inc/comm30.inc'
!Declaration of input variables
      INTEGER        kbinv
      REAL           phi,                     r,
     &               tolinv,                  z
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           br,                      btor,
     &               bz,                      elp,
     &               elt,                     rho,
     &               rp,                      rt,
     &               rx,                      theta,
     &               xiota,                   zp,
     &               zt,                      zx
!Declaration of local variables
      INTEGER        i,                       k,
     &               nmes
      INTEGER        k_vopt(3)
      REAL           value(3)
      REAL           cosk,                    emn0,
     &               phiprm,                  sqrtg,
     &               tau,                     x,
     &               xh
      REAL           z_pi,                    z_precision
      DATA i/        1/
      SAVE           i
      SAVE           phiprm,                  sqrtg
!Initialization
!  Physical and conversion constants
      z_pi=ACOS(-1.0)
      z_precision=2.0e-7
!  Other
      br=0.0
      bz=0.0
      btor=0.0
      k_vopt(1)=1
      k_vopt(2)=0
      k_vopt(3)=0
      IF(kbinv.ne.0) THEN
!       Get flux coordinates rho and theta
        rx=0.0
        rt=0.0
        rp=0.0
        zx=0.0
        zt=0.0
        zp=0.0
        CALL TRAFLX(r,z,phi,tolinv,rho,theta,rx,rt,rp,zx,zt,zp,
     &              iflag,message)
        IF(iflag.ne.0) THEN
          nmes=LEN(message)
          message='TRABFL/'//message(1:nmes-7)
          GOTO 1000
        ENDIF
      ENDIF
!This point is inside solution domain
      elt=0.0
      elp=0.0
      xh=AMAX1(rho,1.0e2*z_precision)
      x=AMIN1(xh,1.0-10.0*z_precision)
      DO k=1,nk3d
        cosk=COS(m3d(k)*theta-n3d(k)*phi)
        CALL W_SPLINE(k_vopt,nx3d,x,xm3d,elmna(1,1,k),i,value)
        emn0=value(1)*x**ABS(m3d(k))
        elt=elt+emn0*m3d(k)*cosk
        elp=elp-emn0*n3d(k)*cosk
      ENDDO   
!Compute magnetic field components
      tau=zx*rt-rx*zt
      sqrtg=r*tau
      phiprm=2.0*rho*phitot
      CALL W_SPLINE(k_vopt,nx3d,x,xm3d,fiota,i,value)
      xiota=value(1)
      br=phiprm*((xiota-elp)*rt+(1.0+elt)*rp)/(2.0*z_pi*sqrtg)
      bz=phiprm*((xiota-elp)*zt+(1.0+elt)*zp)/(2.0*z_pi*sqrtg)
      btor=r*phiprm*(1.0+elt)/(2.0*z_pi*sqrtg)
 1000 RETURN
      END
      SUBROUTINE TRACK(a0,rhoi,nrhoi,rv,zv,av,mv,ir,iz,rint,zint,
     &                 tint,rhoint,thint,sl,nsurf,iflag,message)
!***********************************************************************
!TRACK finds the intersections (rint,zint,tint) of a segmented path,
!  defined by the set of points (rv,zv,av), with the flux surfaces rhoi
!  where (r,z,a) are cylindrical coordinates (R,phi,Z) aligned with R0
!References:
!  Attenberger, Houlberg, Hirshman, J Comp Phys 72 (1987) 435
!  W.A.Houlberg 3/2000
!Input:
!  a0-characteristic minor dimension of plasma (m)
!  rhoi(i)-rho value of i'th flux surface (-)
!  nrhoi-number of flux surfaces (-)
!  rv(j)-major radius of j'th point on the ray (m)
!  zv(j)-vertical coordinate of j'th point on the ray (m)
!  av(j)-toroidal angle of j'th point on the ray (radians)
!  mv-number of points on the ray (-)
!Output:
!  ir(k)-radial node numbers of the intersections (-)
!       =0 if k'th point is a node of the ray
!  iz(k)-cell number for segment k, between intersection k and k+1 (-)
!       =index of current volume element, if k'th point is a node of the ray
!  rint(k),tint(k),zint(k)-ordered coordinates for the intersection
!                        with the flux surfaces rhoi(i)
!                       -first point at rv(1),zv(1),av(1)
!                       -last point at rv(mv),zv(mv),av(mv)
!  rhoint(k)-rho at k'th point (-)
!  thint(k)-theta at k'th point (radians)
!  sl(k)-length along path to intersection k (m)
!  nsurf-number of nodes + intersections (-)
!  iflag-error and warning flag (-)
!       =-1 warning
!       =0 none
!       =1 error
!  message-warning or error message (character)
!Comments:
!  If maxsrf is reached before the end of the ray a return occurs
!  It is assumed that rhoi(i) is increasing with i so that the i'th
!    volume element is bounded by rhoi(i) and rhoi(i+1)
!  The region outside rhoi(nrhoi) is volume element nrhoi
!***********************************************************************
      IMPLICIT NONE
      INCLUDE '../inc/comtr7.inc'
      INCLUDE '../inc/comtr8.inc'
      INCLUDE '../inc/comtr9.inc'
!Declaration of input variables
      INTEGER        mv,                      nrhoi
      REAL           a0
      REAL           av(*),                   rhoi(*),
     &               rv(*),                   zv(*)
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag,                   nsurf
      INTEGER        ir(*),                   iz(*)
      REAL           rhoint(*),               rint(*),
     &               sl(*),                   thint(*),
     &               tint(*),                 zint(*)
!Declaration of local variables
      LOGICAL        l_start
      INTEGER        i,                       i1,
     &               i2,                      ihit,
     &               ilo,                     ilold,
     &               incr,                    j,
     &               jlo,                     k,
     &               l,                       mhalf,
     &               ns,                      nmes
      REAL           di,                      dint,
     &               dist,                    dnode,
     &               drds0,                   drdss,
     &               drdstr,                  drho,
     &               ds,                      dshort,
     &               dtry,                    dx,
     &               dy,                      dz,
     &               p0,                      ps,
     &               rho0,                    rhold,
     &               rhomax,                  rhop,
     &               rhor,                    rhos,
     &               rhoz,                    r0,
     &               rp0,                     rps,
     &               rs,                      rt0,
     &               rts,                     rx0,
     &               rxs,                     step,
     &               stepmn,                  stepmx,
     &               t0,                      tau0,
     &               taus,                    thold,
     &               tol,                     ts,
     &               xs,                      ys,
     &               zp0,                     zps,
     &               zs,                      zt0,
     &               zts,                     zx0,
     &               zxs
      REAL           z_precision
      SAVE           rhold,                   thold
!Declaration of external functions
      REAL           TRAQRH,                  TRAQRS,
     &               U_ZEROIN
      EXTERNAL       TRAQRH,                  TRAQRS
      DATA mhalf/    30/
      DATA rhold/    1.0/
!Initialization
!  Physical and conversion constants
      z_precision=2.0e-7
!  Other
      kwtra=0
      maxsrf=500
      mhalf=30
      maxstp=1000
      rhomax=100.0
      toltra=1.0e-3
!tol-tolerance for U_ZEROIN calculation of flux surface location along chord-m
      tol=toltra*a0
      tolinv=1.0e-2*tol
      stepmx=a0
      stepmn=2.0*tol
!Set starting point and number of intersections
      nrho=nrhoi
      x0=rv(1)*COS(av(1))
      y0=rv(1)*SIN(av(1))
      z0=zv(1)
      r0=rv(1)
      p0=av(1)
      rho0=rhold
      t0=thold
!Get slope in case starting point is far from plasma with bad Jacobian
      dx=rv(2)*COS(av(2))-x0
      dy=rv(2)*SIN(av(2))-y0
      dz=zv(2)-z0
      dist=SQRT(dx**2+dy**2+dz**2)
      gx=dx/dist
      gy=dy/dist
      gz=dz/dist
!  Step along chord until valid starting point is found
c      l_start=.false.
c      i=1
c      DO WHILE((.not.l_start).and.(i.le.40))
c        rho0=1.05
c        iflag=0
c        message=' '
c        CALL TRAFLX(r0,z0,p0,tolinv,rho0,t0,rx0,rt0,rp0,zx0,zt0,zp0,
c     &              iflag,message)
!    Check flag
c        IF((iflag.ne.1).and.(rho0.le.1.1)) THEN
c          l_start=.true.
!         Ignore warnings of solution outside separatrix
c          iflag=0
c          message=' '
c        ELSE
c          i=i+1
c          x0=x0+0.01*dist*gx
c          y0=y0+0.01*dist*gy
c          z0=z0+0.01*dist*gz
c          r0=SQRT(x0**2+y0**2)
c          p0=SIGN(1.0,y0)*ACOS(x0/r0)
c        ENDIF
c      ENDDO
c Reset start of chord
c      rv(1)=r0
c      zv(1)=z0
c      av(1)=p0
        CALL TRAFLX(r0,z0,p0,tolinv,rho0,t0,rx0,rt0,rp0,zx0,zt0,zp0,
     &              iflag,message)
      IF(iflag.eq.1.or.rho0.gt.1.1) THEN
        nmes=LEN(message)
        message='TRACK(1)/'//message(1:nmes-9)
        GOTO 1000
      ENDIF
      rhold=rho0
      thold=t0
      rhos=rho0
      ts=t0
      nsurf=1
      rint(1)=r0
      zint(1)=z0
      tint(1)=p0
      rhoint(1)=rho0
      thint(1)=t0
      ir(1)=0
      sl(1)=0.0
!Search for initial volume element
      iflag=1
      ilo=nrhoi+1
      DO WHILE((iflag.eq.1).and.(ilo.gt.1))
        ilo=ilo-1
        IF(rho0.ge.rhoi(ilo)) iflag=0
      ENDDO
      IF(iflag.eq.1) THEN
        message='TRACK(2)/ERROR:initial volume not found'
        GOTO 1000
      ENDIF
      dnode=0.0
!Loop over segments of path
      DO l=2,mv  ! 60
!  Set parameters for this path segment
        d0=0.0
        dx=rv(l)*COS(av(l))-rv(l-1)*COS(av(l-1))
        dy=rv(l)*SIN(av(l))-rv(l-1)*SIN(av(l-1))
        dz=zv(l)-zv(l-1)
        dist=SQRT(dx**2+dy**2+dz**2)
        gx=dx/dist
        gy=dy/dist
        gz=dz/dist
!  Step along one segment
        DO k=1,maxstp  ! 50
!         Find slope at start of step
          tau0=zx0*rt0-rx0*zt0
          rhor=-zt0/tau0
          rhoz=rt0/tau0
          rhop=(rp0*zt0-rt0*zp0)/tau0
          drds0=(rhor*x0-rhop*y0/r0)*gx/r0+(rhor*y0+rhop*x0/r0)*gy/r0
     &          +rhoz*gz
!         Set step size
          step=stepmx/nrho
          IF(step.gt.stepmx) step=stepmx
          IF(step.lt.stepmn) step=stepmn
          ds=d0+step
!         Do not step past next node
          IF(ds.gt.dist) ds=dist
          step=ds-d0
!  Attempt to take a step
!  Possible halving of step up to mhalf times
          DO j=1,mhalf
            rhos=rho0
            ts=t0
            CALL TRAXYZ(ds,rhos,ts,ps,rs,xs,ys,zs,rxs,rts,rps,zxs,
     &                  zts,zps,iflag,message)
!  Check flag
            IF(iflag.eq.-1) THEN
!             Ignore warnings of solution outside separatrix
              iflag=0
              message=' '
            ELSEIF(iflag.eq.1) THEN
              nmes=LEN(message)
              message='TRACK(3)/'//message(1:nmes-9)
              GOTO 1000
            ENDIF
            taus=zxs*rts-rxs*zts
            rhor=-zts/taus
            rhoz=rts/taus
            rhop=(rps*zts-rts*zps)/taus
            drdss=(rhor*xs-rhop*ys/rs)*gx/rs+(rhor*ys+rhop*xs/rs)*gy/rs
     &            +rhoz*gz
!  Check whether tangency point has been passed
            IF(SIGN(1.0,drds0).ne.SIGN(1.0,drdss)) THEN
!             Shorten step to land near d(rho)/ds = 0
              iflx=0
              rhflx0=rho0
              thflx0=t0
              rhflxs=rhos
              thflxs=ts
              dshort=U_ZEROIN(d0,ds,TRAQRS,tol,iflag,message)
!  Check flag
              IF(iflag.eq.-1) THEN
!               Ignore warnings of solution outside separatrix
                iflag=0
                message=' '
              ELSEIF(iflag.eq.1) THEN
                nmes=LEN(message)
                message='TRACK(4)/'//message(1:nmes-9)
                GOTO 1000
              ENDIF
              ds=dshort+tol+2.0*tolinv
              step=ds-d0
              rhos=rho0
              ts=t0
              CALL TRAXYZ(ds,rhos,ts,ps,rs,xs,ys,zs,rxs,rts,rps,
     &                    zxs,zts,zps,iflag,message)
!  Check flag
              IF(iflag.eq.-1) THEN
!               Ignore warnings of solution outside separatrix
                iflag=0
                message=' '
              ELSEIF(iflag.eq.1) THEN
                nmes=LEN(message)
                message='TRACK(5)/'//message(1:nmes-9)
                GOTO 1000
              ENDIF
            ENDIF
!  Check domain
            IF(rhos.gt.rhomax) THEN
!             rho is out of bounds - shorten step and retry
              step=step*0.5
              ds=d0+step
              GOTO 30
            ENDIF
!  Test whether an intersection has been passed
            IF(rhoi(ilo).lt.rhos.and.
     &      (rhos.lt.rhoi(ilo+1).or.ilo.eq.nrhoi)) THEN
!             Still in the same volume element
              d0=ds
              rho0=rhos
              t0=ts
              p0=ps
              r0=rs
              x0=xs
              y0=ys
              z0=zs
              rx0=rxs
              rt0=rts
              rp0=rps
              zx0=zxs
              zt0=zts
              zp0=zps
              IF((dist-ds).le.a0*10.0*z_precision) THEN
!               Hit a node
                nsurf=nsurf+1
                rint(nsurf)=r0
                zint(nsurf)=z0
                tint(nsurf)=p0
                rhoint(nsurf)=rho0
                thint(nsurf)=t0
                ir(nsurf)=0
                iz(nsurf-1)=ilo
                iz(nsurf)=ilo
                dnode=dnode+dist
                sl(nsurf)=dnode
                GOTO 60
              ENDIF
!  Update drds0
              tau0=zx0*rt0-rx0*zt0
              rhor=-zt0/tau0
              rhoz=rt0/tau0
              rhop=(rp0*zt0-rt0*zp0)/tau0
              drds0=rhoz*gz +
     &        (rhor*x0-rhop*y0/r0)*gx/r0+(rhor*y0+rhop*x0/r0)*gy/r0
!  Step completed
              GOTO 50
            ELSEIF(ilo.gt.1.and.rhoi(ilo-1).lt.rhos
     &        .and.rhos.lt.rhoi(ilo)) THEN
!             Crossed rhoi(ilo)
              ilold=ilo
              ihit=ilo
              ilo=ilo-1
              GOTO 40
            ELSEIF(rhoi(ilo+1).lt.rhos.and.
     &        (ilo+1.eq.nrhoi .or.
     &        (ilo+1.lt.nrhoi.and.rhos.lt.rhoi(ilo+2)))) THEN
!             Crossed rhoi(ilo+1)
              ilold=ilo
              ihit=ilo+1
              ilo=ilo+1
              GOTO 40
            ELSEIF(rhoi(ilo).eq.rhos.or.rhoi(ilo+1).eq.rhos) THEN
!             Stepped onto boundary of present volume element
!             Lengthen step to check for intersection
!             A tangent is treated as a near miss
              step=step+2.0*tolinv
              ds=d0+step
              GOTO 30
            ELSEIF(step.gt.stepmn) THEN
!             Crossed several flux surfaces - shorten step and retry
              step=step*0.5
              ds=d0+step
              GOTO 30
            ELSE
!             Crossed several flux surfaces in one small step
!             Locate present volume element, then interpolate from last
!             step to estimate closely spaced surfaces
              DO i=1,nrhoi
                jlo=nrhoi-(i-1)
                IF(rhos.ge.rhoi(jlo)) GOTO 20
              ENDDO   
   20         IF(jlo.lt.ilo) THEN
                i1=ilo
                i2=jlo+1
                ns=ilo-jlo
                incr=-1
              ELSE
                i1=ilo+1
                i2=jlo
                ns=jlo-ilo
                incr=1
              ENDIF
              dint=step/FLOAT(ns+1)
              di=d0
              DO i=i1,i2,incr
                ilold=ilo
                ihit=i
                ilo=ilo+incr
                di=di+dint
                CALL TRAXYZ(di,rhflx,thflx,p0,r0,x0,y0,z0,rx0,rt0,
     &                      rp0,zx0,zt0,zp0,iflag,message)
!  Check flag
                IF(iflag.eq.-1) THEN
!                 Ignore warnings of solution outside separatrix
                  iflag=0
                  message=' '
                ELSEIF(iflag.eq.1) THEN
                  nmes=LEN(message)
                  message='TRACK(6)/'//message(1:nmes-9)
                  GOTO 1000
                ENDIF
                d0=di
                rho0=rhflx
                t0=thflx
                nsurf=nsurf+1
                sl(nsurf)=dnode+d0
                rint(nsurf)=r0
                zint(nsurf)=z0
                tint(nsurf)=p0
                rhoint(nsurf)=rho0
                thint(nsurf)=t0
                iz(nsurf-1)=ilold
                ir(nsurf)=ihit
              ENDDO   
            ENDIF
          ENDDO
   30     continue
          iflag=1
          message='TRACK(7)/ERROR:halved step too many times'
          GOTO 1000
!  This step spans an intersection
   40     iflx=0
          rhflx0=rho0
          thflx0=t0
          rhflxs=rhos
          thflxs=ts
          rhohit=rhoi(ihit)
          di=U_ZEROIN(d0,ds,TRAQRH,tol,iflag,message)
          IF(iflag.ne.0) THEN
            nmes=LEN(message)
            message='TRACK(8)/'//message(1:nmes-9)
            IF(iflag.eq.1) GOTO 1000
          ENDIF
!  Update rhflx,thflx,x0,y0,and z0
          CALL TRAXYZ(di,rhflx,thflx,p0,r0,x0,y0,z0,rx0,rt0,rp0,
     &                zx0,zt0,zp0,iflag,message)
!  Check flag
          IF(iflag.ne.0) THEN
            nmes=LEN(message)
            message='TRACK(9)/'//message(1:nmes-9)
            IF(iflag.eq.1) GOTO 1000
          ENDIF
          d0=di
          rho0=rhflx
          t0=thflx
          nsurf=nsurf+1
          sl(nsurf)=dnode+d0
          rint(nsurf)=r0
          zint(nsurf)=z0
          tint(nsurf)=p0
          rhoint(nsurf)=rho0
          thint(nsurf)=t0
          iz(nsurf-1)=ilold
          ir(nsurf)=ihit
          IF(nsurf.eq.maxsrf) THEN
!  Maximum number of intersections reached
            iflag=1
            message='TRACK(10)/ERROR:max intersections reached'
            GOTO 1000
          ENDIF
!  Step completed
        ENDDO  ! 50
   50   continue
!  Maximum number of steps in one segment exceeded
        iflag=1
        message='TRACK(11)/ERROR:max steps in segment exceeded'
        GOTO 1000
!Segment completed
      ENDDO  ! 60
   60 continue
 1000 RETURN
      END
      SUBROUTINE TRAFLX(r,z,p,tolinv,x,t,rx,rt,rp,zx,zt,zp,iflag,
     &                  message)
!***********************************************************************
!TRAFLX finds the radial and poloidal flux coordinates of a point
!  specified in cylindrical coordinates
!References:
!  Attenberger, Houlberg, Hirshman J Comp Phys 72 (1987) 435
!  W.A.Houlberg 3/2000
!Input:
!  r-major radius (m)
!  z-vertical coordinate (m)
!  p-toroidal angle, right-hand rule around z-direction (radians)
!  tolinv-convergence tolerance on r and z (m)
!     =1.0e-5*plasma minor radius is suggested (m)
!Output:
!  x-flux surface label (-)
!  t-poloidal coordinate theta (radians)
!  rx-dr/dx (m)
!  rt-dr/dtheta (m/radian)
!  rp-dr/dphi (m/radian)
!  zx-dz/dx (m)
!  zt-dz/dtheta (m/radian)
!  zp-dz/dphi (m/radian)
!  iflag-error and warning flag (-)
!       =-1 warning
!       =0 none
!       =1 error
!  message-warning or error message (character)
!Comments:
!  The input values of x and t are used for the initial guess
!  The basic method is a Newton iteration in 2 dimensions
!  Point a is the previous best guess and point b is a trial point
!  If point b is further from (r,z) than point a, a new trial point is
!    generated by halving the step and using interpolated derivatives
!***********************************************************************
      IMPLICIT NONE
      INCLUDE '../inc/comtr9.inc'
!Declaration of input variables
      REAL           p,                       r,
     &               tolinv,                  z
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           rp,                      rt,
     &               rx,                      t,
     &               x,                       zp,
     &               zt,                      zx
!Declaration of local variables
      INTEGER        it,                      nh
      REAL           dr,                      dra,
     &               drb,                     dt,
     &               dx,                      dz,
     &               dza,                     dzb,
     &               err2a,                   err2b,
     &               rb,                      rpb,
     &               rta,                     rtb,
     &               rxa,                     rxb,
     &               ta,                      tau,
     &               taua,                    taub,
     &               tb,                      xa,
     &               xb,                      zb,
     &               zpb,                     zta,
     &               ztb,                     zxa,
     &               zxb
      REAL           z_pi,                    z_precision
!Initialization
!  Physical and conversion constants
      z_pi=ACOS(-1.0)
      z_precision=2.0e-7
!Set starting points and parameters
      dthmax=0.5
      iter=0
      nh=1
      taua=0.0
      err2a=1.0e35
      xa=x
      ta=t
      xb=xa
      tb=ta
      IF(xb.eq.0.0) xb=1.0e-2
!Iterate to find flux coordinates
      DO it=1,mxit
        iter=it
        CALL TRAGRZ(xb,tb,p,rb,zb,rxb,rtb,rpb,zxb,ztb,zpb)
        rdbg(iter)=rb
        zdbg(iter)=zb
        drb=(r-rb)
        dzb=(z-zb)
        err2b=drb**2+dzb**2
!  Check convergence
        IF(ABS(drb).le.tolinv.and.ABS(dzb).le.tolinv) THEN
!         Successful exit
          IF(xb.gt.1.0) THEN
            iflag=-1
            message='TRAFLX(1)/WARNING:outside MHD solution domain'
          ENDIF
          x=xb
          t=AMOD(tb,2.0*z_pi)
          rx=rxb
          rt=rtb
          rp=rpb
          zx=zxb
          zt=ztb
          zp=zpb
          GOTO 1000
        ENDIF
!  Need improved estimate for rho and theta - get new Jacobian
        IF(xb.eq.0.0) xb=10.0*z_precision
        taub=zxb*rtb-rxb*ztb
        IF(ABS(taub).lt.1.0e-37) THEN
!         No solution exists for the Newton method
!         Assume the solution is at the origin where tau=0
          iflag=-1
          message='TRAFLX(2)/WARNING:solution at origin'
          x=0.0
          t=0.0
        ENDIF
!  Check consistency of sign of Jacobian
        IF(SIGN(1.0,taub).ne.SIGN(1.0,taua).and.taua.ne.0.0) THEN
!         Bad data
          iflag=1
          message='TRAFLX(3)/ERROR:bad Jacobian'
          GOTO 1000
        ENDIF
        IF(err2b.lt.err2a) THEN
!         Converging - set point a = point b
          IF(nh.gt.1) nh=nh/2
          xa=xb
          ta=tb
          err2a=err2b
          rxa=rxb
          zxa=zxb
          rta=rtb
          zta=ztb
          taua=taub
          dra=drb
          dza=dzb
        ELSE
!         Not converging - halve step size
          nh=nh*2
        ENDIF
!  Compute rho and theta steps - a and b may coincide here
        rx=0.75*rxa+0.25*rxb
        zx=0.75*zxa+0.25*zxb
        rt=0.75*rta+0.25*rtb
        zt=0.75*zta+0.25*ztb
        tau=0.75*taua+0.25*taub
        dr=dra/FLOAT(nh)
        dz=dza/FLOAT(nh)
        dx=(rt*dz-zt*dr)/tau
        dt=(zx*dr-rx*dz)/tau
        IF(ABS(dt).gt.dthmax) dt=SIGN(dthmax,dt)
!  Get new point b
        xb=xa+dx
        tb=ta+dt
        IF(xb.lt.0.0) THEN
          xb=-xb
          tb=tb+z_pi-2.0*dt
        ENDIF
      ENDDO      
!Exceeded maximum iterations
      iflag=1
      message='TRAFLX(4)/ERROR:max iterations exceeded'
 1000 RETURN
      END
      SUBROUTINE TRAGRZ(x,t,p,r,z,rx,rt,rp,zx,zt,zp)
!***********************************************************************
!TRAGRZ finds r and z and their derivatives with respect to rho, theta,
!   and phi at the given rho, theta, and phi
!References:
!  Attenberger, Houlberg, Hirshman J Comp Phys 72 (1987) 435
!  W.A.Houlberg, P.I. Strand 1/2000
!Input:
!  x-radial flux coordinate (-)
!  t-poloidal angle coordinate (radians)
!  p-toroidal angle coordinate (radians)
!Output:
!  r-radius of point from major axis (m)
!  z-distance of point from midplane (m)
!  rx-dr/dx (m)
!  rt-dr/dtheta (m/radians)
!  rp-dr/dphi (m/radians)
!  zx-dz/dx (m)
!  zt-dz/dtheta (m/radians)
!  zp-dz/dphi (m/radians)
!Comments:
!  The representation is extended beyond the xm3d grid by linear
!    extrapolation of the m=1, n=0 term in rho, with all other terms
!    held fixed at the edge value. This permits unique representation
!    of all space for flux surfaces which are everywhere convex (like
!    ATF, but not bean-shapes)
!  A factor of rho**m is factored out of the Fourier coefficients prior
!    to spline fitting for increased accuracy. For rho less than
!    xm3d(1), all terms are assumed to vary as rho**m
!***********************************************************************
      IMPLICIT NONE
      INCLUDE '../inc/pamxkm.inc'
      INCLUDE '../inc/pamxnx.inc'
      INCLUDE '../inc/comm30.inc'
!Declaration of input variables
      REAL           p,                       t,
     &               x
!Declaration of output variables
      REAL           r,                       rp,
     &               rt,                      rx,
     &               z,                       zp,
     &               zt,                      zx
!Declaration of local variables
      INTEGER        i,                       k
      INTEGER        k_vopt(3)
      REAL           angl,                    cangl,
     &               ct,                      drhom,
     &               drmnx,                   dzmnx,
     &               rhom,                    rmnx,
     &               sangl,                   st,
     &               xdum,                    xx,
     &               zmnx
      REAL           value(3)
      REAL           z_precision
      DATA i/        1/
      SAVE           i
!Initialization
!  Physical and conversion constants
      z_precision=2.0e-7
!  Other
      r=0.0
      z=0.0
      rx=0.0
      rt=0.0
      rp=0.0
      zx=0.0
      zt=0.0
      zp=0.0
      k_vopt(1)=1
      k_vopt(2)=1
      k_vopt(3)=0
!Check if inside or outside equilibrium solution range
      xx=x
      IF(xx.gt.xm3d(nx3d)) xx=xm3d(nx3d)
      IF(xx.lt.xm3d(1)) xx=xm3d(1)
!Loop over modes
      DO k=1,nk3drz
        angl=m3d(k)*t-n3d(k)*p
        cangl=COS(angl)
        sangl=SIN(angl)
        IF(x.eq.0.0.and.m3d(k).eq.0) THEN
          rhom=1.0
        ELSEIF(x.gt.xm3d(nx3d)) THEN
          rhom=xm3d(nx3d)**ABS(m3d(k))
        ELSE
          rhom=x**ABS(m3d(k))
        ENDIF
        CALL W_SPLINE(k_vopt,nx3d,xx,xm3d,rmna(1,1,k),i,value)
        rmnx=value(1)
        drmnx=value(2)
        CALL W_SPLINE(k_vopt,nx3d,xx,xm3d,zmna(1,1,k),i,value)
        zmnx=value(1)
        dzmnx=value(2)
        r=r+rmnx*cangl*rhom
        z=z+zmnx*sangl*rhom
        rt=rt-m3d(k)*rmnx*sangl*rhom
        rp=rp+n3d(k)*rmnx*sangl*rhom
        zt=zt+m3d(k)*zmnx*cangl*rhom
        zp=zp-n3d(k)*zmnx*cangl*rhom
        IF(x.le.xm3d(nx3d)) THEN
          xdum=AMAX1(x,1.0e2*z_precision)
          drhom=ABS(m3d(k))*xdum**(ABS(m3d(k))-1)
          IF(x.lt.xm3d(1)) THEN
            drmnx=0.0
            dzmnx=0.0
          ENDIF
          rx=rx+drmnx*cangl*rhom+rmnx*cangl*drhom
          zx=zx+dzmnx*sangl*rhom+zmnx*sangl*drhom
        ENDIF
      ENDDO   
      IF(x.gt.xm3d(nx3d)) THEN
!       This point is off the grid toward the wall
        ct=COS(m3d(km1n0)*t)
        st=SIN(m3d(km1n0)*t)
        r=r+(x-xm3d(nx3d))*rmna(1,nx3d,km1n0)*ct
        z=z+(x-xm3d(nx3d))*zmna(1,nx3d,km1n0)*st
        rx=rmna(1,nx3d,km1n0)*ct
        zx=zmna(1,nx3d,km1n0)*st
        rt=rt-(x-xm3d(nx3d))*rmna(1,nx3d,km1n0)*st
        zt=zt+(x-xm3d(nx3d))*zmna(1,nx3d,km1n0)*ct
      ENDIF
      RETURN
      END
      FUNCTION TRAQRH(d,iflag,message)
!***********************************************************************
!TRAQRH returns rho(d)-rhohit
!  W.A.Houlberg 3/2000
!Input:
!  d-distance along chord (m)
!Output:
!  iflag-error and warning flag (-)
!       =-1 warning
!       =0 none
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
      INCLUDE '../inc/comtr7.inc'
!Declaration of input variables
      REAL           d
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           TRAQRH
!Declaration of local variables
      INTEGER        nmes
      REAL           phi,                     r,
     &               rp,                      rt,
     &               rx,                      x,
     &               y,                       z,
     &               zp,                      zt,
     &               zx
      TRAQRH=0.0
!Use existing solutions at endpoints on first two calls
      IF(iflx.eq.0) THEN
        rhflx=rhflx0
        thflx=thflx0
        iflx=1
      ELSEIF(iflx.eq.1) THEN
        rhflx=rhflxs
        thflx=thflxs
        iflx=2
      ENDIF
      CALL TRAXYZ(d,rhflx,thflx,phi,r,x,y,z,rx,rt,rp,zx,zt,zp,
     &            iflag,message)
      IF(iflag.eq.-1) THEN
!       Ignore warnings of solution outside plasma
        iflag=0
        message=' '
      ELSEIF(iflag.eq.1) THEN
        nmes=LEN(message)
        message='TRAQRH/'//message(1:nmes-7)
        GOTO 1000
      ENDIF
      TRAQRH=rhflx-rhohit
 1000 RETURN
      END
      FUNCTION TRAQRS(d,iflag,message)
!***********************************************************************
!TRAQRS returns d(rho)/ds where s is path length
!  W.A.Houlberg 3/2000
!Input:
!  d-distance along chord (m)
!Output:
!  iflag-error and warning flag (-)
!       =-1 warning
!       =0 none
!       =1 error
!  message-warning or error message (character)
!Comments:
!  For efficiency, rhflx and thflx should be a good initial guess for
!    rho and theta on entry
!  TRAQRS modifies rhflx and thflx in COMTR7
!***********************************************************************
      IMPLICIT NONE
      INCLUDE '../inc/comtr7.inc'
      INCLUDE '../inc/comtr8.inc'
!Declaration of input variables
      REAL           d
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           TRAQRS
!Declaration of local variables
      INTEGER        nmes
      REAL           phi,                     r,
     &               rhop,                    rhor,
     &               rhoz,                    rt,
     &               rp,                      rx,
     &               tau,                     x,
     &               y,                       z,
     &               zp,                      zt,
     &               zx
      TRAQRS=0.0
!Use existing solutions at endpoints on first two calls
      IF(iflx.eq.0) THEN
        rhflx=rhflx0
        thflx=thflx0
        iflx=1
      ELSEIF(iflx.eq.1) THEN
        rhflx=rhflxs
        thflx=thflxs
        iflx=2
      ENDIF
      x=x0+gx*(d-d0)
      y=y0+gy*(d-d0)
      z=z0+gz*(d-d0)
      CALL TRARPH(x,y,r,phi)
      CALL TRAFLX(r,z,phi,tolinv,rhflx,thflx,rx,rt,rp,zx,zt,zp,iflag,
     &            message)
      IF(iflag.eq.-1) THEN
!       Ignore warnings of solution outside separatrix
        iflag=0
        message=' '
      ELSEIF(iflag.eq.1) THEN
        nmes=LEN(message)
        message='TRAQRS/'//message(1:nmes-7)
        GOTO 1000
      ENDIF
      tau=zx*rt-rx*zt
      rhor=-zt/tau
      rhoz=rt/tau
      rhop=(rp*zt-rt*zp)/tau
      TRAQRS=(rhor*x-rhop*y/r)*gx/r+(rhor*y+rhop*x/r)*gy/r+rhoz*gz
 1000 RETURN
      END
      SUBROUTINE TRARPH(x,y,r,phi)
!***********************************************************************
!TRARPH finds r and phi given x and y, a Cartesian to cylindrical
!  transformation
!  W.A.Houlberg 3/2000
!Input:
!  x,y-Cartesian coodinates of point (m)
!Output:
!  r,phi-cylindrical coordinates of point (m,radians)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      REAL           x,                       y
!Declaration of output variables
      REAL           phi,                     r
!Declaration of local variables
      REAL           xor
      REAL           z_pi
!Initialization
!  Physical and machine constants
      z_pi=ACOS(-1.0)
      r=SQRT(x**2+y**2)
!Check for round-off error in sqrt
      xor=1.0
      IF(r.gt.0.0.and.x.lt.r) xor=x/r
      phi=ACOS(xor)
      IF(y.ge.0.0) phi=ABS(phi)
      IF(y.lt.0.0) phi=2.0*z_pi-ABS(phi)
      RETURN
      END
      SUBROUTINE TRAXYZ(d,rho,th,phi,r,x,y,z,rx,rt,rp,zx,zt,zp,
     &                  iflag,message)
!***********************************************************************
!TRAXYZ finds (rho,theta,phi), (x,y,z), (r,z,phi) and derivatives of r
!  and z wrt rho and theta, for position d along a chord
!References:
!  Attenberger, Houlberg, Hirshman J Comp Phys 72 (1987) 435
!  W.A.Houlberg 3/2000
!Input:
!  d-distance along chord (m)
!Output:
!  rho-curvilinear radius of point (-)
!  th-curvilinear poloidal angle of point (radians)
!  phi-curvilinear toroidal angle of point (radians)
!  r-radius of point from major axis (m)
!  x,y,z-Cartesian coodinates of point (m,m,m)
!  z-distance of point from midplane (m)
!  rx-dr/dx (m)
!  rt-dr/dtheta (m/radians)
!  rp-dr/dphi (m/radians)
!  zx-dz/dx (m)
!  zt-dz/dtheta (m/radians)
!  zp-dz/dphi (m/radians)
!  iflag-error and warning flag (-)
!       =-1 warning
!       =0 none
!       =1 error
!  message-warning or error message (character)
!Comments:
!  For efficiency, rho and th should be a good initial guess on entry
!***********************************************************************
      IMPLICIT NONE
      INCLUDE '../inc/comtr8.inc'
!Declaration of input variables
      REAL           d
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag,                   nmes
      REAL           phi,                     r,
     &               rho,                     rp,
     &               rt,                      rx,
     &               th,                      x,
     &               y,                       z,
     &               zp,                      zt,
     &               zx
!Declaration of local variables
      x=x0+gx*(d-d0)
      y=y0+gy*(d-d0)
      z=z0+gz*(d-d0)
      CALL TRARPH(x,y,r,phi)
      iflag=0
      message=' '
      CALL TRAFLX(r,z,phi,tolinv,rho,th,rx,rt,rp,zx,zt,zp,iflag,message)
      IF(iflag.eq.1) THEN
        nmes=LEN(message)
        message='TRAXYZ/'//message(1:nmes-7)
      ENDIF
      RETURN
      END
