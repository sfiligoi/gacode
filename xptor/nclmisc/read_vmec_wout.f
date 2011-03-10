      SUBROUTINE READ_VMEC_WOUT(cnin,a0,r0,bt0,f0,iflag,message)
!***********************************************************************
!READ_VMEC_WOUT reads a WOUT file from VMEC
!  W.A. Houlberg 3/2000
!Input:
!  cnin-input file name
!Output:
!  a0-minor radius (half diameter) in midplane (m)
!  r0-major radius (center of outermost surface) in midplane (m)
!  bt0-vacuum toroidal field at r0 (T)
!  f0-poloidal current external to last flux surface (A)
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!Comments:
!  VMEC uses a left-handed coordinate system with the poloidal angle
!    clockwise from the outside midplane on the RHS of the torus
!  To convert to a right-handed system the signs on the poloidal mode
!    numbers, lambdas, toroidal flux and iota-bar (1/q) must be flipped
!***********************************************************************
      IMPLICIT NONE
      INCLUDE '../inc/pamxkm.inc'
      INCLUDE '../inc/pamxmp.inc'
      INCLUDE '../inc/pamxmt.inc'
      INCLUDE '../inc/pamxnx.inc'
      INCLUDE '../inc/comm30.inc'
      INCLUDE '../inc/comm31.inc'
!Declaration of input variables
      CHARACTER*(*)  cnin
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           a0,                      bt0,
     &               f0,                      r0
!Declaration of local variables
      CHARACTER*15   adum
      INTEGER        i,                       jpol,
     &               jtor,                    jtorm,
     &               k,                       l,
     &               mpol,                    nin,
     &               nmes,                    ntor
      REAL           rdum(9)
      REAL           rbt_v(mxnx+1),           x_v(mxnx+1),
     &               xiota_v(mxnx+1)
      REAL           wk(mxnx)
      REAL           rmn(mxnx,mxkm),          zmn(mxnx,mxkm)
      REAL           emn_v(mxnx+1,mxkm)
      REAL           dphi,                    dtheta,
     &               phi,                     rhom,
     &               rin,                     rout,
     &               theta
      REAL           z_mu0,                   z_pi
!Initialization
!  Physical and conversion constants
      z_mu0=1.2566e-06
      z_pi=ACOS(-1.0)
!Open the WOUT file
      nin=20
      OPEN(unit=nin,status='old',file=cnin,form='formatted')
!Read the WOUT file
!  nper-toroidal periods, =1 for tokamaks
!  nx3d-radial grid points
!  mpol-poloidal modes
!  ntor-toroidal modes
!  nk3drz-total modes in R,Z expansions
      READ(nin,'(48x,8i6)') nper,nx3d,mpol,ntor,nk3drz
!  nk3d-total modes in lambda expansion
      nk3d=nk3drz
!     Check maximum number of radial grid points
      IF(nx3d.gt.mxnx) THEN
!       Too many grid points
        iflag=1
        message='READ_VMEC_WOUT(1)/ERROR:x dimension exceeded'
        GOTO 1000
      ENDIF
!     Check maximum number of radial grid points
      IF(nk3drz.eq.mxkm) THEN
!       Too many modes
        iflag=1
        message='READ_VMEC_WOUT(2)/ERROR:mode dimension exceeded'
        GOTO 1000
      ENDIF
!     Check consistency of number of modes
      IF(nk3drz.ne.(ntor+1+(2*ntor+1)*(mpol-1))) THEN
!       Number of modes not consistent
        iflag=1
        message='READ_VMEC_WOUT(3)/ERROR:mode data inconsistent'
        GOTO 1000
      ENDIF
!  One to three lines of data not needed depending on VMEC version
      READ(nin,'(a4)') adum
!      READ(nin,'(a4)') adum
!      READ(nin,'(a4)') adum
!  Expansion coefficients and mode numbers
!  Data on uniformly spaced grid in toroidal flux from axis to edge
      km1n0=0
      DO i=1,nx3d
        k=0
        DO jpol=1,mpol
          IF(jpol.eq.1) THEN
            jtorm=0
          ELSE
            jtorm=-ntor
          ENDIF
          DO jtor=jtorm,ntor
            k=k+1
            IF(i.eq.1) THEN
!             Mode numbers only for first node
              READ(nin,'(2i10)') m3d(k),n3d(k)
!             Flip sign on m for right-handed system
              m3d(k)=-m3d(k)
              IF(ABS(m3d(k)).eq.1.and.n3d(k).eq.0.and.km1n0.eq.0)
     &          km1n0=k
            ENDIF
!           Expansion coefficients
!           R and Z on full mesh starting at origin
!           Lambda on half mesh, first node is a dummy with all zeros
            READ(nin,'(5e20.13)') rmn(i,k),zmn(i,k),emn_v(i,k),
     &                            (rdum(l),l=1,8)
!           Flip sign on lambda for right-handed system
            emn_v(i,k)=-emn_v(i,k)
          ENDDO
        ENDDO
      ENDDO
!  iota-bar and RBt on half mesh starting at first half mesh
!  toroidal flux on full mesh starting at second node
      READ(nin,'(5e20.13)') (xiota_v(i),rdum(1),rdum(2),rdum(3),rdum(4),
     &                       rbt_v(i),phitot,rdum(5),rdum(6),rdum(7),
     &                       rdum(8),rdum(9),i=2,nx3d)
!Set dependent variables
!  Set radial grid proportional to sqrt(tor flux)
!  VMEC grid is uniform in delta(tor flux)
      DO i=1,nx3d
        xm3d(i)=SQRT(FLOAT(i-1)/FLOAT(nx3d-1))
      ENDDO
!  Set half mesh, with extra nodes at origin and edge
      x_v(1)=0.0
      DO i=2,nx3d
        x_v(i)=SQRT(0.5*(xm3d(i)**2+xm3d(i-1)**2))
      ENDDO
      x_v(nx3d+1)=1.0
!  Flip sign of phitot for right-handed system
      phitot=-phitot
!  Normalize lambdas to rho**m for better axial resolution
      DO k=1,nk3d
        DO i=2,nx3d
          rhom=x_v(i)**ABS(m3d(k))
          emn_v(i,k)=emn_v(i,k)/rhom
        ENDDO
!  Shift to full grid
!       Axis - linear extrapolation in toroidal flux
        emn_v(1,k)=1.5*emn_v(2,k)-0.5*emn_v(3,k)
!       Outer boundary - linear extrapolation in toroidal flux
        emn_v(nx3d+1,k)=1.5*emn_v(nx3d,k)-0.5*emn_v(nx3d-1,k)
!       Interpolate normalized lambdas to full grid
        CALL W_LIN_INTERP(nx3d+1,x_v,emn_v(1,k),nx3d,xm3d,wk,iflag,
     &                    message)
        IF(iflag.ne.0) THEN
          nmes=LEN(message)
          message='READ_VMEC_WOUT(4)/'//message(1:nmes-18)
          IF(iflag.eq.1) GOTO 1000
        ENDIF
        DO i=1,nx3d
          elmna(1,i,k)=wk(i)
        ENDDO
!  Set spline coefficients
        CALL V_SPLINE(0,0,nx3d,xm3d,elmna(1,1,k),wk)
      ENDDO
!  Normalize rmn, zmn and set spline coefficients 
      DO k=1,nk3drz
!       Intermediate and boundary points
        DO i=2,nx3d
          rhom=xm3d(i)**ABS(m3d(k))
          rmna(1,i,k)=rmn(i,k)/rhom
          zmna(1,i,k)=zmn(i,k)/rhom
        ENDDO
!       Linearly extrapolate normalized coefficients to axis
        rmna(1,1,k)=rmna(1,2,k)-xm3d(2)**2*(rmna(1,3,k)-rmna(1,2,k))
     &            /(xm3d(3)**2-xm3d(2)**2)
        zmna(1,1,k)=zmna(1,2,k)-xm3d(2)**2*(zmna(1,3,k)-zmna(1,2,k))
     &            /(xm3d(3)**2-xm3d(2)**2)
!       Evaluate spline coefficients
        CALL V_SPLINE(0,0,nx3d,xm3d,rmna(1,1,k),wk)
        CALL V_SPLINE(0,0,nx3d,xm3d,zmna(1,1,k),wk)
      ENDDO
!  Get major radius (geometric center) and minor radius (half diameter)
      CALL TRAGRZ(1.0,0.0,0.0,rout,rdum(1),rdum(2),rdum(3),rdum(4),
     &            rdum(5),rdum(6),rdum(7))
      CALL TRAGRZ(1.0,-z_pi,0.0,rin,rdum(1),rdum(2),rdum(3),rdum(4),
     &            rdum(5),rdum(6),rdum(7))
      a0=ABS(rout-rin)/2.0
      r0=(rout+rin)/2.0
!  Extrapolate xiota_v to center and edge and interpolate to full grid
      xiota_v(1)=xiota_v(2)
      xiota_v(nx3d+1)=1.5*xiota_v(nx3d)-0.5*xiota_v(nx3d-1)
      CALL W_LIN_INTERP(nx3d+1,x_v,xiota_v,nx3d,xm3d,wk,iflag,
     &                  message)
      IF(iflag.ne.0) THEN
        nmes=LEN(message)
        message='READ_VMEC_WOUT(5)/'//message(1:nmes-18)
        IF(iflag.eq.1) GOTO 1000
      ENDIF
!  Flip sign of iota-bar for right-handed system
      DO i=1,nx3d
        fiota(1,i)=-wk(i)
      ENDDO
!  Extrapolate R*Bt to edge to get r0*bt0
      bt0=(1.5*rbt_v(nx3d)-0.5*rbt_v(nx3d-1))/r0
!  Set external poloidal current
      f0=2.0*z_pi*r0*bt0/z_mu0
!Evaluate spline coefficients for iota-bar   
      CALL V_SPLINE(0,0,nx3d,xm3d,fiota,wk)
!Evaluate sines and cosines
      mt3d=mxmt
      dtheta=2.0*z_pi/(mt3d-1)
      IF(nper.gt.1) THEN
!       3-D equilibrium
        mp3d=mxmp
        dphi=2.0*z_pi/(mp3d-1)/FLOAT(nper)
      ELSE
!       2-D axisymmetric equilibrium
        mp3d=1
      ENDIF
!     Loop over modes.
      DO k=1,nk3d
!       Loop over toroidal angles phi
        DO jtor=1,mp3d
          phi=(jtor-1)*dphi
!         Loop over poloidal angles theta
          DO jpol=1,mt3d
            theta=(jpol-1)*dtheta
            cth(jpol,jtor,k)=COS(m3d(k)*theta-n3d(k)*phi)
            sth(jpol,jtor,k)=SIN(m3d(k)*theta-n3d(k)*phi)
          ENDDO   
        ENDDO   
      ENDDO   
!Set up integration weighting arrays for theta
      thfcn(1)=1.0
      DO jpol=2,mt3d-2,2
        thfcn(jpol)=4.0
        thfcn(jpol+1)=2.0
      ENDDO   
      thfcn(mt3d-1)=4.0
      thfcn(mt3d)=1.0
      thnorm=3.0*(mt3d-1)
!Set up integration weighting factors for phi
      phfcn(1)=1.0
      IF(mp3d.eq.1)THEN
        phnorm=1.0
      ELSE
        DO jtor=2,mp3d-2,2
          phfcn(jtor)=4.0
          phfcn(jtor+1)=2.0
        ENDDO   
        phfcn(mp3d-1)=4.0
        phfcn(mp3d)=1.0
        phnorm=3.0*(mp3d-1)
      ENDIF
 1000 CLOSE(unit=nin)
      RETURN
      END
