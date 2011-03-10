      SUBROUTINE GEOM_VMEC(cnin,nr_r,rhot_r,rhop_r,f_r,q_r,rin_r,rout_r,
     &                     elong_r,vol_r,vp_r,phit_r,fm_r,grth_r,gph_r,
     &                     gth_r,b2_r,bm2_r,grho1_r,grho2_r,gr2bm2_r,
     &                     ftrap_r,fhat_r,bpout_r,btout_r,psi_r,a0,r0,
     &                     bt0,iflag,message)
!***********************************************************************
!GEOM_VMEC generates plasma geometry information for FORCEBAL by calling
!  getting VMEC information about the MHD equilibrium and then calling
!  FLUXAV_VMEC to generate additional flux surface quantities 
!  W.A. Houlberg 3/2000
!Input:
!  cnin-name of VMEC 'WOUT' file to be read
!  nr_r-number of radial points (-)
!  rhot_r(i)-normalized toroidal flux grid ~(Phi)**0.5 (-)
!Output:
!  rhop_r-normalized poloidal flux grid ~psi (-)
!  f_r(i)-2*pi*R*B_t/mu0 (A)
!  q_r(i)-safety factor (-)
!  rin_r(i)-major radius grid on inside of torus in axis plane (m)
!  rout_r(i)-major radius grid on outside of torus in axis plane (m)
!  elong_r(i)-elongation (-)
!  vol_r(i)-volume enclosed (m**3)
!  vp_r(i)-d vol_r/d rhot_r/a0 (m**2)
!  phit_r(i)-toroidal flux (Wb)
!  fm_r(3,i)-geometric factor (-)
!  grth_r(i)-<n.grad(theta)> (1/m)
!  gph_r(i)-poloidal flux metric (-)
!  gth_r(i)-toroidal flux metric (-)
!  b2_r(i)-<B**2> (T**2)
!  bm2_r(i)-<1/B**2> (/T**2)
!  grho1_r(i)-a0*<|grad(rhot_r)|>
!  grho2_r(i)-a0**2*<|grad(rhot_r)|**2>
!  gr2bm2_r(i)-a0**2*<|grad(rhot_r)|**2/B**2> (1/T**2)
!  ftrap_r(i)-trapped fraction
!  bpout_r(i)-poloidal field at rout_r(i) (T)
!  btout_r(i)-toroidal field at rout_r(i) (T)
!  psi_r(i)-poloidal flux (Wb/rad)
!  a0-minor radius, half diameter of boundary flux surface (m)
!  r0-major radius, center of boundary flux suface (m)
!  bt0-toroidal field at r0 (T)
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
      INCLUDE '../inc/pamxkm.inc'
      INCLUDE '../inc/pamxnx.inc'
      INCLUDE '../inc/comm30.inc'
!Declaration of input variables
      CHARACTER*(*)  cnin
      INTEGER        nr_r
      REAL           rhot_r(*)
!Declaration of output variables
!  0-D
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           a0,                      r0,
     &               bt0
!  1-D
      REAL           b2_r(*),                 bm2_r(*),
     &               bpout_r(*),              btout_r(*),
     &               elong_r(*),              f_r(*),
     &               fhat_r(*),               fm_r(3,*),
     &               ftrap_r(*),              gph_r(*),
     &               grho1_r(*),              grho2_r(*),
     &               gr2bm2_r(*),             grth_r(*),
     &               gth_r(*),                phit_r(*),
     &               psi_r(*),                q_r(*),
     &               rhop_r(*),               rin_r(*),
     &               rout_r(*),               vol_r(*),
     &               vp_r(*)
!Declaration of local variables
      INTEGER        i,                       nmes
      REAL           f0,                      br,
     &               bz,                      phi,
     &               rp,                      rt,
     &               rx,                      theta,
     &               tolinv,                  zp,                      
     &               z,                       zt,
     &               zx
      REAL           rdum(3)
      REAL           z_mu0,                   z_pi
!  Radial grid-psi
      INTEGER        mxnr_r
      PARAMETER     (mxnr_r=130)
!  1-D psi
      REAL           dvol_r(mxnr_r),          rt_r(mxnr_r),
     &               shift_r(mxnr_r),         triang_r(mxnr_r)
!Initialization
!  Physical and conversion constants
      z_mu0=1.2566e-06
      z_pi=ACOS(-1.0)
!Read WOUT file to get VMEC MHD equilibrium information
      CALL READ_VMEC_WOUT(cnin,a0,r0,bt0,f0,iflag,message)
!  Check warning and error flag
      IF(iflag.ne.0) THEN
        i=LEN(message)
        message='READ_VMEC_WOUT/'//message(1:i-15)
        IF(iflag.eq.1) GOTO 1000
      ENDIF
!Set radial grid for flux surfaces [sqrt(tor flux)] and toroidal flux
      DO i=2,nr_r
        phit_r(i)=phitot*rhot_r(i)**2
      ENDDO
      CALL RARRAY_COPY(nr_r,rhot_r,1,rt_r,1)
      CALL RARRAY_SCALE(nr_r,a0,rt_r,1)
!Call FLUXAV_VMEC to generate metrics from VMEC MHD equilibrium
      CALL FLUXAV_VMEC(nr_r,rt_r,a0,f0,dvol_r,vp_r,shift_r,elong_r,
     &                 triang_r,gth_r,gph_r,b2_r,bm2_r,grho2_r,gr2bm2_r,
     &                 grth_r,ftrap_r,fm_r,f_r,q_r,psi_r,iflag,message)
!  Check warning and error flag
      IF(iflag.ne.0) THEN
        nmes=LEN(message)
        message='GEOM_VMEC/'//message(1:nmes-10)
        IF(iflag.eq.1) GOTO 1000
      ENDIF
!Calculate additional variables
      CALL RARRAY_COPY(nr_r,psi_r,1,rhop_r,1)
      CALL RARRAY_SCALE(nr_r,1.0/psi_r(nr_r),rhop_r,1)
      tolinv=1.0e-5*a0
      DO i=1,nr_r
        IF(i.gt.1) THEN
          grho1_r(i)=SQRT(grho2_r(i))
          vol_r(i)=vol_r(i-1)+dvol_r(i-1)
          fhat_r(i)=q_r(i)/gph_r(i)
        ENDIF
        theta=0.0
        phi=0.0
        CALL TRAGRZ(rhot_r(i),theta,phi,rout_r(i),z,rx,rt,rp,zx,zt,zp)
        IF(i.gt.1) THEN
          iflag=0
          CALL TRABFL(0,rout_r(i),z,phi,tolinv,br,bz,btout_r(i),
     &                rhot_r(i),theta,rx,rt,rp,zx,zt,zp,rdum(1),rdum(2),
     &                rdum(3),iflag,message)
          IF(iflag.ne.0) THEN
            nmes=LEN(message)            
            message='GEOM_VMEC/'//message(1:nmes-10)
            IF(iflag.eq.1) GOTO 1000
          ENDIF
          bpout_r(i)=SQRT(br**2+bz**2)
        ENDIF
        theta=z_pi
        phi=0.0
        CALL TRAGRZ(rhot_r(i),theta,phi,rin_r(i),z,rx,rt,rp,zx,zt,zp)
      ENDDO
      btout_r(1)=btout_r(2)*rout_r(2)/rout_r(1)
      bpout_r(1)=0.0
      grho1_r(1)=0.0
      vol_r(1)=0.0
      fhat_r(1)=fhat_r(2)
 1000 RETURN
      END
