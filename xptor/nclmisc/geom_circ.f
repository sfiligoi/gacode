      SUBROUTINE GEOM_CIRC(cnin,nr_r,a0,r0,bt0,rhot_r,rhop_r,f_r,q_r,
     &                     rin_r,rout_r,elong_r,vol_r,vp_r,phit_r,fm_r,
     &                     grth_r,gph_r,gth_r,b2_r,bm2_r,grho1_r,
     &                     grho2_r,gr2bm2_r,ftrap_r,fhat_r,bpout_r,
     &                     btout_r,psi_r)
!***********************************************************************
!GEOM_CIRC generates plasma geometry information for FORCEBAL using a
!  low beta circular plasma approximation
!  W.A. Houlberg 3/2000
!Input:
!  cnin-name of VMEC 'WOUT' file to be read
!  nr_r-number of radial points (-)
!  a0-minor radius, half diameter of boundary flux surface (m)
!  r0-major radius, center of boundary flux suface (m)
!  bt0-toroidal field at r0 (T)
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
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      CHARACTER*(*)  cnin
      INTEGER        nr_r
      REAL           a0,                      r0,
     &               bt0
      REAL           rhot_r(*)
!Declaration of output variables
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
      CHARACTER*120  message
      INTEGER        i,                       iflag,
     &               m,                       n,
     &               nin
      REAL           eps,                     f0,
     &               phitot,                  scale
      REAL           z_mu0,                   z_pi
      REAL           r(101),                  y(101)
!Initialization
!  Error and warning flag
      iflag=0
!  Physical and conversion constants
      z_mu0=1.2566e-06
      z_pi=ACOS(-1.0)
!  Input unit
      nin=20
!Get q_r
      CALL READ_REC(nin,cnin,n)
      OPEN(unit=nin,status='old',file=cnin,form='formatted')
      DO i=1,n
        READ(nin,*) r(i),y(i)
      ENDDO
      CLOSE(unit=nin)
      scale=1.0/r(n)
      CALL RARRAY_SCALE(n,scale,r,1)
      CALL W_LIN_INTERP(n,r,y,nr_r,rhot_r,q_r,iflag,message)
!Total toroidal flux, poloidal current
      phitot=2.0*z_pi*a0**2*bt0
      f0=2.0*z_pi*r0*bt0/z_mu0
!Set profiles
      DO i=2,nr_r
        eps=rhot_r(i)*a0/r0
        phit_r(i)=phitot*rhot_r(i)**2
        elong_r(i)=1.0
        vp_r(i)=(2.0*z_pi)**2*rhot_r(i)*a0*r0
        vol_r(i)=2.0*(z_pi*rhot_r(i)*a0)**2*r0
        gth_r(i)=eps
        gph_r(i)=eps
        grho2_r(i)=1.0
        grho1_r(i)=1.0
        b2_r(i)=bt0**2*(1.0+0.5*eps**2)
        bm2_r(i)=(1.0+1.5*eps**2)/bt0**2
        gr2bm2_r(i)=1.0/bt0**2
        grth_r(i)=1.0/(q_r(i)*r0)
        ftrap_r(i)=1.46*SQRT(eps)
        fhat_r(i)=q_r(i)/gph_r(i)
        DO m=1,3
          fm_r(m,i)=m*((1.0-SQRT(1.0-eps**2))/eps)**(2.0*m)
     &              *(1.0+m*SQRT(1.0-eps**2))/((1.0-eps**2)**1.5
     &              *(q_r(i)*r0)**2)
        ENDDO   
        f_r(i)=f0
        psi_r(i)=z_mu0*f_r(i)/q_r(i)*eps*rhot_r(i)*a0
        rout_r(i)=r0+rhot_r(i)*a0
        rin_r(i)=r0-rhot_r(i)*a0
        btout_r(i)=r0*bt0/rout_r(i)
        bpout_r(i)=eps*btout_r(i)/q_r(i)
      ENDDO
!Set values at axis
      phit_r(1)=0.0
      elong_r(1)=elong_r(2)
      vp_r(1)=0.0
      vol_r(1)=0.0
      gth_r(1)=0.0
      gph_r(1)=0.0
      grho2_r(1)=grho2_r(2)
      grho1_r(1)=grho1_r(2)
      b2_r(1)=bt0**2
      bm2_r(1)=1.0/bt0**2
      gr2bm2_r(1)=gr2bm2_r(2)
      grth_r(1)=1.0/(q_r(1)*r0)
      ftrap_r(1)=0.0
      fhat_r(1)=fhat_r(2)
      DO m=1,3
        fm_r(m,1)=fm_r(m,2)
      ENDDO
      f_r(1)=f_r(2)
      psi_r(1)=0.0
      rout_r(1)=r0
      rin_r(i)=r0
      btout_r(1)=btout_r(2)*rout_r(2)/rout_r(1)
      bpout_r(1)=0.0
!Poloidal flux grid
      CALL RARRAY_COPY(nr_r,psi_r,1,rhop_r,1)
      CALL RARRAY_SCALE(nr_r,1.0/psi_r(nr_r),rhop_r,1)
      RETURN
      END
