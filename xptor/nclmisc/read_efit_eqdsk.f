      SUBROUTINE READ_EFIT_EQDSK(cnin,mxnx_xy,mxny_xy,mxnr_r,mxn_bdry,
     &                           mxn_lim,nx_xy,ny_xy,nr_r,r0,bt0,rmag,
     &                           zmag,bmag,x_xy,y_xy,psi_xy,psi_r,
     &                           rhop_r,frb_r,ffp_r,p_r,pp_r,q_r,n_bdry,
     &                           x_bdry,y_bdry,n_lim,x_lim,y_lim,iflag,
     &                           message)
!***********************************************************************
!READ_EFIT_EQDSK reads an EQDSK file from EFIT
!  W.A.Houlberg 3/2000
!Input:
!  cnin-input file name
!  mxnx_xy-maximum number of x points on psi(x,y) grid
!  mxny_xy-maximum number of y points on psi(x,y) grid
!  mxnr_r-maximum number of radial points
!  mxn_bdy-maximum number of points on boundary
!  mxn_lim-maximum number of points on limiter surface
!Output:
!  nx_xy-number of x points on psi(x,y) grid
!  ny_xy-number of y points on psi(x,y) grid
!  nr_r-number of radial points
!  r0-major radius = center of boundary flux suface (m)
!  bt0-toroidal field at r0 (T)
!  rmag-major radius of magnetic axis (m)
!  zmag-vertical position of magnetic axis (m)
!  bmag-axial magnetic field (T)
!  x_xy(i)-horizontal grid for 2-D poloidal flux (m)
!  y_xy(i)-vertical grid for 2-D poloidal flux (m)
!  psi_xy(i,j)-poloidal flux/(2*pi) on 2-D grid (Wb/rad)
!  psi_r(i)-poloidal flux/2*pi (Wb/rad)
!  rhop_r-normalized poloidal flux grid proportional to psi
!  frb_r(i)-R*B_t (m*T)
!  ffp_r(i)-frb_r*dfrb_r/dpsi_r (rad*T)
!  p_r(i)-kinetic pressure (N/m**2)
!  pp_r(i)-dp_r/dpsi_r (rad*N/m**2/Wb)
!  q_r(i)-safety factor
!  n_bdry-number of points on plasma boundary
!  x_bdry(i)-major radius of boundary point (m)
!  y_bdry(i)-vertical position of boundary point (m)
!  n_lim-number of points on limiter
!  x_lim(i)-horizontal positions of limiter points (m)
!  y_lim(i)-vertical positions of limiter points (m)
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      CHARACTER*(*)  cnin
      INTEGER        mxnx_xy,                 mxny_xy,                 
     &               mxnr_r,                  mxn_bdry,
     &               mxn_lim
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        nx_xy,                   ny_xy,
     &               nr_r,                    n_bdry,
     &               n_lim,                   iflag
      REAL           r0,                      bt0,
     &               rmag,                    zmag,
     &               bmag
      REAL           x_xy(*),                 y_xy(*),
     &               psi_xy(mxnx_xy,*)
      REAL           psi_r(*),                rhop_r(*),
     &               frb_r(*),                ffp_r(*),      
     &               p_r(*),                  pp_r(*),
     &               q_r(*)
      REAL           x_bdry(*),               y_bdry(*)
      REAL           x_lim(*),                y_lim(*)
!Declaration of local variables
!  EQDSK file
      REAL           rmin,                    rdim,
     &               zmid,                    zdim,
     &               psimag,                  psilim,
     &               current
      INTEGER        mxn_e
      PARAMETER      (mxn_e=130)
      REAL           psi_e(mxn_e),            rhop_e(mxn_e),
     &               f_e(mxn_e),              ffp_e(mxn_e),      
     &               p_e(mxn_e),              pp_e(mxn_e),
     &               q_e(mxn_e)
!  Other
      INTEGER        i,                       j,
     &               nin,                     nr_e
      REAL           dx,                      dy,
     &               dpsi,                    dum
!Open the EQDSK file
      nin=20
      OPEN(unit=nin,status='old',file=cnin,form='formatted')
!Read the EQDSK file
!  Point data
!  All dum values are duplicate information or no longer used
      READ(nin,'(52x,2i4)') nx_xy,ny_xy
!  Same number of points in the radial and horizontal grids
      nr_e=nx_xy
!  Check whether data exceeds maximum dimensions
      IF(nx_xy.gt.mxnx_xy) THEN
!       Horizontal grid points exceed dimensions set by parameters
        iflag=1
        message='READ_EFIT_EQDSK/ERROR:x grid dimension exceeded'
        GOTO 1000
      ENDIF
      IF(ny_xy.gt.mxny_xy) THEN
!       Vertical grid points exceed dimensions set by parameters
        iflag=1
        message='READ_EFIT_EQDSK/ERROR:y grid dimension exceeded'
        GOTO 1000
      ENDIF
      IF(nr_e.gt.mxnr_r) THEN
!       Radial grid points exceed dimensions set by parameters
        iflag=1
        message='READ_EFIT_EQDSK/ERROR:r grid dimension exceeded'
        GOTO 1000
      ENDIF
!  rdim-width of computational region (m)
!  zdim-height of computational region (m)
!  rmin-major radius at inside edge of computational region (m)
!  zmid-vertical center of computational region (m)
      READ(nin,'(5e16.9)') rdim,zdim,r0,rmin,zmid
!  psimag-poloidal flux at magnetic axis/(2*pi) (Wb/rad)
!  psilim-poloidal flux at limiter/(2*pi) (Wb/rad)
      READ(nin,'(5e16.9)') rmag,zmag,psimag,psilim,bt0
!  torcur-toroidal current (A)
      READ(nin,'(5e16.9)') current
      READ(nin,'(5e16.9)') dum
!Profile data
!  The grid is equally spaced in poloidal flux
!  Node 1 is the axis and node nr_r is the edge
      READ(nin,'(5e16.9)') (f_e(i),i=1,nr_e)
      READ(nin,'(5e16.9)') (p_e(i),i=1,nr_e)
      READ(nin,'(5e16.9)') (ffp_e(i),i=1,nr_e)
      READ(nin,'(5e16.9)') (pp_e(i),i=1,nr_e)
      READ(nin,'(5e16.9)') ((psi_xy(i,j),i=1,nx_xy),j=1,ny_xy)
      READ(nin,'(5e16.9)') (q_e(i),i=1,nr_e)
!Boundary and limiter data
      READ(nin,'(2i5)') n_bdry,n_lim
      IF(n_bdry.gt.mxn_bdry) THEN
!       Boundary points exceed dimensions set by parameters
        iflag=1
        message='READ_EFIT_EQDSK/ERROR:bdry grid dim exceeded'
        GOTO 1000
      ENDIF
      IF(n_lim.gt.mxn_lim) THEN
!       Limiter points exceed dimensions set by parameters
        iflag=1
        message='READ_EFIT_EQDSK/ERROR:lim grid dim exceeded'
        GOTO 1000
      ENDIF
      READ(nin,'(5e16.9)') (x_bdry(i),y_bdry(i),i=1,n_bdry)
      READ(nin,'(5e16.9)') (x_lim(i),y_lim(i),i=1,n_lim)
!Derived quantities
!  0-D
      bmag=f_e(1)/rmag
!  1-D grids - map from uniform in Psi to uniform in Psi^0.5
      dpsi=(psilim-psimag)/(nr_e-1)
      psi_e(1)=psimag
      rhop_e(1)=0.0
      psi_e(nr_e)=psilim
      rhop_e(nr_e)=1.0
      DO i=2,nr_e-1
        psi_e(i)=psi_e(i-1)+dpsi 
        rhop_e(i)=(psi_e(i)-psi_e(1))/(psi_e(nr_e)-psi_e(1))
      ENDDO
c     nr_r=41
      nr_r=nx_xy
c      write(*,*) 'MHD grid = ',nr_r
      psi_r(1)=psimag
      rhop_r(1)=0.0
      psi_r(nr_r)=psilim
      rhop_r(nr_r)=1.0
      dpsi=psilim-psimag
      DO i=2,nr_r-1
        psi_r(i)=psi_r(1)+dpsi*(FLOAT(i-1)/FLOAT(nr_r-1))**2
        rhop_r(i)=(psi_r(i)-psi_r(1))/(psi_r(nr_r)-psi_r(1))
      ENDDO
      CALL W_LIN_INTERP(nr_e,rhop_e,f_e,nr_r,rhop_r,frb_r,iflag,message)
      IF(iflag.ne.0) THEN
        i=LEN(message)
        message='READ_EFIT_EQDSK(1)/'//message(1:i-19)
        IF(iflag.eq.1) GOTO 1000
      ENDIF
      CALL W_LIN_INTERP(nr_e,rhop_e,p_e,nr_r,rhop_r,p_r,iflag,message)
      IF(iflag.ne.0) THEN
        i=LEN(message)
        message='READ_EFIT_EQDSK(2)/'//message(1:i-19)
        IF(iflag.eq.1) GOTO 1000
      ENDIF
      CALL W_LIN_INTERP(nr_e,rhop_e,ffp_e,nr_r,rhop_r,ffp_r,iflag,
     &                  message)
      IF(iflag.ne.0) THEN
        i=LEN(message)
        message='READ_EFIT_EQDSK(3)/'//message(1:i-19)
        IF(iflag.eq.1) GOTO 1000
      ENDIF
      CALL W_LIN_INTERP(nr_e,rhop_e,pp_e,nr_r,rhop_r,pp_r,iflag,message)
      IF(iflag.ne.0) THEN
        i=LEN(message)
        message='READ_EFIT_EQDSK(4)/'//message(1:i-19)
        IF(iflag.eq.1) GOTO 1000
      ENDIF
      CALL W_LIN_INTERP(nr_e,rhop_e,q_e,nr_r,rhop_r,q_r,iflag,message)
      IF(iflag.ne.0) THEN
        i=LEN(message)
        message='READ_EFIT_EQDSK(5)/'//message(1:i-19)
        IF(iflag.eq.1) GOTO 1000
      ENDIF
!  2-D grid
      x_xy(1)=rmin
      dx=rdim/(nx_xy-1)
      DO i=2,nx_xy
        x_xy(i)=x_xy(i-1)+dx
      ENDDO
      y_xy(1)=zmid-0.5*zdim
      dy=zdim/(ny_xy-1)
      DO i=2,ny_xy
        y_xy(i)=y_xy(i-1)+dy
      ENDDO
 1000 CLOSE(unit=nin)
      RETURN
      END
