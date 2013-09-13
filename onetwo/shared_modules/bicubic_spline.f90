   MODULE bicubic_spline

     !----------------------------------------------------------------
     ! self contained bicubic spline module (no externals)   
     ! uses coding developed for EFIT
     ! uses the not-a-knot bc (hardwired see sub eknot)
     !-------------------------------------------------------HSJ-------


     USE nrtype,                              ONLY :  DP, I4B

     USE error_handler,                       ONLY :  iomaxerr, lerrno, &
                                                      terminate

     USE io_gcnmp,                            ONLY :  nlog,ncrt

     IMPLICIT NONE
     PRIVATE

     INTEGER(I4B), PARAMETER  ::      kubicx = 4,kubicy = 4 
     INTEGER(I4B), PARAMETER  ::      krord  = 4,kzord  = 4
     INTEGER(I4B), PARAMETER  ::      n00    = 0,n11    = 1,  n22   = 2

     INTEGER(I4B)   lubicx, lubicy,kujunk

     INTEGER(I4B)  nw,nh,lr0,bk_lx,bk_ly
     REAL(DP), ALLOCATABLE,DIMENSION(:,:,:,:)    :: cs


     REAL(DP), ALLOCATABLE :: xknot(:),yknot(:),rknot(:), &
                rgrid(:),zgrid(:),zknot(:),copynew(:,:),           &
                bkx(:),bky(:)

     REAL(DP)  work0(4),work1(4),work2(4)

     PUBLIC ::     set_2dspline_dim
     PUBLIC ::     done_2dspline
     PUBLIC ::     sets2d_spline
     PUBLIC ::     eval_2d_spline
     PUBLIC ::     spline_2d_interp
     PUBLIC ::     spline_2d_interp_setup


     CONTAINS


      SUBROUTINE set_2dspline_dim(nR_grid,nZ_grid)
!-------------------------------------------------------------------------
! --  Call this routine first whenever setting up new 2d spline case
!---------------------------------------------------------HSJ-------------       
       INTEGER(I4B) nR_grid,nZ_grid, alloc_stat

       lubicx = nR_grid - kubicx + 1
       lubicy = nZ_grid - kubicy + 1
       kujunk = kubicx*kubicy*lubicx*lubicy

       nw     = nR_grid
       nh     = nZ_grid
       lr0    = nR_grid - krord+1
       lubicx = nw - kubicx + 1
       lubicy = nh - kubicy + 1

       IF(ALLOCATED(bkx))DEALLOCATE(bkx)
          ALLOCATE(bkx(lubicx+1),stat=alloc_stat)
       IF (alloc_stat .NE. 0) THEN 
        STOP  'Allocation unsuccessful in set_2dspline_dim: bkx'
       ENDIF
     
       IF(ALLOCATED(bky))DEALLOCATE(bky)
          ALLOCATE(bky(lubicy+1),stat=alloc_stat)
       IF (alloc_stat .NE. 0) THEN 
        STOP  'Allocation unsuccessful in set_2dspline_dim: bky'
       ENDIF

       IF(ALLOCATED( cs ))DEALLOCATE(cs)
          ALLOCATE(cs(kubicx,lubicx,kubicy,lubicy),stat=alloc_stat)
       IF (alloc_stat .NE. 0) THEN 
        STOP  'Allocation unsuccessful in set_2dspline_dim: cs'
       ENDIF

       IF(ALLOCATED(xknot))DEALLOCATE(xknot)
          ALLOCATE(xknot(kubicx + nR_grid),stat=alloc_stat)
       IF (alloc_stat .NE. 0) THEN 
        STOP  'Allocation unsuccessful in set_2dspline_dim: xknot'
       ENDIF

       IF(ALLOCATED(yknot))DEALLOCATE(yknot)
          ALLOCATE(yknot(kubicy + nZ_grid),stat=alloc_stat)
       IF (alloc_stat .NE. 0) THEN 
        STOP  'Allocation unsuccessful in set_2dspline_dim: yknot'
       ENDIF

       IF(ALLOCATED(rknot))DEALLOCATE(rknot)
          ALLOCATE(rknot(nR_grid+krord),stat=alloc_stat)
       IF (alloc_stat .NE. 0) THEN 
        STOP  'Allocation unsuccessful in set_2dspline_dim: rknot'
       ENDIF

       IF(ALLOCATED(zknot))DEALLOCATE(zknot)
          ALLOCATE(zknot(nZ_grid+kzord),stat=alloc_stat)
       IF (alloc_stat .NE. 0) THEN 
        STOP  'Allocation unsuccessful in set_2dspline_dim: zknot'
       ENDIF

       IF(ALLOCATED(rgrid))DEALLOCATE(rgrid)
          ALLOCATE(rgrid(nR_grid),stat=alloc_stat)
       IF (alloc_stat .NE. 0) THEN 
        STOP  'Allocation unsuccessful in set_2dspline_dim'
       ENDIF
 
       IF(ALLOCATED(zgrid))DEALLOCATE(zgrid)
          ALLOCATE(zgrid(nZ_grid),stat=alloc_stat)
       IF (alloc_stat .NE. 0) THEN 
        STOP  'Allocation unsuccessful in set_2dspline_dim'
       ENDIF

       IF(ALLOCATED(copynew))DEALLOCATE(copynew)
          ALLOCATE(copynew(nR_grid,nZ_grid),stat=alloc_stat)
       IF (alloc_stat .NE. 0) THEN 
        STOP  'Allocation unsuccessful in set_2dspline_dim'
       ENDIF

!       IF(ALLOCATED(dwork1))DEALLOCATE(dwork1)
!          ALLOCATE(dwork1(nR_grid,nZ_grid))

!       IF(ALLOCATED(dwork2))DEALLOCATE(dwork2)
!          ALLOCATE(dwork2(nR_grid))

!       IF(ALLOCATED(dwork3))DEALLOCATE(dwork3)
!          ALLOCATE(dwork3(nR_grid,2*krord-1))

       RETURN
      END SUBROUTINE set_2dspline_dim




      SUBROUTINE done_2dspline
!-------------------------------------------------------------------------
! --  Call this routine whenever stored spline info is no longer needed
!---------------------------------------------------------HSJ------------- 

          IF(ALLOCATED(xknot))DEALLOCATE(xknot)
          IF(ALLOCATED(yknot))DEALLOCATE(yknot)
          IF(ALLOCATED(rknot))DEALLOCATE(rknot)
          IF(ALLOCATED(rgrid))DEALLOCATE(rgrid)
          IF(ALLOCATED(zgrid))DEALLOCATE(zgrid)
          IF(ALLOCATED(zknot))DEALLOCATE(zknot)
          IF(ALLOCATED(copynew))DEALLOCATE(copynew)
          IF(ALLOCATED( cs ))DEALLOCATE(cs)
          IF(ALLOCATED(bkx))DEALLOCATE(bkx)
          IF(ALLOCATED(bky))DEALLOCATE(bky)
 
      END SUBROUTINE done_2dspline



!      SUBROUTINE seva2d(bkx,lx,bky,ly,cs,xl,yl,fs,ier,icalc)
      SUBROUTINE  eval_2d_spline(xl,yl,fs,icalc)
!-----------------------------------------------------------------------------------
!--  Call this routine after setting up spline coefficients 
! -- ( done in sets2d_spline) HSJ
!--  S.Thompson  92/05/18                                                    --
!--    Bicubic spline routines.                                              --
!--    Put together with routines from E.Solano.                             --
!--  SMWolfe     93/12/17                                                    --
!--    Modifed to avoid over-writing the original array.                     --
!--              94/04/28                                                    --
!--    Updated.                                                              --
!-----------------------------------------------------------------------------------
!  Inputs:
!
!      cs       - array of spline coefficients of dimension (kubicx,
!                 lubicx,kubicy,lubicy) from sets2d.
!                 cs is global to this module HSJ
!
!      bkx, bky - interval coefficients of length lubicx+1 and lubicy+1 from
!                 sets2d. bkx,bky are global to this module
!
!      bk_lx, bk_ly   - number of terms in bkx and bky from sets2d.
!                       global
!
!      xl, yl   - the point at which interpolations are desired.
!
!  Outputs:
!
!      fs       - vector containing results depending on icalc:
!                 icalc              fs
!                   1                f
!                   2                fx
!                   3                fy
!                   4                fxy
!                   5                fxx
!                   6                fyy
!
!      ier      - error parameter.
!
!-------------------------------------------------------------------------------
    
    !  use eparmdud129
    !  use expath
    !  IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h,o-z)

      IMPLICIT NONE

      INTEGER(I4B)  ibk,icalc, jj,lef,mflag,ndummy 



       REAL(DP)  h,xl,yl,fs(6)


!
!  Evaluate function and its partial derivatives at (XL, YL):
!
!  First do all the lookup and interpolation stuff.
!  This is the most time consuming part of the evaluation, so
!  don't do more than needed.
!
      CALL interv(bky,bk_ly,yl,lef,mflag)
      CALL interv(bkx,bk_lx,xl,ibk,ndummy)


      h = xl - bkx(ibk)
      DO 41 jj=1,4
         work0(jj) = ppvalw(cs(1,ibk,jj,lef),h,n00)
         IF (icalc.EQ.1) GOTO 41
         work1(jj) = ppvalw(cs(1,ibk,jj,lef),h,n11)
         IF (icalc.LE.4) GOTO 41
         work2(jj) = ppvalw(cs(1,ibk,jj,lef),h,n22)
 41   CONTINUE
      h = yl - bky(lef)
      fs(1) = ppvalw(work0,h,n00)
      IF (icalc.EQ.1) RETURN
      fs(2) = ppvalw(work1,h,n00)
      IF (icalc.EQ.2) RETURN
      fs(3) = ppvalw(work0,h,n11)
      IF (icalc.EQ.3) RETURN
      fs(4) = ppvalw(work1,h,n11)
      IF (icalc.EQ.4) RETURN
      fs(5) = ppvalw(work2,h,n00)
      IF (icalc.EQ.5) RETURN
      fs(6) = ppvalw(work0,h,n22)
!
      RETURN
      END SUBROUTINE eval_2d_spline



      SUBROUTINE spline_2d_interp_setup(zdata,xdata,ydata)
        !----------------------------------------------------------------------------
        ! Set up and evaluate a 2d bicubic spline with not a knot boundary conditions
        ! Arguments:
        !  zdata = 2D array, having dimensions of (nx,ny), defined on xdata, ydata
        !  xdata = 1D array of length nx
        !  ydata = 1D array of length ny
        !----------------------------------------------------------------------------
        IMPLICIT NONE

        REAL(DP), DIMENSION(:,:), INTENT(IN) :: zdata
        REAL(DP), DIMENSION(:), INTENT(IN) :: xdata, ydata
        REAL(DP), DIMENSION(SIZE(zdata,1)*SIZE(zdata,2)) :: ztemp
        INTEGER(I4B) :: ndatax, ndatay, i, j,zdatax,zdatay
        ndatax = SIZE(xdata)   ; ndatay = SIZE(ydata)
        zdatax = SIZE(zdata,1) ; zdatay = SIZE(zdata,2)
        CALL done_2dspline

        CALL set_2dspline_dim(ndatax,ndatay)

       DO i = 1,  zdatax
           DO j=1,zdatay         ! HSJ 6/15/2012
             ztemp((i-1)*zdatay+j) = zdata(i,j) ! HSJ 6/15/2012
          ENDDO
        ENDDO

        CALL sets2d_spline(ztemp,xdata,ndatax,ydata,ndatay)

      END SUBROUTINE spline_2d_interp_setup

      FUNCTION spline_2d_interp(zdata,xdata,ydata,icalc,xint,yint)
        !-----------------------------------------------------------------------------
        ! Set up and evaluate a 2d bicubic spline with not a knot boundary conditions
        ! Arguments:
        !  zdata = 2D array, having dimensions of (nx,ny), defined on xdata, ydata
        !  xdata = 1D array of length nx
        !  ydata = 1D array of length ny
        !  xint = 1D array of length nxint - the x values at which to interpolate
        !  yint = 1D array of length nyint - the y values at which to interpolate
        ! icalc = 1..6, value or derivative to return
        ! Returns:
        ! zint = 2D array of points interpolated at xint, yint. 
        !---------------------------------------------------------------------------------
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: zdata
        REAL(DP), DIMENSION(:), INTENT(IN) :: xdata, ydata, xint, yint
        REAL(DP), DIMENSION(SIZE(xint),SIZE(yint)) :: spline_2d_interp
        REAL(DP), DIMENSION(SIZE(xdata)*SIZE(YDATA)) :: ztemp
        REAL(DP), DIMENSION(6) :: feval
        INTEGER(I4B), INTENT(IN) :: icalc
        INTEGER(I4B) :: ndatax, ndatay, i, j
        CALL spline_2d_interp_setup(zdata,xdata,ydata)
        spline_2d_interp = eval_2d_spline_array(xint,yint,icalc)
      END FUNCTION spline_2d_interp
      
      FUNCTION eval_2d_spline_array(xint,yint,icalc)
        ! Returns the 2-d evaluation of the current 2-d spline at the points
        ! xint(i),yint(j)
        !
        ! Arguments:
        !  xint - nx sized array - The x grid at which to evaluate the spline
        !  yint - ny sized array - The y grid at which to evaluate the spline
        !  icalc - Only return the specified spline value
        !        = 1 - Value
        !        = 2 - df/dx
        !        = 3 - df/dy
        !        = 4 - d2f/dxdy
        !        = 5 - d2f/dx2
        !        = 6 - d2f/dy2
        !
        ! Returns:
        !  An nx by ny 2d array of interpolated points
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: xint, yint
        INTEGER(I4B), OPTIONAL :: icalc
        REAL(DP), DIMENSION(SIZE(xint),SIZE(yint)) :: eval_2d_spline_array
        INTEGER(I4B) :: nintx, ninty, i, j
        REAL(DP), DIMENSION(6) :: feval
        nintx = SIZE(xint)
        ninty = SIZE(yint)
        IF (.NOT. PRESENT(icalc)) icalc = 1
        DO i = 1, nintx
          DO j = 1, ninty
            CALL eval_2d_spline(xint(i),yint(j),feval,icalc)
            eval_2d_spline_array(i,j) = feval(icalc)
          ENDDO
        ENDDO       
      END FUNCTION eval_2d_spline_array


      REAL(DP) FUNCTION ppvalw (coef, x, jd )
!-----------------------------------------------------------------------       
!  Modified for optimization by S.J. Thompson, 30-Aug-1993
!  Revised to eliminate call to interv by S.M.Wolfe, 17-Dec-1993
!          and to use ASF's for evaluation
!  This routine performs only the innermost guts of the spline evaluation
!  Assumes k=4 (cubic spline only). No other cases considered. 
!  does not call  interv
!  calculates value at  x  of  jd-th derivative of pp fct from pp-repr
!
!******  i n p u t  ****** to PPVALU, on which this is based.
!  break, coef, l, k.....forms the pp-representation of the function  f
!        to be evaluated. specifically, the j-th derivative of  f  is
!        given by
!
!     (d**j)f(x) = coef(j+1,i) + h*(coef(j+2,i) + h*( ... (coef(k-1,i) +
!                             + h*coef(k,i)/(k-j-1))/(k-j-2) ... )/2)/1
!
!        with  h = x - break(i),  and
!
!       i  =  max( 1 , max( j ;  break(j) .le. x , 1 .le. j .le. l ) ).
!
!  x.....the point at which to evaluate.
!        as used here, x is the distance from the break, not the absolute 
!        position. 
!  jd.....integer*4 giving the order of the derivative to be evaluat-
!        ed.  a s s u m e d  to be zero or positive.
!
!******  o u t p u t  ******
!  ppvalw.....the value of the (jd)-th derivative of  f  at  x.
!
!******  m e t h o d  ******
!     the interval index  i , appropriate for  x , is found through a
!  call to  interv . the formula above for the  jd-th derivative
!  of  f  is then evaluated (by nested multipication).
!
!-----------------------------------------------------------------------        
!   Variable declarations.
!-----------------------------------------------------------------------        
!      IMPLICIT INTEGER(I4B)  (i-n)
!      REAL(DP) (a-h, o-z)
      REAL(DP) x,xx,d0,d1,d2
      REAL(DP) coef(4)
      INTEGER(I4B) jd

!----------------------------------------------------------------------
! ASF's may be slightly more efficient than the alternative
!----------------------------------------------------------------------
      d2(xx) = coef(4)*xx + coef(3)
      d1(xx) = (coef(4)*xx/2. + coef(3))*xx + coef(2)
      d0(xx) = ((coef(4)*xx/3. + coef(3))*xx/2. +  &
                 coef(2))*xx + coef(1)
!-----------------------------------------------------------------------        
!   Derivatives of order k or higher are identically zero.
!   Evaluate jd-th derivative of i-th polynomial piece at x .
!-----------------------------------------------------------------------        
      GOTO (1,2,3) jd+1
      ppvalw = 0.
      PRINT *, 'Error (ppvalw): JD must be 0, 1, or 2.'
      PRINT *, 'Execution terminated.'
      RETURN
 1    ppvalw = d0(x)	! k = 4 , jd = 0
      RETURN
 2    ppvalw = d1(x)	! k = 4 , jd = 1
      RETURN
 3    ppvalw = d2(x)	! k = 4 , jd = 2
      RETURN
      END   FUNCTION ppvalw 
!


      SUBROUTINE interv ( xt, lxt, x, left, mflag )
!-------------------------------------------------------------------------
!omputes  left = max( i ; 1 .le. i .le. lxt  .and.  xt(i) .le. x )  .
!
!******  i n p u t  ******
!  xt.....a REAL*8 sequence, of length  lxt , assumed to be nondecreasing
!  lxt.....number of terms in the sequence  xt .
!  x.....the point whose location with respect to the sequence  xt  is
!        to be determined.
!
!******  o u t p u t  ******
!  left, mflag.....both integers, whose value is
!
!   1     -1      if               x .lt.  xt(1)
!   i      0      if   xt(i)  .le. x .lt. xt(i+1)
!  lxt     1      if  xt(lxt) .le. x
!
!        in particular,  mflag = 0 is the 'usual' case.  mflag .ne. 0
!        indicates that  x  lies outside the halfopen interval
!        xt(1) .le. y .lt. xt(lxt) . the asymmetric treatment of the
!        interval is due to the decision to make all pp functions cont-
!        inuous from the right.
!
!******  m e t h o d  ******
!  the program is designed to be efficient in the common situation that
!  it is called repeatedly, with  x  taken from an increasing or decrea-
!  sing sequence. this will happen, e.g., when a pp function is to be
!  graphed. the first guess for  left  is therefore taken to be the val-
!  ue returned at the previous call and stored in the  l o c a l  varia-
!  ble  ilo . a first check ascertains that  ilo .lt. lxt (this is nec-
!  essary since the present call may have nothing to do with the previ-
!  ous call). then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
!  ilo  and are done after just three comparisons.
!     otherwise, we repeatedly double the difference  istep = ihi - ilo
!  while also moving  ilo  and  ihi  in the direction of  x , until
!                      xt(ilo) .le. x .lt. xt(ihi) ,
!  after which we use bisection to get, in addition, ilo+1 = ihi .
!  left = ilo  is then returned.
!
!      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)

      INTEGER(I4B) left,lxt,mflag,ihi,ilo,istep,middle
      REAL(DP) x,xt(lxt)
      DATA ilo /1/
!     save ilo  (a valid fortran statement in the new 1977 standard)
      ihi = ilo + 1
      IF (ihi .LT. lxt)                 go to 20
         IF (x .GE. xt(lxt))            go to 110
         IF (lxt .LE. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt
!
   20 IF (x .GE. xt(ihi))               go to 40
      IF (x .GE. xt(ilo))               go to 100
!
!              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
   30 istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         IF (ilo .LE. 1)                go to 35
         IF (x .GE. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      IF (x .LT. xt(1))                 go to 90
                                        go to 50
!              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         IF (ihi .GE. lxt)              go to 45
         IF (x .LT. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 IF (x .GE. xt(lxt))               go to 110
      ihi = lxt
!
!           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 middle = (ilo + ihi)/2
      IF (middle .EQ. ilo)              go to 100
!     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      IF (x .LT. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50
!**** set output and return.
   90 mflag = -1
      left = 1
                                        RETURN
  100 mflag = 0
      left = ilo
                                        RETURN
  110 mflag = 1
      left = lxt
                                        RETURN
      END SUBROUTINE interv


      SUBROUTINE sets2d_spline(s,x,nx,y,ny)
!------------------------------------------------------------------------------
!--  Call this routine after defining arrays (but before evaluating spline) HSJ
!--  S.Thompson  92/05/18                                                    --
!--    Bicubic spline routines.                                              --
!--    Put together with routines from E.Solano.                             --
!--  SMWolfe     93/12/17                                                    --
!--    Modifed to avoid over-writing the original array.                     --
!--              94/04/28                                                    --
!--    Updated.                                                              --
!------------------------------------------------------------------------------
!  Inputs: 
!
!      s     - A 1-d array of length nx * ny containing the function values at
!             (x,y). This is a 1-d array, such that if the original tabulated
!             function is 2-d f(i,j), then s(k)=f(i,j), k=(i-1)*ny+j.
!
!      x, y  - (x,y) location, arrays of length nx and ny.
!
!  Outputs:
!
!      cs    - array of spline coefficients of dimension (kubicx,
!              lubicx,kubicy,lubicy). 
!              cs is global to this module HSJ
!
!      bkx, bky - interval coefficients of length lubicx+1 and lubicy+1.
!              bkx,bky are global to this module HSJ
!
!      bk_lx, bk_ly -   number of terms in bkx and bky.
!                       global
!
!      ier   - error parameter.
!
!  Work arrays:
!
!      wk    - of dimension at least nx by ny.
!------------------------------------------------------------------------------
!      INCLUDE 'eparmdud129.f90'
!      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h,o-z)
!
       IMPLICIT NONE
       INTEGER(I4B) i,j,k,nx,ny
       REAL(DP) s(nx*ny), x(nx), y(ny)
       REAL(DP)  wk(nx,ny)   ! temporary array
!  Set up knots:
!
      CALL eknot (nx, x, kubicx, xknot)		
      CALL eknot (ny, y, kubicy, yknot)	

		
!
!  Save the original, use the work array
!
    DO 10 i=1,nx
      DO 10 j=1,ny
         k=(i-1)*ny+j
  10     wk(i,j) = s(k)
!
!  Calculate spline coefficients:
!
 
      CALL spl2bc (x, y, xknot, yknot, wk)	
!
!  Coefficients stored in bkx, bky, and c:
!

      CALL spl2pp (xknot, yknot, wk, bkx, bk_lx, bky, bk_ly, cs)
!

!
      RETURN
      END SUBROUTINE sets2d_spline




      SUBROUTINE eknot(n,x,k,xk)
!--------------------------------------------------------------------------
! given the ordered data points x(1)<...<x(n), this subroutine generates
! a knot sequence with not-a-knot end conditions (like BSNAK from IMSL)
! Some of this is discussed in de Boor(1978), page 211.
!--------------------------------------------------------------------------
!      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)

        IMPLICIT NONE
        INTEGER(I4B) i,ii,kh,k,k2,n
        REAL(DP) x(n),xk(n+k)

!
        DO i=1,k
        xk(i)=x(1)
        ii=i+n
        xk(ii)= x(n)+1.e-5_DP
        ENDDO
        kh=k/2
        k2=kh+kh
        IF (k2.EQ.k) THEN
! even k, place knots at data points
        DO i=k+1,n
        xk(i)=x(i-kh)
        ENDDO
        ELSE
! odd k, place knots in between data points
        DO i=k+1,n
        xk(i)=.5*(x(i-kh)+x(i-1-kh))
        ENDDO
        END IF
        RETURN
        END SUBROUTINE eknot


       SUBROUTINE spl2bc(drgrid,dzgrid,drknot,dzknot,dcopynew)
!--------------------------------------------------------------
! -- calculates the b-spline coeficients
!--------------------------------------------------------------------
      IMPLICIT NONE 

       REAL(DP) drgrid(:),dzgrid(:),drknot(:),dzknot(:),dcopynew(:,:)
!
!       REAL(DP), ALLOCATABLE :: twork1(:,:),twork2(:),twork3(:,:)

      INTEGER(I4B) iflag1,iflag2
      !TEMPORARY arrays:
       REAL(DP) twork1(nw,nh),twork2(nh),twork3(nh,2*krord-1)
!      ALLOCATE(twork1(nw,nh),twork2(nh),twork3(nh,2*krord-1))
!------------------------------------------------------------------
!-- change dimension of work2 and work3 from nw to nh            --
!-- to ensure the cases when nh > nw     ll, 93/04/01            --
!------------------------------------------------------------------
!
      CALL spli2d(drgrid,dcopynew,drknot,nw,krord,nh,twork2,twork3,twork1,iflag1)
      IF (iflag1 .NE. 1)  &
           WRITE(ncrt,FMT='(" error in first spli2d, iflag=")'),iflag1

      CALL spli2d(dzgrid,twork1,dzknot,nh,kzord,nw,twork2,twork3,dcopynew,iflag2)
      IF (iflag2 .NE. 1)  &
           WRITE(ncrt,FMT='(" error in second spli2d, iflag=")'),iflag2

       IF(iflag1 .NE. 1 .OR. iflag2 .NE. 1)THEN
         lerrno = iomaxerr+184
         CALL terminate(lerrno,nlog)
          
       ENDIF

!       DEALLOCATE(twork1,twork2,twork3)

      RETURN
      END SUBROUTINE spl2bc


      SUBROUTINE spl2pp(drknot,dzknot,dcopy,breakr,lr,breakz,lz,coef)
!--------------------------------------------------------------------
! translates to pp representation
!      USE eparmdud129,ONLY:nw,nh,lr0,lz0


        IMPLICIT NONE 

        INTEGER(I4B) lr,lz,ndum
!       INTEGER(I4B), PARAMETER :: lr0=nw-krord+1,lz0=nh-kzord+1
        REAL(DP) coef(kubicx,lubicx,kubicy,lubicy)
        REAL(DP) breakr(:),breakz(:)
        REAL(DP) twork4(krord,nw,nh),twork5(nh,krord,lr0), &
                 twork6(kzord,kzord,nw,krord)
        REAL(DP) drknot(:),dzknot(:),dcopy(nw,nh)
!
 
      CALL bspp2d(drknot,dcopy,nw,krord,nh,twork4,breakr,twork5,lr)
      ndum=lr*krord

      CALL bspp2d(dzknot,twork5,nh,kzord,ndum,twork6,breakz,coef,lz)
! 



 
        RETURN
      END  SUBROUTINE spl2pp


      SUBROUTINE spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )
!-------------------------------------------------------------------------
!alls bsplvb, banfac/slv
!  this is an extended version of  splint , for the use in tensor prod-
!  uct interpolation.
!
!   spli2d  produces the b-spline coeff.s  bcoef(j,.)  of the spline of
!   order  k  with knots  t (i), i=1,..., n + k , which takes on the
!   value  gtau (i,j)  at  tau (i), i=1,..., n ; j=1,..., m .
!
!******  i n p u t  ******
!  tau   array of length  n , containing data point abscissae.
!  a s s u m p t i o n . . .  tau  is strictly increasing
!  gtau(.,j)  corresponding array of length  n , containing data point
!        ordinates, j=1,...,m
!  t     knot sequence, of length  n+k
!  n     number of data points and dimension of spline space  s(k,t)
!  k     order of spline
!  m     number of data sets
!
!******  w o r k   a r e a  ******
!  work  a vector of length  n
!
!******  o u t p u t  ******
!  q     array of order  (n,2*k-1), containing the triangular factoriz-
!        ation of the coefficient matrix of the linear system for the b-
!        coefficients of the spline interpolant.
!           the b-coeffs for the interpolant of an additional data set
!        (tau(i),htau(i)), i=1,...,n  with the same data abscissae can
!        be obtained without going through all the calculations in this
!        routine, simply by loading  htau  into  bcoef  and then execut-
!        ing the    call banslv ( q, n, n, 2*k-1, k, bcoef )
!  bcoef the b-coefficients of the interpolant, of length  n
!  iflag an INTEGER*4 indicating success (= 1)  or failure (= 2)
!        the linear system to be solved is (theoretically) invertible if
!        and only if
!              t(i) .lt. tau(i) .lt. tau(i+k),    all i.
!        violation of this condition is certain to lead to  iflag = 2 .
!
!******  m e t h o d  ******
!     the i-th equation of the linear system  a*bcoef = b  for the b-co-
!  effs of the interpolant enforces interpolation at  tau(i), i=1,...,n.
!  hence,  b(i) = gtau(i), all i, and  a  is a band matrix with  2k-1
!  bands (if it is invertible).
!     the matrix  a  is generated row by row and stored, diagonal by di-
!  agonal, in the  c o l u m n s  of the array  q , with the main diag-
!  onal going into column  k .  see comments in the program below.
!     the banded system is then solved by a call to  banfac (which con-
!  structs the triangular factorization for  a  and stores it again in
!   q ), followed by a call to  banslv (which then obtains the solution
!   bcoef  by substitution).
!     banfac  does no pivoting, since the total positivity of the matrix
!  a  makes this unnecessary.
!------------------------------------------------------------------------


!      INTEGER*4 iflag,k,m,n,i,ilp1mx,j,jj,kpkm1,left,np1
!      REAL*8 bcoef(m,n),gtau(n,m),q(n,7),t(n+k),tau(n),work(n),taui

      IMPLICIT NONE

!      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
       INTEGER(I4B) n,k,m,iflag,nn,np1,kpkm1,left,i,iindex,ilp1mx, &
                    j,jj,nnn
       REAL(DP)  bcoef(m,n),gtau(n,m),q(n,2*k-1),t(n+k),tau(n),work(n)
       REAL(DP) taui
!
      nnn=1
      np1 = n + 1
      kpkm1 = 2*k - 1
      left = k
!
!  ***   loop over  i  to construct the  n  interpolation equations
      DO 30 i=1,n
         iindex=i
         taui = tau(iindex)
         ilp1mx = MIN(iindex+k,np1)
!        *** zero out all entries in row  i  of  a (in the 2k-1 bands)
         DO 13 j=1,kpkm1
   13       q(iindex,j) = 0.
!        *** find  left  in the closed interval (i,i+k-1) such that
!                t(left) .le. tau(i) .lt. t(left+1)
!        matrix is singular if this is not possible
         left = MAX(left,i)
         IF (taui .LT. t(left))         go to 998
   15       IF (taui .LT. t(left+1))    go to 16
            left = left + 1
            IF (left .LT. ilp1mx)       go to 15
         left = left - 1
         IF (taui .GT. t(left+1))       go to 998
!        *** the i-th equation enforces interpolation at taui, hence
!        a(i,j) = b(j,k,t)(taui), all j. only the  k  entries with  j =
!        left-k+1,...,left actually might be nonzero. these  k  numbers
!        are returned, in  work  (used for temp.storage here), by the
!        following
   16    CALL bsplvb ( t, k, nnn, taui, left, work )
!        we therefore want  work(j) = b(left-k+j)(taui) to go into
!        a(i,left-k+j), i.e., into  q(i,left-i+j), since the i-th row of
!        a  is so stored in the i-th row of  q  that the (i,i)-entry of
!        a  goes into the  k-th  entry of  q.
         jj = left - iindex
         DO 29 j=1,k
            jj = jj+1
            q(iindex,jj) = work(j)
   29    CONTINUE
   30    CONTINUE
!
!     ***obtain factorization of  a  , stored again in  q.
      CALL banfac ( q, n, n, kpkm1, k, iflag )
                                        go to (40,999), iflag
!     *** solve  a*bcoef = gtau  by backsubstitution
   40 DO 50 j=1,m
         DO 41 i=1,n
   41       work(i) = gtau(i,j)
         CALL banslv ( q, n, n, kpkm1, k, work )
         DO 50 i=1,n
   50    bcoef(j,i) = work(i)
                                        RETURN
  998 iflag = 2
  999 PRINT 699
  699 FORMAT(41h linear system in  splint  not invertible)
                                        RETURN
      END SUBROUTINE spli2d



      SUBROUTINE bspp2d ( t, bcoef, n, k, m, scrtch, break, coef, l )
!---------------------------------------------------------------------
!alls  bsplvb
!  this is an extended version of  bsplpp  for use with tensor products
!
!onverts the b-representation  t, bcoef(.,j), n, k  of some spline into
!  its pp-representation  break, coef(j,.,.), l, k ; j=1, ..., m  .
!
!******  i n p u t  ******
!  t     knot sequence, of length  n+k
!  bcoef(.,j) b-spline coefficient sequence, of length  n ;j=1,...,m
!  n     length of  bcoef  and  dimension of spline space  s(k,t)
!  k     order of the spline
!  m     number of data sets
!
!******  w o r k   a r e a  ******
!  scrtch   of size  (k,k,m), needed to contain bcoeffs of a piece of
!        the spline and its  k-1  derivatives   for each of the m sets
!
!******  o u t p u t  ******
!  break breakpoint sequence, of length  l+1, contains (in increasing
!        order) the distinct points in the sequence  t(k), ..., t(n+1)
!  coef(mm,.,.)  array of size (k,n), with  coef(mm,i,j) = (i-1)st der-
!        ivative of  mm-th  spline at break(j) from the right, mm=1,.,m
!  l     number of polynomial pieces which make up the spline in the
!        interval  (t(k), t(n+1))
!
!******  m e t h o d  ******
!     for each breakpoint interval, the  k  relevant b-coeffs of the
!  spline are found and then differenced repeatedly to get the b-coeffs
!  of all the derivatives of the spline on that interval. the spline and
!  its first  k-1  derivatives are then evaluated at the left end
!  point of that interval, using  bsplvb  repeatedly to obtain the val-
!  ues of all b-splines of the appropriate order at that point.
!---------------------------------------------------------------------
    !  IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
       IMPLICIT NONE 
       INTEGER(I4B), PARAMETER :: kmax  = 4
       INTEGER(I4B)  k,l,m,n,mm,   i,j,jp1,kmj,left,n11,n22
       REAL(DP) bcoef(n,m),break(:),coef(m,k,*),scrtch(k,k,m),t(:), &
           biatx(kmax)
       REAL(DP)  diff,fkmj,sum
!

      n11=1
      n22=2
      l = 0
      break(1) = t(k)
      DO 50 left=k,n
!        find the next nontrivial knot interval.

         IF (t(left+1) .EQ. t(left))    go to 50
         l = l + 1
         break(l+1) = t(left+1)
         IF (k .GT. 1)                  go to 9
         DO 5 mm=1,m
    5       coef(mm,1,l) = bcoef(left,mm)

                                        go to 50
!        store the k b-spline coeff.s relevant to current knot interval
!        in  scrtch(.,1) .
    9    DO 10 i=1,k
            DO 10 mm=1,m
   10          scrtch(i,1,mm) = bcoef(left-k+i,mm)
!        for j=1,...,k-1, compute the k-j b-spline coeff.s relevant to
!        current knot interval for the j-th derivative by differencing
!        those for the (j-1)st derivative, and store in scrtch(.,j+1) .

         DO 20 jp1=2,k
            j = jp1 - 1
            kmj = k - j
            fkmj = float(kmj)
            DO 20 i=1,kmj
               diff = (t(left+i) - t(left+i - kmj))/fkmj
               IF (diff .LE. 0.)         go to 20
               DO 15 mm=1,m
   15             scrtch(i,jp1,mm) = &
                  (scrtch(i+1,j,mm) - scrtch(i,j,mm))/diff
   20          CONTINUE
!        starting with the one b-spline of order 1 not zero at t(left),
!        find the values at t(left) of the j+1 b-splines of order j+1
!        not identically zero there from those of order j, then combine
!        with the b-spline coeff.s found earlier to compute the (k-j)-
!        th derivative at t(left) of the given spline.

         CALL bsplvb ( t, n11, n11, t(left), left, biatx )
         DO 25 mm=1,m
   25       coef(mm,k,l) = scrtch(1  ,k,mm)
         DO 30 jp1=2,k
 
            CALL bsplvb ( t, jp1, n22, t(left), left, biatx )
            kmj = k+1 - jp1
            DO 30 mm=1,m
               sum = 0.
               DO 28 i=1,jp1
   28             sum = biatx(i)*scrtch(i,kmj,mm) + sum
   30          coef(mm,kmj,l) = sum
   50    CONTINUE

                                        RETURN
      END SUBROUTINE bspp2d



      SUBROUTINE bsplvb ( t, jhigh, index, x, left, biatx )
!--------------------------------------------------------------------------
!calculates the value of all possibly nonzero b-splines at  x  of order
!
!               jout  =  max( jhigh , (j+1)*(index-1) )
!
!  with knot sequence  t .
!
!******  i n p u t  ******
!  t.....knot sequence, of length  left + jout  , assumed to be nonde-
!        creasing.  a s s u m p t i o n . . . .
!                       t(left)  .lt.  t(left + 1)   .
!   d i v i s i o n  b y  z e r o  will result if  t(left) = t(left+1)
!  jhigh,
!  index.....integers which determine the order  jout = max(jhigh,
!        (j+1)*(index-1))  of the b-splines whose values at  x  are to
!        be returned.  index  is used to avoid recalculations when seve-
!        ral columns of the triangular array of b-spline values are nee-
!        ded (e.g., in  bvalue  or in  bsplvd ). precisely,
!                     if  index = 1 ,
!        the calculation starts from scratch and the entire triangular
!        array of b-spline values of orders 1,2,...,jhigh  is generated
!        order by order , i.e., column by column .
!                     if  index = 2 ,
!        only the b-spline values of order  j+1, j+2, ..., jout  are ge-
!        nerated, the assumption being that  biatx , j , deltal , deltar
!        are, on entry, as they were on exit at the previous call.
!           in particular, if  jhigh = 0, then  jout = j+1, i.e., just
!        the next column of b-spline values is generated.
!
!  w a r n i n g . . .  the restriction   jout .le. jmax (= 20)  is im-
!        posed arbitrarily by the dimension statement for  deltal  and
!        deltar  below, but is  n o w h e r e  c h e c k e d  for .
!
!  x.....the point at which the b-splines are to be evaluated.
!  left.....an INTEGER*4 chosen (usually) so that
!                  t(left) .le. x .le. t(left+1)  .
!
!******  o u t p u t  ******
!  biatx.....array of length  jout , with  biatx(i)  containing the val-
!        ue at  x  of the polynomial of order  jout  which agrees with
!        the b-spline  b(left-jout+i,jout,t)  on the interval (t(left),
!        t(left+1)) .
!
!******  m e t h o d  ******
!  the recurrence relation
!
!                       x - t(i)              t(i+j+1) - x
!     b(i,j+1)(x)  =  -----------b(i,j)(x) + ---------------b(i+1,j)(x)
!                     t(i+j)-t(i)            t(i+j+1)-t(i+1)
!
!  is used (repeatedly) to generate the (j+1)-vector  b(left-j,j+1)(x),
!  ...,b(left,j+1)(x)  from the j-vector  b(left-j+1,j)(x),...,
!  b(left,j)(x), storing the new values in  biatx  over the old. the
!  facts that
!            b(i,1) = 1  if  t(i) .le. x .lt. t(i+1)
!  and that
!            b(i,j)(x) = 0  unless  t(i) .le. x .lt. t(i+j)
!  are used. the particular organization of the calculations follows al-
!  gorithm  (8)  in chapter x of the text.
!--------------------------------------------------------------------------
!      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)

      IMPLICIT NONE
      INTEGER(I4B),PARAMETER :: jmax = 4
      INTEGER(I4B) index,jhigh,left,   i,j,jp1
      REAL(DP) x,saved,term

      REAL(DP) deltal(jmax),deltar(jmax)
      REAL(DP) biatx(jhigh), t(left+jhigh)
!urrent fortran standard makes it impossible to specify the length of
!  t  and of  biatx  precisely without the introduction of otherwise
!  superfluous additional arguments.
      DATA j/1/
      SAVE j,deltal,deltar  ! (valid in fortran 77)
!
                                        go to (10,20), index
   10 j = 1
      biatx(1) = 1.
      IF (j .GE. jhigh)                 go to 99
!
   20    jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0.
         DO 26 i=1,j
            term = biatx(i)/(deltar(i) + deltal(jp1-i))
            biatx(i) = saved + deltar(i)*term
   26       saved = deltal(jp1-i)*term
         biatx(jp1) = saved
         j = jp1
         IF (j .LT. jhigh)              go to 20
!
   99                                   RETURN
      END SUBROUTINE bsplvb

      SUBROUTINE banfac ( a, nrow, n, ndiag, middle, iflag )
!------------------------------------------------------------
!--
!-------------------------------------------------------------
!      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      IMPLICIT NONE
      INTEGER(i4B)nrow,ndiag,ilo,middle,iflag,i,ihi,jmax,kmax,&
              mmj,n,j,k
      REAL(DP)  a(nrow,ndiag),diag
      iflag = 1
      ilo = middle - 1
      IF (ilo)                          999,10,19
   10 DO 11 i=1,n
         IF(a(i,1) .EQ. 0.)             go to 999
   11    CONTINUE
                                        RETURN
   19 ihi = ndiag - middle
      IF (ihi)                          999,20,29
   20 DO 25 i=1,n
         IF (a(i,middle) .EQ. 0.)       go to 999
         jmax = MIN(ilo,n-i)
         IF (jmax .LT. 1)               go to 25
         DO 23 j=1,jmax
   23       a(i+j,middle-j) = a(i+j,middle-j)/a(i,middle)
   25    CONTINUE
                                        RETURN
   29 DO 50 i=1,n
         diag = a(i,middle)
         IF (diag .EQ. 0.)              go to 999
         jmax = MIN(ilo,n-i)
         IF(jmax .LT. 1)                go to 50
         kmax = MIN(ihi,n-i)
         DO 33 j=1,jmax
            mmj = middle-j
            a(i+j,mmj) = a(i+j,mmj)/diag
            DO 33 k=1,kmax
   33          a(i+j,mmj+k) = a(i+j,mmj+k) - a(i+j,mmj)*a(i,middle+k)
   50    CONTINUE
                                        RETURN
  999 iflag = 2
                                        RETURN
      END       SUBROUTINE banfac 



      SUBROUTINE banslv ( a, nrow, n, ndiag, middle, b )
!----------------------------------------------------------------------------
! --
!----------------------------------------------------------------------------


!      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)

       IMPLICIT NONE
       INTEGER(I4B) n,nrow,ndiag,middle,ilo,jmax,ihi,i,j
       REAL(DP)  a(nrow,ndiag),b(n)
      IF (n .EQ. 1)                     go to 21
      ilo = middle - 1
      IF (ilo .LT. 1)                   go to 21
      DO 19 i=2,n
         jmax = MIN(i-1,ilo)
         DO 19 j=1,jmax
   19       b(i) = b(i) - b(i-j)*a(i,middle-j)
!
   21 ihi = ndiag-middle
      DO 30 i=n,1,-1
         jmax = MIN(n-i,ihi)
         IF (jmax .LT. 1)               go to 30
         DO 25 j=1,jmax
   25       b(i) = b(i) - b(i+j)*a(i,middle+j)
   30    b(i) = b(i)/a(i,middle)
                                        RETURN
      END SUBROUTINE banslv



   END MODULE bicubic_spline
