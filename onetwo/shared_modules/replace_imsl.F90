  MODULE replace_imsl
!-----------------------------------------------------------------------------------
! original imsl routines called from outside this module are prefixed with 'my_***'
! routiens called internally in this module use origianl names
!------------------------------------------------------------HSJ--------------------

     USE nrtype,                               ONLY : DP,I4B

     USE bicube,                               ONLY : n2cspln,nh2,nwork

     USE cubic_spline,                         ONLY : eval_1d_cubic_spline, &
                                                      create_1d_spline_object,&
                                                      done_1dspline, spl1d_inst

     USE error_handler,                        ONLY : lerrno,terminate

     USE io_gcnmp,                             ONLY : ncrt,nlog

     INTEGER(I4B) ioerr
     TYPE(spl1d_inst)spline_objct
     DATA ioerr /-1/
#ifdef TESTING
     include 'imsl.i'
#endif
     CONTAINS

       SUBROUTINE my_usmnmx(array,maxi,stride,min_val,max_val)
       !-----------------------------------------------------------------------------------
       ! find and return the minimum and maximum values
       ! this routine assumes that maxi may be less than actual size of array
       ! the starting value,sti, is given by lower  bound of array
       !-------------------------------------------------------------------HSJ----6/7/2012-
         REAL(DP),DIMENSION(:) ::  array
         REAL(DP) min_val,max_val,min_valt,max_valt
         INTEGER(I4B) j,maxi,stride
         INTEGER(I4B) sti,eti

         sti = LBOUND(array,DIM=1)
         eti = UBOUND(array,DIM=1)


         IF(maxi == eti .AND. stride ==1)THEN
           max_val = MAXVAL(array)
           min_val = MINVAL(array)
         ELSE
            max_val = -HUGE(1._DP)
            min_val =  HUGE(1._DP)
            DO j = sti,maxi,stride
               max_val = MAX(max_val,array(j))
               min_val = MIN(min_val,array(j))
            ENDDO
         ENDIF 
#ifdef TESTING
            CALL usmnmx(array,maxi,stride,min_valt,max_valt)
            IF(ABS(min_valt -min_val) .GT. 1.e-14 .OR. &
                  ABS(max_valt -max_val) .GT. 1.e-14)THEN
               IF(ioerr == -1)CALL set_err_file
               WRITE(ioerr,FMT='("error,my_usmnmx,values =",2(1pe14.6))')min_val,max_val
               WRITE(ioerr,FMT='("error, usmnmx values =",2(1pe14.6))')min_valt,max_valt
               STOP 'error in my_usmnmx,replace_imsl.f90'

            ENDIF
#endif
         RETURN
       END SUBROUTINE my_usmnmx

       SUBROUTINE my_vsrta(array,maxi)
       !-------------------------------------------------------------
       !-- algebraic sort of array, elements from 
       !-------------------------------------------------------HSJ---
         USE m_refsor,                            ONLY : refsor,D_refsor
                                                          
         REAL(DP),DIMENSION(:) ::  array
         INTEGER(I4B) j,maxi,sti,eti
         REAL(DP),DIMENSION(:), ALLOCATABLE :: cpa
         REAL(DP) maxer
         sti = LBOUND(array,DIM=1)
         eti = UBOUND(array,DIM=1)


        CALL D_refsor (array)          ! refsor is overloaded name
                                       ! NO! - pgf90 complains about overload

#ifdef TESTING
            ALLOCATE(cpa(sti:eti))
            cpa(:) = array(:)
            CALL vsrta(cpa,maxi)

           maxer = 0.0_DP
           DO j=1,1,maxi
              maxer = MAX(ABS(cpa(j)-array(j)),maxer)
           ENDDO
           IF(maxer .GT. 1.e-14_DP)THEN
                 IF(ioerr == -1)CALL set_err_file
                 WRITE(ioerr,FMT='("error in my_vsrta,replace_imsl.f90",4(1pe14.6))')max_er
                 STOP 'error in my_vsrta,replace_imsl.f90'
           ENDIF 
           DEALLOCATE(cpa)
#endif          
       END SUBROUTINE my_vsrta
      
       SUBROUTINE my_ibcccu (a2d,rgrid,nr,zgrid,nz,cspln,nw,wnoperm,ier)
       !-------------------------------------------------------------
       !-- 2d spline interpolation - setup
       !-- note that the array a2dt must be of size (nr,nz)
       !-- for the routines in bicubic spline to work properly.
       !-- When called by onetwo it may be that dim a2d is not (nr,nz)
       !-- (but imsl ibcccu  requires that nw is in fact first dimension
       !-- of a2d(nw,*).
       !-- We deal with this by making a copy of a2d with the
       !-- required size,a2dt(nr,nz) and using it in the constructor call
       !-------------------------------------------------------HSJ---

        USE bicubic_spline,                        ONLY : spline_2d_interp_setup, &
                                                          done_2dspline

        IMPLICIT NONE 
        REAL(DP) rgrid(nr),zgrid(nz)
        REAL(DP) a2d(nw,nh2/2)
        REAL(DP),DIMENSION(:,:),ALLOCATABLE :: a2dt
        REAL(DP)  cspln(n2cspln,nw,nh2),wnoperm(nwork)
        INTEGER(I4B) icalc,nr,nz,nw,ier

        ! clear previous instance of 2d spline (non fatal):
         CALL done_2dspline
        ! create  array cs  in bicubic_spline module back end:
        ! spline routine assumes working size of a2d is  same
        ! as dimensioned size
         ALLOCATE(a2dt(nr,nz))
         a2dt(1:nr,1:nz) = a2d(1:nr,1:nz)
        CALL spline_2d_interp_setup(a2dt,rgrid,zgrid) 
        DEALLOCATE(a2dt)
#ifdef TESTING
        ! EXPLICIT RETURN OF CS == CSPLN
        ! testing is done in  my_dbcevl1
        CALL ibcccu (a2d,rgrid,nr,zgrid,nz,cspln,nw,wnoperm,ier)
#endif   

       END SUBROUTINE my_ibcccu
       




       SUBROUTINE my_dbcevl1 (rgrid,nr,zgrid,nz,cspln,nw,reval,zeval, &
                              pds,ier,icalc)
       !-------------------------------------------------------------
       ! -- 2d spline interpolation - evaluate
       ! -- note that the spline coefficient array cs must be set up
       ! -- allready in the back  end of  bicubic_spline
       ! -- this is accomplished by first calling my_ibcccu
       !-------------------------------------------------------HSJ---
        USE bicubic_spline,                        ONLY : eval_2d_spline

        IMPLICIT NONE

        REAL(DP) pds(6),rgrid(nr),zgrid(nz)
        REAL(DP) cspln(n2cspln,nw,nh2),pdst(6)
        REAL(DP) reval,zeval,max_er
        INTEGER(I4B) icalc,ier,nr,nz,nw,j
        ier = 0

          CALL  eval_2d_spline(reval,zeval,pds,icalc)

#ifdef TESTING
        ! EXPLICIT RETURN OF CS == CSPLN 
        CALL dbcevl1 (rgrid,nr,zgrid,nz,cspln,nw,reval,zeval, &
                      pdst,ier,icalc)
        max_er = 0.0_DP
        DO j=1,icalc
           max_er = MAX(max_er,ABS(pds(j)-pdst(j)))
        ENDDO



        IF(max_er .GT. 1.e-5)THEN
            IF(ioerr == -1)CALL set_err_file
            WRITE(ioerr,FMT='("message from",a)')imslmd
            WRITE(ioerr,FMT='("2d interpolation,my_dbcevl1, not good, max_er =",4(1pe14.6))')max_er
            DO j=1,icalc
               WRITE(ioerr,FMT='("2d interpolation,my_dbcevl1,pds,pdst =",2(1pe14.6))')pds(j),pdst(j)
            ENDDO
!call stop('testing',1)
        ENDIF
        pds(:) = pdst(:)
#endif   

       END SUBROUTINE my_dbcevl1





       SUBROUTINE  my_icsicu(xknot,yknot,nknots,bpar,cs,nfitmax,ier)
       !-------------------------------------------------------------
       ! -- 1d spline interpolation - setup
       ! -- note that the spline coefficient array cs must be set up
       ! -- allready in the back  end of  cubic_spline
       ! -- this is accomplished by first calling my_ibcccu
       !-------------------------------------------------------HSJ---
  

        IMPLICIT NONE
        INTEGER(I4B) nfitmax,ier,nknots,ib
        REAL(DP) xknot(nknots),yknot(nknots),bpar(4)
        REAL(DP) cs(nfitmax,3)
        REAL(DP),PARAMETER :: tol = 1.e-12
        REAL(DP) dx,dy,ax,bx


        CALL done_1dspline(spline_objct)

        !---------------------------------------------------------------------
        ! examine bpar to set boundary conditions:
        ! bpar(1) = 1.0   ==> first deriv at xknot(1) is specified
        ! bpar(3) = 1.0   ==>                xknot(nknots)
        ! second derivative
        ! bpar(1:4) = 0.0 ==> natural cubic spline
        ! more generaly Y(1)'',y(nknots)'' need not be zero
        ! bpar(1) = 0.0 and bpar(2) = 2.*y(1)''
        ! bpar(3) = 0.0 and bpar(4) = 
        !---------------------------------------------------------------------
        ib = -1
        IF(ABS(bpar(1) -1._DP) .LT. tol)THEN
           ! first deriv is given at left end
              dx = xknot(2)-xknot(1) ; dy = yknot(2) - yknot(1)
              ax = dy/dx - bpar(2)*dx/6._DP
              IF(ABS(ax) .LT. tol)ax= 0.0_DP
           IF(ABS(bpar(3) -1._DP) .LT. tol)THEN
              ! first drivativwe both ends
              ib = 1 
              dx = xknot(nknots)-xknot(nknots-1) 
              dy = yknot(nknots) - yknot(nknots-1)
              bx = dy/dx -  bpar(4)*dx/6._DP
           ELSEIF(ABS(bpar(3)) .LT. tol)THEN
              ! second deriv at right end is given
              ib = 6
              bx = bpar(4)/2._DP
           ELSE
              IF(ioerr == -1)CALL set_err_file

              WRITE(ioerr,FMT='("my_icsicu, bpar =",4(1pe14.6))')bpar
              STOP "my_icsicu ,bc condition erro"
           ENDIF
        ELSEIF(ABS(bpar(1)) .LT. tol)THEN
              ! second deriv at left end
              ax = bpar(2)/2._DP
           IF(ABS(bpar(3))  .LT. tol)THEN
              ! second derivative at right end
              bx = bpar(4)/2._DP
              ib = 2
           ELSEIF(ABS(bpar(3) -1._DP) .LT. tol)THEN
              ! first deriv at right end
              IF(ioerr == -1)CALL set_err_file
              WRITE(ioerr,FMT='("my_icsicu,bc prob, bpar =",4(1pe14.6))')bpar
              STOP "my_icsicu ,bc condition erro"
              
           ENDIF
        ELSE 
           ! here bpar(1) is not 0.0 nor 1.0
           ! This implies that second derivative is given at xknot(1),xknot(2)
           ! this option is not implemented
           IF(ioerr == -1)CALL set_err_file
           WRITE(ioerr,FMT='("my_icsicu,bc prob, bpar =",4(1pe14.6))')bpar
           WRITE(ioerr,FMT='("bpar(1) =0.0 or bpar(1) =1.0 are implemented")')
           STOP "my_icsicu ,bc condition erro"
        ENDIF
        
        IF(ib == -1)THEN
           IF(ioerr == -1)CALL set_err_file
           WRITE(ioerr,FMT='("my_icsicu, ib = -1, bpar =",4(1pe14.6))')bpar
           STOP "my_icsicu ,bc condition erro"
        ENDIF

        CALL create_1d_spline_object(spline_objct,nknots,xknot,yknot,ax,bx,ib)


#ifdef TESTING
        CALL icsicu(xknot,yknot,nknots,bpar,cs,nfitmax,ier)
        ! testing is done  in sub my_icsevu
#endif   

       END SUBROUTINE  my_icsicu
      

       SUBROUTINE  my_icsccu (xknot,yknot,nknot,cs,kj,ier)
       !-------------------------------------------------------------
       ! -- 1d spline interpolation , not a knot bc conditions
       !-------------------------------------------------------------
        IMPLICIT NONE
        INTEGER(I4B) ier,nknot,ib,kj
        REAL(DP) xknot(nknot),yknot(nknot)
        REAL(DP) cs(kj,3)
        REAL(DP),PARAMETER :: tol = 1.e-12
        REAL(DP) ax,bx

        ib = 5
        CALL create_1d_spline_object(spline_objct,nknot,xknot,yknot,ax,bx,ib)

#ifdef TESTING
        CALL icsccu (xknot,yknot,nknot,cs,kj,ier)
#endif   
       END SUBROUTINE  my_icsccu 

       SUBROUTINE  my_icsevu(xknot,yknot,nknots,cs,nfitmax, &
                      xin,yout,npts,ier)
       !-------------------------------------------------------------
       ! -- 1d spline interpolation - evaluate
       ! -- note that the spline coefficient array cs must be set up
       ! -- allready in the back  end of  cubic_spline
       ! -- this is accomplished by first calling my_ibcccu
       !-------------------------------------------------------HSJ---
 
        IMPLICIT NONE

        INTEGER(I4B) nfitmax,ier,npts,nknots,j


        REAL(DP) xknot(nknots),yknot(nknots),xin(npts),yout(npts),youtt(npts)
        REAL(DP) cs(nfitmax,3)
        REAL(DP) max_er,tol
        DO j=1,npts
           spline_objct%xx = xin(j)
           CALL eval_1d_cubic_spline(spline_objct)
           yout(j) = spline_objct%sp
        ENDDO

#ifdef TESTING
       ier =0
       CALL icsevu(xknot,yknot,nknots,cs,nfitmax, &
                      xin,youtt,npts,ier)

          IF(ier == 33)THEN
             WRITE(6,FMT='("xin(1),xknot(1) =",2(x,1pe16.8))')xin(1),xknot(1)
          ELSEIF(ier == 34 )THEN
             WRITE(6,FMT='("xin(npts),xknot(nknots) =",2(x,1pe16.8))')xin(npts),xknot(nknots)
          ENDIF
       max_er = 0.0_DP ; tol = 1.e-6
       DO j=1,npts
          max_er = MAX(max_er,ABS(youtt(j)-yout(j)))
       ENDDO
       IF(max_er .GT. tol)THEN
          IF(ioerr == -1)CALL set_err_file
          WRITE(ioerr,FMT='("1d interpolation,my_icsevu, not good, max_er =",1pe14.6)')max_er
       ENDIF

#endif

        END SUBROUTINE  my_icsevu

      SUBROUTINE my_icsmou (x, y, nx, dis, sc, maxit, wk, ierr)
!*******************************************************************************
!*                                                                             *
!*    ROUTINE NAME:  ICSMOU                                                    *
!*                                                                             *
!*    PURPOSE:  one-dimensional data smoothing by error detection              *
!*                                                                             *
!*    AUTHOR:   Steve Reid - 20 January 2000                                   *
!*              S & E Computing                                                *
!*                                                                             *
!*    INPUTS:   x     - vector on nx ascending abcissae                        *
!*              y     - vector of nx ordinates; points passing the stopping    *
!*                      criterion (see sc below) are not changed.              *
!*              nx    - number of elements in x and y (must be > 7)            *
!*              dis   - proportion of the distance the ordinate in error is    *
!*                      moved toward the interpolating curve (a number in the  *
!*                      range [0,1], usually input as 1.0)                     *
!*              sc    - stopping criterion - If the absolute distance between  *
!*                      the ordinate and the interpolating curve is less than  *
!*                      or equal to sc * ((x(i + 3) - x(i - 3)) / 6.0), then   *
!*                      the ordinate is not changed for y(i), where i = 3,...  *
!*                      , nx - 3.  Typically, sc is input as 0.0.              *
!*              maxit - maximum number of iterations allowed.  A suggested     *
!*                      input value is the estimated number of ordinates in    *
!*                      error.                                                 *
!*              wk    - work vector having 2 * nx - 12 elements                *
!*                                                                             *
!*    OUTPUTS:  y     - vector of nx smoothed ordinates                        *
!*              maxit - actual number of iterations taken                      *
!*              ierr  - error code (terminal error)                            *
!*                       129 = nx is less than 7                               *
!*                       130 = input abcissae are not in ascending order       *
!*                                                                             *
!*    ROUTINES CALLED:  icsevu, spline, uertst                                 *
!*                                                                             *
!*    LIMITATIONS:  In general, the routine will not work well if more than    *
!*                  25% of the data points are in error.                       *
!*                                                                             *
!*******************************************************************************
      IMPLICIT NONE
      INTEGER(I4B) i, ierr, isp, istart, istop, itlim, j, maxit, nx, nxm3, &
              nxm6, nxm9
      REAL(DP) c(6, 3), const, dis, sc, t(1), wk(2 * nx - 12), x(nx), &
             xpoint(6), y(nx), ydev, ydevmax, ypoint(6)
      CHARACTER*6 isub /'icsmou'/
!
!--Validate the input parameters.
!
      IF (nx .LT. 7) THEN
         ierr = 129
         go to 9000
      ELSE  ! ensure that abcissae are in ascending order
         DO i = 2, nx
            IF (x(i) .LE. x(i - 1)) THEN
               ierr = 130
               go to 9000
            ENDIF
         ENDDO
      ENDIF
      ierr = 0
      nxm3 = nx - 3
      nxm6 = nx - 6
      nxm9 = nx - 9
!
!--Store the stopping criterion for each of the nx-6 interior points.
!
      const = sc / 6.0d0
      DO i = 1, nxm6
         wk(i) = const * (x(i + 6) - x(i))  ! save in work vector
      ENDDO
!
!--Store the deviations between the input functional values and the computed
!--functional values for the nx-6 interior points right after the stopping
!--criteria in the work vector.
!
      itlim  = maxit  ! save iterations limit input by user
      maxit  = 0      ! initialize counter for actual iterations (output
      istart = 4      ! initialize start index
      istop  = nxm3   ! initialize stop index
      isp    = 99     ! spine pointer set to arbitrary value to enter lo
      DO WHILE (isp .GT. 0)  ! there is still a spline point out of toler
         DO i = istart, istop
!
!--Find a spline fit for the ith point using the three preceeding points and
!--the three following points.
!
            isp = i - 4  ! initialize the spline pointer to four points
            DO j = 1, 6
               IF (j .NE. 4) THEN
                  isp = isp + 1
               ELSE  ! skip over the ith point itself
                  isp = isp + 2  ! set spline index to the next point
               ENDIF
               xpoint(j) = x(isp)
               ypoint(j) = y(isp)
            ENDDO
!
!--Get the spline coefficients for the spline segment.
!
            CALL spline (6, xpoint, ypoint, c(1, 1), c(1, 2), c(1, 3))
!
!--Get the functional value at x(i) using these coefficients.
!
            CALL my_icsevu (xpoint, ypoint, 6, c, 6, x(i), t, 1, ierr)
!
!--Store the deviation between the input functional value and the spline value.
!
            wk(nxm9 + i) = y(i) - t(1)     ! save in work vector
         ENDDO
         IF (maxit .GE. itlim) go to 9990  ! exceeded maximum allowed it
!
!--Find the maximum deviation between input and calculated functional values.
!
         isp     = 0  ! initialize the spline index to zero
         ydevmax = 0.0d0
         DO i = 1, nxm6
            ydev = ABS (wk(nxm6 + i))
            IF (ydev .GT. wk(i) .AND. ydev .GT. ydevmax) THEN
               ydevmax = ydev
               isp     = i + 3  ! spline index to maximum unacceptable d
            ENDIF
         ENDDO
         IF (isp .LE. 0) go to 9990  ! stopping criterion has been met
         maxit = maxit + 1           ! increment the iteration count
!
!--Adjust the worst-case point by the specified percentage (dis) of the
!--displacement and recompute points affected by moving this point.
!
         y(isp) = y(isp) - (wk(nxm9 + isp) * dis)
         istart = MAX (isp - 3, 4)     ! start at third point before, if
         istop  = MIN (isp + 3, nxm3)  ! stop at third point after, if p
      ENDDO
!
!--Error exit.
!
 9000 CONTINUE
      CALL my_uertst (ierr, isub)
!
!--Normal exit.
!
 9990 CONTINUE
      RETURN
      END SUBROUTINE my_icsmou



        SUBROUTINE my_icsvku (v_ind,v_dep,nval,xknot,nknots,yknot,csp,nwt,error, &
                              work_icsvku,ier)
        !--------------------------------------------------------------
        ! -- least squares approximation using cubic splines
        ! -- no implementation at this time
        ! -- could use toms642_cubspl.f (/usc-data/p2/linux/tools/math_libs/splinefit_2/cubspl)
        !----------------------------------------------------------------
 
        IMPLICIT NONE

        INTEGER(I4B) nxk,nwt,ier,nval,nknots
        REAL(DP) v_ind(nval),v_dep(nval),xknot(nknots),yknot(nknots) 
        REAL(DP) error
        REAL(DP) csp(nwt,3),work_icsvku(nval*(nknots+6))
#ifdef TESTING
        CALL icsvku (v_ind,v_dep,nval,xknot,nknots,yknot,csp,nwt,error, &
                    work_icsvku,ier)
          IF(ioerr == -1)CALL set_err_file
          WRITE(ioerr,FMT='("least sq spline fit, icsvku used")')
#endif
        END SUBROUTINE my_icsvku



      SUBROUTINE my_icsscv (x, f, nx, y, sc, isc, ijob, wk, ierr)
!*******************************************************************************
!*                                                                             *
!*    ROUTINE NAME:  ICSSCV                                                    *
!*                                                                             *
!*    PURPOSE:  cubic spline smoother (easy-to-use version)                    *
!*                                                                             *
!*    AUTHOR:   Steve Reid - 11 February 2000                                  *
!*              S & E Computing                                                *
!*                                                                             *
!*    INPUTS:   x     - vector on nx ascending abcissae                        *
!*              f     - vector of nx ordinates                                 *
!*              nx    - number of data points (must be > 3)                    *
!*              isc   - row dimension of spline coefficient matrix (sc)        *
!*                      exactly as dimensioned in the calling program          *
!*              ijob  - job selection parameter                                *
!*                      1 = nx <= 20 or unevenly spaced abcissae are used (not *
!*                          implemented)                                       *
!*                      2 = nx > 20 and abcissae are evenly spaced             *
!*              wk    - work vector dimensioned nx * (3 * nx + 5) for ijob = 1 *
!*                      or 6 * nx for ijob = 2.                                *
!*                                                                             *
!*    OUTPUTS:  y     - vector of nx smoothed ordinates                        *
!*              sc    - spline coefficient matrix dimensioned nx -1 by 3.      *
!*                      The value of the spline approximation at t is          *
!*                      s(t) = ((sc(i, 3) * d + c(i, 2)) * d +c(i, 1)) * d     *
!*                             + y(i), where x(i) .le. t .lt. x(i + 1) and     *
!*                                           d = t - x(i).                     *
!*              ierr  - error code (terminal error)                            *
!*                       129 = ic is less than nx - 1                          *
!*                       130 = nx is less than 4                               *
!*                       131 = input abcissae are not in ascending order       *
!*                       132 = failure to converge on a minimum of the         *
!*                             cross-validation function (not implemented)     *
!*                       133 = ijob is 2 and the abcissae are not evenly       *
!*                             spaced                                          *
!*                                                                             *
!*    ROUTINES CALLED:  dcssmo, uertst                                         *
!*                                                                             *
!*    LIMITATIONS:  The number of arithmetic operations is proportional to nx  *
!*                  when ijob = 2.  If the case for ijob = 1 were implemented, *
!*                  the number of arithmetic operations would be proportional  *
!*                  to nx**3.                                                  *
!*                                                                             *
!*******************************************************************************
      IMPLICIT NONE
      INTEGER(I4B) i, ierr, ijob, isc, nx
      REAL(DP) delmax, delmin, deltax, f(*), rho, sc(isc, 3), step_size, &
             wgs(nx), wk(*), x(*), y(*)
      CHARACTER*6 isub /'icsscv'/
!
!--Validate the input parameters.
!
      IF (isc .LT. nx - 1) THEN
         ierr = 129
         go to 9000
      ELSEIF (nx .LT. 4) THEN
         ierr = 130
         go to 9000
      ELSE  ! ensure that abcissae are in ascending order
         DO i = 2, nx
            IF (x(i) .LE. x(i - 1)) THEN
               ierr = 131
               go to 9000
            ENDIF
         ENDDO
      ENDIF
      IF (ijob .NE. 2) THEN  ! unexpected value of ijob input
         PRINT *,  'only ijob = 2 is used in onetwo'
#ifdef  TESTING
         IF(ioerr == -1)CALL set_err_file
         WRITE(ioerr,FMT='("sub my_icsscv only ijob = 2 is used in onetwo")')
#endif
         STOP 'icsscv called with illegal value for ijob'
      ELSE  ! ijob = 2 is the implementation for onetwo
!
!--Find the maximum and minimum abcissa spacing.
!
         step_size = x(2) - x(1)
         delmax    = step_size
         delmin    = step_size
         wgs(1)    = 1.0
         wgs(2)    = 1.0
         DO i = 3, nx
            deltax = x(i) - x(i - 1)
            delmax = MAX (delmax, deltax)
            delmin = MIN (delmin, deltax)
            wgs(i) = 1.0  ! all points given the same weight
         ENDDO
         IF (delmax - delmin .GE. 0.01 * delmax) THEN
            ierr = 133  ! abcissae not considered evenly spaced
            go to 9000
         ENDIF
         rho = 0.5  ! the user might need to adjust this smoothing
                    ! factor (smaller for increased smoothing)
      ENDIF
      CALL dcssmo (step_size, nx, x, f, wgs, rho, y, sc(1, 1), &
                   sc(1, 2), sc(1, 3))
      go to 9990
!
!--Error exit.
!
 9000 CONTINUE
      CALL my_uertst (ierr, isub)
!
!--Normal exit.
!
 9990 CONTINUE
      RETURN
      END      SUBROUTINE my_icsscv


      SUBROUTINE my_ibcieu (f, ifd, x, nx, y, ny, xl, nxl, yl, nyl, &
                         fl, ifld, wk, ierr)
!*******************************************************************************
!*                                                                             *
!*    ROUTINE NAME:  IBCIEU                                                    *
!*                                                                             *
!*    PURPOSE:  bicubic spline two-dimensional interpolator                    *
!*                                                                             *
!*    AUTHOR:   Steve Reid - 09 September 1999                                 *
!*              S & E Computing                                                *
!*                                                                             *
!*    INPUTS:   f    - functional value at (x(i), y(j)),                       *
!*                                               where i = 1, nx; j = 1, ny    *
!*              ifd  - row dimension of array f                                *
!*              x    - vector of nx ascending elements                         *
!*              nx   - number of elements in vector x                          *
!*              y    - vector of ny ascending elements                         *
!*              ny   - number of elements in vector y                          *
!*              xl   - vector of length nxl                                    *
!*              nxl  - number of elements in vector xl                         *
!*              yl   - vector of length nyl                                    *
!*              nyl  - number of elements in vector yl                         *
!*              ifld - row dimension of array fl                               *
!*              wk   - work vector of length                                   *
!*                                    max (3 * (nx - 1), 3 * (ny - 1) + ny)    *
!*                                                                             *
!*    OUTPUTS:  fl   - bicubic spline values at (xl(i), yl(j)),                *
!*                                             where i = 1, nxl; j = 1, nyl    *
!*              ierr - error code (terminal error)                             *
!*                      129 = ifd < nx                                         *
!*                      130 = ifld < nxl                                       *
!*                      131 = nx < 2                                           *
!*                      132 = ny < 2                                           *
!*                      133 = elements of vector x not in ascending order      *
!*                      134 = elements of vector y not in ascending order      *
!*                   - error code (warning error)                              *
!*                       37 = at least one xl(i) is out of range               *
!*                       38 = at least one yl(i) is out of range               *
!*                                                                             *
!*    ROUTINES CALLED:  icsccu, icsevu, uertst                                 *
!*                                                                             *
!*    LIMITATIONS:  The elements of input vectors x and y are assumed to be    *
!*                  pre-sorted into ascending order.  The input ic must be     *
!*                  greater than or equal to nx.                               *
!*                                                                             *
!*******************************************************************************
      IMPLICIT NONE
      INTEGER(I4B) i, ierr, ifd, ifld, iwarn, j, jwarn, kyl, kylp1, nx, nxl, &
              nxm1, ny, nyl, nym1
      REAL(DP) f(ifd, ny), fl(ifld, *), wk(*), x(nx), xl(nxl), y(ny), &
             yl(nyl)
      CHARACTER*6 isub /'ibcieu'/
!
!--Validate input arguments.
!
      ierr = 0
      IF (ifd .LT. nx) THEN
         ierr = 129
      ELSEIF (ifld .LT. nxl) THEN
         ierr = 130
      ELSEIF (nx .LT. 2) THEN
         ierr = 131
      ELSEIF (ny .LT. 2) THEN
         ierr = 132
      ENDIF
      IF (ierr .NE. 0) go to 9000
!
!--Interpolate in the x-direction.
!
      iwarn = 0
      nxm1  = nx - 1
      nym1  = ny - 1
      DO i = 1, ny
         CALL my_icsccu (x, f(1, i), nx, wk(1), nxm1, ierr)
         IF (ierr .NE. 0) THEN
            ierr = 133
            go to 9000
         ENDIF
         CALL my_icsevu (x, f(1, i), nx, wk(1), nxm1, xl, fl(1, i), nxl, &
                      ierr)
         IF (ierr .NE. 0) iwarn = 37
      ENDDO
!
!--Interpolate in the y-direction.
!
      jwarn = 0
      kyl   = nym1 * 3
      kylp1 = kyl + 1
      DO i = 1, nxl
         DO j = 1, ny
            wk(kyl + j) = fl(i, j)
         ENDDO
!
!--Get the spline coefficients.
!
         CALL my_icsccu (y, wk(kylp1), ny, wk(1), nym1, ierr)
         IF (ierr .NE. 0) THEN
            ierr = 134
            go to 9000
         ENDIF
!
!--Evaluate the cubic spline.
!
         DO j = 1, nyl
            CALL my_icsevu (y, wk(kylp1), ny, wk(1), nym1, yl(j), fl(i, j), &
                         1, ierr)
            IF (ierr .NE. 0) jwarn = 38
         ENDDO
      ENDDO
      ierr = MAX (iwarn, jwarn)
      IF (iwarn .NE. 0) CALL my_uertst (iwarn, isub)
      IF (jwarn .NE. 0) CALL my_uertst (jwarn, isub)
      go to 9990
!
!--Error exit.
!
 9000 CONTINUE
      CALL my_uertst (ierr, isub)
!
!--Normal exit.
!
 9990 CONTINUE
      RETURN
      END       SUBROUTINE my_ibcieu 



      SUBROUTINE icsvku (x, f, nx, xk, nxk, y, sc, isc, error, wk, ierr)
!*******************************************************************************
!*                                                                             *
!*    ROUTINE NAME:  ICSVKU                                                    *
!*                                                                             *
!*    PURPOSE:  least squares approximation by cubic splines - variable knots  *
!*                                                                             *
!*    AUTHOR:   Steve Reid - 25 February 2000                                  *
!*              S & E Computing                                                *
!*                                                                             *
!*    INPUTS:   x     - vector on nx ascending abcissae                        *
!*              f     - vector of nx ordinates                                 *
!*              nx    - number of data points                                  *
!*              xk    - vector of nxk ascending estimated knot locations such  *
!*                      that xk(1) .le. x(1) and xk(nxk) .ge. x(nx)            *
!*              nxk   - number of knots (elements of xk vector)                *
!*              isc   - row dimension of spline coefficient matrix (sc)        *
!*              wk    - work vector dimensioned nx * (nxk + 6)                 *
!*                                                                             *
!*    OUTPUTS:  xk    - vector of nx knot locations as determined by icsvku    *
!*              y     - vector of (nxk - 1) spline approximations defined by   *
!*                      s(t) = ((sc(i, 3) * d + sc(i, 2)) * d + sc(i, 1)) * d  *
!*                      + y(i), where xk(i) .le. t .lt. xk(i + 1) and          *
!*                      d = t - xk(i)                                          *
!*              sc    - spline coefficient matrix dimensioned (nxk - 1) by 3   *
!*              error - least squares error of cubic spline approximation      *
!*              ierr  - error code (warning error)                             *
!*                        33 = convergence criteria was not satisfied, so the  *
!*                             best possible results will be returned          *
!*                      error code (terminal error)                            *
!*                       130 = nxk is greater that 28                          *
!*                       131 = input abcissae are not in ascending order       *
!*                       132 = input abcissae do not lie between the boundary  *
!*                             knots                                           *
!*                       133 = input knots are not in ascending order          *
!*                                                                             *
!*    ROUTINES CALLED:  cubgcv, uertst                                         *
!*                                                                             *
!*    LIMITATIONS:  icsvku generally yields low to medium accuracy (two to     *
!*                  five significant digits).                                  *
!*                                                                             *
!*******************************************************************************
      IMPLICIT NONE
      INTEGER(I4B)  i, ierr, isc, job, nx, nxk
      REAL(DP) df(nx), dummy(nx), error, f(nx), sc(isc, 3), var, &
             wk(nx, 1), x(nx), xk(nxk), y(isc)
      CHARACTER*6 isub /'icsvku'/
!
!--Validate the input parameters.
!
      IF (nxk .GT. 28) THEN
         ierr = 130
         go to 9000
      ELSE  ! ensure that abcissae are in ascending order
         DO i = 2, nx
            IF (x(i) .LT. x(i - 1)) THEN
               ierr = 131
               go to 9000
            ENDIF
         ENDDO
      ENDIF
!
!--Ensure that knots are in ascending order.
!
      DO i = 2, nxk
         IF (xk(i) .LE. xk(i - 1)) THEN
            ierr = 133
            go to 9000
         ENDIF
      ENDDO
!
!--Ensure that all abcissae lie between the boundary knots.
!
      DO i = 1, nx
         IF (x(i) .LE. xk(1) .OR. x(I) .GE. xk(nxk)) THEN
            ierr = 132
            go to 9000
         ENDIF
         df(i) = 1  ! initialize the s.d. of each data point
      ENDDO
      var = 0  ! calculate an interpolating natural cubic spline
      job = 0  ! point standard error estimates are not required
      CALL cubgcv (x, f, df, nx, y, sc, isc, var, job, dummy, wk, ierr)
      IF (ierr .NE. 0) THEN
         PRINT *,  'ierr =', ierr
#ifdef TESTING
         IF(ioerr == -1)CALL set_err_file
         WRITE(ioerr,FMT='("sub icsvku,call to cubgv error,ierr = ",i5)')ierr
#endif
         STOP 'Error occurred in cubgcv'
      ELSE
         ierr = 0
      ENDIF
      go to 9990
!
!--Error exit.
!
 9000 CONTINUE
      CALL my_uertst (ierr, isub)
!
!--Normal exit.
!
 9990 CONTINUE
      RETURN
      END SUBROUTINE icsvku



        SUBROUTINE my_eigrf (a, n,lda, ijob, w, z,iz, wk, info)
        !--------------------------------------------------------------
        ! -- eigenval,eigenvectros, real,general matrix, full storage mode
        ! -- ijob = 0 eigenvalues only returned in w(*)
        ! -- ijob =1 eigenvalues and eigenvectors (returned in w(*) 
        ! --         and z(*,*)
        ! -- ijob = 2 same as 1 plus performance index
        ! -- ijob =3 performance index only
        ! -- pay attention to dim handling of array z. it is correct for
        ! -- what cray321u.f wants but is a rediculous,confusing, fudge to
        ! -- make up for real versus complex arrays. This routine is designed to
        ! -- avoid changing the eigen routine  in onetwo,cray321u.f
        !-------------------------------------------------------HSJ-------
        USE nrtype, ONLY : DP,I4B
        IMPLICIT NONE
        REAL(DP),PARAMETER :: tol = 1.e-12_DP  ! smallest  imaginary part not zero
        INTEGER(I4B)  i,j,k,info,n,lda,ijob,iz,iz2,swk,ldvl,ldvr
        REAL(DP) max_er
        REAL(DP) a(lda,*),wk(*),z(2*iz,*),w(*)
        REAL(DP) wr(n),wi(n)
        REAL(DP),ALLOCATABLE,DIMENSION(:) :: wt,wkk
        REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: vl,vr,eigvecr,eigveci, &
                 ac,zt,bc
        CHARACTER*8 routine

        CHARACTER *1, jobvl,jobvr

        IF(ijob == 0)THEN
           swk = n
        ELSEIF(ijob == 1)THEN
           swk =2*n
        ELSEIF(ijob == 2)THEN
           swk = (2+n)*n
        ELSEIF(ijob == 3)THEN
           swk = 1
        ENDIF
        swk = MAX(swk,4*lda) ! works for both imsl and dgeev
        ALLOCATE(wkk(swk))

#ifdef TESTING ! see matching deallocate 
        ALLOCATE(ac(lda,lda),bc(lda,lda),zt(2*iz,lda))
        ac(1:lda,1:lda) = a(1:lda,1:lda)   ! need copy both routines  destroy a
        zt(:,:) = 0.0_DP
        ALLOCATE(wt(2*n))
        wt(:) = 0.0_DP
#endif
 


        
        jobvr = 'V' ; jobvl = 'N' 
        IF(jobvr == 'V')THEN
           ldvr = n
        ELSE
           ldvr =1
        ENDIF
        ALLOCATE(vr(ldvr,n))

        IF(jobvl == 'V')THEN
           ldvl = n
        ELSE
           ldvl =1
        ENDIF
        ALLOCATE(vl(ldvl,n))

        ALLOCATE(eigvecr(n,n),eigveci(n,n))
        eigvecr(:,:) = 0.0_DP ; eigveci(:,:) = 0.0_DP
        routine = 'dgeev'
        CALL Dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,wkk,swk,info)
        IF(info .NE. 0)&
           STOP 'replace_imsl, sub dgeev, info not 0'
        !eigvec extraction is based on eigenvalue structure

        i =0
        DO WHILE( 1 > 0)
           i = i + 1 
           IF(ABS(wi(i)) .LT. tol)THEN   ! eigenvalue is real
              eigvecr(1:n,i) = vr(1:n,i)
           ELSE  ! eigenvalue is complex
              eigvecr(1:n,i)  = vr(1:n,i)
              eigveci(1:n,i)  = vr(1:n,i+1)
              i = i + 1
              eigvecr(1:n,i)  = vr(1:n,i-1)
              eigveci(1:n,i)  = -vr(1:n,i)
           ENDIF
           IF(i == n) EXIT
        ENDDO


        ! put in imsl  form for compatability


        DO i=1,n  ! over eigenvectors,eigenvalues
           k =0
           DO j=1,2*n-1,2   
              k = k+1
              !Real part,element k, eigenvector i
              z(j,i)   = eigvecr(k,i) 
              !imaginary part,element k, eigenvector i
              z(j+1,i) = eigveci(k,i)
              IF(i ==1)THEN
                 w(j)   = wr(k)
                 w(j+1) = wi(k)
              ENDIF
           ENDDO
        ENDDO



#ifdef TESTING
!      check that bc*z = w*z for dgeev solution:
       bc(:,:) = ac(:,:)  ! 88888889999999
       iz2 = 2*iz ! 88888889999999
       CALL check_sol(bc,lda,w,z,iz2,n,routine) ! 88888889999999

!       Now get solution from imsl routine:
        routine = 'eigrf'
        CALL eigrf (ac, n,lda, ijob, wt, zt,iz, wkk, info)
        IF(info .NE. 0)&
           STOP 'replace_imsl, sub eigrf, info not 0'

 ! note eigenvectors and eigenvalues may not be in same order
 ! on output from dgeev and eigrf
 


        iz2 = 2*iz
        CALL check_sol(bc,lda,wt,zt,iz2,n,routine)

        DEALLOCATE(wt,zt,ac,bc)
#endif
        DEALLOCATE(wkk,vr,vl,eigvecr,eigveci)

        END SUBROUTINE my_eigrf

        SUBROUTINE check_sol(m,lda,lam,z,ldz,n,routine)
        !---------------------------------------------------------------
        ! -- eval m*x=lam*x for check
        !----------------------------------------------------------------

          USE nrtype,                          ONLY : DP,I4B
          IMPLICIT NONE
          REAL(DP), PARAMETER :: tol = 1.e-12
          INTEGER(I4B) i,j,k,l,n,lda,ldz
          REAL(DP) m(lda,*), lam(*),z(ldz,*),rhsr(n),rhsi(n),diffr(n),diffi(n)
          REAL(DP) eigvecr(n,n),eigveci(n,n),wr(n),wi(n)
          REAL(DP) sumr,sumi,max_er
          CHARACTER*8 routine
      
!           DO i = 1,n
!                WRITE(6,FMT='(8(1pe12.2,x))')(z(j,i),z(j+1,i),j=1,2*n-1,2)
!           ENDDO

  
          DO i =1,n 
             k =0
             DO j=1,2*n-1,2
                k=k+1
                IF( i == 1)THEN
                   wr(k) = lam(j)
                   wi(k) = lam(j+1)
                ENDIF
                eigvecr(k,i)   = z(j,i)
                eigveci(k,i)   = z(j+1,i)
             ENDDO
          ENDDO
     max_er = 0.0_DP
     DO l=1,n
        DO i =1,n
           sumr = 0.0_DP ; sumi = 0.0_DP
           DO k =1 ,n
              sumr = sumr + m(i,k)*eigvecr(k,l)
              sumi = sumi + m(i,k)*eigveci(k,l)
           ENDDO
           rhsr(i) = wr(l)* eigvecr(i,l) - wi(l)*eigveci(i,l)
           rhsi(i) = wr(l)*eigveci(i,l) +wi(l)* eigvecr(i,l)
           diffr(i) = sumr - rhsr(i)
           diffi(i) =  sumi -rhsi(i)
           max_er = MAX(max_er,ABS(diffi(i)),ABS(diffr(i)))
        ENDDO
     ENDDO


       IF(max_er .GT. tol)THEN
          IF(ioerr == -1)CALL set_err_file
          WRITE(ioerr,FMT='("eigrf problem with eigenvectors,max_er = ",1pe12.4,x,a,x,i5)')max_er,routine,n
        ENDIF


        END SUBROUTINE check_sol




        SUBROUTINE my_lginf (desgnm,maxobsrv,neqn,nfcoil1,tol,fcoilinv, &
                 nfcoil,singval,wdum,ier)
        !--------------------------------------------------------------
        ! --  get the generalized inverse of the (neqn,nfcoil) matrix using
        ! singular value decomposition
        ! -- no implementation at this time
        ! -- 
        !----------------------------------------------------------------
         IMPLICIT NONE
         REAL(DP) desgnm(1),tol,fcoilinv(1),singval(1),wdum(1)
         INTEGER(I4B)nfcoil,ier,maxobsrv,neqn,nfcoil1
#ifdef TESTING
          !CALL lginf (desgnm,maxobsrv,neqn,nfcoil,tol,fcoilinv, &
          !        nfcoil,singval,wdum,ier)
          IF(ioerr == -1)CALL set_err_file
          WRITE(ioerr,FMT='("sub lginf,no imsl replacement")')
          STOP 'my_lginf,imsl replacment not available'
#endif
        END SUBROUTINE my_lginf


     SUBROUTINE spline (n, x, y, b, c, d)
      INTEGER n
      REAL(DP) x(n), y(n), b(n), c(n), d(n)
!
!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
!
      INTEGER nm1, ib, i
      REAL(DP) t
!
      nm1 = n-1
      IF ( n .LT. 2 ) RETURN
      IF ( n .LT. 3 ) go to 50
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      DO 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 CONTINUE
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      IF ( n .EQ. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
!
!  forward elimination
!
   15 DO 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 CONTINUE
!
!  back substitution
!
      c(n) = c(n)/b(n)
      DO 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 CONTINUE
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      DO 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 CONTINUE
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      RETURN
!
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      RETURN
      END      SUBROUTINE spline

       SUBROUTINE set_err_file
       !---------------------------------------------------------------------
       ! --
       !---------------------------------------------------------------------
         ! ioerr unit not defined 
         !CALL  getioun(ioerr,499)
         ioerr = 499
         OPEN (unit = ioerr, file = "rpl_imsl_err.txt", status = 'UNKNOWN')

       END  SUBROUTINE set_err_file




      SUBROUTINE my_mmbsir (arg, order, n, iopt, b, ierr)
! ******************************************************************************
!                                                                              *
!     ROUTINE NAME:  MMBSIR                                                    *
!                                                                              *
!     PURPOSE:  computes the modified Bessel FUNCTION of the First Kind of     *
!               nonnegative REAL order for REAL positive arguments WITH an     *
!               exponential scaling FUNCTION                                   *
!                                                                              *
!     AUTHOR:   Steve Reid - 12 November 1999                                  *
!               S & E Computing                                                *
!                                                                              *
!     INPUTS:   arg   - must be nonnegative and IF no scaling is specified,    *
!                       must not exceed the largest value allowed by the       *
!                       FORTRAN exp FUNCTION.                                  *
!               order - specifies the order of the Bessel FUNCTION (must be    *
!                       positive, but less than one)                           *
!               n     - number of FUNCTION values to be computed               *
!               iopt  - scaling option switch                                  *
!                       1 = calculated Bessel Functions are unscaled           *
!                       2 = calculated Bessel Functions are scaled by          *
!                           EXP (-arg)                                         *
!                                                                              *
!     OUTPUTS:  b     - vector of length n containing the computed FUNCTION    *
!                       values.  b(1) will contain the FUNCTION value for the  *
!                       input order, b(2) will contain the FUNCTION value for  *
!                       order + 1, b(3) for order + 2, etc.                    *
!               ierr  - error code (terminal error)                            *
!                        129 = any one or more of the input arguments is out   *
!                              range and all elements of b are set to machine  *
!                              infinity                                        *
!                        129 + j = the specified n was so much greater than    *
!                                  the largest INTEGER in arg that ONLY the    *
!                                  first j elements of b were calculated       *
!                                  correctly                                   *
!                                                                              *
!     ROUTINES CALLED:  ribesl, my_uertst                                         *
!                                                                              *
!     LIMITATIONS:  NONE                                                       *
!                                                                              *
! ******************************************************************************
      IMPLICIT NONE 
      INTEGER(I4B) ierr, iopt, n, num
      REAL(DP) arg, order, b(*), infinity
      CHARACTER*6 isub /'mmbsir'/
!      DATA infinity /+1.79769313486231e+308/  ! largest real*8 for HP
      DATA infinity /8.0d+307/                !approximate for pgf90  HSJ

      ierr = 0

!--Decompose the matrix A.
!
      CALL ribesl (arg, order, n, iopt, b, ierr)
      IF (ierr .LT. 0) THEN
         ierr = 129    ! at least one of the input arguments is out of range
         num  = MAX (1, n)
         DO WHILE (num .GT. 0)
            b(num) = infinity
            num    = num - 1
         ENDDO
      ELSEIF (ierr .LT. n) THEN
         ierr = 129 + ierr  ! only the first ierr terms were calculated correctly
      ELSE             ! ierr = n and there are no errors
         ierr = 0
      ENDIF
      IF (ierr .EQ. 0) go to 9990
!
!--Error EXIT.
!
      CALL my_uertst (ierr, isub)
!
!--Normal EXIT.
!
 9990 CONTINUE
      RETURN
      END SUBROUTINE my_mmbsir

      SUBROUTINE RIBESL(X, ALPHA, NB, IZE, B, NCALC)                
!-------------------------------------------------------------------
!
!  THIS ROUTINE CALCULATES BESSEL FUNCTIONS I SUB(N+ALPHA) (X)
!  FOR NON-NEGATIVE ARGUMENT X, AND NON-NEGATIVE ORDER N+ALPHA,
!  WITH OR WITHOUT EXPONENTIAL SCALING.
!
!
! EXPLANATION OF VARIABLES IN THE CALLING SEQUENCE
!
! X     - WORKING PRECISION NON-NEGATIVE REAL ARGUMENT FOR WHICH
!         I'S OR EXPONENTIALLY SCALED I'S (I*EXP(-X))
!         ARE TO BE CALCULATED.  IF I'S ARE TO BE CALCULATED,
!         X MUST BE LESS THAN EXPARG (SEE BELOW).
! ALPHA - WORKING PRECISION FRACTIONAL PART OF ORDER FOR WHICH
!         I'S OR EXPONENTIALLY SCALED I'S (I*EXP(-X)) ARE
!         TO BE CALCULATED.  0 .LE. ALPHA .LT. 1.0.
! NB    - INTEGER NUMBER OF FUNCTIONS TO BE CALCULATED, NB .GT. 0.
!         THE FIRST FUNCTION CALCULATED IS OF ORDER ALPHA, AND THE
!         LAST IS OF ORDER (NB - 1 + ALPHA).
! IZE   - INTEGER TYPE.  IZE = 1 IF UNSCALED I'S ARE TO CALCULATED,
!         AND 2 IF EXPONENTIALLY SCALED I'S ARE TO BE CALCULATED.
! B     - WORKING PRECISION OUTPUT VECTOR OF LENGTH NB.  IF THE ROUTINE
!         TERMINATES NORMALLY (NCALC=NB), THE VECTOR B CONTAINS THE
!         FUNCTIONS I/ALPHA/(X) THROUGH I/NB-1+ALPHA/(X), OR THE
!         CORRESPONDING EXPONENTIALLY SCALED FUNCTIONS.
! NCAL! - INTEGER OUTPUT VARIABLE INDICATING POSSIBLE ERRORS.
!         BEFORE USING THE VECTOR B, THE USER SHOULD CHECK THAT
!         NCALC=NB, I.E., ALL ORDERS HAVE BEEN CALCULATED TO
!         THE DESIRED ACCURACY.  SEE ERROR RETURNS BELOW.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! NSIG   - DECIMAL SIGNIFICANCE DESIRED.  SHOULD BE SET TO
!          IFIX(ALOG10(2)*NBIT+1), WHERE NBIT IS THE NUMBER OF
!          BITS IN THE MANTISSA OF A WORKING PRECISION VARIABLE.
!          SETTING NSIG LOWER WILL RESULT IN DECREASED ACCURACY
!          WHILE SETTING NSIG HIGHER WILL INCREASE CPU TIME
!          WITHOUT INCREASING ACCURACY.  THE TRUNCATION ERROR
!          IS LIMITED TO A RELATIVE ERROR OF T=.5*10**(-NSIG).
! ENTEN  - 10.0 ** K, WHERE K IS THE LARGEST INTEGER SUCH THAT
!          ENTEN IS MACHINE-REPRESENTABLE IN WORKING PRECISION.
! ENSIG  - 10.0 ** NSIG.
! RTNSIG - 10.0 ** (-K) FOR THE SMALLEST INTEGER K SUCH THAT
!          K .GE. NSIG/4.
! ENMTEN - THE SMALLEST ABS(X) SUCH THAT X/4 DOES NOT UNDERFLOW.
! XLARGE - UPPER LIMIT ON THE MAGNITUDE OF X WHEN IZE=2.  BEAR
!          IN MIND THAT IF ABS(X)=N, THEN AT LEAST N ITERATIONS
!          OF THE BACKWARD RECURSION WILL BE EXECUTED.
! EXPARG - LARGEST WORKING PRECISION ARGUMENT THAT THE LIBRARY
!          EXP ROUTINE CAN HANDLE AND UPPER LIMIT ON THE
!          MAGNITUDE OF X WHEN IZE=1.
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!         IBM/195    CDC/7600  UNIVAC/1108   VAX 11/780 (UNIX)
!          (D.P.)  (S.P.,RNDG)    (D.P.)     (S.P.)     (D.P.)
!
! NSIG      16         14          18           8        17
! ENTEN   1.0D75     1.0E322     1.0D307     1.0E38    1.0D38
! ENSIG   1.0D16     1.0E14      1.0D18      1.0E8     1.0D17
! RTNSIG  1.0D-4     1.0E-4      1.0D-5      1.0E-2    1.0D-4
! ENMTEN  2.2D-78    1.0E-290    1.2D-308    1.2E-37   1.2D-37
! XLARGE  1.0D4      1.0E4       1.0D4       1.0E4     1.0D4
! EXPARG  174.0D0    740.0E0     709.0D0     88.0E0    88.0D0
!
!*******************************************************************
!*******************************************************************
!
!
! ERROR RETURNS
!
!  IN CASE OF AN ERROR,  NCALC .NE. NB,  AND NOT ALL I'S ARE
!  CALCULATED TO THE DESIRED ACCURACY.
!
!  NCALC .LT. 0:  AN ARGUMENT IS OUT OF RANGE. FOR EXAMPLE,
!     NB .LE. 0, IZE IS NOT 1 OR 2, OR IZE=1 AND ABS(X) .GE. EXPARG.
!     IN THIS CASE, THE B-VECTOR IS NOT CALCULATED, AND NCAL! IS
!     SET TO  MIN0(NB,0)-1  SO THAT NCALC .NE. NB.
!
!  NB .GT. NCALC .GT. 0: NOT ALL REQUESTED FUNCTION VALUES COULD
!     BE CALCULATED ACCURATELY.  THIS USUALLY OCCURS BECAUSE NB IS
!     MUCH LARGER THAN ABS(X).  IN THIS CASE, B(N) IS CALCULATED
!     TO THE DESIRED ACCURACY FOR  N .LE. NCALC,  BUT PRECISION
!     IS LOST FOR NCAL! .LT. N .LE. NB.  IF B(N) DOES NOT VANISH
!     FOR  N .GT. NCAL!  (BECAUSE IT IS TOO SMALL TO BE REPRESENTED),
!     AND  B(N)/B(NCALC) = 10**(-K), THEN ONLY THE FIRST NSIG-K
!     SIGNIFICANT FIGURES OF B(N) CAN BE TRUSTED.
!
!
! OTHER SUBPROGRAMS REQUIRED (SINGLE PRECISION VERSION)
!
!     EXP,GAMMA,AMAX1,SQRT,FLOAT,IFIX,MIN0
!
! OTHER SUBPROGRAMS REQUIRED (DOUBLE PRECISION VERSION)
!
!     DBLE,DEXP,DGAMMA,DMAX1,DSQRT,FLOAT,IFIX,MIN0,SNGL
!
!
! ACKNOWLEDGEMENT
!
!  THIS PROGRAM IS BASED ON A PROGRAM WRITTEN BY DAVID J.
!  SOOKNE (2) THAT COMPUTES VALUES OF THE BESSEL FUNCTIONS J OR
!  I OF REAL ARGUMENT AND INTEGER ORDER.  MODIFICATIONS INCLUDE
!  THE RESTRICTION OF THE COMPUTATION TO THE I BESSEL FUNCTION
!  OF NON-NEGATIVE REAL ARGUMENT, THE EXTENSION OF THE COMPUTATION
!  TO ARBITRARY POSITIVE ORDER, THE INCLUSION OF OPTIONAL
!  EXPONENTIAL SCALING, AND THE ELIMINATION OF MOST UNDERFLOW.
!
! REFERENCES
!
!  (1) OLVER, F. W. J., AND SOOKNE, D. J., "A NOTE ON BACKWARD
!        RECURRENCE ALGORITHMS," MATH. COMP. 26, 1972, PP 941-947.
!
!  (2) SOOKNE, D. J., "BESSEL FUNCTIONS OF REAL ARGUMENT AND
!        INTEGER ORDER," NBS JOUR. OF RES. B. 77B, 1973, PP 125-132.
!
!
!      MODIFIED BY: W. J. CODY
!                   APPLIED MATHEMATICS DIVISION
!                   ARGONNE NATIONAL LABORATORY
!                   ARGONNE, IL  60439
!
!      LATEST MODIFICATION: MAY 18, 1982
!
!-------------------------------------------------------------------
      INTEGER(I4B) IZE, L, MAGX, N, NB, NBMX, NCALC, NEND, NSIG, NSTART
!S    REAL              ALPHA,B,EM,EMPAL,EMP2AL,EN,ENMTEN,ENSIG,
!S   2 ENTEN,EXPARG,GAMMA,HALF,HALFX,ONE,P,PLAST,POLD,PSAVE,PSAVEL,
!S   3 RTNSIG,SUM,TEMPA,TEMPB,TEMPC,TEST,TOVER,TWO,X,XLARGE,ZERO
      REAL(DP) ALPHA, B, EM, EMPAL, EMP2AL, EN, ENMTEN, ENSIG,  &
      ENTEN, EXPARG, HALF, HALFX, ONE, P, PLAST, POLD, PSAVE, &
     PSAVEL, RTNSIG, SUM, TEMPA, TEMPB, TEMPC, TEST, TOVER, TWO, X,  &
      XLARGE, ZERO
      DIMENSION B(NB)
!-------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!-------------------------------------------------------------------
!S    DATA ONE,TWO,ZERO,HALF/1.0E0,2.0E0,0.0E0,0.5E0/
      DATA ONE, TWO, ZERO, HALF /1.0D0,2.0D0,0.0D0,0.5D0/
!-------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!-------------------------------------------------------------------
!S    DATA NSIG,XLARGE,EXPARG / 7,1.0E4,88.0E0/
!S    DATA ENTEN,ENSIG,RTNSIG/1.0E38,1.0E7,1.0E-2/
!S    DATA ENMTEN/1.2E-37/
      DATA NSIG, XLARGE, EXPARG /17,1.0D4,88.0D0/
      DATA ENTEN, ENSIG, RTNSIG /1.0D38,1.0D17,1.0D-4/
      DATA ENMTEN /1.2D-37/
!-------------------------------------------------------------------
!S    MAGX = IFIX(X)
!HSJ      MAGX = IFIX(SNGL(X))
      MAGX = INT(x)
      IF ((NB.GT.0) .AND. (X.GE.ZERO) .AND. (ALPHA.GE.ZERO) .AND. &
      (ALPHA.LT.ONE) .AND. (((IZE.EQ.1) .AND. (X.LE.EXPARG)) .OR. &
      ((IZE.EQ.2) .AND. (X.LE.XLARGE)))) GO TO 10
!-------------------------------------------------------------------
! ERROR RETURN -- X,NB,OR IZE IS OUT OF RANGE
!-------------------------------------------------------------------
      NCALC = MIN0(NB,0) - 1
      RETURN
!-------------------------------------------------------------------
! USE 2-TERM ASCENDING SERIES FOR SMALL X
!-------------------------------------------------------------------
   10 NCALC = NB
      IF (X.LT.RTNSIG) GO TO 210
!-------------------------------------------------------------------
! INITIALIZE THE FORWARD SWEEP, THE P-SEQUENCE OF OLVER
!-------------------------------------------------------------------
      NBMX = NB - MAGX
      N = MAGX + 1
!S    EN = FLOAT(N+N) + (ALPHA+ALPHA)
      EN = DBLE(FLOAT(N+N)) + (ALPHA+ALPHA)
      PLAST = ONE
      P = EN/X
!-------------------------------------------------------------------
! CALCULATE GENERAL SIGNIFICANCE TEST
!-------------------------------------------------------------------
      TEST = ENSIG + ENSIG
!S    IF (2*MAGX .GT. 5*NSIG) TEST = SQRT(TEST*P)
      IF (2*MAGX.GT.5*NSIG) TEST = DSQRT(TEST*P)
!S    IF (2*MAGX .LE. 5*NSIG) TEST = TEST / 1.585E0**MAGX
      IF (2*MAGX.LE.5*NSIG) TEST = TEST/1.585D0**MAGX
      IF (NBMX.LT.3) GO TO 30
!-------------------------------------------------------------------
! CALCULATE P-SEQUENCE UNTIL N = NB-1.  CHECK FOR POSSIBLE OVERFLOW.
!-------------------------------------------------------------------
      TOVER = ENTEN/ENSIG
      NSTART = MAGX + 2
      NEND = NB - 1
      DO 20 N=NSTART,NEND
        EN = EN + TWO
        POLD = PLAST
        PLAST = P
        P = EN*PLAST/X + POLD
        IF (P.GT.TOVER) GO TO 40
   20 CONTINUE
      N = NEND
!S    EN = FLOAT(N+N) + (ALPHA+ALPHA)
      EN = DBLE(FLOAT(N+N)) + (ALPHA+ALPHA)
!-------------------------------------------------------------------
! CALCULATE SPECIAL SIGNIFICANCE TEST FOR NBMX.GT.2.
!-------------------------------------------------------------------
!S    TEST = AMAX1(TEST,SQRT(PLAST*ENSIG)*SQRT(P+P))
      TEST = DMAX1(TEST,DSQRT(PLAST*ENSIG)*DSQRT(P+P))
!-------------------------------------------------------------------
! CALCULATE P-SEQUENCE UNTIL SIGNIFICANCE TEST PASSES
!-------------------------------------------------------------------
   30 N = N + 1
      EN = EN + TWO
      POLD = PLAST
      PLAST = P
      P = EN*PLAST/X + POLD
      IF (P.LT.TEST) GO TO 30
      GO TO 80
!-------------------------------------------------------------------
! TO AVOID OVERFLOW, DIVIDE P-SEQUENCE BY TOVER.  CALCULATE
! P-SEQUENCE UNTIL ABS(P).GT.1.
!-------------------------------------------------------------------
   40 TOVER = ENTEN
      P = P/TOVER
      PLAST = PLAST/TOVER
      PSAVE = P
      PSAVEL = PLAST
      NSTART = N + 1
   50 N = N + 1
      EN = EN + TWO
      POLD = PLAST
      PLAST = P
      P = EN*PLAST/X + POLD
      IF (P.LE.ONE) GO TO 50
      TEMPB = EN/X
!-------------------------------------------------------------------
! CALCULATE BACKWARD TEST, AND FIND NCALC, THE HIGHEST N
! SUCH THAT THE TEST IS PASSED.
!-------------------------------------------------------------------
      TEST = POLD*PLAST*(HALF-HALF/(TEMPB*TEMPB))/ENSIG
      P = PLAST*TOVER
      N = N - 1
      EN = EN - TWO
      NEND = MIN0(NB,N)
      DO 60 L=NSTART,NEND
        NCALC = L
        POLD = PSAVEL
        PSAVEL = PSAVE
        PSAVE = EN*PSAVEL/X + POLD
        IF (PSAVE*PSAVEL.GT.TEST) GO TO 70
   60 CONTINUE
      NCALC = NEND + 1
   70 NCALC = NCALC - 1
!-------------------------------------------------------------------
! INITIALIZE THE BACKWARD RECURSION AND THE NORMALIZATION SUM
!-------------------------------------------------------------------
   80 N = N + 1
      EN = EN + TWO
      TEMPB = ZERO
      TEMPA = ONE/P
!S    EM = FLOAT(N) - ONE
      EM = DBLE(FLOAT(N)) - ONE
      EMPAL = EM + ALPHA
      EMP2AL = (EM-ONE) + (ALPHA+ALPHA)
      SUM = TEMPA*EMPAL*EMP2AL/EM
      NEND = N - NB
      IF (NEND) 130, 110, 90
!-------------------------------------------------------------------
! RECUR BACKWARD VIA DIFFERENCE EQUATION, CALCULATING (BUT
! NOT STORING) B(N), UNTIL N = NB.
!-------------------------------------------------------------------
   90 DO 100 L=1,NEND
        N = N - 1
        EN = EN - TWO
        TEMPC = TEMPB
        TEMPB = TEMPA
        TEMPA = (EN*TEMPB)/X + TEMPC
        EM = EM - ONE
        EMP2AL = EMP2AL - ONE
        IF (N.EQ.1) GO TO 110
        IF (N.EQ.2) EMP2AL = ONE
        EMPAL = EMPAL - ONE
        SUM = (SUM+TEMPA*EMPAL)*EMP2AL/EM
  100 CONTINUE
!-------------------------------------------------------------------
! STORE B(NB)
!-------------------------------------------------------------------
  110 B(N) = TEMPA
      IF (NB.GT.1) GO TO 120
      SUM = (SUM+SUM) + TEMPA
      GO TO 190
!-------------------------------------------------------------------
! CALCULATE AND STORE B(NB-1)
!-------------------------------------------------------------------
  120 N = N - 1
      EN = EN - TWO
      B(N) = (EN*TEMPA)/X + TEMPB
      IF (N.EQ.1) GO TO 180
      EM = EM - ONE
      EMP2AL = EMP2AL - ONE
      IF (N.EQ.2) EMP2AL = ONE
      EMPAL = EMPAL - ONE
      SUM = (SUM+B(N)*EMPAL)*EMP2AL/EM
      GO TO 150
!-------------------------------------------------------------------
! N.LT.NB, SO STORE B(N) AND SET HIGHER ORDERS TO ZERO
!-------------------------------------------------------------------
  130 B(N) = TEMPA
      NEND = -NEND
      DO 140 L=1,NEND
        B(N+L) = ZERO
  140 CONTINUE
  150 NEND = N - 2
      IF (NEND.EQ.0) GO TO 170
!-------------------------------------------------------------------
! CALCULATE VIA DIFFERENCE EQUATION AND STORE B(N), UNTIL N = 2
!-------------------------------------------------------------------
      DO 160 L=1,NEND
        N = N - 1
        EN = EN - TWO
        B(N) = (EN*B(N+1))/X + B(N+2)
        EM = EM - ONE
        EMP2AL = EMP2AL - ONE
        IF (N.EQ.2) EMP2AL = ONE
        EMPAL = EMPAL - ONE
        SUM = (SUM+B(N)*EMPAL)*EMP2AL/EM
  160 CONTINUE
!-------------------------------------------------------------------
! CALCULATE B(1)
!-------------------------------------------------------------------
  170 B(1) = TWO*EMPAL*B(2)/X + B(3)
  180 SUM = (SUM+SUM) + B(1)
!-------------------------------------------------------------------
! NORMALIZE.  DIVIDE ALL B(N) BY SUM.
!-------------------------------------------------------------------
  190 CONTINUE
!S    IF (ALPHA.NE.ZERO)SUM=SUM*GAMMA(ONE+ALPHA)*(X*HALF)**(-ALPHA)
      IF (ALPHA.NE.ZERO) SUM = SUM*DGAMMA(ONE+ALPHA)*(X*HALF)**(-ALPHA)
!S    IF (IZE .EQ. 1) SUM = SUM * EXP(-X)
      IF (IZE.EQ.1) SUM = SUM*DEXP(-X)
      TEMPA = ENMTEN
      IF (SUM.GT.ONE) TEMPA = TEMPA*SUM
      DO 200 N=1,NB
        IF (B(N).LT.TEMPA) B(N) = ZERO
        B(N) = B(N)/SUM
  200 CONTINUE
      RETURN
!-------------------------------------------------------------------
! TWO-TERM ASCENDING SERIES FOR SMALL X
!-------------------------------------------------------------------
  210 TEMPA = ONE
      EMPAL = ONE + ALPHA
      HALFX = ZERO
      IF (X.GT.ENMTEN) HALFX = HALF*X
!S    IF (ALPHA .NE. ZERO) TEMPA = HALFX ** ALPHA / GAMMA(EMPAL)
      IF (ALPHA.NE.ZERO) TEMPA = HALFX**ALPHA/DGAMMA(EMPAL)
!S    IF (IZE .EQ. 2) TEMPA = TEMPA * EXP(-X)
      IF (IZE.EQ.2) TEMPA = TEMPA*DEXP(-X)
      TEMPB = ZERO
      IF ((X+ONE).GT.ONE) TEMPB = HALFX*HALFX
      B(1) = TEMPA + TEMPA*TEMPB/EMPAL
      IF ((X.NE.ZERO) .AND. (B(1).EQ.ZERO)) NCALC = 0
      IF (NB.EQ.1) GO TO 250
      IF (X.GT.ZERO) GO TO 230
      DO 220 N=2,NB
        B(N) = ZERO
  220 CONTINUE
      GO TO 250
!-------------------------------------------------------------------
! CALCULATE HIGHER ORDER FUNCTIONS
!-------------------------------------------------------------------
  230 TEMPC = HALFX
      TOVER = (ENMTEN+ENMTEN)/X
      IF (TEMPB.NE.ZERO) TOVER = ENMTEN/TEMPB
      DO 240 N=2,NB
        TEMPA = TEMPA/EMPAL
        EMPAL = EMPAL + ONE
        TEMPA = TEMPA*TEMPC
        IF (TEMPA.LE.TOVER*EMPAL) TEMPA = ZERO
        B(N) = TEMPA + TEMPA*TEMPB/EMPAL
        IF ((B(N).EQ.ZERO) .AND. (NCALC.GT.N)) NCALC = N - 1
  240 CONTINUE
  250 RETURN
      END SUBROUTINE RIBESL

      SUBROUTINE my_uertst (ierr, isub)
!******************************************************************************
!                                                                             *
!    ROUTINE NAME:  MY_UERTST                                                    *
!                                                                             *
!    PURPOSE:  prints error number and name of calling routine                *
!                                                                             *
!    AUTHOR:   Steve Reid - 9 September 1999                                  *
!              S & E Computing                                                *
!                                                                             *
!    INPUTS:   ierr - error number to be printed                              *
!              isub - name of calling subroutine to be printed                *
!                                                                             *
!    OUTPUTS:  none                                                           *
!                                                                             *
!    ROUTINES CALLED:  none                                                   *
!                                                                             *
!    LIMITATIONS:  none                                                       *
!                                                                             *
!******************************************************************************
      IMPLICIT NONE
      INTEGER(I4B) ierr
      CHARACTER*(*) isub

      IF (ierr .GT. 128) THEN
         WRITE (*, 6000)  'Terminal error', ierr, isub
      ELSEIF (ierr .GT. 64) THEN
         WRITE (*, 6000)  'Warning error with fix', ierr, isub
      ELSEIF (ierr .GT. 32) THEN
         WRITE (*, 6000)  'Warning error', ierr, isub
      ELSE
         WRITE (*, 6000)  'Undefined error', ierr, isub
      ENDIF
 6000 FORMAT (' *** ', a, i4, ' occurred in routine ', a)
      RETURN
      END SUBROUTINE my_uertst

      SUBROUTINE my_uerset (ierr, isub)
!-----------------------------------------------------------------------------
!-- set error level output routine, not used 
!-----------------------------------------------------------------------HSJ---
      IMPLICIT NONE
      INTEGER(I4B) ierr,isub

       RETURN
      END SUBROUTINE my_uerset

      REAL(DP) FUNCTION my_mmdele (iopt, arg, ierr)
!*******************************************************************************
!*                                                                             *
!*    ROUTINE NAME:  MMDELE                                                    *
!*                                                                             *
!*    PURPOSE:  computes the Complete Elliptic Integral of the Second Kind     *
!*                                                                             *
!*    AUTHOR:   Steve Reid - 19 December 1999                                  *
!*              S & E Computing                                                *
!*                                                                             *
!*    INPUTS:   iopt   - input option                                          *
!*                       1 = integral from 0 to pi/2 of SQRT (1 - arg          *
!*                           * (SIN (phi)) ** 2 d(phi) will be evaluated,      *
!*                           WHERE 0.0 <= arg <= 1.0                           *
!*                       2 = integral from 0 to pi/2 of SQRT (1 - arg ** 2     *
!*                           * (SIN (phi)) ** 2) d(phi) will be evaluated,     *
!*                           WHERE |arg| <= 1.0                                *
!*                       3 = integral from 0 to pi/2 of SQRT (1 - (1 - arg)    *
!*                           * (SIN (phi)) ** 2) d(phi) will be evaluated,     *
!*                           WHERE |arg| <= 1.0                                *
!*              arg    - input argument as described above                     *
!*                                                                             *
!*    OUTPUTS:  mmdele - value of the integral                                 *
!*              ierr   - error code (terminal error)                           *
!*                       129 = iopt < 1 or iopt > 3 and mmdele is set to       *
!*                             machine infinity                                *
!*                       130 = arg is out of range and mmdele is set to        *
!*                             machine infinity                                *
!*                                                                             *
!*    ROUTINES CALLED:  my_uertst                                                 *
!*                                                                             *
!*    LIMITATIONS:  NONE                                                       *
!*                                                                             *
!*******************************************************************************
      IMPLICIT NONE 
      INTEGER(I4B) i, ierr, iopt
      REAL(DP) a(9), b(9), suma, sumb
      REAL(DP) arg, arg1, dabsarg, epsilon, infinity, one, zero
      CHARACTER*6 isub /'mmdele'/
      DATA epsilon /0.44408920985006d-15/, one /1.0d0/, zero /0.0d0/
!     DATA infinity /+1.79769313486231d+308/  ! largest real*8 for HP
      DATA infinity /8.0d+307/                !approximate for pgf90  HSJ
      DATA a(1) /0.32519201550639d-3/
      DATA a(2) /0.43025377747931d-2/
      DATA a(3) /0.11785841008734d-1/
      DATA a(4) /0.11841925995501d-1/
      DATA a(5) /0.90355277375409d-2/
      DATA a(6) /0.11716766944658d-1/
      DATA a(7) /0.21836131405487d-1/
      DATA a(8) /0.56805223329308d-1/
      DATA a(9) /0.44314718058337d+0/
      DATA b(1) /0.72031696345716d-4/
      DATA b(2) /0.18645379184063d-2/
      DATA b(3) /0.10087958494375d-1/
      DATA b(4) /0.22660309891604d-1/
      DATA b(5) /0.32811069172721d-1/
      DATA b(6) /0.42672510126592d-1/
      DATA b(7) /0.58592707184265d-1/
      DATA b(8) /0.93749995116367d-1/
      DATA b(9) /0.24999999999746d+0/

      ierr    = 0
      dabsarg = dabs (arg)
      IF (iopt .LT. 1 .OR. iopt .GT. 3) THEN
         ierr = 129
      ELSEIF (iopt .EQ. 3) THEN
         IF (dabsarg .LE. one) THEN
            arg1 = arg
         ELSE
            ierr = 130
         ENDIF
      ELSEIF (iopt .EQ. 2) THEN
         IF (dabsarg .LT. one) THEN
            arg1 = one - arg * arg
         ELSE
            ierr = 130
         ENDIF
      ELSEIF (dabsarg .GE. zero .AND. dabsarg .LE. one) THEN
         arg1 = one - arg
      ELSE
         ierr = 130
      ENDIF
      IF (ierr .NE. 0) go to 9000
      IF (arg1 .LT. epsilon) THEN
         my_mmdele = one
      ELSE
         suma = zero
         sumb = zero
         DO i = 1, 9
            suma = (suma + a(i)) * arg1
            sumb = (sumb + b(i)) * arg1
         ENDDO
         my_mmdele = suma - dlog (arg1) * sumb
         my_mmdele = my_mmdele + one
      ENDIF
      go to 9990

!--Error EXIT.

 9000 CONTINUE
      my_mmdele = infinity
      CALL my_uertst (ierr, isub)

!--Normal EXIT.

 9990 CONTINUE
      RETURN
      END FUNCTION my_mmdele



      REAL(DP) FUNCTION DGAMMA (X)
!***BEGIN PROLOGUE  DGAMMA
!***PURPOSE  Compute the complete Gamma function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      DOUBLE PRECISION (GAMMA-S, DGAMMA-D, CGAMMA-C)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DGAMMA(X) calculates the double precision complete Gamma function
! for double precision argument X.
!
! Series for GAM        on the interval  0.          to  1.00000E+00
!                                        with weighted error   5.79E-32
!                                         log weighted error  31.24
!                               significant figures required  30.00
!                                    decimal places required  32.05
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, D9LGMC, DCSEVL, DGAMLM, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920618  Removed space from variable name.  (RWC, WRB)
!***END PROLOGUE  DGAMMA
      REAL(DP) X, GAMCS(42), DXREL, PI, SINPIY, SQ2PIL, XMAX, &
        XMIN, Y, DCD1MACH
      LOGICAL FIRST
!
      SAVE GAMCS, PI, SQ2PIL, NGAM, XMIN, XMAX, DXREL, FIRST
      DATA GAMCS(  1) / +.8571195590989331421920062399942E-2      /
      DATA GAMCS(  2) / +.4415381324841006757191315771652E-2      /
      DATA GAMCS(  3) / +.5685043681599363378632664588789E-1      /
      DATA GAMCS(  4) / -.4219835396418560501012500186624E-2      /
      DATA GAMCS(  5) / +.1326808181212460220584006796352E-2      /
      DATA GAMCS(  6) / -.1893024529798880432523947023886E-3      /
      DATA GAMCS(  7) / +.3606925327441245256578082217225E-4      /
      DATA GAMCS(  8) / -.6056761904460864218485548290365E-5      /
      DATA GAMCS(  9) / +.1055829546302283344731823509093E-5      /
      DATA GAMCS( 10) / -.1811967365542384048291855891166E-6      /
      DATA GAMCS( 11) / +.3117724964715322277790254593169E-7      /
      DATA GAMCS( 12) / -.5354219639019687140874081024347E-8      /
      DATA GAMCS( 13) / +.9193275519859588946887786825940E-9      /
      DATA GAMCS( 14) / -.1577941280288339761767423273953E-9      /
      DATA GAMCS( 15) / +.2707980622934954543266540433089E-10     /
      DATA GAMCS( 16) / -.4646818653825730144081661058933E-11     /
      DATA GAMCS( 17) / +.7973350192007419656460767175359E-12     /
      DATA GAMCS( 18) / -.1368078209830916025799499172309E-12     /
      DATA GAMCS( 19) / +.2347319486563800657233471771688E-13     /
      DATA GAMCS( 20) / -.4027432614949066932766570534699E-14     /
      DATA GAMCS( 21) / +.6910051747372100912138336975257E-15     /
      DATA GAMCS( 22) / -.1185584500221992907052387126192E-15     /
      DATA GAMCS( 23) / +.2034148542496373955201026051932E-16     /
      DATA GAMCS( 24) / -.3490054341717405849274012949108E-17     /
      DATA GAMCS( 25) / +.5987993856485305567135051066026E-18     /
      DATA GAMCS( 26) / -.1027378057872228074490069778431E-18     /
      DATA GAMCS( 27) / +.1762702816060529824942759660748E-19     /
      DATA GAMCS( 28) / -.3024320653735306260958772112042E-20     /
      DATA GAMCS( 29) / +.5188914660218397839717833550506E-21     /
      DATA GAMCS( 30) / -.8902770842456576692449251601066E-22     /
      DATA GAMCS( 31) / +.1527474068493342602274596891306E-22     /
      DATA GAMCS( 32) / -.2620731256187362900257328332799E-23     /
      DATA GAMCS( 33) / +.4496464047830538670331046570666E-24     /
      DATA GAMCS( 34) / -.7714712731336877911703901525333E-25     /
      DATA GAMCS( 35) / +.1323635453126044036486572714666E-25     /
      DATA GAMCS( 36) / -.2270999412942928816702313813333E-26     /
      DATA GAMCS( 37) / +.3896418998003991449320816639999E-27     /
      DATA GAMCS( 38) / -.6685198115125953327792127999999E-28     /
      DATA GAMCS( 39) / +.1146998663140024384347613866666E-28     /
      DATA GAMCS( 40) / -.1967938586345134677295103999999E-29     /
      DATA GAMCS( 41) / +.3376448816585338090334890666666E-30     /
      DATA GAMCS( 42) / -.5793070335782135784625493333333E-31     /
      DATA PI / 3.14159265358979323846264338327950E0 /
      DATA SQ2PIL / 0.91893853320467274178032973640562E0 /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DGAMMA
      IF (FIRST) THEN
         NGAM = INITDS (GAMCS, 42, 0.1*REAL(D1MACH(3)) )
!
         CALL DGAMLM (XMIN, XMAX)
         DXREL = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
!
      Y = ABS(X)
      IF (Y.GT.10.D0) GO TO 50
!
! COMPUTE GAMMA(X) FOR -XBND .LE. X .LE. XBND.  REDUCE INTERVAL AND FIND
! GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
!
      N = X
      IF (X.LT.0.D0) N = N - 1
      Y = X - N
      N = N - 1
      DGAMMA = 0.9375D0 + DCSEVL (2.D0*Y-1.D0, GAMCS, NGAM)
      IF (N.EQ.0) RETURN
!
      IF (N.GT.0) GO TO 30
!
! COMPUTE GAMMA(X) FOR X .LT. 1.0
!
      N = -N
      IF (X .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DGAMMA', 'X IS 0', 4, 2)
      IF (X .LT. 0.0 .AND. X+N-2 .EQ. 0.D0) CALL XERMSG ('SLATEC', &
         'DGAMMA', 'X IS A NEGATIVE INTEGER', 4, 2)
      IF (X .LT. (-0.5D0) .AND. ABS((X-AINT(X-0.5D0))/X) .LT. DXREL) &
         CALL XERMSG ('SLATEC', 'DGAMMA', &
         'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', &
         1, 1)
!
      DO 20 I=1,N
        DGAMMA = DGAMMA/(X+I-1 )
 20   CONTINUE
      RETURN
!
! GAMMA(X) FOR X .GE. 2.0 AND X .LE. 10.0
!
 30   DO 40 I=1,N
        DGAMMA = (Y+I) * DGAMMA
 40   CONTINUE
      RETURN
!
! GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
!
 50   IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'DGAMMA', &
         'X SO BIG GAMMA OVERFLOWS', 3, 2)
!
      DGAMMA = 0.D0
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'DGAMMA', &
         'X SO SMALL GAMMA UNDERFLOWS', 2, 1)
      IF (X.LT.XMIN) RETURN
!
      DGAMMA = EXP ((Y-0.5D0)*LOG(Y) - Y + SQ2PIL + D9LGMC(Y) )
      IF (X.GT.0.D0) RETURN
!
      IF (ABS((X-AINT(X-0.5D0))/X) .LT. DXREL) CALL XERMSG ('SLATEC', &
         'DGAMMA', &
         'ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER', 1, 1)
!
      SINPIY = SIN (PI*Y)
      IF (SINPIY .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DGAMMA', &
         'X IS A NEGATIVE INTEGER', 4, 2)
!
      DGAMMA = -PI/(Y*SINPIY*DGAMMA)
!
      RETURN
      END  FUNCTION DGAMMA


      SUBROUTINE DGAMLM (XMIN, XMAX)
!***BEGIN PROLOGUE  DGAMLM
!***PURPOSE  Compute the minimum and maximum bounds for the argument in
!            the Gamma function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A, R2
!***TYPE      DOUBLE PRECISION (GAMLIM-S, DGAMLM-D)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Calculate the minimum and maximum legal bounds for X in gamma(X).
! XMIN and XMAX are not the only bounds, but they are the only non-
! trivial ones to calculate.
!
!             Output Arguments --
! XMIN   double precision minimum legal value of X in gamma(X).  Any
!        smaller value of X might result in underflow.
! XMAX   double precision maximum legal value of X in gamma(X).  Any
!        larger value of X might cause overflow.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DGAMLM
      REAL(DP) XMIN, XMAX, ALNBIG, ALNSML, XLN, XOLD, D1MACH
!***FIRST EXECUTABLE STATEMENT  DGAMLM
      ALNSML = LOG(D1MACH(1))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = LOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.5D0)*XLN - XMIN - 0.2258D0 + ALNSML) &
          / (XMIN*XLN+0.5D0)
        IF (ABS(XMIN-XOLD).LT.0.005D0) GO TO 20
 10   CONTINUE
      CALL XERMSG ('SLATEC', 'DGAMLM', 'UNABLE TO FIND XMIN', 1, 2)
!
 20   XMIN = -XMIN + 0.01D0
!
      ALNBIG = LOG (D1MACH(2))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = LOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.5D0)*XLN - XMAX + 0.9189D0 - ALNBIG) &
          / (XMAX*XLN-0.5D0)
        IF (ABS(XMAX-XOLD).LT.0.005D0) GO TO 40
 30   CONTINUE
      CALL XERMSG ('SLATEC', 'DGAMLM', 'UNABLE TO FIND XMAX', 2, 2)
!
 40   XMAX = XMAX - 0.01D0
      XMIN = MAX (XMIN, -XMAX+1.D0)
!
      RETURN
      END       SUBROUTINE DGAMLM


      REAL(DP) FUNCTION DCSEVL (X, CS, N)
!***BEGIN PROLOGUE  DCSEVL
!***PURPOSE  Evaluate a Chebyshev series.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C3A2
!***TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D)
!***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!  Evaluate the N-term Chebyshev series CS at X.  Adapted from
!  a method presented in the paper by Broucke referenced below.
!
!       Input Arguments --
!  X    value at which the series is to be evaluated.
!  CS   array of N terms of a Chebyshev series.  In evaluating
!       CS, only half the first coefficient is summed.
!  N    number of terms in array CS.
!
!***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
!                 Chebyshev series, Algorithm 446, Communications of
!                 the A.C.M. 16, (1973) pp. 254-256.
!               L. Fox and I. B. Parker, Chebyshev Polynomials in
!                 Numerical Analysis, Oxford University Press, 1968,
!                 page 56.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900329  Prologued revised extensively and code rewritten to allow
!           X to be slightly outside interval (-1,+1).  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DCSEVL
      REAL(DP) B0, B1, B2, CS(*), ONEPL, TWOX, X, D1MACH
      LOGICAL FIRST
      SAVE FIRST, ONEPL
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DCSEVL
      IF (FIRST) ONEPL = 1.0D0 + D1MACH(4)
      FIRST = .FALSE.
      IF (N .LT. 1) CALL XERMSG ('SLATEC', 'DCSEVL', &
         'NUMBER OF TERMS .LE. 0', 2, 2)
      IF (N .GT. 1000) CALL XERMSG ('SLATEC', 'DCSEVL', &
         'NUMBER OF TERMS .GT. 1000', 3, 2)
      IF (ABS(X) .GT. ONEPL) CALL XERMSG ('SLATEC', 'DCSEVL', &
         'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1)
!
      B1 = 0.0D0
      B0 = 0.0D0
      TWOX = 2.0D0*X
      DO 10 I = 1,N
         B2 = B1
         B1 = B0
         NI = N + 1 - I
         B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
!
      DCSEVL = 0.5D0*(B0-B2)
!
      RETURN
      END      FUNCTION DCSEVL


 

      SUBROUTINE my_dgear (n, fcn, fcnj, x, h, y, xend, tol, meth, miter, &
                        index, iwk, wk, ierr)
!*******************************************************************************
!*                                                                             *
!*    ROUTINE NAME:  DGEAR                                                     *
!*                                                                             *
!*    PURPOSE:  solve a system of ordinary differential equations by the       *
!*              variable order Adams Predictor Corrector Method or Gears       *
!*              Method                                                         *
!*                                                                             *
!*    AUTHOR:   Steve Reid - 22 November 1999                                  *
!*              S & E Computing                                                *
!*                                                                             *
!*    INPUTS:   n     - number of first-order differential equations           *
!*              fcn   - name of subroutine for evaluating functions            *
!*              fcnj  - name of subroutine for evaluating the n by n Jacobian  *
!*                      Matrix of partial derivitives                          *
!*              x     - initial value of the independent variable              *
!*              h     - initial value of the step size                         *
!*              y     - n initial values of dependent variable (vector)        *
!*              xend  - value of x at which solution is desired next; the      *
!*                      integration will normally go beyond xend and the       *
!*                      routine will interpolate to x = xend.                  *
!*              tol   - tolerance or relative error control parameter (must be *
!*                      > 0).  tol is only used on the first call unless       *
!*                      index = -1.  tol should be atleast an order of magni-  *
!*                      tude larger than the unit roundoff error, but usually  *
!*                      not larger than 0.001.                                 *
!*              meth  - method indicator (used only on the first call, unless  *
!*                      index = -1                                             *
!*                      1 = Adams Method is to be used                         *
!*                      2 = the stiff Methods of Gear or the Backward          *
!*                          Differentiation Formulae are to be used            *
!*              miter - iteration method indicator                             *
!*                      0 = functional iteration is to be used; no partial     *
!*                          derivitives are needed; a dummy fcnj can be used   *
!*                      1 = Chord Method is to be used with an analytic        *
!*                          Jacobian; user-supplied fcnj                       *
!*                      2 = Chord Method is to be used with the Jacobian cal-  *
!*                          culated internally by finite differences; a dummy  *
!*                          fcnj can be used                                   *
!*                      3 = Chord Method is used with the Jacobian replaced by *
!*                          a diagonal approximation based on a directional    *
!*                          derivitive; a dummy fcnj can be used               *
!*                      -1 or -2 = use the same method as for miter = 1 or 2,  *
!*                          respectively, but use a banded Jacobian Matrix.    *
!*                          In either case, bandwidth information is passed    *
!*                          to dgear through a common block (i.e.              *
!*                             common /dband/ nlc, nuc                         *
!*                          where nlc and nuc are the number of lower and up-  *
!*                          per codiagonals, respectively                      *
!*              index - indicates type of call to dgear                        *
!*                      1  = first call for this problem                       *
!*                      0  = not the first call for this problem               *
!*                      -1 = not the first call for this problem and a new     *
!*                           value has been specified for tol                  *
!*                      2  = not the first call for this problem and integra-  *
!*                           tion is to continue until xend is hit exactly     *
!*                           (i.e. no interpolation is to be done).            *
!*                      3  = not the first call for this problem; integration  *
!*                           is to continue and control is to be returned to   *
!*                           the calling program after one step.  xend is      *
!*                           ignored.                                          *
!*              iwk   - work vector of length n (used only if miter = 1 or 2)  *
!*              wk    - work vector of length 4 * n + nmeth + nmiter, where    *
!*                      the value of nmeth depends on the value of meth:       *
!*                      If meth = 1, nmeth = n * 13                            *
!*                      If meth = 2, nmeth = n * 6                             *
!*                      The value of nmiter depends on the value of miter:     *
!*                      If miter = 1 or 2, nmiter = n * 13                     *
!*                      If miter = -1 or -2, nmiter = n * (2 * nlc + nuc + 3), *
!*                      where nlc = number of lower codiagonals                *
!*                            nuc = number of upper codiagonals                *
!*                      If miter = 3, nmiter = n                               *
!*                      If miter = 0, nmiter = 1                               *
!*                      wk must remain unchanged between successive calls      *
!*                      during integration.                                    *
!*                                                                             *
!*    OUTPUTS:  x     - value of the independent variable at which integration *
!*                      has been completed                                     *
!*              h     - step size used last, whether successful or not         *
!*              y     - vector of n dependent values computed at xend          *
!*              index - If the input value of index was 1 or -1 and the inte-  *
!*                      gration was successful, index is set to zero; other-   *
!*                      wise it is left unchanged.                             *
!*              ierr  - error code (warning error) [not implemented]           *
!*                        33 = (x + h) will equal x on the next step.  This    *
!*                             condition does not force the routine to halt,   *
!*                             but it indicates one of two conditions.  Either *
!*                             the user is specifying too much accuracy (i.e.  *
!*                             tol is too small) or the system of differential *
!*                             equations being solved is too stiff (either in  *
!*                             general or over the subinterval of the problem  *
!*                             being solved at the time of the error).  In     *
!*                             this case, the user should make miter non-zero. *
!*                    - error code (warning error with fix) [not implemented]  *
!*                        66 = error test failed; h was decremented by 0.1     *
!*                             until success was achieved.                     *
!*                        67 = corrector convergence could not be achieved; h  *
!*                             was decremented by 0.1 until success was        *
!*                             achieved.                                       *
!*                    - error code (terminal error) [only 135 implemented]     *
!*                       132 = integration was halted after failing to pass    *
!*                             the error test, even after reducing h by a      *
!*                             factor of 1.0e10 from its original value.       *
!*                       133 = integration was halted after failing to achieve *
!*                             corrector convergence, even after reducing h by *
!*                             a factor of 1.0e10 from its initial value.      *
!*                       134 = after some initial success, the integration was *
!*                             halted either by repeated error test failures   *
!*                             or by a test on tol.                            *
!*                       135 = one of the input arguments n, x, h, xend, tol,  *
!*                             meth, miter, or index was specified incorrectly *
!*                       136 = index had an input value of -1, but the desired *
!*                             changes of parameters were not implemented be-  *
!*                             cause xend was not beyond x.  Interpolation to  *
!*                             x = xend was performed.  Retry with index = -1  *
!*                             and a new value for xend.                       *
!*                                                                             *
!*    ROUTINES CALLED:  lsodes, uertst                                         *
!*                                                                             *
!*    LIMITATIONS:  none                                                       *
!*                                                                             *
!*******************************************************************************
      IMPLICIT NONE
      INTEGER(I4B) ierr, index, iopt, istate, itask, itol, iwk(*), lenrat, &
              liw, lrw, meth, mf, miter, n, nnz,neq(1)
      REAL(DP) atolv(1),rtolv(1)
      REAL(DP) atol, infinity, h, rtol, tol, tout, wk(*), x, xend, y(n)
      CHARACTER*6 isub /'dgear'/
!      data infinity /+1.79769313486231e+308/  ! largest real*8 for HP
      DATA infinity /8.0d+307/                !approximate for pgf90  HS
      DATA lenrat /2/     ! number of words used for double precision
      EXTERNAL fcn, fcnj

       ! in case user calls this subrotuine:
          STOP 'replace_imsl,sub dgear not useable'
       ! problem is with calls to xerrwv which is hollerith and not
       ! currently ported
!
!--Validate the input arguments and set up arguments for call to lsodes.
!
      neq(1) =n 
      IF (miter .LT. -2 .OR. miter .GT. 3 .OR. meth .LT. 1 .OR. &
          meth .GT. 2 .OR. tol .LE. 0.0 .OR. n .LT. 0 .OR. &
          (x - xend) * h .GE. 0.0 .OR. index .LT. -1 .OR. &
          index .GT. 3) THEN
         go to 9000
      ENDIF
      tout   = x + h
      itol   = 1            ! indicates atol is a scalar
      rtol   = tol
      atol   = 0.0d0        ! use relative error control
      itask  = 1            ! normal computation of y at t = tout
      iopt   = 0            ! no optional inputs
      liw    = 30           ! minimum size according to isodes prolog
!
!--Set the method flag.
!
      IF (meth .EQ. 1) THEN
         mf  = 10           ! nonstiff method of Adams; no Jacobian used
         lrw = 20 + 16 * n  ! according to isodes prolog
      ELSEIF (meth .EQ. 2) THEN
         mf   = 121   ! stiff method of Gear; user-supplied sparse Jacob
         lrw = 20 + (2 + 1. / lenrat) * nnz + (11 + 9. / lenrat) * n  !
         STOP 'DGEAR is not programmed to handle stiff method of Gear'
      ELSE
         go to 9000
      ENDIF
!
!--Set the integration state indicator.
!
      IF (index .EQ. 1) THEN
         istate = 1  ! first call to dgear for this problem
      ELSE
         istate = 2  ! continuation
      ENDIF
      rtolv(1) = rtol ; atolv(1) = atol ! required for ifort
      CALL lsodes (fcn, neq, y, x, tout, itol, rtolv, atolv, itask, istate, &
                   iopt, wk, lrw, iwk, liw, fcnj, mf)
      IF (istate .EQ. 2) index = 0   ! integration was successful
      go to 9990
!
!--Error exits.
!
 9000 CONTINUE
      ierr = 135
      CALL my_uertst (ierr, isub)
!
!--Normal exit.
!
 9990 CONTINUE
      RETURN
      END       SUBROUTINE my_dgear




      SUBROUTINE DPPFA (AP, N, INFO)
!***BEGIN PROLOGUE  DPPFA
!***PURPOSE  Factor a real symmetric positive definite matrix stored in
!            packed form.
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1B
!***TYPE      DOUBLE PRECISION (SPPFA-S, DPPFA-D, CPPFA-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED,
!             POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DPPFA factors a double precision symmetric positive definite
!     matrix stored in packed form.
!
!     DPPFA is usually called by DPPCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (time for DPPCO) = (1 + 18/N)*(time for DPPFA) .
!
!     On Entry
!
!        AP      DOUBLE PRECISION (N*(N+1)/2)
!                the packed form of a symmetric matrix  A .  The
!                columns of the upper triangle are stored sequentially
!                in a one-dimensional array of length  N*(N+1)/2 .
!                See comments below for details.
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        AP      an upper triangular matrix  R , stored in packed
!                form, so that  A = TRANS(R)*R .
!
!        INFO    INTEGER
!                = 0  for normal return.
!                = K  if the leading minor of order  K  is not
!                     positive definite.
!
!
!     Packed Storage
!
!          The following program segment will pack the upper
!          triangle of a symmetric matrix.
!
!                K = 0
!                DO 20 J = 1, N
!                   DO 10 I = 1, J
!                      K = K + 1
!                      AP(K) = A(I,J)
!             10    CONTINUE
!             20 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DPPFA
      INTEGER(I4B) N,INFO
      REAL(DP) AP(*)
!
      REAL(DP) DDOT,T
      REAL(DP) S
      INTEGER(I4B) J,JJ,JM1,K,KJ,KK
!***FIRST EXECUTABLE STATEMENT  DPPFA
         JJ = 0
         DO 30 J = 1, N
            INFO = J
            S = 0.0D0
            JM1 = J - 1
            KJ = JJ
            KK = 0
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               KJ = KJ + 1
               T = AP(KJ) - DDOT(K-1,AP(KK+1),1,AP(JJ+1),1)
               KK = KK + K
               T = T/AP(KK)
               AP(KJ) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            JJ = JJ + J
            S = AP(JJ) - S
            IF (S .LE. 0.0D0) GO TO 40
            AP(JJ) = SQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END       SUBROUTINE DPPFA 

      SUBROUTINE DPPSL (AP, N, B)
!***BEGIN PROLOGUE  DPPSL
!***PURPOSE  Solve the real symmetric positive definite system using
!            the factors computed by DPPCO or DPPFA.
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1B
!***TYPE      DOUBLE PRECISION (SPPSL-S, DPPSL-D, CPPSL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, PACKED,
!             POSITIVE DEFINITE, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DPPSL solves the double precision symmetric positive definite
!     system A * X = B
!     using the factors computed by DPPCO or DPPFA.
!
!     On Entry
!
!        AP      DOUBLE PRECISION (N*(N+1)/2)
!                the output from DPPCO or DPPFA.
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains
!        a zero on the diagonal.  Technically this indicates
!        singularity, but it is usually caused by improper subroutine
!        arguments.  It will not occur if the subroutines are called
!        correctly and  INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DPPCO(AP,N,RCOND,Z,INFO)
!           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ...
!           DO 10 J = 1, P
!              CALL DPPSL(AP,N,C(1,J))
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DPPSL
      INTEGER(I4B) N
      REAL(DP) AP(*),B(*)
!
      REAL(DP) DDOT,T
      INTEGER(I4B) K,KB,KK
!***FIRST EXECUTABLE STATEMENT  DPPSL
      KK = 0
      DO 10 K = 1, N
         T = DDOT(K-1,AP(KK+1),1,B(1),1)
         KK = KK + K
         B(K) = (B(K) - T)/AP(KK)
   10 CONTINUE
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/AP(KK)
         KK = KK - K
         T = -B(K)
         CALL DAXPY(K-1,T,AP(KK+1),1,B(1),1)
   20 CONTINUE
      RETURN
      END       SUBROUTINE DPPSL



      SUBROUTINE DSCAL (N, DA, DX, INCX)
!***BEGIN PROLOGUE  DSCAL
!***PURPOSE  Multiply a vector by a constant.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A6
!***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DA  double precision scale factor
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!       DX  double precision result (unchanged if N.LE.0)
!
!     Replace double precision DX by double precision DA*DX.
!     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DSCAL
      REAL(DP) DA, DX(*)
      INTEGER(I4B) I, INCX, IX, M, MP1, N
!***FIRST EXECUTABLE STATEMENT  DSCAL
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GOTO 20
!
!     Code for increment not equal to 1.
!
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
!
!     Code for increment equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 5.
!
   20 M = MOD(N,5)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END       SUBROUTINE DSCAL





      REAL(DP) FUNCTION D9LGMC (X)
!***BEGIN PROLOGUE  D9LGMC
!***SUBSIDIARY
!***PURPOSE  Compute the log Gamma correction factor so that
!            LOG(DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-5.)*LOG(X) - X
!            + D9LGMC(X).
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
!             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the log gamma correction factor for X .GE. 10. so that
! LOG (DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + D9lGMC(X)
!
! Series for ALGM       on the interval  0.          to  1.00000E-02
!                                        with weighted error   1.28E-31
!                                         log weighted error  30.89
!                               significant figures required  29.81
!                                    decimal places required  31.48
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  D9LGMC
      REAL(DP) X, ALGMCS(15), XBIG, XMAX, D1MACH
      LOGICAL FIRST
      SAVE ALGMCS, NALGM, XBIG, XMAX, FIRST
      DATA ALGMCS(  1) / +.1666389480451863247205729650822E+0      /
      DATA ALGMCS(  2) / -.1384948176067563840732986059135E-4      /
      DATA ALGMCS(  3) / +.9810825646924729426157171547487E-8      /
      DATA ALGMCS(  4) / -.1809129475572494194263306266719E-10     /
      DATA ALGMCS(  5) / +.6221098041892605227126015543416E-13     /
      DATA ALGMCS(  6) / -.3399615005417721944303330599666E-15     /
      DATA ALGMCS(  7) / +.2683181998482698748957538846666E-17     /
      DATA ALGMCS(  8) / -.2868042435334643284144622399999E-19     /
      DATA ALGMCS(  9) / +.3962837061046434803679306666666E-21     /
      DATA ALGMCS( 10) / -.6831888753985766870111999999999E-23     /
      DATA ALGMCS( 11) / +.1429227355942498147573333333333E-24     /
      DATA ALGMCS( 12) / -.3547598158101070547199999999999E-26     /
      DATA ALGMCS( 13) / +.1025680058010470912000000000000E-27     /
      DATA ALGMCS( 14) / -.3401102254316748799999999999999E-29     /
      DATA ALGMCS( 15) / +.1276642195630062933333333333333E-30     /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  D9LGMC
      IF (FIRST) THEN
         NALGM = INITDS (ALGMCS, 15, REAL(D1MACH(3)) )
         XBIG = 1.0D0/SQRT(D1MACH(3))
         XMAX = EXP (MIN(LOG(D1MACH(2)/12.D0), -LOG(12.D0*D1MACH(1))))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LT. 10.D0) CALL XERMSG ('SLATEC', 'D9LGMC', &
         'X MUST BE GE 10', 1, 2)
      IF (X.GE.XMAX) GO TO 20
!
      D9LGMC = 1.D0/(12.D0*X)
      IF (X.LT.XBIG) D9LGMC = DCSEVL (2.0D0*(10.D0/X)**2-1.D0, ALGMCS, &
        NALGM) / X
      RETURN
!
 20   D9LGMC = 0.D0
      CALL XERMSG ('SLATEC', 'D9LGMC', 'X SO BIG D9LGMC UNDERFLOWS', 2, &
         1)
      RETURN
!
      END FUNCTION D9LGMC 

      INTEGER(I4B) FUNCTION INITDS (OS, NOS, ETA)
!***BEGIN PROLOGUE  INITDS
!***PURPOSE  Determine the number of terms needed in an orthogonal
!            polynomial series so that it meets a specified accuracy.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C3A2
!***TYPE      DOUBLE PRECISION (INITS-S, INITDS-D)
!***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
!             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!  Initialize the orthogonal series, represented by the array OS, so
!  that INITDS is the number of terms needed to insure the error is no
!  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
!  machine precision.
!
!             Input Arguments --
!   OS     double precision array of NOS coefficients in an orthogonal
!          series.
!   NOS    number of coefficients in OS.
!   ETA    single precision scalar containing requested accuracy of
!          series.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891115  Modified error message.  (WRB)
!   891115  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  INITDS
       REAL(DP)  OS(*)
!      real*8, dimension(:) :: OS
!***FIRST EXECUTABLE STATEMENT  INITDS
      IF (NOS .LT. 1) CALL XERMSG ('SLATEC', 'INITDS', &
         'Number of coefficients is less than 1', 2, 1)
!
      ERR = 0.
      DO 10 II = 1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(REAL(OS(I)))
        IF (ERR.GT.ETA) GO TO 20
   10 CONTINUE
!
   20 IF (I .EQ. NOS) CALL XERMSG ('SLATEC', 'INITDS', &
         'Chebyshev series too short for specified accuracy', 1, 1)
      INITDS = I
!
      RETURN
      END FUNCTION INITDS





      SUBROUTINE lsodes (f, neq, y, t, tout, itol, rtol, atol, itask, &
                  istate, iopt, rwork, lrw, iwork, liw, jac, mf)
      EXTERNAL f, jac
      INTEGER neq, itol, itask, istate, iopt, lrw, iwork, liw, mf,iwmint(1)
      REAL(DP) y, t, tout, rtol, atol, rwork
      DIMENSION neq(1), y(1), rtol(1), atol(1), rwork(lrw), iwork(liw)
!-----------------------------------------------------------------------
! this is the march 30, 1987 version of
! lsodes.. livermore solver for ordinary differential equations
!          with general sparse jacobian matrices.
! this version is in double precision.
!
! lsodes solves the initial value problem for stiff or nonstiff
! systems of first order ode-s,
!     dy/dt = f(t,y) ,  or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
! lsodes is a variant of the lsode package, and is intended for
! problems in which the jacobian matrix df/dy has an arbitrary
! sparse structure (when the problem is stiff).
!
! authors..      alan c. hindmarsh,
!                computing and mathematics research division, l-316
!                lawrence livermore national laboratory
!                livermore, ca 94550.
!
! and            andrew h. sherman
!                j. s. nolen and associates
!                houston, tx 77084
!-----------------------------------------------------------------------
! references..
! 1.  alan c. hindmarsh,  odepack, a systematized collection of ode
!     solvers, in scientific computing, r. s. stepleman et al. (eds.),
!     north-holland, amsterdam, 1983, pp. 55-64.
!
! 2.  s. c. eisenstat, m. c. gursky, m. h. schultz, and a. h. sherman,
!     yale sparse matrix package.. i. the symmetric codes,
!     int. j. num. meth. eng., 18 (1982), pp. 1145-1151.
!
! 3.  s. c. eisenstat, m. c. gursky, m. h. schultz, and a. h. sherman,
!     yale sparse matrix package.. ii. the nonsymmetric codes,
!     research report no. 114, dept. of computer sciences, yale
!     university, 1977.
!-----------------------------------------------------------------------
! summary of usage.
!
! communication between the user and the lsodes package, for normal
! situations, is summarized here.  this summary describes only a subset
! of the full set of options available.  see the full description for
! details, including optional communication, nonstandard options,
! and instructions for special situations.  see also the example
! problem (with program and output) following this summary.
!
! a. first provide a subroutine of the form..
!               subroutine f (neq, t, y, ydot)
!               dimension y(neq), ydot(neq)
! which supplies the vector function f by loading ydot(i) with f(i).
!
! b. next determine (or guess) whether or not the problem is stiff.
! stiffness occurs when the jacobian matrix df/dy has an eigenvalue
! whose real part is negative and large in magnitude, compared to the
! reciprocal of the t span of interest.  if the problem is nonstiff,
! use a method flag mf = 10.  if it is stiff, there are two standard
! for the method flag, mf = 121 and mf = 222.  in both cases, lsodes
! requires the jacobian matrix in some form, and it treats this matrix
! in general sparse form, with sparsity structure determined internally.
! (for options where the user supplies the sparsity structure, see
! the full description of mf below.)
!
! c. if the problem is stiff, you are encouraged to supply the jacobian
! directly (mf = 121), but if this is not feasible, lsodes will
! compute it internally by difference quotients (mf = 222).
! if you are supplying the jacobian, provide a subroutine of the form..
!               subroutine jac (neq, t, y, j, ian, jan, pdj)
!               dimension y(1), ian(1), jan(1), pdj(1)
! here neq, t, y, and j are input arguments, and the jac routine is to
! load the array pdj (of length neq) with the j-th column of df/dy.
! i.e., load pdj(i) with df(i)/dy(j) for all relevant values of i.
! the arguments ian and jan should be ignored for normal situations.
! lsodes will call the jac routine with j = 1,2,...,neq.
! only nonzero elements need be loaded.  usually, a crude approximation
! to df/dy, possibly with fewer nonzero elements, will suffice.
!
! d. write a main program which calls subroutine lsodes once for
! each point at which answers are desired.  this should also provide
! for possible use of logical unit 6 for output of error messages
! by lsodes.  on the first call to lsodes, supply arguments as follows..
! f      = name of subroutine for right-hand side vector f.
!          this name must be declared external in calling program.
! neq    = number of first order ode-s.
! y      = array of initial values, of length neq.
! t      = the initial value of the independent variable.
! tout   = first point where output is desired (.ne. t).
! itol   = 1 or 2 according as atol (below) is a scalar or array.
! rtol   = relative tolerance parameter (scalar).
! atol   = absolute tolerance parameter (scalar or array).
!          the estimated local error in y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
!             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.
!          thus the local error test passes if, in each component,
!          either the absolute error is less than atol (or atol(i)),
!          or the relative error is less than rtol.
!          use rtol = 0.0 for pure absolute error control, and
!          use atol = 0.0 (or atol(i) = 0.0) for pure relative error
!          control.  caution.. actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! itask  = 1 for normal computation of output values of y at t = tout.
! istate = integer flag (input and output).  set istate = 1.
! iopt   = 0 to indicate no optional inputs used.
! rwork  = real work array of length at least..
!             20 + 16*neq            for mf = 10,
!             20 + (2 + 1./lenrat)*nnz + (11 + 9./lenrat)*neq
!                                    for mf = 121 or 222,
!          where..
!          nnz    = the number of nonzero elements in the sparse
!                   jacobian (if this is unknown, use an estimate), and
!          lenrat = the real to integer wordlength ratio (usually 1 in
!                   single precision and 2 in double precision).
!          in any case, the required size of rwork cannot generally
!          be predicted in advance if mf = 121 or 222, and the value
!          above is a rough estimate of a crude lower bound.  some
!          experimentation with this size may be necessary.
!          (when known, the correct required length is an optional
!          output, available in iwork(17).)
! lrw    = declared length of rwork (in user-s dimension).
! iwork  = integer work array of length at least 30.
! liw    = declared length of iwork (in user-s dimension).
! jac    = name of subroutine for jacobian matrix (mf = 121).
!          if used, this name must be declared external in calling
!          program.  if not used, pass a dummy name.
! mf     = method flag.  standard values are..
!          10  for nonstiff (adams) method, no jacobian used.
!          121 for stiff (bdf) method, user-supplied sparse jacobian.
!          222 for stiff method, internally generated sparse jacobian.
! note that the main program must declare arrays y, rwork, iwork,
! and possibly atol.
!
! e. the output from the first call (or any call) is..
!      y = array of computed values of y(t) vector.
!      t = corresponding value of independent variable (normally tout).
! istate = 2  if lsodes was successful, negative otherwise.
!          -1 means excess work done on this call (perhaps wrong mf).
!          -2 means excess accuracy requested (tolerances too small).
!          -3 means illegal input detected (see printed message).
!          -4 means repeated error test failures (check all inputs).
!          -5 means repeated convergence failures (perhaps bad jacobian
!             supplied or wrong choice of mf or tolerances).
!          -6 means error weight became zero during problem. (solution
!             component i vanished, and atol or atol(i) = 0.)
!          -7 means a fatal error return flag came from the sparse
!             solver cdrv by way of prjs or slss.  should never happen.
!          a return with istate = -1, -4, or -5 may result from using
!          an inappropriate sparsity structure, one that is quite
!          different from the initial structure.  consider calling
!          lsodes again with istate = 3 to force the structure to be
!          reevaluated.  see the full description of istate below.
!
! f. to continue the integration after a successful return, simply
! reset tout and call lsodes again.  no other parameters need be reset.
!
!-----------------------------------------------------------------------
! example problem.
!
! the following is a simple example problem, with the coding
! needed for its solution by lsodes.  the problem is from chemical
! kinetics, and consists of the following 12 rate equations..
!    dy1/dt  = -rk1*y1
!    dy2/dt  = rk1*y1 + rk11*rk14*y4 + rk19*rk14*y5
!                - rk3*y2*y3 - rk15*y2*y12 - rk2*y2
!    dy3/dt  = rk2*y2 - rk5*y3 - rk3*y2*y3 - rk7*y10*y3
!                + rk11*rk14*y4 + rk12*rk14*y6
!    dy4/dt  = rk3*y2*y3 - rk11*rk14*y4 - rk4*y4
!    dy5/dt  = rk15*y2*y12 - rk19*rk14*y5 - rk16*y5
!    dy6/dt  = rk7*y10*y3 - rk12*rk14*y6 - rk8*y6
!    dy7/dt  = rk17*y10*y12 - rk20*rk14*y7 - rk18*y7
!    dy8/dt  = rk9*y10 - rk13*rk14*y8 - rk10*y8
!    dy9/dt  = rk4*y4 + rk16*y5 + rk8*y6 + rk18*y7
!    dy10/dt = rk5*y3 + rk12*rk14*y6 + rk20*rk14*y7
!                + rk13*rk14*y8 - rk7*y10*y3 - rk17*y10*y12
!                - rk6*y10 - rk9*y10
!    dy11/dt = rk10*y8
!    dy12/dt = rk6*y10 + rk19*rk14*y5 + rk20*rk14*y7
!                - rk15*y2*y12 - rk17*y10*y12
!
! with rk1 = rk5 = 0.1,  rk4 = rk8 = rk16 = rk18 = 2.5,
!      rk10 = 5.0,  rk2 = rk6 = 10.0,  rk14 = 30.0,
!      rk3 = rk7 = rk9 = rk11 = rk12 = rk13 = rk19 = rk20 = 50.0,
!      rk15 = rk17 = 100.0.
!
! the t interval is from 0 to 1000, and the initial conditions
! are y1 = 1, y2 = y3 = ... = y12 = 0.  the problem is stiff.
!
! the following coding solves this problem with lsodes, using mf = 121
! and printing results at t = .1, 1., 10., 100., 1000.  it uses
! itol = 1 and mixed relative/absolute tolerance controls.
! during the run and at the end, statistical quantities of interest
! are printed (see optional outputs in the full description below).
!
!     external fex, jex
!     real*8 atol, rtol, rwork, t, tout, y
!     dimension y(12), rwork(500), iwork(30)
!     data lrw/500/, liw/30/
!     neq = 12
!     do 10 i = 1,neq
! 10    y(i) = 0.0d0
!     y(1) = 1.0d0
!     t = 0.0d0
!     tout = 0.1d0
!     itol = 1
!     rtol = 1.0d-4
!     atol = 1.0d-6
!     itask = 1
!     istate = 1
!     iopt = 0
!     mf = 121
!     do 40 iout = 1,5
!       call lsodes (fex, neq, y, t, tout, itol, rtol, atol,
!    1     itask, istate, iopt, rwork, lrw, iwork, liw, jex, mf)
!       write(6,30)t,iwork(11),rwork(11),(y(i),i=1,neq)
! 30    format(//7h at t =,e11.3,4x,
!    1    12h no. steps =,i5,4x,12h last step =,e11.3/
!    2    13h  y array =  ,4e14.5/13x,4e14.5/13x,4e14.5)
!       if (istate .lt. 0) go to 80
!       tout = tout*10.0d0
! 40    continue
!     lenrw = iwork(17)
!     leniw = iwork(18)
!     nst = iwork(11)
!     nfe = iwork(12)
!     nje = iwork(13)
!     nlu = iwork(21)
!     nnz = iwork(19)
!     nnzlu = iwork(25) + iwork(26) + neq
!     write (6,70) lenrw,leniw,nst,nfe,nje,nlu,nnz,nnzlu
! 70  format(//22h required rwork size =,i4,15h   iwork size =,i4/
!    1   12h no. steps =,i4,12h   no. f-s =,i4,12h   no. j-s =,i4,
!    2   13h   no. lu-s =,i4/23h no. of nonzeros in j =,i5,
!    3   26h   no. of nonzeros in lu =,i5)
!     stop
! 80  write(6,90)istate
! 90  format(///22h error halt.. istate =,i3)
!     stop
!     end
!
!     subroutine fex (neq, t, y, ydot)
!     real*8 t, y, ydot
!     real*8 rk1, rk2, rk3, rk4, rk5, rk6, rk7, rk8, rk9,
!    1   rk10, rk11, rk12, rk13, rk14, rk15, rk16, rk17
!     dimension y(12), ydot(12)
!     data rk1/0.1d0/, rk2/10.0d0/, rk3/50.0d0/, rk4/2.5d0/, rk5/0.1d0/,
!    1   rk6/10.0d0/, rk7/50.0d0/, rk8/2.5d0/, rk9/50.0d0/, rk10/5.0d0/,
!    2   rk11/50.0d0/, rk12/50.0d0/, rk13/50.0d0/, rk14/30.0d0/,
!    3   rk15/100.0d0/, rk16/2.5d0/, rk17/100.0d0/, rk18/2.5d0/,
!    4   rk19/50.0d0/, rk20/50.0d0/
!     ydot(1)  = -rk1*y(1)
!     ydot(2)  = rk1*y(1) + rk11*rk14*y(4) + rk19*rk14*y(5)
!    1           - rk3*y(2)*y(3) - rk15*y(2)*y(12) - rk2*y(2)
!     ydot(3)  = rk2*y(2) - rk5*y(3) - rk3*y(2)*y(3) - rk7*y(10)*y(3)
!    1           + rk11*rk14*y(4) + rk12*rk14*y(6)
!     ydot(4)  = rk3*y(2)*y(3) - rk11*rk14*y(4) - rk4*y(4)
!     ydot(5)  = rk15*y(2)*y(12) - rk19*rk14*y(5) - rk16*y(5)
!     ydot(6)  = rk7*y(10)*y(3) - rk12*rk14*y(6) - rk8*y(6)
!     ydot(7)  = rk17*y(10)*y(12) - rk20*rk14*y(7) - rk18*y(7)
!     ydot(8)  = rk9*y(10) - rk13*rk14*y(8) - rk10*y(8)
!     ydot(9)  = rk4*y(4) + rk16*y(5) + rk8*y(6) + rk18*y(7)
!     ydot(10) = rk5*y(3) + rk12*rk14*y(6) + rk20*rk14*y(7)
!    1           + rk13*rk14*y(8) - rk7*y(10)*y(3) - rk17*y(10)*y(12)
!    2           - rk6*y(10) - rk9*y(10)
!     ydot(11) = rk10*y(8)
!     ydot(12) = rk6*y(10) + rk19*rk14*y(5) + rk20*rk14*y(7)
!    1           - rk15*y(2)*y(12) - rk17*y(10)*y(12)
!     return
!     end
!
!     subroutine jex (neq, t, y, j, ia, ja, pdj)
!     real*8 t, y, pdj
!     real*8 rk1, rk2, rk3, rk4, rk5, rk6, rk7, rk8, rk9,
!    1   rk10, rk11, rk12, rk13, rk14, rk15, rk16, rk17
!     dimension y(1), ia(1), ja(1), pdj(1)
!     data rk1/0.1d0/, rk2/10.0d0/, rk3/50.0d0/, rk4/2.5d0/, rk5/0.1d0/,
!    1   rk6/10.0d0/, rk7/50.0d0/, rk8/2.5d0/, rk9/50.0d0/, rk10/5.0d0/,
!    2   rk11/50.0d0/, rk12/50.0d0/, rk13/50.0d0/, rk14/30.0d0/,
!    3   rk15/100.0d0/, rk16/2.5d0/, rk17/100.0d0/, rk18/2.5d0/,
!    4   rk19/50.0d0/, rk20/50.0d0/
!     go to (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), j
! 1   pdj(1) = -rk1
!     pdj(2) = rk1
!     return
! 2   pdj(2) = -rk3*y(3) - rk15*y(12) - rk2
!     pdj(3) = rk2 - rk3*y(3)
!     pdj(4) = rk3*y(3)
!     pdj(5) = rk15*y(12)
!     pdj(12) = -rk15*y(12)
!     return
! 3   pdj(2) = -rk3*y(2)
!     pdj(3) = -rk5 - rk3*y(2) - rk7*y(10)
!     pdj(4) = rk3*y(2)
!     pdj(6) = rk7*y(10)
!     pdj(10) = rk5 - rk7*y(10)
!     return
! 4   pdj(2) = rk11*rk14
!     pdj(3) = rk11*rk14
!     pdj(4) = -rk11*rk14 - rk4
!     pdj(9) = rk4
!     return
! 5   pdj(2) = rk19*rk14
!     pdj(5) = -rk19*rk14 - rk16
!     pdj(9) = rk16
!     pdj(12) = rk19*rk14
!     return
! 6   pdj(3) = rk12*rk14
!     pdj(6) = -rk12*rk14 - rk8
!     pdj(9) = rk8
!     pdj(10) = rk12*rk14
!     return
! 7   pdj(7) = -rk20*rk14 - rk18
!     pdj(9) = rk18
!     pdj(10) = rk20*rk14
!     pdj(12) = rk20*rk14
!     return
! 8   pdj(8) = -rk13*rk14 - rk10
!     pdj(10) = rk13*rk14
!     pdj(11) = rk10
! 9   return
! 10  pdj(3) = -rk7*y(3)
!     pdj(6) = rk7*y(3)
!     pdj(7) = rk17*y(12)
!     pdj(8) = rk9
!     pdj(10) = -rk7*y(3) - rk17*y(12) - rk6 - rk9
!     pdj(12) = rk6 - rk17*y(12)
! 11  return
! 12  pdj(2) = -rk15*y(2)
!     pdj(5) = rk15*y(2)
!     pdj(7) = rk17*y(10)
!     pdj(10) = -rk17*y(10)
!     pdj(12) = -rk15*y(2) - rk17*y(10)
!     return
!     end
!
! the output of this program (on a cray-1 in single precision)
! is as follows..
!
!
! at t =  1.000e-01     no. steps =   12     last step =  1.515e-02
!  y array =     9.90050e-01   6.28228e-03   3.65313e-03   7.51934e-07
!                1.12167e-09   1.18458e-09   1.77291e-12   3.26476e-07
!                5.46720e-08   9.99500e-06   4.48483e-08   2.76398e-06
!
!
! at t =  1.000e+00     no. steps =   33     last step =  7.880e-02
!  y array =     9.04837e-01   9.13105e-03   8.20622e-02   2.49177e-05
!                1.85055e-06   1.96797e-06   1.46157e-07   2.39557e-05
!                3.26306e-05   7.21621e-04   5.06433e-05   3.05010e-03
!
!
! at t =  1.000e+01     no. steps =   48     last step =  1.239e+00
!  y array =     3.67876e-01   3.68958e-03   3.65133e-01   4.48325e-05
!                6.10798e-05   4.33148e-05   5.90211e-05   1.18449e-04
!                3.15235e-03   3.56531e-03   4.15520e-03   2.48741e-01
!
!
! at t =  1.000e+02     no. steps =   91     last step =  3.764e+00
!  y array =     4.44981e-05   4.42666e-07   4.47273e-04  -3.53257e-11
!                2.81577e-08  -9.67741e-11   2.77615e-07   1.45322e-07
!                1.56230e-02   4.37394e-06   1.60104e-02   9.52246e-01
!
!
! at t =  1.000e+03     no. steps =  111     last step =  4.156e+02
!  y array =    -2.65492e-13   2.60539e-14  -8.59563e-12   6.29355e-14
!               -1.78066e-13   5.71471e-13  -1.47561e-12   4.58078e-15
!                1.56314e-02   1.37878e-13   1.60184e-02   9.52719e-01
!
!
! required rwork size = 442   iwork size =  30
! no. steps = 111   no. f-s = 142   no. j-s =   2   no. lu-s =  20
! no. of nonzeros in j =   44   no. of nonzeros in lu =   50
!-----------------------------------------------------------------------
! full description of user interface to lsodes.
!
! the user interface to lsodes consists of the following parts.
!
! i.   the call sequence to subroutine lsodes, which is a driver
!      routine for the solver.  this includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      following these descriptions is a description of
!      optional inputs available through the call sequence, and then
!      a description of optional outputs (in the work arrays).
!
! ii.  descriptions of other routines in the lsodes package that may be
!      (optionally) called by the user.  these provide the ability to
!      alter error message handling, save and restore the internal
!      common, and obtain specified derivatives of the solution y(t).
!
! iii. descriptions of common blocks to be declared in overlay
!      or similar environments, or to be saved when doing an interrupt
!      of the problem and continued solution later.
!
! iv.  description of two routines in the lsodes package, either of
!      which the user may replace with his own version, if desired.
!      these relate to the measurement of errors.
!
!-----------------------------------------------------------------------
! part i.  call sequence.
!
! the call sequence parameters used for input only are
!     f, neq, tout, itol, rtol, atol, itask, iopt, lrw, liw, jac, mf,
! and those used for both input and output are
!     y, t, istate.
! the work arrays rwork and iwork are also used for conditional and
! optional inputs and optional outputs.  (the term output here refers
! to the return from subroutine lsodes to the user-s calling program.)
!
! the legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by istate = 3 on input.
!
! the descriptions of the call arguments are as follows.
!
! f      = the name of the user-supplied subroutine defining the
!          ode system.  the system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y.  subroutine f is to
!          compute the function f.  it is to have the form
!               subroutine f (neq, t, y, ydot)
!               dimension y(1), ydot(1)
!          where neq, t, and y are input, and the array ydot = f(t,y)
!          is output.  y and ydot are arrays of length neq.
!          (in the dimension statement above, 1 is a dummy
!          dimension.. it can be replaced by any value.)
!          subroutine f should not alter y(1),...,y(neq).
!          f must be declared external in the calling program.
!
!          subroutine f may access user-defined quantities in
!          neq(2),... and/or in y(neq(1)+1),... if neq is an array
!          (dimensioned in f) and/or y has length exceeding neq(1).
!          see the descriptions of neq and y below.
!
!          if quantities computed in the f routine are needed
!          externally to lsodes, an extra call to f should be made
!          for this purpose, for consistent and accurate results.
!          if only the derivative dy/dt is needed, use intdy instead.
!
! neq    = the size of the ode system (number of first order
!          ordinary differential equations).  used only for input.
!          neq may be decreased, but not increased, during the problem.
!          if neq is decreased (with istate = 3 on input), the
!          remaining components of y should be left undisturbed, if
!          these are to be accessed in f and/or jac.
!
!          normally, neq is a scalar, and it is generally referred to
!          as a scalar in this user interface description.  however,
!          neq may be an array, with neq(1) set to the system size.
!          (the lsodes package accesses only neq(1).)  in either case,
!          this parameter is passed as the neq argument in all calls
!          to f and jac.  hence, if it is an array, locations
!          neq(2),... may be used to store other integer data and pass
!          it to f and/or jac.  subroutines f and/or jac must include
!          neq in a dimension statement in that case.
!
! y      = a real array for the vector of dependent variables, of
!          length neq or more.  used for both input and output on the
!          first call (istate = 1), and only for output on other calls.
!          on the first call, y must contain the vector of initial
!          values.  on output, y contains the computed solution vector,
!          evaluated at t.  if desired, the y array may be used
!          for other purposes between calls to the solver.
!
!          this array is passed as the y argument in all calls to
!          f and jac.  hence its length may exceed neq, and locations
!          y(neq+1),... may be used to store other real data and
!          pass it to f and/or jac.  (the lsodes package accesses only
!          y(1),...,y(neq).)
!
! t      = the independent variable.  on input, t is used only on the
!          first call, as the initial point of the integration.
!          on output, after each call, t is the value at which a
!          computed solution y is evaluated (usually the same as tout).
!          on an error return, t is the farthest point reached.
!
! tout   = the next value of t at which a computed solution is desired.
!          used only for input.
!
!          when starting the problem (istate = 1), tout may be equal
!          to t for one call, then should .ne. t for the next call.
!          for the initial t, an input value of tout .ne. t is used
!          in order to determine the direction of the integration
!          (i.e. the algebraic sign of the step sizes) and the rough
!          scale of the problem.  integration in either direction
!          (forward or backward in t) is permitted.
!
!          if itask = 2 or 5 (one-step modes), tout is ignored after
!          the first call (i.e. the first call with tout .ne. t).
!          otherwise, tout is required on every call.
!
!          if itask = 1, 3, or 4, the values of tout need not be
!          monotone, but a value of tout which backs up is limited
!          to the current internal t interval, whose endpoints are
!          tcur - hu and tcur (see optional outputs, below, for
!          tcur and hu).
!
! itol   = an indicator for the type of error control.  see
!          description below under atol.  used only for input.
!
! rtol   = a relative error tolerance parameter, either a scalar or
!          an array of length neq.  see description below under atol.
!          input only.
!
! atol   = an absolute error tolerance parameter, either a scalar or
!          an array of length neq.  input only.
!
!             the input parameters itol, rtol, and atol determine
!          the error control performed by the solver.  the solver will
!          control the vector e = (e(i)) of estimated local errors
!          in y, according to an inequality of the form
!                      rms-norm of ( e(i)/ewt(i) )   .le.   1,
!          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),
!          and the rms-norm (root-mean-square norm) here is
!          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))
!          is a vector of weights which must always be positive, and
!          the values of rtol and atol should all be non-negative.
!          the following table gives the types (scalar/array) of
!          rtol and atol, and the corresponding form of ewt(i).
!
!             itol    rtol       atol          ewt(i)
!              1     scalar     scalar     rtol*abs(y(i)) + atol
!              2     scalar     array      rtol*abs(y(i)) + atol(i)
!              3     array      scalar     rtol(i)*abs(y(i)) + atol
!              4     array      array      rtol(i)*abs(y(i)) + atol(i)
!
!          when either of these parameters is a scalar, it need not
!          be dimensioned in the user-s calling program.
!
!          if none of the above choices (with itol, rtol, and atol
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting
!          user-supplied routines for the setting of ewt and/or for
!          the norm calculation.  see part iv below.
!
!          if global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of rtol and atol (i.e. of ewt) should be scaled
!          down uniformly.
!
! itask  = an index specifying the task to be performed.
!          input only.  itask has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = tout (by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = tout and return.
!          4  means normal computation of output values of y(t) at
!             t = tout but without overshooting t = tcrit.
!             tcrit must be input as rwork(1).  tcrit may be equal to
!             or beyond tout, but not behind it in the direction of
!             integration.  this option is useful if the problem
!             has a singularity at or beyond t = tcrit.
!          5  means take one step, without passing tcrit, and return.
!             tcrit must be input as rwork(1).
!
!          note..  if itask = 4 or 5 and the solver reaches tcrit
!          (within roundoff), it will return t = tcrit (exactly) to
!          indicate this (unless itask = 4 and tout comes before tcrit,
!          in which case answers at t = tout are returned first).
!
! istate = an index used for input and output to specify the
!          the state of the calculation.
!
!          on input, the values of istate are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done).  see note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly tout and itask.
!             (if itol, rtol, and/or atol are changed between calls
!             with istate = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             tout and itask.  changes are allowed in
!             neq, itol, rtol, atol, iopt, lrw, liw, mf,
!             the conditional inputs ia and ja,
!             and any of the optional inputs except h0.
!             in particular, if miter = 1 or 2, a call with istate = 3
!             will cause the sparsity structure of the problem to be
!             recomputed (or reread from ia and ja if moss = 0).
!          note..  a preliminary call with tout = t is not counted
!          as a first call here, as no initialization or checking of
!          input is done.  (such a call is sometimes useful for the
!          purpose of outputting the initial conditions.)
!          thus the first call for which tout .ne. t requires
!          istate = 1 on input.
!
!          on output, istate has the following values and meanings.
!           1  means nothing was done, as tout was equal to t with
!              istate = 1 on input.  (however, an internal counter was
!              set to detect and prevent repeated calls of this type.)
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than mxstep
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as t.  (mxstep is an optional input
!              and is normally 500.)  to continue, the user may
!              simply reset istate to a value .gt. 1 and call again
!              (the excess work step counter will be reset to 0).
!              in addition, the user may increase mxstep to avoid
!              this error return (see below on optional inputs).
!          -2  means too much accuracy was requested for the precision
!              of the machine being used.  this was detected before
!              completing the requested task, but the integration
!              was successful as far as t.  to continue, the tolerance
!              parameters must be reset, and istate must be set
!              to 3.  the optional output tolsf may be used for this
!              purpose.  (note.. if this condition is detected before
!              taking any steps, then an illegal input return
!              (istate = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps.  see written message for details.
!              note..  if the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as t.
!              the problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as t.
!              this may be caused by an inaccurate jacobian matrix,
!              if one is being used.
!          -6  means ewt(i) became zero for some i during the
!              integration.  pure relative error control (atol(i)=0.0)
!              was requested on a variable which has now vanished.
!              the integration was successful as far as t.
!          -7  means a fatal error return flag came from the sparse
!              solver cdrv by way of prjs or slss (numerical
!              factorization or backsolve).  this should never happen.
!              the integration was successful as far as t.
!
!          note.. an error return with istate = -1, -4, or -5 and with
!          miter = 1 or 2 may mean that the sparsity structure of the
!          problem has changed significantly since it was last
!          determined (or input).  in that case, one can attempt to
!          complete the integration by setting istate = 3 on the next
!          call, so that a new structure determination is done.
!
!          note..  since the normal output value of istate is 2,
!          it does not need to be reset for normal continuation.
!          also, since a negative input value of istate will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other inputs, before
!          calling the solver again.
!
! iopt   = an integer flag to specify whether or not any optional
!          inputs are being used on this call.  input only.
!          the optional inputs are listed separately below.
!          iopt = 0 means no optional inputs are being used.
!                   default values will be used in all cases.
!          iopt = 1 means one or more optional inputs are being used.
!
! rwork  = a work array used for a mixture of real (double precision)
!          and integer work space.
!          the length of rwork (in real words) must be at least
!             20 + nyh*(maxord + 1) + 3*neq + lwm    where
!          nyh    = the initial value of neq,
!          maxord = 12 (if meth = 1) or 5 (if meth = 2) (unless a
!                   smaller value is given as an optional input),
!          lwm = 0                                    if miter = 0,
!          lwm = 2*nnz + 2*neq + (nnz+9*neq)/lenrat   if miter = 1,
!          lwm = 2*nnz + 2*neq + (nnz+10*neq)/lenrat  if miter = 2,
!          lwm = neq + 2                              if miter = 3.
!          in the above formulas,
!          nnz    = number of nonzero elements in the jacobian matrix.
!          lenrat = the real to integer wordlength ratio (usually 1 in
!                   single precision and 2 in double precision).
!          (see the mf description for meth and miter.)
!          thus if maxord has its default value and neq is constant,
!          the minimum length of rwork is..
!             20 + 16*neq        for mf = 10,
!             20 + 16*neq + lwm  for mf = 11, 111, 211, 12, 112, 212,
!             22 + 17*neq        for mf = 13,
!             20 +  9*neq        for mf = 20,
!             20 +  9*neq + lwm  for mf = 21, 121, 221, 22, 122, 222,
!             22 + 10*neq        for mf = 23.
!          if miter = 1 or 2, the above formula for lwm is only a
!          crude lower bound.  the required length of rwork cannot
!          be readily predicted in general, as it depends on the
!          sparsity structure of the problem.  some experimentation
!          may be necessary.
!
!          the first 20 words of rwork are reserved for conditional
!          and optional inputs and optional outputs.
!
!          the following word in rwork is a conditional input..
!            rwork(1) = tcrit = critical value of t which the solver
!                       is not to overshoot.  required if itask is
!                       4 or 5, and ignored otherwise.  (see itask.)
!
! lrw    = the length of the array rwork, as declared by the user.
!          (this will be checked by the solver.)
!
! iwork  = an integer work array.  the length of iwork must be at least
!             31 + neq + nnz   if moss = 0 and miter = 1 or 2, or
!             30               otherwise.
!          (nnz is the number of nonzero elements in df/dy.)
!
!          in lsodes, iwork is used only for conditional and
!          optional inputs and optional outputs.
!
!          the following two blocks of words in iwork are conditional
!          inputs, required if moss = 0 and miter = 1 or 2, but not
!          otherwise (see the description of mf for moss).
!            iwork(30+j) = ia(j)     (j=1,...,neq+1)
!            iwork(31+neq+k) = ja(k) (k=1,...,nnz)
!          the two arrays ia and ja describe the sparsity structure
!          to be assumed for the jacobian matrix.  ja contains the row
!          indices where nonzero elements occur, reading in columnwise
!          order, and ia contains the starting locations in ja of the
!          descriptions of columns 1,...,neq, in that order, with
!          ia(1) = 1.  thus, for each column index j = 1,...,neq, the
!          values of the row index i in column j where a nonzero
!          element may occur are given by
!            i = ja(k),  where   ia(j) .le. k .lt. ia(j+1).
!          if nnz is the total number of nonzero locations assumed,
!          then the length of the ja array is nnz, and ia(neq+1) must
!          be nnz + 1.  duplicate entries are not allowed.
!
! liw    = the length of the array iwork, as declared by the user.
!          (this will be checked by the solver.)
!
! note..  the work arrays must not be altered between calls to lsodes
! for the same problem, except possibly for the conditional and
! optional inputs, and except for the last 3*neq words of rwork.
! the latter space is used for internal scratch space, and so is
! available for use by the user outside lsodes between calls, if
! desired (but not for use by f or jac).
!
! jac    = name of user-supplied routine (miter = 1 or moss = 1) to
!          compute the jacobian matrix, df/dy, as a function of
!          the scalar t and the vector y.  it is to have the form
!               subroutine jac (neq, t, y, j, ian, jan, pdj)
!               dimension y(1), ian(1), jan(1), pdj(1)
!          where neq, t, y, j, ian, and jan are input, and the array
!          pdj, of length neq, is to be loaded with column j
!          of the jacobian on output.  thus df(i)/dy(j) is to be
!          loaded into pdj(i) for all relevant values of i.
!          here t and y have the same meaning as in subroutine f,
!          and j is a column index (1 to neq).  ian and jan are
!          undefined in calls to jac for structure determination
!          (moss = 1).  otherwise, ian and jan are structure
!          descriptors, as defined under optional outputs below, and
!          so can be used to determine the relevant row indices i, if
!          desired.  (in the dimension statement above, 1 is a
!          dummy dimension.. it can be replaced by any value.)
!               jac need not provide df/dy exactly.  a crude
!          approximation (possibly with greater sparsity) will do.
!               in any case, pdj is preset to zero by the solver,
!          so that only the nonzero elements need be loaded by jac.
!          calls to jac are made with j = 1,...,neq, in that order, and
!          each such set of calls is preceded by a call to f with the
!          same arguments neq, t, and y.  thus to gain some efficiency,
!          intermediate quantities shared by both calculations may be
!          saved in a user common block by f and not recomputed by jac,
!          if desired.  jac must not alter its input arguments.
!          jac must be declared external in the calling program.
!               subroutine jac may access user-defined quantities in
!          neq(2),... and y(neq(1)+1),... if neq is an array
!          (dimensioned in jac) and y has length exceeding neq(1).
!          see the descriptions of neq and y above.
!
! mf     = the method flag.  used only for input.
!          mf has three decimal digits-- moss, meth, miter--
!             mf = 100*moss + 10*meth + miter.
!          moss indicates the method to be used to obtain the sparsity
!          structure of the jacobian matrix if miter = 1 or 2..
!            moss = 0 means the user has supplied ia and ja
!                     (see descriptions under iwork above).
!            moss = 1 means the user has supplied jac (see below)
!                     and the structure will be obtained from neq
!                     initial calls to jac.
!            moss = 2 means the structure will be obtained from neq+1
!                     initial calls to f.
!          meth indicates the basic linear multistep method..
!            meth = 1 means the implicit adams method.
!            meth = 2 means the method based on backward
!                     differentiation formulas (bdf-s).
!          miter indicates the corrector iteration method..
!            miter = 0 means functional iteration (no jacobian matrix
!                      is involved).
!            miter = 1 means chord iteration with a user-supplied
!                      sparse jacobian, given by subroutine jac.
!            miter = 2 means chord iteration with an internally
!                      generated (difference quotient) sparse jacobian
!                      (using ngp extra calls to f per df/dy value,
!                      where ngp is an optional output described below.)
!            miter = 3 means chord iteration with an internally
!                      generated diagonal jacobian approximation.
!                      (using 1 extra call to f per df/dy evaluation).
!          if miter = 1 or moss = 1, the user must supply a subroutine
!          jac (the name is arbitrary) as described above under jac.
!          otherwise, a dummy argument can be used.
!
!          the standard choices for mf are..
!            mf = 10  for a nonstiff problem,
!            mf = 21 or 22 for a stiff problem with ia/ja supplied
!                     (21 if jac is supplied, 22 if not),
!            mf = 121 for a stiff problem with jac supplied,
!                     but not ia/ja,
!            mf = 222 for a stiff problem with neither ia/ja nor
!                     jac supplied.
!          the sparseness structure can be changed during the
!          problem by making a call to lsodes with istate = 3.
!-----------------------------------------------------------------------
! optional inputs.
!
! the following is a list of the optional inputs provided for in the
! call sequence.  (see also part ii.)  for each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! the use of any of these inputs requires iopt = 1, and in that
! case all of these inputs are examined.  a value of zero for any
! of these optional inputs will cause the default value to be used.
! thus to use a subset of the optional inputs, simply preload
! locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and
! then set those of interest to nonzero values.
!
! name    location      meaning and default value
!
! h0      rwork(5)  the step size to be attempted on the first step.
!                   the default value is determined by the solver.
!
! hmax    rwork(6)  the maximum absolute step size allowed.
!                   the default value is infinite.
!
! hmin    rwork(7)  the minimum absolute step size allowed.
!                   the default value is 0.  (this lower bound is not
!                   enforced on the final step before reaching tcrit
!                   when itask = 4 or 5.)
!
! seth    rwork(8)  the element threshhold for sparsity determination
!                   when moss = 1 or 2.  if the absolute value of
!                   an estimated jacobian element is .le. seth, it
!                   will be assumed to be absent in the structure.
!                   the default value of seth is 0.
!
! maxord  iwork(5)  the maximum order to be allowed.  the default
!                   value is 12 if meth = 1, and 5 if meth = 2.
!                   if maxord exceeds the default value, it will
!                   be reduced to the default value.
!                   if maxord is changed during the problem, it may
!                   cause the current order to be reduced.
!
! mxstep  iwork(6)  maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   the default value is 500.
!
! mxhnil  iwork(7)  maximum number of messages printed (per problem)
!                   warning that t + h = t on a step (h = step size).
!                   this must be positive to result in a non-default
!                   value.  the default value is 10.
!-----------------------------------------------------------------------
! optional outputs.
!
! as optional additional output from lsodes, the variables listed
! below are quantities related to the performance of lsodes
! which are available to the user.  these are communicated by way of
! the work arrays, but also have internal mnemonic names as shown.
! except where stated otherwise, all of these outputs are defined
! on any successful return from lsodes, and on any return with
! istate = -1, -2, -4, -5, or -6.  on an illegal input return
! (istate = -3), they will be unchanged from their existing values
! (if any), except possibly for tolsf, lenrw, and leniw.
! on any error return, outputs relevant to the error will be defined,
! as noted below.
!
! name    location      meaning
!
! hu      rwork(11) the step size in t last used (successfully).
!
! hcur    rwork(12) the step size to be attempted on the next step.
!
! tcur    rwork(13) the current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t.  on output, tcur
!                   will always be at least as far as the argument
!                   t, but may be farther (if interpolation was done).
!
! tolsf   rwork(14) a tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (istate = -3 if detected at the start of
!                   the problem, istate = -2 otherwise).  if itol is
!                   left unaltered but rtol and atol are uniformly
!                   scaled up by a factor of tolsf for the next call,
!                   then the solver is deemed likely to succeed.
!                   (the user may also ignore tolsf and alter the
!                   tolerance parameters in any other way appropriate.)
!
! nst     iwork(11) the number of steps taken for the problem so far.
!
! nfe     iwork(12) the number of f evaluations for the problem so far,
!                   excluding those for structure determination
!                   (moss = 2).
!
! nje     iwork(13) the number of jacobian evaluations for the problem
!                   so far, excluding those for structure determination
!                   (moss = 1).
!
! nqu     iwork(14) the method order last used (successfully).
!
! nqcur   iwork(15) the order to be attempted on the next step.
!
! imxer   iwork(16) the index of the component of largest magnitude in
!                   the weighted local error vector ( e(i)/ewt(i) ),
!                   on an error return with istate = -4 or -5.
!
! lenrw   iwork(17) the length of rwork actually required.
!                   this is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! leniw   iwork(18) the length of iwork actually required.
!                   this is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! nnz     iwork(19) the number of nonzero elements in the jacobian
!                   matrix, including the diagonal (miter = 1 or 2).
!                   (this may differ from that given by ia(neq+1)-1
!                   if moss = 0, because of added diagonal entries.)
!
! ngp     iwork(20) the number of groups of column indices, used in
!                   difference quotient jacobian aproximations if
!                   miter = 2.  this is also the number of extra f
!                   evaluations needed for each jacobian evaluation.
!
! nlu     iwork(21) the number of sparse lu decompositions for the
!                   problem so far.
!
! lyh     iwork(22) the base address in rwork of the history array yh,
!                   described below in this list.
!
! ipian   iwork(23) the base address of the structure descriptor array
!                   ian, described below in this list.
!
! ipjan   iwork(24) the base address of the structure descriptor array
!                   jan, described below in this list.
!
! nzl     iwork(25) the number of nonzero elements in the strict lower
!                   triangle of the lu factorization used in the chord
!                   iteration (miter = 1 or 2).
!
! nzu     iwork(26) the number of nonzero elements in the strict upper
!                   triangle of the lu factorization used in the chord
!                   iteration (miter = 1 or 2).
!                   the total number of nonzeros in the factorization
!                   is therefore nzl + nzu + neq.
!
! the following four arrays are segments of the rwork array which
! may also be of interest to the user as optional outputs.
! for each array, the table below gives its internal name,
! its base address, and its description.
! for yh and acor, the base addresses are in rwork (a real array).
! the integer arrays ian and jan are to be obtained by declaring an
! integer array iwk and identifying iwk(1) with rwork(21), using either
! an equivalence statement or a subroutine call.  then the base
! addresses ipian (of ian) and ipjan (of jan) in iwk are to be obtained
! as optional outputs iwork(23) and iwork(24), respectively.
! thus ian(1) is iwk(ipian), etc.
!
! name    base address      description
!
! ian    ipian (in iwk)  structure descriptor array of size neq + 1.
! jan    ipjan (in iwk)  structure descriptor array of size nnz.
!         (see above)    ian and jan together describe the sparsity
!                        structure of the jacobian matrix, as used by
!                        lsodes when miter = 1 or 2.
!                        jan contains the row indices of the nonzero
!                        locations, reading in columnwise order, and
!                        ian contains the starting locations in jan of
!                        the descriptions of columns 1,...,neq, in
!                        that order, with ian(1) = 1.  thus for each
!                        j = 1,...,neq, the row indices i of the
!                        nonzero locations in column j are
!                        i = jan(k),  ian(j) .le. k .lt. ian(j+1).
!                        note that ian(neq+1) = nnz + 1.
!                        (if moss = 0, ian/jan may differ from the
!                        input ia/ja because of a different ordering
!                        in each column, and added diagonal entries.)
!
! yh      lyh            the nordsieck history array, of size nyh by
!          (optional     (nqcur + 1), where nyh is the initial value
!          output)       of neq.  for j = 0,1,...,nqcur, column j+1
!                        of yh contains hcur**j/factorial(j) times
!                        the j-th derivative of the interpolating
!                        polynomial currently representing the solution,
!                        evaluated at t = tcur.  the base address lyh
!                        is another optional output, listed above.
!
! acor     lenrw-neq+1   array of size neq used for the accumulated
!                        corrections on each step, scaled on output
!                        to represent the estimated local error in y
!                        on the last step.  this is the vector e in
!                        the description of the error control.  it is
!                        defined only on a successful return from
!                        lsodes.
!
!-----------------------------------------------------------------------
! part ii.  other routines callable.
!
! the following are optional calls which the user may make to
! gain additional capabilities in conjunction with lsodes.
! (the routines xsetun and xsetf are designed to conform to the
! slatec error handling package.)
!
!     form of call                  function
!   call xsetun(lun)          set the logical unit number, lun, for
!                             output of messages from lsodes, if
!                             the default is not desired.
!                             the default value of lun is 6.
!
!   call xsetf(mflag)         set a flag to control the printing of
!                             messages by lsodes.
!                             mflag = 0 means do not print. (danger..
!                             this risks losing valuable information.)
!                             mflag = 1 means print (the default).
!
!                             either of the above calls may be made at
!                             any time and will take effect immediately.
!
!   call srcms(rsav,isav,job) saves and restores the contents of
!                             the internal common blocks used by
!                             lsodes (see part iii below).
!                             rsav must be a real array of length 224
!                             or more, and isav must be an integer
!                             array of length 75 or more.
!                             job=1 means save common into rsav/isav.
!                             job=2 means restore common from rsav/isav.
!                                srcms is useful if one is
!                             interrupting a run and restarting
!                             later, or alternating between two or
!                             more problems solved with lsodes.
!
!   call intdy(,,,,,)         provide derivatives of y, of various
!        (see below)          orders, at a specified point t, if
!                             desired.  it may be called only after
!                             a successful return from lsodes.
!
! the detailed instructions for using intdy are as follows.
! the form of the call is..
!
!   lyh = iwork(22)
!   call intdy (t, k, rwork(lyh), nyh, dky, iflag)
!
! the input parameters are..
!
! t         = value of independent variable where answers are desired
!             (normally the same as the t last returned by lsodes).
!             for valid results, t must lie between tcur - hu and tcur.
!             (see optional outputs for tcur and hu.)
! k         = integer order of the derivative desired.  k must satisfy
!             0 .le. k .le. nqcur, where nqcur is the current order
!             (see optional outputs).  the capability corresponding
!             to k = 0, i.e. computing y(t), is already provided
!             by lsodes directly.  since nqcur .ge. 1, the first
!             derivative dy/dt is always available with intdy.
! lyh       = the base address of the history array yh, obtained
!             as an optional output as shown above.
! nyh       = column length of yh, equal to the initial value of neq.
!
! the output parameters are..
!
! dky       = a real array of length neq containing the computed value
!             of the k-th derivative of y(t).
! iflag     = integer flag, returned as 0 if k and t were legal,
!             -1 if k was illegal, and -2 if t was illegal.
!             on an error return, a message is also written.
!-----------------------------------------------------------------------
! part iii.  common blocks.
!
! if lsodes is to be used in an overlay situation, the user
! must declare, in the primary overlay, the variables in..
!   (1) the call sequence to lsodes,
!   (2) the three internal common blocks
!         /ls0001/  of length  257  (218 double precision words
!                         followed by 39 integer words),
!         /lss001/  of length  40    ( 6 double precision words
!                         followed by 34 integer words),
!         /eh0001/  of length  2 (integer words).
!
! if lsodes is used on a system in which the contents of internal
! common blocks are not preserved between calls, the user should
! declare the above three common blocks in his main program to insure
! that their contents are preserved.
!
! if the solution of a given problem by lsodes is to be interrupted
! and then later continued, such as when restarting an interrupted run
! or alternating between two or more problems, the user should save,
! following the return from the last lsodes call prior to the
! interruption, the contents of the call sequence variables and the
! internal common blocks, and later restore these values before the
! next lsodes call for that problem.  to save and restore the common
! blocks, use subroutine srcms (see part ii above).
!
!-----------------------------------------------------------------------
! part iv.  optionally replaceable solver routines.
!
! below are descriptions of two routines in the lsodes package which
! relate to the measurement of errors.  either routine can be
! replaced by a user-supplied version, if desired.  however, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (note.. the means by which the package version of a routine is
! superseded by the user-s version may be system-dependent.)
!
! (a) ewset.
! the following subroutine is called just before each internal
! integration step, and sets the array of error weights, ewt, as
! described under itol/rtol/atol above..
!     subroutine ewset (neq, itol, rtol, atol, ycur, ewt)
! where neq, itol, rtol, and atol are as in the lsodes call sequence,
! ycur contains the current dependent variable vector, and
! ewt is the array of weights set by ewset.
!
! if the user supplies this subroutine, it must return in ewt(i)
! (i = 1,...,neq) a positive quantity suitable for comparing errors
! in y(i) to.  the ewt array returned by ewset is passed to the
! vnorm routine (see below), and also used by lsodes in the computation
! of the optional output imxer, the diagonal jacobian approximation,
! and the increments for difference quotient jacobians.
!
! in the user-supplied version of ewset, it may be desirable to use
! the current values of derivatives of y.  derivatives up to order nq
! are available from the history array yh, described above under
! optional outputs.  in ewset, yh is identical to the ycur array,
! extended to nq + 1 columns with a column length of nyh and scale
! factors of h**j/factorial(j).  on the first call for the problem,
! given by nst = 0, nq is 1 and h is temporarily set to 1.0.
! the quantities nq, nyh, h, and nst can be obtained by including
! in ewset the statements..
!     double precision h, rls
!     common /ls0001/ rls(218),ils(39)
!     nq = ils(35)
!     nyh = ils(14)
!     nst = ils(36)
!     h = rls(212)
! thus, for example, the current value of dy/dt can be obtained as
! ycur(nyh+i)/h  (i=1,...,neq)  (and the division by h is
! unnecessary when nst = 0).
!
! (b) vnorm.
! the following is a real function routine which computes the weighted
! root-mean-square norm of a vector v..
!     d = vnorm (n, v, w)
! where..
!   n = the length of the vector,
!   v = real array of length n containing the vector,
!   w = real array of length n containing weights,
!   d = sqrt( (1/n) * sum(v(i)*w(i))**2 ).
! vnorm is called with n = neq and with w(i) = 1.0/ewt(i), where
! ewt is as set by subroutine ewset.
!
! if the user supplies this function, it should return a non-negative
! value of vnorm suitable for use in the error control in lsodes.
! none of the arguments should be altered by vnorm.
! for example, a user-supplied vnorm routine might..
!   -substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
!   -ignore some components of v in the norm, with the effect of
!    suppressing the error control on those components of y.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! other routines in the lsodes package.
!
! in addition to subroutine lsodes, the lsodes package includes the
! following subroutines and function routines..
!  iprep    acts as an iterface between lsodes and prep, and also does
!           adjusting of work space pointers and work arrays.
!  prep     is called by iprep to compute sparsity and do sparse matrix
!           preprocessing if miter = 1 or 2.
!  jgroup   is called by prep to compute groups of jacobian column
!           indices for use when miter = 2.
!  adjlr    adjusts the length of required sparse matrix work space.
!           it is called by prep.
!  cntnzu   is called by prep and counts the nonzero elements in the
!           strict upper triangle of j + j-transpose, where j = df/dy.
!  intdy    computes an interpolated value of the y vector at t = tout.
!  stode    is the core integrator, which does one step of the
!           integration and the associated error control.
!  cfode    sets all method coefficients and test constants.
!  prjs     computes and preprocesses the jacobian matrix j = df/dy
!           and the newton iteration matrix p = i - h*l0*j.
!  slss     manages solution of linear system in chord iteration.
!  ewset    sets the error weight vector ewt before each step.
!  vnorm    computes the weighted r.m.s. norm of a vector.
!  srcms    is a user-callable routine to save and restore
!           the contents of the internal common blocks.
!  odrv     constructs a reordering of the rows and columns of
!           a matrix by the minimum degree algorithm.  odrv is a
!           driver routine which calls subroutines md, mdi, mdm,
!           mdp, mdu, and sro.  see ref. 2 for details.  (the odrv
!           module has been modified since ref. 2, however.)
!  cdrv     performs reordering, symbolic factorization, numerical
!           factorization, or linear system solution operations,
!           depending on a path argument ipath.  cdrv is a
!           driver routine which calls subroutines nroc, nsfc,
!           nnfc, nnsc, and nntc.  see ref. 3 for details.
!           lsodes uses cdrv to solve linear systems in which the
!           coefficient matrix is  p = i - con*j, where i is the
!           identity, con is a scalar, and j is an approximation to
!           the jacobian df/dy.  because cdrv deals with rowwise
!           sparsity descriptions, cdrv works with p-transpose, not p.
!  d1mach   computes the unit roundoff in a machine-independent manner.
!  xerrwv, xsetun, and xsetf   handle the printing of all error
!           messages and warnings.  xerrwv is machine-dependent.
! note..  vnorm and d1mach are function routines.
! all the others are subroutines.
!
! the intrinsic and external routines used by lsodes are..
! dabs, dmax1, dmin1, dfloat, max0, min0, mod, dsign, dsqrt, and write.
!
! a block data subprogram is also included with the package,
! for loading some of the variables in internal common.
!
!-----------------------------------------------------------------------
! the following card is for optimized compilation on lll compilers.
!lll. optimize
!-----------------------------------------------------------------------
!      EXTERNAL prjs, slss  (not if prjs and slss are in this module)
      INTEGER illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, &
         mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      INTEGER icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
      INTEGER i, i1, i2, iflag, imax, imul, imxer, ipflag, ipgo, irem, &
         j, kgo, lenrat, lenyht, leniw, lenrw, lf0, lia, lja, &
         lrtem, lwtem, lyhd, lyhn, mf1, mord, mxhnl0, mxstp0, ncolm
      REAL(DP) rowns, &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      REAL(DP) con0, conmin, ccmxj, psmall, rbig, seth
      REAL(DP) atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli, &
         tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0, &
         d1mach
      DIMENSION mord(2)
      LOGICAL ihit
!-----------------------------------------------------------------------
! the following two internal common blocks contain
! (a) variables which are local to any subroutine but whose values must
!     be preserved between calls to the routine (own variables), and
! (b) variables which are communicated between subroutines.
! the structure of each block is as follows..  all real variables are
! listed first, followed by all integers.  within each type, the
! variables are grouped with those local to subroutine lsodes first,
! then those local to subroutine stode or subroutine prjs
! (no other routines have own variables), and finally those used
! for communication.  the block ls0001 is declared in subroutines
! lsodes, iprep, prep, intdy, stode, prjs, and slss.  the block lss001
! is declared in subroutines lsodes, iprep, prep, prjs, and slss.
! groups of variables are replaced by dummy arrays in the common
! declarations in routines where those variables are not used.
!-----------------------------------------------------------------------
      COMMON /ls0001/ rowns(209), &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, &
         illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, &
         mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6), &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
!
      COMMON /lss001/ con0, conmin, ccmxj, psmall, rbig, seth, &
         iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
!
      DATA mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/
!-----------------------------------------------------------------------
! in the data statement below, set lenrat equal to the ratio of
! the wordlength for a real number to that for an integer.  usually,
! lenrat = 1 for single precision and 2 for double precision.  if the
! true ratio is not an integer, use the next smaller integer (.ge. 1).
!-----------------------------------------------------------------------
      DATA lenrat/2/
!-----------------------------------------------------------------------
! block a.
! this code block is executed on every call.
! it tests istate and itask for legality and branches appropriately.
! if istate .gt. 1 but the flag init shows that initialization has
! not yet been done, an error return occurs.
! if istate = 1 and tout = t, jump to block g and return immediately.
!-----------------------------------------------------------------------
      IF (istate .LT. 1 .OR. istate .GT. 3) go to 601
      IF (itask .LT. 1 .OR. itask .GT. 5) go to 602
      IF (istate .EQ. 1) go to 10
      IF (init .EQ. 0) go to 603
      IF (istate .EQ. 2) go to 200
      go to 20
 10   init = 0
      IF (tout .EQ. t) go to 430
 20   ntrep = 0
!-----------------------------------------------------------------------
! block b.
! the next code block is executed for the initial call (istate = 1),
! or for a continuation call with parameter changes (istate = 3).
! it contains checking of all inputs and various initializations.
! if istate = 1, the final setting of work space pointers, the matrix
! preprocessing, and other initializations are done in block c.
!
! first check legality of the non-optional inputs neq, itol, iopt,
! mf, ml, and mu.
!-----------------------------------------------------------------------
      IF (neq(1) .LE. 0) go to 604
      IF (istate .EQ. 1) go to 25
      IF (neq(1) .GT. n) go to 605
 25   n = neq(1)
      IF (itol .LT. 1 .OR. itol .GT. 4) go to 606
      IF (iopt .LT. 0 .OR. iopt .GT. 1) go to 607
      moss = mf/100
      mf1 = mf - 100*moss
      meth = mf1/10
      miter = mf1 - 10*meth
      IF (moss .LT. 0 .OR. moss .GT. 2) go to 608
      IF (meth .LT. 1 .OR. meth .GT. 2) go to 608
      IF (miter .LT. 0 .OR. miter .GT. 3) go to 608
      IF (miter .EQ. 0 .OR. miter .EQ. 3) moss = 0
! next process and check the optional inputs. --------------------------
      IF (iopt .EQ. 1) go to 40
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      IF (istate .EQ. 1) h0 = 0.0d0
      hmxi = 0.0d0
      hmin = 0.0d0
      seth = 0.0d0
      go to 60
 40   maxord = iwork(5)
      IF (maxord .LT. 0) go to 611
      IF (maxord .EQ. 0) maxord = 100
      maxord = min0(maxord,mord(meth))
      mxstep = iwork(6)
      IF (mxstep .LT. 0) go to 612
      IF (mxstep .EQ. 0) mxstep = mxstp0
      mxhnil = iwork(7)
      IF (mxhnil .LT. 0) go to 613
      IF (mxhnil .EQ. 0) mxhnil = mxhnl0
      IF (istate .NE. 1) go to 50
      h0 = rwork(5)
      IF ((tout - t)*h0 .LT. 0.0d0) go to 614
 50   hmax = rwork(6)
      IF (hmax .LT. 0.0d0) go to 615
      hmxi = 0.0d0
      IF (hmax .GT. 0.0d0) hmxi = 1.0d0/hmax
      hmin = rwork(7)
      IF (hmin .LT. 0.0d0) go to 616
      seth = rwork(8)
      IF (seth .LT. 0.0d0) go to 609
! check rtol and atol for legality. ------------------------------------
 60   rtoli = rtol(1)
      atoli = atol(1)
      DO 65 i = 1,n
        IF (itol .GE. 3) rtoli = rtol(i)
        IF (itol .EQ. 2 .OR. itol .EQ. 4) atoli = atol(i)
        IF (rtoli .LT. 0.0d0) go to 619
        IF (atoli .LT. 0.0d0) go to 620
 65     CONTINUE
!-----------------------------------------------------------------------
! compute required work array lengths, as far as possible, and test
! these against lrw and liw.  then set tentative pointers for work
! arrays.  pointers to rwork/iwork segments are named by prefixing l to
! the name of the segment.  e.g., the segment yh starts at rwork(lyh).
! segments of rwork (in order) are denoted  wm, yh, savf, ewt, acor.
! if miter = 1 or 2, the required length of the matrix work space wm
! is not yet known, and so a crude minimum value is used for the
! initial tests of lrw and liw, and yh is temporarily stored as far
! to the right in rwork as possible, to leave the maximum amount
! of space for wm for matrix preprocessing.  thus if miter = 1 or 2
! and moss .ne. 2, some of the segments of rwork are temporarily
! omitted, as they are not needed in the preprocessing.  these
! omitted segments are.. acor if istate = 1, ewt and acor if istate = 3
! and moss = 1, and savf, ewt, and acor if istate = 3 and moss = 0.
!-----------------------------------------------------------------------
      lrat = lenrat
      IF (istate .EQ. 1) nyh = n
      lwmin = 0
      IF (miter .EQ. 1) lwmin = 4*n + 10*n/lrat
      IF (miter .EQ. 2) lwmin = 4*n + 11*n/lrat
      IF (miter .EQ. 3) lwmin = n + 2
      lenyh = (maxord+1)*nyh
      lrest = lenyh + 3*n
      lenrw = 20 + lwmin + lrest
      iwork(17) = lenrw
      leniw = 30
      IF (moss .EQ. 0 .AND. miter .NE. 0 .AND. miter .NE. 3) &
         leniw = leniw + n + 1
      iwork(18) = leniw
      IF (lenrw .GT. lrw) go to 617
      IF (leniw .GT. liw) go to 618
      lia = 31
      IF (moss .EQ. 0 .AND. miter .NE. 0 .AND. miter .NE. 3) &
         leniw = leniw + iwork(lia+n) - 1
      iwork(18) = leniw
      IF (leniw .GT. liw) go to 618
      lja = lia + n + 1
      lia = min0(lia,liw)
      lja = min0(lja,liw)
      lwm = 21
      IF (istate .EQ. 1) nq = 1
      ncolm = min0(nq+1,maxord+2)
      lenyhm = ncolm*nyh
      lenyht = lenyh
      IF (miter .EQ. 1 .OR. miter .EQ. 2) lenyht = lenyhm
      imul = 2
      IF (istate .EQ. 3) imul = moss
      IF (moss .EQ. 2) imul = 3
      lrtem = lenyht + imul*n
      lwtem = lwmin
      IF (miter .EQ. 1 .OR. miter .EQ. 2) lwtem = lrw - 20 - lrtem
      lenwk = lwtem
      lyhn = lwm + lwtem
      lsavf = lyhn + lenyht
      lewt = lsavf + n
      lacor = lewt + n
      istatc = istate
      IF (istate .EQ. 1) go to 100
!-----------------------------------------------------------------------
! istate = 3.  move yh to its new location.
! note that only the part of yh needed for the next step, namely
! min(nq+1,maxord+2) columns, is actually moved.
! a temporary error weight array ewt is loaded if moss = 2.
! sparse matrix processing is done in iprep/prep if miter = 1 or 2.
! if maxord was reduced below nq, then the pointers are finally set
! so that savf is identical to yh(*,maxord+2).
!-----------------------------------------------------------------------
      lyhd = lyh - lyhn
      imax = lyhn - 1 + lenyhm
! move yh.  branch for move right, no move, or move left. --------------
      IF (lyhd) 70,80,74
 70   DO 72 i = lyhn,imax
        j = imax + lyhn - i
 72     rwork(j) = rwork(j+lyhd)
      go to 80
 74   DO 76 i = lyhn,imax
 76     rwork(i) = rwork(i+lyhd)
 80   lyh = lyhn
      iwork(22) = lyh
      IF (miter .EQ. 0 .OR. miter .EQ. 3) go to 92
      IF (moss .NE. 2) go to 85
! temporarily load ewt if miter = 1 or 2 and moss = 2. -----------------
      CALL ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      DO 82 i = 1,n
        IF (rwork(i+lewt-1) .LE. 0.0d0) go to 621
 82     rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
 85   CONTINUE
! iprep and prep do sparse matrix preprocessing if miter = 1 or 2. -----
      lsavf = min0(lsavf,lrw)
      lewt = min0(lewt,lrw)
      lacor = min0(lacor,lrw)
      CALL iprep (neq, y, rwork, iwork(lia), iwork(lja), ipflag, f, jac)
      lenrw = lwm - 1 + lenwk + lrest
      iwork(17) = lenrw
      IF (ipflag .NE. -1) iwork(23) = ipian
      IF (ipflag .NE. -1) iwork(24) = ipjan
      ipgo = -ipflag + 1
      go to (90, 628, 629, 630, 631, 632, 633), ipgo
 90   iwork(22) = lyh
      IF (lenrw .GT. lrw) go to 617
! set flag to signal parameter changes to stode. -----------------------
 92   jstart = -1
      IF (n .EQ. nyh) go to 200
! neq was reduced.  zero part of yh to avoid undefined references. -----
      i1 = lyh + l*nyh
      i2 = lyh + (maxord + 1)*nyh - 1
      IF (i1 .GT. i2) go to 200
      DO 95 i = i1,i2
 95     rwork(i) = 0.0d0
      go to 200
!-----------------------------------------------------------------------
! block c.
! the next block is for the initial call only (istate = 1).
! it contains all remaining initializations, the initial call to f,
! the sparse matrix preprocessing (miter = 1 or 2), and the
! calculation of the initial step size.
! the error weights in ewt are inverted after being loaded.
!-----------------------------------------------------------------------
 100  CONTINUE
      lyh = lyhn
      iwork(22) = lyh
      tn = t
      nst = 0
      h = 1.0d0
      nnz = 0
      ngp = 0
      nzl = 0
      nzu = 0
! load the initial value vector in yh. ---------------------------------
      DO 105 i = 1,n
 105    rwork(i+lyh-1) = y(i)
! initial call to f.  (lf0 points to yh(*,2).) -------------------------
      lf0 = lyh + nyh
      CALL f (neq, t, y, rwork(lf0))
      nfe = 1
! load and invert the ewt array.  (h is temporarily set to 1.0.) -------
      CALL ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      DO 110 i = 1,n
        IF (rwork(i+lewt-1) .LE. 0.0d0) go to 621
 110    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
      IF (miter .EQ. 0 .OR. miter .EQ. 3) go to 120
! iprep and prep do sparse matrix preprocessing if miter = 1 or 2. -----
      lacor = min0(lacor,lrw)
      CALL iprep (neq, y, rwork, iwork(lia), iwork(lja), ipflag, f, jac)
      lenrw = lwm - 1 + lenwk + lrest
      iwork(17) = lenrw
      IF (ipflag .NE. -1) iwork(23) = ipian
      IF (ipflag .NE. -1) iwork(24) = ipjan
      ipgo = -ipflag + 1
      go to (115, 628, 629, 630, 631, 632, 633), ipgo
 115  iwork(22) = lyh
      IF (lenrw .GT. lrw) go to 617
! check tcrit for legality (itask = 4 or 5). ---------------------------
 120  CONTINUE
      IF (itask .NE. 4 .AND. itask .NE. 5) go to 125
      tcrit = rwork(1)
      IF ((tcrit - tout)*(tout - t) .LT. 0.0d0) go to 625
      IF (h0 .NE. 0.0d0 .AND. (t + h0 - tcrit)*h0 .GT. 0.0d0) &
         h0 = tcrit - t
! initialize all remaining parameters. ---------------------------------
 125  uround = d1mach(4)
      jstart = 0
      IF (miter .NE. 0) rwork(lwm) = dsqrt(uround)
      msbj = 50
      nslj = 0
      ccmxj = 0.2d0
      psmall = 1000.0d0*uround
      rbig = 0.01d0/psmall
      nhnil = 0
      nje = 0
      nlu = 0
      nslast = 0
      hu = 0.0d0
      nqu = 0
      ccmax = 0.3d0
      maxcor = 3
      msbp = 20
      mxncf = 10
!-----------------------------------------------------------------------
! the coding below computes the step size, h0, to be attempted on the
! first step, unless the user has supplied a value for this.
! first check that tout - t differs significantly from zero.
! a scalar tolerance quantity tol is computed, as max(rtol(i))
! if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
! so as to be between 100*uround and 1.0e-3.
! then the computed value h0 is given by..
!                                      neq
!   h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )
!                                       1
! where   w0     = max ( abs(t), abs(tout) ),
!         f(i)   = i-th component of initial value of f,
!         ywt(i) = ewt(i)/tol  (a weight for y(i)).
! the sign of h0 is inferred from the initial values of tout and t.
!-----------------------------------------------------------------------
      lf0 = lyh + nyh
      IF (h0 .NE. 0.0d0) go to 180
      tdist = dabs(tout - t)
      w0 = dmax1(dabs(t),dabs(tout))
      IF (tdist .LT. 2.0d0*uround*w0) go to 622
      tol = rtol(1)
      IF (itol .LE. 2) go to 140
      DO 130 i = 1,n
 130    tol = dmax1(tol,rtol(i))
 140  IF (tol .GT. 0.0d0) go to 160
      atoli = atol(1)
      DO 150 i = 1,n
        IF (itol .EQ. 2 .OR. itol .EQ. 4) atoli = atol(i)
        ayi = dabs(y(i))
        IF (ayi .NE. 0.0d0) tol = dmax1(tol,atoli/ayi)
 150    CONTINUE
 160  tol = dmax1(tol,100.0d0*uround)
      tol = dmin1(tol,0.001d0)
      sum = vnorm (n, rwork(lf0), rwork(lewt))
      sum = 1.0d0/(tol*w0*w0) + tol*sum**2
      h0 = 1.0d0/dsqrt(sum)
      h0 = dmin1(h0,tdist)
      h0 = dsign(h0,tout-t)
! adjust h0 if necessary to meet hmax bound. ---------------------------
 180  rh = dabs(h0)*hmxi
      IF (rh .GT. 1.0d0) h0 = h0/rh
! load h with h0 and scale yh(*,2) by h0. ------------------------------
      h = h0
      DO 190 i = 1,n
 190    rwork(i+lf0-1) = h0*rwork(i+lf0-1)
      go to 270
!-----------------------------------------------------------------------
! block d.
! the next code block is for continuation calls only (istate = 2 or 3)
! and is to check stop conditions before taking a step.
!-----------------------------------------------------------------------
 200  nslast = nst
      go to (210, 250, 220, 230, 240), itask
 210  IF ((tn - tout)*h .LT. 0.0d0) go to 250
      CALL intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      IF (iflag .NE. 0) go to 627
      t = tout
      go to 420
 220  tp = tn - hu*(1.0d0 + 100.0d0*uround)
      IF ((tp - tout)*h .GT. 0.0d0) go to 623
      IF ((tn - tout)*h .LT. 0.0d0) go to 250
      go to 400
 230  tcrit = rwork(1)
      IF ((tn - tcrit)*h .GT. 0.0d0) go to 624
      IF ((tcrit - tout)*h .LT. 0.0d0) go to 625
      IF ((tn - tout)*h .LT. 0.0d0) go to 245
      CALL intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      IF (iflag .NE. 0) go to 627
      t = tout
      go to 420
 240  tcrit = rwork(1)
      IF ((tn - tcrit)*h .GT. 0.0d0) go to 624
 245  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .LE. 100.0d0*uround*hmx
      IF (ihit) go to 400
      tnext = tn + h*(1.0d0 + 4.0d0*uround)
      IF ((tnext - tcrit)*h .LE. 0.0d0) go to 250
      h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)
      IF (istate .EQ. 2) jstart = -2
!-----------------------------------------------------------------------
! block e.
! the next block is normally executed for all calls and contains
! the call to the one-step core integrator stode.
!
! this is a looping point for the integration steps.
!
! first check for too many steps being taken, update ewt (if not at
! start of problem), check for too much accuracy being requested, and
! check for h below the roundoff level in t.
!-----------------------------------------------------------------------
 250  CONTINUE
      IF ((nst-nslast) .GE. mxstep) go to 500
      CALL ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      DO 260 i = 1,n
        IF (rwork(i+lewt-1) .LE. 0.0d0) go to 510
 260    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
 270  tolsf = uround*vnorm (n, rwork(lyh), rwork(lewt))
      IF (tolsf .LE. 1.0d0) go to 280
      tolsf = tolsf*2.0d0
      IF (nst .EQ. 0) go to 626
      go to 520
 280  IF ((tn + h) .NE. tn) go to 290
      nhnil = nhnil + 1
      IF (nhnil .GT. mxhnil) go to 290
!      CALL xerrwv(50hlsodes-- warning..internal t (=r1) and h (=r2) are, &
!         50, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv( &
!        60h      such that in the machine, t + h = t on the next step  , &
!         60, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv(50h      (h = step size). solver will CONTINUE anyway, &
!         50, 101, 0, 0, 0, 0, 2, tn, h)
      IF (nhnil .LT. mxhnil) go to 290
!      CALL xerrwv(50hlsodes-- above warning has been issued i1 times.  , &
!         50, 102, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv(50h      it will not be issued again for this problem, &
!         50, 102, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
 290  CONTINUE
!-----------------------------------------------------------------------
!    call stode(neq,y,yh,nyh,yh,ewt,savf,acor,wm,wm,f,jac,prjs,slss)
!-----------------------------------------------------------------------
      ! type mismatch here HSJ rwork(lwm) ==> both int and real
      ! see      SUBROUTINE stode (neq, y, yh, nyh, yh1, ewt, savf, acor, &
      !           wm, iwm, f, jac, pjac, slvs)
      iwmint(1) = INT(rwork(lwm))
      !CALL stode (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt), &
      !   rwork(lsavf), rwork(lacor), rwork(lwm), rwork(lwm), &
      !   f, jac, prjs, slss)
      CALL stode (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),    &
                  rwork(lsavf), rwork(lacor), rwork(lwm), iwmint(1),      &
                  f, jac, prjs, slss)
      kgo = 1 - kflag
      go to (300, 530, 540, 550), kgo
!-----------------------------------------------------------------------
! block f.
! the following block handles the case of a successful return from the
! core integrator (kflag = 0).  test for stop conditions.
!-----------------------------------------------------------------------
 300  init = 1
      go to (310, 400, 330, 340, 350), itask
! itask = 1.  if tout has been reached, interpolate. -------------------
 310  IF ((tn - tout)*h .LT. 0.0d0) go to 250
      CALL intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
! itask = 3.  jump to exit if tout was reached. ------------------------
 330  IF ((tn - tout)*h .GE. 0.0d0) go to 400
      go to 250
! itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.
 340  IF ((tn - tout)*h .LT. 0.0d0) go to 345
      CALL intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
 345  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .LE. 100.0d0*uround*hmx
      IF (ihit) go to 400
      tnext = tn + h*(1.0d0 + 4.0d0*uround)
      IF ((tnext - tcrit)*h .LE. 0.0d0) go to 250
      h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)
      jstart = -2
      go to 250
! itask = 5.  see if tcrit was reached and jump to exit. ---------------
 350  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .LE. 100.0d0*uround*hmx
!-----------------------------------------------------------------------
! block g.
! the following block handles all successful returns from lsodes.
! if itask .ne. 1, y is loaded from yh and t is set accordingly.
! istate is set to 2, the illegal input counter is zeroed, and the
! optional outputs are loaded into the work arrays before returning.
! if istate = 1 and tout = t, there is a return with no action taken,
! except that if this has happened repeatedly, the run is terminated.
!-----------------------------------------------------------------------
 400  DO 410 i = 1,n
 410    y(i) = rwork(i+lyh-1)
      t = tn
      IF (itask .NE. 4 .AND. itask .NE. 5) go to 420
      IF (ihit) t = tcrit
 420  istate = 2
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      iwork(19) = nnz
      iwork(20) = ngp
      iwork(21) = nlu
      iwork(25) = nzl
      iwork(26) = nzu
      RETURN
!
 430  ntrep = ntrep + 1
      IF (ntrep .LT. 5) RETURN
 !     CALL xerrwv( &
!        60hlsodes-- repeated calls WITH istate = 1 and tout = t (=r1)  , &
!         60, 301, 0, 0, 0, 0, 1, t, 0.0d0)
      go to 800
!-----------------------------------------------------------------------
! block h.
! the following block handles all unsuccessful returns other than
! those for illegal input.  first the error message routine is called.
! if there was an error test or convergence test failure, imxer is set.
! then y is loaded from yh, t is set to tn, and the illegal input
! counter illin is set to 0.  the optional outputs are loaded into
! the work arrays before returning.
!-----------------------------------------------------------------------
! the maximum number of steps was taken before reaching tout. ----------
 500   CONTINUE
!      CALL xerrwv(50hlsodes-- at current t (=r1), mxstep (=i1) steps   , &
!         50, 201, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv(50h      taken on this CALL before reaching tout     , &
!         50, 201, 0, 1, mxstep, 0, 1, tn, 0.0d0)
      istate = -1
      go to 580
! ewt(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  ewti = rwork(lewt+i-1)
!      CALL xerrwv(50hlsodes-- at t (=r1), ewt(i1) has become r2 .LE. 0., &
!         50, 202, 0, 1, i, 0, 2, tn, ewti)
      istate = -6
      go to 580
! too much accuracy requested for machine precision. -------------------
 520 CONTINUE
!     CALL xerrwv(50hlsodes-- at t (=r1), too much accuracy requested  , &
!         50, 203, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv(50h      for PRECISION of machine..  see tolsf (=r2) , &
!         50, 203, 0, 0, 0, 0, 2, tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      go to 580
! kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----
 530  CONTINUE
!      CALL xerrwv(50hlsodes-- at t(=r1) and step size h(=r2), the error, &
!         50, 204, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv(50h      test failed repeatedly or WITH ABS(h) = hmin, &
!         50, 204, 0, 0, 0, 0, 2, tn, h)
      istate = -4
      go to 560
! kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----
 540  CONTINUE
!      CALL xerrwv(50hlsodes-- at t (=r1) and step size h (=r2), the    , &
!         50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv(50h      corrector convergence failed repeatedly     , &
!         50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv(30h      or WITH ABS(h) = hmin   , &
!         30, 205, 0, 0, 0, 0, 2, tn, h)
      istate = -5
      go to 560
! kflag = -3.  fatal error flag returned by prjs or slss (cdrv). -------
 550  CONTINUE
!      CALL xerrwv(50hlsodes-- at t (=r1) and step size h (=r2), a fatal, &
!         50, 207, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv(50h      error flag was returned by cdrv (by way of  , &
!         50, 207, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv(30h      SUBROUTINE prjs or slss), &
!         30, 207, 0, 0, 0, 0, 2, tn, h)
      istate = -7
      go to 580
! compute imxer if relevant. -------------------------------------------
 560  big = 0.0d0
      imxer = 1
      DO 570 i = 1,n
        size = dabs(rwork(i+lacor-1)*rwork(i+lewt-1))
        IF (big .GE. size) go to 570
        big = size
        imxer = i
 570    CONTINUE
      iwork(16) = imxer
! set y vector, t, illin, and optional outputs. ------------------------
 580  DO 590 i = 1,n
 590    y(i) = rwork(i+lyh-1)
      t = tn
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      iwork(19) = nnz
      iwork(20) = ngp
      iwork(21) = nlu
      iwork(25) = nzl
      iwork(26) = nzu
      RETURN
!-----------------------------------------------------------------------
! block i.
! the following block handles all error returns due to illegal input
! (istate = -3), as detected before calling the core integrator.
! first the error message routine is called.  then if there have been
! 5 consecutive such returns just before this call to the solver,
! the run is halted.
!-----------------------------------------------------------------------
 601  CONTINUE
!      CALL xerrwv(30hlsodes-- istate (=i1) illegal , &
!         30, 1, 0, 1, istate, 0, 0, 0.0d0, 0.0d0)
      go to 700
 602  CONTINUE
!      CALL xerrwv(30hlsodes-- itask (=i1) illegal  , &
!         30, 2, 0, 1, itask, 0, 0, 0.0d0, 0.0d0)
      go to 700
 603 CONTINUE
!       CALL xerrwv(50hlsodes-- istate .GT. 1 but lsodes not initialized , &
!         50, 3, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      go to 700
 604  CONTINUE
!     CALL xerrwv(30hlsodes-- neq (=i1) .LT. 1     , &
!         30, 4, 0, 1, neq(1), 0, 0, 0.0d0, 0.0d0)
      go to 700
 605  CONTINUE
!     CALL xerrwv(50hlsodes-- istate = 3 and neq increased (i1 to i2)  , &
!         50, 5, 0, 2, n, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 606  CONTINUE
!     CALL xerrwv(30hlsodes-- itol (=i1) illegal   , &
!         30, 6, 0, 1, itol, 0, 0, 0.0d0, 0.0d0)
      go to 700
 607  CONTINUE
!     CALL xerrwv(30hlsodes-- iopt (=i1) illegal   , &
!         30, 7, 0, 1, iopt, 0, 0, 0.0d0, 0.0d0)
      go to 700
 608  CONTINUE
!      CALL xerrwv(30hlsodes-- mf (=i1) illegal     , &
!         30, 8, 0, 1, mf, 0, 0, 0.0d0, 0.0d0)
      go to 700
 609  CONTINUE
!     CALL xerrwv(30hlsodes-- seth (=r1) .LT. 0.0  , &
!         30, 9, 0, 0, 0, 0, 1, seth, 0.0d0)
      go to 700
 611  CONTINUE
!     CALL xerrwv(30hlsodes-- maxord (=i1) .LT. 0  , &
!         30, 11, 0, 1, maxord, 0, 0, 0.0d0, 0.0d0)
      go to 700
 612  CONTINUE
!      CALL xerrwv(30hlsodes-- mxstep (=i1) .LT. 0  , &
!         30, 12, 0, 1, mxstep, 0, 0, 0.0d0, 0.0d0)
      go to 700
 613  CONTINUE
!     CALL xerrwv(30hlsodes-- mxhnil (=i1) .LT. 0  , &
!         30, 13, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
      go to 700
 614  CONTINUE
!     CALL xerrwv(40hlsodes-- tout (=r1) behind t (=r2)      , &
!         40, 14, 0, 0, 0, 0, 2, tout, t)
!      CALL xerrwv(50h      integration direction is given by h0 (=r1)  , &
!         50, 14, 0, 0, 0, 0, 1, h0, 0.0d0)
      go to 700
 615 CONTINUE
!       CALL xerrwv(30hlsodes-- hmax (=r1) .LT. 0.0  , &
!         30, 15, 0, 0, 0, 0, 1, hmax, 0.0d0)
      go to 700
 616  CONTINUE
!     CALL xerrwv(30hlsodes-- hmin (=r1) .LT. 0.0  , &
!         30, 16, 0, 0, 0, 0, 1, hmin, 0.0d0)
      go to 700
 617  CONTINUE
!     CALL xerrwv(50hlsodes-- rwork length is insufficient to proceed. , &
!         50, 17, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv( &
!        60h        length needed is .GE. lenrw (=i1), exceeds lrw (=i2), &
!         60, 17, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 618  CONTINUE
!      CALL xerrwv(50hlsodes-- iwork length is insufficient to proceed. , &
!         50, 18, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv( &
!        60h        length needed is .GE. leniw (=i1), exceeds liw (=i2), &
!         60, 18, 0, 2, leniw, liw, 0, 0.0d0, 0.0d0)
      go to 700
 619  CONTINUE
!     CALL xerrwv(40hlsodes-- rtol(i1) is r1 .LT. 0.0        , &
!         40, 19, 0, 1, i, 0, 1, rtoli, 0.0d0)
      go to 700
 620  CONTINUE
!      CALL xerrwv(40hlsodes-- atol(i1) is r1 .LT. 0.0        , &
!         40, 20, 0, 1, i, 0, 1, atoli, 0.0d0)
      go to 700
 621  ewti = rwork(lewt+i-1)
!      CALL xerrwv(40hlsodes-- ewt(i1) is r1 .LE. 0.0         , &
!         40, 21, 0, 1, i, 0, 1, ewti, 0.0d0)
      go to 700
 622  CONTINUE
!     CALL xerrwv( &
!        60hlsodes-- tout (=r1) too CLOSE to t(=r2) to start integration, &
!         60, 22, 0, 0, 0, 0, 2, tout, t)
      go to 700
 623  CONTINUE
!       CALL xerrwv( &
!       60hlsodes-- itask = i1 and tout (=r1) behind tcur - hu (= r2)  , &
!         60, 23, 0, 1, itask, 0, 2, tout, tp)
      go to 700
 624 CONTINUE
! CALL xerrwv( &
!        60hlsodes-- itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   , &
!         60, 24, 0, 0, 0, 0, 2, tcrit, tn)
      go to 700
 625  CONTINUE
! CALL xerrwv( &
!       60hlsodes-- itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   , &
!         60, 25, 0, 0, 0, 0, 2, tcrit, tout)
      go to 700
 626  CONTINUE
! CALL xerrwv(50hlsodes-- at start of problem, too much accuracy   , &
!        50, 26, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!    CALL xerrwv( &
!     60h      requested for PRECISION of machine..  see tolsf (=r1) , &
!         60, 26, 0, 0, 0, 0, 1, tolsf, 0.0d0)
      rwork(14) = tolsf
      go to 700
 627 CONTINUE
!  CALL xerrwv(50hlsodes-- trouble from intdy. itask = i1, tout = r1, &
!        50, 27, 0, 1, itask, 0, 1, tout, 0.0d0)
      go to 700
 628 CONTINUE
!  CALL xerrwv( &
!        60hlsodes-- rwork length insufficient (for SUBROUTINE prep).   , &
!         60, 28, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv( &
!        60h        length needed is .GE. lenrw (=i1), exceeds lrw (=i2), &
!         60, 28, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 629 CONTINUE
!  CALL xerrwv( &
!        60hlsodes-- rwork length insufficient (for SUBROUTINE jgroup). , &
!         60, 29, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv( &
!        60h        length needed is .GE. lenrw (=i1), exceeds lrw (=i2), &
!         60, 29, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 630 CONTINUE
!  CALL xerrwv( &
!       60hlsodes-- rwork length insufficient (for SUBROUTINE odrv).   , &
!         60, 30, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv( &
!        60h        length needed is .GE. lenrw (=i1), exceeds lrw (=i2), &
!         60, 30, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 631 CONTINUE
!  CALL xerrwv( &
!        60hlsodes-- error from odrv in yale sparse matrix package      , &
!         60, 31, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      imul = (iys - 1)/n
      irem = iys - imul*n
!      CALL xerrwv( &
!        60h      at t (=r1), odrv returned error flag = i1*neq + i2.   , &
!         60, 31, 0, 2, imul, irem, 1, tn, 0.0d0)
      go to 700
 632 CONTINUE
!  CALL xerrwv( &
!        60hlsodes-- rwork length insufficient (for SUBROUTINE cdrv).   , &
!         60, 32, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      CALL xerrwv( &
!        60h        length needed is .GE. lenrw (=i1), exceeds lrw (=i2), &
!         60, 32, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 633 CONTINUE
!  CALL xerrwv( &
!        60hlsodes-- error from cdrv in yale sparse matrix package      , &
!         60, 33, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      imul = (iys - 1)/n
      irem = iys - imul*n
!      CALL xerrwv( &
!        60h      at t (=r1), cdrv returned error flag = i1*neq + i2.   , &
!         60, 33, 0, 2, imul, irem, 1, tn, 0.0d0)
!      IF (imul .EQ. 2) CALL xerrwv( &
!        60h        duplicate ENTRY in sparsity structure descriptors   , &
!         60, 33, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!      IF (imul .EQ. 3 .OR. imul .EQ. 6) CALL xerrwv( &
!        60h        insufficient storage for nsfc (called by cdrv)      , &
!         60, 33, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!
 700  IF (illin .EQ. 5) go to 710
      illin = illin + 1
      istate = -3
      RETURN
 710 CONTINUE
! CALL xerrwv(50hlsodes-- repeated occurrences of illegal input    , &
!         50, 302, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
!
 800 CONTINUE
!  CALL xerrwv(50hlsodes-- run aborted.. apparent infinite loop     , &
!         50, 303, 2, 0, 0, 0, 0, 0.0d0, 0.0d0)
      RETURN
!----------------------- end of subroutine lsodes ----------------------
      END SUBROUTINE lsodes






      SUBROUTINE my_leqt1p (a, m, n, b, ib, idgt, d1, d2, ierr)
!*******************************************************************************
!*                                                                             *
!*    ROUTINE NAME:  LEQT1P                                                    *
!*                                                                             *
!*    PURPOSE:  solves linear equations returning the solution matrix in b     *
!*              and components of the determinant of matrix A in d1 and d2     *
!*                                                                             *
!*    AUTHOR:   Steve Reid - 10 November 1999                                  *
!*              S & E Computing                                                *
!*                                                                             *
!*    INPUTS:   a     - vector of length n(n + 1) / 2 containing the n * n     *
!*                      coefficient matrix of the equation AX = B.  Matrix A   *
!*                      is a positive definite symmetric matrix stored in      *
!*                      symmetric storage mode.                                *
!*              m     - number of right-hand sides (columns in B)              *
!*              n     - order of matrix a (number of rows in B)                *
!*              b     - n * m matrix containing the right-hand sides of the    *
!*                      equation AX = B                                        *
!*              ib    - row dimension of matrix b exactly as specified in the  *
!*                      calling program                                        *
!*              idgt  - decimal digit precision for elements of matrix a       *
!*                      (currently not used)                                   *
!*                                                                             *
!*    OUTPUTS:  a     - lower triangular matrix L, where A = L * L-transpose   *
!*                      L is stored in symmetric storage mode with the diagon- *
!*                      al elements of L in reciprocal form.                   *
!*              b     - contains the n * m solution matrix B                   *
!*              d1,d2 - components of the determinant of matrix A, such that   *
!*                      det(A) = d1 * 2. ** d2                                 *
!*              ierr  - error code (terminal error)                            *
!*                       129 = input matrix A is algorithmically not positive  *
!*                             definite                                        *
!*                                                                             *
!*    ROUTINES CALLED:  dppdi, dppfa, dppsl, uertst                            *
!*                                                                             *
!*    LIMITATIONS:  none                                                       *
!*                                                                             *
!*******************************************************************************
      IMPLICIT NONE
      INTEGER(I4B) ib, idgt, ierr, ijob, m, n
      REAL(DP) a(*), b(ib, *), d1, d2, det(2)
      CHARACTER*6 isub /'leqt1p'/
!
      ierr = 0
!
!--Decompose the matrix A.
!
      CALL dppfa (a, n, ierr)
      IF (ierr .NE. 0) THEN
         ierr = 129    ! matrix A is not positive definite
         GOTO 9000
      ENDIF
!
!--Compute the determinant components of matrix A.
!
      ijob = 10     ! only compute the determinant
      CALL dppdi (a, n, det, ijob)
      d1 = det(1)
      d2 = det(2) / dlog (2.d0)
      IF (ierr .EQ. 0) go to 9990
!
!--Solve the real symmetric positive definite system.
!
      CALL dppsl (a, n, b)
      go to 9990
!
!--Error exit.
!
 9000 CONTINUE
      CALL my_uertst (ierr, isub)
!
!--Normal exit.
!
 9990 CONTINUE
      RETURN
      END       SUBROUTINE my_leqt1p



!
      SUBROUTINE CUBGCV(X,F,DF,N,Y,C,IC,VAR,JOB,SE,WK,IER)
!     ALGORITHM 642 COLLECTED ALGORITHMS FROM ACM.
!     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.12, NO. 2,
!     JUN., 1986, P. 150.
!   SUBROUTINE NAME     - CUBGCV
!
!--------------------------------------------------------------------------
!
!   COMPUTER            - VAX/DOUBLE
!
!   AUTHOR              - M.F.HUTCHINSON
!                         CSIRO DIVISION OF MATHEMATICS AND STATISTICS
!                         P.O. BOX 1965
!                         CANBERRA, ACT 2601
!                         AUSTRALIA
!
!   LATEST REVISION     - 15 AUGUST 1985
!
!   PURPOSE             - CUBIC SPLINE DATA SMOOTHER
!
!   USAGE               - CALL CUBGCV (X,F,DF,N,Y,C,IC,VAR,JOB,SE,WK,IER)
!
!   ARGUMENTS    X      - VECTOR OF LENGTH N CONTAINING THE
!                           ABSCISSAE OF THE N DATA POINTS
!                           (X(I),F(I)) I=1..N. (INPUT) X
!                           MUST BE ORDERED SO THAT
!                           X(I) .LT. X(I+1).
!                F      - VECTOR OF LENGTH N CONTAINING THE
!                           ORDINATES (OR FUNCTION VALUES)
!                           OF THE N DATA POINTS (INPUT).
!                DF     - VECTOR OF LENGTH N. (INPUT/OUTPUT)
!                           DF(I) IS THE RELATIVE STANDARD DEVIATION
!                           OF THE ERROR ASSOCIATED WITH DATA POINT I.
!                           EACH DF(I) MUST BE POSITIVE.  THE VALUES IN
!                           DF ARE SCALED BY THE SUBROUTINE SO THAT
!                           THEIR MEAN SQUARE VALUE IS 1, AND UNSCALED
!                           AGAIN ON NORMAL EXIT.
!                           THE MEAN SQUARE VALUE OF THE DF(I) IS RETURNED
!                           IN WK(7) ON NORMAL EXIT.
!                           IF THE ABSOLUTE STANDARD DEVIATIONS ARE KNOWN,
!                           THESE SHOULD BE PROVIDED IN DF AND THE ERROR
!                           VARIANCE PARAMETER VAR (SEE BELOW) SHOULD THEN
!                           BE SET TO 1.
!                           IF THE RELATIVE STANDARD DEVIATIONS ARE UNKNOWN,
!                           SET EACH DF(I)=1.
!                N      - NUMBER OF DATA POINTS (INPUT).
!                           N MUST BE .GE. 3.
!                Y,C    - SPLINE COEFFICIENTS. (OUTPUT) Y
!                           IS A VECTOR OF LENGTH N. C IS
!                           AN N-1 BY 3 MATRIX. THE VALUE
!                           OF THE SPLINE APPROXIMATION AT T IS
!                           S(T)=((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
!                           WHERE X(I).LE.T.LT.X(I+1) AND
!                           D = T-X(I).
!                IC     - ROW DIMENSION OF MATRIX C EXACTLY
!                           AS SPECIFIED IN THE DIMENSION
!                           STATEMENT IN THE CALLING PROGRAM. (INPUT)
!                VAR    - ERROR VARIANCE. (INPUT/OUTPUT)
!                           IF VAR IS NEGATIVE (I.E. UNKNOWN) THEN
!                           THE SMOOTHING PARAMETER IS DETERMINED
!                           BY MINIMIZING THE GENERALIZED CROSS VALIDATION
!                           AND AN ESTIMATE OF THE ERROR VARIANCE IS
!                           RETURNED IN VAR.
!                           IF VAR IS NON-NEGATIVE (I.E. KNOWN) THEN THE
!                           SMOOTHING PARAMETER IS DETERMINED TO MINIMIZE
!                           AN ESTIMATE, WHICH DEPENDS ON VAR, OF THE TRUE
!                           MEAN SQUARE ERROR, AND VAR IS UNCHANGED.
!                           IN PARTICULAR, IF VAR IS ZERO, THEN AN
!                           INTERPOLATING NATURAL CUBIC SPLINE IS CALCULATED.
!                           VAR SHOULD BE SET TO 1 IF ABSOLUTE STANDARD
!                           DEVIATIONS HAVE BEEN PROVIDED IN DF (SEE ABOVE).
!                JOB    - JOB SELECTION PARAMETER. (INPUT)
!                         JOB = 0 SHOULD BE SELECTED IF POINT STANDARD ERROR
!                           ESTIMATES ARE NOT REQUIRED IN SE.
!                         JOB = 1 SHOULD BE SELECTED IF POINT STANDARD ERROR
!                           ESTIMATES ARE REQUIRED IN SE.
!                SE     - VECTOR OF LENGTH N CONTAINING BAYESIAN STANDARD
!                           ERROR ESTIMATES OF THE FITTED SPLINE VALUES IN Y.
!                           SE IS NOT REFERENCED IF JOB=0. (OUTPUT)
!                WK     - WORK VECTOR OF LENGTH 7*(N + 2). ON NORMAL EXIT THE
!                           FIRST 7 VALUES OF WK ARE ASSIGNED AS FOLLOWS:-
!
!                           WK(1) = SMOOTHING PARAMETER (= RHO/(RHO + 1))
!                           WK(2) = ESTIMATE OF THE NUMBER OF DEGREES OF
!                                   FREEDOM OF THE RESIDUAL SUM OF SQUARES
!                           WK(3) = GENERALIZED CROSS VALIDATION
!                           WK(4) = MEAN SQUARE RESIDUAL
!                           WK(5) = ESTIMATE OF THE TRUE MEAN SQUARE ERROR
!                                   AT THE DATA POINTS
!                           WK(6) = ESTIMATE OF THE ERROR VARIANCE
!                           WK(7) = MEAN SQUARE VALUE OF THE DF(I)
!
!                           IF WK(1)=0 (RHO=0) AN INTERPOLATING NATURAL CUBIC
!                           SPLINE HAS BEEN CALCULATED.
!                           IF WK(1)=1 (RHO=INFINITE) A LEAST SQUARES
!                           REGRESSION LINE HAS BEEN CALCULATED.
!                           WK(2) IS AN ESTIMATE OF THE NUMBER OF DEGREES OF
!                           FREEDOM OF THE RESIDUAL WHICH REDUCES TO THE
!                           USUAL VALUE OF N-2 WHEN A LEAST SQUARES REGRESSION
!                           LINE IS CALCULATED.
!                           WK(3),WK(4),WK(5) ARE CALCULATED WITH THE DF(I)
!                           SCALED TO HAVE MEAN SQUARE VALUE 1.  THE
!                           UNSCALED VALUES OF WK(3),WK(4),WK(5) MAY BE
!                           CALCULATED BY DIVIDING BY WK(7).
!                           WK(6) COINCIDES WITH THE OUTPUT VALUE OF VAR IF
!                           VAR IS NEGATIVE ON INPUT.  IT IS CALCULATED WITH
!                           THE UNSCALED VALUES OF THE DF(I) TO FACILITATE
!                           COMPARISONS WITH A PRIORI VARIANCE ESTIMATES.
!
!                IER    - ERROR PARAMETER. (OUTPUT)
!                         TERMINAL ERROR
!                           IER = 129, IC IS LESS THAN N-1.
!                           IER = 130, N IS LESS THAN 3.
!                           IER = 131, INPUT ABSCISSAE ARE NOT
!                             ORDERED SO THAT X(I).LT.X(I+1).
!                           IER = 132, DF(I) IS NOT POSITIVE FOR SOME I.
!                           IER = 133, JOB IS NOT 0 OR 1.
!
!   PRECISION/HARDWARE  - DOUBLE
!
!   REQUIRED ROUTINES   - SPINT1,SPFIT1,SPCOF1,SPERR1
!
!   REMARKS      THE NUMBER OF ARITHMETIC OPERATIONS REQUIRED BY THE
!                SUBROUTINE IS PROPORTIONAL TO N.  THE SUBROUTINE
!                USES AN ALGORITHM DEVELOPED BY M.F. HUTCHINSON AND
!                F.R. DE HOOG, 'SMOOTHING NOISY DATA WITH SPLINE
!                FUNCTIONS', NUMER. MATH. (IN PRESS)
!
!--------------------------------------------------------------
!
!---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER N,IC,JOB,IER
      REAL(DP) X(N),F(N),DF(N),Y(N),C(IC,3),SE(N),VAR, &
                       WK(0:N+1,7)
!
!---SPECIFICATIONS FOR LOCAL VARIABLES---
      REAL(DP) DELTA,ERR,GF1,GF2,GF3,GF4,R1,R2,R3,R4,TAU,RATIO, &
                       AVH,AVDF,AVAR,ZERO,ONE,STAT(6),P,Q
!
      DATA RATIO/2.0D0/
      DATA TAU/1.618033989D0/
      DATA ZERO,ONE/0.0D0,1.0D0/
!
!---INITIALIZE---
      IER = 133
      IF (JOB.LT.0 .OR. JOB.GT.1) GO TO 140
      CALL SPINT1(X,AVH,F,DF,AVDF,N,Y,C,IC,WK,WK(0,4),IER)
      IF (IER.NE.0) GO TO 140
      AVAR = VAR
      IF (VAR.GT.ZERO) AVAR = VAR*AVDF*AVDF
!
!---CHECK FOR ZERO VARIANCE---
      IF (VAR.NE.ZERO) GO TO 10
      R1 = ZERO
      GO TO 90
!
!---FIND LOCAL MINIMUM OF GCV OR THE EXPECTED MEAN SQUARE ERROR---
   10 R1 = ONE
      R2 = RATIO*R1
      CALL SPFIT1(X,AVH,DF,N,R2,P,Q,GF2,AVAR,STAT,Y,C,IC,WK,WK(0,4), &
                  WK(0,6),WK(0,7))
   20 CALL SPFIT1(X,AVH,DF,N,R1,P,Q,GF1,AVAR,STAT,Y,C,IC,WK,WK(0,4), &
                  WK(0,6),WK(0,7))
      IF (GF1.GT.GF2) GO TO 30
!
!---EXIT IF P ZERO---
      IF (P.LE.ZERO) GO TO 100
      R2 = R1
      GF2 = GF1
      R1 = R1/RATIO
      GO TO 20

   30 R3 = RATIO*R2
   40 CALL SPFIT1(X,AVH,DF,N,R3,P,Q,GF3,AVAR,STAT,Y,C,IC,WK,WK(0,4), &
                  WK(0,6),WK(0,7))
      IF (GF3.GT.GF2) GO TO 50
!
!---EXIT IF Q ZERO---
      IF (Q.LE.ZERO) GO TO 100
      R2 = R3
      GF2 = GF3
      R3 = RATIO*R3
      GO TO 40

   50 R2 = R3
      GF2 = GF3
      DELTA = (R2-R1)/TAU
      R4 = R1 + DELTA
      R3 = R2 - DELTA
      CALL SPFIT1(X,AVH,DF,N,R3,P,Q,GF3,AVAR,STAT,Y,C,IC,WK,WK(0,4), &
                  WK(0,6),WK(0,7))
      CALL SPFIT1(X,AVH,DF,N,R4,P,Q,GF4,AVAR,STAT,Y,C,IC,WK,WK(0,4), &
                  WK(0,6),WK(0,7))
!
!---GOLDEN SECTION SEARCH FOR LOCAL MINIMUM---
   60 IF (GF3.GT.GF4) GO TO 70
      R2 = R4
      GF2 = GF4
      R4 = R3
      GF4 = GF3
      DELTA = DELTA/TAU
      R3 = R2 - DELTA
      CALL SPFIT1(X,AVH,DF,N,R3,P,Q,GF3,AVAR,STAT,Y,C,IC,WK,WK(0,4), &
                  WK(0,6),WK(0,7))
      GO TO 80

   70 R1 = R3
      GF1 = GF3
      R3 = R4
      GF3 = GF4
      DELTA = DELTA/TAU
      R4 = R1 + DELTA
      CALL SPFIT1(X,AVH,DF,N,R4,P,Q,GF4,AVAR,STAT,Y,C,IC,WK,WK(0,4), &
                  WK(0,6),WK(0,7))
   80 ERR = (R2-R1)/ (R1+R2)
      IF (ERR*ERR+ONE.GT.ONE .AND. ERR.GT.1.0D-6) GO TO 60
      R1 = (R1+R2)*0.5D0
!
!---CALCULATE SPLINE COEFFICIENTS---
   90 CALL SPFIT1(X,AVH,DF,N,R1,P,Q,GF1,AVAR,STAT,Y,C,IC,WK,WK(0,4), &
                  WK(0,6),WK(0,7))
  100 CALL SPCOF1(X,AVH,F,DF,N,P,Q,Y,C,IC,WK(0,6),WK(0,7))
!
!---OPTIONALLY CALCULATE STANDARD ERROR ESTIMATES---
      IF (VAR.GE.ZERO) GO TO 110
      AVAR = STAT(6)
      VAR = AVAR/ (AVDF*AVDF)
  110 IF (JOB.EQ.1) CALL SPERR1(X,AVH,DF,N,WK,P,AVAR,SE)
!
!---UNSCALE DF---
      DO 120 I = 1,N
         DF(I) = DF(I)*AVDF
  120 CONTINUE
!
!--PUT STATISTICS IN WK---
      DO 130 I = 0,5
         WK(I,1) = STAT(I+1)
  130 CONTINUE
      WK(5,1) = STAT(6)/ (AVDF*AVDF)
      WK(6,1) = AVDF*AVDF
      GO TO 150
!
!---CHECK FOR ERROR CONDITION---
  140 CONTINUE
!     IF (IER.NE.0) CONTINUE
  150 RETURN
      END       SUBROUTINE CUBGCV




      SUBROUTINE SPCOF1(X,AVH,Y,DY,N,P,Q,A,C,IC,U,V)
!
! CALCULATES COEFFICIENTS OF A CUBIC SMOOTHING SPLINE FROM
! PARAMETERS CALCULATED BY SUBROUTINE SPFIT1.
!
!---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER IC,N
      REAL(DP) X(N),Y(N),DY(N),P,Q,A(N),C(IC,3),U(0:N+1), &
                       V(0:N+1),AVH
!
!---SPECIFICATIONS FOR LOCAL VARIABLES---
      INTEGER I
      REAL(DP) H,QH
!
!---CALCULATE A---
      QH = Q/ (AVH*AVH)
      DO 10 I = 1,N
         A(I) = Y(I) - P*DY(I)*V(I)
         U(I) = QH*U(I)
   10 CONTINUE
!
!---CALCULATE C---
      DO 20 I = 1,N - 1
         H = X(I+1) - X(I)
         C(I,3) = (U(I+1)-U(I))/ (3.0D0*H)
         C(I,1) = (A(I+1)-A(I))/H - (H*C(I,3)+U(I))*H
         C(I,2) = U(I)
   20 CONTINUE
      RETURN
      END       SUBROUTINE SPCOF1



      SUBROUTINE SPERR1(X,AVH,DY,N,R,P,VAR,SE)
!
! CALCULATES BAYESIAN ESTIMATES OF THE STANDARD ERRORS OF THE FITTED
! VALUES OF A CUBIC SMOOTHING SPLINE BY CALCULATING THE DIAGONAL ELEMENTS
! OF THE INFLUENCE MATRIX.
!
!---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER N
      REAL(DP) X(N),DY(N),R(0:N+1,3),SE(N),AVH,P,VAR
!
!---SPECIFICATIONS FOR LOCAL VARIABLES---
      INTEGER I
      REAL(DP) F,G,H,F1,G1,H1,ZERO,ONE
      DATA ZERO,ONE/0.0D0,1.0D0/
!
!---INITIALIZE---
      H = AVH/ (X(2)-X(1))
      SE(1) = ONE - P*DY(1)*DY(1)*H*H*R(2,1)
      R(1,1) = ZERO
      R(1,2) = ZERO
      R(1,3) = ZERO
!
!---CALCULATE DIAGONAL ELEMENTS---
      DO 10 I = 2,N - 1
         F = H
         H = AVH/ (X(I+1)-X(I))
         G = -F - H
         F1 = F*R(I-1,1) + G*R(I-1,2) + H*R(I-1,3)
         G1 = F*R(I-1,2) + G*R(I,1) + H*R(I,2)
         H1 = F*R(I-1,3) + G*R(I,2) + H*R(I+1,1)
         SE(I) = ONE - P*DY(I)*DY(I)* (F*F1+G*G1+H*H1)
   10 CONTINUE
      SE(N) = ONE - P*DY(N)*DY(N)*H*H*R(N-1,1)
!
!---CALCULATE STANDARD ERROR ESTIMATES---
      DO 20 I = 1,N
         SE(I) = DSQRT(DMAX1(SE(I)*VAR,ZERO))*DY(I)
   20 CONTINUE
      RETURN
      END       SUBROUTINE SPERR1


      SUBROUTINE SPFIT1(X,AVH,DY,N,RHO,P,Q,FUN,VAR,STAT,A,C,IC,R,T,U,V)
!
! FITS A CUBIC SMOOTHING SPLINE TO DATA WITH RELATIVE
! WEIGHTING DY FOR A GIVEN VALUE OF THE SMOOTHING PARAMETER
! RHO USING AN ALGORITHM BASED ON THAT OF C.H. REINSCH (1967),
! NUMER. MATH. 10, 177-183.
!
! THE TRACE OF THE INFLUENCE MATRIX IS CALCULATED USING AN
! ALGORITHM DEVELOPED BY M.F.HUTCHINSON AND F.R.DE HOOG (NUMER.
! MATH., IN PRESS), ENABLING THE GENERALIZED CROSS VALIDATION
! AND RELATED STATISTICS TO BE CALCULATED IN ORDER N OPERATIONS.
!
! THE ARRAYS A, C, R AND T ARE ASSUMED TO HAVE BEEN INITIALIZED
! BY THE SUBROUTINE SPINT1.  OVERFLOW AND UNDERFLOW PROBLEMS ARE
! AVOIDED BY USING P=RHO/(1 + RHO) AND Q=1/(1 + RHO) INSTEAD OF
! RHO AND BY SCALING THE DIFFERENCES X(I+1) - X(I) BY AVH.
!
! THE VALUES IN DF ARE ASSUMED TO HAVE BEEN SCALED SO THAT THE
! SUM OF THEIR SQUARED VALUES IS N.  THE VALUE IN VAR, WHEN IT IS
! NON-NEGATIVE, IS ASSUMED TO HAVE BEEN SCALED TO COMPENSATE FOR
! THE SCALING OF THE VALUES IN DF.
!
! THE VALUE RETURNED IN FUN IS AN ESTIMATE OF THE TRUE MEAN SQUARE
! WHEN VAR IS NON-NEGATIVE, AND IS THE GENERALIZED CROSS VALIDATION
! WHEN VAR IS NEGATIVE.
!
!---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER IC,N
      REAL(DP) X(N),DY(N),RHO,STAT(6),A(N),C(IC,3),R(0:N+1,3), &
                       T(0:N+1,2),U(0:N+1),V(0:N+1),FUN,VAR,AVH,P,Q
!
!---LOCAL VARIABLES---
      INTEGER I
      REAL(DP) E,F,G,H,ZERO,ONE,TWO,RHO1
      DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/
!
!---USE P AND Q INSTEAD OF RHO TO PREVENT OVERFLOW OR UNDERFLOW---
      RHO1 = ONE + RHO
      P = RHO/RHO1
      Q = ONE/RHO1
      IF (RHO1.EQ.ONE) P = ZERO
      IF (RHO1.EQ.RHO) Q = ZERO
!
!---RATIONAL CHOLESKY DECOMPOSITION OF P*C + Q*T---
      F = ZERO
      G = ZERO
      H = ZERO
      DO 10 I = 0,1
         R(I,1) = ZERO
   10 CONTINUE
      DO 20 I = 2,N - 1
         R(I-2,3) = G*R(I-2,1)
         R(I-1,2) = F*R(I-1,1)
         R(I,1) = ONE/ (P*C(I,1)+Q*T(I,1)-F*R(I-1,2)-G*R(I-2,3))
         F = P*C(I,2) + Q*T(I,2) - H*R(I-1,2)
         G = H
         H = P*C(I,3)
   20 CONTINUE
!
!---SOLVE FOR U---
      U(0) = ZERO
      U(1) = ZERO
      DO 30 I = 2,N - 1
         U(I) = A(I) - R(I-1,2)*U(I-1) - R(I-2,3)*U(I-2)
   30 CONTINUE
      U(N) = ZERO
      U(N+1) = ZERO
      DO 40 I = N - 1,2,-1
         U(I) = R(I,1)*U(I) - R(I,2)*U(I+1) - R(I,3)*U(I+2)
   40 CONTINUE
!
!---CALCULATE RESIDUAL VECTOR V---
      E = ZERO
      H = ZERO
      DO 50 I = 1,N - 1
         G = H
         H = (U(I+1)-U(I))/ ((X(I+1)-X(I))/AVH)
         V(I) = DY(I)* (H-G)
         E = E + V(I)*V(I)
   50 CONTINUE
      V(N) = DY(N)* (-H)
      E = E + V(N)*V(N)
!
!---CALCULATE UPPER THREE BANDS OF INVERSE MATRIX---
      R(N,1) = ZERO
      R(N,2) = ZERO
      R(N+1,1) = ZERO
      DO 60 I = N - 1,2,-1
         G = R(I,2)
         H = R(I,3)
         R(I,2) = -G*R(I+1,1) - H*R(I+1,2)
         R(I,3) = -G*R(I+1,2) - H*R(I+2,1)
         R(I,1) = R(I,1) - G*R(I,2) - H*R(I,3)
   60 CONTINUE
!
!---CALCULATE TRACE---
      F = ZERO
      G = ZERO
      H = ZERO
      DO 70 I = 2,N - 1
         F = F + R(I,1)*C(I,1)
         G = G + R(I,2)*C(I,2)
         H = H + R(I,3)*C(I,3)
   70 CONTINUE
      F = F + TWO* (G+H)
!
!---CALCULATE STATISTICS---
      STAT(1) = P
      STAT(2) = F*P
      STAT(3) = N*E/ (F*F)
      STAT(4) = E*P*P/N
      STAT(6) = E*P/F
      IF (VAR.GE.ZERO) GO TO 80
      STAT(5) = STAT(6) - STAT(4)
      FUN = STAT(3)
      GO TO 90

   80 STAT(5) = DMAX1(STAT(4)-TWO*VAR*STAT(2)/N+VAR,ZERO)
      FUN = STAT(5)
   90 RETURN
      END       SUBROUTINE SPFIT1


      SUBROUTINE SPINT1(X,AVH,Y,DY,AVDY,N,A,C,IC,R,T,IER)
!
! INITIALIZES THE ARRAYS C, R AND T FOR ONE DIMENSIONAL CUBIC
! SMOOTHING SPLINE FITTING BY SUBROUTINE SPFIT1.  THE VALUES
! DF(I) ARE SCALED SO THAT THE SUM OF THEIR SQUARES IS N
! AND THE AVERAGE OF THE DIFFERENCES X(I+1) - X(I) IS CALCULATED
! IN AVH IN ORDER TO AVOID UNDERFLOW AND OVERFLOW PROBLEMS IN
! SPFIT1.
!
! SUBROUTINE SETS IER IF ELEMENTS OF X ARE NON-INCREASING,
! IF N IS LESS THAN 3, IF IC IS LESS THAN N-1 OR IF DY(I) IS
! NOT POSITIVE FOR SOME I.
!
!---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER N,IC,IER
      REAL(DP) X(N),Y(N),DY(N),A(N),C(IC,3),R(0:N+1,3), &
                       T(0:N+1,2),AVH,AVDY
!
!---SPECIFICATIONS FOR LOCAL VARIABLES---
      INTEGER I
      REAL(DP) E,F,G,H,ZERO
      DATA ZERO/0.0D0/
!
!---INITIALIZATION AND INPUT CHECKING---
      IER = 0
      IF (N.LT.3) GO TO 60
      IF (IC.LT.N-1) GO TO 70
!
!---GET AVERAGE X SPACING IN AVH---
      G = ZERO
      DO 10 I = 1,N - 1
         H = X(I+1) - X(I)
         IF (H.LE.ZERO) GO TO 80
         G = G + H
   10 CONTINUE
      AVH = G/ (N-1)
!
!---SCALE RELATIVE WEIGHTS---
      G = ZERO
      DO 20 I = 1,N
         IF (DY(I).LE.ZERO) GO TO 90
         G = G + DY(I)*DY(I)
   20 CONTINUE
      AVDY = DSQRT(G/N)
!
      DO 30 I = 1,N
         DY(I) = DY(I)/AVDY
   30 CONTINUE
!
!---INITIALIZE H,F---
      H = (X(2)-X(1))/AVH
      F = (Y(2)-Y(1))/H
!
!---CALCULATE A,T,R---
      DO 40 I = 2,N - 1
         G = H
         H = (X(I+1)-X(I))/AVH
         E = F
         F = (Y(I+1)-Y(I))/H
         A(I) = F - E
         T(I,1) = 2.0D0* (G+H)/3.0D0
         T(I,2) = H/3.0D0
         R(I,3) = DY(I-1)/G
         R(I,1) = DY(I+1)/H
         R(I,2) = -DY(I)/G - DY(I)/H
   40 CONTINUE
!
!---CALCULATE C = R'*R---
      R(N,2) = ZERO
      R(N,3) = ZERO
      R(N+1,3) = ZERO
      DO 50 I = 2,N - 1
         C(I,1) = R(I,1)*R(I,1) + R(I,2)*R(I,2) + R(I,3)*R(I,3)
         C(I,2) = R(I,1)*R(I+1,2) + R(I,2)*R(I+1,3)
         C(I,3) = R(I,1)*R(I+2,3)
   50 CONTINUE
      RETURN
!
!---ERROR CONDITIONS---
   60 IER = 130
      RETURN

   70 IER = 129
      RETURN

   80 IER = 131
      RETURN

   90 IER = 132
      RETURN
      END       SUBROUTINE SPINT1




      SUBROUTINE DCSSMO(H, N, TNODE, G, WGS, RHO, GSMO, B, C, D)
      IMPLICIT REAL(DP) (a-h, o-z)  ! added by S. Reid for onetwo 2/17/00
!
!  THIS SUBROUTINE COMPUTES THE DISCRETE NATURAL CUBIC
!  SPLINE DEFINED ON THE INTERVAL (TNODE(1),TNODE(N)) WHICH
!  SMOOTHS THROUGH THE DATA (TNODE(I),G(I)),I=1,2,...,N.
!  N MUST BE 2 OR GREATER. THE NODES MUST SATISFY TNODE(I)
!  .LT.TNODE(I+1). THE SOLUTION S(T) FOR T IN THE INTERVAL
!  (TNODE(I),TNODE(I+1)) IS GIVEN BY
!
!     S(T)=GSMO(I)+B(I)*(T-TNODE(I))+
!             C(I)*(T-TNODE(I))**2+D(I)*(T-TNODE(I))**3
!
      DIMENSION TNODE(N), G(N), WGS(N), GSMO(N), B(N), C(N), D(N)
!
!  INPUT  PARAMETERS(NONE OF THE INPUT PARAMETERS ARE CHANGED
!         BY THIS SUBROUTINE)
!
!  H     - THE STEP SIZE USED FOR THE DISCRETE CUBIC SPLINE
!  N     - NUMBER OF NODES (TNODE) AND DATA VALUES(G)
!  TNODE - REAL ARRAY CONTAINING THE NODES (TNODE(I).LT.
!          TNODE(I+1)).
!  G     - REAL ARRAY CONTAINING THE DATA VALUES.
!  WGS   - REAL ARRAY CONTAINING THE WEIGHTS WGS(I)
!          CORRESPONDING TO THE DATA (TNODE(I),G(I)).
!  RHO   - SIMPLE REAL VARIABLE CONTAINING THE POSITIVE
!          PARAMETER FOR VARYING THE SMOOTHNESS OF THE FIT.
!          IF RHO IS SMALL SMOOTHNESS IS EMPHASIZED.
!          IF RHO IS LARGE DATA FITTING IS EMPHASIZED.
!
!  OUTPUT PARAMETERS
!
!  GSMO  - REAL ARRAY CONTAINING THE SMOOTHED VALUES OF
!          THE DATA G(I),I=1,2,....,N.
!  B     - REAL ARRAY CONTAINING THE COEFFICIENTS B(I) FOR
!          THE TERMS (T-TNODE(I)).
!  C     - REAL ARRAY CONTAINING THE COEFFICIENTS C(I) FOR
!          THE TERMS (T-TNODE(I))**2.
!  D     - REAL ARRAY CONTAINING THE COEFFICIENTS D(I) FOR
!          THE TERMS (T-TNODE(I))**3.
!
      IF (N.EQ.2) GO TO 180
      N1 = N - 1
      N2 = N1 - 1
      N3 = N2 - 1
!  THE RIGHT HAND SIDE OF THE LINEAR SYSTEM FOR THE
!  C(I)'S WILL NOW BE CONSTRUCTED.
      DO 10 I=1,N
        C(I) = G(I)
   10 CONTINUE
      DO 20 I=1,N1
        C(I) = (C(I+1)-C(I))/(TNODE(I+1)-TNODE(I))
   20 CONTINUE
      DO 30 I=1,N2
        C(I) = 3.0*(C(I+1)-C(I))
   30 CONTINUE
!  THE RIGHT HAND SIDE IS NOW IN ARRAY C.
!
!  THE P.D. 5 BANDED SYMMETRIC MATRIX WILL NOW BE CONSTRUCTED.
!  THE THREE NEEDED DIAGONALS WILL BE STORED IN ARRAYS
!  GSMO,B,D.
      H2 = H*H
      H3 = H2*H
      R6 = 6.0*H3/RHO
      HI3 = TNODE(2) - TNODE(1)
      HI4 = TNODE(3) - TNODE(2)
      ETA3 = HI3 + HI3 + H2/HI3
      BETA3 = R6/(WGS(1)*HI3)
      BETA4 = R6/(WGS(2)*HI3*HI4)
      EPS3 = (BETA3+BETA4*HI4)/HI3
      H2DHI = H2/HI4
      ETA4 = HI4 + HI4 + H2DHI
      IF (N.EQ.3) GO TO 60
      HI5 = TNODE(4) - TNODE(3)
      BETA5 = R6/(WGS(3)*HI4*HI5)
      EPS4 = (BETA4*HI3+BETA5*HI5)/HI4
      GSMO(1) = ETA3 + ETA4 + BETA4 + BETA4 + EPS3 + EPS4
      P = H2DHI + BETA4 + BETA5 + EPS4
      B(1) = HI4 - P
      IF (N.EQ.4) GO TO 50
      DO 40 I=2,N3
        HI3 = HI4
        HI4 = HI5
        HI5 = TNODE(I+3) - TNODE(I+2)
        ETA3 = ETA4
        H2DHI = H2/HI4
        ETA4 = HI4 + HI4 + H2DHI
        BETA3 = BETA4
        BETA4 = BETA5
        BETA5 = R6/(WGS(I+2)*HI4*HI5)
        EPS3 = EPS4
        EPS4 = (BETA4*HI3+BETA5*HI5)/HI4
        D(I-1) = BETA4
        P = H2DHI + BETA4 + BETA5 + EPS4
        B(I) = HI4 - P
        GSMO(I) = ETA3 + ETA4 + BETA4 + BETA4 + EPS3 + EPS4
   40 CONTINUE
   50 HI3 = HI4
      HI4 = HI5
      ETA3 = ETA4
      ETA4 = HI4 + HI4 + H2/HI4
      BETA4 = BETA5
      EPS3 = EPS4
   60 BETA5 = R6/(WGS(N)*HI4)
      EPS4 = (BETA4*HI3+BETA5)/HI4
      GSMO(N2) = ETA3 + ETA4 + BETA4 + BETA4 + EPS3 + EPS4
!  THE P.D. 5 BANDED SYMMETRIC MATRIX IS COMPLETE.
!  THE SYSTEM OF LINEAR EQUATION WILL NOW BE SOLVED FOR THE
!  C(I)'S.
      IF (N.GT.3) GO TO 70
      C(1) = C(1)/GSMO(1)
      GO TO 150
   70 IF (N.GT.4) GO TO 80
      C(1) = (C(1)*GSMO(2)-C(2)*B(1))/(GSMO(1)*GSMO(2)-B(1)**2)
      C(2) = (C(2)-C(1)*B(1))/GSMO(2)
      GO TO 150
!  THIS SOLVE THE 5 BANDED SYSTEM WHEN K=N-2.GT.3.
   80 K = N2
      K1 = K - 1
      K2 = K1 - 1
      K3 = K2 - 1
!  THE 5 BANDED MATRIX WILL NOW BE FACTORED.
      B(1) = B(1)/GSMO(1)
      D(1) = D(1)/GSMO(1)
      P = GSMO(1)*B(1)
      GSMO(2) = GSMO(2) - P*B(1)
      B(2) = (B(2)-P*D(1))/GSMO(2)
      IF (K.EQ.3) GO TO 110
      D(2) = D(2)/GSMO(2)
      IF (K.EQ.4) GO TO 100
      DO 90 I=3,K2
        I1 = I - 1
        I2 = I1 - 1
        P = GSMO(I1)*B(I1)
        GSMO(I) = GSMO(I) - GSMO(I2)*(D(I2)**2) - P*B(I1)
        B(I) = (B(I)-P*D(I1))/GSMO(I)
        D(I) = D(I)/GSMO(I)
   90 CONTINUE
  100 P = GSMO(K2)*B(K2)
      GSMO(K1) = GSMO(K1) - GSMO(K3)*(D(K3)**2) - P*B(K2)
      B(K1) = (B(K1)-P*D(K2))/GSMO(K1)
  110 GSMO(K) = GSMO(K) - GSMO(K2)*(D(K2)**2) - GSMO(K1)*(B(K1)**2)
!  FACTORIZATION COMPLETE.
!  CARRY OUT FORWARD  AND BACKWARD SUBSTITUTION.
      C(2) = C(2) - B(1)*C(1)
      DO 120 I=3,K
        I1 = I - 1
        I2 = I - 2
        C(I) = C(I) - B(I1)*C(I1) - D(I2)*C(I2)
  120 CONTINUE
      DO 130 I=1,K
        C(I) = C(I)/GSMO(I)
  130 CONTINUE
      C(K1) = C(K1) - B(K1)*C(K)
      DO 140 I=2,K1
        J = K - I
        C(J) = C(J) - B(J)*C(J+1) - D(J)*C(J+2)
  140 CONTINUE
!  THE 5 BANDED SYSTEM HAS BEEN SOLVED.THE SOLUTION IS IN
!  ARRAY C. THE COEFFICIENTS GSMO, B, C, AND D WILL NOW BE
!  SET UP.
  150 C(N) = 0.0
      D(N) = 0.0
      C(N1) = C(N2)
      HK1 = TNODE(N) - TNODE(N1)
      D(N1) = -C(N1)/(3.0*HK1)
      GSMO(N) = G(N) + R6*D(N1)/WGS(N)
      IF (N.EQ.3) GO TO 170
      DO 160 I=2,N2
        K = N - I
        K1 = K + 1
        HK2 = HK1
        HK1 = TNODE(K1) - TNODE(K)
        C(K) = C(K-1)
        D(K) = (C(K1)-C(K))/(3.0*HK1)
        GSMO(K1) = G(K1) - R6*(D(K1)-D(K))/WGS(K1)
        B(K1) = (GSMO(K1+1)-GSMO(K1))/HK2 - HK2*(C(K1)+C(K1)+C(K1+1))/ &
          3.0
  160 CONTINUE
  170 C(1) = 0.0
      HK2 = HK1
      HK1 = TNODE(2) - TNODE(1)
      D(1) = (C(2)-C(1))/(3.0*HK1)
      GSMO(2) = G(2) - R6*(D(2)-D(1))/WGS(2)
      GSMO(1) = G(1) - R6*D(1)/WGS(1)
      B(2) = (GSMO(3)-GSMO(2))/HK2 - HK2*(C(2)+C(2)+C(3))/3.0
      B(1) = (GSMO(2)-GSMO(1))/HK1 - HK1*(C(1)+C(1)+C(2))/3.0
!  THE DISCRETE CUBIC SMOOTHING SPLINE IS NOW COMPLETE.
      RETURN
!  THE TRIVIAL CASE WHEN N=2 IS HANDLED HERE.
  180 GSMO(1) = G(1)
      GSMO(2) = G(2)
      B(1) = (G(2)-G(1))/(TNODE(2)-TNODE(1))
      C(1) = 0.0
      D(1) = 0.0
      RETURN
      END       SUBROUTINE DCSSMO



      SUBROUTINE DPPDI (AP, N, DET, JOB)
!***BEGIN PROLOGUE  DPPDI
!***PURPOSE  Compute the determinant and inverse of a real symmetric
!            positive definite matrix using factors from DPPCO or DPPFA.
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1B, D3B1B
!***TYPE      DOUBLE PRECISION (SPPDI-S, DPPDI-D, CPPDI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
!             PACKED, POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DPPDI computes the determinant and inverse
!     of a double precision symmetric positive definite matrix
!     using the factors computed by DPPCO or DPPFA .
!
!     On Entry
!
!        AP      DOUBLE PRECISION (N*(N+1)/2)
!                the output from DPPCO or DPPFA.
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        JOB     INTEGER
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     On Return
!
!        AP      the upper triangular half of the inverse .
!                The strict lower triangle is unaltered.
!
!        DET     DOUBLE PRECISION(2)
!                determinant of original matrix if requested.
!                Otherwise not referenced.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                with  1.0 .LE. DET(1) .LT. 10.0
!                or  DET(1) .EQ. 0.0 .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        It will not occur if the subroutines are called correctly
!        and if DPOCO or DPOFA has set INFO .EQ. 0 .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DSCAL
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DPPDI
      INTEGER(I4B) N,JOB
      REAL(DP) AP(*)
      REAL(DP) DET(2)
!
      REAL(DP) T
      REAL(DP) S
      INTEGER(I4B) I,II,J,JJ,JM1,J1,K,KJ,KK,KP1,K1
!***FIRST EXECUTABLE STATEMENT  DPPDI
!
!     COMPUTE DETERMINANT
!
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         S = 10.0D0
         II = 0
         DO 50 I = 1, N
            II = II + I
            DET(1) = AP(II)**2*DET(1)
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (DET(1) .GE. 1.0D0) GO TO 20
               DET(1) = S*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (DET(1) .LT. S) GO TO 40
               DET(1) = DET(1)/S
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
!
!     COMPUTE INVERSE(R)
!
      IF (MOD(JOB,10) .EQ. 0) GO TO 140
         KK = 0
         DO 100 K = 1, N
            K1 = KK + 1
            KK = KK + K
            AP(KK) = 1.0D0/AP(KK)
            T = -AP(KK)
            CALL DSCAL(K-1,T,AP(K1),1)
            KP1 = K + 1
            J1 = KK + 1
            KJ = KK + K
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = AP(KJ)
               AP(KJ) = 0.0D0
               CALL DAXPY(K,T,AP(K1),1,AP(J1),1)
               J1 = J1 + J
               KJ = KJ + J
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        FORM  INVERSE(R) * TRANS(INVERSE(R))
!
         JJ = 0
         DO 130 J = 1, N
            J1 = JJ + 1
            JJ = JJ + J
            JM1 = J - 1
            K1 = 1
            KJ = J1
            IF (JM1 .LT. 1) GO TO 120
            DO 110 K = 1, JM1
               T = AP(KJ)
               CALL DAXPY(K,T,AP(J1),1,AP(K1),1)
               K1 = K1 + K
               KJ = KJ + 1
  110       CONTINUE
  120       CONTINUE
            T = AP(JJ)
            CALL DSCAL(J,T,AP(J1),1)
  130    CONTINUE
  140 CONTINUE
      RETURN
      END      SUBROUTINE DPPDI



      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
!***BEGIN PROLOGUE  XERMSG
!***PURPOSE  Process error messages for SLATEC and other libraries.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERMSG-A)
!***KEYWORDS  ERROR MESSAGE, XERROR
!***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!***DESCRIPTION
!
!   XERMSG processes a diagnostic message in a manner determined by the
!   value of LEVEL and the current value of the library error control
!   flag, KONTRL.  See subroutine XSETF for details.
!
!    LIBRAR   A character constant (or character variable) with the name
!             of the library.  This will be 'SLATEC' for the SLATEC
!             Common Math Library.  The error handling package is
!             general enough to be used by many libraries
!             simultaneously, so it is desirable for the routine that
!             detects and reports an error to identify the library name
!             as well as the routine name.
!
!    SUBROU   A character constant (or character variable) with the name
!             of the routine that detected the error.  Usually it is the
!             name of the routine that is calling XERMSG.  There are
!             some instances where a user callable library routine calls
!             lower level subsidiary routines where the error is
!             detected.  In such cases it may be more informative to
!             supply the name of the routine the user called rather than
!             the name of the subsidiary routine that detected the
!             error.
!
!    MESSG    A character constant (or character variable) with the text
!             of the error or warning message.  In the example below,
!             the message is a character constant that contains a
!             generic message.
!
!                   CALL XERMSG ('SLATEC', 'MMPY',
!                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
!                  *3, 1)
!
!             It is possible (and is sometimes desirable) to generate a
!             specific message--e.g., one that contains actual numeric
!             values.  Specific numeric values can be converted into
!             character strings using formatted WRITE statements into
!             character variables.  This is called standard Fortran
!             internal file I/O and is exemplified in the first three
!             lines of the following example.  You can also catenate
!             substrings of characters to construct the error message.
!             Here is an example showing the use of both writing to
!             an internal file and catenating character strings.
!
!                   CHARACTER*5 CHARN, CHARL
!                   WRITE (CHARN,10) N
!                   WRITE (CHARL,10) LDA
!                10 FORMAT(I5)
!                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
!                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
!                  *   CHARL, 3, 1)
!
!             There are two subtleties worth mentioning.  One is that
!             the // for character catenation is used to construct the
!             error message so that no single character constant is
!             continued to the next line.  This avoids confusion as to
!             whether there are trailing blanks at the end of the line.
!             The second is that by catenating the parts of the message
!             as an actual argument rather than encoding the entire
!             message into one large character variable, we avoid
!             having to know how long the message will be in order to
!             declare an adequate length for that large character
!             variable.  XERMSG calls XERPRN to print the message using
!             multiple lines if necessary.  If the message is very long,
!             XERPRN will break it into pieces of 72 characters (as
!             requested by XERMSG) for printing on multiple lines.
!             Also, XERMSG asks XERPRN to prefix each line with ' *  '
!             so that the total line length could be 76 characters.
!             Note also that XERPRN scans the error message backwards
!             to ignore trailing blanks.  Another feature is that
!             the substring '$$' is treated as a new line sentinel
!             by XERPRN.  If you want to construct a multiline
!             message without having to count out multiples of 72
!             characters, just use '$$' as a separator.  '$$'
!             obviously must occur within 72 characters of the
!             start of each line to have its intended effect since
!             XERPRN is asked to wrap around at 72 characters in
!             addition to looking for '$$'.
!
!    NERR     An integer value that is chosen by the library routine's
!             author.  It must be in the range -99 to 999 (three
!             printable digits).  Each distinct error should have its
!             own error number.  These error numbers should be described
!             in the machine readable documentation for the routine.
!             The error numbers need be unique only within each routine,
!             so it is reasonable for each routine to start enumerating
!             errors from 1 and proceeding to the next integer.
!
!    LEVEL    An integer value in the range 0 to 2 that indicates the
!             level (severity) of the error.  Their meanings are
!
!            -1  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.  An attempt is made to only print this
!                message once.
!
!             0  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.
!
!             1  A recoverable error.  This is used even if the error is
!                so serious that the routine cannot return any useful
!                answer.  If the user has told the error package to
!                return after recoverable errors, then XERMSG will
!                return to the Library routine which can then return to
!                the user's routine.  The user may also permit the error
!                package to terminate the program upon encountering a
!                recoverable error.
!
!             2  A fatal error.  XERMSG will not return to its caller
!                after it receives a fatal error.  This level should
!                hardly ever be used; it is much better to allow the
!                user a chance to recover.  An example of one of the few
!                cases in which it is permissible to declare a level 2
!                error is a reverse communication Library routine that
!                is likely to be called repeatedly until it integrates
!                across some interval.  If there is a serious error in
!                the input such that another step cannot be taken and
!                the Library routine is called again without the input
!                error having been corrected by the caller, the Library
!                routine will probably be called forever with improper
!                input.  In this case, it is reasonable to declare the
!                error to be fatal.
!
!    Each of the arguments to XERMSG is input; none will be modified by
!    XERMSG.  A routine may make multiple calls to XERMSG with warning
!    level messages; however, after a call to XERMSG with a recoverable
!    error, the routine should return to the user.  Do not try to call
!    XERMSG with a second recoverable error after the first recoverable
!    error because the error package saves the error number.  The user
!    can retrieve this error number by calling another entry point in
!    the error handling package and then clear the error number when
!    recovering from the error.  Calling XERMSG in succession causes the
!    old error number to be overwritten by the latest error number.
!    This is considered harmless for error numbers associated with
!    warning messages but must not be done for error numbers of serious
!    errors.  After a call to XERMSG with a recoverable error, the user
!    must be given a chance to call NUMXER or XERCLR to retrieve or
!    clear the error number.
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
!***REVISION HISTORY  (YYMMDD)
!   880101  DATE WRITTEN
!   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
!           THERE ARE TWO BASIC CHANGES.
!           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
!               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
!               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
!               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
!               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
!               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
!               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
!               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
!           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
!               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
!               OF LOWER CASE.
!   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
!           THE PRINCIPAL CHANGES ARE
!           1.  CLARIFY COMMENTS IN THE PROLOGUES
!           2.  RENAME XRPRNT TO XERPRN
!           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
!               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
!               CHARACTER FOR NEW RECORDS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           CLEAN UP THE CODING.
!   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
!           PREFIX.
!   891013  REVISED TO CORRECT COMMENTS.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
!           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
!           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
!           XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERMSG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
!***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
!
!       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
!       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
!          SHOULD BE PRINTED.
!
!       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
!          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
!          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
!
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR. &
         LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' // &
            'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '// &
            'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
!
!       RECORD THE MESSAGE.
!
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
!
!       HANDLE PRINT-ONCE WARNING MESSAGES.
!
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
!
!       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
!
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
!
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
!
!       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
!       ZERO AND THE ERROR IS NOT FATAL.
!
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
!
!       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
!       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
!       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
!       IS NOT ZERO.
!
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
!
!       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
!       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
!       FROM EACH OF THE FOLLOWING THREE OPTIONS.
!       1.  LEVEL OF THE MESSAGE
!              'INFORMATIVE MESSAGE'
!              'POTENTIALLY RECOVERABLE ERROR'
!              'FATAL ERROR'
!       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
!              'PROG CONTINUES'
!              'PROG ABORTED'
!       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
!           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
!           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
!              'TRACEBACK REQUESTED'
!              'TRACEBACK NOT REQUESTED'
!       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
!       EXCEED 74 CHARACTERS.
!       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
!
      IF (LKNTRL .GT. 0) THEN
!
!       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
!
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
!
!       THEN WHETHER THE PROGRAM WILL CONTINUE.
!
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR. &
             (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
!
!       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
!
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
!
!       NOW SEND OUT THE MESSAGE.
!
      CALL XERPRN (' *  ', -1, MESSG, 72)
!
!       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
!          TRACEBACK.
!
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
!
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
!
!       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
!
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
!
!       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
!       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
!
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
!
!       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
!       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
!       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
!
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN &
               (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END       SUBROUTINE XERMSG





      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
!***BEGIN PROLOGUE  XERPRN
!***SUBSIDIARY
!***PURPOSE  Print error messages processed by XERMSG.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERPRN-A)
!***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
!***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!***DESCRIPTION
!
! This routine sends one or more lines to each of the (up to five)
! logical units to which error messages are to be sent.  This routine
! is called several times by XERMSG, sometimes with a single line to
! print and sometimes with a (potentially very long) message that may
! wrap around into multiple lines.
!
! PREFIX  Input argument of type CHARACTER.  This argument contains
!         characters to be put at the beginning of each line before
!         the body of the message.  No more than 16 characters of
!         PREFIX will be used.
!
! NPREF   Input argument of type INTEGER.  This argument is the number
!         of characters to use from PREFIX.  If it is negative, the
!         intrinsic function LEN is used to determine its length.  If
!         it is zero, PREFIX is not used.  If it exceeds 16 or if
!         LEN(PREFIX) exceeds 16, only the first 16 characters will be
!         used.  If NPREF is positive and the length of PREFIX is less
!         than NPREF, a copy of PREFIX extended with blanks to length
!         NPREF will be used.
!
! MESSG   Input argument of type CHARACTER.  This is the text of a
!         message to be printed.  If it is a long message, it will be
!         broken into pieces for printing on multiple lines.  Each line
!         will start with the appropriate prefix and be followed by a
!         piece of the message.  NWRAP is the number of characters per
!         piece; that is, after each NWRAP characters, we break and
!         start a new line.  In addition the characters '$$' embedded
!         in MESSG are a sentinel for a new line.  The counting of
!         characters up to NWRAP starts over for each new line.  The
!         value of NWRAP typically used by XERMSG is 72 since many
!         older error messages in the SLATEC Library are laid out to
!         rely on wrap-around every 72 characters.
!
! NWRAP   Input argument of type INTEGER.  This gives the maximum size
!         piece into which to break MESSG for printing on multiple
!         lines.  An embedded '$$' ends a line, and the count restarts
!         at the following character.  If a line break does not occur
!         on a blank (it would split a word) that word is moved to the
!         next line.  Values of NWRAP less than 16 will be treated as
!         16.  Values of NWRAP greater than 132 will be treated as 132.
!         The actual line length will be NPREF + NWRAP after NPREF has
!         been adjusted to fall between 0 and 16 and NWRAP has been
!         adjusted to fall between 16 and 132.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  I1MACH, XGETUA
!***REVISION HISTORY  (YYMMDD)
!   880621  DATE WRITTEN
!   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
!           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
!           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
!           SLASH CHARACTER IN FORMAT STATEMENTS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
!           LINES TO BE PRINTED.
!   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
!           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
!   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Added code to break messages between words.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER(I4B) NPREF, NWRAP
      CHARACTER*148 CBUFF
      INTEGER(I4B) IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
!***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
!
!       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
!       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
!       ERROR MESSAGE UNIT.
!
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
!
!       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
!       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
!       THE REST OF THIS ROUTINE.
!
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
!
!       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
!       TIME FROM MESSG TO PRINT ON ONE LINE.
!
      LWRAP = MAX(16, MIN(132, NWRAP))
!
!       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
!
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
!
!       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
!
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
!
!       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
!       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
!       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
!       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
!
!       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
!       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
!       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
!       OF THE SECOND ARGUMENT.
!
!       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
!       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
!       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
!       POSITION NEXTC.
!
!       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
!                       REMAINDER OF THE CHARACTER STRING.  LPIECE
!                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
!                       WHICHEVER IS LESS.
!
!       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
!                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
!                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
!                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
!                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
!                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
!                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
!                       SHOULD BE INCREMENTED BY 2.
!
!       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
!
!       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
!                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
!                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
!                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
!                       AT THE END OF A LINE.
!
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
!
!       THERE WAS NO NEW LINE SENTINEL FOUND.
!
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
!
!       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
!       DON'T PRINT A BLANK LINE.
!
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
!
!       LPIECE SHOULD BE SET DOWN TO LWRAP.
!
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
!
!       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
!       WE SHOULD DECREMENT LPIECE BY ONE.
!
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
!
!       PRINT
!
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
!
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
      END       SUBROUTINE XERPRN



      SUBROUTINE XERHLT (MESSG)
!***BEGIN PROLOGUE  XERHLT
!***SUBSIDIARY
!***PURPOSE  Abort program execution and print error message.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERHLT-A)
!***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        ***Note*** machine dependent routine
!        XERHLT aborts the execution of the program.
!        The error message causing the abort is given in the calling
!        sequence, in case one needs it for printing on a dayfile,
!        for example.
!
!     Description of Parameters
!        MESSG is as in XERMSG.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900206  Routine changed from user-callable to subsidiary.  (WRB)
!   900510  Changed calling sequence to delete length of character
!           and changed routine name from XERABT to XERHLT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
      END       SUBROUTINE XERHLT


      SUBROUTINE xerrwv (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
      INTEGER(I4B) msg, nmes, nerr, level, ni, i1, i2, nr, &
         i, lun, lunit, mesflg, ncpw, nch, nwds
      REAL(DP) r1, r2
      DIMENSION msg(nmes)
!-----------------------------------------------------------------------
! subroutines xerrwv, xsetf, and xsetun, as given here, constitute
! a simplified version of the slatec error handling package.
! written by a. c. hindmarsh at llnl.  version of march 30, 1987.
! this version is in double precision.
!
! all arguments are input arguments.
!
! msg    = the message (hollerith literal or integer array).
! nmes   = the length of msg (number of characters).
! nerr   = the error number (not used).
! level  = the error level..
!          0 or 1 means recoverable (control returns to caller).
!          2 means fatal (run is aborted--see note below).
! ni     = number of integers (0, 1, or 2) to be printed with message.
! i1,i2  = integers to be printed, depending on ni.
! nr     = number of reals (0, 1, or 2) to be printed with message.
! r1,r2  = reals to be printed, depending on nr.
!
! note..  this routine is machine-dependent and specialized for use
! in limited context, in the following ways..
! 1. the number of hollerith characters stored per word, denoted
!    by ncpw below, is a data-loaded constant.
! 2. the value of nmes is assumed to be at most 60.
!    (multi-line messages are generated by repeated calls.)
! 3. if level = 2, control passes to the statement   stop
!    to abort the run.  this statement may be machine-dependent.
! 4. r1 and r2 are assumed to be in double precision and are printed
!    in d21.13 format.
! 5. the common block /eh0001/ below is data-loaded (a machine-
!    dependent feature) with default values.
!    this block is needed for proper retention of parameters used by
!    this routine which the user can reset by calling xsetf or xsetun.
!    the variables in this block are as follows..
!       mesflg = print control flag..
!                1 means print all messages (the default).
!                0 means no printing.
!       lunit  = logical unit number for messages.
!                the default is 6 (machine-dependent).
!-----------------------------------------------------------------------
! the following are instructions for installing this routine
! in different machine environments.
!
! to change the default output unit, change the data statement
! in the block data subprogram below.
!
! for a different number of characters per word, change the
! data statement setting ncpw below, and format 10.  alternatives for
! various computers are shown in comment cards.
!
! for a different run-abort command, change the statement following
! statement 100 at the end.
!-----------------------------------------------------------------------
      COMMON /eh0001/ mesflg, lunit
!-----------------------------------------------------------------------
! the following data-loaded value of ncpw is valid for the cdc-6600
! and cdc-7600 computers.
!     data ncpw/10/
! the following is valid for the cray-1 computer.
!     data ncpw/8/
! the following is valid for the burroughs 6700 and 7800 computers.
!     data ncpw/6/
! the following is valid for the pdp-10 computer.
!     data ncpw/5/
! the following is valid for the vax computer with 4 bytes per integer,
! and for the ibm-360, ibm-370, ibm-303x, and ibm-43xx computers.
      DATA ncpw/4/
! the following is valid for the pdp-11, or vax with 2-byte integers.
!     data ncpw/2/
!-----------------------------------------------------------------------
      IF (mesflg .EQ. 0) go to 100
! get logical unit number. ---------------------------------------------
      lun = lunit
! get number of words in message. --------------------------------------
      nch = min0(nmes,60)
      nwds = nch/ncpw
      IF (nch .NE. nwds*ncpw) nwds = nwds + 1
! write the message. ---------------------------------------------------
      WRITE (lun, 10) (msg(i),i=1,nwds)
!-----------------------------------------------------------------------
! the following format statement is to have the form
! 10  format(1x,mmann)
! where nn = ncpw and mm is the smallest integer .ge. 60/ncpw.
! the following is valid for ncpw = 10.
! 10  format(1x,6a10)
! the following is valid for ncpw = 8.
! 10  format(1x,8a8)
! the following is valid for ncpw = 6.
! 10  format(1x,10a6)
! the following is valid for ncpw = 5.
! 10  format(1x,12a5)
! the following is valid for ncpw = 4.
  10  FORMAT(1x,15a4)
! the following is valid for ncpw = 2.
! 10  format(1x,30a2)
!-----------------------------------------------------------------------
      IF (ni .EQ. 1) WRITE (lun, 20) i1
 20   FORMAT(6x,23hin above message,  i1 =,i10)
      IF (ni .EQ. 2) WRITE (lun, 30) i1,i2
 30   FORMAT(6x,23hin above message,  i1 =,i10,3x,4hi2 =,i10)
      IF (nr .EQ. 1) WRITE (lun, 40) r1
 40   FORMAT(6x,23hin above message,  r1 =,d21.13)
      IF (nr .EQ. 2) WRITE (lun, 50) r1,r2
 50   FORMAT(6x,15hin above,  r1 =,d21.13,3x,4hr2 =,d21.13)
! abort the run if level = 2. ------------------------------------------
 100  IF (level .NE. 2) RETURN
      STOP
!----------------------- end of subroutine xerrwv ----------------------
      END       SUBROUTINE xerrwv

      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
!***BEGIN PROLOGUE  XERCNT
!***SUBSIDIARY
!***PURPOSE  Allow user control over handling of errors.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERCNT-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        Allows user control over handling of individual errors.
!        Just after each message is recorded, but before it is
!        processed any further (i.e., before it is printed or
!        a decision to abort is made), a call is made to XERCNT.
!        If the user has provided his own version of XERCNT, he
!        can then override the value of KONTROL used in processing
!        this message by redefining its value.
!        KONTRL may be set to any value from -2 to 2.
!        The meanings for KONTRL are the same as in XSETF, except
!        that the value of KONTRL changes only for this message.
!        If KONTRL is set to a value outside the range from -2 to 2,
!        it will be moved back into that range.
!
!     Description of Parameters
!
!      --Input--
!        LIBRAR - the library that the routine is in.
!        SUBROU - the subroutine that XERMSG is being called from
!        MESSG  - the first 20 characters of the error message.
!        NERR   - same as in the call to XERMSG.
!        LEVEL  - same as in the call to XERMSG.
!        KONTRL - the current value of the control flag as set
!                 by a call to XSETF.
!
!      --Output--
!        KONTRL - the new value of KONTRL.  If KONTRL is not
!                 defined, it will remain at its original value.
!                 This changed value of control affects only
!                 the current occurrence of the current message.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900206  Routine changed from user-callable to subsidiary.  (WRB)
!   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
!           names, changed routine name from XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERCNT
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
!***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END  SUBROUTINE XERCNT 





      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, &
         ICOUNT)
!***BEGIN PROLOGUE  XERSVE
!***SUBSIDIARY
!***PURPOSE  Record that an error has occurred.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (XERSVE-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
! *Usage:
!
!        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
!        CHARACTER * (len) LIBRAR, SUBROU, MESSG
!
!        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
!
! *Arguments:
!
!        LIBRAR :IN    is the library that the message is from.
!        SUBROU :IN    is the subroutine that the message is from.
!        MESSG  :IN    is the message to be saved.
!        KFLAG  :IN    indicates the action to be performed.
!                      when KFLAG > 0, the message in MESSG is saved.
!                      when KFLAG=0 the tables will be dumped and
!                      cleared.
!                      when KFLAG < 0, the tables will be dumped and
!                      not cleared.
!        NERR   :IN    is the error number.
!        LEVEL  :IN    is the error severity.
!        ICOUNT :OUT   the number of times this message has been seen,
!                      or zero if the table has overflowed and does not
!                      contain this message specifically.  When KFLAG=0,
!                      ICOUNT will not be altered.
!
! *Description:
!
!   Record that this error occurred and possibly dump and clear the
!   tables.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  I1MACH, XGETUA
!***REVISION HISTORY  (YYMMDD)
!   800319  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900413  Routine modified to remove reference to KFLAG.  (WRB)
!   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
!           sequence, use IF-THEN-ELSE, make number of saved entries
!           easily changeable, changed routine name from XERSAV to
!           XERSVE.  (RWC)
!   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERSVE
      PARAMETER (LENTAB=10)
      INTEGER(I4B) LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
!***FIRST EXECUTABLE STATEMENT  XERSVE
!
      IF (KFLAG.LE.0) THEN
!
!        Dump the table.
!
         IF (NMSG.EQ.0) RETURN
!
!        Print to each unit.
!
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
!
!           Print the table header.
!
            WRITE (IUNIT,9000)
!
!           Print body of table.
!
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I), &
                  NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
!
!           Print number of other errors.
!
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
!
!        Clear the error tables.
!
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
!
!        PROCESS A MESSAGE...
!        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
!
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND. &
               MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND. &
               LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
!
         IF (NMSG.LT.LENTAB) THEN
!
!           Empty slot found for new message.
!
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
!
!           Table is full.
!
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
!
!     Formats.
!
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' / &
         ' LIBRARY    SUBROUTINE MESSAGE START             NERR', &
         '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END       SUBROUTINE XERSVE




      SUBROUTINE XGETUA (IUNITA, N)
!***BEGIN PROLOGUE  XGETUA
!***PURPOSE  Return unit number(s) to which error messages are being
!            sent.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XGETUA-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        XGETUA may be called to determine the unit number or numbers
!        to which error messages are being sent.
!        These unit numbers may have been set by a call to XSETUN,
!        or a call to XSETUA, or may be a default value.
!
!     Description of Parameters
!      --Output--
!        IUNIT - an array of one to five unit numbers, depending
!                on the value of N.  A value of zero refers to the
!                default unit, as defined by the I1MACH machine
!                constant routine.  Only IUNIT(1),...,IUNIT(N) are
!                defined by XGETUA.  The values of IUNIT(N+1),...,
!                IUNIT(5) are not defined (for N .LT. 5) or altered
!                in any way by XGETUA.
!        N     - the number of units to which copies of the
!                error messages are being sent.  N will be in the
!                range from 1 to 5.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  J4SAVE
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
!***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END SUBROUTINE XGETUA


      INTEGER(I4B)  FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
!***BEGIN PROLOGUE  J4SAVE
!***SUBSIDIARY
!***PURPOSE  Save or recall global variables needed by error
!            handling routines.
!***LIBRARY   SLATEC (XERROR)
!***TYPE      INTEGER (J4SAVE-I)
!***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        J4SAVE saves and recalls several global variables needed
!        by the library error handling routines.
!
!     Description of Parameters
!      --Input--
!        IWHICH - Index of item desired.
!                = 1 Refers to current error number.
!                = 2 Refers to current error control flag.
!                = 3 Refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                = 4 Refers to the maximum number of times any
!                     message is to be printed (as set by XERMAX).
!                = 5 Refers to the total number of units to which
!                     each error message is to be written.
!                = 6 Refers to the 2nd unit for error messages
!                = 7 Refers to the 3rd unit for error messages
!                = 8 Refers to the 4th unit for error messages
!                = 9 Refers to the 5th unit for error messages
!        IVALUE - The value to be set for the IWHICH-th parameter,
!                 if ISET is .TRUE. .
!        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
!                 given the value, IVALUE.  If ISET=.FALSE., the
!                 IWHICH-th parameter will be unchanged, and IVALUE
!                 is a dummy parameter.
!      --Output--
!        The (old) value of the IWHICH-th parameter will be returned
!        in the function value, J4SAVE.
!
!***SEE ALSO  XERMSG
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900205  Minor modifications to prologue.  (WRB)
!   900402  Added TYPE section.  (WRB)
!   910411  Added KEYWORDS section.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
!***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END FUNCTION j4save


      SUBROUTINE FDUMP
!***BEGIN PROLOGUE  FDUMP
!***PURPOSE  Symbolic dump (should be locally written).
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (FDUMP-A)
!***KEYWORDS  ERROR, XERMSG
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!        ***Note*** Machine Dependent Routine
!        FDUMP is intended to be replaced by a locally written
!        version which produces a symbolic dump.  Failing this,
!        it should be replaced by a version which prints the
!        subprogram nesting list.  Note that this dump must be
!        printed on each of up to five files, as indicated by the
!        XGETUA routine.  See XSETUA and XGETUA for details.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  FDUMP
!***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END       SUBROUTINE FDUMP

      SUBROUTINE ewset (n, itol, rtol, atol, ycur, ewt)
!lll. optimize
!-----------------------------------------------------------------------
! this subroutine sets the error weight vector ewt according to
!     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
! with the subscript on rtol and/or atol possibly replaced by 1 above,
! depending on the value of itol.
!-----------------------------------------------------------------------
      INTEGER n, itol
      INTEGER i
      REAL(DP) rtol, atol, ycur, ewt
      DIMENSION rtol(1), atol(1), ycur(n), ewt(n)
!
      go to (10, 20, 30, 40), itol
 10   CONTINUE
      DO 15 i = 1,n
 15     ewt(i) = rtol(1)*dabs(ycur(i)) + atol(1)
      RETURN
 20   CONTINUE
      DO 25 i = 1,n
 25     ewt(i) = rtol(1)*dabs(ycur(i)) + atol(i)
      RETURN
 30   CONTINUE
      DO 35 i = 1,n
 35     ewt(i) = rtol(i)*dabs(ycur(i)) + atol(1)
      RETURN
 40   CONTINUE
      DO 45 i = 1,n
 45     ewt(i) = rtol(i)*dabs(ycur(i)) + atol(i)
      RETURN
!----------------------- end of subroutine ewset -----------------------
      END SUBROUTINE ewset 




      SUBROUTINE iprep (neq, y, rwork, ia, ja, ipflag, f, jac)
!lll. optimize
      EXTERNAL f, jac
      INTEGER neq, ia, ja, ipflag
      INTEGER illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, &
         mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      INTEGER icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
      INTEGER i, imax, lewtn, lyhd, lyhn ,iwkint(1)
      REAL(DP) y, rwork
      REAL(DP) rowns, &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      REAL(DP) rlss
      DIMENSION neq(1), y(1), rwork(1), ia(1), ja(1)
      COMMON /ls0001/ rowns(209), &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, &
         illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, &
         mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6), &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      COMMON /lss001/ rlss(6), &
         iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
!-----------------------------------------------------------------------
! this routine serves as an interface between the driver and
! subroutine prep.  it is called only if miter is 1 or 2.
! tasks performed here are..
!  * call prep,
!  * reset the required wm segment length lenwk,
!  * move yh back to its final location (following wm in rwork),
!  * reset pointers for yh, savf, ewt, and acor, and
!  * move ewt to its new position if istate = 1.
! ipflag is an output error indication flag.  ipflag = 0 if there was
! no trouble, and ipflag is the value of the prep error flag ipper
! if there was trouble in subroutine prep.
!-----------------------------------------------------------------------
      ipflag = 0
! call prep to do matrix preprocessing operations. ---------------------
      iwkint(1) = INT(rwork(lwm))
      !CALL prep (neq, y, rwork(lyh), rwork(lsavf), rwork(lewt), &
      !   rwork(lacor), ia, ja, rwork(lwm), rwork(lwm), ipflag, f, jac)
      CALL prep (neq, y, rwork(lyh), rwork(lsavf), rwork(lewt), &
         rwork(lacor), ia, ja, rwork(lwm), iwkint(1), ipflag, f, jac)
      lenwk = max0(lreq,lwmin)
      IF (ipflag .LT. 0) RETURN
! if prep was successful, move yh to end of required space for wm. -----
      lyhn = lwm + lenwk
      IF (lyhn .GT. lyh) RETURN
      lyhd = lyh - lyhn
      IF (lyhd .EQ. 0) go to 20
      imax = lyhn - 1 + lenyhm
      DO 10 i = lyhn,imax
 10     rwork(i) = rwork(i+lyhd)
      lyh = lyhn
! reset pointers for savf, ewt, and acor. ------------------------------
 20   lsavf = lyh + lenyh
      lewtn = lsavf + n
      lacor = lewtn + n
      IF (istatc .EQ. 3) go to 40
! if istate = 1, move ewt (left) to its new position. ------------------
      IF (lewtn .GT. lewt) RETURN
      DO 30 i = 1,n
 30     rwork(i+lewtn-1) = rwork(i+lewt-1)
 40   lewt = lewtn
      RETURN
!----------------------- end of subroutine iprep -----------------------
      END       SUBROUTINE iprep

      REAL(DP) FUNCTION vnorm (n, v, w)
!lll. optimize
!-----------------------------------------------------------------------
! this function routine computes the weighted root-mean-square norm
! of the vector of length n contained in the array v, with weights
! contained in the array w of length n..
!   vnorm = sqrt( (1/n) * sum( v(i)*w(i) )**2 )
!-----------------------------------------------------------------------
      INTEGER n,   i
      REAL(DP) v, w,   sum
      DIMENSION v(n), w(n)
      sum = 0.0d0
      DO 10 i = 1,n
 10     sum = sum + (v(i)*w(i))**2
      vnorm = dsqrt(sum/dfloat(n))
      RETURN
!----------------------- end of function vnorm -------------------------
      END FUNCTION vnorm 

      SUBROUTINE intdy (t, k, yh, nyh, dky, iflag)
!lll. optimize
      INTEGER k, nyh, iflag
      INTEGER iownd, iowns, &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER i, ic, j, jb, jb2, jj, jj1, jp1
      REAL(DP) t, yh, dky
      REAL(DP) rowns, &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      REAL(DP) c, r, s, tp
      DIMENSION yh(nyh,1), dky(1)
      COMMON /ls0001/ rowns(209), &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, &
         iownd(14), iowns(6), &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
!-----------------------------------------------------------------------
! intdy computes interpolated values of the k-th derivative of the
! dependent variable vector y, and stores it in dky.  this routine
! is called within the package with k = 0 and t = tout, but may
! also be called by the user for any k up to the current order.
! (see detailed instructions in the usage documentation.)
!-----------------------------------------------------------------------
! the computed values in dky are gotten by interpolation using the
! nordsieck history array yh.  this array corresponds uniquely to a
! vector-valued polynomial of degree nqcur or less, and dky is set
! to the k-th derivative of this polynomial at t.
! the formula for dky is..
!              q
!  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
!             j=k
! where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
! the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are
! communicated by common.  the above sum is done in reverse order.
! iflag is returned negative if either k or t is out of bounds.
!-----------------------------------------------------------------------
      iflag = 0
      IF (k .LT. 0 .OR. k .GT. nq) go to 80
      tp = tn - hu -  100.0d0*uround*(tn + hu)
      IF ((t-tp)*(t-tn) .GT. 0.0d0) go to 90
!
      s = (t - tn)/h
      ic = 1
      IF (k .EQ. 0) go to 15
      jj1 = l - k
      DO 10 jj = jj1,nq
 10     ic = ic*jj
 15   c = dfloat(ic)
      DO 20 i = 1,n
 20     dky(i) = c*yh(i,l)
      IF (k .EQ. nq) go to 55
      jb2 = nq - k
      DO 50 jb = 1,jb2
        j = nq - jb
        jp1 = j + 1
        ic = 1
        IF (k .EQ. 0) go to 35
        jj1 = jp1 - k
        DO 30 jj = jj1,j
 30       ic = ic*jj
 35     c = dfloat(ic)
        DO 40 i = 1,n
 40       dky(i) = c*yh(i,jp1) + s*dky(i)
 50     CONTINUE
      IF (k .EQ. 0) RETURN
 55   r = h**(-k)
      DO 60 i = 1,n
 60     dky(i) = r*dky(i)
      RETURN
!
 80    CONTINUE
!       call xerrwv(30hintdy--  k (=i1) illegal      , &
!         30, 51, 0, 1, k, 0, 0, 0.0d0, 0.0d0)
      iflag = -1
      RETURN
 90   CONTINUE
!4324
!     CALL xerrwv(30hintdy--  t (=r1) illegal      , &
!         30, 52, 0, 0, 0, 0, 1, t, 0.0d0)
!      call xerrwv( &
!        60h      t not in interval tcur - hu (= r1) to tcur (=r2)      , &
!         60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2
      RETURN
!----------------------- end of subroutine intdy -----------------------
      END SUBROUTINE intdy

      SUBROUTINE prjs (neq,y,yh,nyh,ewt,ftem,savf,wk,iwk,f,jac)
!lll. optimize
      EXTERNAL f,jac
      INTEGER neq, nyh, iwk
      INTEGER iownd, iowns, &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
      INTEGER i, imul, j, jj, jok, jmax, jmin, k, kmax, kmin, ng
      REAL(DP) y, yh, ewt, ftem, savf, wk
      REAL(DP) rowns, &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      REAL(DP) con0, conmin, ccmxj, psmall, rbig, seth
      REAL(DP) con, di, fac, hl0, pij, r, r0, rcon, rcont, &
         srur
!     dimension neq(1), y(1), yh(nyh,1), ewt(1), ftem(1), savf(1),
      DIMENSION neq(1), y(1), yh(nyh,*), ewt(1), ftem(1), savf(1), &
         wk(*), iwk(1)
!    1   wk(1), iwk(1)
      COMMON /ls0001/ rowns(209), &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, &
         iownd(14), iowns(6), &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      COMMON /lss001/ con0, conmin, ccmxj, psmall, rbig, seth, &
         iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
!-----------------------------------------------------------------------
! prjs is called to compute and process the matrix
! p = i - h*el(1)*j , where j is an approximation to the jacobian.
! j is computed by columns, either by the user-supplied routine jac
! if miter = 1, or by finite differencing if miter = 2.
! if miter = 3, a diagonal approximation to j is used.
! if miter = 1 or 2, and if the existing value of the jacobian
! (as contained in p) is considered acceptable, then a new value of
! p is reconstructed from the old value.  in any case, when miter
! is 1 or 2, the p matrix is subjected to lu decomposition in cdrv.
! p and its lu decomposition are stored (separately) in wk.
!
! in addition to variables described previously, communication
! with prjs uses the following..
! y     = array containing predicted values on entry.
! ftem  = work array of length n (acor in stode).
! savf  = array containing f evaluated at predicted y.
! wk    = real work space for matrices.  on output it contains the
!         inverse diagonal matrix if miter = 3, and p and its sparse
!         lu decomposition if miter is 1 or 2.
!         storage of matrix elements starts at wk(3).
!         wk also contains the following matrix-related data..
!         wk(1) = sqrt(uround), used in numerical jacobian increments.
!         wk(2) = h*el0, saved for later use if miter = 3.
! iwk   = integer work space for matrix-related data, assumed to
!         be equivalenced to wk.  in addition, wk(iprsp) and iwk(ipisp)
!         are assumed to have identical locations.
! el0   = el(1) (input).
! ierpj = output error flag (in common).
!       = 0 if no error.
!       = 1  if zero pivot found in cdrv.
!       = 2  if a singular matrix arose with miter = 3.
!       = -1 if insufficient storage for cdrv (should not occur here).
!       = -2 if other error found in cdrv (should not occur here).
! jcur  = output flag = 1 to indicate that the jacobian matrix
!         (or approximation) is now current.
! this routine also uses other variables in common.
!-----------------------------------------------------------------------
      hl0 = h*el0
      con = -hl0
      IF (miter .EQ. 3) go to 300
! see whether j should be reevaluated (jok = 0) or not (jok = 1). ------
      jok = 1
      IF (nst .EQ. 0 .OR. nst .GE. nslj+msbj) jok = 0
      IF (icf .EQ. 1 .AND. dabs(rc - 1.0d0) .LT. ccmxj) jok = 0
      IF (icf .EQ. 2) jok = 0
      IF (jok .EQ. 1) go to 250
!
! miter = 1 or 2, and the jacobian is to be reevaluated. ---------------
 20   jcur = 1
      nje = nje + 1
      nslj = nst
      iplost = 0
      conmin = dabs(con)
      go to (100, 200), miter
!
! if miter = 1, call jac, multiply by scalar, and add identity. --------
 100  CONTINUE
      kmin = iwk(ipian)
      DO 130 j = 1, n
        kmax = iwk(ipian+j) - 1
        DO 110 i = 1,n
 110      ftem(i) = 0.0d0
        CALL jac (neq, tn, y, j, iwk(ipian), iwk(ipjan), ftem)
        DO 120 k = kmin, kmax
          i = iwk(ibjan+k)
          wk(iba+k) = ftem(i)*con
          IF (i .EQ. j) wk(iba+k) = wk(iba+k) + 1.0d0
 120      CONTINUE
        kmin = kmax + 1
 130    CONTINUE
      go to 290
!
! if miter = 2, make ngp calls to f to approximate j and p. ------------
 200  CONTINUE
      fac = vnorm(n, savf, ewt)
      r0 = 1000.0d0 * dabs(h) * uround * dfloat(n) * fac
      IF (r0 .EQ. 0.0d0) r0 = 1.0d0
      srur = wk(1)
      jmin = iwk(ipigp)
      DO 240 ng = 1,ngp
        jmax = iwk(ipigp+ng) - 1
        DO 210 j = jmin,jmax
          jj = iwk(ibjgp+j)
          r = dmax1(srur*dabs(y(jj)),r0/ewt(jj))
 210      y(jj) = y(jj) + r
        CALL f (neq, tn, y, ftem)
        DO 230 j = jmin,jmax
          jj = iwk(ibjgp+j)
          y(jj) = yh(jj,1)
          r = dmax1(srur*dabs(y(jj)),r0/ewt(jj))
          fac = -hl0/r
          kmin =iwk(ibian+jj)
          kmax =iwk(ibian+jj+1) - 1
          DO 220 k = kmin,kmax
            i = iwk(ibjan+k)
            wk(iba+k) = (ftem(i) - savf(i))*fac
            IF (i .EQ. jj) wk(iba+k) = wk(iba+k) + 1.0d0
 220        CONTINUE
 230      CONTINUE
        jmin = jmax + 1
 240    CONTINUE
      nfe = nfe + ngp
      go to 290
!
! if jok = 1, reconstruct new p from old p. ----------------------------
 250  jcur = 0
      rcon = con/con0
      rcont = dabs(con)/conmin
      IF (rcont .GT. rbig .AND. iplost .EQ. 1) go to 20
      kmin = iwk(ipian)
      DO 275 j = 1,n
        kmax = iwk(ipian+j) - 1
        DO 270 k = kmin,kmax
          i = iwk(ibjan+k)
          pij = wk(iba+k)
          IF (i .NE. j) go to 260
          pij = pij - 1.0d0
          IF (dabs(pij) .GE. psmall) go to 260
            iplost = 1
            conmin = dmin1(dabs(con0),conmin)
 260      pij = pij*rcon
          IF (i .EQ. j) pij = pij + 1.0d0
          wk(iba+k) = pij
 270      CONTINUE
        kmin = kmax + 1
 275    CONTINUE
!
! do numerical factorization of p matrix. ------------------------------
 290  nlu = nlu + 1
      con0 = con
      ierpj = 0
      DO 295 i = 1,n
 295    ftem(i) = 0.0d0
      CALL cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan), &
         wk(ipa),ftem,ftem,nsp,iwk(ipisp),wk(iprsp),iesp,2,iys)
      IF (iys .EQ. 0) RETURN
      imul = (iys - 1)/n
      ierpj = -2
      IF (imul .EQ. 8) ierpj = 1
      IF (imul .EQ. 10) ierpj = -1
      RETURN
!
! if miter = 3, construct a diagonal approximation to j and p. ---------
 300  CONTINUE
      jcur = 1
      nje = nje + 1
      wk(2) = hl0
      ierpj = 0
      r = el0*0.1d0
      DO 310 i = 1,n
 310    y(i) = y(i) + r*(h*savf(i) - yh(i,2))
      CALL f (neq, tn, y, wk(3))
      nfe = nfe + 1
      DO 320 i = 1,n
        r0 = h*savf(i) - yh(i,2)
        di = 0.1d0*r0 - h*(wk(i+2) - savf(i))
        wk(i+2) = 1.0d0
        IF (dabs(r0) .LT. uround/ewt(i)) go to 320
        IF (dabs(di) .EQ. 0.0d0) go to 330
        wk(i+2) = 0.1d0*r0/di
 320    CONTINUE
      RETURN
 330  ierpj = 2
      RETURN
!----------------------- end of subroutine prjs ------------------------
      END SUBROUTINE prjs

      SUBROUTINE slss (wk, iwk, x, tem)
!lll. optimize
      INTEGER iwk
      INTEGER iownd, iowns, &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
      INTEGER i
      DOUBLE PRECISION wk, x, tem
      REAL(DP) rowns, &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      REAL(DP) rlss
      REAL(DP) di, hl0, phl0, r
!     dimension wk(1), iwk(1), x(1), tem(1)
      DIMENSION wk(*), iwk(*), x(*), tem(*)
      COMMON /ls0001/ rowns(209), &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, &
         iownd(14), iowns(6), &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      COMMON /lss001/ rlss(6), &
         iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
!-----------------------------------------------------------------------
! this routine manages the solution of the linear system arising from
! a chord iteration.  it is called if miter .ne. 0.
! if miter is 1 or 2, it calls cdrv to accomplish this.
! if miter = 3 it updates the coefficient h*el0 in the diagonal
! matrix, and then computes the solution.
! communication with slss uses the following variables..
! wk    = real work space containing the inverse diagonal matrix if
!         miter = 3 and the lu decomposition of the matrix otherwise.
!         storage of matrix elements starts at wk(3).
!         wk also contains the following matrix-related data..
!         wk(1) = sqrt(uround) (not used here),
!         wk(2) = hl0, the previous value of h*el0, used if miter = 3.
! iwk   = integer work space for matrix-related data, assumed to
!         be equivalenced to wk.  in addition, wk(iprsp) and iwk(ipisp)
!         are assumed to have identical locations.
! x     = the right-hand side vector on input, and the solution vector
!         on output, of length n.
! tem   = vector of work space of length n, not used in this version.
! iersl = output flag (in common).
!         iersl = 0  if no trouble occurred.
!         iersl = -1 if cdrv returned an error flag (miter = 1 or 2).
!                    this should never occur and is considered fatal.
!         iersl = 1  if a singular matrix arose with miter = 3.
! this routine also uses other variables in common.
!-----------------------------------------------------------------------
      iersl = 0
      go to (100, 100, 300), miter
 100  CALL cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan), &
         wk(ipa),x,x,nsp,iwk(ipisp),wk(iprsp),iesp,4,iersl)
      IF (iersl .NE. 0) iersl = -1
      RETURN
!
 300  phl0 = wk(2)
      hl0 = h*el0
      wk(2) = hl0
      IF (hl0 .EQ. phl0) go to 330
      r = hl0/phl0
      DO 320 i = 1,n
        di = 1.0d0 - r*(1.0d0 - 1.0d0/wk(i+2))
        IF (dabs(di) .EQ. 0.0d0) go to 390
 320    wk(i+2) = 1.0d0/di
 330  DO 340 i = 1,n
 340    x(i) = wk(i+2)*x(i)
      RETURN
 390  iersl = 1
      RETURN
!
!----------------------- end of subroutine slss ------------------------
      END SUBROUTINE slss

      SUBROUTINE stode (neq, y, yh, nyh, yh1, ewt, savf, acor, &
         wm, iwm, f, jac, pjac, slvs)
!lll. optimize
      EXTERNAL f, jac, pjac, slvs
      INTEGER neq, nyh, iwm
      INTEGER iownd, ialth, ipup, lmax, meo, nqnyh, nslp, &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER i, i1, iredo, iret, j, jb, m, ncf, newq
      REAL(DP) y, yh, yh1, ewt, savf, acor, wm
      REAL(DP) conit, crate, el, elco, hold, rmax, tesco, &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      REAL(DP) dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup, &
         r, rh, rhdn, rhsm, rhup, told
      DIMENSION neq(1), y(1), yh(nyh,*), yh1(1), ewt(1), savf(1), &
         acor(1), wm(1), iwm(1)
      COMMON /ls0001/ conit, crate, el(13), elco(13,12), &
         hold, rmax, tesco(3,12), &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(14), &
         ialth, ipup, lmax, meo, nqnyh, nslp, &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
!-----------------------------------------------------------------------
! stode performs one step of the integration of an initial value
! problem for a system of ordinary differential equations.
! note.. stode is independent of the value of the iteration method
! indicator miter, when this is .ne. 0, and hence is independent
! of the type of chord method used, or the jacobian structure.
! communication with stode is done with the following variables..
!
! neq    = integer array containing problem size in neq(1), and
!          passed as the neq argument in all calls to f and jac.
! y      = an array of length .ge. n used as the y argument in
!          all calls to f and jac.
! yh     = an nyh by lmax array containing the dependent variables
!          and their approximate scaled derivatives, where
!          lmax = maxord + 1.  yh(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by h**j/factorial(j)
!          (j = 0,1,...,nq).  on entry for the first step, the first
!          two columns of yh must be set from the initial values.
! nyh    = a constant integer .ge. n, the first dimension of yh.
! yh1    = a one-dimensional array occupying the same space as yh.
! ewt    = an array of length n containing multiplicative weights
!          for local error measurements.  local errors in y(i) are
!          compared to 1.0/ewt(i) in various error tests.
! savf   = an array of working storage, of length n.
!          also used for input of yh(*,maxord+2) when jstart = -1
!          and maxord .lt. the current order nq.
! acor   = a work array of length n, used for the accumulated
!          corrections.  on a successful return, acor(i) contains
!          the estimated one-step local error in y(i).
! wm,iwm = real and integer work arrays associated with matrix
!          operations in chord iteration (miter .ne. 0).
! pjac   = name of routine to evaluate and preprocess jacobian matrix
!          and p = i - h*el0*jac, if a chord method is being used.
! slvs   = name of routine to solve linear system in chord iteration.
! ccmax  = maximum relative change in h*el0 before pjac is called.
! h      = the step size to be attempted on the next step.
!          h is altered by the error control algorithm during the
!          problem.  h can be either positive or negative, but its
!          sign must remain constant throughout the problem.
! hmin   = the minimum absolute value of the step size h to be used.
! hmxi   = inverse of the maximum absolute value of h to be used.
!          hmxi = 0.0 is allowed and corresponds to an infinite hmax.
!          hmin and hmxi may be changed at any time, but will not
!          take effect until the next change of h is considered.
! tn     = the independent variable. tn is updated on each step taken.
! jstart = an integer used for input only, with the following
!          values and meanings..
!               0  perform the first step.
!           .gt.0  take a new step continuing from the last.
!              -1  take the next step with a new value of h, maxord,
!                    n, meth, miter, and/or matrix parameters.
!              -2  take the next step with a new value of h,
!                    but with other inputs unchanged.
!          on return, jstart is set to 1 to facilitate continuation.
! kflag  = a completion code with the following meanings..
!               0  the step was succesful.
!              -1  the requested error could not be achieved.
!              -2  corrector convergence could not be achieved.
!              -3  fatal error in pjac or slvs.
!          a return with kflag = -1 or -2 means either
!          abs(h) = hmin or 10 consecutive failures occurred.
!          on a return with kflag negative, the values of tn and
!          the yh array are as of the beginning of the last
!          step, and h is the last step size attempted.
! maxord = the maximum order of integration method to be allowed.
! maxcor = the maximum number of corrector iterations allowed.
! msbp   = maximum number of steps between pjac calls (miter .gt. 0).
! mxncf  = maximum number of convergence failures allowed.
! meth/miter = the method flags.  see description in driver.
! n      = the number of first-order differential equations.
!-----------------------------------------------------------------------
      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0d0
      IF (jstart .GT. 0) go to 200
      IF (jstart .EQ. -1) go to 100
      IF (jstart .EQ. -2) go to 160
!-----------------------------------------------------------------------
! on the first call, the order is set to 1, and other variables are
! initialized.  rmax is the maximum ratio by which h can be increased
! in a single step.  it is initially 1.e4 to compensate for the small
! initial h, but then is normally equal to 10.  if a failure
! occurs (in corrector convergence or error test), rmax is set at 2
! for the next increase.
!-----------------------------------------------------------------------
      lmax = maxord + 1
      nq = 1
      l = 2
      ialth = 2
      rmax = 10000.0d0
      rc = 0.0d0
      el0 = 1.0d0
      crate = 0.7d0
      hold = h
      meo = meth
      nslp = 0
      ipup = miter
      iret = 3
      go to 140
!-----------------------------------------------------------------------
! the following block handles preliminaries needed when jstart = -1.
! ipup is set to miter to force a matrix update.
! if an order increase is about to be considered (ialth = 1),
! ialth is reset to 2 to postpone consideration one more step.
! if the caller has changed meth, cfode is called to reset
! the coefficients of the method.
! if the caller has changed maxord to a value less than the current
! order nq, nq is reduced to maxord, and a new h chosen accordingly.
! if h is to be changed, yh must be rescaled.
! if h or meth is being changed, ialth is reset to l = nq + 1
! to prevent further changes in h for that many steps.
!-----------------------------------------------------------------------
 100  ipup = miter
      lmax = maxord + 1
      IF (ialth .EQ. 1) ialth = 2
      IF (meth .EQ. meo) go to 110
      CALL cfode (meth, elco, tesco)
      meo = meth
      IF (nq .GT. maxord) go to 120
      ialth = l
      iret = 1
      go to 150
 110  IF (nq .LE. maxord) go to 160
 120  nq = maxord
      l = lmax
      DO 125 i = 1,l
 125    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5d0/dfloat(nq+2)
      ddn = vnorm (n, savf, ewt)/tesco(1,l)
      exdn = 1.0d0/dfloat(l)
      rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
      rh = dmin1(rhdn,1.0d0)
      iredo = 3
      IF (h .EQ. hold) go to 170
      rh = dmin1(rh,dabs(h/hold))
      h = hold
      go to 175
!-----------------------------------------------------------------------
! cfode is called to get all the integration coefficients for the
! current meth.  then the el vector and related constants are reset
! whenever the order nq is changed, or at the start of the problem.
!-----------------------------------------------------------------------
 140  CALL cfode (meth, elco, tesco)
 150  DO 155 i = 1,l
 155    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5d0/dfloat(nq+2)
      go to (160, 170, 200), iret
!-----------------------------------------------------------------------
! if h is being changed, the h ratio rh is checked against
! rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to
! l = nq + 1 to prevent a change of h for that many steps, unless
! forced by a convergence or error test failure.
!-----------------------------------------------------------------------
 160  IF (h .EQ. hold) go to 200
      rh = h/hold
      h = hold
      iredo = 3
      go to 175
 170  rh = dmax1(rh,hmin/dabs(h))
 175  rh = dmin1(rh,rmax)
      rh = rh/dmax1(1.0d0,dabs(h)*hmxi*rh)
      r = 1.0d0
      DO 180 j = 2,l
        r = r*rh
        DO 180 i = 1,n
 180      yh(i,j) = yh(i,j)*r
      h = h*rh
      rc = rc*rh
      ialth = l
      IF (iredo .EQ. 0) go to 690
!-----------------------------------------------------------------------
! this section computes the predicted values by effectively
! multiplying the yh array by the pascal triangle matrix.
! rc is the ratio of new to old values of the coefficient  h*el(1).
! when rc differs from 1 by more than ccmax, ipup is set to miter
! to force pjac to be called, if a jacobian is involved.
! in any case, pjac is called at least every msbp steps.
!-----------------------------------------------------------------------
 200  IF (dabs(rc-1.0d0) .GT. ccmax) ipup = miter
      IF (nst .GE. nslp+msbp) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      DO 215 jb = 1,nq
        i1 = i1 - nyh
!dir$ ivdep
        DO 210 i = i1,nqnyh
 210      yh1(i) = yh1(i) + yh1(i+nyh)
 215    CONTINUE
!-----------------------------------------------------------------------
! up to maxcor corrector iterations are taken.  a convergence test is
! made on the r.m.s. norm of each correction, weighted by the error
! weight vector ewt.  the sum of the corrections is accumulated in the
! vector acor(i).  the yh array is not altered in the corrector loop.
!-----------------------------------------------------------------------
 220  m = 0
      DO 230 i = 1,n
 230    y(i) = yh(i,1)
      CALL f (neq, tn, y, savf)
      nfe = nfe + 1
      IF (ipup .LE. 0) go to 250
!-----------------------------------------------------------------------
! if indicated, the matrix p = i - h*el(1)*j is reevaluated and
! preprocessed before starting the corrector iteration.  ipup is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
      CALL pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac)
      ipup = 0
      rc = 1.0d0
      nslp = nst
      crate = 0.7d0
      IF (ierpj .NE. 0) go to 430
 250  DO 260 i = 1,n
 260    acor(i) = 0.0d0
 270  IF (miter .NE. 0) go to 350
!-----------------------------------------------------------------------
! in the case of functional iteration, update y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
      DO 290 i = 1,n
        savf(i) = h*savf(i) - yh(i,2)
 290    y(i) = savf(i) - acor(i)
      del = vnorm (n, y, ewt)
      DO 300 i = 1,n
        y(i) = yh(i,1) + el(1)*savf(i)
 300    acor(i) = savf(i)
      go to 400
!-----------------------------------------------------------------------
! in the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! p as coefficient matrix.
!-----------------------------------------------------------------------
 350  DO 360 i = 1,n
 360    y(i) = h*savf(i) - (yh(i,2) + acor(i))
      CALL slvs (wm, iwm, y, savf)
      IF (iersl .LT. 0) go to 430
      IF (iersl .GT. 0) go to 410
      del = vnorm (n, y, ewt)
      DO 380 i = 1,n
        acor(i) = acor(i) + y(i)
 380    y(i) = yh(i,1) + el(1)*acor(i)
!-----------------------------------------------------------------------
! test for convergence.  if m.gt.0, an estimate of the convergence
! rate constant is stored in crate, and this is used in the test.
!-----------------------------------------------------------------------
 400  IF (m .NE. 0) crate = dmax1(0.2d0*crate,del/delp)
      dcon = del*dmin1(1.0d0,1.5d0*crate)/(tesco(2,nq)*conit)
      IF (dcon .LE. 1.0d0) go to 450
      m = m + 1
      IF (m .EQ. maxcor) go to 410
      IF (m .GE. 2 .AND. del .GT. 2.0d0*delp) go to 410
      delp = del
      CALL f (neq, tn, y, savf)
      nfe = nfe + 1
      go to 270
!-----------------------------------------------------------------------
! the corrector iteration failed to converge.
! if miter .ne. 0 and the jacobian is out of date, pjac is called for
! the next try.  otherwise the yh array is retracted to its values
! before prediction, and h is reduced, if possible.  if h cannot be
! reduced or mxncf failures have occurred, exit with kflag = -2.
!-----------------------------------------------------------------------
 410  IF (miter .EQ. 0 .OR. jcur .EQ. 1) go to 430
      icf = 1
      ipup = miter
      go to 220
 430  icf = 2
      ncf = ncf + 1
      rmax = 2.0d0
      tn = told
      i1 = nqnyh + 1
      DO 445 jb = 1,nq
        i1 = i1 - nyh
!dir$ ivdep
        DO 440 i = i1,nqnyh
 440      yh1(i) = yh1(i) - yh1(i+nyh)
 445    CONTINUE
      IF (ierpj .LT. 0 .OR. iersl .LT. 0) go to 680
      IF (dabs(h) .LE. hmin*1.00001d0) go to 670
      IF (ncf .EQ. mxncf) go to 670
      rh = 0.25d0
      ipup = miter
      iredo = 1
      go to 170
!-----------------------------------------------------------------------
! the corrector has converged.  jcur is set to 0
! to signal that the jacobian involved may need updating later.
! the local error test is made and control passes to statement 500
! if it fails.
!-----------------------------------------------------------------------
 450  jcur = 0
      IF (m .EQ. 0) dsm = del/tesco(2,nq)
      IF (m .GT. 0) dsm = vnorm (n, acor, ewt)/tesco(2,nq)
      IF (dsm .GT. 1.0d0) go to 500
!-----------------------------------------------------------------------
! after a successful step, update the yh array.
! consider changing h if ialth = 1.  otherwise decrease ialth by 1.
! if ialth is then 1 and nq .lt. maxord, then acor is saved for
! use in a possible order increase on the next step.
! if a change in h is considered, an increase or decrease in order
! by one is considered also.  a change in h is made only if it is by a
! factor of at least 1.1.  if not, ialth is set to 3 to prevent
! testing for that many steps.
!-----------------------------------------------------------------------
      kflag = 0
      iredo = 0
      nst = nst + 1
      hu = h
      nqu = nq
      DO 470 j = 1,l
        DO 470 i = 1,n
 470      yh(i,j) = yh(i,j) + el(j)*acor(i)
      ialth = ialth - 1
      IF (ialth .EQ. 0) go to 520
      IF (ialth .GT. 1) go to 700
      IF (l .EQ. lmax) go to 700
      DO 490 i = 1,n
 490    yh(i,lmax) = acor(i)
      go to 700
!-----------------------------------------------------------------------
! the error test failed.  kflag keeps track of multiple failures.
! restore tn and the yh array to their previous values, and prepare
! to try the step again.  compute the optimum step size for this or
! one lower order.  after 2 or more failures, h is forced to decrease
! by a factor of 0.2 or less.
!-----------------------------------------------------------------------
 500  kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      DO 515 jb = 1,nq
        i1 = i1 - nyh
!dir$ ivdep
        DO 510 i = i1,nqnyh
 510      yh1(i) = yh1(i) - yh1(i+nyh)
 515    CONTINUE
      rmax = 2.0d0
      IF (dabs(h) .LE. hmin*1.00001d0) go to 660
      IF (kflag .LE. -3) go to 640
      iredo = 2
      rhup = 0.0d0
      go to 540
!-----------------------------------------------------------------------
! regardless of the success or failure of the step, factors
! rhdn, rhsm, and rhup are computed, by which h could be multiplied
! at order nq - 1, order nq, or order nq + 1, respectively.
! in the case of failure, rhup = 0.0 to avoid an order increase.
! the largest of these is determined and the new order chosen
! accordingly.  if the order is to be increased, we compute one
! additional scaled derivative.
!-----------------------------------------------------------------------
 520  rhup = 0.0d0
      IF (l .EQ. lmax) go to 540
      DO 530 i = 1,n
 530    savf(i) = acor(i) - yh(i,lmax)
      dup = vnorm (n, savf, ewt)/tesco(3,nq)
      exup = 1.0d0/dfloat(l+1)
      rhup = 1.0d0/(1.4d0*dup**exup + 0.0000014d0)
 540  exsm = 1.0d0/dfloat(l)
      rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
      rhdn = 0.0d0
      IF (nq .EQ. 1) go to 560
      ddn = vnorm (n, yh(1,l), ewt)/tesco(1,nq)
      exdn = 1.0d0/dfloat(nq)
      rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
 560  IF (rhsm .GE. rhup) go to 570
      IF (rhup .GT. rhdn) go to 590
      go to 580
 570  IF (rhsm .LT. rhdn) go to 580
      newq = nq
      rh = rhsm
      go to 620
 580  newq = nq - 1
      rh = rhdn
      IF (kflag .LT. 0 .AND. rh .GT. 1.0d0) rh = 1.0d0
      go to 620
 590  newq = l
      rh = rhup
      IF (rh .LT. 1.1d0) go to 610
      r = el(l)/dfloat(l)
      DO 600 i = 1,n
 600    yh(i,newq+1) = acor(i)*r
      go to 630
 610  ialth = 3
      go to 700
 620  IF ((kflag .EQ. 0) .AND. (rh .LT. 1.1d0)) go to 610
      IF (kflag .LE. -2) rh = dmin1(rh,0.2d0)
!-----------------------------------------------------------------------
! if there is a change of order, reset nq, l, and the coefficients.
! in any case h is reset according to rh and the yh array is rescaled.
! then exit from 690 if the step was ok, or redo the step otherwise.
!-----------------------------------------------------------------------
      IF (newq .EQ. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 150
!-----------------------------------------------------------------------
! control reaches this section if 3 or more failures have occured.
! if 10 failures have occurred, exit with kflag = -1.
! it is assumed that the derivatives that have accumulated in the
! yh array have errors of the wrong order.  hence the first
! derivative is recomputed, and the order is set to 1.  then
! h is reduced by a factor of 10, and the step is retried,
! until it succeeds or h reaches hmin.
!-----------------------------------------------------------------------
 640  IF (kflag .EQ. -10) go to 660
      rh = 0.1d0
      rh = dmax1(hmin/dabs(h),rh)
      h = h*rh
      DO 645 i = 1,n
 645    y(i) = yh(i,1)
      CALL f (neq, tn, y, savf)
      nfe = nfe + 1
      DO 650 i = 1,n
 650    yh(i,2) = h*savf(i)
      ipup = miter
      ialth = 5
      IF (nq .EQ. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 150
!-----------------------------------------------------------------------
! all returns are made through this section.  h is saved in hold
! to allow the caller to change h on the next step.
!-----------------------------------------------------------------------
 660  kflag = -1
      go to 720
 670  kflag = -2
      go to 720
 680  kflag = -3
      go to 720
 690  rmax = 10.0d0
 700  r = 1.0d0/tesco(2,nqu)
      DO 710 i = 1,n
 710    acor(i) = acor(i)*r
 720  hold = h
      jstart = 1
      RETURN
!----------------------- end of subroutine stode -----------------------
      END SUBROUTINE stode




      SUBROUTINE cfode (meth, elco, tesco)
!lll. optimize
      INTEGER meth
      INTEGER i, ib, nq, nqm1, nqp1
      REAL*8 elco, tesco
      REAL*8 agamq, fnq, fnqm1, pc, pint, ragq, &
         rqfac, rq1fac, tsign, xpin
      DIMENSION elco(13,12), tesco(3,12)
!-----------------------------------------------------------------------
! cfode is called by the integrator routine to set coefficients
! needed there.  the coefficients for the current method, as
! given by the value of meth, are set for all orders and saved.
! the maximum order assumed here is 12 if meth = 1 and 5 if meth = 2.
! (a smaller value of the maximum order is also allowed.)
! cfode is called once at the beginning of the problem,
! and is not called again unless and until meth is changed.
!
! the elco array contains the basic method coefficients.
! the coefficients el(i), 1 .le. i .le. nq+1, for the method of
! order nq are stored in elco(i,nq).  they are given by a genetrating
! polynomial, i.e.,
!     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
! for the implicit adams methods, l(x) is given by
!     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
! for the bdf methods, l(x) is given by
!     l(x) = (x+1)*(x+2)* ... *(x+nq)/k,
! where         k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
!
! the tesco array contains test constants used for the
! local error test and the selection of step size and/or order.
! at order nq, tesco(k,nq) is used for the selection of step
! size at order nq - 1 if k = 1, at order nq if k = 2, and at order
! nq + 1 if k = 3.
!-----------------------------------------------------------------------
      DIMENSION pc(12)
!
      go to (100, 200), meth
!
 100  elco(1,1) = 1.0d0
      elco(2,1) = 1.0d0
      tesco(1,1) = 0.0d0
      tesco(2,1) = 2.0d0
      tesco(1,2) = 1.0d0
      tesco(3,12) = 0.0d0
      pc(1) = 1.0d0
      rqfac = 1.0d0
      DO 140 nq = 2,12
!-----------------------------------------------------------------------
! the pc array will contain the coefficients of the polynomial
!     p(x) = (x+1)*(x+2)*...*(x+nq-1).
! initially, p(x) = 1.
!-----------------------------------------------------------------------
        rq1fac = rqfac
        rqfac = rqfac/dfloat(nq)
        nqm1 = nq - 1
        fnqm1 = dfloat(nqm1)
        nqp1 = nq + 1
! form coefficients of p(x)*(x+nq-1). ----------------------------------
        pc(nq) = 0.0d0
        DO 110 ib = 1,nqm1
          i = nqp1 - ib
 110      pc(i) = pc(i-1) + fnqm1*pc(i)
        pc(1) = fnqm1*pc(1)
! compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        pint = pc(1)
        xpin = pc(1)/2.0d0
        tsign = 1.0d0
        DO 120 i = 2,nq
          tsign = -tsign
          pint = pint + tsign*pc(i)/dfloat(i)
 120      xpin = xpin + tsign*pc(i)/dfloat(i+1)
! store coefficients in elco and tesco. --------------------------------
        elco(1,nq) = pint*rq1fac
        elco(2,nq) = 1.0d0
        DO 130 i = 2,nq
 130      elco(i+1,nq) = rq1fac*pc(i)/dfloat(i)
        agamq = rqfac*xpin
        ragq = 1.0d0/agamq
        tesco(2,nq) = ragq
        IF (nq .LT. 12) tesco(1,nqp1) = ragq*rqfac/dfloat(nqp1)
        tesco(3,nqm1) = ragq
 140    CONTINUE
      RETURN
!
 200  pc(1) = 1.0d0
      rq1fac = 1.0d0
      DO 230 nq = 1,5
!-----------------------------------------------------------------------
! the pc array will contain the coefficients of the polynomial
!     p(x) = (x+1)*(x+2)*...*(x+nq).
! initially, p(x) = 1.
!-----------------------------------------------------------------------
        fnq = dfloat(nq)
        nqp1 = nq + 1
! form coefficients of p(x)*(x+nq). ------------------------------------
        pc(nqp1) = 0.0d0
        DO 210 ib = 1,nq
          i = nq + 2 - ib
 210      pc(i) = pc(i-1) + fnq*pc(i)
        pc(1) = fnq*pc(1)
! store coefficients in elco and tesco. --------------------------------
        DO 220 i = 1,nqp1
 220      elco(i,nq) = pc(i)/pc(2)
        elco(2,nq) = 1.0d0
        tesco(1,nq) = rq1fac
        tesco(2,nq) = dfloat(nqp1)/elco(1,nq)
        tesco(3,nq) = dfloat(nq+2)/elco(1,nq)
        rq1fac = rq1fac/fnq
 230    CONTINUE
      RETURN
!----------------------- end of subroutine cfode -----------------------
      END SUBROUTINE cfode




      SUBROUTINE prep (neq, y, yh, savf, ewt, ftem, ia, ja, &
                           wk, iwk, ipper, f, jac)
!lll. optimize
      EXTERNAL f,jac
      INTEGER neq, ia, ja, iwk, ipper
      INTEGER iownd, iowns, &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
      INTEGER i, ibr, ier, ipil, ipiu, iptt1, iptt2, j, jfound, k, &
         knew, kmax, kmin, ldif, lenigp, liwk, maxg, np1, nzsut
      REAL*8 y, yh, savf, ewt, ftem, wk
      REAL*8 rowns, &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      REAL*8 con0, conmin, ccmxj, psmall, rbig, seth
      REAL*8 dq, dyj, erwt, fac, yj
      DIMENSION neq(1), y(1), yh(1), savf(1), ewt(1), ftem(1), &
         ia(1), ja(1), wk(1), iwk(1)
      COMMON /ls0001/ rowns(209), &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, &
         iownd(14), iowns(6), &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      COMMON /lss001/ con0, conmin, ccmxj, psmall, rbig, seth, &
         iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
!-----------------------------------------------------------------------
! this routine performs preprocessing related to the sparse linear
! systems that must be solved if miter = 1 or 2.
! the operations that are performed here are..
!  * compute sparseness structure of jacobian according to moss,
!  * compute grouping of column indices (miter = 2),
!  * compute a new ordering of rows and columns of the matrix,
!  * reorder ja corresponding to the new ordering,
!  * perform a symbolic lu factorization of the matrix, and
!  * set pointers for segments of the iwk/wk array.
! in addition to variables described previously, prep uses the
! following for communication..
! yh     = the history array.  only the first column, containing the
!          current y vector, is used.  used only if moss .ne. 0.
! savf   = a work array of length neq, used only if moss .ne. 0.
! ewt    = array of length neq containing (inverted) error weights.
!          used only if moss = 2 or if istate = moss = 1.
! ftem   = a work array of length neq, identical to acor in the driver,
!          used only if moss = 2.
! wk     = a real work array of length lenwk, identical to wm in
!          the driver.
! iwk    = integer work array, assumed to occupy the same space as wk.
! lenwk  = the length of the work arrays wk and iwk.
! istatc = a copy of the driver input argument istate (= 1 on the
!          first call, = 3 on a continuation call).
! iys    = flag value from odrv or cdrv.
! ipper  = output error flag with the following values and meanings..
!          0  no error.
!         -1  insufficient storage for internal structure pointers.
!         -2  insufficient storage for jgroup.
!         -3  insufficient storage for odrv.
!         -4  other error flag from odrv (should never occur).
!         -5  insufficient storage for cdrv.
!         -6  other error flag from cdrv.
!-----------------------------------------------------------------------
      ibian = lrat*2
      ipian = ibian + 1
      np1 = n + 1
      ipjan = ipian + np1
      ibjan = ipjan - 1
      liwk = lenwk*lrat
      IF (ipjan+n-1 .GT. liwk) go to 210
      IF (moss .EQ. 0) go to 30
!
      IF (istatc .EQ. 3) go to 20
! istate = 1 and moss .ne. 0.  perturb y for structure determination. --
      DO 10 i = 1,n
        erwt = 1.0d0/ewt(i)
        fac = 1.0d0 + 1.0d0/(dfloat(i)+1.0d0)
        y(i) = y(i) + fac*dsign(erwt,y(i))
 10     CONTINUE
      go to (70, 100), moss
!
 20   CONTINUE
! istate = 3 and moss .ne. 0.  load y from yh(*,1). --------------------
      DO 25 i = 1,n
 25     y(i) = yh(i)
      go to (70, 100), moss
!
! moss = 0.  process user-s ia,ja.  add diagonal entries if necessary. -
 30   knew = ipjan
      kmin = ia(1)
      iwk(ipian) = 1
      DO 60 j = 1,n
        jfound = 0
        kmax = ia(j+1) - 1
        IF (kmin .GT. kmax) go to 45
        DO 40 k = kmin,kmax
          i = ja(k)
          IF (i .EQ. j) jfound = 1
          IF (knew .GT. liwk) go to 210
          iwk(knew) = i
          knew = knew + 1
 40       CONTINUE
        IF (jfound .EQ. 1) go to 50
 45     IF (knew .GT. liwk) go to 210
        iwk(knew) = j
        knew = knew + 1
 50     iwk(ipian+j) = knew + 1 - ipjan
        kmin = kmax + 1
 60     CONTINUE
      go to 140
!
! moss = 1.  compute structure from user-supplied jacobian routine jac.
 70   CONTINUE
! a dummy call to f allows user to create temporaries for use in jac. --
      CALL f (neq, tn, y, savf)
      k = ipjan
      iwk(ipian) = 1
      DO 90 j = 1,n
        IF (k .GT. liwk) go to 210
        iwk(k) = j
        k = k + 1
        DO 75 i = 1,n
 75       savf(i) = 0.0d0
        CALL jac (neq, tn, y, j, iwk(ipian), iwk(ipjan), savf)
        DO 80 i = 1,n
          IF (dabs(savf(i)) .LE. seth) go to 80
          IF (i .EQ. j) go to 80
          IF (k .GT. liwk) go to 210
          iwk(k) = i
          k = k + 1
 80       CONTINUE
        iwk(ipian+j) = k + 1 - ipjan
 90     CONTINUE
      go to 140
!
! moss = 2.  compute structure from results of n + 1 calls to f. -------
 100  k = ipjan
      iwk(ipian) = 1
      CALL f (neq, tn, y, savf)
      DO 120 j = 1,n
        IF (k .GT. liwk) go to 210
        iwk(k) = j
        k = k + 1
        yj = y(j)
        erwt = 1.0d0/ewt(j)
        dyj = dsign(erwt,yj)
        y(j) = yj + dyj
        CALL f (neq, tn, y, ftem)
        y(j) = yj
        DO 110 i = 1,n
          dq = (ftem(i) - savf(i))/dyj
          IF (dabs(dq) .LE. seth) go to 110
          IF (i .EQ. j) go to 110
          IF (k .GT. liwk) go to 210
          iwk(k) = i
          k = k + 1
 110      CONTINUE
        iwk(ipian+j) = k + 1 - ipjan
 120    CONTINUE
!
 140  CONTINUE
      IF (moss .EQ. 0 .OR. istatc .NE. 1) go to 150
! if istate = 1 and moss .ne. 0, restore y from yh. --------------------
      DO 145 i = 1,n
 145    y(i) = yh(i)
 150  nnz = iwk(ipian+n) - 1
      lenigp = 0
      ipigp = ipjan + nnz
      IF (miter .NE. 2) go to 160
!
! compute grouping of column indices (miter = 2). ----------------------
      maxg = np1
      ipjgp = ipjan + nnz
      ibjgp = ipjgp - 1
      ipigp = ipjgp + n
      iptt1 = ipigp + np1
      iptt2 = iptt1 + n
      lreq = iptt2 + n - 1
      IF (lreq .GT. liwk) go to 220
      CALL jgroup (n, iwk(ipian), iwk(ipjan), maxg, ngp, iwk(ipigp), &
         iwk(ipjgp), iwk(iptt1), iwk(iptt2), ier)
      IF (ier .NE. 0) go to 220
      lenigp = ngp + 1
!
! compute new ordering of rows/columns of jacobian. --------------------
 160  ipr = ipigp + lenigp
      ipc = ipr
      ipic = ipc + n
      ipisp = ipic + n
      iprsp = (ipisp - 2)/lrat + 2
      iesp = lenwk + 1 - iprsp
      IF (iesp .LT. 0) go to 230
      ibr = ipr - 1
      DO 170 i = 1,n
 170    iwk(ibr+i) = i
      nsp = liwk + 1 - ipisp
      CALL odrv (n, iwk(ipian), iwk(ipjan), wk, iwk(ipr), iwk(ipic), &
         nsp, iwk(ipisp), 1, iys)
      IF (iys .EQ. 11*n+1) go to 240
      IF (iys .NE. 0) go to 230
!
! reorder jan and do symbolic lu factorization of matrix. --------------
      ipa = lenwk + 1 - nnz
      nsp = ipa - iprsp
      lreq = max0(12*n/lrat, 6*n/lrat+2*n+nnz) + 3
      lreq = lreq + iprsp - 1 + nnz
      IF (lreq .GT. lenwk) go to 250
      iba = ipa - 1
      DO 180 i = 1,nnz
 180    wk(iba+i) = 0.0d0
      ipisp = lrat*(iprsp - 1) + 1
      CALL cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan), &
         wk(ipa),wk(ipa),wk(ipa),nsp,iwk(ipisp),wk(iprsp),iesp,5,iys)
      lreq = lenwk - iesp
      IF (iys .EQ. 10*n+1) go to 250
      IF (iys .NE. 0) go to 260
      ipil = ipisp
      ipiu = ipil + 2*n + 1
      nzu = iwk(ipil+n) - iwk(ipil)
      nzl = iwk(ipiu+n) - iwk(ipiu)
      IF (lrat .GT. 1) go to 190
      CALL adjlr (n, iwk(ipisp), ldif)
      lreq = lreq + ldif
 190  CONTINUE
      IF (lrat .EQ. 2 .AND. nnz .EQ. n) lreq = lreq + 1
      nsp = nsp + lreq - lenwk
      ipa = lreq + 1 - nnz
      iba = ipa - 1
      ipper = 0
      RETURN
!
 210  ipper = -1
      lreq = 2 + (2*n + 1)/lrat
      lreq = max0(lenwk+1,lreq)
      RETURN
!
 220  ipper = -2
      lreq = (lreq - 1)/lrat + 1
      RETURN
!
 230  ipper = -3
      CALL cntnzu (n, iwk(ipian), iwk(ipjan), nzsut)
      lreq = lenwk - iesp + (3*n + 4*nzsut - 1)/lrat + 1
      RETURN
!
 240  ipper = -4
      RETURN
!
 250  ipper = -5
      RETURN
!
 260  ipper = -6
      lreq = lenwk
      RETURN
!----------------------- end of subroutine prep ------------------------
      END SUBROUTINE prep
  
   END MODULE replace_imsl
