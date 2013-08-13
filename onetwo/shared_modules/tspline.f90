     MODULE tension_spline
       ! note this  module also has non tension spline items in it !!

       USE nrtype,                       ONLY : DP,I4B
       USE io_gcnmp,                     ONLY : ncrt,nlog

       USE error_handler
       USE io_gcnmp,                     ONLY : nlog,ncrt


       IMPLICIT NONE
       
       !yp,sigma,wk are used to pass spline fit data from sub 
       !t716_TSPSI  to sub t716_TSVAL1 in this module
       REAL(DP),ALLOCATABLE,DIMENSION(:),PRIVATE :: yp,sigma,wk,w
       REAL(DP),ALLOCATABLE, DIMENSION(:) :: s_coef
       REAL(DP) bc_0,bc_1,rho_w,sp,dsp,d2sp
       INTEGER(I4B) s_bc
       LOGICAL get_scoef




       CONTAINS 

       SUBROUTINE tspline90 (x, y, nx, bpar, rgrid, tspl, npts,tmax,t)
!-------------------------------------------------------------------------
!
! --- tspline calculates the second derivatives of the spline at the knots.
! --- these coefficients define the spline and are stored
! --- in array cs(i,j),i = 1,#knots,j=1,3
! --- subroutine EVTSPLN is used to evaluate the spline,
!
! --- input
!
!      x      vector,length nx,of knot locations (must be in ascending order)
!      y          function values at the knots
!      n          #points in x and y
!    rgrid(i)     i = 1,2...npts values of x at which evaluation is desired
!      bpar       array used to specifiy boundary conditions
!      ic         exact row dimension of matrix c
!      t          tension parameter,(optional)
!
!      a,b,c,     temporary work vectors of length nx
!      fppr,dx
!
! --- output
!
!      cs         array of spline coefficients (see above)
! --- the value of the spline for x in the interval (x(i),x(i+1)) is
! ---           (the range on i is 1 to nx-1 )
!          f(x) = cs(i,1)*cs(i,2) * SINH (t*(x(i+1)-x))+
!               (y(i)*t*t-cs(i,1))*cs(i,3)*(x(i+1)-x)
!               +cs(i+1,1)*cs(i,2) * SINH (t*(x-x(i)))+
!               (y(i+1)*t*t-cs(i+1,1))*cs(j,3)*(x-x(i))
!       tmax   max allowed value of tension parameter
!       tspl(i)      i = 1,2..npts values of y at points  rgrid(i)
!
! ------------------------------------------------------------------ HSJ
!
      INTEGER(I4b),INTENT(IN)             :: npts,nx
      REAL(DP),INTENT(inout),DIMENSION(:) :: rgrid,x,y,bpar
      REAL(DP),INTENT(inout),OPTIONAL     :: t
      REAL(DP),INTENT(out)                :: tmax
      REAL(DP),INTENT(out),DIMENSION(:)   :: tspl

      LOGICAL  reverseset
      REAL(DP),DIMENSION(:,:)             :: cs(SIZE(x),3)
      REAL(DP),DIMENSION(:)               :: a(SIZE(x)),b(SIZE(x)),c(SIZE(x)),    &
                                             fpp(SIZE(x)),r(SIZE(x)),dx(SIZE(x))
      REAL(DP) dxmin,bp1,bp2,bp3,bp4,tsq,a1,b1,an,bn,const,         &
               bcl,xhold,yhold, ac,pv,thn,th2
      INTEGER(I4b) ier,nm1,j,i
      const = .0194444444_DP
!
! --- check input
!
      cs(:,:) = 0.0_DP
      reverseset = .FALSE.
      dxmin = ABS (x(nx)-x(1))
      nm1 = nx-1
      ier = 0
      IF (nx .LT. 2)  ier = 130
      IF (ier .NE. 0)  go to 1000
      DO 15 j=1,nm1
        dxmin = MIN (dxmin, ABS (x(j+1)-x(j)))
        IF (x(j) .LT. x(j+1))  go to 15
        ier = 131
        go to 1000
   15 CONTINUE
!
! --- some initialization
!
   16 bp1   = bpar(1)
      bp2   = bpar(2)
      bp3   = bpar(3)
      bp4   = bpar(4)
      tmax  = 30.0_DP / dxmin    ! max tension, avoids overflow
      IF(PRESENT(t))THEN
         t     = MIN (t, tmax)
      ELSE
         t     = 0.0_DP           !default to zero tension
      ENDIF
      tsq   = t*t
      nm1   = nx-1
      dx(1) = x(2) - x(1)

! --- set up vectors a,b,c
! --- vector a is lower diagonal
! --- vector b is diagonal
! --- vector c is upper diagonal
!
      IF (t .GT. 0.0_DP)THEN    
         DO 20 j=2,nm1
            dx(j) = x(j+1)-x(j)
            a(j) = (1.0_DP/dx(j-1)-t / SINH (t*dx(j-1)))/tsq
            b(j) = t * COSH (t*dx(j-1)) / SINH (t*dx(j-1))
            b(j) = b(j)-1.0_DP/dx(j-1)-1.0_DP/dx(j)
            b(j) = (b(j)+t * COSH (t*dx(j)) / SINH (t*dx(j)))/tsq
            c(j) = (1.0_DP/dx(j)-t / SINH (t*dx(j)))/tsq
20          r(j) = (y(j+1)-y(j))/dx(j)-(y(j)-y(j-1))/dx(j-1)
            !
            ! --- boundary conditions on spline
            !
            th2 = t*dx(1)
            thn = t*dx(nm1)
            b1 = (1.0_DP/t)*(1.0_DP / SINH (th2)-1.0_DP/th2)
            a1 = -1.0_DP/TANH(th2)+1.0_DP/th2
            a1 = a1/(1.0_DP / SINH (th2)-1.0_DP/th2)
            an = (1.0_DP/t)*(-1.0_DP / SINH (thn)+1.0_DP/thn)
            bn = 1.0_DP / TANH(thn)-1.0_DP/thn
            bn = bn/(-1.0_DP / SINH (thn)+1.0_DP/thn)
            IF (bpar(1) .NE. -1.0D30)  go to 30
            bp1 = 1.0_DP
            !***    bp2 = (-1.0_DP/(dx(1)*b1))*(dx(1)*bpar(2)+y(1)-y(2))
            bp2 = (+1.0_DP/(dx(1)*b1))*(dx(1)*bpar(2)+y(1)-y(2))
30          IF (bpar(3) .NE. -1.0D30)  go to 40
            bp3 = 1.0_DP
            bp4 = (1.0_DP/(dx(nm1)*an))*(dx(nm1)*bpar(4)+y(nm1)-y(nx))
40          b(1) = a1
            c(1) = bp1
            a(nx) = bp3
            b(nx) = bn
            r(1) = bp2
            r(nx) = bp4
!
! --- for tension sufficiently small go over into cubic spline
!
      ELSE ! (here t =0.0)
         DO  j=2,nm1
            dx(j) = x(j+1)-x(j)
            a(j) = dx(j-1)/6.0_DP
            b(j) = (dx(j)+dx(j-1))/3.0_DP
            c(j) = dx(j)/6.0_DP
            r(j) = (y(j+1)-y(j))/dx(j)-(y(j)-y(j-1))/dx(j-1)
         ENDDO
            !
            ! --- boundary conditions
            !
            IF (bpar(1) .NE. -1.0D30)  go to 60
            bp1 = 1.0_DP
            bcl = dx(1)*(-0.166666667_DP+const*dx(1)**2*tsq)
            bp2 = (1.0_DP/(dx(1)*bcl))*(y(1)-y(2)+dx(1)*bpar(2))
60          IF (bpar(3) .NE. -1.0D30)  go to 70
            bp3 = 1.0_DP
            ac = dx(nm1)*(0.166666667_DP-1.94444444*dx(nm1)**2*tsq)
            bp4 = (1.0_DP/(dx(nm1)*ac))*(dx(nm1)*bpar(4)+y(nm1)-y(nx))
70          a(nx) = bp3
            b( 1) = 2.0_DP + dx(1  )**2*tsq / 10.0_DP
            c( 1) = bp1
            b(nx) = 2.0_DP + dx(nm1)**2*tsq / 10.0_DP
            r( 1) = bp2
            r(nx) = bp4
            !
            ! --- note a(1) and c(nx) are not used
            ! --- r is vector of rhs
            !
            ! --- solve the tridiagonal system
            ! --- forward elimination
            !
      ENDIF
      DO j=2,nx
        pv = a(j)/b(j-1)
        b(j) = b(j)-c(j-1)*pv
        r(j) = r(j)-r(j-1)*pv
      END DO
!
! --- back substitution
!
      fpp(nx) = r(nx)/b(nx)
      DO j=1,nm1
        i = nx-j
        fpp(i) = (r(i)-c(i)*fpp(i+1))/b(i)
      END DO
!
! --- now have vector of second derivatives fpp. set up the vector cs(i,j)
! --- to be used in evaluating the tension spline.
! --- convention used for natural spline will not work here.
! --- however we still pack all the necessary information into array cs.
!
      DO 150 j=1,nx
        cs(j,1) = fpp(j)
        IF (j .EQ. nx )  go to 150
        IF (t .EQ. 0.0_DP)  go to 150
        cs(j,2) = 1.0_DP / (SINH (t*dx(j))*t*t)
        cs(j,3) = 1.0_DP / (dx(j)*t*t)
  150 CONTINUE
!
      CALL evtspl90 (x, y, nx, cs, rgrid, tspl, npts, t)
      IF (.NOT. reverseset)  RETURN
!
 1000 DO j=1,nx/2
        xhold     = x(j)
        x(j)      = x(nx-j+1)
        x(nx-j+1) = xhold
        yhold     = y(j)
        y(j)      = y(nx-j+1)
        y(nx-j+1) = yhold
      END DO
!
      reverseset  = .NOT. reverseset
      IF (reverseset)  go to 16
      RETURN
!
      END       SUBROUTINE tspline90



      SUBROUTINE evtspl90 (x, y, nx, cs, r, tspl, npts, t)
! ----------------------------------------------------------------------
! --- evaluate the tension spline at npts points of r, return results in tspl
! ----------------------------------------------------------------------
!
      REAL(DP), INTENT(IN) ,DIMENSION(:)   ::   x,y,r
      REAL(DP), INTENT(IN) ,DIMENSION(:,:) ::   cs(SIZE(x),3)
      REAL(DP), INTENT(IN)                 ::   t
      INTEGER , INTENT(IN)                 ::   nx,npts
      
      REAL(DP), INTENT(OUT),DIMENSION(:)   ::   tspl(SIZE(x))

      INTEGER(I4B) j,i
      REAL(DP) z,dxl,dxr,delx


      DO j=1,npts
        z = r(j)
        DO i=1,nx-1
          IF (z .LE. x(i+1))  go to 30
        END DO
        IF (ABS ((z-x(nx))/x(nx)) .LT. 1.0e-8)  go to 30
!
        STOP 'subroutine EVTSPL: unspecified problem'
!
!       evaluate the spline
!
   30   dxl = z-x(i)
        dxr = x(i+1)-z
        IF (t .NE. 0.0_DP) THEN
             tspl(j) = cs(i,1)*cs(i,2) * SINH (t*dxr)+ &
                       (y(i)*t*t-cs(i,1))*cs(i,3)*dxr &
                      +cs(i+1,1)*cs(i,2) * SINH (t*dxl)+ &
                       (y(i+1)*t*t-cs(i+1,1))*cs(i,3)*dxl
        ELSE
             delx    = x(i+1)-x(i)
             tspl(j) = (dxr*(-cs(i  ,1) * dxl*(dxr+delx)+6.0_DP * y(i  )) &
                      + dxl*(-cs(i+1,1) * dxr*(dxl+delx)+6.0_DP * y(i+1))) &
                       /(6.0_DP * delx)
        END IF
      END DO
      RETURN
!
      END SUBROUTINE evtspl90 

      
      SUBROUTINE t716_TSPSI (n,x, y, ncd ,iendc,per,UNIFRM,deriv_bc,t)
!-------------------------------------------------------------------------
! --- this subroutine  is a wrapper for toms 716 tspline package
! --- routine TSPSI  (all toms tspline routines are in splinelib.a)
!
!
! --- input
!
! On input:
!
!       N = Number of data points.  N .GE. 2 and N .GE. 3 if
!           PER = TRUE.
!
!       X = Array of length N containing a strictly in-
!           creasing sequence of abscissae:  X(I) < X(I+1)
!           for I = 1,...,N-1.
!
!       Y = Array of length N containing data values asso-
!           ciated with the abscissae.  H(X(I)) = Y(I) for
!           I = 1,...,N.  If NCD = 1 and PER = TRUE, Y(N)
!           is set to Y(1).
!
!       NCD = Number of continuous derivatives at the knots.
!             NCD = 1 or NCD = 2.  If NCD = 1, the YP values
!             are computed by local monotonicity-constrained
!             quadratic fits.  Otherwise, a linear system is
!             solved for the derivative values which result
!             in second derivative continuity.  Unless
!             UNIFRM = TRUE, this requires iterating on
!             calls to YPC2 or YPC2P and calls to SIGS, and
!             generally results in more nonzero tension
!             factors (hence more expensive evaluation).
!
!       IENDC = End condition indicator for NCD = 2 and PER
!               = FALSE (or dummy parameter otherwise):
!               IENDC = 0 if YP(1) and YP(N) are to be com-
!                         puted by monotonicity-constrained
!                         parabolic fits to the first three
!                         and last three points, respective-
!                         ly.  This is identical to the
!                         values computed by YPC1.
!               IENDC = 1 if the first derivatives of H at
!                         X(1) and X(N) are user-specified
!                         in YP(1) and YP(N), respectively.
!               IENDC = 2 if the second derivatives of H at
!                         X(1) and X(N) are user-specified
!                         in YP(1) and YP(N), respectively.
!               IENDC = 3 if the end conditions are to be
!                         computed by Subroutine ENDSLP and
!                         vary with SIGMA(1) and SIGMA(N-1).
!
!       PER = Logical variable with value TRUE if and only
!             H(x) is to be a periodic function with period
!             X(N)-X(1).  It is assumed without a test that
!             Y(N) = Y(1) in this case.  On output, YP(N) =
!             YP(1).  If H(x) is one of the components of a
!             parametric curve, this option may be used to
!             obtained a closed curve.
!
!       UNIFRM = Logical variable with value TRUE if and
!                only if constant (uniform) tension is to be
!                used.  The tension factor must be input in
!                SIGMA(1) in this case and must be in the
!                range 0 to 85.  If SIGMA(1) = 0, H(x) is
!                piecewise cubic (a cubic spline if NCD =
!                2), and as SIGMA increases, H(x) approaches
!                the piecewise linear interpolant.  If
!                UNIFRM = FALSE, tension factors are chosen
!                (by SIGS) to preserve local monotonicity
!                and convexity of the data.  This often
!                improves the appearance of the curve over
!                the piecewise cubic fit.
!
!       LWK = Length of work space WK:  no work space is
!             needed if NCD = 1; at least N-1 locations
!             are required if NCD = 2; another N-1 locations
!             are required if PER = TRUE; and an additional
!             N-1 locations are required for the convergence
!             test if SIGS is called (UNIFRM = FALSE):
!
!             LWK GE 0    if NCD=1
!             LWK GE N-1  if NCD=2, PER=FALSE, UNIFRM=TRUE
!             LWK GE 2N-2 if NCD=2, PER=TRUE,  UNIFRM=TRUE
!             LWK GE 2N-2 if NCD=2, PER=FALSE, UNIFRM=FALSE
!             LWK GE 3N-3 if NCD=2, PER=TRUE,  UNIFRM=FALSE
!
!       t   = optional argument,tension parameter , between 0.0 and 85.
!   The above parameters, except possibly Y(N), are not
! altered by this routine.
!
!       WK = Array of length at least LWK to be used as
!            temporary work space.
!
!       YP = Array of length .GE. N containing end condition
!            values in positions 1 and N if NCD = 2 and
!            IENDC = 1 or IENDC = 2.
!
!       SIGMA = Array of length .GE. N-1 containing a ten-
!               sion factor (0 to 85) in the first position
!               if UNIFRM = TRUE.
!
! On output:
!
!       WK = Array containing convergence parameters in the
!            first two locations if IER > 0 (NCD = 2 and
!            UNIFRM = FALSE):
!            WK(1) = Maximum relative change in a component
!                    of YP on the last iteration.
!            WK(2) = Maximum relative change in a component
!                    of SIGMA on the last iteration.
!
!       YP = Array containing derivatives of H at the
!            abscissae.  YP is not altered if -4 < IER < 0,
!            and YP is only partially defined if IER = -4.
!
!       SIGMA = Array containing tension factors.  SIGMA(I)
!               is associated with interval (X(I),X(I+1))
!               for I = 1,...,N-1.  SIGMA is not altered if
!               -4 < IER < 0 (unless IENDC is invalid), and
!               SIGMA is constant (not optimal) if IER = -4
!               or IENDC (if used) is invalid.
!
!       IER = Error indicator or iteration count:
!             IER = IC .GE. 0 if no errors were encountered
!                      and IC calls to SIGS and IC+1 calls
!                      to YPC1, YPC1P, YPC2 or YPC2P were
!                      employed.  (IC = 0 if NCD = 1).
!             IER = -1 if N, NCD, or IENDC is outside its
!                      valid range.
!             IER = -2 if LWK is too small.
!             IER = -3 if UNIFRM = TRUE and SIGMA(1) is out-
!                      side its valid range.
!             IER = -4 if the abscissae X are not strictly
!                      increasing.
!
! Modules required by TSPSI:  ENDSLP, SIGS, SNHCSH, STORE,
!                               YPCOEF, YPC1, YPC1P, YPC2,
!                               YPC2P
!
! Intrinsic functions called by TSPSI:  ABS, MAX
!


! ------------------------------------------------------------------ HSJ
!
      INTEGER(I4b),INTENT(IN)             :: n,ncd,iendc
      REAL(DP),INTENT(inout),DIMENSION(:) :: x,y,deriv_bc
      REAL(DP),INTENT(in),OPTIONAL        :: t

      INTEGER(I4b) lwk,ier

      LOGICAL per,unifrm
 


       IF(ALLOCATED(wk))DEALLOCATE(wk)
       IF(ALLOCATED(yp))DEALLOCATE(yp)
       IF(ALLOCATED(sigma))DEALLOCATE(sigma)
       lwk = 1
       IF( ncd == 2 .AND. .NOT. per  .AND. UNIFRM )lwk = 2*n-2
       IF( ncd == 2 .AND. per   .AND. UNIFRM )lwk = 2*n-2
       IF( ncd == 2 .AND. .NOT. per  .AND.  .NOT. UNIFRM )lwk = 2*n-2
       IF( ncd == 2 .AND. per  .AND. .NOT. UNIFRM)lwk = 3*n-3

       ALLOCATE(wk(lwk))
       ALLOCATE(yp(n))
       ALLOCATE(sigma(n))
       sigma(:)= 0.0_DP
       IF(unifrm .AND. PRESENT(t))sigma(1) = t
       IF(iendc == 1)THEN
          yp(1) = deriv_bc(1)
          yp(n) = deriv_bc(2)
       ENDIF
       CALL TSPSI (N,X,Y,NCD,IENDC,PER,UNIFRM,LWK, WK, &
                        YP,SIGMA, IER)
       IF(ier == -1_I4B)THEN

          WRITE(ncrt,1)N, NCD,IENDC
          WRITE(nlog,1)
1         FORMAT(2x,'N, NCD, or IENDC out of range :',/,&
                 2x,'N, NCD, or IENDC =',3(2x,i3))      

       ELSEIF(ier == -2_I4B)THEN
          WRITE(ncrt,2)lwk
          WRITE(nlog,2)lwk
2         FORMAT(2x,'dimension lwk too small,lwk,n =',2(2X,I3))
       ELSEIF(ier == -3_I4B)THEN
          WRITE(ncrt,3)sigma(1)
          WRITE(nlog,3)sigma(1)
3         FORMAT(2x,"Selected uniform tension which is out of range",/,&
                 2x,'tension = ',f10.3)
       ELSEIF(ier == -4_I4B)THEN
          WRITE(ncrt,4)x(1:n)
          WRITE(nlog,4)x(1:n)
4         FORMAT(2x,'The abscissa values are not strictly increasing',/,&
                 (5(2x,1pe14.8)))
       ENDIF

       IF(ier .LT. 0_I4B)STOP 'sub t716_tspsi problem'

      RETURN
      END SUBROUTINE t716_TSPSI 

      
      SUBROUTINE t716_TSVAL1 (N,X,Y,IFLAG,NE,TE,V)
! -----------------------------------------------------------------
!   this is a  wrapper for TSVAL1.
!   TSVAL1  evaluates a Hermite interpolatory ten-
!   sion spline H or its first or second derivative at a set
!   of points TE. (Can be used with a number of routiens, including
!   TSPSI).

!   Note that a large tension factor in SIGMA may cause
! underflow.  The result is assumed to be zero.  If not the
! default, this may be specified by either a compiler option
! or operating system option.
!
! On input:
!
!       N = Number of data points.  N .GE. 2.
!
!       X = Array of length N containing the abscissae.
!           These must be in strictly increasing order:
!           X(I) < X(I+1) for I = 1,...,N-1.
!
!       Y = Array of length N containing data values or
!           function values returned by Subroutine SMCRV.
!           Y(I) = H(X(I)) for I = 1,...,N.
 
!
!       YP = Array of length N containing first deriva-
!            tives.  YP(I) = HP(X(I)) for I = 1,...,N, where
!            HP denotes the derivative of H.
!
!       SIGMA = Array of length N-1 containing tension fac-
!               tors whose absolute values determine the
!               balance between cubic and linear in each
!               interval.  SIGMA(I) is associated with int-
!               erval (I,I+1) for I = 1,...,N-1.
!
!       IFLAG = Output option indicator:
!               IFLAG = 0 if values of H are to be computed.
!               IFLAG = 1 if first derivative values are to
!                         be computed.
!               IFLAG = 2 if second derivative values are to
!                         be computed.
!
!       NE = Number of evaluation points.  NE > 0.
!
!       TE = Array of length NE containing the evaluation
!            points.  The sequence should be strictly in-
!            creasing for maximum efficiency.  Extrapolation
!            is performed if a point is not in the interval
!            [X(1),X(N)].
!
! The above parameters are not altered by this routine.
!
!       V = Array of length at least NE.
!
! On output:
!
!       V = Array of function, first derivative, or second
!           derivative values at the evaluation points un-
!           less IER < 0.  If IER = -1, V is not altered.
!           If IER = -2, V may be only partially defined.
!
!       IER = Error indicator:
!             IER = 0  if no errors were encountered and
!                      no extrapolation occurred.
!             IER > 0  if no errors were encountered but
!                      extrapolation was required at IER
!                      points.
!             IER = -1 if N < 2, IFLAG < 0, IFLAG > 2, or
!                      NE < 1.
!             IER = -2 if the abscissae are not in strictly
!                      increasing order.  (This error will
!                      not necessarily be detected.)
!
! Modules required by TSVAL1:  HPPVAL, HPVAL, HVAL, INTRVL,
!                                SNHCSH

!--------------------------------------------------------------------




      INTEGER(I4B)  N, IFLAG, NE, IER

      REAL(DP),INTENT(inout),DIMENSION(:)   :: x,y,te,v


      CALL TSVAL1 (N,X,Y,YP,SIGMA,IFLAG,NE,TE, V, &
                         IER)
       IF(ier == -1_I4B)THEN
          WRITE(ncrt,1)N,IFLAG,NE
          WRITE(nlog,1)
1         FORMAT(2x,'N, IFLAG, or NE  out of range :',/,&
                 2x,'N, IFLAG, or NE  =',3(2x,i3))      
       ELSEIF(ier == -2_I4B)THEN
          WRITE(ncrt,4)x(1:n)
          WRITE(nlog,4)x(1:n)
4         FORMAT(2x,'The abscissa values are not strictly increasing',/,&
                 (5(2x,1pe14.8)))
       ENDIF

       IF(ier .LT. 0_I4B) STOP 'sub t716_tsval1 problem'

      RETURN
      END SUBROUTINE t716_TSVAL1 






      SUBROUTINE t716_TSPSS(N,X,Y,PER,UNIFRM,SM,SMTOL,YS,smlevel,t)
!-------------------------------------------------------------------------
! --- this subroutine  is a wrapper for toms 716 tspline package
! --- routine TSPSS   (all toms tspline routines are in splinelib.a)
!
! This subroutine computes a set of parameter values which
! define a smoothing tension spline H(x).  The parameters
! consist of knot values YS and derivatives YP computed
! by Subroutine SMCRV, and tension factors SIGMA computed by
! Subroutine SIGS (unless UNIFRM = TRUE).  The Hermite
! interpolatory tension spline H(x) defined by the knot
! values and derivatives has two continuous derivatives and
! satisfies either natural or periodic end conditions.
!
!   The tension spline may be evaluated by Subroutine TSVAL1
! or Functions HVAL (values), HPVAL (first derivatives),
! HPPVAL (second derivatives), and TSINTL (integrals).
!
! On input:
!
!       N = Number of data points.  N .GE. 2 and N .GE. 3 if
!           PER = TRUE.
!
!       X = Array of length N containing a strictly in-
!           creasing sequence of abscissae:  X(I) < X(I+1)
!           for I = 1,...,N-1.
!
!       Y = Array of length N containing data values asso-
!           ciated with the abscissae.  If PER = TRUE, it is
!           assumed that Y(N) = Y(1).
!
!       PER = Logical variable with value TRUE if and only
!             H(x) is to be a periodic function with period
!             X(N)-X(1).  It is assumed without a test that
!             Y(N) = Y(1) in this case.  On output, YP(N) =
!             YP(1) and, more generally, the values and
!             first two derivatives of H at X(1) agree with
!             those at X(N).  If H(x) is one of the compo-
!             nents of a parametric curve, this option may
!             be used to obtained a closed curve.  If PER =
!             FALSE, H satisfies natural end conditions:
!             zero second derivatives at X(1) and X(N).
!
!       UNIFRM = Logical variable with value TRUE if and
!                only if constant (uniform) tension is to be
!                used.  The tension factor must be input in
!                SIGMA(1) in this case and must be in the
!                range 0 to 85.  If SIGMA(1) = 0, H(x) is
!                a cubic spline, and as SIGMA increases,
!                H(x) approaches piecewise linear.  If
!                UNIFRM = FALSE, tension factors are chosen
!                (by SIGS) to preserve local monotonicity
!                and convexity of the data.  This may re-
!                sult in a better fit than the case of
!                uniform tension, but requires an iteration
!                on calls to SMCRV and SIGS.
!
!       W = Array of length N containing positive weights
!           associated with the data values.  The recommend-
!           ed value of W(I) is 1/DY**2, where DY is the
!           standard deviation associated with Y(I).  If
!           nothing is known about the errors in Y, a con-
!           stant (estimated value) should be used for DY.
!           If PER = TRUE, it is assumed that W(N) = W(1).
!
!       SM = Positive parameter specifying an upper bound on
!            Q2(YS), where Q2(YS) is the weighted sum of
!            squares of deviations from the data (differ-
!            ences between YS and Y).  H(x) is linear (and
!            Q2 is minimized) if SM is sufficiently large
!            that the constraint is not active.  It is
!            recommended that SM satisfy N-SQRT(2N) .LE. SM
!            .LE. N+SQRT(2N) and SM = N is reasonable if
!            W(I) = 1/DY**2.
!
!       SMTOL = Parameter in the range (0,1) specifying the
!               relative error allowed in satisfying the
!               constraint:  the constraint is assumed to
!               be satisfied if SM*(1-SMTOL) .LE. Q2 .LE.
!               SM*(1+SMTOL).  A reasonable value for SMTOL
!               is SQRT(2/N).
!
!       LWK = Length of work space WK:
!             LWK .GE. 6N   if PER=FALSE  and  UNIFRM=TRUE
!             LWK .GE. 7N   if PER=FALSE  and  UNIFRM=FALSE
!             LWK .GE. 10N  if PER=TRUE   and  UNIFRM=TRUE
!             LWK .GE. 11N  if PER=TRUE   and  UNIFRM=FALSE
!
!      smlevel weight factor. Here it is a single vlaue
!              but more general use is possible .
!              Set smlevel to yield desired smoothness >=0.0
!
! The above parameters are not altered by this routine.
!
!       WK = Array of length at least LWK to be used as
!            temporary work space.
!
!       SIGMA = Array of length .GE. N-1 containing a ten-
!               sion factor (0 to 85) in the first position
!               if UNIFRM = TRUE.
!
!       YS = Array of length .GE. N.
!
!       YP = Array of length .GE. N.
!
! On output:
!
!       WK = Array containing convergence parameters in the
!            first two locations if NIT > 0:
!            WK(1) = Maximum relative change in a component
!                    of YS on the last iteration.
!            WK(2) = Maximum relative change in a component
!                    of SIGMA on the last iteration.
!
!       SIGMA = Array containing tension factors.  SIGMA(I)
!               is associated with interval (X(I),X(I+1))
!               for I = 1,...,N-1.  SIGMA is not altered if
!               N is invalid or -4 < IER < -1, and SIGMA is
!               constant if IER = -1 (and N is valid) or
!               IER = -4.
!
!
!       YS = Array of length N containing values of H at the
!            abscissae.  YS(N) = YS(1) if PER = TRUE.  YS is
!            not altered if IER < 0.
!
!       YP = Array of length N containing first derivative
!            values of H at the abscissae.  YP(N) = YP(1)
!            if PER = TRUE.  YP is not altered if IER < 0.
!
!       NIT = Number of iterations (calls to SIGS).  NIT = 0
!             if IER < 0 or UNIFRM = TRUE.  If NIT > 0,
!             NIT+1 calls to SMCRV were employed.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered and the
!                     constraint is active:  Q2(YS) is ap-
!                     proximately equal to SM.
!             IER = 1 if no errors were encountered but the
!                     constraint is not active:  YS and YP
!                     are the values and derivatives of the
!                     linear function (constant function if
!                     PERIOD = TRUE) which minimizes Q2, and
!                     Q1 = 0 (refer to SMCRV).
!
!             IER = -1 if N, W, SM, or SMTOL is outside its
!                      valid range.
!             IER = -2 if LWK is too small.
!             IER = -3 if UNIFRM = TRUE and SIGMA(1) is out-
!                      side its valid range.
!             IER = -4 if the abscissae X are not strictly
!
!-------------------------------------------------------------------------

      INTEGER(I4b),INTENT(IN)             :: n
      REAL(DP),INTENT(inout),DIMENSION(:) :: x,y,ys
      REAL(DP),INTENT(in),OPTIONAL        :: t
      INTEGER      LWK, NIT, IER
      LOGICAL      PER, UNIFRM
      REAL(DP)     SM, SMTOL,smlevel

 
       IF(per .eqv. .FALSE. .AND. unifrm .eqv. .TRUE. )lwk =6*n
       IF(per .eqv. .FALSE. .AND. unifrm .eqv. .FALSE.)lwk =7*n
       IF(per .eqv. .TRUE.  .AND. unifrm .eqv. .TRUE. )lwk =10*n
       IF(per .eqv. .TRUE.  .AND. unifrm .eqv. .FALSE.)lwk =11*n

       IF(ALLOCATED(wk))DEALLOCATE(wk)
       IF(ALLOCATED(w))DEALLOCATE(w)
       IF(ALLOCATED(yp))DEALLOCATE(yp)
       IF(ALLOCATED(sigma))DEALLOCATE(sigma)
       ALLOCATE(wk(lwk))
       ALLOCATE(w(n))
       ALLOCATE(yp(n))
       ALLOCATE(sigma(n))


       sigma(:)= 0.0_DP 
       IF(unifrm .AND. PRESENT(t))sigma(1) = t
       w(:)     = smlevel
       CALL TSPSS (N,X,Y,PER,UNIFRM,W,SM,SMTOL,LWK, WK, &
                        SIGMA,YS,YP, NIT,IER) 

       IF(ier == -1_I4B)THEN

          WRITE(ncrt,1)N,W,SM,SMTOL
          WRITE(nlog,1)
1         FORMAT(2x,'N,W,SM,or SMTOL  out of range :',/,&
                 2x,'N,W,SM ,SMTOL = ',i3,/,5(2x,1pe12.4))      
       ELSEIF(ier == 1_I4B)THEN
          WRITE(ncrt,5)
          WRITE(nlog,5)
5         FORMAT(2x,'t716_tspss constraints not active in smoothing =',2(2X,I3))
       ELSEIF(ier == -2_I4B)THEN
          WRITE(ncrt,2)lwk
          WRITE(nlog,2)lwk
2         FORMAT(2x,'dimension lwk too small,lwk,n =',2(2X,I3))
       ELSEIF(ier == -3_I4B)THEN
          WRITE(ncrt,3)sigma(1)
          WRITE(nlog,3)sigma(1)
3         FORMAT(2x,"Selected uniform tension which is out of range",/,&
                 2x,'tension = ',f10.3)
       ELSEIF(ier == -4_I4B)THEN
          WRITE(ncrt,4)x(1:n)
          WRITE(nlog,4)x(1:n)
4         FORMAT(2x,'The abscissa values are not strictly increasing',/,&
                 (5(2x,1pe14.8)))
       ENDIF

       IF(ier .LT. 0_I4B)STOP 'sub t716_tspss problem '

      RETURN
      END SUBROUTINE t716_TSPSS



      SUBROUTINE smooth_values(ydata,xdata,ndata,smthdata,dsmthdatadx,d2smthdatadx2,  &
                               weight,no_weight,sdevfctr,ibc,derivl,derivr)
!------------------------------------------------------------------------------------------
!   INPUT
!   ydata(i)          noisy data , i =1,ndata
!   xdata(i)          abscissa associated with ydata
!   weights(i)        the weighting of each data point. 
!                     small weights lead to an interpolatory spline
!   no_weight = 0     means weights are given,eg. weight is input
!   no_weight = 1     means weights are internally 
!                     calculated in this routine (using the variance of the data)
!                     weight(i) is then defined in this routine  and becomes an
!                     output
!   sdevfctr          A multiplier for the weights,used if no_weight = 1.
!                     hence the actual weight becomes some (possibly fractional)
!                     multiple of the standard deviation
!   ibc               boundary condition specifier:
!                      ibc =1 input first deriv at ends,derivl at xdata(1)
!                                                       derivr at x(ndata)
!                          =2 means second deriv = 0 at ends, (no input required)
!                          =3 means periodic at ends,         (no input required)

!   OUPUT
!   smthdata(i)        the smoothed data
!   dsmthdatadx(I)     first deriv
!   d2smthdatadx2(I)   second derivative 
!   weight(i)          see above
! -------------------------------------------------------------------HSJ----------
       IMPLICIT NONE

       REAL(DP),DIMENSION(:) :: ydata,xdata,smthdata,weight,dsmthdatadx(ndata),           &
                                d2smthdatadx2(ndata)    
       !temporary arrays:
       REAL(DP),DIMENSION(:) :: a(ndata),b(ndata),c(ndata),d(ndata),e(ndata),g(ndata),    &  
                                zm(ndata),aa(ndata),dd(ndata),                            &
                                u(ndata+2),v(ndata+2),w(ndata+2),p(ndata+2),q(ndata+2),   &
                                r(ndata+2),s(ndata+2),t(ndata+2) 


       REAL(DP)              :: derivl,derivr,ave,sdev,adev,var,skew,curt,                &
                                deriv1,deriv2,sdevfctr
       INTEGER(I4B)          :: ndata, no_weight,ibc,ind,i


       a(:)=0.0_DP; b(:)=0.0_DP; c(:)=0.0_DP ;d(:)=0.0_DP ;e (:)=0.0_DP
       g(:)=0.0_DP ; zm(:) =0.0_DP ; aa(:) =0.0_DP ;  dd(:) =0.0_DP  
       u(:) =0.0_DP ;  v(:) =0.0_DP ;  w(:) =0.0_DP 
       p(:) =0.0_DP ;  q(:) =0.0_DP ;  r(:) =0.0_DP 
       s(:) =0.0_DP ;  t(:) =0.0_DP
       
       IF(no_weight == 1)THEN
          CALL moment(ydata,ndata,ave,adev,sdev,var,skew,curt)
          weight(:) = sdevfctr*sdev
       ENDIF


       ind = 0 ! determine spline coefficients 
       CALL Smspline(ndata,ibc,ind,xdata,ydata,a,b,c,d,e,g,weight,        &
                      derivl,derivl,u,v,w,p,q,r,s,t,aa,dd,zm,xdata(1),    &
                      smthdata(1),deriv1,deriv2)







       ind = 1  ! evaluate spline and first and seconde derivative point xx
       DO i = 1,ndata
            CALL Smspline(ndata,ibc,ind,xdata,ydata,a,b,c,d,e,g,weight,    &
                          derivr,derivl,u,v,w,p,q,r,s,t,aa,dd,zm,xdata(i),  &
                          smthdata(i),deriv1,deriv2)
       ENDDO




       RETURN

       END SUBROUTINE smooth_values







        SUBROUTINE Progon5(n,a,b,c,d,e,g,u,v,w, &
                            p,q,r,s,t,aa,dd,z)

!
!***********************************************************
!    Program Progon5 for solving linear system of the      *
!    form                                                  *
! Ú                                          ¿ Ú  ¿ Ú  ¿   *
! ³a1   b1   c1    0   ..........0   e1   d1 ³ ³z1³ ³g1³   *
! ³d2   a2   b2   c2    0 ........    0   e2 ³ ³z2³ ³g2³   *
! ³e3   d3   a3   b3   c3    0 ...    0    0 ³ ³z3³ ³g3³   *
! ³O    e4   d4   a4   b4   c4   .         0 ³ ³z4³=³g4³   *
! ³..........................................³ ³..³ ³  ³   *
! ³                                          ³ ³  ³ ³  ³   *
! ³bn   cn   0     0  .....  0   en  dn   an ³ ³zn³ ³gn³   *
! À                                          Ù À  Ù À  Ù   *
!  with pentadiagonal matrix by sweep method               *
!  n - number of equations                                 *
!  a,b,c,d,e - arrays (of size n) for elements of          *
!              the matrix's diagonal                       *
!  g - array (of size n) for elements of the right parts   *
!                                                          *
!  u,v,w,p,q,                                              *
!  r,s,t,aa,dd - work arrays (of size n+2)                 *
!                                                          *
!        Output:                                           *
!  z - array (of size n) containing solution of the        *
!      system                                              *
!***********************************************************
!
       REAL(DP)  a(:),b(:),c(:),d(:),e(:),g(:),u(:),aa(:),r(:)
       REAL(DP)  v(:),w(:),p(:),q(:),t(:),s(:),z(:),dd(:)
       REAL(DP)  a11,a12,a21,a22,b1,b2
       INTEGER(I4B) n,i
!
              v(1) = 0.
              v(2) = 0.
              w(1) = 0.
              w(2) = 0.
              u(1) = 0.
              u(2) = 0.
!
        DO i = 1,n
         dd(i)  =  d(i) + e(i)*v(i)
         aa(i)  =  a(i) + dd(i)*v(i+1) +e(i)*w(i)
         u(i+2) = (g(i) - dd(i)*u(i+1) - e(i)*u(i))/aa(i)
         v(i+2) = -(b(i)+dd(i)*w(i+1))/aa(i)
         w(i+2) = -c(i)/aa(i)
      ENDDO
!
            p(1) = 0.
            q(2) = 0.
            p(2) = 1.
            q(1) = 1.
!
           DO i = 1,n
         p(i+2) = -(dd(i)*p(i+1) + e(i)*p(i))/aa(i)
         q(i+2) = -(dd(i)*q(i+1) + e(i)*q(i))/aa(i)
      ENDDO
!
                t(n-1) = 0.
                s(n-1) = 0.
                t(n)   = 0.
                r(n)   = 0.
             r(n-1) = 1.
             s(n)   = 1.
!
       DO i = n-2,1,-1
          t(i) = v(i+2)*t(i+1) +w(i+2)*t(i+2) + u(i+2)
          s(i) = v(i+2)*s(i+1) +w(i+2)*s(i+2) + p(i+2)
          r(i) = v(i+2)*r(i+1) +w(i+2)*r(i+2) + q(i+2)
       ENDDO
!
       a11 =  1. - q(n+1) - w(n+1)*r(1)
       a12 = -(p(n+1) + v(n+1) + w(n+1)*s(1))
       a21 = -(v(n+2)*r(1) + w(n+2)*r(2) + q(n+2))
       a22 =  1. - p(n+2) - v(n+2)*s(1) - w(n+2)*s(2)
       b1  = w(n+1)*t(1) + u(n+1)
       b2  = v(n+2)*t(1) + w(n+2)*t(2) + u(n+2)
       z(n-1) = ( b1*a22 - b2*a12)/(a11*a22-a12*a21)
       z(n)   = (-b1*a21 + b2*a11)/(a11*a22-a12*a21)
!
       DO i = 1,n-2
          z(i) = t(i) + s(i)*z(n) + r(i)*z(n-1)
       ENDDO
!
       RETURN
       END SUBROUTINE progon5




       SUBROUTINE Smspline(n,ib,ind,x,y,a,b,c,d,e,g,rho, &
          ax,bx,u,v,w,p,q,r,s,t,aa,dd,zm,xx,sp,dsp,d2sp)
!
!************************************************************
!       Program Smspline constructs the smoothing spline    *
!                for given array of points                  *
!                                                           *
!  n   - number of points in array                          *
!  x   - array (of size n) containing x-coordinates         *
!  y   - array (of size n) containing y-coordinates         *
!  rho - array (of size n) containing weight coefficients   *
!                                                           *
!  ind - mode type:                                         *
!        ind = 0 - program calculates parameters of spline  *
!        ind = 1 - the spline parameters are known          *
!  ib  - type of boundary conditions                        *
!         ib = i (i=1-3) - conditions of ith type           *
!                                                           *
!  a,b,c,d,e,g,zm,aa,dd - work arrays (of size n)           *
!  u,v,w,p                                                  *
!  q,r,s,t - work arrays (of size n+2)                      *
!  ax, bx - parameter values in boundary conditions         *
!           of the first type                               *
!                                                           *
!         Output:                                           *
!  sp   - value of spline at point xx                       *
!  dsp  - value of its first derivative at point xx         *
!  d2sp - value of its second derivative at point xx        *
!                                                           *
!************************************************************
!
       REAL(DP)  a(:),b(:),c(:),d(:),e(:), &
            g(:),rho(:),aa(:),r(:)
       REAL(DP)  u(:),v(:),w(:),p(:),q(:),t(:),s(:),zm(:),dd(:)
       REAL(DP)  x(:), y(:),ax,bx,xx,sp,dsp,d2sp,h1,h2,h3,h4,di,h,tt
       INTEGER(I4B) ib,ind,n,nn,i,ni
!
        nn = n
        IF (ind.EQ.0) THEN
        SELECT CASE (ib)
!
        CASE (1)
         h1     = x(2)-x(1)
         a(1)   = h1/3.+(rho(1)+rho(2))/h1**2
         h2     = x(3) - x(2)
         b(1)   = h1/6.-(1./h1 +1./h2)*rho(2)/h1 &
                  - rho(1)/h1**2
         c(1)   = rho(2)/(h1*h2)
         g(1)   = (y(2) - y(1))/h1 - ax
         e(3)   = c(1)
         h1     = x(n) - x(n-1)
         h2     = x(n-1) - x(n-2)
         a(n)   = h1/3. +(rho(n-1)+rho(n))/h1**2
         b(n-1) = h1/6.-(1./h1 +1./h2)*rho(n-1)/h1 &
                 -rho(n)/h1**2
         c(n-2) = rho(n-1)/(h1*h2)
         e(1) = 0.
         e(2) = 0.
         d(1) = 0.
         d(2) = b(1)
         c(n-1) = 0.
         c(n) = 0.
         b(n) = 0.
         d(n) = b(n-1)
         g(n) = bx - (y(n) - y(n-1))/h1
         e(n) = c(n-2)
!
        CASE (2)
         a(1)   = 1.
         a(n)   = 1.
         b(1)   = 0.
         b(n-1) = 0.
         c(1)   = 0.
         c(n-2) = 0.
         g(1)   = 0.
         g(n)   = 0.
         e(1)   = 0.
         e(2)   = 0.
         e(3)   = 0.
         d(1)   = 0.
         b(n)   = 0.
         c(n)   = 0.
         c(n-1) = 0.
!
        CASE (3)
         h1 = x(2) - x(1)
         h2 = x(n) - x(n-1)
         h3 = x(3) - x(2)
         a(1) = (h1 + h2)/3. &
               + rho(n-1)/h2**2 + rho(2)/h1**2 + &
                (1./h1 + 1./h2)**2*rho(1)
         b(1) = h1/6. -((1./h2 + 1./h1)*rho(1)+ &
                        (1./h2 + 1./h3)*rho(2))/h2
         c(1) = rho(1)/(h1*h3)
         g(1) = (y(2) - y(1))/h1 - (y(1) - y(n-1))/h2
         d(2) = b(1)
         e(3) = c(1)
         h4   = x(n-1) - x(n-2)
         c(n-2) = rho(n-1)/(h2*h4)
         c(n-1) = rho(n)/(h2*h1)
         b(n-1) = h2/3. - ((1./h4 + 1./h2)*rho(n-1) + &
                           (1./h2 + 1./h1)*rho(n))/h2
         d(1) = b(n-1)
         e(1) = c(n-2)
         e(2) = c(n-1)
         g(n-1) = (y(n) - y(n-1))/h2 - (y(n-1) - y(n-2))/h4
         nn = n - 1
       END SELECT
!
      DO i = 2, n-1
         h1   = x(i)   - x(i-1)
       h2   = x(i+1) - x(i)
         a(i) = (h1+h2)/3. &
               + rho(i-1)/h1**2 + rho(i+1)/h2**2 + &
                 (1./h1 + 1./h2)**2*rho(i)
       g(i) = (y(i+1)-y(i))/h2 - (y(i)-y(i-1))/h1
            IF (i.LE.(n-2)) THEN
                h3 = (x(i+2)-x(i+1))
              b(i) = h2/6. - ((1./h1 + 1./h2)*rho(i) &
                               +(1./h2 + 1./h3)*rho(i+1))/h2
              d(i+1) = b(i)
            ENDIF
         IF (i.LE.(n-3)) THEN
               c(i) = rho(i+1)/(h2*h3)
               e(i+2) = c(i)
           ENDIF
      ENDDO
!
        CALL Progon5 &
             (nn,a,b,c,d,e,g,u,v,w,p,q,r,s,t,aa,dd,zm)
        aa(1) = y(1) - rho(1)*(zm(2)-zm(1))/(x(2)-x(1))
      aa(n) = y(n) + rho(n)*(zm(n)-zm(n-1))/(x(n)-x(n-1))

        IF (ib.EQ.3) THEN
            zm(n) = zm(1)
            di = (zm(2)-zm(1))/(x(2)-x(1))
            di = di - (zm(n)-zm(n-1))/(x(n)-x(n-1))
            aa(1) = y(1) - rho(1)*di
            aa(n) = aa(1)
        ENDIF
!
      DO i = 2, n-1
            di = (zm(i+1)-zm(i))/(x(i+1)-x(i))
            di = di - (zm(i)-zm(i-1))/(x(i)-x(i-1))
            aa(i) = y(i) - rho(i)*di
      ENDDO
!
      ELSE
!
        DO  i = 2,n
            ni= i
            IF(x(i).GT.xx) EXIT
        ENDDO
!
        i  = ni - 1
        h  = x(i+1) - x(i)
        tt = (xx - x(i))/h
        sp = aa(i)*(1.-tt) + aa(i+1)*tt &
            - h**2/6.*tt*(1.-tt)* &
             ((2.-tt)*zm(i) + (1.+tt)*zm(i+1))
        dsp= (aa(i+1) - aa(i))/h - &
              h/6.*((2. - 6.*tt + 3.*tt**2)*zm(i) &
            +(1. - 3.*tt**2)*zm(i+1))
        d2sp = (1. - tt)*zm(i) + tt*zm(i+1)
!
      ENDIF
      RETURN

      END       SUBROUTINE Smspline



        
      SUBROUTINE moment(DATA,n,ave,adev,sdev,var,skew,curt)
! ORIGINALlY FROM NUMERICAL RECEPIES
      IMPLICIT NONE
      INTEGER(I4B) n,j
      REAL(DP)  adev,ave,curt,sdev,skew,var,DATA(n)
      REAL(DP)  p,s,ep

      IF(n .LE. 1) THEN
         PRINT *,'Moments not calculated'
         AVE = DATA(1)
         ADEV =DATA(1)
         SDEV =DATA(1)
         VAR =SDEV*SDEV
         SKEW =0
         CURT =0
         RETURN
      ENDIF
      s=0.
      DO 11 j=1,n
        s=s+DATA(j)
11    CONTINUE
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      DO 12 j=1,n
        s=DATA(j)-ave
        ep=ep+s
        adev=adev+ABS(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
12    CONTINUE
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=SQRT(var)

      IF(var.NE.0.)THEN
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3.
      ELSE
        PRINT *, 'no skew or kurtosis when zero variance in moment'
        SKEW =0.0
        CURT =0.0
      ENDIF

      RETURN

      END       SUBROUTINE moment



   SUBROUTINE Progon3 (a,b,c,d,x,n)
     !***********************************************************
     !        Program Progon3 for solving linear system         *
     !        of the form                                       *
     !                                                          *
     !        a(1)*x(1) + b(1)*x(2)      + c(1)*x(n) = d(1)     *
     !        ............................................      *
     !        c(i)*x(i-1) + a(i)*x(i)  + b(i)*x(i+1) = d(i)     *
     !        .............................................     *
     !        b(n)*x(1)     + a(n)*x(n-1) + b(n)*x(n)= d(n)     *
     !        with tridiagonal matrix by sweep method           *
     !    n - number of equations                               *
     !  a,b,                                                    *
     !  c,d - arrays (of size n) containing elements of the     *
     !        matrix' diagonals and of the right parts          *
     !  w,s,t,                                                  *
     !  u,v - work arrays (of size n+1)                         *
     !        Output:                                           *
     !    x - array (of size n) containing solution of the      *
     !        system                                            *
     !***********************************************************
     USE nrtype,                         ONLY : DP,I4B
     IMPLICIT NONE
     INTEGER(I4B) i,n,i1
     REAL(DP)  a(n),b(n),c(n),d(n),x(n)
     !        tEMPORARY ARRAYS:
     REAL(DP) u(n+1),v(n+1),w(n+1),s(n+1),t(n+1)
     REAL(DP) z

     u(1) = 0.
     v(1) = 0.
     w(1) = 1.
     !
     DO  i = 1,n
        i1 = i + 1
        z  = 1._DP/(a(i) + c(i)*v(i))
        v(i1) = -b(i)*z
        u(i1) = (-c(i)*u(i) + d(i))*z
        w(i1) = - c(i)*w(i)*z
     ENDDO
     !
     s(n) = 1._DP
     t(n) = 0._DP
     DO  i = n-1,1, -1
        s(i) = v(i+1)*s(i+1) + w(i+1)
        t(i) = v(i+1)*t(i+1) + u(i+1)
     ENDDO
     !
     x(n) = (d(n) - b(n)*t(1) - c(n)*t(n-1))/ &
          (a(n) + b(n)*s(1) + c(n)*s(n-1))

     DO i = 1, n-1
        x(i) = s(i)*x(n) + t(i)
     ENDDO
     !
     RETURN
   END SUBROUTINE Progon3


   SUBROUTINE Spline_intrp (n,x,y,z,get_scoef,ib,ax,bx,xx,sp,dsp,d2sp)
     !***********************************************************
     !    Program Spline constructs the interpolating cubic     *
     !    spline for the function given in tabulated form       *
     !                                                          *
     !     n - number of knots                                  *
     !     x - array (of size n) containing the interpolation   *
     !         knots                                            *
     !     y - array (of size n) containing the function values *
     !         at knots                                         *
     !     z - array (of size n) containing the spline          *
     !         parameters                                       *
     !     get_scoef :                                          *
     !         true - program calculates the spline parameters  *
     !         false  - the spline parameters are known         *
     !    xx - point where values of splines and its            *
     !         derivatives are calculated                       *
     !    ib - type of boundary conditions                      *
     !               ind = i (i=1-4) - conditions of ith type   *
     !    ax, bx - parameter values in boundary conditions      *
     !         of the first or second type                      *
     !         Output:                                          *
     !    sp   - the spline value at point xx                   *
     !    dsp  - value of its first derivative at point xx      *
     !    d2sp - value of its second derivative at point xx     *

     !   INTERNAL (temporary arrays):
     !   a,b,                                                   *
     !   c,d - work arrays (of size n)                          *
     !***********************************************************
     USE nrtype,                        ONLY : DP,I4B
     IMPLICIT NONE
     INTEGER(I4B) n,ip
     REAL(DP)  x(n),y(n)
     REAL(DP)  a(n),b(n),c(n),d(n),z(n)
     REAL(DP) sp,dsp,d2sp,h1,xx,ax,bx,h2,am,al,g0,cc,gn,h,tt,rp,aa,bb
     INTEGER(I4b) ind,ib,ne,ns,nf,nj,j
     LOGICAL get_scoef
     !
     ind = 1
     IF(get_scoef) ind =0
     ip = 1
     IF (ind.EQ.0) THEN
        a(1) =  2.
        ne   =  n
        ns   =  2
        nf   =  n - 1
        SELECT CASE (ib)
           !
        CASE (1)
           b(1 ) = 0.
           c(1) =  0.
           d(1) =  2.*ax
           a(n) =  2.
           b(n) =  0.
           c(n) =  0.
           d(n) =  2.*bx
           !
        CASE (2)
           b(1) =  1.
           c(1) =  0.
           h1   =  x(2)-x(1)
           d(1) =  3.*(y(2)-y(1))/h1 - 0.5*h1*ax
           a(n) =  2.
           b(n) =  0.
           c(n) =  1.
           h1   =  x(n)-x(n-1)
           d(n) =  3.*(y(n)-y(n-1))/h1 + 0.5*h1*bx
           !
        CASE (3)
           h1   = x(2) - x(1)
           h2   = x(n) - x(n-1)
           am   = h2/(h1 + h2)
           al   = 1. - am
           b(1) = am
           c(1) = al
           d(1) = 3.*(am*(y(2) - y(1))/h1 + &
                al*(y(1) - y(n-1))/h2)
           h1   = x(n-1) - x(n-2)
           h2   = x(n) - x(n-1)
           am   = h1/(h1 + h2)
           al   = 1. - am
           a(n-1) = 2.
           b(n-1) = am
           c(n-1) = al
           d(n-1) = 3.*(am*(y(n) - y(n-1))/h2 + &
                al*(y(n-1) - y(n-2))/h1)
           nf = n - 2
           ne = n - 1
           !
        CASE (4)
           h1   = x(2)-x(1)
           h2   = x(3)-x(2)
           g0   = h1/h2
           a(2) = 1. + g0
           b(2) = g0
           c(2) = 0.
           am   = h1/(h1+h2)
           al   = 1.- am
           cc   = am*(y(3)-y(2))/h2 + al*(y(2)-y(1))/h1
           d(2) = cc + 2.*g0*(y(3)-y(2))/h2
           h2   = x(n)-x(n-1)
           h1   = x(n-1)-x(n-2)
           gn   = h1/h2
           a(n-1) = 1. + gn
           b(n-1) = 0.
           c(n-1) = gn
           am = h1/(h1+h2)
           al = 1. - am
           cc = am*(y(n)-y(n-1))/h2 +al*(y(n-1)-y(n-2))/h1
           d(n-1) = cc + 2.*gn*(y(n-1)-y(n-2))/h1
           ns = 3
           nf = n - 2
           ne = n - 2
           ip = 2
        END SELECT
        !
        DO  j  = ns ,nf
           h1   = x(j + 1) - x(j)
           h2  = x(j) - x(j-1)
           am   = h2/(h2 + h1)
           al   = 1. - am
           c(j) = al
           a(j) = 2.
           b(j) = am
           d(j) = 3.*(am*(y(j+1) - y(j))/h1 + &
                al*(y(j) - y(j-1))/h2)
        ENDDO
        !
        CALL Progon3 (a(ip),b(ip),c(ip),d(ip), &
             z(ip),ne)
        !
        SELECT CASE (ib)
           !
        CASE (3)
           z(n) = z(1)
           !
        CASE (4)
           z(1) = g0**2*z(3)+(g0**2-1.)*z(2)+ &
                2.*((y(2)-y(1))/(x(2)-x(1)) &
                -g0**2*(y(3)-y(2))/(x(3)-x(2)))

           z(n) = gn**2*z(n-2)+(gn**2-1.)*z(n-1)+ &
                2.*((y(n)-y(n-1))/(x(n)-x(n-1)) &
                -gn**2*(y(n-1)-y(n-2)) &
                /(x(n-1)-x(n-2)))
        END SELECT


        RETURN
     ENDIF

     
     !
     DO  j = 2,n
        nj= j
        IF(x(j).GT.xx) go to 1
     ENDDO
     !
1    j  = nj - 1
     h  = x(j+1) - x(j)
     tt = (xx - x(j))/h
     rp = (y(j+1) - y(j))/h
     aa = -2.*rp + z(j) + z(j+1)
     bb = -aa +rp - z(j)
     sp = y(j) + (xx - x(j))*(z(j) + tt*(bb + tt*aa))
     dsp= z(j) + tt*(bb + aa*tt) + tt*(bb + 2.*aa*tt)
     d2sp = (2.*bb + 6.*aa*tt)/h
     RETURN
   END SUBROUTINE Spline_intrp




      SUBROUTINE spline_coef (n, x, y, b, c, d)
!-------------------------------------------------------------------
!  the coefficients b(i), c(i), and d(i), i = 1,2,...,n are computed
!  for a cubic interpolating spline

!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3

!    for  x(i) .le. x .le. x(i+1)

!  input..

!    n = the number of data points or knots (n .ge. 2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots

!  output..

!    b, c, d  = arrays of spline coefficients as defined above.

!  using  p  to denote differentiation,

!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)

!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!-------------------------------------------------------------------

      USE nrtype,                                       ONLY : Dp,I4B

                                  
!      IMPLICIT  INTEGER (i-n), REAL*8 (a-h, o-z)
      IMPLICIT NONE
!
      REAL(DP)  x(n), y(n), b(n), c(n), d(n)
      INTEGER(I4B) nm1, ib, i,n
      REAL(DP)  t

      nm1 = n - 1
      IF (n .LT. 2)  RETURN
      IF (n .LT. 3)  go to 50

!  set up tridiagonal system

!  b = diagonal, d = offdiagonal, ! = right hand side.

      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      DO i=2, nm1
         d(i)   = x(i+1) - x(i)
         b(i)   = 2.0 * (d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i)   = c(i+1) - c(i)
      END DO

!  end conditions.
!  third derivatives at x(1) and x(n) obtained from divided differences

      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.0
      c(n) = 0.0
      IF (n .EQ. 3)  go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))

!  forward elimination

   15 DO i=2, n
        t    = d(i-1)/b(i-1)
        b(i) = b(i) - t*d(i-1)
        c(i) = c(i) - t*c(i-1)
      END DO

!  back substitution

      c(n) = c(n)/b(n)
      DO ib=1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
      END DO

!  c(i) is now the sigma(i) of the text

!  compute polynomial coefficients

      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.0*c(n))
      DO i=1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.0*c(i)
      END DO
      c(n) = 3.0*c(n)
      d(n) = d(n-1)
      RETURN

   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.0
      d(1) = 0.0
      b(2) = b(1)
      c(2) = 0.0
      d(2) = 0.0
      RETURN

      END SUBROUTINE spline_coef



      FUNCTION seval (n, u, x, y, b, c, d)

! this subroutine evaluates the cubic spline function

!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3

!    where  x(i) < u < x(i+1), using Horner's rule

!  if  u .lt. x(1)  then  i = 1  is used.
!  if  u .ge. x(n)  then  i = n  is used.

!  ---------------------------HSJ Modification 1/30/98------------------
!  The above two rules are not acceptible in the context in which
!  this routine is used in the ONETWO code. (actually they are not
!  acceptable in any context in my opinion.
!  Interpolating outside the cubic spline fit is totally bogus.)
!  I have changed the rules to read
!   if u .gt. x(n) then use Y(n) as the most reasonable
!   approximation if you dont actually want to stop the code.
!   similarly if u .lt. x(1) then use y(1) as the returned value.
!  using the first or last point may not be much better in some circumstances
!  but at least we will get physically relevant values !!!!!!!!!!!HSJ
!  Note that the switch extend_seval must be set to 1 to get this
!  behavior. Otherwise an error stop is taken when trying to interpoalte
!  outside the tables.
! -----------------------------------------------------------------------

!  input..

!    n     = the number of data points
!    u     = the abscissa at which the spline is to be evaluated
!    x,y   = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline

!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.

!      USE nrtype,                                 ONLY : DP,I4B
!      USE param
!      USE io
!      IMPLICIT  INTEGER (i-n), REAL*8 (a-h, o-z)

       IMPLICIT NONE

      REAL(DP) seval
      INTEGER(I4B)  n,i, j, k,extend_seval
      REAL(DP)  u, x(n), y(n), b(n), c(n), d(n)
      REAL(DP)  dx

      DATA    i,extend_seval /1,0/


!  added  1/29/98  HSJ
 1    FORMAT('Subroutine Seval,a spline fit evaluater,',/,                &
            ' has detected that out of bounds interpolation',/,           &
            ' is occuring. The value at which the cubic spline',/,        &
            ' is to be evaluated is ',1pe18.8,/,                          &
            ' the interpolating table extends from ',1pe18.8,             &
            ' to ',1pe18.8,/,                                             &
            ' You can set extend_seval=1 in the first namelist',/,        &
            ' in inone if you want to ignore this problem')
      IF (u .LT. x(1)) THEN
         IF (extend_seval .EQ. 1) THEN
           seval = y(1)
           RETURN
         ELSE
           WRITE (ncrt, 1)  u, x(1), x(n)
!           CALL STOP ('subroutine SEVAL: bad interpolation #1', 267)
           lerrno = 224 + iomaxerr
           CALL terminate(lerrno,nlog)
         END IF
      ELSE IF (u .GT. x(n)) THEN
         IF (extend_seval .EQ. 1) THEN
           seval = y(n)
           RETURN
         ELSE
           WRITE (ncrt, 1)  u, x(1), x(n)
           lerrno = 224 + iomaxerr
           CALL terminate(lerrno,nlog)
!           CALL STOP ('subroutine SEVAL: bad interpolation #2', 268)
         END IF
      END IF


      IF (i .GE. n     )  i = 1
      IF (u .LT. x(i  ))  go to 10
      IF (u .LE. x(i+1))  go to 30

!  binary search

   10 i = 1
      j = n+1
   20 k = (i+j)/2
      IF (u .LT. x(k))  j = k
      IF (u .GE. x(k))  i = k
      IF (j .GT. i+1 )  go to 20

!  evaluate spline

   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      RETURN

      END FUNCTION seval


  END MODULE tension_spline
