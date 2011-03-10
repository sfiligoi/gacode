      FUNCTION U_ZEROIN(ax,bx,F,tol,iflag,message)
!***********************************************************************
!U_ZEROIN finds a zero of the function F(x) in the interval (ax,bx)
!References:
!  Brent, Algorithms for Minimization Without Derivatives, Prentice-Hall
!    Inc.(1973)
!  W.A.Houlberg 3/2000
!Input:
!  ax-left endpoint of initial interval
!  bx-right endpoint of initial interval
!  F-function subprogram which evaluates F(x) for any x in (ax,bx)
!  tol-desired length of the final interval of uncertainty
!Output:
!  U_ZEROIN-abcissa approximating a zero of F in the (ax,bx)
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!Comments:
!  U_ZEROIN returns a zero x in (ax,bx) to within a tolerance of
!    4*eps*abs(x)+tol, where eps is the relative machine precision
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      REAL           ax,                      bx,
     &               tol
      REAL           F
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           U_ZEROIN
!Declaration of local variables
      INTEGER        idone,                   nbis,
     &               nmes,                    nstep
      REAL           a,                       b,
     &               c,                       d,
     &               e,                       fa,
     &               fb,                      fc,
     &               p,                       q,
     &               r,                       s,
     &               toltst,                  xm
      REAL           z_precision
!Initialization
!  Physical and conversion constants
      z_precision=2.0e-7
!  Other
      nbis=0
      nstep=0
      a=ax
      b=bx
      fa=F(a,iflag,message)
      IF(iflag.eq.-1) THEN
        iflag=0
        message=' '
      ELSEIF(iflag.eq.1) THEN
!       Error evaluating function at a
        nmes=LEN(message)
        message='U_ZEROIN(1)/'//message(1:nmes-11)
        GOTO 1000
      ENDIF
      fb=F(b,iflag,message)
      IF(iflag.eq.-1) THEN
        iflag=0
        message=' '
      ELSEIF(iflag.eq.1) THEN
!       Error evaluating function at b
        nmes=LEN(message)
        message='U_ZEROIN(2)/'//message(1:nmes-11)
        IF(iflag.eq.1) GOTO 1000
      ENDIF
!Check whether F changes sign over the interval
      IF(ABS(SIGN(1.0,fa)+SIGN(1.0,fb)).gt.10.0*z_precision) THEN
        iflag=-1
        message='U_ZEROIN(3)/WARNING:no or multiple solutions'
        GOTO 1000
      ENDIF
!Initialize
      c=a
      fc=fa
      d=b-a
      e=d
      xm=b-a
      U_ZEROIN=b
      idone=0
!Begin step
      DO WHILE(idone.eq.0)
        IF(ABS(fc).lt.ABS(fb)) THEN
!         Reorder points a, b, and c   
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        ENDIF
!       Convergence test
        toltst=2.0*z_precision*ABS(b)+0.5*tol
        xm=0.5*(c-b)
        IF((ABS(xm).le.toltst).or.(fb.eq.0)) THEN
!         Converged within tolerance at x=b
          U_ZEROIN=b
          idone=1
        ELSE
!       Check if bisection is necessary
          IF((ABS(e).lt.toltst).or.(ABS(fa).le.ABS(fb))) THEN
!           Bisection
            d=xm
            e=d
            nbis=nbis+1
          ELSE
!           Check if quadratic interpolation is possible
            IF(a.ne.c) THEN
!             Inverse quadratic interpolation
              q=fa/fc
              r=fb/fc
              s=fb/fa
              p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
              q=(q-1.0)*(r-1.0)*(s-1.0)
            ELSE
!             Linear interpolation
              s=fb/fa
              p=2.0*xm*s
              q=1.0-s
            ENDIF
!           Adjust signs
            IF(p.gt.0.0) q=-q
            p=ABS(p)
!           Check if interpolation is acceptable
            IF(((2.0*p).ge.(3.0*xm*q-ABS(toltst*q)))
     &        .or.(p.ge.ABS(0.5*e*q))) THEN
!             Bisection
              d=xm
              e=d
              nbis=nbis+1
            ELSE
!             Interpolation
              e=d
              d=p/q
            ENDIF
          ENDIF
!Complete step
          a=b
          fa=fb
          IF(ABS(d).gt.toltst) b=b+d
          IF(ABS(d).le.toltst) b=b+SIGN(toltst,xm)
          fb=F(b,iflag,message)
          IF(iflag.gt.1) THEN
            nmes=LEN(message)
            message='U_ZEROIN(4)/'//message(1:nmes-11)
            GOTO 1000
          ENDIF
          nstep=nstep+1
          IF((fb*(fc/ABS(fc))).gt.0.0) THEN
            c=a
            fc=fa
            d=b-a
            e=d
          ENDIF
        ENDIF
      ENDDO
 1000 RETURN
      END
