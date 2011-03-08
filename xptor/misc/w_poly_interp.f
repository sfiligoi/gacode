      FUNCTION W_POLY_INTERP(x,y,npts,k_order,xin)
!***********************************************************************
!W_POLY_INTERP is a polynomial interpolating routine
!References:
!  Bevington, Data Reduction and Error Analysis for the Physical
!    Sciences, McGraw-Hill, Inc., 1969
!  W.A.Houlberg 1/1999
!Input: 
!  x(i)-array of abscissas 
!  y(i)-array of ordinates 
!  npts-number of data points in x and y arrays 
!  k_order-order of interpolation-(1<nterms<11) 
!         =2 linear 
!         =3 quadratic 
!         =4-10 etc 
!  xin-abscissa at which the ordinate value is desired 
!Output:
!  W_POLY_INTERP-value of ordinate at xin
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        npts,                    k_order
      REAL           x(*),                    y(*)
      REAL           xin
!Declaration of output variables
      REAL           W_POLY_INTERP
!Declaration of local variables
      INTEGER        i,                       ix,
     &               i1,                      i2,
     &               imax,                    ixmax,
     &               j,                       k,
     &               ntermt
      REAL           deltax,                  denom,
     &               prod,                    sum
      REAL           delta(10),               a(10)
!Search for appropriate value of x(1)
      ntermt=k_order
      IF(ntermt.lt.2) ntermt=2
      IF(ntermt.gt.npts) ntermt=npts
      DO i=1,npts
        IF(xin.eq.x(i)) THEN
!         Solution is at node i
          W_POLY_INTERP=y(i)
          GOTO 1000
        ELSEIF(xin.lt.x(i)) THEN
!         Bracketed solution
          i1=i-ntermt/2
          IF(i1.le.0) i1=1
          GOTO 10
        ENDIF
      ENDDO
!xin > x(npts)
      i1=npts-ntermt+1
   10 i2=i1+ntermt-1
      IF(npts.lt.i2) THEN
        i2=npts
        i1=i2-ntermt+1
        IF(i1.le.0) THEN
          i1=1
          ntermt=i2-i1+1
        ENDIF
      ENDIF
!Evaluate deviations delta
      denom=x(i1+1)-x(i1)
      deltax=(xin-x(i1))/denom
      DO i=1,ntermt
        ix=i1+i-1
        delta(i)=(x(ix)-x(i1))/denom
      ENDDO   
!Accumulate coefficients a
      a(1)=y(i1)
      DO k=2,ntermt
        prod=1.0
        sum=0.0
        imax=k-1
        ixmax=i1+imax
        DO i=1,imax
          j=k-i
          prod=prod*(delta(k)-delta(j))
          sum=sum-a(j)/prod
        ENDDO   
        a(k)=sum+y(ixmax)/prod
      ENDDO   
!Accumulate sum of expansion
      sum=a(1)
      DO j=2,ntermt
        prod=1.0
        imax=j-1
        DO i=1,imax
          prod=prod*(deltax-delta(i))
        ENDDO   
        sum=sum+a(j)*prod
      ENDDO   
      W_POLY_INTERP=sum
 1000 RETURN
      END
      