      SUBROUTINE W_SPLINE(k_vopt,n,u,x,f,i,value)
!***********************************************************************
!W_SPLINE evaluates the cubic spline function and its derivatives
!  W.A.Houlberg, D.McCune 3/2000
!Input:
!  k_vopt(1)-calculate the function, 0=off, else=on
!  k_vopt(2)-calculate the first derivative, 0=off, else=on
!  k_vopt(3)-calculate the second derivative, 0=off, else=on
!  n-number of data points
!  u-abscissa at which the spline is to be evaluated
!  x-array containing the data abcissas
!  f(4,n)-array containing the data ordinates
!Input/Output:
!  i-guess for target lower bound input if 1<i<n
!  i-target lower bound on output
!Output:
!  value(i)-ordering is a subset of the sequence given under k_vopt
!  value(1)=function value or lowest order derivative requested
!  value(2)=next order derivative requested
!  value(3)=second derivative if all k_vopt are non-zero
!Comments:
!  s=f(1,i)+f(2,i)*dx+f(3,i)*dx**2/2!+f(4,i)*dx**3/3!
!  s'=f(2,i)+f(3,i)*dx+f(4,i)*dx**2/2!
!  s''=f(3,i)+f(4,i)*dx
!           where dx=u-x(i) and x(i).lt.u.lt.x(i+1)
!  If u.le.x(1) then i=1 is used
!  If u.ge.x(n) then i=n is used
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        k_vopt(*)
      INTEGER        n
      REAL           u
      REAL           f(4,*),                  x(*)
!Declaration of input/output variables
      INTEGER        i 
!Declaration of output variables
      REAL           value(*)
!Declaration of local variables
      INTEGER        j
      REAL           dx
      CALL RARRAY_SEARCH(x,n,u,i)
      IF(i.le.0) THEN
        i=1
        dx=0.0
      ELSEIF(i.ge.n) THEN
        i=n
        dx=0.0
      ELSE
        dx=u-x(i)
      ENDIF
!Evaluate spline
      j=0
      IF(k_vopt(1).ne.0) THEN
        j=j+1
        value(j)=f(1,i)+dx*(f(2,i)+0.5*dx*(f(3,i)+dx*f(4,i)/3.0))
      ENDIF
      IF(k_vopt(2).ne.0) THEN
        j=j+1
        value(j)=f(2,i)+dx*(f(3,i)+0.5*dx*f(4,i))
      ENDIF
      IF(k_vopt(3).ne.0) THEN
        j=j+1
        value(j)=f(3,i)+dx*f(4,i)
      ENDIF
      RETURN
      END
