      SUBROUTINE W_LIN_INTEG(x,y,n,sum)
!***********************************************************************
!W_LIN_INTEG integrates y with respect to x from x(1) to x(n) assuming
!  linear variation in y between nodes
!References:
!  W.A.Houlberg 4/1997
!Input:
!  x(i)-abscissas
!  y(i)-ordinates
!  n-number of nodes
!Output:
!  sum-integral
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables:
      INTEGER        n
      REAL           x(*),                    y(*)
!Declaration of output variables:
      REAL           sum
!Declaration of local variables
      INTEGER        i
      sum=0.0
      DO i=2,n
        sum=sum+(x(i)-x(i-1))*(y(i)+y(i-1))
      ENDDO
      sum=0.5*sum
      RETURN
      END
