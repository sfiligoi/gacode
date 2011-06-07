      SUBROUTINE RARRAY_ZERO(n,x)
!***********************************************************************
!RARRAY_ZERO sets the elements of array x to 0.0
!References:
!  W.A.Houlberg 12/1998
!Input:
!  n-number of elements to be zeroed
!Input/Output:
!  x-array to be zeroed
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        n
!Declaration of input/output variables
      REAL           x(*)
!Declaration of local variables
      INTEGER        i
      DO i=1,n
        x(i)=0.0
      ENDDO   
      RETURN
      END 
