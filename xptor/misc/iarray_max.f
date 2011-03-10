      FUNCTION IARRAY_MAX(n,x,incx)
!***********************************************************************
!IARRAY_MAX is the index of the largest element of the x array
!References:
!  W.A.Houlberg 12/1998
!Input:
!  n-number of elements to be checked
!  x-array to be checked
!  incx-increment in x index
!Output:
!  IARRAY_XMAX-index of the largest element of the x array
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        incx,                    n
      REAL           x(*)
!Declaration of output variables
      INTEGER        IARRAY_MAX
!Declaration of local variables
      INTEGER        i
      IARRAY_MAX=1
      DO i=1,n,incx
        IF(x(i).gt.x(iarray_max)) IARRAY_MAX=i
      ENDDO   
      RETURN
      END
