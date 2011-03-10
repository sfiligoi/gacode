      SUBROUTINE RARRAY_SCALE(n,a,x,incx)
!***********************************************************************
!SSCAL scales the array x by the factor a
!References:
!  W.A.Houlberg 4/1997
!Input:
!  n-number of elements to be scaled
!  a-scale factor
!  incx-increment in sx index
!Input/Output:
!  x(i)-array to be scaled/scaled array
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        incx,                    n
      REAL           a
!Declaration of input/output variables
      REAL           x(*)
!Declaration of local variables
      INTEGER        i,                       ix
      ix=1
      DO i=1,n
        x(ix)=a*x(ix)
        ix=ix+incx
      ENDDO   
      RETURN
      END
