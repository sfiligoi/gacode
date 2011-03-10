      FUNCTION RARRAY_XDOTY(n,x,incx,y,incy)
!***********************************************************************
!RARRAY_XDOTY is the dot product of the arrays x and y
!References:
!  W.A.Houlberg 12/1998
!Input: 
!  n-number of elements to be summed
!  x-array to be dotted
!  incx-increment in x index
!  y-array to be dotted
!  incy-increment in y index
!Output:
!  RARRAY_XDOTY-dot product of the arrays x and y
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        incx,                    incy,
     &               n
      REAL           x(*),                    y(*)
!Declaration of output variables
      REAL           RARRAY_XDOTY
!Declaration of local variables
      INTEGER        i,                       ix,
     &               iy
      RARRAY_XDOTY=0.0
      ix=1
      iy=1
      DO i=1,n
        RARRAY_XDOTY=RARRAY_XDOTY+x(ix)*y(iy)
        ix=ix+incx
        iy=iy+incy
      ENDDO   
      RETURN
      END
