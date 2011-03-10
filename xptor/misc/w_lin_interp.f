      SUBROUTINE W_LIN_INTERP(n1,x1,y1,n2,x2,y2,iflag,message)
!***********************************************************************
!W_LIN_INTERP does a linear interpolation to obtain n2 values for the
!  target array y2 on the grid x2 from the n1 values of y1 on the grid
!  x1
!References:
!  W.A.Houlberg 3/2000
!Input:
!  n1-number of abscissas and values in the source arrays
!  x1-array of source abscissas
!  y1-array of source values
!  n2-number of target values to be found
!  x2-array of target abscissas
!Output:
!  y2-array of target values
!  iflag-error and warning flag
!       =-1 warning
!       =0 no warnings or errors
!       =1 error
!  message-warning or error message (character)
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        n1,                      n2
      REAL           x1(*),                   y1(*),
     &               x2(*)
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           y2(*)
!Declaration of local variables
      INTEGER        i,                       il
      IF(n1.lt.2) THEN
        iflag=1
        message='W_LIN_INTERP/ERROR:less than 2 points in source array'
      ELSE
        il=1
        DO i=1,n2
   10     IF(x2(i).lt.x1(1)) THEN
!           Use innermost data value
            y2(i)=y1(1)
            iflag=-1
            message='W_LIN_INTERP(1)/WARNING:x<x(1), use end point'
          ELSEIF(x2(i).eq.x1(1)) THEN
            y2(i)=y1(1)
          ELSEIF(x2(i).gt.x1(il+1)) THEN
            IF(il.lt.n1-1) THEN
!             Step x1 grid forward and loop
              il=il+1
              GOTO 10
            ELSE
!             Set to last value
              y2(i)=y1(n1)
              iflag=-1
              message='W_LIN_INTERP(2)/WARNING:x>x(n1), use end point'
            ENDIF
          ELSE
!           Interpolate
            y2(i)=y1(il)
     &            +(y1(il+1)-y1(il))*(x2(i)-x1(il))/(x1(il+1)-x1(il))
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END
