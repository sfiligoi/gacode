      SUBROUTINE RARRAY_SEARCH(x,n,xl,jlo)
!***********************************************************************
!RARRAY_SEARCH is a correlated table search routine to find the indices
!  of the array x(j) that bound xl
!References:
!  W.A.Houlberg 3/2000
!Input:
!  x(j)-monotonically increasing array
!  n-number of elements in x
!  xl-target value
!Input/Output:
!  jlo-input starting lower index
!     <1     binary search
!     =1,n-1 use value
!     >n-1   binary search
!  jlo-output starting lower index
!     =0     xl.lt.x(1) 
!     =1     x(1).le.xl.le.x(2)
!     =2,n-1 x(jlo).lt.xl.le.x(jlo+1)
!     =n     x(jlo).gt.x(n)
!Comments:
!  This is a modified version of the Numerical Recipes routine HUNT
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        n
      REAL           xl,                      x(*)
!Declaration of input/output variables
      INTEGER        jlo
!Declaration of local variables
      INTEGER        inc,                     jhi,
     &               jmid
      IF(xl.lt.x(1)) THEN
!       xl is out of range, below x(1)
        jlo=0                     
      ELSEIF(xl.le.x(2)) THEN
!       x(1).le.xl.le.x(2)
        jlo=1                     
      ELSEIF(xl.le.x(n)) THEN
!       Search for x(2).lt.xl.le.x(n)
!       Check if jlo from previous call is usable
        IF(jlo.lt.1.or.jlo.gt.(n-1)) THEN
          jlo=2
          jhi=n
        ELSE
!         Bracket xl
          inc=1
          IF(xl.gt.x(jlo)) THEN
!           Search up
            jhi=jlo+1
            DO WHILE(xl.gt.x(jhi))
              inc=inc+inc
              jlo=jhi
              jhi=MIN0(jlo+inc,n)
            ENDDO
          ELSE
!           Search down
            jhi=jlo
            jlo=jlo-1
            DO WHILE(xl.le.x(jlo))
              inc=inc+inc
              jhi=jlo
              jlo=MAX0(jlo-inc,1)
            ENDDO
          ENDIF
        ENDIF
!       Bisection
        DO WHILE(jhi-jlo.gt.1)
          jmid=(jhi+jlo)/2
          IF(xl.gt.x(jmid)) THEN
            jlo=jmid
          ELSE
            jhi=jmid
          ENDIF
        ENDDO
      ELSE
!       xl.gt.x(n)
        jlo=n
      ENDIF
      RETURN
      END
