      MODULE ENAME
!
!
      IMPLICIT NONE
      CHARACTER(len =64)    eqdskfilename, eqdskoldname, eqfile, &
                            intfl, eqdsk_tdem
!
!
!
      DATA  eqdskfilename,  eqdskoldname,  eqfile               &
         / 'none'        , 'none'       , 'none' /
!
!
!


      CONTAINS  

    SUBROUTINE eqdsk_name(eqdname)

! ----------------------------------------------------------------------
! construct an eqdsk file name
! ----------------------------------------------------------------------
!

      USE solcon,   only : steady_state,time
      USE yoka,     only :  ishot,itime
      character(len = *),  intent(OUT) :: eqdname
      character  itchar*6, ischar*7
      character (len = 24) intfl
      INTEGER itms,iashot




      IF(ABS(steady_state) .LT. 1.e-5)THEN
            itms = 99999
      ELSE
            itms = ABS (1000 * time)
      ENDIF
      IF (itms .LT. 10) THEN
            WRITE (intfl, 9000) itms
            READ  (intfl, 8025) itchar
          ELSE IF (itms .LT. 100   ) THEN
            WRITE (intfl, 9020) itms
            READ  (intfl, 8025) itchar
          ELSE IF (itms .LT. 1000  ) THEN
            WRITE (intfl, 9030) itms
            READ  (intfl, 8025) itchar
          ELSE IF (itms .LT. 10000 ) THEN
            WRITE (intfl, 9040) itms
            READ  (intfl, 8025) itchar
          ELSE IF (itms .LT. 100000) THEN
            WRITE (intfl, 9050) itms
            READ  (intfl, 8025) itchar
          ELSE
            itms = itms / 10
  200       IF (itms .GE. 100000) THEN
              itms = itms / 10
              go to 200
            END IF
            WRITE (intfl, 9050) itms
            READ  (intfl, 8025) itchar
      END IF
!
! --- itchar is now a 6-CHARACTER symbol; DO the same WITH the shot number
!
      iashot = IABS (ishot)
      IF (iashot .LT. 10) THEN
            WRITE (intfl, 9500) iashot
            READ  (intfl, 8025) ischar
          ELSE IF (iashot .LT. 100    ) THEN
            WRITE (intfl, 9520) iashot
            READ  (intfl, 8025) ischar
          ELSE IF (iashot .LT. 1000   ) THEN
            WRITE (intfl, 9530) iashot
            READ  (intfl, 8025) ischar
          ELSE IF (iashot .LT. 10000  ) THEN
            WRITE (intfl, 9540) iashot
            READ  (intfl, 8025) ischar
          ELSE IF (iashot .LT. 100000 ) THEN
            WRITE (intfl, 9550) iashot
            READ  (intfl, 8025) ischar
          ELSE IF (iashot .LT. 1000000) THEN
            WRITE (intfl, 9560) iashot
            READ  (intfl, 8025) ischar
          ELSE
            iashot = iashot / 10
  300       IF (iashot .GE. 1000000) THEN
              iashot = iashot / 10
              go to 300
            END IF
            WRITE (intfl, 9550) iashot
            READ  (intfl, 8025) ischar
      END IF
!
! --- ischar is now a 7-CHARACTER symbol
! --- form the file name (on CRAY we must truncate to 8 characters)
!
      WRITE (intfl, 9700)  itchar
      READ  (intfl, 8025)  eqdname
      eqdname = ADJUSTL(eqdname)

 8025 format (a)
 9000 format ('.0000' , i1)
 9020 format ('.000'  , i2)
 9030 format ('.00'   , i3)
 9040 format ('.0'    , i4)
 9050 format ('.'     , i5)
 9500 format ('g00000', i1)
 9520 format ('g0000' , i2)
 9530 format ('g000'  , i3)
 9540 format ('g00'   , i4)
 9550 format ('g0'    , i5)
 9560 format ('g'     , i6)
 9700 format ('g0'    ,  a)

      END SUBROUTINE  eqdsk_name


       END MODULE ENAME
