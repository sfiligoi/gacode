

      SUBROUTINE strip_comments (input, output)
! ------------------------------------------------------------------------------------------
! remove comments from namelists 
! input is io unit number of file to be read
! output is io unit number of new file to be written wit comments stripped out
! ------------------------------------------------------------------------------HSJ---------
      USE nrtype,                             ONLY : DP,I4B,I2B
      USE error_handler,                      ONLY :  lerrno,iomaxerr,terminate 
      USE io_gcnmp,                            ONLY : nlog


      IMPLICIT NONE
!
      INTEGER(I4B)   INDEX, nchars, last, i,                           &
                     tab_loc,number_of_delimiters,maxchr
      INTEGER(I2B)   input,output
      PARAMETER                  ( maxchr = 256 )
      PARAMETER                  (number_of_delimiters = 3)
      CHARACTER comment_delimiter(number_of_delimiters)*1
      CHARACTER(LEN = maxchr) line
      CHARACTER  TAB
      DATA      comment_delimiter / ';', '!', '#' /
      TAB = CHAR ( 9 )
!
      DO WHILE (.TRUE.)
        READ (input, '(a)', END = 10) line
        nchars = LEN_TRIM(line)
        line   = 'X' // line(1:nchars)
        nchars = nchars+1
        IF (nchars .GT. maxchr) THEN
          lerrno = 10
          CALL terminate(lerrno,nlog)
        ELSE IF (nchars .GT. 1) THEN
          last   = nchars
          i      = 0
          DO WHILE (i .LT. number_of_delimiters .AND. last .GT. 1)
            i    = i + 1
            last = INDEX (line, comment_delimiter(i)) - 1
            IF (last .EQ. -1) &
            last = nchars
            line = line(1:last)
          END DO
!         replace  tabs with single blanks:
          tab_loc = 1
          DO WHILE(tab_loc .GT. 0)
             tab_loc = INDEX(line,TAB)
             IF(tab_loc .NE. 0)line(tab_loc:tab_loc) = ' '
          ENDDO
          WRITE (output, '(a)') line(2:last)
        END IF
      END DO
   10 RETURN
!
      END


#ifndef GCNMP

      SUBROUTINE get_beam_id (input,name_beam)
! ------------------------------------------------------------------------------------------
! Beam id line is a coment in nbwave generated namelists
! we need to get this identifier so we can group it witht he appropriate beam.
! input is io unit number of file to be read
! ------------------------------------------------------------------------------HSJ---------
      USE nrtype,                             ONLY : DP,I4B,I2B
      USE error_handler,                      ONLY :  lerrno,iomaxerr,terminate 
      USE io_gcnmp,                           ONLY : nlog
      USE beam_structure,                     ONLY : neutral_beam


      IMPLICIT NONE
!
      INTEGER(I4B)   INDEX,fc,lc,no_bm,nchars
      INTEGER(I2B)   input
      INTEGER(I4B),PARAMETER  ::   maxchr = 256 
      CHARACTER(LEN = maxchr) line
      TYPE(neutral_beam)name_beam

      no_bm = 0
      REWIND(input)

      DO WHILE (.TRUE.)
        READ (input, '(a)', END = 10) line
        nchars = LEN_TRIM(line)
        IF (nchars .GT. maxchr) THEN
          lerrno = 10
          CALL terminate(lerrno,nlog)
        ELSE IF (nchars .GT. 1) THEN
          fc  = INDEX (line,'BEAMLINE' ) -1
          IF(fc .gt. -1)THEN
             fc = fc +1
             no_bm = no_bm +1
             lc = INDEX (line,'%NBWAVE' ) -1
             DO WHILE(.TRUE.) ! strip trailing blanks
                IF(line(lc:lc) .NE. ' ')EXIT
                lc = lc-1
             ENDDO
             name_beam%beam_id(no_bm) = line(fc:lc)
          END IF
        END IF
      END DO

   10 REWIND(input)

      RETURN
!
      END

#endif
