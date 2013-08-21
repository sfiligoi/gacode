    MODULE string_util
      USE nrtype,                        ONLY : DP,I4B,SP
#ifndef ONETWO
      USE  error_handler,                ONLY : lerrno,terminate,iomaxerr
      USE io_gcnmp,                             ONLY : ncrt,nlog
#endif
     CONTAINS

     SUBROUTINE strip_path(path,prefix,extension,filename,task,ishot,limit,ncrt)!JMP
!--------------------------------------------------------------------------
! return path, shot number and filename (with path removed)
! task   input values
!        if task   =0 then try to get shot number,
!        if task   =1 return after determining path and file names.
!------------------------------------------------------------HSJ-11/18/03--
   IMPLICIT NONE
   INTEGER j,k,l,ll,lll,task
   INTEGER, INTENT (out) :: ishot
   CHARACTER*(*), INTENT(inout):: filename
   CHARACTER*(*), INTENT(out):: path,prefix,extension !JMP
   INTEGER,  INTENT(in):: limit      !limit search to this many characters
                                     !ufiles has 64 character limit
   INTEGER ,INTENT(in) :: ncrt
   CHARACTER(len=6):: shot,buf


   l = LEN_TRIM(filename)
   IF(l .LE. limit)THEN
     DO j=l-1,2,-1
        k = j
        IF(filename(j:j) .EQ. '/')EXIT
     ENDDO
     IF(k .EQ. 2)THEN
        !no path specification found set to local path:
        path = "./"
     ELSE  ! found '/' in filename, assume it is path specifier
        path = filename(1:k)
        !rotate filename until it is at the begining of the string:
        filename(1:k) =' '
        filename = ADJUSTL(filename) 
        l = LEN_TRIM(filename)
     ENDIF 

     If(task == 1) RETURN ! return with just path and filename set


     ll =SCAN(filename(1:l),'0123456789')
     prefix = filename(1:ll-1)
     lll = SCAN(filename(ll:l),'.')
     extension = filename(lll+2:l) !JMP
     IF (lll-ll .LE. 6)THEN
            shot = filename(ll:lll)
            !convert to integer
            WRITE(buf,FMT='(a)')shot
            READ(buf,FMT='(i6)')ishot

     ELSE
        WRITE(ncrt,2)filename(1:l)
 2      FORMAT(2x,'could not determine shot no. for',a)
#ifdef ONETWO
        CALL STOP('strip_path',1)
#else
        lerrno = 265 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
     ENDIF
!     print *,'shot =',shot
!     print *,'ishot =',ishot
!     print *,'prefix =',prefix
!     print *,'in strip_path, path =',path(1:len_trim(path))
!     print *,'in strip_path,filename =',filename(1:l)
   ELSE
     WRITE(ncrt,1)filename(1:LEN_TRIM(filename))
 1   FORMAT(2x,'Error, ufile routines require that ',a,/, &
            2x,'be 64 characters or less in length')
#ifdef ONETWO
     CALL STOP('strip_path',2)
#else
        lerrno = 266 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
   ENDIF


   RETURN
   END SUBROUTINE strip_path



    FUNCTION to_upper_case(string) RESULT (UC_string)
! -------------------------------------------------------------------
!
      
      CHARACTER*(*), INTENT (IN) :: string
      CHARACTER(LEN(string))  UC_string
      UC_string = string
      DO  l =1,LEN_TRIM(string),1
         IF(LGE (string(l:l),'a') .AND. &
              LLE(string(l:l),'z'))     &
              UC_string(l:l) = ACHAR(IACHAR(string(l:l)) - 32)
      ENDDO
      END FUNCTION to_upper_case

      SUBROUTINE to_upper_case1(string)
! -------------------------------------------------------------------
!
      CHARACTER*(*), INTENT (INOUT) :: string
      INTEGER l
      DO  l =1,LEN_TRIM(string),1
         IF(LGE (string(l:l),'a') .AND. &
              LLE(string(l:l),'z'))     &
              string(l:l) = ACHAR(IACHAR(string(l:l)) - 32)
      ENDDO
      RETURN
      END SUBROUTINE to_upper_case1


   FUNCTION to_lower_case(string) RESULT (LC_string)
! -------------------------------------------------------------------
!
      
      CHARACTER*(*), INTENT (IN) :: string
      CHARACTER(LEN(string)) LC_string
      LC_string = string
      DO  l =1,LEN_TRIM(string),1
         IF(LGE (string(l:l),'A') .AND. &
              LLE(string(l:l),'Z'))     &
              LC_string(l:l) = ACHAR(IACHAR(string(l:l)) + 32)
      ENDDO
      END FUNCTION to_lower_case

      SUBROUTINE  to_lower_case1(string)
! -------------------------------------------------------------------
!     FUNCTION FORM NOR ACCEPTED BY PGF90  HSJ
      
      CHARACTER*(*), INTENT (INOUT) :: string
      INTEGER l
      DO  l =1,LEN_TRIM(string),1
         IF(LGE (string(l:l),'A') .AND. &
              LLE(string(l:l),'Z'))     &
              string(l:l) = ACHAR(IACHAR(string(l:l)) + 32)
      ENDDO
      RETURN
      END SUBROUTINE to_lower_case1



      SUBROUTINE  set_r4_1darray_value(k,kl,line1,value)
! -----------------------------------------------------------------------
! --------------------------------------------------------------HSJ------
       IMPLICIT NONE
       CHARACTER*(*) line1
       CHARACTER (len = 256) line
       CHARACTER(len =4) cha
       CHARACTER(len =1) cha1,chs
       CHARACTER (len =24)intfile,char
       REAL (SP) ,DIMENSION(:) :: value
       INTEGER kbl,kbr,kbe,k,kk,kbs,kbb,inta,nvalu,nvall,j,JJ,KNE,intr,kl,kn
       REAL (SP) DATA

       nvalu = UBOUND(value,1)
       nvall = LBOUND(value,1)
       kn = k-kl+1
       line(1:kn) = line1(kl:k)
       kn=kn+1
       line(kn:kn) =' '
       k=kn

       kbl = INDEX(line(1:k),'(')
       kbr = INDEX(line(1:k),')')
       kbe = INDEX(line(1:k),'=')
       intr = 0
       IF(kbl .GE. kbr)THEN
                        !no parentheses found. assume array input without indecies
            IF(kbe .GT. 0 .AND. kbe .LT. k)THEN
               inta = 0
               !determine number of values present and setup indecies
               kbe = kbe+1
               chs = ' '
               kbs =kbe
 30              DO j= kbs,k !scan until we hit a single digit or + or - :
                  kk = j
                  cha1 = line(j:j)
                  IF(cha1 .EQ. '-' .OR. cha1 .EQ. '+')THEN
                     chs = cha
                  ELSE IF(cha1 .EQ. '0' .OR. cha1 .EQ.'1' .OR.   &
                         cha1 .EQ.'2' .OR. cha1 .EQ.'3'  .OR.    &
                         cha1 .EQ.'4' .OR. cha1 .EQ.'5'  .OR.    &
                         cha1 .EQ.'6' .OR. cha1 .EQ.'7'  .OR.    &
                         cha1 .EQ.'8' .OR. cha1 .EQ.'9')THEN
                         kbe = j
                         !number starts at position kbe, find where it ends
                         ! by searching for comma or blank or *
                         DO jj = kbe+1,k
                            cha1 = line(jj:jj)
                            IF(cha1 .EQ. ' ' .OR. cha1 .EQ. ',' )THEN
                               kne = jj-1
                               go to 10
                            ENDIF
                            IF(cha1 .EQ. '*')THEN  !repeat field
                               kne = jj -1
                               cha = LINE(KBE:KNE)
                               !print *,'repeat field =',cha
                               WRITE(intfile,FMT='(A)')cha
                               READ(intfile,FMT='(i4)')intr
                               go to 20
                            ENDIF
                         ENDDO
                         kne = k+1 !number extends to end of line
 10                      CONTINUE !number extends from kbe to kne-1
                         char = chs//LINE(KBE:KNE)
                         WRITE(intfile,FMT='(A)')char
                         !print *,'line =',line(1:k)
                         !print *,'char =',char
                         READ(intfile,FMT='(e16.8)')DATA
                         IF(intr .EQ. 0)THEN
                            inta = inta+1
                            IF(inta .GE. nvall .AND. inta .LE. nvalu)THEN
                               value(inta)= DATA
                               kbs = kne+1
                               go to 25
                            ELSE
#ifdef ONETWO
                               CALL STOP('set_r4_1darray 10',1)
#else
        lerrno = 267 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
                            ENDIF
                         ELSE
                            IF(1 .LE. nvall .AND. intr .LE. nvalu)THEN
                               DO jj =1,intr
                                  VALUE(JJ) =DATA
                               ENDDO
                               RETURN
                            ELSE
#ifdef ONETWO
                               CALL STOP('set_r4_1darray 11',1)
#else
       lerrno = 268 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
                            ENDIF
                         ENDIF
                         chs = ' '
                         !print *,'value(inta),inta =',value(inta),inta
                  ENDIF
 20               ENDDO
 25               CONTINUE
                  !print *,'kk,k,j,kbs =',kk,k,j,kbs
                  !print *,'value(inta),inta =',value(inta),inta
                 IF(kk .LT. k .AND. kbs .LT. k) go to 30
            ELSE
#ifdef ONETWO
               CALL STOP('set_r4_1darray_value a',1)
#else
        lerrno = 269 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
            ENDIF
       ELSE
          !PRINT *,LINE(1:K)
          cha =line(kbl+1:kbr-1)
          WRITE(intfile,FMT='(A)')cha
          READ(intfile,FMT='(i4)')inta
          IF(inta .GE. nvall .AND. inta .LE. nvalu)THEN
            kbe = INDEX(line(1:k),'=')
            IF(kbe .GT. 0)THEN
               DO j = kbe+1,k
                  IF(line(j:j) .NE. ' ')THEN !skip blanks after '='
                    jj = j
                    go to 15
                  ENDIF
               ENDDO
 15            kbb = INDEX(line(jj:k),' ') + jj-1
               !PRINT *,'KBB,KBE =',KBB,KBE
               IF(kbb == 0 .OR. kbb .GT. kbe+2)THEN !kbb=0 ==> no blanks
                   IF(kbb == 0)kbb = k+1
                   char = line(kbe+1:kbb-1)
                   WRITE(intfile,FMT='(A)')char
                   READ(intfile,FMT='(e16.8)')DATA
                   value(inta) = DATA
                   !PRINT *,'VALUE,INTA =',VALUE,INTA
               ELSE
                   RETURN
               ENDIF
            ELSE
              RETURN
            ENDIF
          ELSE
              PRINT *,'error, index out of range :'
              PRINT *,'lower bound =',nvall
              PRINT *,'upper bound =',nvalu
              PRINT *,'index = ',inta
#ifdef ONETWO
              CALL STOP('set_r4_1darray_value',2)
#else
        lerrno = 270 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
          ENDIF
       ENDIF 
       RETURN
       END       SUBROUTINE  set_r4_1darray_value



      SUBROUTINE  set_i4_1darray_value(k,kl,line,value)
! -----------------------------------------------------------------------
! --------------------------------------------------------------HSJ------
       IMPLICIT NONE
       CHARACTER*(*) line
       CHARACTER(len =4) cha
       CHARACTER(len =1) cha1,chs
       CHARACTER (len =24)intfile,char
       INTEGER ,DIMENSION(:) :: value
       INTEGER kbl,kbr,kbe,k,kbb,inta,nvalu,nvall,j,JJ,KNE,intr,kl
       INTEGER DATA

       nvalu = UBOUND(value,1)
       nvall = LBOUND(value,1)
       kbl = INDEX(line(1:k),'(')
       kbr = INDEX(line(1:k),')')
       kbe = INDEX(line(1:k),'=')
       intr = 0
       IF(kbl .GE. kbr)THEN
                        !no parentheses found. assume array input without indecies
            IF(kbe .GT. 0 .AND. kbe .LT. k)THEN
               inta = 0
               !determine number of values present and setup indecies
               kbe = kbe+1
               chs = ' '
               DO j= kbe,k !scan until we hit a single digit or + or - :
                  cha1 = line(j:j)
                  IF(cha1 .EQ. '-' .OR. cha1 .EQ. '+')THEN
                     chs = cha
                  ELSE IF(cha1 .EQ. '0' .OR. cha1 .EQ.'1' .OR.   &
                         cha1 .EQ.'2' .OR. cha1 .EQ.'3'  .OR.    &
                         cha1 .EQ.'4' .OR. cha1 .EQ.'5'  .OR.    &
                         cha1 .EQ.'6' .OR. cha1 .EQ.'7'  .OR.    &
                         cha1 .EQ.'8' .OR. cha1 .EQ.'9')THEN
                         kbe = j
                         !number starts at position kbe, find where it ends
                         ! by searching for comma or blanck or *
                         DO jj = kbe+1,k
                            cha1 = line(jj:jj)
                            IF(cha1 .EQ. ' ' .OR. cha1 .EQ. ',' )THEN
                               kne = jj-1
                               go to 10
                            ENDIF
                            IF(cha1 .EQ. '*')THEN  !repeat field
                               kne = jj -1
                               cha = LINE(KBE:KNE)
                               !print *,'repeat field =',cha
                               WRITE(intfile,FMT='(A)')cha
                               READ(intfile,FMT='(i4)')intr
                               go to 20
                            ENDIF
                         ENDDO
                         kne = k+1 !number extends to end of line
 10                      CONTINUE !number extends from kbe to kne-1
                         char = chs//LINE(KBE:KNE)
                         WRITE(intfile,FMT='(A)')char
                         !print *,'line =',line(1:k)
                         !print *,'char =',char
                         READ(intfile,FMT='(i8)')DATA
                         IF(intr .EQ. 0)THEN
                            inta = inta+1
                            value(inta)= DATA
                         ELSE
                            DO jj =1,intr
                               VALUE(JJ) =DATA
                            ENDDO
                            RETURN
                         ENDIF
                         chs = ' '
                         !print *,'value(inta),inta =',value(inta),inta
                  ENDIF
 20               ENDDO
            ELSE
#ifdef ONETWO
               CALL STOP('set_r8_1darray_value a',1)
#else
        lerrno = 271 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
            ENDIF
       ELSE
          cha =line(kbl+1:kbr-1)
          WRITE(intfile,FMT='(A)')cha
          READ(intfile,FMT='(i4)')inta
          IF(inta .GE. nvall .AND. inta .LE. nvalu)THEN
            kbe = INDEX(line(1:k),'=')
            IF(kbe .GT. 0)THEN
               kbb = INDEX(line(kbe+1:k),' ',BACK=.TRUE.)
               IF(kbb == 0 .OR. kbb .GT. kbe+2)THEN !kbb=0 ==> no blanks
                   IF(kbb == 0)kbb = k+1
                   char = line(kbe+1:kbb-1)
                   WRITE(intfile,FMT='(A)')char
                   READ(intfile,FMT='(i8)')DATA
                   value(inta) = DATA
               ELSE
                   RETURN
               ENDIF
            ELSE
              RETURN
            ENDIF
          ELSE
              PRINT *,'error, index out of range :'
              PRINT *,'lower bound =',nvall
              PRINT *,'upper bound =',nvalu
              PRINT *,'index = ',inta
#ifdef ONETWO
              !call STOP('set_i4_1darray_value',1)
#else
        lerrno = 272 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
          ENDIF
       ENDIF 

       RETURN
       END        SUBROUTINE  set_i4_1darray_value


      SUBROUTINE  set_L_1darray_value(k,kl,line,value)
! -----------------------------------------------------------------------
! --------------------------------------------------------------HSJ------
       IMPLICIT NONE
       CHARACTER*(*) line
       CHARACTER(len =4) cha
       CHARACTER(len =1) cha1,chs
       CHARACTER (len =24)intfile,char
       LOGICAL,DIMENSION(:) :: value
       INTEGER kbl,kbr,kbe,k,kbb,inta,nvalu,nvall,j,JJ,KNE,intr,kl,KBB2
       LOGICAL * 1  DATA

       nvalu = UBOUND(value,1)
       nvall = LBOUND(value,1)
       kbl = INDEX(line(1:k),'(')
       kbr = INDEX(line(1:k),')')
       kbe = INDEX(line(1:k),'=')
       intr = 0
       IF(kbl .GE. kbr)THEN
            !no parentheses found. assume array input without indecies
            IF(kbe .GT. 0 .AND. kbe .LT. k)THEN
               inta = 0
               !determine number of values present and setup indecies
               kbe = kbe+1
               chs = ' '
               DO j= kbe,k !scan until we hit a single digit or + or .
                  cha1 = line(j:j)
                  IF(cha1 .EQ. '.' .OR. cha1 .EQ. '+')THEN
                     chs ='.'
                  ELSE IF(cha1 .EQ. '0' .OR. cha1 .EQ.'1' .OR.   &
                         cha1 .EQ.'2' .OR. cha1 .EQ.'3'  .OR.    &
                         cha1 .EQ.'4' .OR. cha1 .EQ.'5'  .OR.    &
                         cha1 .EQ.'6' .OR. cha1 .EQ.'7'  .OR.    &
                         cha1 .EQ.'8' .OR. cha1 .EQ.'9')THEN
                         kbe = j
                         !number starts at position kbe, find where it ends
                         ! by searching for comma or blanck or *
                         DO jj = kbe+1,k
                            cha1 = line(jj:jj)
                            IF(cha1 .EQ. ' ' .OR. cha1 .EQ. ',' )THEN
                               kne = jj-1
                               go to 10
                            ENDIF
                            IF(cha1 .EQ. '*')THEN  !repeat field
                               kne = jj -1
                               cha = LINE(KBE:KNE)
                               !print *,'repeat field =',cha
                               WRITE(intfile,FMT='(A)')cha
                               READ(intfile,FMT='(i4)')intr
                               go to 20
                            ENDIF
                         ENDDO
                         kne = k+1 !number extends to end of line
 10                      CONTINUE !number extends from kbe to kne-1
                         char = chs//LINE(KBE:KNE)
                         WRITE(intfile,FMT='(A)')char
                         !print *,'line =',line(1:k)
                         !print *,'char =',char
                         READ(intfile,FMT='(l)')DATA
                         IF(intr .EQ. 0)THEN
                            inta = inta+1
                            value(inta)= DATA
                         ELSE
                            DO jj =1,intr
                               VALUE(JJ) =DATA
                            ENDDO
#ifdef ONETWO
                            CALL STOP('set_L_1darray_value a',1)
#else
        lerrno = 273 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
                            RETURN
                         ENDIF
                         chs = ' '
                         !print *,'value(inta),inta =',value(inta),inta
                     ELSE IF(cha1 .EQ. 'T' .OR. cha1 .EQ.'t') THEN
                        CHAR(1:4) = line(j:j+3)
                        char = to_upper_case(char)
                        IF(CHAR(1:4) .EQ. 'TRUE')THEN
                          DATA = .TRUE.
                          IF(intr .EQ. 0)THEN
                            inta = inta+1
                            value(inta)= DATA
                          ELSE
                            DO jj =1,intr
                               VALUE(JJ) =DATA
                            ENDDO
                            !print *,'value(1:intr) =',value(1:intr)
                            RETURN
                          ENDIF
                        ENDIF
                     ELSE IF(cha1 .EQ. 'F' .OR. cha1 .EQ.'f') THEN
                        char = line(j:j+3)
                        char = to_upper_case(char)
                        IF(CHAR(1:5) .EQ. 'FALSE')THEN
                           DATA = .FALSE.
                           IF(intr .EQ. 0)THEN
                            inta = inta+1
                            value(inta)= DATA
                           ELSE
                            DO jj =1,intr
                               VALUE(JJ) =DATA
                            ENDDO
                            !print *,'value(1:intr) =',value(1:intr)
                            RETURN
                           ENDIF
                        ENDIF
                     ENDIF
 20               ENDDO
            ELSE
#ifdef ONETWO
               CALL STOP('set_L_1darray_value a',2)
#else
        lerrno = 274 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
            ENDIF
       ELSE
          cha =line(kbl+1:kbr-1)
          WRITE(intfile,FMT='(A)')cha
          READ(intfile,FMT='(i4)')inta
          IF(inta .GE. nvall .AND. inta .LE. nvalu)THEN
            kbe = INDEX(line(1:k),'=')
            IF(kbe .GT. 0)THEN
               kbb = INDEX(line(kbe+1:k),'.') !kbb is relative to kbe
               IF(kbb == 0 ) THEN
                   PRINT*,'parse error in set_l_1darray_value'
                   PRINT *,'current line :'
                   PRINT *,line(1:k)
               ENDIF
               IF( kbb .ne. 0) kbb =kbb + kbe  !make kbb count from start of line
               kbb2 = INDEX(line(kbb+1:k),'.')
               IF(kbb2 .EQ.0)THEN
                      PRINT *,'parse error in set_l_1darray_value 1'
                      PRINT *,'current line :'
                      PRINT *,line(1:k)
               ENDIF
               if(kbb2 .ne. 0) kbb2 = kbb + kbb2  !make kbb2  count from start of line
               IF ( kbb2 .NE. 0 .AND. kbb .NE. 0)THEN
                  CHAR =LINE(KBB:KBB2+1)
                  WRITE(intfile,FMT='(A)')char
                  READ(intfile,FMT='(l)')DATA
                  value(inta) = DATA
               ELSE
#ifdef ONETWO
                  CALL STOP('parse error',1)
#else
        lerrno = 275 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
               ENDIF
               
            ELSE
              RETURN
            ENDIF
          ELSE
              PRINT *,'error, index out of range :'
              PRINT *,'lower bound =',nvall
              PRINT *,'upper bound =',nvalu
              PRINT *,'index = ',inta
              !call STOP('set_L_1darray_value',1)
          ENDIF
       ENDIF 

       RETURN
       END        SUBROUTINE  set_L_1darray_value



      SUBROUTINE  set_r8_1darray_value(k,kl,line1,value)
! -----------------------------------------------------------------------
! --------------------------------------------------------------HSJ------
       IMPLICIT NONE
       CHARACTER*(*) line1
       CHARACTER (len = 256) line
       CHARACTER(len =4) cha
       CHARACTER(len =1) cha1,chs
       CHARACTER (len =24)intfile,char
       REAL (DP) ,DIMENSION(:) :: value
       INTEGER kbl,kbr,kbe,k,kk,kbs,kbb,inta,nvalu,nvall,j,JJ,KNE,intr,kl,kn
       REAL (DP) DATA

       nvalu = UBOUND(value,1)
       nvall = LBOUND(value,1)

       kn = k-kl+1
       line(1:kn) = line1(kl:k)
       kn=kn+1
       line(kn:kn) =' '
       k=kn
       kbl = INDEX(line(1:k),'(')
       kbr = INDEX(line(1:k),')')
       kbe = INDEX(line(1:k),'=')
       intr = 0
       IF(kbl .GE. kbr)THEN
            !no parentheses found. assume array input without indecies
            IF(kbe .GT. 0 .AND. kbe .LT. k)THEN
               inta = 0
               !determine number of values present and setup indecies
               kbe = kbe+1
               chs = ' '
               kbs =kbe
 30              DO j= kbs,k !scan until we hit a single digit or + or - :
                  kk = j
                  cha1 = line(j:j)
                  IF(cha1 .EQ. '-' .OR. cha1 .EQ. '+')THEN
                     chs = cha
                  ELSE IF(cha1 .EQ. '0' .OR. cha1 .EQ.'1' .OR.   &
                         cha1 .EQ.'2' .OR. cha1 .EQ.'3'  .OR.    &
                         cha1 .EQ.'4' .OR. cha1 .EQ.'5'  .OR.    &
                         cha1 .EQ.'6' .OR. cha1 .EQ.'7'  .OR.    &
                         cha1 .EQ.'8' .OR. cha1 .EQ.'9')THEN
                         kbe = j
                         !number starts at position kbe, find where it ends
                         ! by searching for comma or blanck or *
                         DO jj = kbe+1,k
                            cha1 = line(jj:jj)
                            IF(cha1 .EQ. ' ' .OR. cha1 .EQ. ',' )THEN
                               kne = jj-1
                               go to 10
                            ENDIF
                            IF(cha1 .EQ. '*')THEN  !repeat field
                               kne = jj -1
                               cha = LINE(KBE:KNE)
                               !print *,'repeat field =',cha
                               WRITE(intfile,FMT='(A)')cha
                               READ(intfile,FMT='(i4)')intr
                               go to 20
                            ENDIF
                         ENDDO
                         kne = k+1 !number extends to end of line
 10                      CONTINUE !number extends from kbe to kne-1
                         char = chs//LINE(KBE:KNE)
                         WRITE(intfile,FMT='(A)')char
                         !print *,'line =',line(1:k)
                         !print *,'char =',char
                         READ(intfile,FMT='(e16.8)')DATA
                         IF(intr .EQ. 0)THEN
                            inta = inta+1
                            IF(inta .GE. nvall .AND. inta .LE. nvalu)THEN
                               value(inta)= DATA
                               kbs = kne+1
                               go to 25
                            ELSE
                               PRINT *,'error in array dims'
                               PRINT *,'current index =',inta
                               PRINT *,'array range =',nvall,nvalu
#ifdef ONETWO
                               CALL STOP('set_r8_1darray 10',3)
#else
        lerrno = 276 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
                            ENDIF
                         ELSE
                            IF(1 .LE. nvall .AND. intr .LE. nvalu)THEN
                               DO jj =1,intr
                                  VALUE(JJ) =DATA
                               ENDDO
                               RETURN
                            ELSE
#ifdef ONETWO
                               CALL STOP('set_r8_1darray 11',4)
#else
       lerrno = 277 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
                            ENDIF
                         ENDIF
                         chs = ' '
                         !print *,'value(inta),inta =',value(inta),inta
                  ENDIF
 20               ENDDO
 25               CONTINUE
                  !print *,'kk,k,j,kbs =',kk,k,j,kbs
                  !print *,'value(inta),inta =',value(inta),inta
                 IF(kk .LT. k .AND. kbs .LT. k) go to 30
            ELSE
#ifdef ONETWO
               CALL STOP('set_r8_1darray_value a',1)
#else
        lerrno = 278 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
            ENDIF
       ELSE
          !PRINT *,LINE(1:K)
          cha =line(kbl+1:kbr-1)
          WRITE(intfile,FMT='(A)')cha
          READ(intfile,FMT='(i4)')inta
          IF(inta .GE. nvall .AND. inta .LE. nvalu)THEN
            kbe = INDEX(line(1:k),'=')
            IF(kbe .GT. 0)THEN
               DO j = kbe+1,k
                  IF(line(j:j) .NE. ' ')THEN !skip blanks after '='
                    jj = j
                    go to 15
                  ENDIF
               ENDDO
 15            kbb = INDEX(line(jj:k),' ') + jj-1
               !PRINT *,'KBB,KBE =',KBB,KBE
               IF(kbb == 0 .OR. kbb .GT. kbe+2)THEN !kbb=0 ==> no blanks
                   IF(kbb == 0)kbb = k+1
                   char = line(kbe+1:kbb-1)
                   WRITE(intfile,FMT='(A)')char
                   READ(intfile,FMT='(e16.8)')DATA
                   value(inta) = DATA
                   !PRINT *,'VALUE,INTA =',VALUE,INTA
               ELSE
                   RETURN
               ENDIF
            ELSE
              RETURN
            ENDIF
          ELSE
              PRINT *,'error, index out of range :'
              PRINT *,'lower bound =',nvall
              PRINT *,'upper bound =',nvalu
              PRINT *,'index = ',inta
#ifdef ONETWO
              CALL STOP('set_r8_1darray_value',2)
#else
        lerrno = 279 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
          ENDIF
       ENDIF 

       RETURN
       END        SUBROUTINE  set_r8_1darray_value



      SUBROUTINE  set_r8_value(k,kl,line,value)
! -----------------------------------------------------------------------
! (your'e right I dont believe in function or operator  overloading.)
! --------------------------------------------------------------HSJ------
       IMPLICIT NONE
       INTEGER k,kbe,kbb,kl
       CHARACTER*(*) line
       CHARACTER(len =4) cha
       CHARACTER (len =24)intfile,char
       REAL (DP) value,DATA

            kbe = INDEX(line(1:k),'=')
            IF(kbe .GT. 0)THEN
               kbb = INDEX(line(kbe+1:k),' ',BACK=.TRUE.)
               IF(kbb .GT. kbe )THEN ! some blanks after '='
                   char = line(kbb:k)
                   WRITE(intfile,FMT='(A)')char
                   READ(intfile,FMT='(e16.8)')DATA
                   value = DATA
               ELSE  ! no blanks after '='
                   char = line(kbe+1:k)
                   WRITE(intfile,FMT='(A)')char
                   READ(intfile,FMT='(e16.8)')DATA
                   value = DATA
                   RETURN
               ENDIF
            ELSE
#ifdef ONETWO
              CALL STOP('set_r8_value failed',1)
#else
        lerrno = 280 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
              RETURN
            ENDIF

        RETURN
        END   SUBROUTINE  set_r8_value


      SUBROUTINE  set_r4_value(k,kl,line,value)
! -----------------------------------------------------------------------
! (your'e right I dont believe in function or operator  overloading.)
! --------------------------------------------------------------HSJ------
       IMPLICIT NONE
       INTEGER k,kbe,kbb,kl
       CHARACTER*(*) line
       CHARACTER(len =4) cha
       CHARACTER (len =24)intfile,char
       REAL (SP) value,DATA

            kbe = INDEX(line(1:k),'=')
            IF(kbe .GT. 0)THEN
               kbb = INDEX(line(kbe+1:k),' ',BACK=.TRUE.)
               IF(kbb == 0 .OR. kbb .GT. kbe+2)THEN !kbb=0 ==> no blanks
                   IF(kbb == 0)kbb = k+1
                   char = line(kbe+1:kbb-1)
                   WRITE(intfile,FMT='(A)')char
                   READ(intfile,FMT='(e16.8)')DATA
                   value = DATA
               ELSE
                    RETURN
               ENDIF
            ELSE
              RETURN
            ENDIF

        RETURN
        END        SUBROUTINE  set_r4_value



      SUBROUTINE  set_integer_valueo(k,kl,line,value)
! -----------------------------------------------------------------------
! (your'e right I dont believe in function or operator  overloading.)
! --------------------------------------------------------------HSJ------
       IMPLICIT NONE
       INTEGER k,kbe,kbb,kl
       CHARACTER*(*) line
       CHARACTER (len =24)intfile,char
       INTEGER value,DATA
       LOGICAL back

            kbe = INDEX(line(1:k),'=')
            IF(kbe .GT. 0)THEN
               kbb = INDEX(line(kbe+1:k),' ',BACK=.TRUE.)
               kbb =kbe+kbb       !kbb is returned as relative to kbe+1 above
               IF(kbb == 0 .OR. kbb .GT. kbe+2)THEN !kbb=0 ==> no blanks
                   IF(kbb == 0)kbb = k+1
                   char = line(kbe+1:kbb-1)
                   WRITE(intfile,FMT='(A)')char
                   PRINT *,'k,kl,kbe,kbb =',k,kl,kbe,kbb
                   PRINT *,'char =',char
                   READ(intfile,FMT='(i10)')DATA
                   value = DATA
               ELSE
                   PRINT *,'k,kbe,kbb =',k,kbe,kbb
                   PRINT *,'line =',line(1:k)
                   RETURN
               ENDIF
            ELSE
               PRINT *,'k,kbe =',k,kbe
               PRINT *,'line =',line(1:k)
              RETURN
            ENDIF

        RETURN
        END       SUBROUTINE  set_integer_valueo


      SUBROUTINE  set_integer_value(k,kl,line,value)
! -----------------------------------------------------------------------
! (your'e right I dont believe in function or operator  overloading.)
! --------------------------------------------------------------HSJ------
       IMPLICIT NONE
       INTEGER k,kbe,kbb,kl
       CHARACTER*(*) line
       CHARACTER (len =24)intfile,char
       INTEGER value,DATA,j
       LOGICAL back

            kbe = INDEX(line(1:k),'=')
            IF(kbe .GT. 0)THEN
               !search for first non blank after = sign:
               kbb=0
               DO j=kbe+1,k
                  IF(line(j:j) .NE. ' ')kbb=j
                  IF(kbb .NE. 0)EXIT
               ENDDO
               IF(kbb ==0)THEN
                  PRINT *,'line =',line
#ifdef ONETWO
                  CALL STOP("set_integer_value",1) 
#else
                  lerrno = 281 + iomaxerr
                  CALL terminate(lerrno,nlog)
#endif
               ENDIF
               !search for firt blank after kbb
               kbe = 0
               DO j =kbb+1,k
                  IF(line(j:j) .EQ. ' ')kbe = j
                  IF(kbe .NE. 0)EXIT
               ENDDO
               IF(kbe .EQ. 0)kbe =k
               char = line(kbb:kbe)
                   WRITE(intfile,FMT='(A)')char
                   READ(intfile,FMT='(i10)')DATA
                   value = DATA
                   RETURN
            ELSE
              RETURN
            ENDIF

        RETURN
        END       SUBROUTINE  set_integer_value

      SUBROUTINE  set_L_value(k,kl,line,value)
! -----------------------------------------------------------------------
! (your'e right I dont believe in function or operator  overloading.)
! --------------------------------------------------------------HSJ------
       IMPLICIT NONE
       INTEGER k,kbe,kbb,kl
       CHARACTER*(*) line
       CHARACTER(len =4) cha
       CHARACTER (len =24)intfile,char
       LOGICAL value,DATA

            kbe = INDEX(line(1:k),'=')
            IF(kbe .GT. 0)THEN
               kbb = INDEX(line(kbe+1:k),' ',BACK=.TRUE.)
               IF(kbb == 0 .OR. kbb .GT. kbe+2)THEN !kbb=0 ==> no blanks
                   IF(kbb == 0)kbb = k+1
                   char = line(kbe+1:kbb-1)
                   WRITE(intfile,FMT='(A)')char
                   READ(intfile,FMT='(l)')DATA
                   value = DATA
               ELSE
                    RETURN
               ENDIF
            ELSE
              RETURN
            ENDIF

        RETURN
        END         SUBROUTINE  set_L_value






   SUBROUTINE readblock(sstr,rdlbl,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1, &
                             i2,i3,i4,skip_summary,swap_lbl,aout,value_intg)
! -----------------------------------------------------------HSJ
! INPUT 
!  sstr is the string to find  in line
!  rdlbl = false  if sstr does not have a label associated with it
!  rdlbl = true  if sstr is a quantity defined further by a label
!  routine assumes that line is the curret line to search for
!  sstr. On return line is blank!
!  line is assumed tobe in all uppercase on input
!  io is fortran unit number to read next line from
!  nj size of rho grid in ONETWO
!  type ="V" or type ="A" means data to be read is on zone boundaries
!      zbuf(1) = r(2) ,... zbuf(nj-1) = r(nj) 
!  type = anything else means  data to be read  is  on zone centers
!  cvrt conversion factor to be applied to aout
!  zthlbls thermal species labels
!  skip_summary not all outputs have summaries
!  swap_lbl reverse some lables so they can be read with existing code
! OUTPUT
!  i2 = -1 if sstr is not contained in line
!  i2 =0   if rdlbl = .false.
!  i2 =    reaction index if rdlbl = true and lable ="REACTION**"
!  i2 =    zflbl index  if rdlbl = true and  lable = "**_MCBEAM"
!  i1 = 0  if not a beam energy indexed quantity
!  i1 = 1 if indexed by first beam energy component
!  i1 = 2 if indexed by second beam energy component
!  i1 = 3 if indexed by third beam energy component
!  i3 = 
!  aout(1:kj)    values on [0,1]  UNIFORM !!!!! rho grid 
! -----------------------------------------------------------------
!       USE param,ONLY : kj  HSJ 3/29/11
       IMPLICIT NONE
       LOGICAL,INTENT(in)  :: rdlbl,skip_summary,swap_lbl
       CHARACTER(len =*) sstr,line
       CHARACTER(len = 2) chrf
       CHARACTER(len = 3) chrf1
       CHARACTER(len = 1) TYPE
       CHARACTER(len = 132) intfile,lable,numbr,lbl1,lbl2
       CHARACTER(len =*),DIMENSION(:)   :: zflbls,zthlbls
       INTEGER io,k,kk,kr,njm1,nj,i2,l,i1,i3,i4,kc,kcf,kcb
!       REAL (DP)   aout(*),zbuf(kj),cvrt,value_intg
       REAL (DP)   aout(*),cvrt,value_intg
       REAL(DP) zbuf(nj) ! TEMPORARY ARRAY   HSJ 3/29/11


       i1 = 0 ;i2 =0 ; i3 =0
       njm1  = nj-1     ! grid size used in nubeam for zbuf array values
       kk = 0
       k = LEN_TRIM(line) 
       kk= INDEX(line(1:k),sstr(1:LEN(sstr)))
       IF(kk .EQ. 0)THEN
          i2 = -1                        !line does not have sstr in it
          RETURN
       ENDIF
      IF(rdlbl)THEN                       !sstr requires additional info
         kk = INDEX(line(1:k),'(')       !look for (
         kr = INDEX(line(1:k),')')
         IF( kr + kk - 2 .LE. 0)THEN
            PRINT *,'kk,kr =',kk,kr
            PRINT *,'line =',line(1:132)
#ifdef ONETWO
            CALL STOP('readblock',1) 
#else
        lerrno = 282 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
         ENDIF
         lable = ADJUSTL(line(kk+1:kr-1))
         l=LEN_TRIM(lable)
         IF(swap_lbl)THEN
            kcf = INDEX(lable,',')
            lbl1 = lable(1:kcf-1)
            lbl2 = lable(kcf+1:)
            lable = lbl2(1:LEN_TRIM(lbl2))//','//lbl1(1:LEN_TRIM(lbl1))
         ENDIF
         IF(INDEX(lable,'REACTION_') .GT. 0) THEN !10,11 comes from here
             chrf = lable(10:11)                  ! chrf contains 2 digit number
             WRITE(intfile,FMT='(A)')chrf
             READ(intfile,FMT='(i2)')i2
         ELSEIF (INDEX(lable,'_MCBEAM') .GT. 0) THEN
             kc =0 ;kcf =0;kcb =0
             kcb = INDEX(lable,',',BACK=.TRUE.)
             kcf = INDEX(lable,',')               
             IF(kcf == kcb .AND.  kcb .GT. 0)kc =1    !kc =# of commas in lable
             IF(kcf .NE. kcb)kc =2
             l=LEN_TRIM(lable)
             IF(kc == 0)THEN !lable is just ?_MCBEAM
               i2 = 0
               DO  k =1,SIZE(zflbls)
                IF(zflbls(k)(1:l) == lable(1:l)) THEN
                   i2 = k
                   EXIT
                ENDIF        
               ENDDO
             ELSEIF(kc == 1)THEN
                i1 = 0 ; i2 =0 ;i3 = 0

                 kc = INDEX(lable,',',BACK=.TRUE.)
                 DO  k =1,SIZE(zflbls)
                   IF(zflbls(k)(1:l) == lable(kc+1:l)) THEN
                     i2 = k
                     EXIT
                   ENDIF        
                 ENDDO
                IF(lable(1:2) == 'E1')I1 = 1
                IF(lable(1:2) == 'E2')I1 = 2
                IF(lable(1:2) == 'E3')I1 = 3

                IF(lable(1:1) == "H")i3=1
                IF(lable(1:1) == "D")i3=2
                IF(lable(1:1) == "T")i3=3
                IF(lable(1:3) == "He3")i3=3
                IF(lable(1:3) == "He4")i3=4
                IF(i1 .EQ. 0  .AND. i2 .EQ. 0 .AND. i3 .EQ. 0)THEN
                   PRINT *,'i1,i3,lable =',i1,i3,lable(1:l)
#ifdef ONETWO
                   CALL STOP('i1,i3 problem',2)
#else
        lerrno = 283 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
                ENDIF

             ELSEIF(kc == 2)THEN
                 kc = INDEX(lable,',')      !get loc of first comma
                                            !return thermal species index 
                 IF(lable(1:1) == "H")i3=1
                 IF(lable(1:1) == "D")i3=2
                 IF(lable(1:1) == "T")i3=3
                 IF(lable(1:3) == "He3")i3=3
                 IF(lable(1:3) == "He4")i3=4
                 IF(lable(kc+1:kc+2) == 'E1')i1=1
                 IF(lable(kc+1:kc+2) == 'E2')i1=2
                 IF(lable(kc+1:kc+2) == 'E3')i1=3
                 kc = INDEX(lable,',',BACK=.TRUE.)
                 DO  k =1,SIZE(zflbls)
                   IF(zflbls(k)(1:l) == lable(kc+1:l)) THEN
                     i2 = k
                     EXIT
                   ENDIF        
                 ENDDO
             ELSE
#ifdef ONETWO
                CALL STOP('kc problem',3)
#else
        lerrno = 284 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
             ENDIF
         ELSE IF(LEN_TRIM(lable) ==  2)THEN  
            IF(lable(1:2) == 'E1')i1 =1
            IF(lable(1:2) == 'E2')i1 =2
            IF(lable(1:2) == 'E3')i1 =3
         ELSE IF(lable(1:5) == 'BEAM_')THEN
            i1 =0 ;i4 =0
            chrf1 = lable(6:8)
            WRITE(intfile,FMT='(a)')chrf1
            READ(intfile,FMT='(I3)')i4
            IF(lable(10:11) == 'E1')i1 =1
            IF(lable(10:11) == 'E2')i1 =2
            IF(lable(10:11) == 'E3')i1 =3
#ifdef ONETWO
            IF(i1 == 0)CALL STOP('i1,i4,error',4)
#else
        lerrno = 285 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
         ELSE IF(INDEX(lable,'_TH') .GT. 0)THEN
              kc = INDEX(lable,'_TH')
              DO k=1,SIZE(zthlbls)
                 IF(zthlbls(k)(1:LEN_TRIM(zthlbls(k))) == lable(1:kc-1)) THEN
                   i3 = k
                   EXIT
                 ENDIF
              ENDDO      
         ELSE
            PRINT *,'lable = ',lable
            PRINT *,'sstr =',sstr
!            CALL STOP('readblock',10)
            RETURN ! tried to read input that is not present in file
         ENDIF
      ENDIF
!     nubeam writes zone data for nj-1 zones :
      READ(io,1015,END = 20,ERR= 20)(zbuf(k),k=1,nj-1)  
      CALL convrt_12(zbuf,nj,njm1,TYPE,cvrt)      ! convert units and grid

      aout(1:nj) = zbuf(1:nj)



      !read summary line:
      IF(.NOT. skip_summary)THEN
         READ(io,FMT='(A)',ERR = 20,END=20) line
         kk  = LEN_TRIM(line)  
         k   = INDEX(line,':')           !looking for : in  "integrated value:"
         IF(k .EQ. 0)THEN
            PRINT *,'line=',line(1:kk)
#ifdef ONETWO
            CALL STOP('readblock   read file error 6',6)
#else
        lerrno = 286 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
         ENDIF
         numbr =ADJUSTL(line(k+1:kk))           !remainder of line contains number
         WRITE(intfile,FMT='(A)')numbr
         READ(intfile,FMT='(1pd12.5)') value_intg              !integrated value
         !print *,'value_intg =',value_intg
      ENDIF
      !read blank line at end of data block
      READ(io,FMT='(A)',ERR = 20,END=20) line     !should be a blank line

      !position file for next call to this  routine:
      READ(io,FMT='(A)',ERR = 20,END=20) line
      CALL to_upper_case1(line)
      RETURN

20     PRINT *,'ERROR in readblock, expected blank line'
#ifdef ONETWO
       CALL STOP('readblock',7)
#else
        lerrno = 287 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
1015  FORMAT(5(1x,1pd12.5))
       END  SUBROUTINE readblock






      SUBROUTINE read_scalar(sstr,line,io,scal_val)
      IMPLICIT NONE
      CHARACTER(len =*) line,sstr
      CHARACTER(len =12)numbr
      CHARACTER(len = 24)intfile
      REAL (DP),INTENT(OUT) :: scal_val
      INTEGER k,kk,io
      IF(INDEX(line, sstr ) .GT. 0) THEN
         k=INDEX(line,':',BACK =.TRUE.)  ! note a blank follows ':' 
         IF(k .EQ. 0)THEN
           PRINT *,'k,sstr =',k,sstr
#ifdef ONETWO
           CALL STOP('read_scalar',1)
#else
        lerrno = 288 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
         ELSE
           numbr = line(k+2:k+14)       !format in nubeam is 1pd12.5
           WRITE(intfile,FMT='(a)')numbr
           READ(intfile,FMT='(1pd12.5)')scal_val
           !print *,'scal_val =',scal_val
         ENDIF
      ELSE
         PRINT *,'sstr =',sstr
#ifdef ONETWO
         CALL STOP('read_scalar',2)
#else
        lerrno = 289 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
      ENDIF
      !read blank line at end of data block
       READ(io,FMT='(A)',ERR = 20,END=21) line 
      !position file for next call to this  routine:
       READ(io,FMT='(A)',ERR = 20,END=20) line
       CALL to_upper_case1(line)
21       RETURN

20     PRINT *,'ERROR in read_scalar, expected blank line'
#ifdef ONETWO
       CALL STOP('read_scalar',3)
#else
        lerrno = 290 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
      END SUBROUTINE read_scalar



      SUBROUTINE read_scalar1(sstr,line,io,i2,zflbls,scal_val)
      IMPLICIT NONE
      CHARACTER(len =*) line,sstr
      CHARACTER(len =12)numbr
      CHARACTER(len = 24)intfile
      CHARACTER(len= 3) chra
      CHARACTER(len =*),DIMENSION(:)   :: zflbls
      REAL (DP),INTENT(OUT) :: scal_val
      INTEGER kc,k,kk,kl,io,i2
      IF(INDEX(line, sstr ) .GT. 0) THEN
         k=INDEX(line,':',BACK =.TRUE.)  ! note a blank follows ':' 
         IF(k .EQ. 0)THEN
           PRINT *,'k,sstr =',k,sstr
#ifdef ONETWO
           CALL STOP('read_scalar1',1)
#else
        lerrno = 291 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
         ELSE
           numbr = line(k+2:k+14)       !format in nubeam is 1pd12.5
           WRITE(intfile,FMT='(a)')numbr
           READ(intfile,FMT='(1pd12.5)')scal_val
           IF(INDEX(line,'REACTION_') .GT. 0) THEN !
             kk = INDEX(line,'_')
             chra = line(kk+1:kk+3)                  ! chrf contains 2 digit number
             WRITE(intfile,FMT='(A)')chra
             READ(intfile,FMT='(i2)')i2
           ELSE IF(INDEX(line,'_MCBEAM') .GT. 0) THEN
             kk = INDEX(line,'(')
             kc = INDEX(line,')')
             kl = kc -1 -kk 
             line(1:kl) = line(kk+1:kc-1)
             i2 = 0
             DO  k =1,SIZE(zflbls)
               IF(zflbls(k)(1:kl) == line(1:kl)) THEN
                 i2 = k
                 EXIT
               ENDIF        
             ENDDO
           ELSE IF(INDEX(line,'BEAM_') .GT. 0) THEN
             kk = INDEX(line,'_',BACK =.TRUE.)
             chra = line(kk+1:kk+3)   
             WRITE(intfile,FMT='(A)')chra
             READ(intfile,FMT='(i3)')i2
           ELSE
             PRINT *,'read_scalar1 failed'
             PRINT *,'line =',line
#ifdef ONETWO
             CALL STOP('parse error',2)
#else
        lerrno = 292 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
           ENDIF
         ENDIF
      ELSE
         PRINT *,'sstr =',sstr
         RETURN
      ENDIF
      !read blank line at end of data block
       READ(io,FMT='(A)',ERR = 20,END=20) line 
      !position file for next call to this  routine:
       READ(io,FMT='(A)',ERR = 20,END=20) line
       CALL to_upper_case1(line)
       RETURN

20     PRINT *,'ERROR in read_scalar1, expected blank line'
#ifdef ONETWO
       CALL STOP('readblock',3)
#else
        lerrno = 293 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
      END SUBROUTINE read_scalar1



      SUBROUTINE read_scalar2(sstr,line,io,i1,i4,scal_val)
      IMPLICIT NONE
      CHARACTER(len =*) line,sstr
      CHARACTER(len =12)numbr
      CHARACTER(len = 24)intfile
      CHARACTER(len= 3) chra
      REAL (DP),INTENT(OUT) :: scal_val
      INTEGER j,kc,k,kk,ku,io,i1,i4,nbeams
      IF(INDEX(line, sstr ) .GT. 0) THEN
         k=INDEX(line,':',BACK =.TRUE.)  ! note a blank follows ':' 
         IF(k .EQ. 0)THEN
           PRINT *,'k,sstr =',k,sstr
#ifdef ONETWO
           CALL STOP('read_scalar2',1)
#else
        lerrno = 294 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
         ELSE
           numbr = line(k+2:k+14)       !format in nubeam is 1pd12.5
           WRITE(intfile,FMT='(a)')numbr
           READ(intfile,FMT='(1pd12.5)')scal_val
           k = INDEX(line,'(')
           IF(line(k+1:k+2) == 'E1')i1= 1
           IF(line(k+1:k+2) == 'E2')i1= 2
           IF(line(k+1:k+2) == 'E3')i1= 3
           kk = INDEX(line,')')
           kc = INDEX(line,',')
           ku = INDEX(line,'_',BACK=.TRUE.)
           chra = line(ku+1:kk-1)
           WRITE(intfile,FMT='(a)')chra
           READ(intfile,FMT='(i3)')i4
         ENDIF
      ELSE
         PRINT *,'sstr =',sstr
         RETURN
      ENDIF
      !read blank line at end of data block
       READ(io,FMT='(A)',ERR = 20,END=20) line 
      !position file for next call to this  routine:
       READ(io,FMT='(A)',ERR = 20,END=20) line
       CALL to_upper_case1(line)
       RETURN

20     PRINT *,'ERROR in read_scalar, expected blank line'
#ifdef ONETWO
       CALL STOP('read_scalar2',2)
#else
        lerrno = 295 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
      END SUBROUTINE read_scalar2


      SUBROUTINE  convrt_12(zbuf,nj,njm1,TYPE,cvrt) 
#ifdef ONETWO
      USE mesh,ONLY : roa,imesh
#else
      USE grid_class, ONLY: roa,uniform_mesh
#endif

      IMPLICIT NONE
      REAL(DP) zbuf(*)   !size of zbuf must be njm1 +1
      REAL(DP) c,SAVE,SAVE2,cvrt,roa1,roa2
      INTEGER nj,njm1,j
      CHARACTER(len =1) TYPE
#ifndef ONETWO
      INTEGER imesh
      IF(uniform_mesh)THEN
         imesh =0
      ELSE
         imesh = 1
      ENDIF
#endif


      IF(imesh .NE. 0)THEN
        PRINT *,'convrt_12 is programed only for uniform mesh'
#ifdef ONETWO
        CALL STOP('string_util.f90,convrt_12',1)
#else
        lerrno = 296 + iomaxerr
        CALL terminate(lerrno,nlog)
#endif
      ENDIF
      IF(TYPE .EQ. 'A' .OR. TYPE .EQ. 'V')THEN
         !move zone bundaries up by one index to correspond
         !with onetwo rho grid starting at rho =0.0
         DO j= njm1,1, -1
            zbuf(j+1) = zbuf(j)
         ENDDO
         !Etrapolate to rho = 0 with zero gradient condtion
         !zbuf = a + c*r**2, evaluated at r= 0:
         c = (zbuf(3) -zbuf(2))/(roa(3)**2-roa(2)**2)
         zbuf(1)  = zbuf(2) -c*roa(2)**2
      ELSE
         !zbuf is on zone centers, convert to roa grid:
         !FIRST EXTRAPOLATE TO ROA =0
         !zbuf = a + c*r**2, evaluated at r= 0:
         ROA2 = (ROA(2)+ROA(3))*0.5
         ROA1 = (ROA(1)+ROA(2))*0.5
         c = (zbuf(3) -zbuf(2))/(roa2**2-roa1**2)
         SAVE   = zbuf(2) -c*roa1**2
         SAVE2 = ZBUF(NJM1)
         DO j = NJM1,2,-1
            ZBUF(J) = (ZBUF(J) + ZBUF(J-1))*0.5
         ENDDO
         ZBUF(1) = SAVE
         ZBUF(NJM1+1) =  SAVE2  !CAN DO BETTER, CONSERVE INTEGRALS !!
      ENDIF
         zbuf(1:nj) = zbuf(1:nj)*cvrt
      RETURN
      END    SUBROUTINE  convrt_12

    END MODULE string_util
