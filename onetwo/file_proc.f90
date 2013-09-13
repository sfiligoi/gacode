 MODULE file_proc
    INTEGER,PARAMETER :: string_length = 256
    CHARACTER (len = string_length) int_io,new_name


    CONTAINS 

    SUBROUTINE append_time_to_filename (time,filename)
!----------------------------------------------------------
! -- loads new_name
!----------------------------------------------------------
    USE nrtype,                             ONLY : DP
    USE echdat_module,                      ONLY : model_globl
    IMPLICIT NONE

    REAL(DP) time
    CHARACTER (len=*) filename
    CHARACTER *4 id
    INTEGER j

       IF(LEN_TRIM(filename) .GT. string_length)THEN
          PRINT *,'error in file_proc, string length too short'
          PRINT *,'string length = ',string_length
          PRINT *,'require at least :', LEN_TRIM(filename)
          CALL STOP ('subroutine append_time_to_filename :', 1)
       ENDIF


       !find location of .nc or .cdf  (if netcdf file)
       new_name = ADJUSTL(filename)

       j = 0
       j = INDEX(new_name,'.nc',BACK = .TRUE. )
       IF( j .NE. 0)id ='.nc '
       IF(j == 0)THEN
          j = INDEX(new_name,'.cdf',BACK = .TRUE. )
          id = '.cdf'
       ENDIF
       IF(j ==0)THEN
          j = INDEX(new_name,'.txt',BACK = .TRUE. )
          id = '.txt'
       ENDIF
       IF(j == 0)THEN
          j = LEN_TRIM(new_name)
          id =''
       ENDIF
 
!       WRITE(int_io,'(A,A,1pe12.6,A,i1,A,A)')new_name(1:j-1),'_',time,'_',model_globl,'_',id
       WRITE(int_io,'(A,A,1pe12.6,A)')new_name(1:j-1),'_',time,id
       READ (int_io,'(A)')new_name
 

    RETURN
    END     SUBROUTINE append_time_to_filename


   SUBROUTINE append_gyro_id_to_filename (time,filename)
!----------------------------------------------------------
! -- loads new_name
!----------------------------------------------------------
    USE nrtype,                             ONLY : DP
    USE echdat_module,                      ONLY : model_globl
    IMPLICIT NONE

    REAL(DP) time
    CHARACTER (len=*) filename
    CHARACTER *4 id
    INTEGER j

       IF(LEN_TRIM(filename) .GT. string_length)THEN
          PRINT *,'error in file_proc, string length too short'
          PRINT *,'string length = ',string_length
          PRINT *,'require at least :', LEN_TRIM(filename)
          CALL STOP ('subroutine append_time_to_filename :', 1)
       ENDIF


       !find location of .nc or .cdf  (if netcdf file)
       new_name = ADJUSTL(filename)
 
       j = 0
       j = INDEX(new_name,'.nc',BACK = .TRUE. )
       IF( j .NE. 0)id ='.nc '
       IF(j == 0)THEN
          j = INDEX(new_name,'.cdf',BACK = .TRUE. )
          id = '.cdf'
       ENDIF
       IF(j ==0)THEN
          j = INDEX(new_name,'.txt',BACK = .TRUE. )
          id = '.txt'
       ENDIF
       IF(j == 0)THEN
          j = LEN_TRIM(new_name)
          id =''
       ENDIF
 
       WRITE(int_io,'(A,A,1pe12.6,A,i1,A,A)')new_name(1:j-1),'_',time,'_',model_globl,'_',id
 
       READ (int_io,'(A)')new_name

 
    RETURN
    END     SUBROUTINE append_gyro_id_to_filename

    SUBROUTINE rename_disk_file(orig_name)
!----------------------------------------------------------------------
! -- rename file in current working directory fromorig_name to new_name
! -- Note that new_name is assumed loaded in this module by a previous
! -- call to append_time_to_filename
!-----------------------------------------------------------------------

     CHARACTER(len=*) orig_name

        call system("mv " // trim(orig_name) // " " // trim(new_name))

     RETURN
     END SUBROUTINE rename_disk_file

    SUBROUTINE  temp_file_test
!---------------------------------------------------------------------
! --- just a tester to  rename disk files
! --- NOT CALLED IN NORMAL USE OF CODE
!-------------------------------------------------------HSJ------------
      USE nrtype,                                      ONLY : DP,I4B
      USE io,                                          ONLY : prtlst
      USE solcon,                                      ONLY : time,time0,timmax 
      IMPLICIT NONE
      INTEGER j 
      REAL(DP) ltime
      CHARACTER(len=256) filename

      filename ='toray1.nc'
      DO j=1,SIZE(prtlst)
         if(prtlst(j) .lt. time0)CYCLE
         IF(prtlst(j) .gt. timmax)EXIT
         CALL SYSTEM('touch ' // trim(filename) )
         ltime = prtlst(j)
         CALL append_time_to_filename (ltime,filename)
         CALL rename_disk_file(filename)
      ENDDO
      RETURN
    END SUBROUTINE   temp_file_test




      SUBROUTINE delete_file(filename)
! ----------------------------------------------------------------------
! unlink not standard ?? So use this instead - HSJ
! _____________________________________________________________________

        USE nrtype,                      ONLY : I4B
        USE io,                          ONLY : io_temp
        IMPLICIT NONE
        CHARACTER(LEN = *) filename
        LOGICAL exists,opnd
        INTEGER(I4B) iounit                         


        INQUIRE(FILE=filename,EXIST = exists, OPENED = opnd, &
                                           NUMBER = iounit)
        IF(exists)THEN
           IF( .NOT. opnd)THEN
               !CALL getioun(io_temp) does nothing
               iounit= io_temp
               OPEN(UNIT=iounit,FILE=filename)
           ENDIF
           CLOSE(UNIT=iounit,STATUS = "DELETE")
        ENDIF
        RETURN
      END SUBROUTINE delete_file


 END MODULE file_proc
