
FUNCTION get_next_io_unit () RESULT (next)
! -------------------------------------------------------------
! find a unit number available for i/o action
! NOTE: this routien is here for gcnmp compatibility
! Onetwo tradiationally uses getioun (unitv, unitreq)(see cray101.f)
! --------------------------------------------------------------
  USE nrtype,    ONLY : I2B

  USE io_gcnmp,   ONLY : niterdb,nout,nlog,ncrt,ioplot

  USE error_handler,  ONLY : lerrno, terminate

  USE MPI_data,         ONLY : myid,mpiierr,master

#if defined (USEMPI)
      USE mpi
#endif

  IMPLICIT NONE
  INTEGER(I2B) :: next   ! the next available unit number 

  INTEGER(I2B), PARAMETER :: min_unit = 100_I2B, max_unit = 999_I2B
  INTEGER(I2B), SAVE      :: last_unit = ioplot  ! start  with unit ioplot+1
  ! (unit 6 is ncrt above)
  INTEGER(I2B)           :: count                ! number of failures
  LOGICAL            :: OPEN                     ! file status

#if defined USEMPI  
  IF(myid .NE. master)THEN
    PRINT *,'myid =',myid , ' is calling get_next_io_unit'
    PRINT *,'only master(process # ',master,' )' 
    PRINT *,' is allowed to call get_next_io_unit'
    CALL MPI_ABORT(MPI_COMM_WORLD,lerrno,mpiierr)
    CALL EXIT(1)
  ENDIF
#endif

  count = 0_I2B ; next = min_unit - 1_I2B
  IF ( last_unit > 0_I2B ) THEN ! check next in line
     next = last_unit + 1_I2B
     INQUIRE (unit=next, opened=OPEN)
     IF ( .NOT. OPEN ) last_unit = next ! found it
     RETURN
  ELSE ! loop through allowed units
     DO ! forever
        next = next + 1_I2B
        INQUIRE (unit=next, opened=OPEN)
        IF ( .NOT. OPEN ) THEN 
           last_unit = next     ! found it
           EXIT ! the unit loop
        END IF
        IF ( next == max_unit ) THEN ! attempt reset 3 times
           last_unit = 0
           count     = count + 1_I2B
           IF ( count <= 3 ) next = min_unit - 1_I2B
        END IF ! reset try
        IF ( next > max_unit ) THEN ! abort
           lerrno =6
           CALL terminate(lerrno,nlog)
        END IF ! abort
     END DO ! over unit numbers
  END IF ! last_unit
END FUNCTION get_next_io_unit



