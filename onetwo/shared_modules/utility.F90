MODULE utility
! The purpose of this module is to collect subroutines and functions, which may
! be of general utility.  These subprograms should not depend on other modules,
! except perhaps nrtype.  To make any routine MPI compatible, then they should
! use preprocessor statements of the form 
! #if defined (USEMPI)
!   ...
! #endif
! 
CONTAINS
  FUNCTION get_executable_name(name,return_size,ierr)
    !
    ! NAME: 
    !   get_executable_name
    ! PURPOSE: 
    !   Get a valid path to the executable given by name.  This is
    !     accomplished by first testing to see if name in its given form
    !     exists.  If it exists then it is returned.  If it doesn't
    !     exist, then the shell command  `which` is run, writing the
    !     results to the file .program_name, which is then read and
    !     destroyed.

    
    ! ARGUMENTS:
    !   name - A string giving the executable name for which a valid
    !     path will  be found
    !   return_size - The size of the character string which will
    !     receive the return value of this function.  On return, it
    !     holds the length of the executable string (thus giving the
    !     necessary length of the return string)
    !   ierr - The error return code variable.  It can take the 
    !       following values
    !     0 - No error
    !     1 - The return_size is not large enough to hold the valid
    !       path  to the executable.  See the return value of
    !       return_size to see how large it needed to be. 
    !     2 - A valid path could not be found.
    ! INTERNAL FILE UNITS:
    !   For communication with the 'which' command output file, the unit
    !     10 is  used.
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, INTENT(INOUT) :: return_size
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: name_len, io_unit = 10
    CHARACTER(len=return_size) :: get_executable_name
    CHARACTER(len=*), PARAMETER :: out_name = '.program_name'
    CHARACTER(len=5000) :: long_name = ' '
    CHARACTER(len=256) :: command
    LOGICAL :: ex
    get_executable_name = ' '
    INQUIRE(file=name,exist=ex)
    !WRITE (*,*) 'Looking for executable',name(1:len_trim(name)),ex
    IF (.not.ex) THEN
      command = 'which '//name(1:LEN_TRIM(name))//' > '//out_name
      CALL SYSTEM(command)
      OPEN (file=out_name,unit=io_unit,status='old')
      READ (io_unit,fmt='(a)') long_name
      CLOSE (io_unit)
      command = 'rm '//out_name
      CALL SYSTEM(command)
    ELSE
      long_name(1:LEN_TRIM(name)) = name(1:LEN_TRIM(name))
      long_name(LEN_TRIM(name)+1:LEN(long_name)) = ' '
    ENDIF
    name_len = LEN_TRIM(long_name)
    !write (*,*) name_len,return_size
    
    IF (name_len <= return_size) THEN
      get_executable_name(1:name_len) = long_name(1:name_len)
      ierr = 0
      INQUIRE(file=get_executable_name,exist=ex)
      if (.not.ex) THEN
        ierr = 2
      ENDIF
    ELSE
      get_executable_name(1:return_size) = long_name(1:return_size)
      ierr = 1
    ENDIF
    return_size = name_len
    !write(*,*) get_executable_name,ierr
  END FUNCTION get_executable_name

  FUNCTION get_executable_creation_time(ierr)
    IMPLICIT NONE
    CHARACTER*24 get_executable_creation_time
    INTEGER*8 :: time
    INTEGER :: input_args, getarg
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER*24 :: ctime
    CHARACTER*256 :: name
    INTEGER*4 :: lstat, statb(13), name_len
    input_args = getarg(0,name)
    name_len = len(name)
    name = get_executable_name(name,name_len,ierr)
    if (ierr .ne. 0) return 
    ierr = lstat ( name, statb )
    if ( ierr .ne. 0 ) return
    get_executable_creation_time = ctime(statb(9))
    !write (*,*) get_executable_creation_time
  END FUNCTION get_executable_creation_time
END MODULE
