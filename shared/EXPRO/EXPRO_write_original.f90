!--------------------------------------------------------
! EXPRO_write_original.f90
!
! PURPOSE:
!  Rewrite original quantities to input.profiles
!
!  NOTE:
!   q is kept positive definite.  All information
!   about sign of q is retained in IPCCW and BTCCW.  
!--------------------------------------------------------

subroutine EXPRO_write_original(tag)

  use EXPRO_globals
  use EXPRO_interface

  implicit none

  integer :: ierr

  character(len=*), intent(in) :: tag
  character (len=80) :: line

  
  open(unit=1,file=trim(path)//'input.profiles',status='old')
  open(unit=2,file=trim(path)//'input.profiles.new',status='replace')

  write(2,'(a)') '# '//tag
  write(2,'(a)') '# '
  do 
     read(1,'(a)',iostat=ierr) line
     if (ierr < 0) exit
     if (line(2:4) /= 'rho') then
        write(2,'(a)') line
     else
        exit
     endif
  enddo
  call EXPRO_write(2)
  
  close(1)
  close(2)

end subroutine EXPRO_write_original
