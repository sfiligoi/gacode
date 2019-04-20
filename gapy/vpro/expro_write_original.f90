!--------------------------------------------------------
! expro_write_original.f90
!
! PURPOSE:
!  Rewrite original quantities to input.profiles
!
!  NOTE:
!   q is kept positive definite.  All information
!   about sign of q is retained in IPCCW and BTCCW.  
!--------------------------------------------------------

subroutine expro_write_original(datafile1,datafile2,tag)

  use vpro

  implicit none

  character (len=*), intent(in) :: datafile1,datafile2
  character (len=*), intent(in) :: tag

  integer :: ierr
  character (len=80) :: line

  open(unit=2,file=trim(datafile1),status='old')
  open(unit=1,file=trim(datafile2),status='replace')

  write(1,'(a)') '# '//tag
  write(1,'(a)') '# '
  do 
     read(2,'(a)',iostat=ierr) line
     if (ierr < 0) exit
     if (line(2:4) /= 'rho') then
        write(1,'(a)') line
     else
        exit
     endif
  enddo
  close(2)
  
  ! NOTE: vpro will use and close unit 1
  call vpro_write
  
end subroutine expro_write_original
