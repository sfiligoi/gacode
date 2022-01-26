!--------------------------------------------------------------
! prgen_read_genf.f90
!
! PURPOSE:
!  Read profiles from General Fusion Statefile
!--------------------------------------------------------------

subroutine prgen_read_genf

  use prgen_globals

  implicit none

  integer :: ios,i
  character(len=16) :: a
  real, dimension(:), allocatable :: dummy

  open(unit=1,file=file_state,status='old')
  nx = -1
  do
     read(1,*,iostat=ios) a
     if (ios < 0) exit
     nx = nx+1
  enddo
  close(1)

  allocate(dummy(nx))

  open(unit=1,file=file_state,status='old')
  read(1,*)
  do i=1,nx
     read(1,*) &
          dpsi(i),&
          dummy(i),&
          dummy(i),&
          p_tot(i),&
          q(i),&
          te_kev(i),&
          ti_kev(i)
  enddo

  call prgen_allocate

end subroutine prgen_read_genf
