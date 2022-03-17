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

  call prgen_allocate
  allocate(dummy(nx))

  open(unit=1,file=file_state,status='old')
  read(1,*)
  do i=1,nx
     read(1,*) &
          dummy(i),&
          dpsi(i),&
          dummy(i),&
          p_tot(i),&
          q(i),&
          te_kev(i),&
          ti_kev(i),&
          ne_e19m3(i),&
          ni_e19m3(i)
  enddo

  dpsi = (-dpsi+dpsi(1))/(2*pi)
 
  ! Units conversions
  te_kev = te_kev*1e-3
  ti_kev = ti_kev*1e-3
  
  ne_e19m3 = ne_e19m3*1e-19
  ni_e19m3 = ni_e19m3*1e-19

  deallocate(dummy)
  
end subroutine prgen_read_genf
