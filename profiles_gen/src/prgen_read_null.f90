!--------------------------------------------------------------
! prgen_read_null.f90
!
! PURPOSE:
!  Minimal setup for gfile reading.
!--------------------------------------------------------------

subroutine prgen_read_null

  use prgen_globals

  implicit none

  integer :: i
  
  nx = n_null

  call prgen_allocate

  do i=1,nx
     rho(i) = (i-1)/(nx-1.0)
  enddo
 
end subroutine prgen_read_null
