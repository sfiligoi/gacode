!---------------------------------------------------------------
! GEO_interface.f90
!
! PURPOSE:
!  Routine to write geometry coefficients to file for plotting.
!---------------------------------------------------------------

subroutine GEO_write(datafile,io)

  use GEO_interface

  implicit none

  !-------------------------------------------------------
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !
  integer :: i
  integer :: n_theta
  !
  real, parameter :: pi=3.141592653589793
  !-------------------------------------------------------

  n_theta = size(GEOV_b)

  open(unit=io,file=datafile,status='replace')

  do i=1,n_theta
     write(io,10) GEOV_theta(i),&
          GEOV_b(i),&
          GEOV_dbdt(i),&
          GEOV_dbdt2(i),&
          GEOV_gsin(i),&
          GEOV_bp(i),&
          GEOV_bigr(i),&
          GEOV_grad_r(i),&
          GEOV_jac_r(i),&
          GEOV_captheta(i),&
          GEOV_l_t(i)

  enddo

  close(io)

10 format(20(es11.4,1x))

end subroutine GEO_write
