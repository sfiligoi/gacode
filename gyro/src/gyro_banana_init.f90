!-----------------------------------------------------
! gyro_banana_init.f90
!
! PURPOSE:
!  Simple routine to initialize internal variables.
!-----------------------------------------------------

subroutine gyro_banana_init(n_IN)

  use gyro_banana_private
  use GEO_interface

  !-----------------------------------------------------------
  implicit none
  !
  integer :: i
  integer, intent(in) :: n_IN 
  !
  real, dimension(:), allocatable :: f
  real :: d_theta
  real :: theta_0
  !-----------------------------------------------------------

  n = n_IN

  allocate(f(n))

  !----------------------------------------------
  ! Determine maximum lambda, and lambda at
  ! the trapped-passing boundary:
  !
  call GEO_interp(0.0)
  !
  lambda_max = 1.0/GEO_b
  !
  call GEO_interp(pi)
  !
  lambda_tp = 1.0/GEO_b
  !-----------------------------------------------

  d_theta = 2.0*pi/(n-1)

  do i=1,n
     theta_0 = -pi+(i-1)*d_theta
     call GEO_interp(theta_0)
     f(i) = GEO_g_theta/GEO_b
  enddo

  fluxave = (0.5*(f(1)+f(n))+sum(f(2:n-1)))*d_theta

  deallocate(f)

end subroutine gyro_banana_init
