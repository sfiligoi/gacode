!-----------------------------------------------------
! ORB_do.f90
!
! PURPOSE:
!  Find iterative solution for 
!
!  d(theta) = d(tau) sqrt[ 1-lambda B(theta) ]
!-----------------------------------------------------

subroutine gyro_banana_uniform_taugrid(lambda,n_tau,n_sub,theta,tau)

  use gyro_banana_private
  use GEO_interface

  !-----------------------------------------------------------
  implicit none
  !
  integer, intent(in) :: n_tau
  integer, intent(in) :: n_sub
  integer :: i
  integer :: j
  integer :: ni
  integer, parameter :: n_iterate=20
  !
  real, dimension(n_tau) :: theta
  real, dimension(n_tau) :: tau
  !
  real, intent(in) :: lambda
  real, dimension(:), allocatable :: f
  real, dimension(:), allocatable :: x
  real :: fsum
  real :: dt
  real :: theta_bm
  real :: theta_bp
  !-----------------------------------------------------------


  if (lambda < lambda_tp) then

     ! Passing particles

     theta_bm = -pi
     theta_bp = pi

  else

     ! Trapped particles

     call gyro_bounce_points(lambda,theta_bp,theta_bm)

  endif

  !----------------------------------------
  ! Compute number of integration intervals
  !
  ni = 1+(n_tau-1)*(n_sub+1)
  !
  allocate(x(ni))
  allocate(f(ni))
  !----------------------------------------

  do i=1,ni
     x(i) = theta_bm+(i-1)*(theta_bp-theta_bm)/(ni-1)
  enddo

  do j=1,n_iterate

     do i=1,ni
        call GEO_interp(x(i))
        f(i) = sqrt(abs(1.0-lambda*GEO_b))/GEO_g_theta
     enddo

     fsum = 0.5*(f(1)+f(ni))+sum(f(2:ni-1))
     dt = (theta_bp-theta_bm)/fsum

     x(1) = theta_bm
     do i=2,ni
        x(i) = x(i-1)+0.5*dt*(f(i)+f(i-1))
     enddo

  enddo ! j

  do j=1,n_tau
     i = 1+(j-1)*(n_sub+1)
     theta(j) = x(i)
     tau(j)   = (i-1)*dt
  enddo

  deallocate(x)
  deallocate(f)

end subroutine gyro_banana_uniform_taugrid
