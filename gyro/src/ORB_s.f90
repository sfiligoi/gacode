!-----------------------------------------------------
! ORB_s.f90
!
! PURPOSE:
!  Compute s(lambda) = Int[0,lambda] tau(lambda)
!  numerically.
!---------------------------------------------

subroutine ORB_s(lambda,s)

  use gyro_banana_private
  use GEO_interface

  !-----------------------------------------------------------
  implicit none
  !
  integer :: i
  !
  real, intent(in) :: lambda
  real, intent(inout) :: s
  real, dimension(:), allocatable :: f
  !
  real :: d_theta
  real :: f0
  real :: theta_bm
  real :: theta_bp
  real :: theta_0
  !
  real, parameter :: zero_tol = 1e-14
  !-----------------------------------------------------------

  allocate(f(n))

  ! Added '=' in test to take care of
  ! calls at the tp boundary.

  if (lambda <= lambda_tp) then

     ! Passing particle
     !
     ! Limits of integration: [-pi,pi]

     d_theta = 2.0*pi/(n-1)

     do i=1,n

        ! The zero-check here will pick out the 
        ! pinch orbit (lambda=ORB_lambda_tp) at 
        ! theta_0=-pi, and set f=0 to maintain
        ! accuracy.

        theta_0 = -pi+(i-1)*d_theta
        call GEO_interp(theta_0)
        f0 = abs(1.0-lambda*GEO_b)
        if (f0 < zero_tol) f0 = 0.0
        f(i) = (1.0-sqrt(f0))*GEO_g_theta/GEO_b

     enddo

     s = d_theta*(0.5*(f(1)+f(n))+sum(f(2:n-1)))/fluxave

  else

     ! Trapped particle
     !
     ! Limits of integration: [theta_bm,theta_bp]

     call gyro_bounce_points(lambda,theta_bp,theta_bm)

     d_theta = (theta_bp-theta_bm)/(n-1)

     do i=1,n
        theta_0 = theta_bm+(i-1)*d_theta
        call GEO_interp(theta_0)
        f(i) = sqrt(abs(1.0-lambda*GEO_b))*GEO_g_theta/GEO_b
     enddo

     s = 1.0-d_theta*(0.5*(f(1)+f(n))+sum(f(2:n-1)))/fluxave

  endif

  deallocate(f)

end subroutine ORB_s
