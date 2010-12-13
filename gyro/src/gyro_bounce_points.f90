!-----------------------------------------------------
! gyro_bounce_points.f90
!
! PURPOSE:
!  For a given lambda, return the upper and lower
!  bounce points theta_bp,theta_bm.  These solve 
!  
!         B(theta_b) = 1/lambda 
! 
!  If there are no bounce points, return 
!  theta_b=pi or theta_b=0.0.
!
! NOTES:
!  This routine is valid for arbitrary flux
!  surface shape.  Also, speed is not essential, 
!  so determine roots via bisection for stability.
!----------------------------------------------------

subroutine gyro_bounce_points(lambda,theta_bp,theta_bm)

  use math_constants
  use GEO_interface

  !---------------------------------
  implicit none
  !
  real, parameter :: eps_theta = 1e-15
  !
  real, intent(in)    :: lambda
  real, intent(inout) :: theta_bp
  real, intent(inout) :: theta_bm
  !
  real :: lambda_tp
  real :: lambda_max
  real :: theta_b
  !
  real :: b0_1
  real :: b0_2
  real :: residual
  real :: x
  real :: dtheta
  !  
  real :: b0
  !---------------------------------
  
  theta_b = 0.0
  call GEO_interp(theta_b) ; b0 = GEO_b
  lambda_max = 1.0/b0

  ! Test for lambda too large:
  if (lambda >= lambda_max) then
     print *,'lambda >= lambda_max'
     return
  endif

  theta_b = pi
  call GEO_interp(theta_b) ; b0 = GEO_b
  lambda_tp = 1.0/b0

  ! Test for lambda too small:
  if (lambda <= lambda_tp) then
     print *,'lambda <= lambda_tp'
     return
  endif

  ! Upper bounce point

  residual = 1.0

  dtheta  = pi/2.0-eps_theta
  theta_b = pi/2.0

  x = 1.0/lambda

  do while (residual > eps_theta)

     call GEO_interp(theta_b) ; b0_1 = GEO_b
     call GEO_interp(theta_b+dtheta) ; b0_2 = GEO_b

     dtheta = 0.5*dtheta

     if (b0_1 <= x .and. b0_2 >= x) then
        theta_b = theta_b+dtheta
     else
        theta_b = theta_b-dtheta
     endif

     residual = abs(b0_1-b0_2)

  enddo

  theta_bp = theta_b

  ! Lower bounce point

  residual = 1.0

  dtheta  = -pi/2.0+eps_theta
  theta_b = -pi/2.0

  x = 1.0/lambda

  do while (residual > eps_theta)

     call GEO_interp(theta_b) ; b0_1 = GEO_b
     call GEO_interp(theta_b+dtheta) ; b0_2 = GEO_b

     dtheta = 0.5*dtheta

     if (b0_1 <= x .and. b0_2 >= x) then
        theta_b = theta_b+dtheta
     else
        theta_b = theta_b-dtheta
     endif

     residual = abs(b0_1-b0_2)

  enddo

  theta_bm = theta_b

end subroutine gyro_bounce_points
