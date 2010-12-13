!-----------------------------------------------------
! gyro_ave_rot.f90
!
! PURPOSE:
!  Subroutine for construction of gyro-operator
!  needed for computation of momentum flux.
! 
! NOTES:
!
!  Eikonal for H(x-rho):
!   
!    exp[ -i kx rho cos(a) -i ky rho sin(a) ]
!  
!  a = |grad(r)|/L 
!
!  u = k_theta*Gq*Theta
!
!  v = k_theta*Gq 
!
!  rho = v_perp/Omega_c
!---------------------------------------------

subroutine gyro_ave_rot(rho,a,u,v,g)

  use gyro_globals, only : z_gyro, i_gyro, m_gyro, n_x
  use math_constants

  !-----------------------------------------
  implicit none
  !
  real, intent(in) :: a
  real, intent(in) :: u
  real, intent(in) :: v
  real, intent(in) :: rho
  !
  integer :: p
  integer :: p0
  integer :: m
  integer :: ierr
  !
  complex, intent(inout) :: g(-m_gyro:m_gyro-i_gyro)
  !
  real :: x
  real :: bessel(0:2)
  real :: func(-n_x/2:n_x/2-1)
  complex :: g0
  !-----------------------------------------

  p0 = n_x/2

  do p=-p0,p0-1

     x = rho*sqrt((pi_2*p*a+u)**2+v**2)
     call RJBESL(x,0.0,3,bessel,ierr)

     ! func -> -(1/2)*k_x*rho*[ J0(z)+J2(z) ]

     func(p) = -0.5*(pi_2*p*a+u)*rho*(bessel(0)+bessel(2))/n_x

  enddo

  do m=-m_gyro,m_gyro-i_gyro
     g0 = (0.0,0.0)
     do p=-p0,p0-1
        ! Add an additional factor of i here:
        g0 = g0 + z_gyro(m,p)*func(p)*i_c
     enddo
     g(m) = g0
  enddo

end subroutine gyro_ave_rot
