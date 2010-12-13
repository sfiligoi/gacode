!-----------------------------------------------------
! gyro_bdoubleave.f90
!
! PURPOSE:
!  Subroutine for construction of v-integrated 
!  double gyroaverage.
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
!  rho = rhos_unit*sqrt(T)/(mu*z*b0)
!---------------------------------------------

subroutine gyro_bdoubleave(rho,a,u,v,g)

  use gyro_globals, only : z_gyro, i_gyro, m_gyro , n_x
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
  !
  complex, intent(inout) :: g(-m_gyro:m_gyro-i_gyro)
  !
  real :: rho2
  real :: func(-n_x/2:n_x/2-1)
  complex :: g0
  !
  real, external :: BESEI0
  !-----------------------------------------

  p0 = n_x/2

  rho2 = rho*rho

  do p=-p0,p0-1
     func(p) = BESEI0(rho2*((pi_2*p*a+u)**2+v**2))/n_x
  enddo

  do m=-m_gyro,m_gyro-i_gyro
     g0 = (0.0,0.0)
     do p=-p0,p0-1
        g0 = g0 + z_gyro(m,p)*func(p)
     end do
     g(m) = g0
  end do

end subroutine gyro_bdoubleave
