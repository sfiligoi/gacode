!-----------------------------------------------------
! gyro_ave.f90
!
! PURPOSE:
!  Subroutine for construction of standard 
!  gyro-averaging operator.
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

subroutine gyro_ave(rho,a,u,v,g,p_gyro)

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
  integer :: p_gyro
  integer :: p
  integer :: p0
  integer :: m
  !
  complex, intent(inout) :: g(-m_gyro:m_gyro-i_gyro)
  !
  real :: func(-n_x/2:n_x/2-1)
  complex :: g0
  !
  real, external :: BESJ0
  !-----------------------------------------

  p0 = n_x/2

  do p=-p0,p0-1

     ! kx = pi_2*p*a + u
     ! ky = v

     func(p)=BESJ0(rho*sqrt((pi_2*p*a+u)**2+v**2))**p_gyro/n_x

  enddo

  do m=-m_gyro,m_gyro-i_gyro
     g0 = (0.0,0.0)
     do p=-p0,p0-1
        g0 = g0 + z_gyro(m,p)*func(p)
     enddo
     g(m) = g0
  enddo

end subroutine gyro_ave
