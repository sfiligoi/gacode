!-----------------------------------------------------
! gyro_ave_aperp.f90
!
! PURPOSE:
!  Subroutine for construction of gyro-operator
!  needed for A_perp field equation.
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

subroutine gyro_ave_aperp(rho,a,u,v,g,g_type)

  use gyro_globals, only : z_gyro, i_gyro, m_gyro, n_x
  use math_constants

  !-----------------------------------------
  implicit none
  !
  real, intent(in) :: a
  real, intent(in) :: u
  real, intent(in) :: v
  real, intent(in) :: rho
  integer, intent(in) :: g_type
  !
  integer :: p
  integer :: p0
  integer :: m
  integer :: ierr
  !
  complex, intent(inout) :: g(-m_gyro:m_gyro-i_gyro)
  !
  real :: x, y
  real :: bessel(0:2)
  real :: func(-n_x/2:n_x/2-1)
  complex :: g0
  !-----------------------------------------

  p0 = n_x/2

  do p=-p0,p0-1

     x = rho*sqrt((pi_2*p*a+u)**2+v**2)
     call RJBESL(x,0.0,3,bessel,ierr)

     ! G_perp -> (1/2)*[ J0(z)+J2(z) ]
     y = 0.5*(bessel(0)+bessel(2))
     if(g_type == 2) then
        ! G_perp^2
        y = y**2
     else if(g_type == 3) then
        ! G_perp * G
        y = y * bessel(0)
     endif

     func(p) = y/n_x

  enddo

  do m=-m_gyro,m_gyro-i_gyro
     g0 = (0.0,0.0)
     do p=-p0,p0-1
        g0 = g0 + z_gyro(m,p)*func(p)
     enddo
     g(m) = g0
  enddo

end subroutine gyro_ave_aperp
