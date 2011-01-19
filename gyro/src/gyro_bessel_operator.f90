!--------------------------------------------------------------
! gyro_bessel_operator.f90
!
! PURPOSE:
!  Compute all the Bessel-function kernel (averaging) stencils
!  to the required bandwidth (m_gyro).  In the limit of maximal
!  m_gyro (n_x/2), these are the pseudospectral operators.
!
! NOTES:
!
! itype=1:  J_0(x) 
!       2:  J_0^2(x)
!       3:  -(i/2)*k_x*rho*[ J_0(x)+J_2(x) ]
!       4:  G_perp = (1/2)*[ J_0(x)+J_2(x) ] 
!       5:  G_perp^2(x) 
!       6:  G_perp(x)*J_0(x)
!       7:  I_0(x^2)
!
! x = rho*sqrt((2*pi*p*a+u)^2+v^2)
!
! kx -> pi_2*p*a + u
! ky -> v
! 
! INPUT:
!
!  a = |grad(r)|/L 
!  u = k_theta*Gq*Theta
!  v = k_theta*Gq 
!
!  rho = v_perp/Omega_c               (itype=1,..,6)
!  rho = rhos_unit*sqrt(T)/(mu*z*b0)  (itype=7)
!---------------------------------------------------------------

subroutine gyro_bessel_operator(rho,a,u,v,g,itype)

  use gyro_globals, only : z_gyro, i_gyro, m_gyro, n_x
  use math_constants

  !-----------------------------------------
  implicit none
  !
  real, intent(in) :: a
  real, intent(in) :: u
  real, intent(in) :: v
  real, intent(in) :: rho
  integer, intent(in) :: itype
  complex, intent(inout) :: g(-m_gyro:m_gyro-i_gyro)
  !
  integer :: p
  integer :: p0
  integer :: m
  integer :: ierr
  !
  real :: x
  real :: bessel(0:2)
  real :: func(-n_x/2:n_x/2-1)
  complex :: g0
  !
  real, external :: BESJ0
  real, external :: BESEI0
  !-----------------------------------------

  p0 = n_x/2


  select case (itype)

  case (1)

     ! J_0

     do p=-p0,p0-1
        x = rho*sqrt((pi_2*p*a+u)**2+v**2)
        func(p)=BESJ0(x)/n_x
     enddo

  case (2)

     ! J_0^2

     do p=-p0,p0-1
        x = rho*sqrt((pi_2*p*a+u)**2+v**2)
        func(p)=BESJ0(x)**2/n_x
     enddo

  case (3)

     ! POSSIBLE ISSUE: the accuracy of this operator when 
     ! banded has not been studied.

     ! -(i/2)*k_x*rho*[ J0(z)+J2(z) ]

     ! The factor -(i/2) will be applied outside this loop

     do p=-p0,p0-1
        x = rho*sqrt((pi_2*p*a+u)**2+v**2)
        call RJBESL(x,0.0,3,bessel,ierr)
        func(p) = (pi_2*p*a+u)*rho*(bessel(0)+bessel(2))/n_x
     enddo

  case (4)

     ! G = (1/2)*[ J_0(z)+J_2(z) ] 

     do p=-p0,p0-1
        x = rho*sqrt((pi_2*p*a+u)**2+v**2)
        call RJBESL(x,0.0,3,bessel,ierr)
        func(p) = 0.5*(bessel(0)+bessel(2))/n_x
     enddo

  case (5)

     ! G^2 

     do p=-p0,p0-1
        x = rho*sqrt((pi_2*p*a+u)**2+v**2)
        call RJBESL(x,0.0,3,bessel,ierr)
        func(p) = (0.5*(bessel(0)+bessel(2)))**2/n_x
     enddo

  case (6)

     ! G * J_0

     do p=-p0,p0-1
        x = rho*sqrt((pi_2*p*a+u)**2+v**2)
        call RJBESL(x,0.0,3,bessel,ierr)
        func(p) = 0.5*(bessel(0)+bessel(2))*bessel(0)/n_x
     enddo

  case (7)

     ! I_0

     do p=-p0,p0-1
        x = rho*sqrt((pi_2*p*a+u)**2+v**2)
        func(p) = BESEI0(x*x)/n_x
     enddo

  end select


  do m=-m_gyro,m_gyro-i_gyro
     g0 = (0.0,0.0)
     do p=-p0,p0-1
        g0 = g0 + z_gyro(m,p)*func(p)
     enddo
     g(m) = g0
  enddo

  ! Add final factor of i for case 3
  if (itype == 3) g = -(i_c/2.0)*g

end subroutine gyro_bessel_operator
