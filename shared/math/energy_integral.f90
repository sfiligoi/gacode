!---------------------------------------------------------
! energy_integral.f90
!
! PURPOSE:
!  Construct energy grid and integration weights.
!
! NOTES:
!  For one energy gridpoint, the rule is trivial.
!
!  For n_energy > 1, we break the integration into
!  two regions: 
!
!  region 1: (0,energy_max)
!  region 2: (energy_max,inf)
!
!  Region 1:
!  --------
!
!  To do the region 1 integral, we change variable to s, 
!  where:
!                          x
!                   2      /          -t
!        s(x) = --------   | sqrt(t) e   dt
!               sqrt(pi)   /
!                          0
!
!  Clearly, s(inf) = 2/sqrt(pi) * Gamma(3/2) = 1.0
!
!  s(x) = P(3/2,x)  [p. 262, Abramowitz and Stegun]
!
!  P(3/2,x) = P(1/2,x) - sqrt(x) exp(-x) / Gamma(3/2)
!  P(1/2,x) = Erf(sqrt(x))
!
!  We use Gauss-Legendre weights on (0,s(energy_max)) to 
!  generate n_energy-1 points:
!
!     "call GaussLegendre(0.0,s0,sn,wn,n_energy-1)"
! 
!  Then, we invert s(energy(i)) = sn(i) for the 
!  abscissae energy(i).
!
!
!  Region 2:
!  --------
!
!  Method 1:
!
!   The region 2 integral (which ought to be exponentially
!   small) is done with one weight at position 
!
!            energy(n_energy) = energy_max.
!
!   The weight is just 1-s0.
!---------------------------------------------------------

subroutine energy_integral(n_energy,energy_max,energy,w_energy)

  use math_constants

  !------------------------------------------
  implicit none
  !
  integer, intent(in) :: n_energy
  real, intent(in) :: energy_max
  ! 
  real, intent(inout) :: energy(n_energy)
  real, intent(inout) :: w_energy(n_energy)
  !
  integer :: i
  !
  real :: s_max
  real :: s1
  real :: d_e
  !
  real, dimension(:), allocatable :: sn
  real, dimension(:), allocatable :: wn
  !
  real, external :: p32
  !------------------------------------------

  allocate(sn(n_energy-1))
  allocate(wn(n_energy-1))

  d_e = energy_max/n_energy

  if (n_energy == 1) then

     energy(1)   = energy_max
     w_energy(1) = 1.0

  else 

     ! p32(x) == P(3/2,x)

     s_max = p32(energy_max)
     s1 = 1.0-s_max

     call gauss_legendre(0.0,s_max,sn,wn,n_energy-1)

     ! Map abscissae to energy:

     do i=1,n_energy-1
        w_energy(i) = wn(i)
        call invert_p32(sn(i),energy(i),energy_max)
     enddo

     ! Add remainder

     energy(n_energy) = energy_max

     w_energy(n_energy) = s1 

  endif

  deallocate(sn)
  deallocate(wn)

end subroutine energy_integral
