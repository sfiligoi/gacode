!-----------------------------------------------------------
! sigv.f90
!
! PURPOSE:
!  Compute D-T fusion reactivity <sigma*v> in cm^3/s.
!----------------------------------------------------------

real function sigv(ti,type)

  implicit none

  real, intent(in) :: ti
  character(len=*), intent(in) :: type

  real :: c1,c2,c3,c4,c5,c6,c7
  real :: r0,theta,xi,bg,er

  select case (trim(type))

  case ('hively')

     ! L.M. Hively, Nucl. Fusion 17 (1977) 873.

     ! Table 1: (S5)
     r0  = 0.2935
     c1 = -21.377692
     c2 = -25.204054
     c3 = -7.1013427e-2
     c4 = 1.9375451e-4
     c5 = 4.9246592e-6
     c6 = -3.9836572e-8

     ! Eq. (5)
     sigv = exp(c1/ti**r0+c2+ti*(c3+ti*(c4+ti*(c5+ti*c6))))

  case ('bosch')

     ! H.-S. Bosch and G.M. Hale, Nucl. Fusion 32 (1992) 611.

     ! Table VII:
     c1 = 1.17302e-9
     c2 = 1.51361e-2
     c3 = 7.51886e-2
     c4 = 4.60643e-3
     c5 = 1.3500e-2
     c6 = -1.06750e-4
     c7 = 1.36600e-5
     bg = 34.3827
     er = 1.124656e6

     ! Eq. (12) 
     r0    = ti*(c2+ti*(c4+ti*c6))/(1.0+ti*(c3+ti*(c5+ti*c7)))
     theta = ti/(1.0-r0)
     xi    = (bg**2/(4.0*theta))**(1.0/3.0)

     sigv = c1*theta*sqrt(xi/(er*ti**3))*exp(-3.0*xi)

  case default

     print '(a)','ERROR: (sigv) Unknown reactivity type.'
     sigv = 0.0
     stop

  end select

end function sigv

!--------------------------------------------------------
! sivukhin.f90
!
! PURPOSE:
!  Compute a low-accuracy but fast approximation to the 
!  ion-alpha heating fraction.
!
!               x
!            1  /     dy
!    F(x) = --- | -----------
!            x  /  1+y^(3/2)
!               0
!          
!  Here, F is the fraction of the alpha energy transferred
!  to ions (at common temperature Ti) by collisions, and
!
!               x = E_alpha/E_crit  
!
!  Details are given in Stix, Plasma Phys. 14 (1972) 367.
!  The function F is derived from Sivukhin's energy loss
!  equation and so that is the rationale for the name.
!---------------------------------------------------------

real function sivukhin(x)

  use math_constants

  implicit none

  real, intent(in) :: x
  integer :: i
  real :: f,dy,yi

  integer, parameter :: n=12


  if (x > 0.1) then

     if (x > 4.0) then

        ! Large-x asymptotic formula

        f = (2*pi/3)/sin(2*pi/3)-2.0/sqrt(x)+0.5/(x*x)

     else

        ! Numerical integration

        dy = x/(n-1)
        f  = 0.0
        do i=1,n
           yi = (i-1)*dy
           if (i == 1 .or. i == n) then
              f = f+0.5/(1.0+yi**1.5) 
           else
              f = f+1.0/(1.0+yi**1.5) 
           endif
        enddo
        f = f*dy

     endif

     sivukhin = f/x

  else

     ! Small-x asymptotic series

     sivukhin = 1.0-0.4*x**1.5

  endif

end function sivukhin

subroutine rad_alpha(ne,ni,te,ti,s_alpha_he,s_alpha_i,s_alpha_e,frac_ai,e_cross,n,nion)

  use tgyro_globals, only : &
       pi,&
       me,&
       mi,&
       malpha,&
       therm_flag,&
       zi_vec,&
       tgyro_dt_method,&
       e_alpha,&
       k,&
       tgyro_input_fusion_scale

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: nion
  real, intent(in) :: ni(nion,n)
  real, intent(in) :: ne(n)
  real, intent(in) :: ti(nion,n)
  real, intent(in) :: te(n)
  real, intent(inout) :: s_alpha_he(n)
  real, intent(inout) :: s_alpha_i(n)
  real, intent(inout) :: s_alpha_e(n)
  real, intent(inout) :: frac_ai(n)
  real, intent(inout) :: e_cross(n)

  real, external :: sivukhin
  real, external :: sigv

  real :: x_a
  real :: n_d,n_t
  real :: s_alpha
  real, dimension(:), allocatable :: c_a

  integer :: i

  allocate(c_a(n))

  ! Alpha heating coefficients [Stix, Plasma Phys. 14 (1972) 367] 
  ! See in particular Eqs. 15 and 17.
  c_a(:) = 0.0
  do i=1,nion
     if (therm_flag(i) == 1) then
        c_a(:) = c_a(:)+(ni(i,:)/ne(:))*zi_vec(i)**2/(mi(i)/malpha)
     endif
  enddo

  do i=1,n

     e_cross(i) = k*te(i)*(4*sqrt(me/malpha)/(3*sqrt(pi)*c_a(i)))**(-2.0/3.0)
     x_a = e_alpha/e_cross(i)
     frac_ai(i) = sivukhin(x_a)

     if (tgyro_dt_method == 1) then
        ! Assume D and T given by ion 1 and ion 2 (order doesn't matter)
        n_d = ni(1,i)
        n_t = ni(2,i)
     else
        ! Assume ion 1 is DT hybrid.
        n_d = 0.5*ni(1,i)
        n_t = 0.5*ni(1,i)
     endif

     ! Alpha particle source and power 
     !  - Can use 'hively' or 'bosch' formulae.
     !  - sigv in cm^3/s
     s_alpha_he(i) = n_d*n_t*sigv(ti(1,i)/1e3,'bosch') * tgyro_input_fusion_scale

     s_alpha      = s_alpha_he(i)*e_alpha
     s_alpha_i(i) = s_alpha*frac_ai(i)
     s_alpha_e(i) = s_alpha*(1-frac_ai(i))

  enddo

  deallocate(c_a)

end subroutine rad_alpha
