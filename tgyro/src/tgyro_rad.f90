!---------------------------------------------------------
! Portable subroutines for calculation of radiation:
!
! rad_sync -> synchrotron 
! rad_ion  -> Bremsstrahlung and impurity line 
!---------------------------------------------------------

subroutine rad_sync(aspect_rat,r_min,b_ref,ne,te,qsync,n)

  ! IN:
  !
  ! aspect_rat = R/a [-]
  ! r_min      = a   [cm]
  ! b_ref      = B   [G]
  !
  ! ne [1/cm^3]
  ! te [keV]
  !
  ! OUT:
  !
  ! qsync [erg/cm^3/s]

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: aspect_rat,r_min
  real, intent(in), dimension(n) :: b_ref,ne,te
  real, intent(inout), dimension(n) :: qsync

  integer :: i
  real :: pi,e,k,me,c
  real :: g,phi,wpe,wce

  ! Reflection coefficient (Rosenbluth)
  real, parameter :: r_coeff=0.8

  !------------------------------------------------------
  ! PHYSICAL CONSTANTS
  !
  pi = 4.0*atan(1.0)
  e  = 4.8032e-10 ! statcoul
  k  = 1.6022e-12 ! erg/eV
  me = 9.1094e-28 ! g
  c  = 2.9979e10  ! cm/s
  !------------------------------------------------------

  !-------------------------------------------------------
  ! Synchrotron radiation
  ! - Trubnikov, JETP Lett. 16 (1972) 25.
  do i=1,n   
     wpe = sqrt(4*pi*ne(i)*e**2/me)
     wce = e*abs(b_ref(i))/(me*c)
     g   = k*te(i)/(me*c**2)
     phi = 60*g**1.5*sqrt((1.0-r_coeff)*(1+1/aspect_rat/sqrt(g))/(r_min*wpe**2/c/wce))

     qsync(i) = me/(3*pi*c)*g*(wpe*wce)**2*phi
  enddo
  !-------------------------------------------------------

end subroutine rad_sync

subroutine rad_ion(&
     te,&
     ne,&
     ni,&
     z,&
     name,&
     qbrem,&
     qline,&
     nion,&
     n)

  ! IN:
  !
  ! te    [keV]
  ! ne    [1/cm^3]
  ! ni    [1/cm^3]
  ! z (charge) [-]
  ! name (ion names) [string] 
  ! nion (number of ions)
  ! n    (number of radii)
  !
  ! OUT:
  !
  ! qbrem  [erg/cm^3/s]
  ! qline  [erg/cm^3/s]

  implicit none

  integer, intent(in) :: nion,n
  real, intent(in), dimension(n) :: te,ne
  real, intent(in), dimension(nion,n) :: ni
  real, intent(in), dimension(nion) :: z
  character(len=3), intent(in), dimension(nion) :: name
  real, intent(inout), dimension(n) :: qbrem,qline

  ! Internals

  integer :: i
  real, dimension(n) :: lz
  real, dimension(nion,n) :: qbremi,qpost

  qbrem = 0.0
  qline = 0.0

  do i=1,nion

     ! Bremsstrahlung radiation [erg/cm^3/s]
     ! - From NRL formulary
     ! - 1 W/cm^3 = 1e7 erg/cm^3/s

     qbremi(i,:) = 1e7*1.69e-32*ne*ni(i,:)*z(i)**2*sqrt(te)

     if (z(i) > 1.0) then
        ! lz is Post 1977 cooling rate in erg cm^3/s
        call post77(te/1e3,name(i),lz,n)
     else
        qpost(i,:) = 0.0
     endif
     qpost(i,:) = ne*ni(i,:)*lz

  enddo

  ! Compute totals by summing over ions
  !
  ! NOTE: Post appears (?) to give sum of Bremsstrahlung and Line
  !       (thus requires subtraction).
  do i=1,nion
     qbrem = qbrem+qbremi(i,:)
     if (z(i) > 1.0) then
        qline = qline+qpost(i,:)-qbremi(i,:)          
     endif
  enddo
 
end subroutine rad_ion

subroutine post77(te,name,lz,n)

  implicit none

  integer, intent(in) :: n

  character(len=3), intent(in) :: name
  real, intent(in), dimension(n) :: te
  real, intent(inout), dimension(n) :: lz
  real, dimension(5,8) :: c
  real, dimension(6) :: a

  integer :: i,j,success
  real :: x,t0,s

  ! Set coefficients

  c = 0.0

  if (name == 'He4') then
     ! Helium 4 (A=4,Z=2) [Agrees with Post for Te > 20]
     c(1,:) = (/2.000e-03,1.000e-02, 3.84322e+03, 8.93072e+03,8.17947e+03, 3.71287e+03, 8.35739e+02, 7.46792e+01/)
     c(2,:) = (/1.000e-02,2.000e-01,-2.25831e+01, 1.15730e-02,-8.32355e-01,-1.17916e+00,-4.74033e-01,-8.64483e-02/)
     c(3,:) = (/2.000e-01,2.000e+00,-2.25560e+01, 3.28276e-01,1.39226e-01,-1.22085e-01,-2.76602e-01,-2.90494e-01/)
     c(4,:) = (/2.000e+00,2.000e+01,-2.25798e+01, 4.57017e-01,-1.59274e-01, 2.71952e-01,-1.71810e-01, 3.86609e-02/)
     c(5,:) = (/2.000e+01,1.000e+02,-1.73046e+01,-1.62462e+01,2.10079e+01,-1.31207e+01, 4.06935e+00,-5.00944e-01/)
  else if (name == 'he4') then
     ! original (for testing)
     c(1,:) = (/2.000e-03,2.000e-02, 1.441278e2,2.940867e2,1.761164e2,3.387430e1,-3.075936e0,-1.204179e0/)
     c(2,:) = (/1.000e-02,2.000e-01,-2.274210e1,-7.402954e-1,-2.177691e0,-2.426768e0,-1.026211e0,-1.798547e-1/)
     c(3,:) = (/2.000e-01,2.000e+00,-2.254156e1,3.503190e-1,1.210755e-1,-1.171573e-1,8.237547e-2,4.361719e-2/)
     c(4,:) = (/2.000e+00,2.000e+01,-2.258311e1,6.858961e-1,-8.628176e-1,1.205242e0,-7.386310e-1,1.686530e-1/)
     c(5,:) = (/2.000e+01,1.000e+02,-1.73046e+01,-1.62462e+01,2.10079e+01,-1.31207e+01, 4.06935e+00,-5.00944e-01/)
  else if (name == 'C') then
     ! Carbon (A=12,Z=6) [Agrees with Post for Te > 20]
     c(1,:) = (/2.000e-03,2.000e-02, 2.26196e+03, 5.27071e+03,4.80768e+03, 2.16700e+03, 4.83472e+02, 4.27841e+01/)
     c(2,:) = (/2.000e-02,2.000e-01, 5.01266e+01, 3.32680e+02,6.00316e+02, 5.17230e+02, 2.13035e+02, 3.36380e+01/)
     c(3,:) = (/2.000e-01,2.000e+00,-2.12371e+01,-3.13175e-01,7.60470e-01,-3.01657e-01, 7.63153e-02, 1.40114e-01/)
     c(4,:) = (/2.000e+00,2.000e+01,-2.12367e+01,-3.44789e-01,1.03546e+00,-1.01249e+00, 5.92047e-01,-1.43559e-01/)
     c(5,:) = (/2.00000e+01,1.00000e+04,-2.47680e+01,9.40818e+00,-9.65745e+00,4.99916e+00,-1.23738e+00,1.16061e-01/)
  else if (name == 'O') then
     ! Oxygen (A=16,Z=8) [Agrees with Post for Te > 20]
     c(1,:) = (/2.000e-03,2.000e-02,-4.52975e+02,-9.70819e+02,-8.54356e+02,-3.70311e+02,-7.89495e+01,-6.59517e+00/)
     c(2,:) = (/2.000e-02,2.000e-01,-5.85047e+01,-1.59971e+02,-2.41395e+02,-1.60654e+02,-4.43559e+01,-3.33047e+00/)
     c(3,:) = (/2.000e-01,2.000e+00,-2.07261e+01,-7.67836e-01,8.75731e-01,-1.08357e+00, 2.41126e+00, 5.75616e+00/)
     c(4,:) = (/2.000e+00,2.000e+01,-2.06891e+01,-1.09982e+00,2.00613e+00,-1.58614e+00, 6.69091e-01,-1.11969e-01/)
     c(5,:) = (/2.000e+01,1.000e+02,-2.78060e+01, 2.14606e+01,-2.66591e+01, 1.67083e+01,-5.19194e+00, 6.41029e-01/)
  else if (name == 'Ar') then
     ! Argon (A=40,Z=18) [Agrees with Post]
     c(1,:) = (/3.000e-02,2.000e-01,-2.05304e+01,-2.83429e+00,1.50690e+01, 3.51718e+01, 2.40012e+01, 5.07272e+00/)
     c(2,:) = (/2.000e-01,2.000e+00,-1.96520e+01,-1.17276e-01,7.83322e+00,-6.35158e+00,-3.05885e+01,-1.52853e+01/)
     c(3,:) = (/2.000e+00,2.000e+01,-1.97488e+01, 2.96484e+00,-8.82939e+00, 9.79100e+00,-4.96002e+00, 9.82003e-01/)
     c(4,:) = (/2.000e+01,1.000e+02,-2.11794e+01, 5.19148e+00,-7.43972e+00, 4.96902e+00,-1.55318e+00, 1.87705e-01/)
  else if (name == 'W') then
     ! Tungsten (A=184,Z=74) [Agrees with Post]
     c(1,:) = (/1.000e-01,2.000e-01, 5.34083e+00, 1.56088e+02,4.17170e+02, 5.50258e+02, 3.56758e+02, 9.04279e+01/)
     c(2,:) = (/2.000e-01,2.000e+00,-1.72389e+01, 5.42375e-02,-1.22107e+00, 4.41181e-01,-4.48582e+00,-7.83614e+00/)
     c(3,:) = (/2.000e+00,2.000e+01,-1.47488e+01,-1.43954e+01,2.10585e+01,-4.39475e+00,-1.10601e+01, 5.61699e+00/)
     c(4,:) = (/2.000e+01,1.000e+02,-2.62426e+02, 7.12559e+02,-8.25017e+02, 4.74241e+02,-1.35517e+02, 1.54189e+01/)
  else
     print *,'ERROR: '//name//'not found'
     stop
  endif

  do j=1,n

     t0 = te(j)
     x = log10(t0)

     ! Determine range
     success = 0
     do i=1,5
        if (t0 >= c(i,1) .and. t0 <= c(i,2)) then
           ! Shift coefficients
           a(1:6) = c(i,3:8)
           success = 1
           exit
        endif
     enddo

     ! Error
     if (success == 0) then
        print *,'ERROR: (post77) Temperature out of bounds for '//name
        stop
     endif

     ! Compute fit1
     s = 0.0
     do i=1,6
        s = s+a(i)*x**(i-1)
     enddo

     ! Lz in units of erg cm^3/s
     lz(j) = 10**s

  enddo

end subroutine post77
