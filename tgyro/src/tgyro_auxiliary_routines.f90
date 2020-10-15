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

subroutine rad_brem(ne,te,zeff,s_brem,n)

  implicit none
  
  integer, intent(in) :: n
  real, intent(in) :: ne(n)
  real, intent(in) :: te(n)
  real, intent(in) :: zeff(n)
  real, intent(inout) :: s_brem(n)
  
  ! Bremsstrahlung radiation [erg/cm^3/s]
  ! - From NRL formulary
  ! - 1 W/cm^3 = 1e7 erg/cm^3/s
  
  s_brem = 1e7*1.69e-32*ne**2*sqrt(te)*zeff
  
end subroutine rad_brem

subroutine rad_sync(b_ref,ne,te,s_sync,n)

  use tgyro_globals , only : pi,e,me,c,k,aspect_rat,r_min

  implicit none
  
  integer, intent(in) :: n
  real, intent(in) :: b_ref(n)
  real, intent(in) :: ne(n)
  real, intent(in) :: te(n)
  real, intent(inout) :: s_sync(n)

  real :: g,phi,wpe,wce
  ! Reflection coefficient (Rosenbluth)
  real, parameter :: r_coeff=0.8

  integer :: i

  !-------------------------------------------------------
  ! Synchrotron radiation
  ! - Trubnikov, JETP Lett. 16 (1972) 25.
  do i=1,n   
     wpe = sqrt(4*pi*ne(i)*e**2/me)
     wce = e*abs(b_ref(i))/(me*c)
     g   = k*te(i)/(me*c**2)
     phi = 60*g**1.5*sqrt((1.0-r_coeff)*(1+1/aspect_rat/sqrt(g))/(r_min*wpe**2/c/wce))

     s_sync(i) = me/(3*pi*c)*g*(wpe*wce)**2*phi
  enddo
  !-------------------------------------------------------

end subroutine rad_sync

subroutine rad_alpha(ne,ni,te,ti,s_alpha_i,s_alpha_e,frac_ai,e_cross,n,nion)

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
  real, intent(inout) :: s_alpha_i(n)
  real, intent(inout) :: s_alpha_e(n)
  real, intent(inout) :: frac_ai(n)
  real, intent(inout) :: e_cross(n)

  real, external :: sivukhin
  real, external :: sigv

  real :: x_a
  real :: n_d,n_t
  real :: s_alpha,sn_alpha
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
     sn_alpha = n_d*n_t*sigv(ti(1,i)/1e3,'bosch') * tgyro_input_fusion_scale

     s_alpha      = sn_alpha*e_alpha
     s_alpha_i(i) = s_alpha*frac_ai(i)
     s_alpha_e(i) = s_alpha*(1-frac_ai(i))

  enddo

  deallocate(c_a)

end subroutine rad_alpha

subroutine radiation(&
     te,&
     ne,&
     ni,&
     z,&
     r,&
     rmajor,&
     mhd_btor,&
     name,&
     nmax,&
     nr,&
     qbrem,&
     qsync,&
     qline)

  implicit none

  integer, intent(in) :: nmax,nr
  character(len=*), intent(in), dimension(nmax) :: name

  real, intent(in), dimension(nr) :: r,te,ne
  real, intent(in), dimension(nr,nmax) :: ni
  real, intent(in), dimension(nmax) :: z
  real :: mhd_btor 
  real :: rmajor 

  real, intent(inout), dimension(nr) :: qbrem,qsync,qline

  ! Internals

  integer :: nprim,nimp
  character(len=3), dimension(:), allocatable :: namep
  character(len=3), dimension(:), allocatable :: namei

  real, dimension(:,:), allocatable :: brems_nions
  real, dimension(:), allocatable :: prad_imp

  real, external :: fradhe,pwr_brem_factr

  real, parameter :: pi=3.1415926535
  real, parameter :: cyclotron_reflection = 0.0
  real, parameter :: kevpjou = 6.2415097e15      ! Kev/jou

  integer :: j,i,k

  real :: cyclo,teerg,wpsq,chi,phi,vdotsq,wbx
  real :: rm,cee,xmasse,charge,dene,btor

  ! Need to map from GACODE to ONETWO species variables
  nprim = 0
  nimp = 0

  do i=1,nmax
     if (z(i) > 1.0) then
        nimp = nimp + 1
     else
        nprim = nprim + 1
     endif
  enddo

  allocate(namep(nprim))
  allocate(namei(nimp))
  allocate(prad_imp(nr))
  allocate(brems_nions(nr,nprim))
  
  namep(1:nprim) = name(1:nprim)
  namei(1:nimp) = name(nprim+1:nprim+nimp)

  prad_imp(:) = 0.0
  brems_nions(:,:) = 0.0

  qbrem = 0.0
  do j=1,nr
     do k=1,nprim
        if (namep(k) == 'He') then
           brems_nions(j,k) = fradhe(te(j))*ni(j,k)    
        else ! ok if namep ='dt'
           brems_nions(j,k) = pwr_brem_factr(te(j),z(k))*ni(j,k)
        endif
        qbrem(j) = qbrem(j) + brems_nions(j,k)
     enddo
  enddo

  if (nimp > 0) then
     do i=1,nimp
        call nradfit(te,ni(:,nprim+i),namei(i),prad_imp,nr)
        qline = qline+prad_imp
     enddo
  endif

  qbrem = 10*ne*qbrem
  qline = 10*ne*qline

  ! ----------------------------------------------------------------------
  !  add radiative energy loss due to cyclotron (synchrotron) radiation.
  !  a realistic treatment of this term would require knowledge of the
  !  emission and absoption characteristics of the plasma.  the present
  !  calculation is a simplified model which will probably only give
  !  good global results. calculation is in CGS (Gaussian)
  !  then converted to keV/cm**3-s.
  !  ref: trubnikov: jetp letters 16, 25 (1972)
  !
  !  cyclotron_reflection == wall reflection coefficient
  !  If cyclotron_reflection =1 then
  !  all the synchrotron radiation is assumed reabsorbed so
  !  no addition to qrad is made:
  ! ----------------------------------------------------------------------

  cyclo = 0.0
  qsync = 0.0

  if (cyclotron_reflection < 1.0) then
     charge = 4.8032067e-10        ! electron charge (esu)
     xmasse = 9.1093897e-28        ! electron mass (g)
     cee    = 2.99792458e10        ! speed of light (cm/s)
     rm     = r(nr)/rmajor            ! 
     btor   = mhd_btor*1e4            ! gauss
     wbx     = charge*abs(btor)/(xmasse*cee)
     do j=1,nr
        dene    = ne(j)*1e-6
        teerg   = te(j)*1.6e-9
        wpsq    = 4.0*pi*dene*charge**2/xmasse
        chi     = rm*sqrt(xmasse*cee**2/teerg)
        phi     = 60.0*(teerg/(xmasse*cee**2))**1.5 &
             * sqrt(cee*wbx/(r(nr)*100.0*wpsq)*(1.0-cyclotron_reflection)*(1.0+chi))
        vdotsq  = wbx**2*2.0*teerg/xmasse
        cyclo   = dene/1.5*charge**2/cee**3*vdotsq*phi/1.6e-9
        qsync(j) = cyclo*1.e6   ! kev/(m^3 sec)
     enddo
  endif

  qsync = qsync*10/kevpjou

end subroutine radiation

subroutine nradfit(te,ni,namei,qline,nr)

  implicit none

  integer, intent(in) :: nr

  character(len=3), intent(in) :: namei
  real, intent(in), dimension(nr) :: te,ni
  real, intent(inout), dimension(nr) :: qline
  real, dimension(5,8) :: c
  real, dimension(6) :: a

  integer :: i,j,success
  real :: x,t0,s

  c = 0.0
  success = 0

  ! Set coefficients

  if (namei == 'Ar') then
     ! Argon
     c(1,:) = (/3.000e-02,2.000e-01,-2.05304e+01,-2.83429e+00,1.50690e+01, 3.51718e+01, 2.40012e+01, 5.07272e+00/)
     c(2,:) = (/2.000e-01,2.000e+00,-1.96520e+01,-1.17276e-01,7.83322e+00,-6.35158e+00,-3.05885e+01,-1.52853e+01/)
     c(3,:) = (/2.000e+00,2.000e+01,-1.97488e+01, 2.96484e+00,-8.82939e+00, 9.79100e+00,-4.96002e+00, 9.82003e-01/)
     c(4,:) = (/2.000e+01,1.000e+02,-2.11794e+01, 5.19148e+00,-7.43972e+00, 4.96902e+00,-1.55318e+00, 1.87705e-01/)
  else if (namei == 'C') then
     ! Carbon
     c(1,:) = (/2.000e-03,2.000e-02, 2.26196e+03, 5.27071e+03,4.80768e+03, 2.16700e+03, 4.83472e+02, 4.27841e+01/)
     c(2,:) = (/2.000e-02,2.000e-01, 5.01266e+01, 3.32680e+02,6.00316e+02, 5.17230e+02, 2.13035e+02, 3.36380e+01/)
     c(3,:) = (/2.000e-01,2.000e+00,-2.12371e+01,-3.13175e-01,7.60470e-01,-3.01657e-01, 7.63153e-02, 1.40114e-01/)
     c(4,:) = (/2.000e+00,2.000e+01,-2.12367e+01,-3.44789e-01,1.03546e+00,-1.01249e+00, 5.92047e-01,-1.43559e-01/)
     c(5,:) = (/2.00000e+01,1.00000e+04,-2.47680e+01,9.40818e+00,-9.65745e+00,4.99916e+00,-1.23738e+00,1.16061e-01/)
  else if (namei == 'W') then
     ! Tungsten
     c(1,:) = (/1.000e-01,2.000e-01, 5.34083e+00, 1.56088e+02,4.17170e+02, 5.50258e+02, 3.56758e+02, 9.04279e+01/)
     c(2,:) = (/2.000e-01,2.000e+00,-1.72389e+01, 5.42375e-02,-1.22107e+00, 4.41181e-01,-4.48582e+00,-7.83614e+00/)
     c(3,:) = (/2.000e+00,2.000e+01,-1.47488e+01,-1.43954e+01,2.10585e+01,-4.39475e+00,-1.10601e+01, 5.61699e+00/)
     c(4,:) = (/2.000e+01,1.000e+02,-2.62426e+02, 7.12559e+02,-8.25017e+02, 4.74241e+02,-1.35517e+02, 1.54189e+01/)
  else if (namei == 'He4') then
     ! Helium 4 (A=4,Z=2)
     c(1,:) = (/2.000e-03,1.000e-02, 3.84322e+03, 8.93072e+03,8.17947e+03, 3.71287e+03, 8.35739e+02, 7.46792e+01/)
     c(2,:) = (/1.000e-02,2.000e-01,-2.25831e+01, 1.15730e-02,-8.32355e-01,-1.17916e+00,-4.74033e-01,-8.64483e-02/)
     c(3,:) = (/2.000e-01,2.000e+00,-2.25560e+01, 3.28276e-01,1.39226e-01,-1.22085e-01,-2.76602e-01,-2.90494e-01/)
     c(4,:) = (/2.000e+00,2.000e+01,-2.25798e+01, 4.57017e-01,-1.59274e-01, 2.71952e-01,-1.71810e-01, 3.86609e-02/)
     c(5,:) = (/2.000e+01,1.000e+02,-1.73046e+01,-1.62462e+01,2.10079e+01,-1.31207e+01, 4.06935e+00,-5.00944e-01/)
  else
     print *,'ERROR: '//namei//'not found'
     stop
  endif

  do j=1,nr

     t0 = te(j)
     x = log10(t0)

     ! Determine range
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
        print *,'Temperature out of bounds in nradfit' ; stop
     endif

     ! Compute fit1
     s = 0.0
     do i=1,6
        s = s+a(i)*x**(i-1)
     enddo

     qline(j) = 10**s

  enddo

  ! W/m^3
  qline = qline*ni*1e-13

end subroutine nradfit

real function fradhe(te)

  !------------------------------------------------------------------------
  !     this function calculates the radiation losses from helium as a
  !     function of elec. temp. (keV) using a formula derived by
  !     George Hopkins and John M. Rawls (ga-a14141 impurity radiations)
  !
  !     temkev = temperature of electrons in keV
  !     value returned is in units of  (watts-m**3)
  !-------------------------------------------------------------------------

  implicit none
  
  real :: te,x,y,sqrtem

  x = te
  x = min(100.0,te)
  x = max(0.01,x)
  if (x > 0.4) goto 10
  x = log10(x)
  y = 10.0**(-35.632+x*(0.473114+x*0.521255))
  goto 20
10 sqrtem = sqrt(x)
  y = 2.136e-36*sqrtem+2.08e-37/sqrtem+4.8e-38/(x**1.34)
20 fradhe  = y
 
end function fradhe

real function pwr_brem_factr(te,z)

  !----------------------------------------------------------------------------------
  !-- this is power radiated in bremsstrhalung divided by electron and ion density
  !-- eg bremssthralung power(watts/m^3)  = pwr_brem_factr *ne*ni
  !-- peculiar way of doing this is due to onetwo compatibility
  !-- See Wesson,Tokamaks,pg 100 eq 4.9.2, (1987 )
  !--------------------------------------------------------------------------HJ------ 

  implicit none

  real :: te,z  ! Te is in kev in this formula

  pwr_brem_factr = 5.35e-37*sqrt(te)*z**2

end function pwr_brem_factr
