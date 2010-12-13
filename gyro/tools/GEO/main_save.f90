program main

  use GEO_interface

  !--------------------------------
  implicit none
  !
  integer :: i
  !
  real :: x_IN(10)
  real :: x,t,gsin,gcos,q,s,kxoky,a,beta
  !--------------------------------

  !-----------------------------------------------------------
  ! Equilibrium parameters:
  !
  !   1. rmin      = x_IN(1)
  !   2. rmaj      = x_IN(2)
  !   3. q         = x_IN(3)
  !   4. delta     = x_IN(4)
  !   5. kappa     = x_IN(5)
  !   6. shift     = x_IN(6)
  !   7. shat      = x_IN(7)
  !   8. s_delta   = x_IN(8)
  !   9. s_kappa   = x_IN(9)
  !  10. beta_star = x_IN(10)
  !
  ! where beta_star = beta_unit*dlnpdr > 0
  !-----------------------------------------------------------

  call GEO_alloc(51,1)

  open(unit=1,file='geo.dat',status='old')
  read(1,*) x_IN(1) 
  read(1,*) x_IN(2) 
  read(1,*) x_IN(3) 
  read(1,*) x_IN(4) 
  read(1,*) x_IN(5) 
  read(1,*) x_IN(6) 
  read(1,*) x_IN(7) 
  read(1,*) x_IN(8) 
  read(1,*) x_IN(9) 
  read(1,*) x_IN(10) 
  close(1)  

  call GEO_do(x_IN)
  call GEO_write('geo.out',1)

  x = x_IN(1)/x_IN(2)
  q = x_IN(3)
  s = x_IN(7)
  a = x_IN(6)/x
  beta = x_IN(10)*x_IN(2)

  open(unit=1,file='geov.out',status='replace')
  do i=1,size(GEOV_gsin)
     t = GEOV_theta(i)
     write(1,10) t,&
          GEOV_gsin(i)-sin(t),&
          GEOV_gcos1(i)-cos(t),&
          GEOV_kxoky(i)-s*t+q*q*beta*sin(t)
  enddo
  close(1)

  open(unit=1,file='circlev.out',status='replace')
  do i=1,size(GEOV_gsin)
     t = GEOV_theta(i)
     gsin = -x*sin(t)*cos(t)
     gcos = -x*(q**2*cos(t)**2-1)/q**2
     kxoky = -x*(((-3*a+1)*s+2*a+2)*sin(t)+(2*a-1)*s*t*cos(t))- &
          ((6*a-1)*q**2*beta*sin(2*t)+((-16*a+8)*q*q*beta*cos(t))*sin(t))*x/8

     write(1,10) t,&
          gsin,&
          gcos,&
          kxoky
  enddo
  close(1)

10 format(4(1pe12.5,1x))

end program main
