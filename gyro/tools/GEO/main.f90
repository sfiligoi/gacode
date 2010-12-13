program main

  use GEO_interface

  !--------------------------------
  implicit none
  !
  integer :: n_theta,j
  !
  real :: x_IN(14)
  real :: y_IN(8,1)
  real :: dtheta,theta,pi
  !--------------------------------

  !-----------------------------------------------------------
  ! Equilibrium parameters:
  !
  !   1. rmin      = x_IN(1)
  !   2. rmaj      = x_IN(2)
  !   3. drmaj     = x_IN(3)
  !   4. zmag      = x_IN(4)
  !   5. dzmag     = x_IN(5)
  !   6. q         = x_IN(6)
  !   7. s         = x_IN(7)
  !   8. kappa     = x_IN(8)
  !   9. s_kappa   = x_IN(9)
  !  10. delta     = x_IN(10)
  !  11. s_delta   = x_IN(11)
  !  12. zeta      = x_IN(12)
  !  13. s_zeta    = x_IN(13)
  !  14. beta_star = x_IN(14)
  !
  ! where beta_star = beta_unit*dlnpdr > 0
  !-----------------------------------------------------------

  call GEO_alloc(1001,1)
  n_theta=128

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
  read(1,*) x_IN(11) 
  read(1,*) x_IN(12) 
  read(1,*) x_IN(13) 
  read(1,*) x_IN(14) 
  close(1)  

  y_IN = 0.0

  call GEO_do(x_IN,y_IN,1,0)
  call GEO_write('geo.out',1)

  open(unit=1,file='geov.out',status='replace')

  pi = 4.0*atan(1.0)
  dtheta = 2*pi/(n_theta-1)

  do j=1,n_theta
     theta = -pi+(j-1)*dtheta
     call GEO_interp(theta)  
     print 10,theta,GEO_gsin,GEO_usin,GEO_gcos1,GEO_ucos
     write(1,10) theta,GEO_gsin,GEO_usin,GEO_gcos1,GEO_ucos
  enddo

  close(1)

10 format(5(1pe12.5,1x))

end program main
