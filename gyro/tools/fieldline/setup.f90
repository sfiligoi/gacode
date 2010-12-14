subroutine setup(dir)

  use fieldline_input_data

  implicit none

  character (len=*) :: dir

  integer :: p
  integer :: i 
  integer :: j
  integer :: i_n

  ! Sanity check
  if (n_x_decimate > n_x-2) then
     print *,'ERROR: n_x_decimate > n_x-2'
     stop
  endif
  if (n_n_decimate > n_n-1) then
     print *,'ERROR: n_n_decimate > n_n-1'
     stop
  endif

  ! Assign central parameter values based on profile data:
  r0 = r(n_x/2+1)
  q0 = q(n_x/2+1)
  s0 = shat_s(n_x/2+1)
  rmaj0 = aspect_s(n_x/2+1)*r_s(n_x/2+1)

  ! Toroidal mode spacing
  if (n_n == 1) d_n = n0

  ! Grid spacing
  dr     = (r(n_x)-r(1))/(n_x-1)
  dx     = 2*pi/n_x
  dtheta = 2*pi/n_theta_plot

  ! Box lengths
  r_length  = n_x*dr
  y_length  = 2*pi*r0/(d_n*q0)
  r_natural = r0/(d_n*q0*s0)
  n_box     = r_length/r_natural

  ! Stochasticity parameter (epsilon is an artifical knob)
  lambda = eps_lambda*&
       (2*pi)**2*q0*rmaj0*(rhos_norm/r_length)*(rhos_norm/y_length)

  ! Array allocation
  allocate(rp(-n_x/2:n_x/2-1))
  allocate(xp(-n_x/2:n_x/2-1))
  allocate(qp(-n_x/2:n_x/2-1))
  allocate(n(0:n_n-1))
  allocate(phase(-n_x/2:n_x/2-1,0:n_n-1))
  allocate(theta(n_theta_plot+1))
  if (n_n == 1) then
     allocate(anp(n_theta_plot+1,-n_x/2:n_x/2-1,0:1))
  else
     allocate(anp(n_theta_plot+1,-n_x/2:n_x/2-1,-n_n+1:n_n-1))
  endif

  do j=1,n_theta_plot+1
     theta(j) = -pi+(j-1)*dtheta
  enddo

  if (n_n == 1) then
     ! Special case for single linear mode 
     n(0) = n0
  else
     do i_n=0,n_n-1
        n(i_n) = i_n*d_n
     enddo
  endif

  do i=-n_x/2,n_x/2-1
     rp(i) = i*dr
     xp(i) = i*dx
     qp(i) = q0+q0*s0*rp(i)/r0
     phase(i,:) = exp(-2*pi*i_c*n(:)*qp(i))
  enddo

  ! Add boundary element to input data
  a(n_theta_plot+1,:,:) = a(1,:,:)*phase(:,:)

  ! Add boundary element to G_theta - periodic, no phase factor
  g_theta(n_theta_plot+1,:) = g_theta(1,:) ! WG

  ! Now, construct Fourier coefficients (making sure to rescale A
  ! using gyroBohm factor rhos_norm).

  anp = (0.0,0.0)
  do p=-n_x/2,n_x/2-1
     do i=-n_x/2,n_x/2-1
        anp(:,p,0) = anp(:,p,0) + exp(-i_c*p*xp(i))*a(:,i,0)/n_x/rhos_norm
        if (n_n == 1) then
           anp(:,p,1) = anp(:,p,1) + exp(-i_c*p*xp(i))*conjg(a(:,i,0))/n_x/rhos_norm
        endif
        do i_n=1,n_n-1
           anp(:,p,i_n) = anp(:,p,i_n) + &
                exp(-i_c*p*xp(i))*a(:,i,i_n)/n_x/rhos_norm
           anp(:,p,-i_n) = anp(:,p,-i_n) + &
                exp(-i_c*p*xp(i))*conjg(a(:,i,i_n))/n_x/rhos_norm
        enddo
     enddo
  enddo

  !--------------------------------------------------
  ! Write diagnostics to screen
  !
  print *,'GYRO Poincare Surface Generation'
  print *
  print '(t2,a,a)','DATA: ',dir
  print * 
  print *,'Simulation parameters'
  print *
  print 30,'rmaj0',rmaj0
  print 30,'q0',q0
  print 30,'s0',s0
  print 30,'r_length/rhos',r_length/rhos_norm
  print 30,'y_length/rhos',y_length/rhos_norm
  print 30,'lambda',lambda
  print *
  print *,'Grid dimensions'
  print *
  print 10,'n_time_slice',n_time_slice
  print 10,'n_x',n_x
  print 10,'n_theta_plot',n_theta_plot
  print 10,'n_pass',n_pass
  print 10,'n_trap',n_trap
  print 10,'n_energy',n_energy
  print 10,'n_field',n_field
  print 10,'n_kinetic',n_kinetic
  print 10,'n_ion',n_ion
  print 10,'n_spec',n_spec
  print 10,'n_n',n_n
  print 10,'n_box',n_box
  print *
  print *,'Mode spectrum'
  print *
  if (n_n <= 4) then
     print "(t2,a,t20,': ',4(f4.2,1x))",'ky',ky
  else
     print "(t2,a,t20,': ',3(f4.2,1x),'... ',f4.2)",'ky',ky(0:2),ky(n_n-1)
  endif
  print *
  print *,'Integration parameters'
  print *
  print 15,'n_turn',n_turn
  print 15,'n_sub', n_sub
  print 15,'n_ic', n_ic
  print 15,'points',n_ic*n_turn
  print *
  print *,'Modes retained'
  print *
  print 10,'n_xd',n_x-n_x_decimate
  print 10,'n_nd',n_n-n_n_decimate
  print *
  !--------------------------------------------------

  deallocate(a)

10 format(t2,a,t20,': ',i3)
15 format(t2,a,t20,': ',i6)
  !20 format(t2,a,t20,': ',1pe12.5)
30 format(t2,a,t20,':',f11.6)

end subroutine setup
