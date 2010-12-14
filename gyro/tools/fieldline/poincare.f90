program poincare

  use fieldline_input_data
  use runge_kutta

  implicit none

  integer :: i
  integer :: n_plot
  integer :: i_turn
  integer :: i_sub
  integer :: i_ic
  integer, parameter :: plot=1

  real :: h
  real :: z

  real :: r_in
  real :: q_in
  real :: xi_initial

  real, allocatable, dimension(:) :: y0
  real, allocatable, dimension(:) :: y
  real, allocatable, dimension(:) :: yp

  character (len=100) :: simdir

  !---------------------------------------------------
  ! INPUT data
  !
  open(unit=1,file='poincare.dat',status='old')
  read(1,*) simdir
  read(1,*) n_turn
  read(1,*) n_sub
  read(1,*) n_ic
  read(1,*) n_time_slice
  read(1,*) n_x_decimate
  read(1,*) n_n_decimate
  read(1,*) eps_lambda
  close(1)
  !---------------------------------------------------

  ! Read GYRO data
  call read_input(trim(simdir))

  ! Compute parameters, work arrays, etc.
  call setup(trim(simdir))

  !-----------------------------------------------------
  ! y(1) -> zeta  (radial phase)
  ! y(2) -> alpha (toroidal phase) 
  ! y(3) -> time
  !
  allocate(y(n_eq))
  allocate(y0(n_eq))
  allocate(yp(n_eq))
  !-----------------------------------------------------

  !----------------------------------------------------
  ! Plot fields as functions of (r,theta) as a sanity check.
  !
  if (plot == 1) then

     n_plot = 256

     open(unit=1,file='xi.out',status='replace')

     write(1,*) n_plot
     z = pi
     do i=1,n_plot
        y(1) = -z+(i-1)*2*z/(n_plot-1)
        y(2) = pi/2
        y(3) = pi/2
        call func(y,yp,n_eq)
        write(1,11) y(1)/(2*pi),yp(1:2)
     enddo
     close(1)

     open(unit=1,file='alpha.out',status='replace')

     write(1,*) n_plot
     z = pi
     do i=1,n_plot
        y(1) = pi/2
        y(2) = -z+(i-1)*2*z/(n_plot-1)
        y(3) = pi/2
        call func(y,yp,n_eq)
        write(1,11) y(2)/(2*pi),yp(1:2)
     enddo
     close(1)

     ! Inspect G_theta
     open(unit=1,file='gtheta.out',status='replace') ! WG
     do i=1,n_theta_plot+1 ! WG
        write(1,12) g_theta(i,:) ! WG
     enddo ! WG
     close(1) ! WG

  endif

11 format(t2,8(1pe12.5,1x))
12 format(t2,32(1pe12.5,1x)) ! WG
  !----------------------------------------------------

  h = 2*pi/n_sub

  open(unit=1,file='surface.out',status='replace')
  write(1,*) n_ic

  do i_ic=1,n_ic

     print '(t2,a,i2,a,i2)','Working on orbit ',i_ic,' of ',n_ic 

     ! New initial conditions

     y0(1) = -pi+(i_ic-1)*2*pi/n_ic+pi/n_ic
     y0(2) = -pi+(i_ic-1)*2*pi/n_ic+pi/n_ic
     y0(3) = -pi

     xi_initial = y0(1)

     do i_turn=1,n_turn

        ! Integrate for n_turn loops in theta (tau)

        do i_sub=1,n_sub

           ! Use n_sub intervals to integrate from -pi to pi:

           call rk4(y,y0,h,n_eq)

           !y(2)  = -pi+modulo(y(2)+pi,2*pi)
           y0(:) = y(:)

        enddo

        ! r = xi L/(2pi)
        r_in  = y0(1)*r_length/(2*pi)
        q_in  = q0+q0*s0*r_in/r0

        ! Jump in Clebsch coordinate: 2 pi q Delta_n 
        y0(2) = y0(2)+2*pi*q_in*d_n

        ! RHS periodic in alpha=y0(2)
        y0(2) = -pi+modulo(y0(2)+pi,2*pi)

        ! Move theta from pi back to -pi.
        y0(3) = -pi

        write(1,11) y0(1:2)/(2*pi),(y0(1)-xi_initial)**2/(2*n_turn)

     enddo

  enddo ! i_ic

  print *,'Done.'

  close(1)
  !-----------------------------------------------------

  print *,'Wrote: xi.out, alpha.out, surface.out'

end program poincare

