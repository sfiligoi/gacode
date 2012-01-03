subroutine read_all

  use gyro2gkv_globals

  implicit none

  !---------------------------
  real :: x(2)
  integer :: ierr,i,j,is
  !---------------------------

  pi = 4*atan(1.0)

  !------------------------------------------
  ! Read some tags
  !
  open(unit=1,file='temp',status='old')
  read(1,'(a)') tag(1)
  read(1,'(a)') tag(2)
  read(1,'(a)') tag(3)
  close(1)
  !
  ! Set additional tags
  !
  if (boundary_method == 1) then
     tag(4) = "periodic"
  else
     tag(4) = "nonperiodic"
  endif

  lt(:) = len_trim(tag(:))
  !------------------------------------------

  !------------------------------------------
  ! Pull data out of profile_vugyro.out
  !
  open(unit=1,file=trim(tag(2))//'/out.gyro.profile',status='old')

  ! Essential scalars
  read(1,*) n_x
  read(1,*) n_theta_section
  read(1,*) n_pass
  read(1,*) n_trap
  read(1,*) n_energy
  read(1,*) n_theta_plot
  read(1,*) n0
  read(1,*) n_n
  read(1,*) d_n
  read(1,*) n_explicit_damp
  read(1,*) nonlinear_flag
  read(1,*) electron_method
  read(1,*) n_field
  read(1,*) n_ion
  read(1,*) n_kinetic
  read(1,*) n_spec
  read(1,*) field_r0_flag
  read(1,*) field_r0_grid
  read(1,*) n_grid_exp
  read(1,*) boundary_method

  allocate(r(n_x))
  allocate(q(n_x))

  ! Basic profile data
  read(1,*) r(:)
  read(1,*) q(:)

  close(1)
  !------------------------------------------

  !------------------------------------------
  ! Check for file existence
  !
  open(unit=1,file='exists_u',status='old')
  read(1,*) exists_u
  close(1)
  open(unit=1,file='exists_n',status='old')
  read(1,*) exists_n
  close(1)
  open(unit=1,file='exists_e',status='old')
  read(1,*) exists_e
  close(1)
  open(unit=1,file='exists_diff_i',status='old')
  read(1,*) exists_diff_i
  close(1)
  !------------------------------------------

  !------------------------------------------
  ! Check for gyro2gkv flags
  !
  open(unit=1,file='f_flag',status='old')
  read(1,*) f_flag
  close(1)
  !------------------------------------------

  !------------------------------------------
  ! (1) Create the refined theta grid, and the 
  !     toroidal mode grid.
  !
  allocate(theta_plot(n_theta_plot))
  allocate(n(n_n))

  if (n_theta_plot > 1) then

     ! Span [-pi,pi]

     do j=1,n_theta_plot
        theta_plot(j) = -pi+(j-1)*2*pi/n_theta_plot
     enddo
  else

     ! Special case for one point:
     !  only theta=0 (outboard midplane)

     theta_plot(1) = 0.0

  endif

  do i=1,n_n
     n(i) = (i-1)*d_n
  enddo
  !------------------------------------------

  !------------------------------------------
  ! (1a) Create the refined theta-r0 grid, and the 
  !     toroidal mode grid.
  !
  allocate(theta_r0_plot(field_r0_grid))

  if (field_r0_grid > 1) then

     ! Span [-pi,pi]

     do j=1,field_r0_grid
        theta_r0_plot(j) = -pi+(j-1)*2*pi/field_r0_grid
     enddo

  else

     ! Special case for one point:
     !  only theta=0 (outboard midplane)

     theta_r0_plot(1) = 0.0

  endif
  !------------------------------------------

  !------------------------------------------
  ! (2) Read time vector:
  ! 
  open(unit=1,file=trim(tag(2))//'/t.out',status='old')
  n_time = 0
  do 
     read(1,*,end=10) x(:)
     n_time = n_time+1 
  enddo
10 close(1)
  allocate(t(n_time))
  open(unit=1,file=trim(tag(2))//'/t.out',status='old')

  do i=1,n_time
     read(1,*) x(:)
     t(i) = x(2)
  enddo
  close(1)
  !------------------------------------------

  !------------------------------------------
  ! (3) Some IO
  !
  print *,'------------------'
  print 20,'n_n',n_n
  print 20,'n_x',n_x
  print 30,'t_max',t(n_time)
  print *,'------------------'

20 format(t2,a,t8,':',1x,i4)
30 format(t2,a,t8,':',1x,f6.1)
  !------------------------------------------

  !------------------------------------------
  ! (4) Read complex field (phi,a_parallel):
  !
  if (exists_u == 1) then
     if (f_flag == 0) then
     print *,'FOUND: u.out'

     open(unit=1,file=trim(tag(2))//'/out.gyro.moment_u',status='old')
     allocate(phi_in(2,n_theta_plot,n_x,n_field))
     allocate(phi0(2,n_x,n_n,n_time))
     if (n_field == 2) allocate(a_parallel0(2,n_x,n_n,n_time))

     print *,'u.out:' 
     do i_time=1,n_time
        if (i_time ==          1) print *,'[0/3 Done]'
        if (i_time ==   n_time/3) print *,'[1/3 Done]'
        if (i_time == 2*n_time/3) print *,'[2/3 Done]'
        if (i_time ==     n_time) print *,'[3/3 Done]'
        do i_n=1,n_n
           read(1,*) phi_in(:,:,:,:)
           phi0(:,:,i_n,i_time) = phi_in(:,n_theta_plot/2+1,:,1)
           if (n_field == 2) then 
              a_parallel0(:,:,i_n,i_time) = phi_in(:,n_theta_plot/2+1,:,2)
           endif
        enddo
     enddo
     close(1)
     endif
  endif
  !------------------------------------------

  !------------------------------------------
  ! (4a) Read full potential (including parallel modes)
  !
if (exists_u == 1) then  
   if (f_flag == 1) then
      open(unit=1,file=trim(tag(2))//'/out.gyro.moment_u',status='old')
      allocate(phi_in_full(2,n_theta_plot,n_x,n_field))
      allocate(phi0_full(2,n_theta_plot,n_x,n_n,n_time))
      allocate(a_parallel0_full(2,n_theta_plot,n_x,n_n,n_time))

      print *,'u.out: full flag on, will output full potential'
      print *,'n_theta_plot=', n_theta_plot
      
      do i_time=1,n_time
         
         if (i_time ==          1) print *,'[0/3 Done]'
         if (i_time ==   n_time/3) print *,'[1/3 Done]'
         if (i_time == 2*n_time/3) print *,'[2/3 Done]'
         if (i_time ==     n_time) print *,'[3/3 Done]'
         
         do i_n=1,n_n                      
            read(1,*) phi_in_full(:,:,:,:)
            if (n_field == 1) then
               phi0_full(:,:,:,i_n,i_time) = phi_in_full(:,:,:,1)
            else 
               phi0_full(:,:,:,i_n,i_time) = phi_in_full(:,:,:,1)
               a_parallel0_full(:,:,:,i_n,i_time) = phi_in_full(:,:,:,2)
            endif
         enddo
      enddo            
      close(1)
   endif
endif
!------------------------------------------

  !------------------------------------------
  ! (5) Read complex field at r0 (phi,a_parallel):
  !
  if (field_r0_flag == 1) then

     print *,'FOUND: field_r0.out'

     open(unit=1,file=trim(tag(2))//'/out.gyro.field_r0',status='old')

     allocate(phi_r0_in(2,field_r0_grid,n_field))
     allocate(phi_r0(2,field_r0_grid,n_n,n_time))

     print *,'field_r0.out:' 
     do i_time=1,n_time
        if (i_time ==          1) print *,'[0/3 Done]'
        if (i_time ==   n_time/3) print *,'[1/3 Done]'
        if (i_time == 2*n_time/3) print *,'[2/3 Done]'
        if (i_time ==     n_time) print *,'[3/3 Done]'
        do i_n=1,n_n
           read(1,*) phi_r0_in(:,:,:)
           phi_r0(:,:,i_n,i_time) = phi_r0_in(:,:,1)
        enddo
     enddo
     close(1)
  endif
  !------------------------------------------

  !------------------------------------------
  ! (6) Read moment_n for ions:
  !
  if (exists_n == 1) then
     if (f_flag == 0) then
     print *,'FOUND: moment_n.out'

     open(unit=1,file=trim(tag(2))//'/out.gyro.moment_n',status='old')

     allocate(moment_n_in(2,n_theta_plot,n_x,n_kinetic))
     allocate(moment_n(2,n_x,n_n,n_time))

     print *,'moment_n.out:' 
     do i_time=1,n_time
        if (i_time ==          1) print *,'[0/3 Done]'
        if (i_time ==   n_time/3) print *,'[1/3 Done]'
        if (i_time == 2*n_time/3) print *,'[2/3 Done]'
        if (i_time ==     n_time) print *,'[3/3 Done]'
        do i_n=1,n_n
           read(1,*) moment_n_in(:,:,:,:)
           moment_n(:,:,i_n,i_time) = moment_n_in(:,n_theta_plot/2+1,:,1)
        enddo
     enddo
     close(1)
     endif
  endif
  !------------------------------------------

  !------------------------------------------
  ! (6a) Read full moment_n (with parallel modes)
  !
  if (exists_n == 1) then
     if (f_flag == 1) then
        print *,'moment_n.out: full flag on, will output full moment_n'
     open(unit=1,file=trim(tag(2))//'/out.gyro.moment_n',status='old')
     allocate(moment_n_full_in(2,n_theta_plot,n_x,n_kinetic))
     allocate(moment_n_full(2,n_theta_plot,n_x,n_n,n_time))
     do i_time=1,n_time
        if (i_time ==          1) print *,'[0/3 Done]'
        if (i_time ==   n_time/3) print *,'[1/3 Done]'
        if (i_time == 2*n_time/3) print *,'[2/3 Done]'
        if (i_time ==     n_time) print *,'[3/3 Done]'
        do i_n=1,n_n
           read(1,*) moment_n_full_in(:,:,:,:)
           moment_n_full(:,:,:,i_n,i_time) = moment_n_full_in(:,:,:,1)
        enddo
     enddo
     close(1)
     endif 
  endif
  !------------------------------------------

  !------------------------------------------
  ! (7) Read moment_e for ions:
  !
  if (exists_e == 1) then
     if (f_flag == 0) then
     print *,'FOUND: moment_e.out'

     open(unit=1,file=trim(tag(2))//'/out.gyro.moment_e',status='old')

     allocate(moment_e_in(2,n_theta_plot,n_x,n_kinetic))
     allocate(moment_e(2,n_x,n_n,n_time))

     print *,'out.gyro.moment_e:' 
     do i_time=1,n_time
        if (i_time ==          1) print *,'[0/3 Done]'
        if (i_time ==   n_time/3) print *,'[1/3 Done]'
        if (i_time == 2*n_time/3) print *,'[2/3 Done]'
        if (i_time ==     n_time) print *,'[3/3 Done]'
        do i_n=1,n_n
           read(1,*) moment_e_in(:,:,:,:)
           moment_e(:,:,i_n,i_time) = moment_e_in(:,n_theta_plot/2+1,:,1)
        enddo
     enddo
     close(1)
     endif
  endif
  !------------------------------------------

  !------------------------------------------
  ! (7a) Read full moment_e (with parallel modes)
  !
  if (exists_e == 1) then
     if (f_flag == 1) then
        print *,'out.gyro.moment_e: full flag on, will output full moment_e'
     open(unit=1,file=trim(tag(2))//'/out.gyro.moment_e',status='old')
     allocate(moment_e_full_in(2,n_theta_plot,n_x,n_kinetic))
     allocate(moment_e_full(2,n_theta_plot,n_x,n_n,n_time))
     do i_time=1,n_time
        if (i_time ==          1) print *,'[0/3 Done]'
        if (i_time ==   n_time/3) print *,'[1/3 Done]'
        if (i_time == 2*n_time/3) print *,'[2/3 Done]'
        if (i_time ==     n_time) print *,'[3/3 Done]'
        do i_n=1,n_n
           read(1,*) moment_e_full_in(:,:,:,:)
           moment_e_full(:,:,:,i_n,i_time) = moment_e_full_in(:,:,:,1)
        enddo
     enddo
     close(1)
     endif 
  endif

  !------------------------------------------

  !------------------------------------------
  ! (8) Read diff_i:
  !
  if (exists_diff_i == 1) then

     print *,'FOUND: out.gyro.diff_i'

     open(unit=1,file=trim(tag(2))//'/out.gyro.diff_i',status='old')

     allocate(diff0(n_kinetic,n_field,2))
     allocate(diff_density(n_x,n_kinetic,n_time))
     allocate(diff_energy(n_x,n_kinetic,n_time))
     if (n_field == 2) then
        allocate(diff_density_em(n_x,n_kinetic,n_time))
        allocate(diff_energy_em(n_x,n_kinetic,n_time))
     endif
     
     do i_time=1,n_time !n_time ! 2 replaced from n_time        
        do i=1,n_x           
           read(1,*,end=35) diff0(:,:,:)
           !print *, "i_time = ", i_time, " and i = ", i
           diff_density(i,:,i_time) = diff0(:,1,1)
           diff_energy(i,:,i_time) = diff0(:,1,2)
           if (n_field == 2) then
              diff_density_em(i,:,i_time) = diff0(:,2,1)
              diff_energy_em(i,:,i_time) = diff0(:,2,2)
           endif
        enddo        
     enddo
35   close(1)
     
     print 40,'Ave(D)  :',(sum(diff_density(:,is,n_time))/n_x,is=1,n_kinetic)
     print 40,'Ave(Chi):',(sum(diff_energy(:,is,n_time))/n_x,is=1,n_kinetic)
  endif

40 format(t2,a,3(1x,1pe11.4))
  !------------------------------------------

end subroutine read_all
