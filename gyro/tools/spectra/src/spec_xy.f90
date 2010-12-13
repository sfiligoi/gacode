program spec_xy

  use input_data

  !-------------------------------------------------------
  implicit none
  !
  integer :: i_k,n_k
  integer :: i_moment
  integer :: i_r,i_n,i_t,i0,i1
  integer :: p,i_s

  real :: pi
  real ::  d1,d2
  real :: kx_min,kx_max
  real :: ky_min,ky_max
  real :: s0,s1
  real, dimension(:,:), allocatable :: f2

  complex :: tpi
  complex, dimension(:,:), allocatable :: f
  complex, dimension(:,:), allocatable :: fp
  complex, dimension(:,:), allocatable :: cr

  character (len=50) :: dir 
  character (len=20) :: name
  !-------------------------------------------------------

  pi  = 4.0*datan(1.0)
  tpi = 2.0*pi*(0.0,1.0)

  !-------------------------------------------------------
  ! Read input/configuration file
  !
  open(unit=1,file='spec_xy.dat',status='old')
  read(1,*) dir
  read(1,*) name
  read(1,*) i_moment
  read(1,*) i_k
  read(1,*) i0
  read(1,*) i1
  close(1)

  if (n_field == 2) then
     print *,"n_field must be 1"
     stop
  endif
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Read simulation input data
  !
  print *,"DIR: ", dir

  call read_input(dir)

  if (i_moment == 1) then
     open(unit=1,file=trim(dir)//'u.out',status='old')
     n_k = 1 ! Number of field indices
     print *,"Reading u.out."
  else
     open(unit=1,file=trim(dir)//'moment_n.out',status='old')
     n_k = n_kinetic ! Number of field indices
     print *,"Reading moment_n.out"
  endif

  print *,"rho_s:", rho
  print *,"n_k:", n_k
  print *,"i_k:", i_k
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Allocate large arrays and read sim data
  !
  allocate(f(0:n_r-1,0:n_n-1))
  allocate(fp(0:n_r/2-1,0:n_n-1))
  allocate(f2(0:n_r/2-1,0:n_n-1))
  allocate(cr(0:n_r-1,0:n_r-1))

  do p=0,n_r-1
     do i_r=0,n_r-1
        cr(p,i_r) = exp(p*(tpi*i_r)/n_r)/n_r
     enddo
  enddo

  f2 = 0.0

  do i_t=1,i1

     do i_n=0,n_n-1
        do i_s=1,n_k
           do i_r=0,n_r-1
              read(1,*) d1,d2
              if (i_s == i_k) f(i_r,i_n) = d1+(0.0,1.0)*d2
           enddo
        enddo
     enddo

     if (i_t > i0) then

        fp(:,:) = (0.0,0.0)
        do i_n=0,n_n-1
           do i_r=0,n_r-1
              do p=0,n_r/2-1
                 fp(p,i_n) = fp(p,i_n)+cr(p,i_r)*f(i_r,i_n)
              enddo
           enddo
        enddo

        f2(:,:) = f2(:,:)+real(fp(:,:)*conjg(fp(:,:)))/rho**2

     endif

     !     print *,i_t

  enddo

  f2(:,:) = f2(:,:)/(i1-i0)

  close(1)
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Write complete spectral data
  !
  open(unit=1,file=trim(dir)//'ky.out',status='replace')
  write(1,*) n_n
  do i_n=0,n_n-1
     write(1,*) ky(i_n)
  enddo
  close(1)

  open(unit=1,file=trim(dir)//'kx.out',status='replace')
  write(1,*) n_r/2
  do p=0,n_r/2-1
     write(1,*) kx(p)
  enddo
  close(1)

  open(unit=1,file=trim(dir)//'k.out',status='replace')
  do i_n=0,n_n-1
     do p=0,n_r/2-1
        write(1,*) f2(p,i_n)
     enddo
  enddo
  close(1)
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Analysis
  !
  print *,'kx_max:  ',kx(n_r/2-1)
  print *,'ky_max:  ',ky(n_n-1)
  print *,'sqrt(sum(f2)): ',sqrt(sum(f2))
  print *,'sum(sqrt(f2)): ',sum(sqrt(f2))

  kx_min = 4.0
  kx_max = 8.0
  ky_min = 0.1
  ky_max = 0.2

  s0 = 0.0
  s1 = 0.0
  do p=0,n_r/2-1
     if (kx(p) >= kx_min .and. kx(p) <= kx_max) then
        do i_n=0,n_n-1
           if (ky(i_n) >= ky_min .and. ky(i_n) <= ky_max) then
              s0 = s0+f2(p,i_n)
              s1 = s1+sqrt(f2(p,i_n))
           endif
        enddo
     endif
  enddo
  print *,'sqrt(I(f2)): ',s0
  print *,'I(sqrt(f2)): ',s1
  !-------------------------------------------------------

end program spec_xy
