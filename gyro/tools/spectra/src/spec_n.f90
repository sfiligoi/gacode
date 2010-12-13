program spec_n

  use input_data

  implicit none

  integer :: i_r,i_n,i_t,i_s
  integer :: i0,i1,i_k
  integer :: n_s,n_k
  real ::  d1,d2
  complex, dimension(:,:), allocatable :: f
  complex, dimension(:,:), allocatable :: fp
  complex, dimension(:,:), allocatable :: fpp
  complex, dimension(:), allocatable :: temp
  real, dimension(:), allocatable :: z_i
  real, dimension(:), allocatable :: z_s
  real, dimension(:), allocatable :: z_n
  real :: pi
  complex :: tpi

  character (len=60) :: dir 
  character (len=30) :: file

  pi  = 4.0*datan(1.0)
  tpi = 2*pi*(0.0,1.0)

  !-------------------------------------------------------
  file = "c64x64.Bnoi.m20" 
  dir = "/home/candy/SIM/etg-ki/"//trim(file)//"/" 

  n_s = 3 ! Number of n_species
  n_k = 1 ! Number of field indices
  i_k = 1 ! desired field index

  call read_input(dir,n_s)

  i0 = 700
  i1 = 1378

  !-------------------------------------------------------

  open(unit=1,file=trim(dir)//'u.out',status='old')

  allocate(f(0:n_r-1,0:n_n-1))
  allocate(fp(0:n_r-1,0:n_n-1))
  allocate(fpp(0:n_r-1,0:n_n-1))
  allocate(temp(0:n_r-1))
  allocate(z_i(0:n_n-1))
  allocate(z_s(0:n_n-1))
  allocate(z_n(0:n_n-1))

  z_i = 0.0
  z_s = 0.0
  z_n = 0.0

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

        do i_n=0,n_n-1

           temp(:) = f(:,i_n)
           z_i(i_n) = z_i(i_n)+sum(real(temp(:)*conjg(temp(:))))/n_r

           call per_deriv(fp(:,i_n),f(:,i_n),r,n_r)
           call per_deriv(fpp(:,i_n),fp(:,i_n),r,n_r)

           temp(:) = (ky(i_n)/rho)**2*f(:,i_n)+fpp(:,i_n)
           z_s(i_n) = z_s(i_n)+sum(real(temp(:)*conjg(temp(:))))/n_r

           temp(:) = ky(i_n)/rho*fp(:,i_n)
           z_n(i_n) = z_n(i_n)+sum(real(temp(:)*conjg(temp(:))))/n_r

        enddo

     endif

     print *,i_t

  enddo

  close(1)

  open(unit=1,file=trim(file)//".spec_n",status='replace')
  do i_n=0,n_n-1
     write(1,10) ky(i_n),&
          z_i(i_n)/(i1-i0)/rho**2, &
          z_s(i_n)/(i1-i0)*rho**2, &
          z_n(i_n)/(i1-i0)*rho**2

  enddo
  close(1)

10 format(4(1pe12.5,1x))

end program spec_n
