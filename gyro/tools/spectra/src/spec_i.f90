program spec_i

  use input_data

  implicit none

  integer :: i_r,i_n,i_t,i_s
  integer :: i0,i1,i_k
  integer :: n_s,n_k
  real ::  d1,d2
  complex, dimension(:,:), allocatable :: f
  real, dimension(:), allocatable :: z
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
  allocate(z(0:n_r-1))

  z = 0.0

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

           if (ky(i_n) > 1.0) then

              z(:) = z(:)+real(f(:,i_n)*conjg(f(:,i_n)))

           endif

        enddo

     endif

     print *,i_t

  enddo

  close(1)

  open(unit=1,file=trim(file)//".spec_i",status='replace')
  do i_r=0,n_r-1
     write(1,10) r(i_r+1),z(i_r)/(i1-i0)/rho**2
  enddo
  close(1)

10 format(4(1pe12.5,1x))

end program spec_i
