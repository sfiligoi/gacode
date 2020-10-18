program curves

  implicit none

  character(len=3) :: name
  real, dimension(:), allocatable :: te
  real, dimension(:), allocatable :: lz

  integer, parameter :: n=128
  integer :: i
  real :: tmin,tmax,r,dx
  
  allocate(te(n))
  allocate(lz(n))

  tmin = 2.1e-3
  tmax = 99.0
  r = tmax/tmin
  dx = log(r)/(n-1)

  do i=1,n
     te(i) = tmin*exp((i-1)*dx)
  enddo
  
  ! Bremsstrahlung and line radiation (Post 1977) 
  name = 'He4'
  call post77(te,name,lz,n)

  open(unit=1,file=name//'.txt',status='replace')
  do i=1,n
     write(1,*) te(i),lz(i)
  enddo
  
  
end program curves
