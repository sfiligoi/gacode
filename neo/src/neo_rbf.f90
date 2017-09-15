subroutine neo_rbf(x0,c0)

  implicit none

  real, dimension(6), intent (in) :: x0
  real, dimension(6), intent (out) :: c0

  integer :: p,q,k
  integer, dimension(9) :: n

  integer :: ntot
  real, dimension(:,:), allocatable :: indata
  real, dimension(:,:), allocatable :: ingeodata
  real, dimension(:,:), allocatable :: outdata
  real, dimension(:,:), allocatable :: x
  real, dimension(:), allocatable :: w
  real, dimension(16) :: dum
  real, dimension(6) :: xmin,xmax,xs0
  integer :: stat

  real, dimension(:,:), allocatable :: a
  real, dimension(:,:), allocatable :: b
  integer, dimension(:), allocatable :: ipiv
  integer :: info

  open(unit=1,file='input.dat',status='old',iostat=stat)
  if (stat /= 0) then
     call neo_error('ERROR: (NEO) RBF files not available')
     return
  endif
  do k=1,9
     read(1,*) n(k)
     read(1,*) dum(1:n(k))
  enddo
  close(1)

  ntot = n(1)*n(2)*n(3)*n(4)*n(5)*n(6)*n(7)*n(8)*n(9)
  allocate(indata(9,ntot))
  allocate(ingeodata(12,ntot))
  allocate(outdata(6,ntot))
  allocate(w(ntot))
  allocate(x(6,ntot))

  open(unit=1,file='indata.dat',status='old',iostat=stat)
  if (stat /= 0) then
     call neo_error('ERROR: (NEO) RBF files not available')
     return
  endif
  open(unit=2,file='ingeodata.dat',status='old',iostat=stat)
  if (stat /= 0) then
     call neo_error('ERROR: (NEO) RBF files not available')
     return
  endif
  open(unit=3,file='outdata.dat',status='old',iostat=stat)
  if (stat /= 0) then
     call neo_error('ERROR: (NEO) RBF files not available')
     return
  endif

  do p=1,ntot

     do k=1,9
        read(1,*) indata(k,p)
     enddo
     do k=1,6
        read(3,*) outdata(k,p)
     enddo
     do k=1,12
        read(2,*) ingeodata(k,p)
     enddo

  enddo

  close(1)
  close(2)
  close(3)

  x(1:5,:) = indata(1:5,:)
  x(6,:)   = ingeodata(2,:)

  ! Magic geo is eps
  n(6) = n(1)

  do k=1,6
     xmin(k) = minval(x(k,:))
     xmax(k) = maxval(x(k,:))
  enddo

  ! test point
  !x0(1)=0.3391239
  !x0(2)=4.4351
  !x0(3)=log10(0.798065)
  !x0(4)=0.741015
  !x0(5)=1.564186
  !x0(6)=0.757284526

  print *,n(1:6)

  ! Rescale data
  do k=1,6
     x(k,:) = (x(k,:)-xmin(k))/(xmax(k)-xmin(k))*n(k)
  enddo
  do k=1,6
     xs0(k) = (x0(k)-xmin(k))/(xmax(k)-xmin(k))*n(k)
  enddo

  allocate(ipiv(ntot))
  allocate(b(ntot,6))
  allocate(a(ntot,ntot))

  do p=1,ntot
     do q=1,ntot
        a(p,q) = sqrt( sum((x(:,q)-x(:,p))**2) )**3
     enddo
  enddo
  b(:,:) = transpose(outdata(:,:))
  call DGESV(ntot,6,a,ntot,ipiv,b,ntot,info) 

  do k=1,6

     c0(k) = 0.0
     do p=1,ntot
        c0(k) = c0(k)+b(p,k)*sqrt( sum((xs0(:)-x(:,p))**2) )**3
     enddo

  enddo

10 format(1pe12.5)

end subroutine neo_rbf
