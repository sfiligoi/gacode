subroutine neo_rbf(x0,c0)

  implicit none

  real, dimension(6), intent (in) :: x0
  real, dimension(6), intent (out) :: c0
  
  integer :: p,k
  integer, dimension(9) :: n

  integer :: ntot
  real, dimension(:,:), allocatable :: indata
  real, dimension(:,:), allocatable :: ingeodata
  real, dimension(:,:), allocatable :: outdata
  real, dimension(:,:), allocatable :: x
  real, dimension(:), allocatable :: w
  real, dimension(16) :: dum
  real, dimension(6) :: xmin,xmax,xs0
  character(len=10), dimension(6) :: file
  integer :: stat

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

  file(1) = 'node_cne'
  file(2) = 'node_cte'
  file(3) = 'node_cni1'
  file(4) = 'node_cti1'
  file(5) = 'node_cni2'
  file(6) = 'node_cti2'

  do k=1,6
     open(unit=1,file=trim(file(k)),status='old',iostat=stat)
     if (stat /= 0) then
        call neo_error('ERROR: (NEO) RBF files not available')
        return
     endif
     do p=1,ntot
        read(1,*) w(p)
     enddo
     close(1) 

     c0(k) = 0.0
     do p=1,ntot
        c0(k) = c0(k)+w(p)*sqrt( sum((xs0(:)-x(:,p))**2) )**3
     enddo
     !print *,file(k),c0(k)
  enddo

10 format(1pe12.5)

end subroutine neo_rbf
