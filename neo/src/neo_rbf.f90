subroutine neo_rbf(x0,c0)

  implicit none

  real, dimension(6), intent(in) :: x0
  real, dimension(6), intent(out) :: c0
  real, dimension(6) :: xs0

  integer :: p,k
  integer :: stat

  integer :: ntot
  real, dimension(:,:), allocatable :: x
  real, dimension(:,:), allocatable :: b

  real, dimension(2,6) :: s
 
  allocate(x(6,ntot))
  allocate(b(ntot,6))

  open(unit=1,file='out.pneo.scale',status='old',iostat=stat)
  if (stat /= 0) then
     call neo_error('ERROR: (neo_rbf) out.pneo.scale not available')
     return
  endif
  read(1,*) ntot
  do k=1,6
     read(1,*) s(:,k)  
  enddo
  close(1)

  open(unit=1,file='out.pneo.rbf',status='old',iostat=stat)
  if (stat /= 0) then
     call neo_error('ERROR: (neo_rbf) out.pneo.rbf not available')
     return
  endif
  read(1,*) x(:,:)
  read(1,*) b(:,:)
  close(1)

  ! Rescale data
  do k=1,6
     xs0(k) = s(1,k)*x0(k)+s(2,k)
  enddo

 do k=1,6
     c0(k) = 0.0
     do p=1,ntot
        c0(k) = c0(k)+b(p,k)*sqrt( sum((xs0(:)-x(:,p))**2) )**3
     enddo
  enddo

end subroutine neo_rbf
