subroutine gauss_integ(xmin,xmax,func,order,n_subdiv,answer)

  implicit none

  !--------------------------------------------
  ! Subroutine arguments
  !
  real, intent(in) :: xmin
  real, intent(in) :: xmax
  real :: func
  integer, intent(in) :: order
  integer, intent(in) :: n_subdiv
  !
  real, intent(inout) :: answer
  !--------------------------------------------

  !--------------------------------------------
  ! Internal variables
  !
  integer :: n_node
  integer :: i
  integer :: j
  integer :: p
  !
  real :: dx
  !
  real, dimension(:), allocatable :: x
  real, dimension(:), allocatable :: w
  real, dimension(:), allocatable :: x0
  real, dimension(:), allocatable :: w0
  !--------------------------------------------

  allocate(x0(order))
  allocate(w0(order))
  call gauss_legendre(0.0,1.0,x0,w0,order)     

  dx = (xmax-xmin)/n_subdiv

  n_node = n_subdiv*order

  allocate(x(n_node))
  allocate(w(n_node))

  do i=1,n_subdiv
     do j=1,order
        p = (i-1)*order+j
        x(p) = xmin+((i-1)+x0(j))*dx
        w(p) = w0(j)*dx
     enddo
  enddo

  answer = 0.0
  do p=1,n_node
     answer = answer+w(p)*func(x(p))
  enddo

  deallocate(x)
  deallocate(w)
  deallocate(x0)
  deallocate(w0)

end subroutine  gauss_integ
