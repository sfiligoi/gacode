subroutine vgen_newton(ia,ib,w0_newton)

  use vgen_globals
  use EXPRO_interface

  implicit none

  integer, intent(in) :: ia,ib
  real, intent(inout), dimension(ib-ia+1) :: w0_newton

  integer :: ierr
  integer :: n
  integer :: i,j
  integer :: i_newt
  integer, parameter :: n_newt_max=20
  real,parameter     :: tol_newt=5.0e-3
  integer :: ii
  integer, dimension(:), allocatable :: ipiv
  real, dimension(:), allocatable :: x,x0
  real, dimension(:), allocatable :: xp,xp0
  real, dimension(:), allocatable :: f,f0
  real, dimension(:), allocatable :: r
  real, dimension(:), allocatable :: b
  real, dimension(:,:), allocatable :: jac
  real :: dx 
  real :: er=0.0

  n = ib-ia+1
  allocate(ipiv(n))
  allocate(r(n))
  allocate(x0(n))
  allocate(x(n))
  allocate(xp0(n))
  allocate(xp(n))
  allocate(f0(n))
  allocate(f(n))
  allocate(b(n))
  allocate(jac(n,n))

  open(unit=1,file='vgenerror.out',status='replace')
  
  ! Get the initial guess x0 -- use the weak rotation limit (rotation_model=1)
  do i=1,n
     ii = i+ia-1
     r(i) = EXPRO_rmin(ii)
     call vgen_compute_neo(ii,f0(i),1,er,x0(i),xp0(i)) 
     x0(i) = f0(i)/(EXPRO_rmaj(ii)+EXPRO_rmin(ii))
  enddo
  
  print *, 'iteration=', 0
  !print *,x0(:)

  dx = 0.01*maxval(abs(x0))

  i_newt = 1
  do

     call bound_deriv(xp0,x0,r,n)
     do i=1,n
        call vgen_compute_neo(i+ia-1,f0(i),2,er,x0(i),xp0(i))
     enddo

     jac(:,:) = 0.0
     do j=1,n

        ! Initial guess plus tweak at i=j:
        x(:) = x0(:)
        x(j) = x(j)+dx
        call bound_deriv(xp,x,r,n)

        f(:) = 0.0
        do i=1,n
           if (abs(xp(i)-xp0(i)) > 1e-6) then
              call vgen_compute_neo(i+ia-1,f(i),2,er,x(i),xp(i))
              ! J_ij = df_i/dx_j
              jac(i,j) = (f(i)-f0(i))/dx
           endif
        enddo ! i
     enddo ! j

     ! LAPACK matrix factorization into L/U components
     call DGETRF(n,n,jac,n,ipiv,ierr)

     ! LAPACK matrix solve jac*delta_x=b
     b(:) = -f0(:)
     call DGETRS('N',n,1,jac,n,ipiv,b,n,ierr)

     ! Update guess and check convergence
     x0(:) = x0(:)+b(:)
     print *, 'iteration=',i_newt, 'error=', maxval(abs(b))
     !print *, x0(:)
     open(unit=1,file='vgenerror.out',status='old',position='append')
     write(1,10) i_newt,maxval(abs(b)), x0(:)
     close(1)
     if(maxval(abs(b)) < tol_newt) exit
     
     i_newt = i_newt + 1
     if(i_newt > n_newt_max) exit

  enddo
  
  if(i_newt > n_newt_max) then
     if(i_proc == 0) then
        print *, 'vgen_newton failed to converge'
     endif
     call MPI_finalize(i_err)
     stop
  endif

  w0_newton = x0
  
  deallocate(ipiv)
  deallocate(r)
  deallocate(x0)
  deallocate(x)
  deallocate(xp0)
  deallocate(xp)
  deallocate(f0)
  deallocate(f)
  deallocate(b)
  deallocate(jac)
  
10 format(t2,i2,10(1pe16.8,1x))
  
end subroutine vgen_newton
