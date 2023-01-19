!--------------------------------------------------
! Generate Legendre pseudospectral data on x=(-1,1) 
!
! x(n)    -> nodes 
! w(n)    -> weights with Sum(w)=2.
! d1(n,n) -> d/dx
! dl(n,n) -> d/dx (1-x^2) d/dx
!
!       1
!       /
! Note: | dx f(x) = Sum w(i) f(i)
!       /
!      -1
!---------------------------------------------------

subroutine pseudo_legendre(n,x,w,d1,dl)

  implicit none

  integer, intent(in) :: n

  real, intent(out) :: x(n)
  real, intent(out) :: w(n)
  real, intent(out) :: d1(n,n)
  real, intent(out) :: dl(n,n)

  integer :: i,j
  integer :: lwork,info
  integer, dimension(:), allocatable :: ipiv
  real, dimension(:), allocatable :: work
  real, dimension(:,:), allocatable :: c
  real, dimension(:,:), allocatable :: cp
  real, dimension(:,:), allocatable :: cpp
  real, dimension(:,:), allocatable :: cl

  integer, parameter :: print_flag=0

  call gauss_legendre(-1.0,1.0,x,w,n)

  lwork = 2*n
  allocate(ipiv(n))
  allocate(work(lwork))
  allocate(c(n,n))
  allocate(cp(n,n))
  allocate(cpp(n,n))
  allocate(cl(n,n))

  do i=1,n
     do j=1,n
        ! d/dxi
        call pseudo_rec_legendre(j-1,x(i),c(i,j),cp(i,j),cpp(i,j))
        ! L
        cl(i,j) = -(j-1)*j*c(i,j)
     enddo
  enddo

  call DGETRF(n,n,c,n,ipiv,info)
  call DGETRI(n,c,n,ipiv,work,lwork,info)
  ! c -> C^(-1)
  call DGEMM('N','N',n,n,n,1.0,cp,n,c,n,0.0,d1,n)
  ! d1 -> (Cp C^(-1))
  call DGEMM('N','N',n,n,n,1.0,cl,n,c,n,0.0,dl,n)
  ! dl -> (L C^(-1))

  if (print_flag == 1) then
     print *
     print *,'Integration weights (x,w,e,w_e)'
     print *
     do i=1,n
        print '(i2,2x,10(1pe14.7,1x))',i,x(i),w(i)
     enddo
     print *
     print *,'Sum(w) = ',sum(w)
     print *
     print *,'Pseudospectral 1st derivative error on exp(-x)' 
     print *
     do i=1,n
        print '(i2,2x,10(1pe14.7,1x))',i,sum(d1(i,:)*exp(-x(:)))/(-exp(-x(i)))-1
     enddo
     print *
     print *,'Pseudospectral L-derivative error on exp(-x)' 
     print *
     do i=1,n
        print '(i2,2x,10(1pe14.7,1x))',i,sum(dl(i,:)*exp(-x(:)))/(exp(-x(i))*(1.0+2.0*x(i)-x(i)**2))-1
     enddo

  endif

  deallocate(ipiv)
  deallocate(work)
  deallocate(c)
  deallocate(cp)
  deallocate(cpp)
  deallocate(cl)

end subroutine pseudo_legendre

!-------------------------------------------------------
! p = P(n,x) ; q = P'(n,x) ; r = P''(n,x)
!-------------------------------------------------------

subroutine pseudo_rec_legendre(n,x,p,q,r)

  implicit none

  integer, intent (in) :: n
  real, intent (in) :: x
  real, intent(out) :: p,q,r
  integer :: j
  real :: a,b
  real :: pmm,pm
  real :: qmm,qm
  real :: rmm,rm

  if (n == 0) then
     p = 1.0
     q = 0.0
     r = 0.0
  else
     pmm = 0.0 ; pm = 1.0
     qmm = 0.0 ; qm = 0.0
     rmm = 0.0 ; rm = 0.0
     
     do j=0,n-1

        a = (2*j+1.0)/(j+1.0)
        b = j/(j+1.0)

        ! P(n,x)
        p = a*x*pm-b*pmm
        pmm = pm ; pm = p

        ! P'(n,x)
        q = a*x*qm-b*qmm + a*pmm
        qmm = qm ; qm = q

        ! P''(n,x)
        r = a*x*rm-b*rmm + 2*a*qmm
        rmm = rm ; rm = r
        
     enddo
  endif

end subroutine pseudo_rec_legendre

real function p_legendre(n,x)

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: x
  real :: p,q,r

  call pseudo_rec_legendre(n,x,p,q,r)

  p_legendre = p

end function p_legendre
