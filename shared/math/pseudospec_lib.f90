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

subroutine pseudospec_legendre(n,x,w,d1,dl)

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
  real, dimension(:,:), allocatable :: cl

  call gauss_legendre(-1.0,1.0,x,w,n)

  print *
  print *,'Integration weights (x,w,e,w_e)'
  print *

  do i=1,n
     print '(i2,2x,10(1pe14.7,1x))',i,x(i),w(i)
  enddo
  print *
  print *,'Sum(w) = ',sum(w)

  lwork = 2*n
  allocate(ipiv(n))
  allocate(work(lwork))

  allocate(c(n,n))
  allocate(cp(n,n))
  allocate(cl(n,n))

  do i=1,n
     do j=1,n
        call legendre(j-1,x(i),c(i,j),cp(i,j))
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

  print *
  print *,'Pseudospectral 1st derivative error on exp(-x)' 
  print *
  do i=1,n
     print '(i2,2x,10(1pe14.7,1x))',i,sum(d1(i,:)*exp(-x(:)))/(-exp(-x(i)))-1
  enddo
  print *
  print *,'Pseudospectral L error on exp(-x)' 
  print *
  do i=1,n
     print '(i2,2x,10(1pe14.7,1x))',i,sum(dl(i,:)*exp(-x(:)))/(exp(-x(i))*(1.0+2.0*x(i)-x(i)**2))-1
  enddo

end subroutine pseudospec_legendre


!---------------------------------------------------
! Generate Maxwell pseudospectral data on e=(0,emax) 
!
! e(n)    -> nodes 
! w(n)    -> weights with Sum(w)=1.0
! d1(n,n) -> d/dx
! d2(n,n) -> d^2/dx^2
!
! where e = emax x^2.
!
!               emax
!           2    /
! Note: -------- | de sqrt(e) exp(-e) h(x) = Sum w(i) h(i)
!       sqrt(pi) /
!                0
!----------------------------------------------------

subroutine pseudospec_maxwell(n,emax,e,w,d1,d2)

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: emax

  real, intent(out) :: e(n)
  real, intent(out) :: w(n)
  real, intent(out) :: d1(n,n)
  real, intent(out) :: d2(n,n)

  integer :: i,j
  real, dimension(:), allocatable :: x0
  real, dimension(:), allocatable :: w0
  real, dimension(:), allocatable :: mu
  real, dimension(:), allocatable :: nu
  real, dimension(:), allocatable :: alpha
  real, dimension(:), allocatable :: beta
  real, dimension(:), allocatable :: am
  real, dimension(:), allocatable :: bm
  real, dimension(:), allocatable :: am0
  real, dimension(:), allocatable :: bm0
  real, dimension(:), allocatable :: pn
  real, dimension(:), allocatable :: pnp

  integer :: lwork,info
  integer, dimension(:), allocatable :: ipiv
  real, dimension(:), allocatable :: work

  real, dimension(:,:), allocatable :: c
  real, dimension(:,:), allocatable :: cp
  real, dimension(:,:), allocatable :: cpp

  real, external :: pythag
  real :: pi,s,xp,dum

  integer :: n_plot=256

  pi = 3.14159265358979323846264338

  allocate(x0(n))
  allocate(w0(n))

  lwork = 2*n
  allocate(ipiv(n))
  allocate(work(lwork))
  allocate(c(n,n))
  allocate(cp(n,n))
  allocate(cpp(n,n))

  allocate(mu(0:2*n-1))
  allocate(nu(0:2*n-1))
  allocate(alpha(0:2*n-2))
  allocate(beta(0:2*n-2))
  allocate(am(n))
  allocate(bm(n))
  allocate(am0(n))
  allocate(bm0(n))
  allocate(pn(0:n-1))
  allocate(pnp(0:n-1))

  print *,'Pseudospectral Node-Weight-Derivative Generator'
  print *
  print '(t2,a,f3.1)','emax = ',emax

  select case (nint(emax))

  case(1)
     open(unit=1,file='nu1.txt',status='old')
  case(2)
     open(unit=1,file='nu2.txt',status='old')
  case(3)
     open(unit=1,file='nu3.txt',status='old')
  case(4)
     open(unit=1,file='nu4.txt',status='old')
  case(5)
     open(unit=1,file='nu5.txt',status='old')
  case(6)
     open(unit=1,file='nu6.txt',status='old')
  case(7)
     open(unit=1,file='nu7.txt',status='old')
  case(8)
     open(unit=1,file='nu8.txt',status='old')
  case(9)
     open(unit=1,file='nu9.txt',status='old')
  case(10)
     open(unit=1,file='nu10.txt',status='old')
  case default
     print '(a)','ERROR: (pseudospec_lib) emax out of bounds.'
  end select

  do j=0,2*n-1
     read(1,*) nu(j),mu(j)
  enddo
  close(1)

  alpha(0) = 0.5
  beta(0)  = 0.0
  do j=1,2*n-2
     alpha(j) = 0.5
     beta(j) = 1.0/(4.0*(4.0-1.0/j**2))
  enddo

  print *
  print *,'Modified moments from Maxima (nu)'
  print *

  do j=0,2*n-1
     print '(i2,2x,10(1pe23.16,1x))',j,nu(j)
  enddo

  print *
  print *,'Recursion Coefficients (a,b)'
  print *

  call orthog(n,nu,alpha,beta,am,bm)

  do j=1,n
     print '(i2,2x,10(1pe23.16,1x))',j-1,am(j),bm(j)
  enddo

  bm0 = bm
  am0 = am
  call gaucof(n,am0,bm0,nu(0),x0,w0)

  print *
  print *,'Integration weights (x,w,e,w_e)'
  print *

  ! Energy nodes and weights
  e = emax*x0(i)**2
  w = w0*emax**1.5*4/sqrt(pi)*x0**2

  do i=1,n
     print '(i2,2x,10(1pe14.7,1x))',i,x0(i),w0(i),e(i),w(i)
  enddo
  print *
  print *,'Sum(w) = ',sum(w)

  print *
  print *,'Discrete sum error check (integrals of monomials)'
  print *
  do j=0,2*n-2
     print '(i2,2x,1pe14.7)',j,sum(w0*x0**j)/mu(j)-1.0
  enddo

  open(unit=1,file='out',status='replace')
  open(unit=2,file='outp',status='replace')
  do i=1,n_plot
     xp = (i-1.0)/(n_plot-1.0)
     !s = nu(0)     
     do j=0,n-1
        call newpoly(j,xp,pn(j),pnp(j),dum,am,bm,n)
        !pn(j) = pn(j)/sqrt(s)/(j+1)
        !if (j < n-1) s = bm(j+2)*s
     enddo
     write(1,"(16(1pe14.7,1x))") xp,pn(:)
     write(2,"(16(1pe14.7,1x))") xp,pnp(:)
  enddo
  close(1)
  close(2)

  ! Derivative weights

  do j=1,n
     do i=1,n
        call newpoly(j-1,x0(i),c(i,j),cp(i,j),cpp(i,j),am,bm,n)
     enddo
  enddo

  call DGETRF(n,n,c,n,ipiv,info)
  call DGETRI(n,c,n,ipiv,work,lwork,info)
  ! c -> C^(-1)
  call DGEMM('N','N',n,n,n,1.0,cp,n,c,n,0.0,d1,n)
  ! d -> (Cp C^(-1))
  call DGEMM('N','N',n,n,n,1.0,cpp,n,c,n,0.0,d2,n)
  ! d2 -> (Cpp C^(-1))

  ! First derivative
  print *
  print *,'Pseudospectral 1st derivative error on exp(-x)' 
  print *
  do i=1,n
     print '(i2,2x,10(1pe14.7,1x))',i,sum(d1(i,:)*exp(-x0(:)))/(-exp(-x0(i)))-1
  enddo
  print *
  print *,'Pseudospectral 2nd derivative error on exp(-x)' 
  print *
  do i=1,n
     print '(i2,2x,10(1pe14.7,1x))',i,sum(d2(i,:)*exp(-x0(:)))/(exp(-x0(i)))-1
  enddo

end subroutine pseudospec_maxwell

subroutine orthog(n,nu,alpha,beta,a,b)  

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: alpha(2*n-1),beta(2*n-1),nu(2*n) 
  real, intent(inout) :: a(n),b(n)  
  integer :: k,l  
  real, dimension(:,:), allocatable :: sig

  allocate(sig(2*n+1,2*n+1))  

  do l=3,2*n  
     sig(1,l) = 0.0  
  enddo
  do l=2,2*n+1  
     sig(2,l) = nu(l-1)  
  enddo
  a(1) = alpha(1)+nu(2)/nu(1)  
  b(1) = 0.0  
  do k=3,n+1  
     do l=k,2*n-k+3  
        sig(k,l) = sig(k-1,l+1)+(alpha(l-1)-a(k-2))*sig(k-1,l)  & 
             -b(k-2)*sig(k-2,l)+beta(l-1)*sig(k-1,l-1)  
     enddo
     a(k-1) = alpha(k-1)+sig(k,k+1)/sig(k,k)-sig(k-1,k)/sig(k-1,k-1)  
     b(k-1) = sig(k,k)/sig(k-1,k-1)  
  enddo

  deallocate(sig)

end subroutine orthog

subroutine gaucof(n,a,b,mu0,x,w)  

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: mu0
  real, intent(inout) :: a(n),b(n),w(n),x(n)  
  real, dimension(:,:), allocatable :: z

  integer :: i,j  

  allocate(z(n,n))

  do i=1,n  
     if (i /= 1) b(i) = sqrt(b(i))  
     do j=1,n  
        if (i == j)then  
           z(i,j) = 1.0  
        else  
           z(i,j) = 0.0  
        endif
     enddo
  enddo

  call tqli(a,b,n,z)  
  call eigsrt(a,z,n)  

  do i=1,n  
     x(i) = a(i)  
     w(i) = mu0*z(1,i)**2  
  enddo

  deallocate(z)

end subroutine gaucof

real function pythag(a,b)  

  implicit none

  real a,b  
  real absa,absb  
  absa=abs(a)  
  absb=abs(b)  
  if(absa.gt.absb)then  
     pythag=absa*sqrt(1.0+(absb/absa)**2)  
  else  
     if(absb.eq.0.)then  
        pythag=0.  
     else  
        pythag=absb*sqrt(1.0+(absa/absb)**2)  
     endif
  endif

end function pythag

subroutine eigsrt(d,v,n)  

  implicit none

  integer, intent(in) :: n
  real, intent(inout) :: d(n)
  real, intent(inout) :: v(n,n)  
  integer :: i,j,k  
  real :: p

  do i=1,n-1  
     k = i  
     p = d(i)  
     do j=i+1,n  
        if (d(j) <= p)then  
!        if (d(j) >= p)then  
           k = j  
           p = d(j)  
        endif
     enddo
     if (k /= i) then  
        d(k) = d(i)  
        d(i) = p  
        do j=1,n  
           p      = v(j,i)  
           v(j,i) = v(j,k)  
           v(j,k) = p  
        enddo
     endif
  enddo

end subroutine eigsrt

subroutine tqli(d,e,n,z) 

  implicit none

  integer, intent(in) :: n  
  real, intent(inout) :: d(n),e(n)
  real, intent(inout) :: z(n,n) 
   
  integer :: i,iter,k,l,m  
  real :: b,c,dd,f,g,p,r,s,pythag

  do i=2,n  
     e(i-1)=e(i)  
  enddo
  e(n)=0.0  
  do l=1,n  
     iter=0  
1    do m=l,n-1  
        dd=abs(d(m))+abs(d(m+1))  
        if (abs(e(m))+dd.eq.dd) goto 2  
     enddo
     m=n  

2    if (m /= l) then  
        if (iter == 30) then
           print *,'too many iterations in tqli'
           stop
        endif
        iter=iter+1  
        g=(d(l+1)-d(l))/(2.0*e(l))  
        r=pythag(g,1.)  
        g=d(m)-d(l)+e(l)/(g+sign(r,g))  
        s=1.  
        c=1.  
        p=0.  
        do i=m-1,l,-1  
           f=s*e(i)  
           b=c*e(i)  
           r=pythag(f,g)  
           e(i+1)=r  
           if(r.eq.0.)then  
              d(i+1)=d(i+1)-p  
              e(m)=0.  
              goto 1  
           endif
           s=f/r  
           c=g/r  
           g=d(i+1)-p  
           r=(d(i)-g)*s+2.0*c*b  
           p=s*r  
           d(i+1)=g+p  
           g=c*r-b  
           !     Omit lines from here ...  
           do k=1,n  
              f=z(k,i+1)  
              z(k,i+1)=s*z(k,i)+c*f  
              z(k,i)=c*z(k,i)-s*f  
           enddo
           !     ... to here when finding only eigenvalues.  
        enddo
        d(l) = d(l)-p  
        e(l) = g  
        e(m) = 0.0  
        goto 1  
     endif
  enddo

end subroutine tqli

!-------------------------------------------------------
! p -> P(n,x)
! q -> P(n,x)'
!-------------------------------------------------------

subroutine legendre(n,x,p,q)

  implicit none

  integer, intent (in) :: n
  real, intent (in) :: x
  real, intent(out) :: p,q
  integer :: j
  real :: a,b
  real :: pmm,pm
  real :: qmm,qm

  if (n == 0) then
     p  = 1.0
     q  = 0.0
  else
     pmm = 0.0
     pm  = 1.0

     qmm = 0.0
     qm  = 0.0

     do j=0,n-1

        a = (2*j+1.0)/(j+1.0)
        b = j/(j+1.0)

        ! P(n,x)
        p = a*x*pm-b*pmm
        pmm = pm
        pm = p

        ! P'(n,x)
        q = a*x*qm-b*qmm + a*pmm
        qmm = qm
        qm  = q
     enddo
  endif

end subroutine legendre


!-------------------------------------------------------
! p -> Poly(n,x)
! q -> Poly(n,x)'
! r -> Poly(n,x)''
!-------------------------------------------------------
subroutine newpoly(n,x,p,q,r,am,bm,nm)

  implicit none

  integer, intent (in) :: n,nm
  real, intent (in) :: x
  real, intent (in) :: am(0:nm-1),bm(0:nm-1)
  real, intent(out) :: p,q,r
  integer :: j
  real :: pmm,pm
  real :: qmm,qm
  real :: rmm,rm

  if (n == 0) then
     p  = 1.0
     q  = 0.0
     r  = 0.0
  else
     pmm = 0.0
     pm  = 1.0
 
     qmm = 0.0
     qm  = 0.0

     rmm = 0.0
     rm  = 0.0

     do j=0,n-1

        ! P(n,x)
        p = (x-am(j))*pm-bm(j)*pmm
        pmm = pm
        pm = p

        ! P'(n,x)
        q = (x-am(j))*qm-bm(j)*qmm + pmm
        qmm = qm
        qm  = q

        ! P''(n,x)
        r = (x-am(j))*rm-bm(j)*rmm + 2*qmm
        rmm = rm
        rm  = r
     enddo
  endif

end subroutine newpoly

