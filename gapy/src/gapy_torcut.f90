subroutine torcut(nr,nth,nn,nx,ny,c,f)

  implicit none

  double precision :: q=2.0
  integer :: m=8
  integer, intent(in) :: nr,nn,nth,nx,ny
  double complex, intent(in) :: c(0:nr-1,0:nth-1,0:nn-1)
  double precision, intent(inout) :: f(0:nx-1,0:ny-1)

  double complex, dimension(:,:), allocatable :: epx,eny,cj
  double complex, dimension(:,:,:), allocatable :: cx
  double precision, dimension(:), allocatable :: xi,yj,th

  integer :: i,j,n,p,jc,pp
  double precision :: pi,fsum
  double complex :: ic

  ! f2py intent(in) nr
  ! f2py intent(in) nth
  ! f2py intent(in) nn
  ! f2py intent(in) nx
  ! f2py intent(in) ny
  ! f2py intent(in) c
  ! f2py intent(in,out) f
  ! f2py intent(in) q
  ! f2py intent(in) m

  allocate(epx(0:nr-1,0:nx-1))
  allocate(eny(0:nn-1,0:ny-1))
  allocate(xi(0:nx-1))
  allocate(yj(0:ny-1))
  allocate(th(0:nth))
  allocate(cj(0:nr-1,0:nn-1))
  allocate(cx(0:nr-1,0:nth,0:nn-1))

  ic = (0d0,1d0)
  pi = atan(1d0)*4d0

  do j=0,nth
     th(j) = j*2*pi/nth
  enddo

  do i=0,nx-1
     xi(i) = i*2*pi/(nx-1)
     do p=0,nr-1    
        epx(p,i) = exp(ic*(p-nr/2)*xi(i))
     enddo
  enddo

  do j=0,ny-1
     yj(j) = j*2*pi/(ny-1)
     do n=0,nn-1    
        eny(n,j) = exp(ic*n*q*yj(j))
     enddo
  enddo

  ! Extended coefficients
  cx(:,0:nth-1,:) = c(:,:,:) 
  do n=0,nn-1
     do p=0,nr-1
        pp = p+n*m
        if (pp > nr-1) pp = pp-nr
        cx(p,nth,n) = c(pp,0,n)*exp(-2*pi*ic*n*q)
     enddo
  enddo

  ! factor of 1/2 for n=0
  eny(0,:) = 0.5*eny(0,:)

  f = 0.0

!!$omp parallel do private(j,i,fsum,n,p)
  do j=0,ny-1

     do jc=0,nth-1
        if (yj(j) >= th(jc) .and. yj(j) <= th(jc+1)) then 
           cj(:,:) = (cx(:,jc+1,:)*(yj(j)-th(jc))+cx(:,jc,:)*(th(jc+1)-yj(j)))/(th(1)-th(0))
           if (jc+1 > nth) print *,'ERROR1' 
           exit
        endif
     enddo

     do i=0,nx-1
        fsum = 0.0
        do n=0,nn-1
           do p=0,nr-1
              fsum = fsum+real(cj(p,n)*epx(p,i)*eny(n,j)*exp(ic*n*m*xi(i)*yj(j)/(2*pi)))
           enddo
        enddo
        f(i,j) = fsum
     enddo

  enddo

  deallocate(epx,eny,xi,yj,th,cj)

end subroutine torcut
