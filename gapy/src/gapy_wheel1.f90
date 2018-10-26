subroutine wheel1(nr,nth,nn,nx,ny,c,f)

  implicit none

  integer, intent(in) :: nr,nth,nn,nx,ny
  double complex, intent(in) :: c(0:nr-1,0:nth-1,0:nn-1)
  double precision, intent(inout) :: f(0:nx-1,0:ny-1)

  double complex, dimension(:,:), allocatable :: epx,cj
  double precision, dimension(:), allocatable :: eny
  double precision, dimension(:), allocatable :: yj,th

  integer :: i,j,n,p,jc
  double precision :: pi,fsum,xi
  double complex :: ic

  ! f2py intent(in) nr
  ! f2py intent(in) nth  
  ! f2py intent(in) nn
  ! f2py intent(in) nx
  ! f2py intent(in) ny
  ! f2py intent(in) c
  ! f2py intent(in,out) f

  allocate(epx(0:nr-1,0:nx-1))
  allocate(eny(0:nn-1))
  allocate(yj(0:ny-1))
  allocate(th(0:nth))
  allocate(cj(0:nr-1,0:nn-1))

  ic = (0d0,1d0)
  pi = atan(1d0)*4d0

  do i=0,nx-1
     xi = i*2*pi/(nx-1)
     do p=0,nr-1    
        epx(p,i) = exp(ic*(p-nr/2)*xi)
     enddo
  enddo

  do j=0,nth
     th(j) = j*2*pi/nth-pi
  enddo
  do j=0,ny-1
     yj(j) = j*pi/(ny-1)-pi
  enddo

  ! factor of 1/2 for n=0
  eny(:) = 1.0
  eny(0) = 0.5 
  f = 0.0

  do j=0,ny-1
     do jc=0,nth-1
        if (yj(j) >= th(jc) .and. yj(j) <= th(jc+1)) then 
           cj(:,:) = ( c(:,jc+1,:)*(yj(j)-th(jc)) &
                + c(:,jc,:)*(th(jc+1)-yj(j)) )/(th(1)-th(0))
           exit
        endif
     enddo
     do i=0,nx-1
        fsum = 0.0
        do n=0,nn-1
           do p=0,nr-1
              fsum = fsum+real(cj(p,n)*epx(p,i)*eny(n))
           enddo
        enddo
        f(i,j) = fsum
     enddo
  enddo

  deallocate(epx,eny,yj,th,cj)

end subroutine wheel1
