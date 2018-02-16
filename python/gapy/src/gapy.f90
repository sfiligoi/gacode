subroutine p_cub_spline_deriv(x,y,n,yp)

  integer, intent(in) :: n
  double precision, intent(in), dimension(n) :: x,y
  double precision, intent(inout), dimension(n) :: yp
  
!f2py intent(in) n
!f2py intent(in) x
!f2py intent(in) y
!f2py intent(out) yp

  call cub_spline_deriv(x,y,n,yp)
  
end subroutine p_cub_spline_deriv

subroutine p_bound_deriv(df,f,r,n)

  integer, intent(in) :: n

  double precision, intent(inout), dimension(n) :: df
  double precision, intent(in), dimension(n) :: f
  double precision, intent(in),dimension(n) :: r

!f2py intent(in) n
!f2py intent(in) f
!f2py intent(in) r
!f2py intent(out) df

  call bound_deriv(df,f,r,n)
  
end subroutine p_bound_deriv

subroutine realfluct(nr,nn,nx,ny,c,f)

  implicit none
  
  integer, intent(in) :: nr,nn,nx,ny
  double complex, intent(in) :: c(0:nr-1,0:nn-1)
  double precision, intent(inout) :: f(0:nx-1,0:ny-1)

  double complex :: epx(0:nx-1,0:nr-1)
  double complex :: eny(0:ny-1,0:nn-1)

  integer :: i,j,n,p
  double precision :: pi,xi,yj
  double complex :: ic

  ! f2py intent(in) nr
  ! f2py intent(in) nn
  ! f2py intent(in) nx
  ! f2py intent(in) ny
  ! f2py intent(in) c
  ! f2py intent(in,out) f

  ic = (0d0,1d0)
  pi = atan(1d0)*4d0

  do i=0,nx-1
     xi = i*2*pi/(nx-1)
     do p=0,nr-1    
        epx(i,p) = exp(ic*(p-nr/2)*xi)
     enddo
  enddo

  do j=0,ny-1
     yj = j*2*pi/(ny-1)
     do n=0,nn-1    
        eny(j,n) = exp(-ic*n*yj)
     enddo
  enddo

  ! factor of 1/2 for n=0
  eny(:,0) = 0.5*eny(:,0)

  f = 0.0
  do n=0,nn-1
     do p=0,nr-1 
        do j=0,ny-1
           do i=0,nx-1
              f(i,j) = f(i,j)+real(c(p,n)*epx(i,p)*eny(j,n))
           enddo
        enddo
     enddo
  enddo

end subroutine realfluct
