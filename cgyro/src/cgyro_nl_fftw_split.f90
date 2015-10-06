!-----------------------------------------------------------------
! cgyro_nl_fftw_split.f90
!
! PURPOSE:
!  Evaluate nonlinear bracket with dealiased FFT.  It is natural 
!  to use the FFTW complex-to-real (c2r) transform to move to real 
!  space, compute the convolution (uv), and transform back to the 
!  spectral form using a real-to-complex transform (r2c).
!  
!  NOTE: Need to be careful with (p=-nr/2,n=0) component.
!-----------------------------------------------------------------

subroutine cgyro_nl_fftw_split(ij)

  use timer_lib
  use parallel_lib

  use cgyro_globals

  integer, intent(in) :: ij

  integer :: nx,ny
  integer :: nx0,ny0
  integer :: j,p,iexch

  complex :: f0,g0

  real, dimension(:,:), allocatable :: ux
  real, dimension(:,:), allocatable :: uy
  real, dimension(:,:), allocatable :: vx
  real, dimension(:,:), allocatable :: vy
  real, dimension(:,:), allocatable :: uv
  complex, dimension(:,:),allocatable :: fx
  complex, dimension(:,:),allocatable :: fy
  complex, dimension(:,:),allocatable :: gx
  complex, dimension(:,:),allocatable :: gy

  complex, dimension(:,:), allocatable :: fpack
  complex, dimension(:,:), allocatable :: gpack

  include 'fftw3.f03'


  ! 2D FFT lengths 
  nx0 = n_radial
  ny0 = 2*n_toroidal-1

  ! 3/2-rule for dealiasing the nonlinear product
  nx = (3*nx0)/2
  ny = (3*ny0)/2

  call timer_lib_in('nl_comm')
  allocate(fpack(n_radial,nv_loc*n_theta))
  allocate(gpack(n_radial,nv_loc*n_theta))
  fpack = 0.0
  gpack = 0.0
  iexch = 0
  do iv_loc=1,nv_loc
     do it=1,n_theta
        iexch = iexch+1
        do ir=1,n_radial
           fpack(ir,iexch) = h_x(ic_c(ir,it),iv_loc)
           gpack(ir,iexch) = psi(ic_c(ir,it),iv_loc)
        enddo
     enddo
  enddo
  call parallel_slib_f(fpack,f_nl)
  call parallel_slib_f(gpack,g_nl)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

  ! Allocate and deallocate these every time.
  allocate(fx(0:ny/2,0:nx-1))
  allocate(gx(0:ny/2,0:nx-1))
  allocate(fy(0:ny/2,0:nx-1))
  allocate(gy(0:ny/2,0:nx-1))

  allocate(ux(0:ny-1,0:nx-1))
  allocate(vx(0:ny-1,0:nx-1))
  allocate(uy(0:ny-1,0:nx-1))
  allocate(vy(0:ny-1,0:nx-1))
  allocate(uv(0:ny-1,0:nx-1))

  ! To avoid memory leak, create this plan every time :-(
  plan_c2r = fftw_plan_dft_c2r_2d(nx,ny,fx,ux,FFTW_MEASURE)
  plan_r2c = fftw_plan_dft_r2c_2d(nx,ny,ux,fx,FFTW_MEASURE)

  do j=1,nsplit

     fx = 0.0
     gx = 0.0
     fy = 0.0
     gy = 0.0

     ! Array mapping
     do in=1,n_toroidal
        iy = in-1
        do ir=1,n_radial
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx  
           f0 = f_nl(ir,j,in)
           g0 = g_nl(ir,j,in)
           fx(iy,ix) = i_c*p*f0
           gx(iy,ix) = i_c*p*g0
           fy(iy,ix) = i_c*iy*f0
           gy(iy,ix) = i_c*iy*g0
        enddo
     enddo

     if (kxfilter_flag == 1) then
        fx(:,-nx0/2+nx) = 0.0
        fy(:,-nx0/2+nx) = 0.0
        gx(:,-nx0/2+nx) = 0.0
        gy(:,-nx0/2+nx) = 0.0
     endif

     call fftw_execute_dft_c2r(plan_c2r,fx,ux)
     call fftw_execute_dft_c2r(plan_c2r,fy,uy)
     call fftw_execute_dft_c2r(plan_c2r,gx,vx)
     call fftw_execute_dft_c2r(plan_c2r,gy,vy)

     ! Poisson bracket in real space

     uv = (ux*vy-uy*vx)/(nx*ny)

     call fftw_execute_dft_r2c(plan_r2c,uv,fx)

     ! Must annhilate n=0,p=-nr/2
     fx(0,-nx0/2+nx) = 0.0       
     fx(0,0)         = 0.0       

     do ir=1,n_radial 
        ix = ir-1-nx0/2
        if (ix < 0) ix = ix+nx
        do in=1,n_toroidal
           iy = in-1
           g_nl(ir,j,in) = fx(iy,ix)
        enddo
     enddo

  enddo ! j
  call fftw_destroy_plan(plan_c2r)
  call fftw_destroy_plan(plan_r2c)

  deallocate(fx)
  deallocate(gx)
  deallocate(fy)
  deallocate(gy)

  deallocate(ux)
  deallocate(vx)
  deallocate(uy)
  deallocate(vy)
  deallocate(uv)

  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  call parallel_slib_r(g_nl,gpack)
  iexch = 0
  do iv_loc=1,nv_loc
     do it=1,n_theta
        iexch = iexch+1
        do ir=1,n_radial
           psi(ic_c(ir,it),iv_loc) = gpack(ir,iexch) 
        enddo
     enddo
  enddo
  deallocate(fpack)
  deallocate(gpack)
  call timer_lib_out('nl_comm')

  ! RHS -> -[f,g] = [f,g]_{r,-alpha}

  rhs(ij,:,:) = rhs(ij,:,:)+(q*rho/rmin)*(2*pi/length)*psi(:,:)

end subroutine cgyro_nl_fftw_split
