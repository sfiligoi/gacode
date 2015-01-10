!--------------------------------------------------------------
! Evaluate nonlinear bracket with direct dealiased convolution
!--------------------------------------------------------------

subroutine cgyro_nl_direct(ij)

  use timer_lib
  use parallel_lib

  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer, intent(in) :: ij
  integer :: nx,ny,nx0,ny0
  integer :: ix,ixp,iy,iyp
  integer :: ixpp,iypp
  integer :: ir,it,j,in
  complex, dimension(:,:), allocatable :: f
  complex, dimension(:,:), allocatable :: g
  complex, dimension(:,:), allocatable :: fg


  ny0 = n_toroidal-1
  nx0 = n_radial/2
  ny = int(1.5*ny0)+1
  nx = int(1.5*nx0)+1

  if (.not.allocated(pcyc)) then
     allocate(pcyc(-3*nx:3*nx-1))
     allocate(ncyc(-3*ny-1:3*ny+1))
     do ix=-nx,nx-1
        pcyc(ix-2*nx) = ix
        pcyc(ix) = ix
        pcyc(ix+2*nx) = ix
     enddo
     do iy=-ny,ny
        ncyc(iy-(2*ny+1)) = iy
        ncyc(iy) = iy
        ncyc(iy+(2*ny+1)) = iy
     enddo
  endif

  allocate( f(-nx:nx,-ny:ny) )
  allocate( g(-nx:nx,-ny:ny) )
  allocate(fg(-nx:nx,-ny:ny) )

  call timer_lib_in('comm_nl')
  call parallel_slib_f(h_x,f_nl)
  call parallel_slib_f(psi,g_nl)
  call timer_lib_out('comm_nl')

  call timer_lib_in('rhs_nl')
  do j=1,nsplit
     do it=1,n_theta 

        f = 0.0
        g = 0.0

        ! Array mapping
        do in=1,n_toroidal
           iy = in-1
           do ir=1,n_radial
              ic = ic_c(ir,it) 
              ix = ir-1-nx0
              f(ix,iy) = f_nl(ic,j,in)
              g(ix,iy) = g_nl(ic,j,in)
              f(-ix,-iy) = conjg(f(ix,iy))
              g(-ix,-iy) = conjg(g(ix,iy)) 
           enddo
           f(nx0,-iy) = f(-nx0,iy)
           g(nx0,-iy) = g(-nx0,iy) 
        enddo

        ! Reality
        do ix=-nx0,nx0
           do iy=1,ny0
              f(-ix,-iy) = conjg(f(ix,iy))
              g(-ix,-iy) = conjg(g(ix,iy))
           enddo
           if (it == 7 .and. j == 2 .and. i_proc == 0) then
              print '(10(1pe10.3,1x,1pe10.3,2x))',f(ix,-ny0:ny0)
           endif
        enddo

        fg = (0.0,0.0)
        do ix=-nx0,nx0-1
           do ixp=-nx,nx-1
              ixpp = pcyc(ix-ixp)
              do iy=0,ny0
                 do iyp=-ny+iy,ny
                    iypp = ncyc(iy-iyp)
                    fg(ix,iy) = fg(ix,iy)+f(ixpp,iypp)*g(ixp,iyp)*(iypp*ixp-iyp*ixpp)
                 enddo
              enddo
           enddo
        enddo

        do ir=1,n_radial
           ic = ic_c(ir,it) 
           ix = ir-1-nx0
           do in=1,n_toroidal
              iy = in-1
              g_nl(ic,j,in) = fg(ix,iy)
           enddo
        enddo

        if (it == 7 .and. j == 2 .and. i_proc == 0) then
           do ir=1,n_radial
              ix = ir-1-nx0
              print '(10(1pe11.4,1x,1pe11.4,3x))',fg(ix,0:ny0)
           enddo
        endif

     enddo ! it
  enddo ! j
  stop
  call timer_lib_out('rhs_nl')

  call timer_lib_in('comm_nl')
  call parallel_slib_r(g_nl,psi)
  call timer_lib_out('comm_nl')

  deallocate( f)
  deallocate( g)
  deallocate(fg)

  ! RHS -> -[f,g] = (n'' p' - n' p'') f'' g' 

  rhs(ij,:,:) = rhs(ij,:,:)+(q*rho/rmin)*(2*pi/length)*psi(:,:)

end subroutine cgyro_nl_direct

!--------------------------------------------------------------
! Evaluate nonlinear bracket with dealiased FFT 
!--------------------------------------------------------------

subroutine cgyro_nl_fftw(ij)

  use, intrinsic :: iso_c_binding

  use timer_lib
  use parallel_lib

  use cgyro_globals

  integer, intent(in) :: ij

  integer :: nx,ny
  integer :: nx0,ny0
  integer :: i,j,p

  complex :: f0,g0

  real(C_DOUBLE), dimension(:,:), allocatable :: ux
  real(C_DOUBLE), dimension(:,:), allocatable :: uy
  real(C_DOUBLE), dimension(:,:), allocatable :: vx
  real(C_DOUBLE), dimension(:,:), allocatable :: vy
  real(C_DOUBLE), dimension(:,:), allocatable :: uv
  complex(C_DOUBLE_COMPLEX), dimension(:,:),allocatable :: fx
  complex(C_DOUBLE_COMPLEX), dimension(:,:),allocatable :: fy
  complex(C_DOUBLE_COMPLEX), dimension(:,:),allocatable :: gx
  complex(C_DOUBLE_COMPLEX), dimension(:,:),allocatable :: gy

  include 'fftw3.f03'

  type(C_PTR) :: plan_r2c
  type(C_PTR) :: plan_c2r

  nx0 = n_radial
  ny0 = 2*n_toroidal-1
  ! 3/2-rule for dealiasing
  nx = (3*nx0)/2+3
  ny = (3*ny0)/2+3

  allocate(fx(0:ny/2,0:nx-1))
  allocate(gx(0:ny/2,0:nx-1))
  allocate(fy(0:ny/2,0:nx-1))
  allocate(gy(0:ny/2,0:nx-1))

  allocate(ux(0:ny-1,0:nx-1))
  allocate(vx(0:ny-1,0:nx-1))
  allocate(uy(0:ny-1,0:nx-1))
  allocate(vy(0:ny-1,0:nx-1))

  allocate(uv(0:ny-1,0:nx-1))

  plan_c2r = fftw_plan_dft_c2r_2d(nx,ny,fx,ux,FFTW_ESTIMATE)
  plan_r2c = fftw_plan_dft_r2c_2d(nx,ny,uv,fx,FFTW_ESTIMATE)

  call timer_lib_in('comm_nl')
  call parallel_slib_f(h_x,f_nl)
  call parallel_slib_f(psi,g_nl)
  call timer_lib_out('comm_nl')

  call timer_lib_in('rhs_nl')
  do j=1,nsplit
     do it=1,n_theta 

        fx = 0.0
        gx = 0.0
        fy = 0.0
        gy = 0.0

        ! Array mapping
        do in=1,n_toroidal
           iy = in-1
           do ir=1,n_radial
              ic = ic_c(ir,it) 
              p  = ir-1-nx0/2
              ix = -p
              if (ix < 0) ix = ix+nx  
              f0 = conjg(f_nl(ic,j,in))
              g0 = conjg(g_nl(ic,j,in))
              if (it == 7 .and. j == 2 .and. i_proc == 0) then
                 if (iy == 0 .and. p == 0) print *,f0
              endif
              fx(iy,ix) = i_c*p*f0
              gx(iy,ix) = i_c*p*g0
              fy(iy,ix) = i_c*iy*f0
              gy(iy,ix) = i_c*iy*g0
           enddo
        enddo

        !do i=1,nx
        !   print '(10(1pe12.5,2x))', fr(i,:)
        !enddo

        call fftw_execute_dft_c2r(plan_c2r,fx,ux)
        call fftw_execute_dft_c2r(plan_c2r,fy,uy)
        call fftw_execute_dft_c2r(plan_c2r,gx,vx)
        call fftw_execute_dft_c2r(plan_c2r,gy,vy)

        ! Poisson bracket in real space

        uv = -(ux*vy-uy*vx)/(nx*ny)

        fx = 0.0
        call fftw_execute_dft_r2c(plan_r2c,uv,fx)


        if (it == 7 .and. j == 2 .and. i_proc == 0) then
           print *
           do ir=1,n_radial
              ix = -(ir-1-nx0/2)
              if (ix < 0) ix=ix+nx
              print '(10(1pe11.4,1x,1pe11.4,3x))',conjg(fx(0:ny0/2,ix))
           enddo
        endif

        do ir=1,n_radial
           ic = ic_c(ir,it) 
           p = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx
           do in=1,n_toroidal
              iy = in-1
              g_nl(ic,j,in) = fx(iy,ix)
           enddo
        enddo


     enddo ! it
  enddo ! j
  stop
  call timer_lib_out('rhs_nl')

  call timer_lib_in('comm_nl')
  call parallel_slib_r(g_nl,psi)
  call timer_lib_out('comm_nl')

  call fftw_destroy_plan(plan_c2r)
  call fftw_destroy_plan(plan_r2c)

  deallocate(fx)
  deallocate(fy)
  deallocate(gx)
  deallocate(gy)
  deallocate(ux)
  deallocate(uy)
  deallocate(vx)
  deallocate(vy)
  deallocate(uv)

  ! RHS -> -[f,g] = (n'' p' - n' p'') f'' g' 

  rhs(ij,:,:) = rhs(ij,:,:)+(q*rho/rmin)*(2*pi/length)*psi(:,:)

end subroutine cgyro_nl_fftw
