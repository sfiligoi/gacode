!-----------------------------------------------------------------
! cgyro_nl_fftw.gpu.f90
!
! PURPOSE:
!  Evaluate nonlinear bracket with dealiased FFT.  It is natural 
!  to use the FFTW complex-to-real (c2r) transform to move to real 
!  space, compute the convolution (uv), and transform back to the 
!  spectral form using a real-to-complex transform (r2c).
!  
!  NOTE: Need to be careful with (p=-nr/2,n=0) component.
!-----------------------------------------------------------------

subroutine cgyro_nl_fftw(ij)

  use timer_lib
  use parallel_lib

  use cgyro_globals
  implicit none

  integer, intent(in) :: ij
  integer :: j,p,iexch
  integer :: it,ir,in,ix,iy

  complex :: f0,g0
  complex, dimension(:,:), allocatable :: fpack
  complex, dimension(:,:), allocatable :: gpack

  include 'fftw3.f03'

  call timer_lib_in('nl_comm')

  allocate(fpack(n_radial,nv_loc*n_theta))
  allocate(gpack(n_radial,nv_loc*n_theta))
  ! fpack = 0.0
  ! gpack = 0.0
  iexch = 0
!$omp  parallel do collapse(2) &
!$omp& private(iv_loc,it,ir,iexch,ic_loc)
  do iv_loc=1,nv_loc
     do it=1,n_theta
        do ir=1,n_radial
        ! iexch = iexch+1
        iexch = it + (iv_loc-1)*n_theta
 
           ! ic_loc = ic_c(ir,it)
           ic_loc = it + (ir-1)*n_theta
           fpack(ir,iexch) = h_x(ic_loc,iv_loc)
           gpack(ir,iexch) = psi(ic_loc,iv_loc)
        enddo
     enddo
  enddo
  call parallel_slib_f(fpack,f_nl)
  call parallel_slib_f(gpack,g_nl)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

!$omp parallel private(fx,gx,fy,gy,in,iy,ir,p,ix,f0,g0,uv,ux,uy,vx,vy)
!$omp do
  do j=1,nsplit

     fx = 0.0
     gx = 0.0
     fy = 0.0
     gy = 0.0

     ! Array mapping
     do ir=1,n_radial
        p  = ir-1-nx0/2
        ix = p
        if (ix < 0) ix = ix+nx  
        do in=1,n_toroidal
           iy = in-1
           f0 = i_c*f_nl(ir,j,in)
           g0 = i_c*g_nl(ir,j,in)
           fx(iy,ix) = p*f0
           gx(iy,ix) = p*g0
           fy(iy,ix) = iy*f0
           gy(iy,ix) = iy*g0
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

     ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
     ! that will be filtered in the main time-stepping loop

     do ir=1,n_radial 
        ix = ir-1-nx0/2
        if (ix < 0) ix = ix+nx
        do in=1,n_toroidal
           iy = in-1
           g_nl(ir,j,in) = fx(iy,ix)
        enddo
     enddo

  enddo ! j
!$omp end do
!$omp end parallel

  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  call parallel_slib_r(g_nl,gpack)
  iexch = 0
!$omp  parallel do collapse(2) &
!$omp& private(iv_loc,it,ir,iexch,ic_loc)
  do iv_loc=1,nv_loc
     do it=1,n_theta
        do ir=1,n_radial
        ! iexch = iexch+1
        iexch = it + (iv_loc-1)*n_theta
           ! ic_loc = ic_c(ir,it)
           ic_loc = it + (ir-1)*n_theta
           psi(ic_loc,iv_loc) = gpack(ir,iexch) 
        enddo
     enddo
  enddo
  deallocate(fpack)
  deallocate(gpack)
  call timer_lib_out('nl_comm')

  ! RHS -> -[f,g] = [f,g]_{r,-alpha}

  rhs(ij,:,:) = rhs(ij,:,:)+(q*rho/rmin)*(2*pi/length)*psi(:,:)

end subroutine cgyro_nl_fftw
