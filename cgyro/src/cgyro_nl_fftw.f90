!-----------------------------------------------------------------
! cgyro_nl_fftw.f90
!
! PURPOSE:
!  Evaluate nonlinear bracket with dealiased FFT.  It is natural 
!  to use the FFTW complex-to-real (c2r) transform to move to real 
!  space, compute the convolution (uv), and transform back to the 
!  spectral form using a real-to-complex transform (r2c).
!  
! NOTE: Need to be careful with (p=-nr/2,n=0) component.
!-----------------------------------------------------------------

subroutine cgyro_nl_fftw(ij)

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  integer :: ix,iy
  integer :: ir,it,in
  integer :: j,p,iexch

  complex :: f0,g0

  include 'fftw3.f03'

  call timer_lib_in('nl_comm')

!$omp parallel do private(iv_loc,it,iexch,ir)
  do iv_loc=1,nv_loc
     iexch = (iv_loc-1)*n_theta
     do it=1,n_theta
        iexch = iexch+1
        do ir=1,n_radial
           fpack(ir,iexch) = h_x(ic_c(ir,it),iv_loc)
        enddo
     enddo
  enddo

  do iexch=nv_loc*n_theta+1,nsplit*n_toroidal
     fpack(:,iexch) = (0.0,0.0)
  enddo

  call parallel_slib_f_i_start(fpack,f_nl,f_com_req)

!$omp parallel do private(iv_loc,it,iexch,ir)
  do iv_loc=1,nv_loc
     iexch = (iv_loc-1)*n_theta
     do it=1,n_theta
        iexch = iexch+1
        do ir=1,n_radial
           gpack(ir,iexch) = psi(ic_c(ir,it),iv_loc)
        enddo
     enddo
  enddo

  do iexch=nv_loc*n_theta+1,nsplit*n_toroidal
     gpack(:,iexch) = (0.0,0.0)
  enddo


  call parallel_slib_f_i_start(gpack,g_nl,g_com_req)
  call parallel_slib_f_i_complete2(f_com_req,g_com_req)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

!$omp parallel private(fx,gx,fy,gy,in,iy,ir,p,ix,f0,g0,uv,ux,uy,vx,vy)
!$omp do schedule(dynamic,1)
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
  call parallel_slib_r_nc(g_nl,gpack)
!$omp parallel do private(iv_loc,it,iexch,ir)
  do iv_loc=1,nv_loc
     iexch = (iv_loc-1)*n_theta
     do it=1,n_theta
        iexch = iexch+1
        do ir=1,n_radial
           psi(ic_c(ir,it),iv_loc) = gpack(ir,iexch) 
        enddo
     enddo
  enddo
  call timer_lib_out('nl_comm')

  ! RHS -> -[f,g] = [f,g]_{r,-alpha}

  rhs(:,:,ij) = rhs(:,:,ij)+(q*rho/rmin)*(2*pi/length)*psi(:,:)

end subroutine cgyro_nl_fftw
