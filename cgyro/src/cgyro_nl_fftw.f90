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

subroutine cgyro_nl_fftw_comm1
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m
  integer :: iexch

!$omp parallel do private(iv_loc,it,iexch,ir)
  do iv_loc_m=1,nv_loc
     iexch = (iv_loc_m-1)*n_theta
     do it=1,n_theta
        iexch = iexch+1
        do ir=1,n_radial
           fpack(ir,iexch) = h_x(ic_c(ir,it),iv_loc_m)
        enddo
     enddo
  enddo

  do iexch=nv_loc*n_theta+1,nsplit*n_toroidal
     fpack(:,iexch) = (0.0,0.0)
  enddo

  call parallel_slib_f_nc(fpack,f_nl)
end subroutine cgyro_nl_fftw_comm1

subroutine cgyro_nl_fftw_comm2
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m
  integer :: iexch

!$omp parallel do private(iv_loc,it,iexch,ir)
  do iv_loc_m=1,nv_loc
     iexch = (iv_loc_m-1)*n_theta
     do it=1,n_theta
        iexch = iexch+1
        do ir=1,n_radial
           gpack(ir,iexch) = psi(ic_c(ir,it),iv_loc_m)
        enddo
     enddo
  enddo

  do iexch=nv_loc*n_theta+1,nsplit*n_toroidal
     gpack(:,iexch) = (0.0,0.0)
  enddo

  call parallel_slib_f_nc(gpack,g_nl)
end subroutine cgyro_nl_fftw_comm2

! NOTE: call cgyro_nl_fftw_comm1 before cgyro_nl_fftw
subroutine cgyro_nl_fftw(ij)

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  integer :: ix,iy
  integer :: ir,it,in
  integer :: j,p,iexch
  integer :: i_omp

  complex :: f0,g0

  integer, external :: omp_get_thread_num

  include 'fftw3.f03'

  if (is_staggered_comm_2) then ! stagger comm2, to load ballance network traffic
    call timer_lib_in('nl_comm')
    call cgyro_nl_fftw_comm2
    call timer_lib_out('nl_comm')
  endif

  call timer_lib_in('nl')

!$omp parallel private(in,iy,ir,p,ix,f0,i_omp,j)
!$omp do schedule(dynamic,1)
  do j=1,nsplit
     i_omp = omp_get_thread_num()+1

     fx(:,:,i_omp) = 0.0
     fy(:,:,i_omp) = 0.0

     ! Array mapping
     do ir=1,n_radial
        p  = ir-1-nx0/2
        ix = p
        if (ix < 0) ix = ix+nx
        do in=1,n_toroidal
           iy = in-1
           f0 = i_c*f_nl(ir,j,in)
           fx(iy,ix,i_omp) = p*f0
           fy(iy,ix,i_omp) = iy*f0
        enddo
     enddo

     call fftw_execute_dft_c2r(plan_c2r,fx(:,:,i_omp),uxmany(:,:,j))
     call fftw_execute_dft_c2r(plan_c2r,fy(:,:,i_omp),uymany(:,:,j))
  enddo ! j
!$omp end do
!$omp end parallel

  call timer_lib_out('nl')

  if (.not. is_staggered_comm_2) then ! stagger comm2, to load ballance network traffic
    call timer_lib_in('nl_comm')
    call cgyro_nl_fftw_comm2
    call timer_lib_out('nl_comm')
  endif

  call timer_lib_in('nl')

!$omp parallel private(in,iy,ir,p,ix,g0,i_omp,j)
!$omp do schedule(dynamic,1)
  do j=1,nsplit
     i_omp = omp_get_thread_num()+1

     gx(:,:,i_omp) = 0.0
     gy(:,:,i_omp) = 0.0

     ! Array mapping
     do ir=1,n_radial
        p  = ir-1-nx0/2
        ix = p
        if (ix < 0) ix = ix+nx  
        do in=1,n_toroidal
           iy = in-1
           g0 = i_c*g_nl(ir,j,in)
           gx(iy,ix,i_omp) = p*g0
           gy(iy,ix,i_omp) = iy*g0
        enddo
     enddo

     call fftw_execute_dft_c2r(plan_c2r,gx(:,:,i_omp),vx(:,:,i_omp))
     call fftw_execute_dft_c2r(plan_c2r,gy(:,:,i_omp),vy(:,:,i_omp))

     ! Poisson bracket in real space

     uv(:,:,i_omp) = (uxmany(:,:,j)*vy(:,:,i_omp)-uymany(:,:,j)*vx(:,:,i_omp))/(nx*ny)

     call fftw_execute_dft_r2c(plan_r2c,uv(:,:,i_omp),fx(:,:,i_omp))

     ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
     ! that will be filtered in the main time-stepping loop

     do ir=1,n_radial 
        ix = ir-1-nx0/2
        if (ix < 0) ix = ix+nx
        do in=1,n_toroidal
           iy = in-1
           g_nl(ir,j,in) = fx(iy,ix,i_omp)
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
