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

subroutine cgyro_nl_fftw_stepr(j, i_omp)

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: j, i_omp
  integer :: ix,iy
  integer :: ir,in

  include 'fftw3.f03'

  ! Poisson bracket in real space

!$omp parallel workshare
  uv(:,:,i_omp) = (uxmany(:,:,j)*vy(:,:,i_omp)-uymany(:,:,j)*vx(:,:,i_omp))/(nx*ny)
!$omp end parallel workshare

  call fftw_execute_dft_r2c(plan_r2c,uv(:,:,i_omp),fx(:,:,i_omp))

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

  ! this should really be accounted against nl_mem, but hard to do with OMP
  do ir=1,n_radial
     ix = ir-1-nx0/2
     if (ix < 0) ix = ix+nx
     do in=1,n_toroidal
        iy = in-1
        f_nl(ir,j,in) = fx(iy,ix,i_omp)
     enddo
  enddo

end subroutine cgyro_nl_fftw_stepr

! NOTE: call cgyro_nl_fftw_comm1 before cgyro_nl_fftw
subroutine cgyro_nl_fftw(ij)

  use timer_lib
  use parallel_lib
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: ij
  !-----------------------------------
  integer :: ix,iy
  integer :: ir,it,in
  integer :: it_loc
  integer :: j,p,iexch
  integer :: i_omp
  logical :: force_early_comm2, one_pass_fft
  integer :: o,num_one_pass

  integer :: nsplit_chunks, c,c1,c2

  complex :: f0,g0

  integer, external :: omp_get_thread_num

  include 'fftw3.f03'
 
  nsplit_chunks = nsplit/n_omp
  if (modulo(nsplit,n_omp)/=0) nsplit_chunks = nsplit_chunks + 1

  ! no staggering, will need all threads for FFTW in one pass
  force_early_comm2 = (n_omp>=(4*nsplit)) 
  one_pass_fft = force_early_comm2

  ! time to wait for the F_nl to become avaialble
  call timer_lib_in('nl_comm')
  call parallel_slib_f_nc_wait(fpack,f_nl,f_req)
  ! make sure g_req progresses
  call parallel_slib_test(g_req)
  call timer_lib_out('nl_comm')

  if (force_early_comm2) then
     ! time to wait for the g_nl to become avaialble
     call timer_lib_in('nl_comm')
     call parallel_slib_f_fd_wait(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)
     call timer_lib_out('nl_comm')
  endif

  call timer_lib_in('nl')

  ! run in chunks to minimize temp storage costs
  do c=1,nsplit_chunks
    c1 = ((c-1)*n_omp)+1
    c2 = c*n_omp
    if (c2>nsplit) c2=nsplit

!$omp parallel do schedule(static,1) private(in,iy,ir,p,ix,f0,i_omp,j)
    do j=c1,c2
        i_omp = j-c1+1

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

        call cleanx(fx(0,:,i_omp),nx)
        call cleanx(fy(0,:,i_omp),nx)
    enddo
    
    ! run these sequentially, the ffts themselves are parallel
    do j=c1,c2
        i_omp = j-c1+1
        call fftw_execute_dft_c2r(plan_c2r,fx(:,:,i_omp),uxmany(:,:,j))
        call fftw_execute_dft_c2r(plan_c2r,fy(:,:,i_omp),uymany(:,:,j))
    enddo

  enddo

  call timer_lib_out('nl')
  
  ! stagger comm2, to load ballance network traffic
  if (.not. force_early_comm2) then
     ! time to wait for the g_nl to become avaialble
     call timer_lib_in('nl_comm')
     call parallel_slib_f_fd_wait(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)
     call timer_lib_out('nl_comm')
  endif

  call timer_lib_in('nl')

  ! run in chunks to minimize temp storage costs
  do c=1,nsplit_chunks
    c1 = ((c-1)*n_omp)+1
    c2 = c*n_omp
    if (c2>nsplit) c2=nsplit

!$omp parallel do schedule(static,1) private(in,iy,ir,p,ix,g0,i_omp,j,it,iv_loc,it_loc)
    do j=c1,c2
        i_omp = j-c1+1

        gx(:,:,i_omp) = 0.0
        gy(:,:,i_omp) = 0.0

        ! Array mapping
        do ir=1,n_radial
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx  
           do in=1,n_toroidal
              it = it_j(j,in)
              iv_loc =iv_j(j,in)

              iy = in-1
              if (iv_loc == 0) then
                 g0 = (0.0,0.0)
              else
                 it_loc = it-jtheta_min+1
                 g0 = i_c*sum( jvec_c_nl(:,ir,it_loc,iv_loc,in)*g_nl(:,ir,it_loc,in))
              endif
              gx(iy,ix,i_omp) = p*g0
              gy(iy,ix,i_omp) = iy*g0
           enddo
        enddo

        call cleanx(gx(0,:,i_omp),nx)
        call cleanx(gy(0,:,i_omp),nx)
    enddo

    ! run these sequentially, the ffts themselves are parallel
    do j=c1,c2
        i_omp = j-c1+1
        call fftw_execute_dft_c2r(plan_c2r,gx(:,:,i_omp),vx(:,:,i_omp))
        call fftw_execute_dft_c2r(plan_c2r,gy(:,:,i_omp),vy(:,:,i_omp))
        call cgyro_nl_fftw_stepr(j, i_omp)
    enddo ! j

  enddo

  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc(f_nl,fpack)
  call timer_lib_out('nl_comm')

  call cgyro_nl_fftw_comm1_r(ij)

end subroutine cgyro_nl_fftw

subroutine cleanx(f,n)

  implicit none

  integer, intent(in) :: n
  complex, intent(inout) :: f(0:n-1)

  integer :: i

  ! Average elements so as to ensure
  !
  !   f(kx,ky=0) = f(-kx,jy=0)^*
  !
  ! This symmetry is required for complex input to the FFTW
  ! c2r transform 
  
  do i=1,n/2-1
     f(i)   = 0.5*( f(i)+conjg(f(n-i)) )
     f(n-i) = conjg(f(i)) 
  enddo

end subroutine cleanx
