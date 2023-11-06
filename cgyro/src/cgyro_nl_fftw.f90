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
  integer :: ir,itm,itl,itor

  include 'fftw3.f03'

  ! Poisson bracket in real space

  uv(:,:,i_omp) = (ux(:,:,i_omp)*vymany(:,:,j)-uy(:,:,i_omp)*vxmany(:,:,j))/(nx*ny)

  call fftw_execute_dft_r2c(plan_r2c,uv(:,:,i_omp),fx(:,:,i_omp))

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

  ! this should really be accounted against nl_mem, but hard to do with OMP
  do itm=1,n_toroidal_procs
   do itl=1,nt_loc
    itor=itl + (itm-1)*nt_loc
    do ir=1,n_radial
     ix = ir-1-nx0/2
     if (ix < 0) ix = ix+nx
     iy = itor-1
     f_nl(ir,itl,j,itm) = fx(iy,ix,i_omp)
    enddo
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
  integer :: ir,it,itm,itl,it_loc
  integer :: itor,mytm
  integer :: j,p
  integer :: i_omp
  integer :: jtheta_min

  complex :: f0,g0

  integer, external :: omp_get_thread_num

  include 'fftw3.f03'
  
  ! time to wait for the g_nl to become avaialble
  call timer_lib_in('nl_comm')
  call parallel_slib_f_fd_wait(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)
  ! make sure f_req progresses
  call parallel_slib_test(f_req)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

! g_nl      is (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! jcev_c_nl is (n_field,n_radial,n_jtheta,nv_loc,nt_loc,n_toroidal_procs)
!$omp parallel do schedule(dynamic,1) &
!$omp& private(itor,mytm,itm,itl,iy,ir,p,ix,g0,i_omp,j,it,iv_loc,it_loc,jtheta_min)
  do j=1,nsplit
        i_omp = omp_get_thread_num()+1

        ! zero elements not otherwise set below
        gx(0:ny2,nx2:nx0-1,i_omp) = 0.0
        gy(0:ny2,nx2:nx0-1,i_omp) = 0.0

        ! Array mapping
        do ir=1,n_radial
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx
           do itm=1,n_toroidal_procs
            do itl=1,nt_loc
              itor = itl + (itm-1)*nt_loc
              mytm = 1 + nt1/nt_loc !my toroidal proc number
              it = 1+((mytm-1)*nsplit+j-1)/nv_loc
              iv_loc = 1+modulo((mytm-1)*nsplit+j-1,nv_loc)
              jtheta_min = 1+((mytm-1)*nsplit)/nv_loc
              it_loc = it-jtheta_min+1

              iy = itor-1
              if (it > n_theta) then
                 g0 = (0.0,0.0)
              else
                 g0 = i_c*sum( jvec_c_nl(:,ir,it_loc,iv_loc,itor)*g_nl(:,ir,it_loc,itor))
              endif
              gx(iy,ix,i_omp) = p*g0
              gy(iy,ix,i_omp) = iy*g0
            enddo
           enddo
           if ((ix/=0) .and. (ix<(nx/2))) then ! happens after ix>nx/2
              ! Average elements so as to ensure
              !   g(kx,ky=0) = g(-kx,ky=0)^*
              ! This symmetry is required for complex input to c2r
              g0 = 0.5*( gx(0,ix,i_omp)+conjg(gx(0,nx-ix,i_omp)) )
              gx(0,ix   ,i_omp) = g0
              gx(0,nx-ix,i_omp) = conjg(g0)
           endif
           gx(n_toroidal:ny2,ix,i_omp) = 0.0
           gy(n_toroidal:ny2,ix,i_omp) = 0.0
        enddo

        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call parallel_slib_test(f_req)
        endif

        call fftw_execute_dft_c2r(plan_c2r,gx(:,:,i_omp),vxmany(:,:,j))

        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call parallel_slib_test(f_req)
        endif

        call fftw_execute_dft_c2r(plan_c2r,gy(:,:,i_omp),vymany(:,:,j))
  enddo ! j

  call timer_lib_out('nl')

  ! time to wait for the f_nl to become avaialble
  call timer_lib_in('nl_comm')
  call parallel_slib_f_nc_wait(fpack,f_nl,f_req)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc
!$omp parallel do private(itm,itl,itor,iy,ir,p,ix,f0,i_omp,j)
  do j=1,nsplit
        i_omp = omp_get_thread_num()+1

        ! zero elements not otherwise set below
        fx(0:ny2,nx2:nx0-1,i_omp) = 0.0
        fy(0:ny2,nx2:nx0-1,i_omp) = 0.0

        ! Array mapping
        do ir=1,n_radial
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx
           do itm=1,n_toroidal_procs
            do itl=1,nt_loc
              itor=itl + (itm-1)*nt_loc
              iy = itor-1
              f0 = i_c*f_nl(ir,itl,j,itm)
              fx(iy,ix,i_omp) = p*f0
              fy(iy,ix,i_omp) = iy*f0
            enddo
           enddo
           if ((ix/=0) .and. (ix<(nx/2))) then ! happens after ix>nx/2
             ! Average elements so as to ensure
             !   f(kx,ky=0) = f(-kx,ky=0)^*
             ! This symmetry is required for complex input to c2r
             f0 = 0.5*( fx(0,ix,i_omp)+conjg(fx(0,nx-ix,i_omp)) )
             fx(0,ix   ,i_omp) = f0
             fx(0,nx-ix,i_omp) = conjg(f0)
           endif
           fx(n_toroidal:ny2,ix,i_omp) = 0.0
           fy(n_toroidal:ny2,ix,i_omp) = 0.0
        enddo

        call fftw_execute_dft_c2r(plan_c2r,fx(:,:,i_omp),ux(:,:,i_omp))
        call fftw_execute_dft_c2r(plan_c2r,fy(:,:,i_omp),uy(:,:,i_omp))

        call cgyro_nl_fftw_stepr(j, i_omp)
  enddo ! j

  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc(f_nl,fpack)
  call timer_lib_out('nl_comm')

  call cgyro_nl_fftw_comm1_r(ij)

end subroutine cgyro_nl_fftw
