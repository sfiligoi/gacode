!-----------------------------------------------------------------
! cgyro_nl_fftw.gpu.f90 [GPU (acc-cuFFT/hipFFT) version]
!
! PURPOSE:
!  Evaluate nonlinear bracket with dealiased FFT.  It is natural 
!  to use the FFTW complex-to-real (c2r) transform to move to real 
!  space, compute the convolution (uv), and transform back to the 
!  spectral form using a real-to-complex transform (r2c).
!  
! NOTE: Need to be careful with (p=-nr/2,n=0) component.
!-----------------------------------------------------------------

subroutine cgyro_nl_fftw_zero4(sz,v1,v2,v3,v4)
  implicit none

  integer, intent(in) :: sz
  complex, dimension(*), intent(inout) :: v1,v2,v3,v4

  integer :: i

!$acc parallel loop independent gang vector &
!$acc&         present(v1,v2,v3,v4) private(i)
  do i=1,sz
    v1(i) = 0.0
    v2(i) = 0.0
    v3(i) = 0.0
    v4(i) = 0.0
  enddo

end subroutine

subroutine cgyro_nl_fftw_mul(sz,uvm,uxm,vym,uym,vxm,inv_nxny)
  implicit none

  integer, intent(in) :: sz
  real, dimension(*),intent(out) :: uvm
  real, dimension(*),intent(in) :: uxm,vym,uym,vxm
  real, intent(in) :: inv_nxny

  integer :: i

!$acc parallel loop independent gang vector &
!$acc&         present(uvm,uxm,vym,uym,vxm) private(i)
  do i=1,sz
    uvm(i) = (uxm(i)*vym(i)-uym(i)*vxm(i))*inv_nxny
  enddo

end subroutine

subroutine cgyro_nl_fftw(ij)

#ifdef HIPGPU
  use hipfort_hipfft
#else
  use cufft
#endif
  use timer_lib
  use parallel_lib
  use cgyro_nl_comm
  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: ij
  !-----------------------------------
  integer :: j,p,iexch
  integer :: it,ir,itm,itl,ix,iy
  integer :: itor,mytm
  integer :: i1,i2
  integer :: it_loc
  integer :: ierr
  integer :: rc
  complex :: f0,g0
  integer :: jtheta_min

  real :: inv_nxny


  ! time to wait for the F_nl to become avaialble
  call timer_lib_in('nl_comm')
  call parallel_slib_f_nc_wait(fpack,f_nl,f_req)
  ! make sure g_req progresses
  call parallel_slib_test(g_req)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl_mem')
!$acc  data present(f_nl)  &
!$acc&      present(fxmany,fymany,gxmany,gymany) &
!$acc&      present(uxmany,uymany,vxmany,vymany) &
!$acc&      present(uvmany)

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl')

  call cgyro_nl_fftw_zero4(size(fxmany,1)*size(fxmany,2)*size(fxmany,3), &
                           fxmany,fymany,gxmany,gymany)

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc
!$acc parallel loop gang vector independent collapse(4) private(j,ir,p,ix,itor,iy,f0,g0) async
  do j=1,nsplit
     do ir=1,n_radial
       do itm=1,n_toroidal_procs
         do itl=1,nt_loc
           itor=itl + (itm-1)*nt_loc
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx

           iy = itor-1
           f0 = i_c*f_nl(ir,itl,j,itm)
           fxmany(iy,ix,j) = p*f0
           fymany(iy,ix,j) = iy*f0
         enddo
       enddo
     enddo
  enddo

  ! make sure g_req progresses
  call parallel_slib_test(g_req)

!$acc wait
  ! Average elements so as to ensure
  !   f(kx,ky=0) = f(-kx,ky=0)^*
  ! This symmetry is required for complex input to c2r
!$acc parallel loop gang vector independent collapse(2) async(2) &
!$acc&         private(j,ix,f0) &
!$acc&         present(fxmany) &
!$acc&         present(nsplit,nx)
  do j=1,nsplit
    do ix=1,nx/2-1
      f0 = 0.5*( fxmany(0,ix,j)+conjg(fxmany(0,nx-ix,j)) )
      fxmany(0,ix   ,j) = f0
      fxmany(0,nx-ix,j) = conjg(f0)
    enddo
  enddo

  ! make sure g_req progresses
  call parallel_slib_test(g_req)

     ! --------------------------------------
     ! perform many Fourier Transforms at once
     ! --------------------------------------

!$acc  host_data &
!$acc& use_device(fymany) &
!$acc& use_device(uymany)
#ifdef HIPGPU
  rc = hipfftExecZ2D(hip_plan_c2r_many,c_loc(fymany),c_loc(uymany))
#else
  rc = cufftExecZ2D(cu_plan_c2r_many,fymany,uymany)
#endif

!$acc end host_data

  ! make sure g_req progresses
  call parallel_slib_test(g_req)

!$acc wait
  ! fxmany is complete now

!$acc  host_data &
!$acc& use_device(fxmany) &
!$acc& use_device(uxmany)

#ifdef HIPGPU
  rc = hipfftExecZ2D(hip_plan_c2r_many,c_loc(fxmany),c_loc(uxmany))
#else
  rc = cufftExecZ2D(cu_plan_c2r_many,fxmany,uxmany)
#endif

!$acc end host_data
  call timer_lib_out('nl')

  ! time to wait for the g_nl to become avaialble
  call timer_lib_in('nl_comm')
  call parallel_slib_f_fd_wait(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

! g_nl      is (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! jcev_c_nl is (n_field,n_radial,n_jtheta,nv_loc,nt_loc,n_toroidal_procs)

!$acc parallel loop gang vector independent collapse(4) &
!$acc&         private(j,ir,p,ix,itor,mytm,iy,g0,it,iv_loc,it_loc,jtheta_min) &
!$acc&         present(g_nl,jvec_c_nl) &
!$acc&         present(nsplit,n_radial,n_toroidal_procs,nt_loc,nt1,n_theta,nv_loc,nx0)
  do j=1,nsplit
     do ir=1,n_radial
       do itm=1,n_toroidal_procs
         do itl=1,nt_loc
           itor = itl + (itm-1)*nt_loc
           mytm = 1 + nt1/nt_loc !my toroidal proc number
           it = 1+((mytm-1)*nsplit+j-1)/nv_loc
           iv_loc = 1+modulo((mytm-1)*nsplit+j-1,nv_loc)
           jtheta_min = 1+((mytm-1)*nsplit)/nv_loc
           it_loc = it-jtheta_min+1

           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx

           iy = itor-1
           if (it > n_theta) then
              g0 = (0.0,0.0)
           else
              g0 = i_c*sum( jvec_c_nl(1:n_field,ir,it_loc,iv_loc,itor)*g_nl(1:n_field,ir,it_loc,itor))
           endif
           gxmany(iy,ix,j) = p*g0
           gymany(iy,ix,j) = iy*g0
         enddo
       enddo
     enddo
  enddo

  ! Average elements so as to ensure
  !   g(kx,ky=0) = g(-kx,ky=0)^*
  ! This symmetry is required for complex input to c2r
!$acc parallel loop gang vector independent collapse(2) async(2) &
!$acc&         private(j,ix,g0) &
!$acc&         present(gxmany) &
!$acc&         present(nsplit,nx)
  do j=1,nsplit
    do ix=1,nx/2-1
      g0 = 0.5*( gxmany(0,ix,j)+conjg(gxmany(0,nx-ix,j)) )
      gxmany(0,ix   ,j) = g0
      gxmany(0,nx-ix,j) = conjg(g0)
    enddo
  enddo

!$acc  host_data &
!$acc& use_device(gymany) &
!$acc& use_device(vymany)

#ifdef HIPGPU
  rc = hipfftExecZ2D(hip_plan_c2r_many,c_loc(gymany),c_loc(vymany))
#else
  rc = cufftExecZ2D(cu_plan_c2r_many,gymany,vymany)
#endif

!$acc end host_data

!$acc wait
  ! gxmany is complete now

!$acc  host_data &
!$acc& use_device(gxmany) &
!$acc& use_device(vxmany)

#ifdef HIPGPU
  rc = hipfftExecZ2D(hip_plan_c2r_many,c_loc(gxmany),c_loc(vxmany))
#else
  rc = cufftExecZ2D(cu_plan_c2r_many,gxmany,vxmany)
#endif

!$acc wait
!$acc end host_data

  ! Poisson bracket in real space
  ! uv = (ux*vy-uy*vx)/(nx*ny)

  inv_nxny = 1.0/(nx*ny)

  call cgyro_nl_fftw_mul(size(uvmany,1)*size(uvmany,2)*size(uvmany,3), &
                         uvmany,uxmany,vymany,uymany,vxmany,inv_nxny)

  ! ------------------
  ! Transform uv to fx
  ! ------------------

!$acc wait
!$acc host_data use_device(uvmany,fxmany)
#ifdef HIPGPU
  rc = hipfftExecD2Z(hip_plan_r2c_many,c_loc(uvmany),c_loc(fxmany))
#else
  rc = cufftExecD2Z(cu_plan_r2c_many,uvmany,fxmany)
#endif
!$acc wait
!$acc end host_data
!$acc wait

  call timer_lib_out('nl')
  call timer_lib_in('nl_mem')

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

!$acc parallel loop independent collapse(4) gang vector &
!$acc&         private(itor,ix,iy) present(f_nl,fxmany)
  do itm=1,n_toroidal_procs
     do itl=1,nt_loc
       do j=1,nsplit
         do ir=1,n_radial 
           itor=itl + (itm-1)*nt_loc
           ix = ir-1-nx0/2
           if (ix < 0) ix = ix+nx

           iy = itor-1
           f_nl(ir,itl,j,itm) = fxmany(iy,ix,j)
         enddo
       enddo
    enddo
  enddo

  ! end data f_nl
!$acc end data

  call timer_lib_out('nl_mem')

  call cgyro_nl_fftw_comm1_r(ij)

end subroutine cgyro_nl_fftw
