!-----------------------------------------------------------------
! cgyro_nl_fftw.gpu.f90 [GPU (acc-cuFFT) version]
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

!$acc parallel loop independent present(v1,v2,v3,v4) private(i)
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

!$acc parallel loop independent present(uvm,uxm,vym,uym,vxm) private(i)
  do i=1,sz
    uvm(i) = (uxm(i)*vym(i)-uym(i)*vxm(i))*inv_nxny
  enddo

end subroutine

subroutine cgyro_nl_fftw(ij)

  use cufft
  use timer_lib
  use parallel_lib
  use cgyro_nl_comm
  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: ij
  !-----------------------------------
  integer :: j,p,iexch
  integer :: it,ir,in,ix,iy
  integer :: i1,i2
  integer :: it_loc
  integer :: ierr
  integer :: rc
  complex :: f0,g0

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

!$acc parallel loop gang vector independent collapse(3) private(j,ir,p,ix,in,iy,f0,g0) async
  do j=1,nsplit
     do ir=1,n_radial
        do in=1,n_toroidal
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx

           iy = in-1
           f0 = i_c*f_nl(ir,j,in)
           fxmany(iy,ix,j) = p*f0
           fymany(iy,ix,j) = iy*f0
        enddo
     enddo
  enddo

  ! make sure g_req progresses
  call parallel_slib_test(g_req)

     ! --------------------------------------
     ! perform many Fourier Transforms at once
     ! --------------------------------------
!$acc wait
!$acc  host_data &
!$acc& use_device(fxmany,fymany) &
!$acc& use_device(uxmany,uymany)

  rc = cufftExecZ2D(cu_plan_c2r_many,fxmany,uxmany)
  ! make sure g_req progresses
  call parallel_slib_test(g_req)
  rc = cufftExecZ2D(cu_plan_c2r_many,fymany,uymany)

!$acc wait
!$acc end host_data
  call timer_lib_out('nl')

  ! time to wait for the g_nl to become avaialble
  call timer_lib_in('nl_comm')
  call parallel_slib_f_fd_wait(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')


!$acc parallel loop gang vector independent collapse(3) private(j,ir,p,ix,in,iy,f0,g0,it,iv_loc,it_loc) &
!$acc&         present(g_nl)
  do j=1,nsplit
     do ir=1,n_radial
        do in=1,n_toroidal
           it = it_j(j,in)
           iv_loc =iv_j(j,in)

           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx

           iy = in-1
           if (iv_loc == 0) then
              g0 = (0.0,0.0)
           else
              it_loc = it-jtheta_min+1
              g0 = i_c*sum( jvec_c_nl(1:n_field,ir,it_loc,iv_loc,in)*g_nl(1:n_field,ir,it_loc,in))
           endif
           gxmany(iy,ix,j) = p*g0
           gymany(iy,ix,j) = iy*g0
        enddo
     enddo
  enddo

!$acc  host_data &
!$acc& use_device(gxmany,gymany) &
!$acc& use_device(vxmany,vymany)

  rc = cufftExecZ2D(cu_plan_c2r_many,gxmany,vxmany)
  rc = cufftExecZ2D(cu_plan_c2r_many,gymany,vymany)

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
  rc = cufftExecD2Z(cu_plan_r2c_many,uvmany,fxmany)
!$acc wait
!$acc end host_data
!$acc wait

  call timer_lib_out('nl')
  call timer_lib_in('nl_mem')

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

!$acc parallel loop independent collapse(3) private(j,ir,in,ix,iy)
  do j=1,nsplit
     do ir=1,n_radial 
        do in=1,n_toroidal
           ix = ir-1-nx0/2
           if (ix < 0) ix = ix+nx

           iy = in-1
           f_nl(ir,j,in) = fxmany(iy,ix,j)
        enddo
     enddo
  enddo

  ! end data f_nl
!$acc end data

  call timer_lib_out('nl_mem')

  call cgyro_nl_fftw_comm1_r(ij)

end subroutine cgyro_nl_fftw
