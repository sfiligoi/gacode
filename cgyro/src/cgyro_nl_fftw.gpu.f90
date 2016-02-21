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
#ifdef _OPENACC
  use precision_m
  use cufft_m
#endif

  use timer_lib
  use parallel_lib

  use cgyro_globals
  implicit none

  integer, intent(in) :: ij
  integer :: j,p,iexch
  integer :: it,ir,in,ix,iy
  integer :: i1,i2

  complex*16 :: f0,g0
  complex*16, dimension(:,:), allocatable :: fpack
  complex*16, dimension(:,:), allocatable :: gpack
  real*8 :: inv_nxny

  logical, parameter :: use_cufft = .true.
  logical, parameter :: use_acc = .true.


  include 'fftw3.f03'

  call timer_lib_in('nl_comm')

  allocate(fpack(n_radial,nv_loc*n_theta))
  allocate(gpack(n_radial,nv_loc*n_theta))
!$omp workshare
   fpack = 0.0
   gpack = 0.0
!$omp end workshare
  iexch = 0
!$omp  parallel do  &
!$omp& private(iv_loc,it,ir,iexch,ic_loc)
  do iv_loc=1,nv_loc
     do it=1,n_theta
        ! iexch = iexch+1
        iexch = it + (iv_loc-1)*n_theta
        do ir=1,n_radial
 
           ic_loc = ic_c(ir,it)
           ! ic_loc = it + (ir-1)*n_theta
           fpack(ir,iexch) = h_x(ic_loc,iv_loc)
           gpack(ir,iexch) = psi(ic_loc,iv_loc)
        enddo
     enddo
  enddo
  call parallel_slib_f(fpack,f_nl)
  call parallel_slib_f(gpack,g_nl)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

#ifdef _OPENACC
!$acc  data pcopyin(f_nl) pcopy(g_nl)   &
!$acc& pcreate(fxmany,fymany,gxmany,gymany) &
!$acc& pcreate(uxmany,uymany,vxmany,vymany) &
!$acc& pcreate(uvmany)
#endif

!$acc kernels
    fxmany(:,:,:) = 0
    fymany(:,:,:) = 0
    gxmany(:,:,:) = 0
    gymany(:,:,:) = 0
!$acc end kernels



  if (use_acc) then
!$acc kernels 
  do j=1,nsplit

     ! Array mapping
     do ir=1,n_radial

        p  = ir-1-nx0/2
        ix = p
        if (ix < 0) ix = ix+nx  
        do in=1,n_toroidal
           iy = in-1
           f0 = i_c*f_nl(ir,j,in)
           g0 = i_c*g_nl(ir,j,in)
           fxmany(iy,ix,j) = p*f0
           gxmany(iy,ix,j) = p*g0
           fymany(iy,ix,j) = iy*f0
           gymany(iy,ix,j) = iy*g0
      enddo
      enddo
    enddo
!$acc end kernels
  else
!$acc update host(f_nl,g_nl) &
!$acc& host(fxmany,gxmany,fymany,gymany)
!$acc wait
  do j=1,nsplit

     ! Array mapping
     do ir=1,n_radial

        p  = ir-1-nx0/2
        ix = p
        if (ix < 0) ix = ix+nx  
        do in=1,n_toroidal
           iy = in-1
           f0 = i_c*f_nl(ir,j,in)
           g0 = i_c*g_nl(ir,j,in)
           fxmany(iy,ix,j) = p*f0
           gxmany(iy,ix,j) = p*g0
           fymany(iy,ix,j) = iy*f0
           gymany(iy,ix,j) = iy*g0
      enddo
      enddo
    enddo
!$acc update device(fxmany,gxmany,fymany,gymany)
!$acc wait
  endif





     if (kxfilter_flag == 1) then
!$acc  kernels
!$acc  loop independent gang
       do j=1,nsplit
        fxmany(:,-nx0/2+nx,j) = 0.0
        fymany(:,-nx0/2+nx,j) = 0.0
        gxmany(:,-nx0/2+nx,j) = 0.0
        gymany(:,-nx0/2+nx,j) = 0.0
       enddo
!$acc  end kernels
     endif

   if (use_cufft) then
! --------------------------------------
! perform many Fourier Transforms at once
! --------------------------------------
!$acc wait
!$acc  host_data &
!$acc& use_device(fxmany,fymany,gxmany,gymany) &
!$acc& use_device(uxmany,uymany,vxmany,vymany)

  if (kind(uxmany).eq.singlePrecision) then
    call cufftExecC2R(cu_plan_c2r_many,fxmany,uxmany)
    call cufftExecC2R(cu_plan_c2r_many,fymany,uymany)
    call cufftExecC2R(cu_plan_c2r_many,gxmany,vxmany)
    call cufftExecC2R(cu_plan_c2r_many,gymany,vymany)
  else
    call cufftExecZ2D(cu_plan_c2r_many,fxmany,uxmany)
    call cufftExecZ2D(cu_plan_c2r_many,fymany,uymany)
    call cufftExecZ2D(cu_plan_c2r_many,gxmany,vxmany)
    call cufftExecZ2D(cu_plan_c2r_many,gymany,vymany)
  endif
!$acc wait
!$acc end host_data
!$acc wait
   else
!$acc wait
!$acc update host(fxmany,fymany,gxmany,gymany)
!$acc wait

!$omp  parallel do &
!$omp& private(fx,fy,gx,gy,ux,uy,vx,vy)
   do j=1,nsplit
     fx(:,:) = fxmany(:,:,j)
     fy(:,:) = fymany(:,:,j)
     gx(:,:) = gxmany(:,:,j)
     gy(:,:) = gymany(:,:,j)

     call fftw_execute_dft_c2r(plan_c2r, fx,ux)
     call fftw_execute_dft_c2r(plan_c2r, fy,uy)
     call fftw_execute_dft_c2r(plan_c2r, gx,vx)
     call fftw_execute_dft_c2r(plan_c2r, gy,vy)

     uxmany(:,:,j) = ux(:,:)
     uymany(:,:,j) = uy(:,:)
     vxmany(:,:,j) = vx(:,:)
     vymany(:,:,j) = vy(:,:)
    enddo
!$acc update device(uxmany,uymany,vxmany,vymany)
!$acc wait
   endif



     ! Poisson bracket in real space
     ! uv = (ux*vy-uy*vx)/(nx*ny)

     inv_nxny = dble(1)/dble(nx*ny)

   if (use_acc) then
!$acc  kernels 
!$acc loop independent gang
   do j=1,nsplit
   do ix=lbound(uvmany,2),ubound(uvmany,2)
   do iy=lbound(uvmany,1),ubound(uvmany,1)
    uvmany(iy,ix,j) = (uxmany(iy,ix,j)*vymany(iy,ix,j)- &
                       uymany(iy,ix,j)*vxmany(iy,ix,j))*inv_nxny
   enddo
   enddo
   enddo
!$acc  end kernels
   else
!$acc update host(uxmany,uymany,vxmany,vymany)
!$acc wait

!$omp workshare
  uvmany = (uxmany*vymany-uymany*vxmany)/(nx*ny)
!$omp end workshare

!$acc update device(uvmany)
!$acc wait
   endif





! ------------------
! Transform uv to fx
! ------------------

       
   if (use_cufft) then
!$acc wait
!$acc host_data use_device(uvmany,fxmany)
   if (kind(uvmany).eq.singlePrecision) then
     call cufftExecR2C(cu_plan_r2c_many,uvmany,fxmany)
   else
     call cufftExecD2Z(cu_plan_r2c_many,uvmany,fxmany)
   endif
!$acc wait
!$acc end host_data
!$acc wait
   else
!$acc wait
!$acc update host(uvmany)
!$acc wait
!$omp parallel do private(j,uv,fx)
     do j=1,nsplit
       uv(:,:) = uvmany(:,:,j)
       call fftw_execute_dft_r2c(plan_r2c,uv,fx)
       fxmany(:,:,j) = fx(:,:)
     enddo
!$acc update device(fxmany)
!$acc wait
   endif





     ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
     ! that will be filtered in the main time-stepping loop

     
   if (use_acc) then
!$acc kernels  
     do j=1,nsplit
     do ir=1,n_radial 
        ix = ir-1-nx0/2

        if (ix < 0) ix = ix+nx

        do in=1,n_toroidal
           iy = in-1
           g_nl(ir,j,in) = fxmany(iy,ix,j)
        enddo
     enddo
     enddo
!$acc end kernels
!$acc wait
  else
!$acc update host(fxmany,g_nl)
!$acc wait
     do j=1,nsplit
     do ir=1,n_radial 
        ix = ir-1-nx0/2

        if (ix < 0) ix = ix+nx

        do in=1,n_toroidal
           iy = in-1
           g_nl(ir,j,in) = fxmany(iy,ix,j)
        enddo
     enddo
     enddo
!$acc update device(g_nl)
!$acc wait
  endif



#ifdef _OPENACC
!$acc wait
!$acc end data
!$acc wait
#endif
  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  call parallel_slib_r(g_nl,gpack)
  iexch = 0
!$omp  parallel do &
!$omp& private(iv_loc,it,ir,iexch,ic_loc)
  do iv_loc=1,nv_loc
     do it=1,n_theta
        ! iexch = iexch+1
        iexch = it + (iv_loc-1)*n_theta
        do ir=1,n_radial
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

!$omp workshare
  rhs(:,:,ij) = rhs(:,:,ij)+((q*rho/rmin)*(2*pi/length))*psi(:,:)
!$omp end workshare

end subroutine cgyro_nl_fftw
