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

  complex :: f0,g0
  complex, dimension(:,:), allocatable :: fpack
  complex, dimension(:,:), allocatable :: gpack
  real :: inv_nxny
  real :: r_ux,r_uy,r_vx,r_vy,r_uv



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

#ifdef _OPENACC
!$acc  data pcopyin(f_nl) pcopy(g_nl)   &
!$acc& pcreate(fxmany,fymany,gxmany,gymany) &
!$acc& pcreate(uxmany,uymany,vxmany,vymany) &
!$acc& pcreate(uvmany)

!$acc parallel
!$acc loop gang 
#else
!$omp parallel do private(fx,gx,fy,gy,in,iy,ir,p,ix,f0,g0)
#endif
  do j=1,nsplit

!$acc loop worker
    do ix=lbound(fxmany,2),ubound(fxmany,2)
!$acc loop vector
    do iy=lbound(fxmany,1),ubound(fxmany,1)
     fxmany(iy,ix,j) = 0.0
     gxmany(iy,ix,j) = 0.0
     fymany(iy,ix,j) = 0.0
     gymany(iy,ix,j) = 0.0
    enddo
    enddo

     ! Array mapping
!$acc loop worker private(ir,p,ix)
     do ir=1,n_radial
        p  = ir-1-nx0/2
        ix = p
        if (ix < 0) ix = ix+nx  
!$acc   loop  vector private(in,iy,f0,g0)
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

     if (kxfilter_flag == 1) then
!$acc  loop vector private(iy)
       do iy=lbound(fxmany,1),ubound(fxmany,1)
        fxmany(iy,-nx0/2+nx,j) = 0.0
        fymany(iy,-nx0/2+nx,j) = 0.0
        gxmany(iy,-nx0/2+nx,j) = 0.0
        gymany(iy,-nx0/2+nx,j) = 0.0
       enddo
     endif
   enddo ! end do j
!$acc end parallel

#ifdef _OPENACC
! --------------------------------------
! perform many Fourier Transforms at once
! --------------------------------------
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
!$acc end host_data
#else
!$omp  parallel do   
   do j=1,nsplit
     call fftw_execute_dft_c2r(plan_c2r, fxmany(:,:,j),uxmany(:,:,j))
     call fftw_execute_dft_c2r(plan_c2r, fymany(:,:,j),uymany(:,:,j))
     call fftw_execute_dft_c2r(plan_c2r, gxmany(:,:,j),vxmany(:,:,j))
     call fftw_execute_dft_c2r(plan_c2r, gymany(:,:,j),vymany(:,:,j))
    enddo
#endif



     ! Poisson bracket in real space
     ! uv = (ux*vy-uy*vx)/(nx*ny)

     inv_nxny = dble(1)/dble(nx*ny)

#ifdef _OPENACC
!$acc  parallel 
!$acc  loop gang 
#else
!$omp parallel do default(none) collapse(2) &
!$omp& private(j,ix,iy) &
!$omp& private(r_ux,r_uy,r_vx,r_vy,r_uv) &
!$omp& shared(nsplit,inv_nxny) &
!$omp& shared(uxmany,uymany,vxmany,vymany) &
!$omp& shared(uvmany)
#endif
     do j=1,nsplit
!$acc loop worker
     do ix=lbound(ux,2),ubound(ux,2)
!$acc loop vector
     do iy=lbound(ux,1),ubound(ux,1)


       r_ux = uxmany(iy,ix,j)
       r_uy = uymany(iy,ix,j)
       r_vx = vxmany(iy,ix,j)
       r_vy = vymany(iy,ix,j)

       r_uv = (r_ux*r_vy-r_uy*r_vx)*inv_nxny

       uvmany(iy,ix,j) = r_uv

     enddo
     enddo
     enddo
!$acc end parallel



! ------------------
! Transform uv to fx
! ------------------

       
#ifdef _OPENACC
!$acc host_data use_device(uvmany,fxmany)
   if (kind(uvmany).eq.singlePrecision) then
     call cufftExecR2C(cu_plan_r2c_many,uvmany,fxmany)
   else
     call cufftExecD2Z(cu_plan_r2c_many,uvmany,fxmany)
   endif
!$acc end host_data
#else
!$omp parallel do private(j)
     do j=1,nsplit
       call fftw_execute_dft_r2c(plan_r2c,uvmany(:,:,j),fxmany(:,:,j))
     enddo
#endif





     ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
     ! that will be filtered in the main time-stepping loop

#ifdef _OPENACC
!$acc parallel 
!$acc loop gang
#else
!$omp parallel do  &
!$omp& private(ir,in,ix,iy)
#endif
     do j=1,nsplit
!$acc loop worker
     do in=1,n_toroidal
        ix = ir-1-nx0/2
        if (ix < 0) ix = ix+nx
!$acc   loop vector
        do ir=1,n_radial 
           iy = in-1
           g_nl(ir,j,in) = fxmany(iy,ix,j)
        enddo
     enddo
     enddo
!$acc end parallel


#ifdef _OPENACC
!$acc end data
#endif
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
