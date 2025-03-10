!-----------------------------------------------------------------
! cgyro_upwind.f90
!
! PURPOSE:
!  Compute corrected distribution used in conservative advection scheme:
!
!                      /
!              J0 |vp| | dv J0 |vp| g 
!                      /
!  |vp| g -  --------------------------
!                   /
!                   | dv J0^2 |vp|
!                   /
!
! NOTE: |vp| prefactor in conserving term added on Jan 4 2023
!-----------------------------------------------------------------

subroutine cgyro_upwind_prepare_async_r64

  use cgyro_globals
  use parallel_lib
  use timer_lib

  implicit none

  integer :: is,itor
  complex :: res_loc
  complex :: g_val

  call timer_lib_in('str')

#if defined(OMPGPU) || defined(_OPENACC)

  ! running on the GPU
  ! The vector/SIMT nature of the GPU allows us to compute effeciently both together
#if defined(OMPGPU)
  ! No async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&         private(res_loc,iv,iv_loc,g_val) 
#else
!$acc parallel loop collapse(3) gang vector independent async(1) &
!$acc&         private(res_loc,iv,iv_loc,g_val) &
!$acc&         present(g_x,h_x,z,temp,jvec_c,field,upfac1,is_v,upwind_res_loc) &
!$acc&         present(nt1,nt2,ns1,ns2,nc,nv1,nv2,n_field) default(none)
#endif
  do itor=nt1,nt2
   do is=ns1,ns2
     do ic=1,nc
       res_loc = (0.0,0.0)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             g_val = h_x(ic,iv_loc,itor)
             if (n_field > 1) then
                g_val = g_val + & 
                   (z(is)/temp(is))*jvec_c(2,ic,iv_loc,itor)*field(2,ic,itor)
             endif
             res_loc = res_loc+upfac1(ic,iv_loc,itor)*g_val
             g_x(ic,iv_loc,itor) = g_val
          endif
       enddo
       upwind_res_loc(ic,is,itor) = res_loc
    enddo
   enddo
  enddo

#else

  ! running on the CPU
  ! Split, to avoid false sharing
!$omp parallel do collapse(3) &
!$omp&         private(g_val,iv,iv_loc,is) 
  do itor=nt1,nt2
    do iv=nv1,nv2
      do ic=1,nc
          is = is_v(iv)
          iv_loc = iv-nv1+1
          g_val = h_x(ic,iv_loc,itor)
          if (n_field > 1) then
             g_val = g_val + & 
                   (z(is)/temp(is))*jvec_c(2,ic,iv_loc,itor)*field(2,ic,itor)
          endif
          g_x(ic,iv_loc,itor) = g_val
      enddo
    enddo
  enddo
!$omp parallel do collapse(3) &
!$omp&         private(res_loc,iv,iv_loc,is) 
  do itor=nt1,nt2
   do is=ns1,ns2
     do ic=1,nc
       res_loc = (0.0,0.0)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             res_loc = res_loc+upfac1(ic,iv_loc,itor)*g_x(ic,iv_loc,itor)
          endif
       enddo
       upwind_res_loc(ic,is,itor) = res_loc
    enddo
   enddo
  enddo

#endif

  call timer_lib_out('str')

end subroutine cgyro_upwind_prepare_async_r64

subroutine cgyro_upwind_complete_r64

  use cgyro_globals
  use parallel_lib
  use timer_lib

  implicit none

  integer :: is,ie,ix,itor

  call timer_lib_in('str_comm')

  call parallel_clib_sum_upwind(upwind_res_loc,upwind_res)

  call timer_lib_out('str_comm')

  call timer_lib_in('str')

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&         private(iv_loc,is,ix,ie,ic)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) independent gang vector &
!$acc&         present(is_v,ix_v,ie_v,xi,vel,upfac2,g_x,upwind_res) &
!$acc&         private(iv_loc,is,ix,ie,ic) present(nt1,nt2,nv1,nv2,nc) default(none)
#else
!$omp parallel do collapse(2) private(iv_loc,is,ix,ie,ic)
#endif
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        g_x(ic,iv_loc,itor) = abs(xi(ix))*vel(ie)*g_x(ic,iv_loc,itor) &
             -upfac2(ic,iv_loc,itor)*upwind_res(ic,is,itor)
     enddo
   enddo
  enddo

  call timer_lib_out('str')

end subroutine cgyro_upwind_complete_r64

subroutine cgyro_upwind_prepare_async_r32

  use, intrinsic :: iso_fortran_env
  use cgyro_globals
  use parallel_lib
  use timer_lib

  implicit none

  integer :: is,itor
  complex(KIND=REAL32) :: res_loc
  complex :: g_val

  call timer_lib_in('str')

#if defined(OMPGPU) || defined(_OPENACC)

  ! running on the GPU
  ! The vector/SIMT nature of the GPU allows us to compute effeciently both together
#if defined(OMPGPU)
  ! no sync for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&         private(res_loc,iv,iv_loc,g_val) 
#else
!$acc parallel loop collapse(3) gang vector independent async(1) &
!$acc&         private(res_loc,iv,iv_loc,g_val) &
!$acc&         present(g_x,h_x,z,temp,jvec_c,field,upfac1,is_v,upwind32_res_loc) &
!$acc&         present(nt1,nt2,ns1,ns2,nc,nv1,nv2,n_field) default(none)
#endif
  do itor=nt1,nt2
   do is=ns1,ns2
     do ic=1,nc
       res_loc = (0.0,0.0)

       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             g_val = h_x(ic,iv_loc,itor)
             if (n_field > 1) then
                g_val = g_val + & 
                   (z(is)/temp(is))*jvec_c(2,ic,iv_loc,itor)*field(2,ic,itor)
             endif
             res_loc = res_loc+upfac1(ic,iv_loc,itor)*g_val
             g_x(ic,iv_loc,itor) = g_val
          endif
       enddo
       upwind32_res_loc(ic,is,itor) = res_loc
    enddo
   enddo
  enddo

#else

  ! running on the CPU
  ! Split, to avoid false sharing
!$omp parallel do collapse(3) &
!$omp&         private(g_val,iv,iv_loc,is) 
  do itor=nt1,nt2
    do iv=nv1,nv2
      do ic=1,nc
          is = is_v(iv)
          iv_loc = iv-nv1+1
          g_val = h_x(ic,iv_loc,itor)
          if (n_field > 1) then
             g_val = g_val + & 
                   (z(is)/temp(is))*jvec_c(2,ic,iv_loc,itor)*field(2,ic,itor)
          endif
          g_x(ic,iv_loc,itor) = g_val
      enddo
    enddo
  enddo
!$omp parallel do collapse(3) &
!$omp&         private(res_loc,iv,iv_loc,is) 
  do itor=nt1,nt2
   do is=ns1,ns2
     do ic=1,nc
       res_loc = (0.0,0.0)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             res_loc = res_loc+upfac1(ic,iv_loc,itor)*g_x(ic,iv_loc,itor)
          endif
       enddo
       upwind32_res_loc(ic,is,itor) = res_loc
    enddo
   enddo
  enddo

#endif

  call timer_lib_out('str')

end subroutine cgyro_upwind_prepare_async_r32

subroutine cgyro_upwind_complete_r32

  use cgyro_globals
  use parallel_lib
  use timer_lib

  implicit none

  integer :: is,ie,ix,itor

  call timer_lib_in('str_comm')

  call parallel_clib_sum_upwind32(upwind32_res_loc,upwind32_res)

  call timer_lib_out('str_comm')

  call timer_lib_in('str')

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&         private(iv_loc,is,ix,ie,ic)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) independent gang vector &
!$acc&         present(is_v,ix_v,ie_v,xi,vel,upfac2,g_x,upwind32_res) &
!$acc&         private(iv_loc,is,ix,ie,ic) &
!$acc&         present(nt1,nt2,nv1,nv2,nc) default(none)
#else
!$omp parallel do collapse(2) private(iv_loc,is,ix,ie,ic)
#endif
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        g_x(ic,iv_loc,itor) = abs(xi(ix))*vel(ie)*g_x(ic,iv_loc,itor) &
             -upfac2(ic,iv_loc,itor)*upwind32_res(ic,is,itor)
     enddo
   enddo
  enddo

  call timer_lib_out('str')

end subroutine cgyro_upwind_complete_r32

subroutine cgyro_upwind_prepare_async

  use cgyro_globals

  implicit none

  if (upwind_single_flag == 0) then
    call cgyro_upwind_prepare_async_r64
  else
    call cgyro_upwind_prepare_async_r32
  endif

end subroutine cgyro_upwind_prepare_async

subroutine cgyro_upwind_complete

  use cgyro_globals

  implicit none

  if (upwind_single_flag == 0) then
    call cgyro_upwind_complete_r64
  else
    call cgyro_upwind_complete_r32
  endif

end subroutine cgyro_upwind_complete

