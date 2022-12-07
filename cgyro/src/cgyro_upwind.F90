!-----------------------------------------------------------------
! cgyro_upwind.f90
!
! PURPOSE:
!  Compute corrected distribution used in conservative advection scheme:
!
!                  /
!              J0  | dv J0 |vp| g 
!                  /
!  |vp| g -  ----------------------
!                   /
!                   | dv J0^2 
!                   /
!
!-----------------------------------------------------------------

subroutine cgyro_upwind_r64

  use cgyro_globals
  use mpi
  use timer_lib

  implicit none

  integer :: is,ie,ix
#ifdef _OPENACC
  complex :: res_loc_one, res_loc_two
#endif

  call timer_lib_in('str')

#ifdef _OPENACC
!$acc parallel loop collapse(2) gang &
!$acc&         private(res_loc_one,iv) &
!$acc&         present(g_x,upfac1,is_v,upwind_res_loc) default(none)
  do is=ns1,ns2
     do ic=1,nc
       res_loc_one = (0.0,0.0)
       res_loc_two = (0.0,0.0)

!$acc loop vector private(iv_loc) reduction(+:res_loc_one)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             res_loc_one = res_loc_one+upfac1(ic,iv_loc,my_toroidal,1)*g_x(ic,iv_loc)
          endif
       enddo
       upwind_res_loc(ic,is,1) = res_loc_one

!$acc loop vector private(iv_loc) reduction(+:res_loc_two)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             res_loc_two = res_loc_two+upfac1(ic,iv_loc,my_toroidal,2)*g_x(ic,iv_loc)
          endif
       enddo
       upwind_res_loc(ic,is,2) = res_loc_two
    enddo

  enddo
#else
  upwind_res_loc(:,:,:) = (0.0,0.0)

!$omp parallel private(iv_loc,ic,is)
!$omp do reduction(+:upwind_res_loc)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        upwind_res_loc(ic,is,:) = upwind_res_loc(ic,is,:)+upfac1(ic,iv_loc,my_toroidal,:)*g_x(ic,iv_loc)
     enddo
  enddo
!$omp end do
!$omp end parallel
#endif

  call timer_lib_out('str')

  call timer_lib_in('str_comm')

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update host(upwind_res_loc)
#else
!$acc host_data use_device(upwind_res_loc,upwind_res)
#endif

  call MPI_ALLREDUCE(upwind_res_loc(:,:,:),&
       upwind_res(:,:,:),&
       size(upwind_res(:,:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_3,&
       i_err)

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update device(upwind_res)
#else
!$acc end host_data
#endif

  call timer_lib_out('str_comm')

  call timer_lib_in('str')

#ifdef _OPENACC
!$acc parallel loop collapse(2) independent &
!$acc&         present(is_v,ix_v,ie_v,xi,vel,upfac2,g_x,upwind_res,up_cutoff) &
!$acc&         private(iv_loc,is,ix,ie) default(none)
#else
!$omp parallel do private(iv_loc,is,ix,ie,ic)
#endif
  do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        g_x(ic,iv_loc) = abs(xi(ix))*vel(ie)*g_x(ic,iv_loc) &
             -upfac2(ic,iv_loc,my_toroidal,1)*upwind_res(ic,is,1) &
             -upfac2(ic,iv_loc,my_toroidal,2)*upwind_res(ic,is,2)*up_cutoff
     enddo
  enddo

  call timer_lib_out('str')

end subroutine cgyro_upwind_r64

subroutine cgyro_upwind_r32

  use, intrinsic :: iso_fortran_env
  use cgyro_globals
  use mpi
  use timer_lib

  implicit none

  integer :: is,ie,ix
#ifdef _OPENACC
  complex(KIND=REAL32) :: res_loc_one, res_loc_two
#endif

  call timer_lib_in('str')

#ifdef _OPENACC
!$acc parallel loop collapse(2) gang &
!$acc&         private(res_loc_one,iv) &
!$acc&         present(g_x,upfac1,is_v,upwind32_res_loc) default(none)
  do is=ns1,ns2
     do ic=1,nc
       res_loc_one = (0.0,0.0)
       res_loc_two = (0.0,0.0)

!$acc loop vector private(iv_loc) reduction(+:res_loc_one)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             res_loc_one = res_loc_one+upfac1(ic,iv_loc,my_toroidal,1)*g_x(ic,iv_loc)
          endif
       enddo
       upwind32_res_loc(ic,is,1) = res_loc_one

!$acc loop vector private(iv_loc) reduction(+:res_loc_two)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             res_loc_two = res_loc_two+upfac1(ic,iv_loc,my_toroidal,2)*g_x(ic,iv_loc)
          endif
       enddo
       upwind32_res_loc(ic,is,2) = res_loc_two
    enddo

  enddo
#else
  upwind32_res_loc(:,:,:) = (0.0,0.0)

!$omp parallel private(iv_loc,ic,is)
!$omp do reduction(+:upwind32_res_loc)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        upwind32_res_loc(ic,is,:) = upwind32_res_loc(ic,is,:)+upfac1(ic,iv_loc,my_toroidal,:)*g_x(ic,iv_loc)
     enddo
  enddo
!$omp end do
!$omp end parallel
#endif

  call timer_lib_out('str')

  call timer_lib_in('str_comm')

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update host(upwind32_res_loc)
#else
!$acc host_data use_device(upwind32_res_loc,upwind32_res)
#endif

  call MPI_ALLREDUCE(upwind32_res_loc(:,:,:),&
       upwind32_res(:,:,:),&
       size(upwind32_res(:,:,:)),&
       MPI_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_3,&
       i_err)

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update device(upwind32_res)
#else
!$acc end host_data
#endif

  call timer_lib_out('str_comm')

  call timer_lib_in('str')

#ifdef _OPENACC
!$acc parallel loop collapse(2) independent &
!$acc&         present(is_v,ix_v,ie_v,xi,vel,upfac2,g_x,upwind32_res,up_cutoff) &
!$acc&         private(iv_loc,is,ix,ie) default(none)
#else
!$omp parallel do private(iv_loc,is,ix,ie,ic)
#endif
  do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        g_x(ic,iv_loc) = abs(xi(ix))*vel(ie)*g_x(ic,iv_loc) &
             -upfac2(ic,iv_loc,my_toroidal,1)*upwind32_res(ic,is,1) &
             -upfac2(ic,iv_loc,my_toroidal,2)*upwind32_res(ic,is,2)*up_cutoff
     enddo
  enddo

  call timer_lib_out('str')

end subroutine cgyro_upwind_r32

subroutine cgyro_upwind

  use cgyro_globals

  implicit none

  if (upwind_single_flag == 0) then
    call cgyro_upwind_r64
  else
    call cgyro_upwind_r32
  endif

end subroutine cgyro_upwind

