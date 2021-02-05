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
  complex :: res_loc_one, res_loc_two

  call timer_lib_in('str')

!$acc parallel loop collapse(2) gang &
!$acc&         private(res_loc_one,iv) &
!$acc&         present(g_x,upfac1,is_v,upwind_res_loc) default(none)
  do is=1,n_species
     do ic=1,nc
       res_loc_one = (0.0,0.0)
       res_loc_two = (0.0,0.0)

!$acc loop vector private(iv_loc) reduction(+:res_loc_one)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             res_loc_one = res_loc_one+upfac1(ic,iv_loc,1)*g_x(ic,iv_loc)
          endif
       enddo
       upwind_res_loc(ic,is,1) = res_loc_one

!$acc loop vector private(iv_loc) reduction(+:res_loc_two)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             res_loc_two = res_loc_two+upfac1(ic,iv_loc,2)*g_x(ic,iv_loc)
          endif
       enddo
       upwind_res_loc(ic,is,2) = res_loc_two
    enddo

  enddo

  call timer_lib_out('str')

  call timer_lib_in('str_comm')

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update host(upwind_res_loc)
#else
!$acc host_data use_device(upwind_res_loc,upwind_res)
#endif

#ifdef SUMMIT
  call MPI_ALLREDUCE(upwind_res_loc(:,:,:),&
       upwind_res(:,:,:),&
       2*size(upwind_res(:,:,:)),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
#else
  call MPI_ALLREDUCE(upwind_res_loc(:,:,:),&
       upwind_res(:,:,:),&
       size(upwind_res(:,:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
#endif

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update device(upwind_res)
#else
!$acc end host_data
#endif

  call timer_lib_out('str_comm')

  call timer_lib_in('str')

!$acc parallel loop collapse(2) independent &
!$acc&         present(is_v,ix_v,ie_v,xi,vel,upfac2,g_x,upwind_res) &
!$acc&         private(iv_loc,is,ix,ie) default(none)
  do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        g_x(ic,iv_loc) = abs(xi(ix))*vel(ie)*g_x(ic,iv_loc) &
             -upfac2(ic,iv_loc,1)*upwind_res(ic,is,1) &
             -upfac2(ic,iv_loc,2)*upwind_res(ic,is,2) 
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
  complex(KIND=REAL32) :: res_loc_one, res_loc_two

  call timer_lib_in('str')

!$acc parallel loop collapse(2) gang &
!$acc&         private(res_loc_one,iv) &
!$acc&         present(g_x,upfac1,is_v,upwind32_res_loc) default(none)
  do is=1,n_species
     do ic=1,nc
       res_loc_one = (0.0,0.0)
       res_loc_two = (0.0,0.0)

!$acc loop vector private(iv_loc) reduction(+:res_loc_one)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             res_loc_one = res_loc_one+upfac1(ic,iv_loc,1)*g_x(ic,iv_loc)
          endif
       enddo
       upwind32_res_loc(ic,is,1) = res_loc_one

!$acc loop vector private(iv_loc) reduction(+:res_loc_two)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             res_loc_two = res_loc_two+upfac1(ic,iv_loc,2)*g_x(ic,iv_loc)
          endif
       enddo
       upwind32_res_loc(ic,is,2) = res_loc_two
    enddo

  enddo

  call timer_lib_out('str')

  call timer_lib_in('str_comm')

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update host(upwind32_res_loc)
#else
!$acc host_data use_device(upwind32_res_loc,upwind32_res)
#endif

#ifdef SUMMIT
  call MPI_ALLREDUCE(upwind32_res_loc(:,:,:),&
       upwind32_res(:,:,:),&
       2*size(upwind32_res(:,:,:)),&
       MPI_REAL,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
#else
  call MPI_ALLREDUCE(upwind32_res_loc(:,:,:),&
       upwind32_res(:,:,:),&
       size(upwind32_res(:,:,:)),&
       MPI_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
#endif

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update device(upwind32_res)
#else
!$acc end host_data
#endif

  call timer_lib_out('str_comm')

  call timer_lib_in('str')

!$acc parallel loop collapse(2) independent &
!$acc&         present(is_v,ix_v,ie_v,xi,vel,upfac2,g_x,upwind32_res) &
!$acc&         private(iv_loc,is,ix,ie) default(none)
  do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        g_x(ic,iv_loc) = abs(xi(ix))*vel(ie)*g_x(ic,iv_loc) &
             -upfac2(ic,iv_loc,1)*upwind32_res(ic,is,1) &
             -upfac2(ic,iv_loc,2)*upwind32_res(ic,is,2)
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

