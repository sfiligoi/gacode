!-----------------------------------------------------------------
! cgyro_nl_comm.F90
!
! PURPOSE:
!  Nonlinear communication routines
!-----------------------------------------------------------------

!
! Comm is a transpose
! First half of the transpose is done locally
!  from (theta,radial,nv_loc) -> (radial, theta, nv_lov)
! Then AlltoAll finishes the transpose
!  from (radial, theta, nv_loc_1, nv_loc_2) x toroidal -> (radial, theta, nv_loc_1 , toroidal) x nv_loc_2
! Implies nv_loc_2 == toroidal
!

! NOTE: call cgyro_nl_fftw_comm1/2_async before cgyro_nl_fftw
subroutine cgyro_nl_fftw_comm1_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m
  integer :: iexch

  call timer_lib_in('nl_mem')

#ifdef _OPENACC
!$acc parallel loop gang independent private(it,iv_loc_m) &
!$acc&         present(iv_e,it_e,ic_c,h_x,fpack) default(none)
#else
!$omp parallel do private(iv_loc_m,it,ir)
#endif
  do iexch=1,nsplit*n_toroidal
     it = it_e(iexch)
     iv_loc_m = iv_e(iexch)
     if (iv_loc_m == 0) then
        ! padding
        fpack(1:n_radial,iexch) = (0.0,0.0)
     else
!$acc loop vector
        do ir=1,n_radial
           fpack(ir,iexch) = h_x(ic_c(ir,it),iv_loc_m)
        enddo
     endif
  enddo

  call parallel_slib_f_nc_async(fpack,f_nl,f_req)

  call timer_lib_out('nl_mem')

end subroutine cgyro_nl_fftw_comm1_async

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_nl_fftw_comm1_test
  use parallel_lib
  use cgyro_globals

  implicit none

  call parallel_slib_test(f_req)

end subroutine cgyro_nl_fftw_comm1_test

subroutine cgyro_nl_fftw_comm1_r
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m,ic_loc_m
  integer :: iexch
  complex :: val

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc(f_nl,fpack)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl_mem')

#ifdef _OPENACC
!$acc parallel loop gang independent private(it,iv_loc_m) &
!$acc&         present(iv_e,it_e,ic_c,px,psi,fpack) default(none)
#else
!$omp parallel do private(iv_loc_m,it,ir,ic_loc_m,val)
#endif
  do iexch=1,nsplit*n_toroidal
     it = it_e(iexch)
     iv_loc_m = iv_e(iexch)
     if (iv_loc_m /= 0 ) then ! else it is padding and can be ignored
!$acc loop vector private(ic_loc_m,val)
        do ir=1,n_radial
           ic_loc_m = ic_c(ir,it)
           if ( (my_toroidal == 0) .and.  (ir == 1 .or. px(ir) == 0) ) then
              ! filter
              val = (0.0,0.0)
           else
              val = fpack(ir,iexch)
           endif
           psi(ic_loc_m,iv_loc_m) = val
        enddo
     endif
  enddo

  call timer_lib_out('nl_mem')
end subroutine cgyro_nl_fftw_comm1_r

subroutine cgyro_nl_fftw_comm2_async

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m
  integer :: iexch

  call timer_lib_in('nl_mem')

#ifdef _OPENACC
!$acc parallel loop gang independent private(it,iv_loc_m) &
!$acc&         present(iv_e,it_e,ic_c,jvec_c,field,gpack) default(none)
#else
!$omp parallel do private(iv_loc_m,it,ir,ic)
#endif
  do iexch=1,nsplit*n_toroidal
     it = it_e(iexch)
     iv_loc_m = iv_e(iexch)
     if (iv_loc_m == 0) then
        ! padding
        gpack(1:n_radial,iexch) = (0.0,0.0)
     else
!$acc loop vector private(ic)
        do ir=1,n_radial
           ic = ic_c(ir,it)
           gpack(ir,iexch) = sum( jvec_c(:,ic,iv_loc_m)*field(:,ic))
        enddo
     endif
  enddo

  call parallel_slib_f_nc_async(gpack,g_nl,g_req)

  call timer_lib_out('nl_mem')

end subroutine cgyro_nl_fftw_comm2_async

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_nl_fftw_comm2_test
  use parallel_lib
  use cgyro_globals

  implicit none

  call parallel_slib_test(g_req)

end subroutine cgyro_nl_fftw_comm2_test

