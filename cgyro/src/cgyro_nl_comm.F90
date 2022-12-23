!-----------------------------------------------------------------
! cgyro_nl_comm.F90
!
! PURPOSE:
!  Nonlinear communication routines
!-----------------------------------------------------------------

module cgyro_nl_comm

  implicit none

contains

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
!$acc parallel loop collapse(3) gang vector independent private(iexch) &
!$acc&         present(ic_c,h_x,fpack) default(none)
#else
!$omp parallel do collapse(3) private(iexch)
#endif
  do it=1,n_theta
   do iv_loc_m=1,nv_loc
     do ir=1,n_radial
       iexch = iv_loc_m + (it-1)*nv_loc
       fpack(ir,iexch) = h_x(ic_c(ir,it),iv_loc_m,my_toroidal)
     enddo
   enddo
  enddo

#ifdef _OPENACC
!$acc parallel loop independent present(fpack) if ( (nv_loc*n_theta) < (nsplit*n_toroidal) )
#endif
  do iexch=nv_loc*n_theta+1,nsplit*n_toroidal
      fpack(1:n_radial,iexch) = (0.0,0.0)
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

subroutine cgyro_nl_fftw_comm1_r(ij)
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: ij
  !-----------------------------------

  integer :: ir,it,iv_loc_m,ic_loc_m
  integer :: iexch
  complex :: my_psi
  real :: psi_mul

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc(f_nl,fpack)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

  psi_mul = ((q*rho/rmin)*(2*pi/length))

#ifdef _OPENACC
!$acc parallel loop collapse(3) gang vector independent private(iexch,ic_loc_m,my_psi) &
!$acc&         present(ic_c,px,rhs,fpack) default(none)
#else
!$omp parallel do collapse(3) private(iexch,ic_loc_m,my_psi)
#endif
  do it=1,n_theta
     do iv_loc_m=1,nv_loc
        do ir=1,n_radial
           iexch = iv_loc_m + (it-1)*nv_loc
           ic_loc_m = ic_c(ir,it)
           if ( (my_toroidal == 0) .and.  (ir == 1 .or. px(ir) == 0) ) then
              ! filter
              my_psi = (0.0,0.0)
           else
              my_psi = fpack(ir,iexch)
           endif
           ! RHS -> -[f,g] = [f,g]_{r,-alpha}
           rhs(ic_loc_m,iv_loc_m,my_toroidal,ij) = rhs(ic_loc_m,iv_loc_m,my_toroidal,ij)+psi_mul*my_psi
        enddo
     enddo
  enddo

  call timer_lib_out('nl')
end subroutine cgyro_nl_fftw_comm1_r

!
! Comm2 is a transpose
! First half of the transpose is done locally
!  from (field,:,n_radial) -> (field,n_radial, :)
! Then AlltoAll finishes the transpose
! 

subroutine cgyro_nl_fftw_comm2_async
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,il,it_loc
  integer :: iexch

  call timer_lib_in('nl_mem')

#ifdef _OPENACC
!$acc parallel loop gang collapse(2) independent private(it) &
!$acc&         present(ic_c,it_f,field,gpack) default(none)
#else
!$omp parallel do collapse(2) private(it,ir)
#endif
  do il=1,n_toroidal
   do it_loc=1,n_jtheta
     it = it_f(it_loc,il)
     if (it == 0) then
        ! just padding
        gpack(1:n_field,1:n_radial,it_loc,il) = (0.0,0.0)
     else
!$acc loop vector
        do ir=1,n_radial
           gpack(1:n_field,ir,it_loc,il) = field(1:n_field,ic_c(ir,it),my_toroidal)
        enddo
     endif
   enddo
  enddo

  call parallel_slib_f_fd_async(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)

  call timer_lib_out('nl_mem')

end subroutine cgyro_nl_fftw_comm2_async

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_nl_fftw_comm2_test
  use parallel_lib
  use cgyro_globals

  implicit none

  call parallel_slib_test(g_req)

end subroutine cgyro_nl_fftw_comm2_test

end module cgyro_nl_comm

