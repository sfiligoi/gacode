!------------------------------------------------------------------------------
! cgyro_error_estimate.f90
!
! PURPOSE:
!  Compute two time-integration error estimates based on 
!  (1) difference between computed field and quadratic extrapolation, 
!  (2) 3rd-order estimate of collisionless error based on RK algebra
!------------------------------------------------------------------------------

subroutine cgyro_error_estimate

  use mpi
  use cgyro_globals
  use cgyro_io
  use timer_lib

  implicit none

  integer :: i_f
  real, dimension(2) :: norm_loc,norm
  real, dimension(2) :: pair_loc,pair
  real, dimension(2) :: error_loc

  real :: norm_loc_s,error_loc_s,h_s,r_s

#ifdef _OPENACC
  ! launch Estimate of collisionless error via 3rd-order linear estimate async ahead of time on GPU
  ! CPU-only code will work on it later
  h_s=0.0
  r_s=0.0
!$acc parallel loop collapse(2) independent present(h_x,rhs(:,:,1)) reduction(+:h_s,r_s) async(2)
  do iv_loc=1,nv_loc
     do ic=1,nc
        h_s = h_s + abs(h_x(ic,iv_loc))
        r_s = r_s + abs(rhs(ic,iv_loc,1))
     enddo
  enddo
#endif

  call timer_lib_in('field')

  norm_loc_s = 0.0
  error_loc_s = 0.0

  ! field_olds are always only in system memory... too expensive to keep in GPU memory
  ! assuming field was already synched to system memory
!$omp parallel do collapse(2) reduction(+:norm_loc_s,error_loc_s)
  do ic=1,nc
     do i_f=1,n_field

        ! 1. Estimate of total (field) error via quadratic interpolation

        field_loc(i_f,ic) = 3*field_old(i_f,ic)-3*field_old2(i_f,ic)+field_old3(i_f,ic)
        field_dot(i_f,ic) = (3*field(i_f,ic)-4*field_old(i_f,ic)+field_old2(i_f,ic))/(2*delta_t)

        ! Define norm and error for each mode number n
        norm_loc_s  = norm_loc_s  + abs(field(i_f,ic))
        error_loc_s = error_loc_s + abs(field(i_f,ic)-field_loc(i_f,ic))

        ! save old values for next iteration
        field_old3(i_f,ic) = field_old2(i_f,ic)
        field_old2(i_f,ic) = field_old(i_f,ic)
        field_old(i_f,ic)  = field(i_f,ic)
     enddo
  enddo

  norm_loc(1)  = norm_loc_s
  error_loc(1) = error_loc_s

  ! JC: Optimize?
  cap_h_c_dot = (3*cap_h_c-4*cap_h_c_old+cap_h_c_old2)/(2*delta_t)
  cap_h_c_old = cap_h_c
  cap_h_c_old2 = cap_h_c_old

  call timer_lib_out('field')

  ! 2. Estimate of collisionless error via 3rd-order linear estimate

  call timer_lib_in('str')

#ifdef _OPENACC
  ! wait for the async GPU compute to be completed
!$acc wait(2)
#else
  h_s=0.0
  r_s=0.0
!$omp parallel do collapse(2) reduction(+:h_s,r_s)
  do iv_loc=1,nv_loc
     do ic=1,nc
        h_s = h_s + abs(h_x(ic,iv_loc))
        r_s = r_s + abs(rhs(ic,iv_loc,1))
     enddo
  enddo
#endif

  pair_loc(1) = h_s
  pair_loc(2) = r_s

  call timer_lib_out('str')

  call timer_lib_in('str_comm')

  ! sum over velocity space
  call MPI_ALLREDUCE(pair_loc,&
       pair,&
       2,&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  norm_loc(2) = pair(1)
  error_loc(2) = pair(2)

  ! Get sum of all errors
  call MPI_ALLREDUCE(error_loc, &
       integration_error, &
       2, &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  ! Get sum of all norms
  call MPI_ALLREDUCE(norm_loc, &
       norm, &
       2, &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  call timer_lib_out('str_comm')

  ! 1: total (field) error
  ! 2: collisionless error
  integration_error = integration_error/norm

  ! Trigger code shutdown on large collisionless error
  if (integration_error(2) > 5e2*error_tol .and. i_time > 2) then
     call cgyro_error('Integration error exceeded limit.')
  endif

end subroutine cgyro_error_estimate
