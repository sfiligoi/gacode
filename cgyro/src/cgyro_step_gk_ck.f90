! Cash-Karp 6:5(4) adaptive integrator  |  multithreaded version

subroutine cgyro_step_gk_ck

  use mpi
  use timer_lib
  use cgyro_globals
  use cgyro_io
  
  implicit none

  !
  !           z e             vpar            z e  vperp^2
  !  h = H - ----- G0 ( phi - ----- Apar ) + ----- ---------- Gperp Bpar
  !            T               c               T   omega_a c
  !
  ! After time advance, we will have 
  !
  ! h    -> h_x
  ! H    -> cap_h_c
  ! phi  -> field(1)
  ! Apar -> field(2)
  ! Bpar -> field(3)

  real, parameter :: EPS = 2e-12
  
  real :: deltah2, orig_delta_t
  real :: delta_x_min, delta_x_max, local_max_error
  real :: tol_ck, total_delta_step
  real :: error_x(5), error_sum(5)
  real :: tau_ck, delta_t_old, delta_t_gk_old, delta_t_last, delta_t_last_step
  real :: last_total_error, rel_error, var_error
  real :: deltah2_min, deltah2_max, scale_x, scale_old
  
  integer :: conv, iiter, rk_count, rkMAXITER, converged, err_x
  
  orig_delta_t = delta_t       ! keep track of delta_t for exiting of subroutine  
  tol_ck = delta_t_tol
  delta_t_gk_old = delta_t_gk

  scale_old = 0.d0
  scale_x = 0.d0
  deltah2 = delta_t_gk
  delta_t_last_step = 0.d0

  delta_x_min = orig_delta_t*1.D-10
  delta_x_max = delta_t
  delta_t_old = delta_t_gk

  total_delta_step = 0.d0
  total_local_error = 0.d0
  last_total_error = 0.d0
  var_error = 0.d0
  conv = 0
  iiter = 0
  rk_count = 0
  rkMAXITER = 1000
  converged = 0
  deltah2_min = 1.d10
  deltah2_max = -1.d0

  call timer_lib_in('str_mem')
!$omp parallel do collapse(2)
  do iv_loc=1,nv_loc
     do ic_loc=1,nc
        h0_old(ic_loc,iv_loc) = h_x(ic_loc,iv_loc)
     enddo
  enddo
  call timer_lib_out('str_mem')
        
  delta_t_gk = 0.d0
  delta_t_last = deltah2
  local_max_error=0.d0
  conv = 1

  do while (total_delta_step .lt. orig_delta_t )
     
     call timer_lib_in('str')
     if ( total_delta_step + deltah2 .gt. orig_delta_t ) then
        deltah2 = orig_delta_t - total_delta_step
        delta_t_last_step = deltah2
     else
        delta_t_last = deltah2
        deltah2_min = min(deltah2, deltah2_min)
        deltah2_max = max(deltah2, deltah2_max)
     endif
     scale_old = scale_x

     if (conv == 1) then
!$omp parallel do collapse(2)
        do iv_loc=1,nv_loc
           do ic_loc=1,nc
              h0_x(ic_loc,iv_loc) = h_x(ic_loc,iv_loc)
           enddo
        enddo
     else
!$omp parallel do collapse(2)
        do iv_loc=1,nv_loc
           do ic_loc=1,nc
              h0_x(ic_loc,iv_loc) = h0_old(ic_loc,iv_loc)
              h_x(ic_loc,iv_loc) = h0_old(ic_loc,iv_loc)
           enddo
        enddo
     endif
     call timer_lib_out('str')
     
     call cgyro_field_c
     call cgyro_rhs(1)

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc, iv_loc) &
                + 0.2d0*deltah2*rhs(ic_loc, iv_loc, 1)
        enddo
     enddo
     call timer_lib_out('str')
     
     call cgyro_field_c
     call cgyro_rhs(2)

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + 1.d0/40.d0*deltah2*(3.d0*rhs(ic_loc, iv_loc, 1) &
                + 9.d0*rhs(ic_loc, iv_loc, 2))
        enddo
     enddo
     call timer_lib_out('str')

     call cgyro_field_c
     call cgyro_rhs(3)

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*( 3.d0/10.d0*rhs(ic_loc, iv_loc, 1) &
                - 9.d0/10.d0*rhs(ic_loc, iv_loc, 2) &
                + 6.d0/5.d0*rhs(ic_loc, iv_loc, 3))
        enddo
     enddo
     call timer_lib_out('str')
     
     call cgyro_field_c
     call cgyro_rhs(4) 

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*(-11.d0/54.d0 *rhs(ic_loc, iv_loc, 1) &
                + 5.d0/2.d0*rhs(ic_loc, iv_loc, 2) &
                - 70.d0/27.d0*rhs(ic_loc, iv_loc, 3) &
                + 35.d0/27.d0*rhs(ic_loc, iv_loc, 4))
        enddo
     enddo
     call timer_lib_out('str')
     
     call cgyro_field_c
     call cgyro_rhs(5)

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*(1631.d0/55296.d0*rhs(ic_loc, iv_loc, 1) &
                + 175.d0/512.d0*rhs(ic_loc, iv_loc, 2) &
                + 575.d0/13824.d0*rhs(ic_loc, iv_loc, 3) &
                + 44275.d0/110592.d0*rhs(ic_loc, iv_loc, 4) &
                + 253.d0/4096.d0*rhs(ic_loc, iv_loc, 5))
        enddo
     enddo
     call timer_lib_out('str')
     
     call cgyro_field_c
     call cgyro_rhs(6)
     
     ! Solution of order 5

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*( 37.d0/378.d0*rhs(ic_loc, iv_loc, 1) &
                + 250.d0/621.d0 * rhs(ic_loc, iv_loc, 3) &
                + 125.d0/594.d0*rhs(ic_loc, iv_loc, 4) &
                + 512.d0/1771.d0*rhs(ic_loc, iv_loc, 6))
        enddo
     enddo
     call timer_lib_out('str')


     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           rhs(ic_loc, iv_loc, 1) = deltah2*(& 
                (37.d0/378.d0-2825.d0/27648.d0)*rhs(ic_loc, iv_loc, 1) &
                + (250.d0/621.d0-18575.d0/48384.d0)*rhs(ic_loc, iv_loc, 3) &
                + (125.d0/594.d0-13525.d0/55296.d0)*rhs(ic_loc, iv_loc, 4) &
                + (-277.d0/14336.d0)*rhs(ic_loc, iv_loc, 5) &
                + (512.d0/1771.d0-1.d0/4.d0)*rhs(ic_loc, iv_loc, 6))
        enddo
     enddo

     error_sum = 0.0
     error_x(1) = sum(abs(rhs(:,:,1)))
     error_x(2) = sum(abs(h_x))
     call timer_lib_out('str')

     call timer_lib_in('str_comm')
     call MPI_ALLREDUCE(error_x, error_sum, 2, &
          MPI_DOUBLE_PRECISION,&
          MPI_SUM, MPI_COMM_WORLD, err_x)
     call timer_lib_out('str_comm')

     error_x = error_sum
     
     tau_ck = tol_ck*max(error_x(2), 1.0d0)
     
     rel_error = error_x(1)/(error_x(2) +1.d-12)

     var_error = sqrt(total_local_error + rel_error*rel_error)
     
     if ( var_error .lt. tol_ck ) then
        call cgyro_field_c

        total_delta_step = total_delta_step + deltah2
        total_local_error = total_local_error + rel_error*rel_error

        scale_x = max((tol_ck/(error_x(1) + EPS)*1./delta_t)**(.2d0), &
             (tol_ck/(error_x(1) + EPS)*1.d0/delta_t)**(.25d0))

        deltah2 = deltah2*max(1.d0, min(6.d0,scale_x))
        
        local_max_error = max(local_max_error, rel_error)
        conv = 1
        converged=converged+1

        call timer_lib_in('str_mem')
!$omp parallel do collapse(2)
        do iv_loc=1,nv_loc
           do ic_loc=1,nc
              h0_old(ic_loc,iv_loc) = h0_x(ic_loc, iv_loc)
           enddo
        enddo
        call timer_lib_out('str_mem')

     else
        deltah2 = .5d0*deltah2
        conv = 0
     endif

     deltah2 = min(deltah2, delta_x_max)
     deltah2 = max(delta_x_min, deltah2)

     iiter = iiter + 1
    
     if ( iiter .gt. rkMAXITER) then
        call cgyro_error('Cash-Carp step exceeded max iteration count')
        return
     endif

  enddo

  delta_t_gk = max(delta_t_last, deltah2)

  if (delta_t_last_step == 0.0) delta_t_last_step = delta_t_last

  if ( delta_t_last_step < 0.1*delta_t_gk ) then
     delta_t_gk = delta_t_last + delta_t_last_step
  else
     if ( delta_t_last_step/iiter < 0.1*delta_t_gk ) then
        delta_t_gk = delta_t_gk + delta_t_last_step/iiter
     endif
  endif
  delta_t_gk = min(delta_t, delta_t_gk)
  total_local_error = var_error

end subroutine cgyro_step_gk_ck
