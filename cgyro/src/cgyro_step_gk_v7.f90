subroutine cgyro_step_gk_v7
  !!
  !! verner embedded 6-7 method
  !! 

  use timer_lib
  use mpi
  use cgyro_globals

  implicit none

  real deltah2, orig_delta_t
  real total_delta_step, delta_t_last_step, delta_t_last
  real var_error, rel_error
  real tol
  integer converged, conv, rk_MAX , iiter
  real max_scale_factor, min_scale_factor, scale_x
  real delta_t_min, delta_t_max
  real delta_x, tau
  real deltah2_min, deltah2_max
  real delta2_h_min, delta2_h_max

  real, dimension(3):: error_sum, error_x
  complex, dimension(:,:), allocatable :: h0_old

  !
  ! embedded time-advance for the distribution , verner 7(6) method
  ! 143 stages
  !
  !           z e             vpar            z e  vperp^2
  !  h = H - ----- G0 ( phi - ----- Apar ) + ----- ---------- Gperp Bpar
  !            T               c               T   omega_a c
  !

  real, parameter :: c1 = 7.0/90.0
  real, parameter :: c4 = 32.0/90.0
  real, parameter :: c5 = 32.0/90.0
  real, parameter :: c7 = 12.0/90.0
  real, parameter :: c8 = 7.0/90.0
  
  real, parameter :: a2 = 1.0/12.0
  real, parameter :: a3 = 1.0/6.0
  real, parameter :: a4 = 1.0/4.0
  real, parameter :: a5 = 3.0/4.0
  real, parameter :: a6 = 16.0/17.0
  real, parameter :: a7 = 1.0/2.0
  real, parameter :: a9 = 2.0/3.0
  
  real, parameter :: b41 = 1.0/16.0
  real, parameter :: b43 = 3.0/16.0
  real, parameter :: b51 = 21.0/16.0
  real, parameter :: b53 = -81.0/16.0
  real, parameter :: b54 = 72.0/16.0
  real, parameter :: b61 = 1344688.0/250563.0
  real, parameter :: b63 = -5127552.0/250563.0
  real, parameter :: b64 = 4096896.0/250563.0
  real, parameter :: b65 = -78208.0/250563.0
  real, parameter :: b71 = -341549.0/234624.0
  real, parameter :: b73 = 1407744.0/234624.0
  real, parameter :: b74 = -1018368.0/234624.0
  real, parameter :: b75 = 84224.0/234624.0
  real, parameter :: b76 = -14739.0/234624.0
  real, parameter :: b81 = -381875.0/136864.0
  real, parameter :: b83 = 1642368.0 / 136864.0
  real, parameter :: b84 = -1327872.0 / 136864.0
  real, parameter :: b85 = 72192.0 / 136864.0
  real, parameter :: b86 = 14739.0 / 136864.0
  real, parameter :: b87 = 117312.0 / 136864.0
  real, parameter :: b91 = -2070757.0 / 16755336.0
  real, parameter :: b93 = 9929088.0 / 16755336.0
  real, parameter :: b94 = 584064.0 / 16755336.0
  real, parameter :: b95 = 3023488.0 / 16755336.0
  real, parameter :: b96 = -447083.0 / 16755336.0
  real, parameter :: b97 = 151424.0 / 16755336.0
  
  real, parameter :: b10_1 = 130521209.0/10743824.0
  real, parameter :: b10_3 = -499279872.0/10743824.0
  real, parameter :: b10_4 = 391267968.0/10743824.0
  real, parameter :: b10_5 = 13012608.0/10743824.0
  real, parameter :: b10_6 = -3522621.0/10743824.0
  real, parameter :: b10_7 = 9033024.0/10743824.0
  real, parameter :: b10_9 = -30288492.0/10743824.0
  
  real, parameter :: e1 = -1090635.0/172448640.0
  real, parameter :: e4 = 9504768.0/172448640.0
  real, parameter :: e5 = -171816960.0/172448640.0
  real, parameter :: e6 =  72412707.0/172448640.0
  real, parameter :: e7 = -55840512.0/172448640.0
  real, parameter :: e8 = -13412672.0/172448640.0
  real, parameter :: e9 =  181730952.0/172448640.0
  real, parameter :: e10 = -21487648.0/172448640.0
  
  real, parameter :: EPS = 1.e-12
  
  allocate(h0_old(nc,nv_loc))
  
  call timer_lib_in('str')
  
  orig_delta_t = delta_t
  deltah2 = delta_t_gk
  
  deltah2_min = 1.e10
  deltah2_max = 0.
  
  iiter = 0
  min_scale_factor = 0.125
  max_scale_factor = 4.0
  
  tol = delta_t_tol
  delta_t_min = orig_delta_t*1.e-10
  delta_t_max = orig_delta_t
  delta_t_last = deltah2  ! just a dummy initializer
  delta2_h_min = 1.0
  delta2_h_max = -1.0

  
  converged = 0
  conv = 0
  rk_MAX = 10000
  
  scale_x = 1.0
  iiter = 0
  total_delta_step = 0.
  total_local_error = 0.

!$omp parallel workshare
  h0_old = h_x
!$omp end parallel workshare

  do while (total_delta_step .lt. (orig_delta_t) .and. iiter .le. rk_MAX )
     
     if ( total_delta_step + deltah2 .gt. orig_delta_t ) then
        !!
        !! last step in this interval; unless backup
        !!
        deltah2 = orig_delta_t - total_delta_step
        delta_t_last_step = deltah2
     else
        delta_t_last = deltah2
        deltah2_min = min(deltah2, deltah2_min)
        deltah2_max = max(deltah2, deltah2_max)
     endif

     if (( conv .eq. 0 ) .and. (iiter .ge. 1)) then
        !!
        !! last converged state, backing up
        !!
!$omp parallel workshare
        h0_x = h0_old
        h_x = h0_old
!$omp end parallel workshare
     else
!$omp parallel workshare
        h0_x = h_x
!$omp end parallel workshare
     endif
     
     !! if (i_proc == 0 ) write(*,*) iiter, " current time step size ", deltah2
     
     call cgyro_field_c     
     call cgyro_rhs(1)
     
     !!
     !!
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2 * a2 * rhs(ic_loc,iv_loc,1)
        enddo
     enddo
     
     call cgyro_field_c
     call cgyro_rhs(2)
     
     !!
     !! k3 = (*f)(x0+h6, *y + h6 * k2 )
     !!

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*a3*rhs(ic_loc,iv_loc,2)
        enddo
     enddo
     call cgyro_field_c
     call cgyro_rhs(3)
     

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2 * (b41*rhs(ic_loc,iv_loc,1) &
                + b43 * rhs(ic_loc,iv_loc,3))
        enddo
     enddo
     
      call cgyro_field_c
      call cgyro_rhs(4)
      
      !!
      !!
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2 * (b51*rhs(ic_loc,iv_loc,1) &
                + b53*rhs(ic_loc,iv_loc,3) &
                + b54 * rhs(ic_loc,iv_loc,4))
              enddo
     enddo

      call cgyro_field_c
      call cgyro_rhs(5)
      
      !!

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2 * ( b61*rhs(ic_loc,iv_loc,1) &
                + b63*rhs(ic_loc,iv_loc,3) &
                + b64*rhs(ic_loc,iv_loc,4) &
                + b65*rhs(ic_loc, iv_loc, 5))
        enddo
     enddo

      call cgyro_field_c
      call cgyro_rhs(6)

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2 * (b71*rhs(ic_loc,iv_loc,1) &
                + b73*rhs(ic_loc,iv_loc,3) &
                + b74*rhs(ic_loc,iv_loc,4) &
                + b75*rhs(ic_loc, iv_loc, 5) &
                + b76*rhs(ic_loc,iv_loc,6))
        enddo
     enddo

      call cgyro_field_c
      call cgyro_rhs(7)
      
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2 *(b81*rhs(ic_loc,iv_loc,1) &
                + b83*rhs(ic_loc,iv_loc,3) &
                + b84*rhs(ic_loc,iv_loc,4) &
                + b85*rhs(ic_loc, iv_loc, 5) &
                + b86*rhs(ic_loc,iv_loc,6) &
                + b87*rhs(ic_loc,iv_loc,7))
        enddo
     enddo
     call cgyro_field_c
     call cgyro_rhs(8)

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2 * ( b91*rhs(ic_loc,iv_loc,1) &
                + b93*rhs(ic_loc,iv_loc,3) &
                + b94*rhs(ic_loc,iv_loc,4) &
                + b95*rhs(ic_loc, iv_loc, 5) &
                + b96*rhs(ic_loc,iv_loc,6) &
                + b97*rhs(ic_loc,iv_loc,7))
              enddo
     enddo

     call cgyro_field_c
     call cgyro_rhs(9)
     
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2 * ( b10_1*rhs(ic_loc,iv_loc,1) &
                + b10_3*rhs(ic_loc,iv_loc,3) &
                + b10_4*rhs(ic_loc,iv_loc,4) &
                + b10_5*rhs(ic_loc, iv_loc, 5) &
                + b10_6*rhs(ic_loc,iv_loc,6) &
                + b10_7*rhs(ic_loc,iv_loc,7) &
                + b10_9*rhs(ic_loc,iv_loc,9))
              enddo
     enddo
      
      call cgyro_field_c
      call cgyro_rhs(10)
      
      !!
      !! soln order 7 
      !!

!$omp parallel do collapse(2)
      do iv_loc=1,nv_loc
         do ic_loc=1,nc
            h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                 +  deltah2*(c1*rhs(ic_loc,iv_loc,1) &
                 + c4*rhs(ic_loc,iv_loc,4) &
                 + c5*rhs(ic_loc, iv_loc, 5) &
                 + c7*rhs(ic_loc,iv_loc,7) &
                 + c8*rhs(ic_loc,iv_loc,8))
         enddo
     enddo

      call cgyro_field_c

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           rhs(ic_loc,iv_loc,1) = deltah2*(e1*rhs(ic_loc,iv_loc,1) &
                + e4*rhs(ic_loc,iv_loc,4) &
                + e5* rhs(ic_loc, iv_loc, 5) &
                + e6*rhs(ic_loc,iv_loc,6) &
                + e7*rhs(ic_loc,iv_loc,7) &
                + e8*rhs(ic_loc,iv_loc,8) &
                + e9*rhs(ic_loc,iv_loc,9) &
                + e10*rhs(ic_loc,iv_loc,10))
        enddo
     enddo

     error_sum = 0.
     error_x(1) = sum(abs(rhs(:,:,1)))
     error_x(2) = sum(abs(h_x))
      
      call MPI_ALLREDUCE(error_x, error_sum, 2, MPI_DOUBLE_PRECISION, &
           MPI_SUM, MPI_COMM_WORLD, i_err)
      
      error_x = error_sum

      delta_x = error_x(1)
      tau = tol*max(error_x(2), 1.0)

      rel_error = error_x(1)/(error_x(2)+1.d-12)
      var_error = sqrt(total_local_error + rel_error*rel_error)

      ! if mode is var_error
      !
      !     if ( var_error .lt. tol ) then
      !            if (i_proc == 0 ) &
      ! write(*,*) "after me = ", i_proc, " var error ", var_error
      !
      
      if ( error_x(1) .lt. tau ) then
         
!!         if (i_proc == 0 ) &
!!              write(*,*) "V7 deltat", deltah2, &
!!              " rel_error ", rel_error, " var error ", var_error

         total_local_error = total_local_error + rel_error*rel_error

!$omp parallel workshare
         h0_old = h0_x
!$omp end parallel workshare

         total_delta_step = total_delta_step + deltah2

         !! scale_x = max(0.95*(tol/(error_x(1) + EPS))**(1./6.), &
         !! .95*(tol/(error_x(1) + EPS))**(1./7.))

         scale_x = max((tol/(error_x(1) + EPS))**(1./6.), &
              (tol/(error_x(1) + EPS))**(1./7.))

         
         deltah2 = deltah2*max(scale_x, 1.0)
         
!!         if (( scale_x .gt. 1.0 ) .and. (i_proc == 0 )) then
!!            write(*,*) iiter, " new deltah2 ", deltah2
!!         endif
         
         converged = converged + 1
         conv = 1
      else
         conv = 0
         deltah2 = .5*deltah2
         if (i_proc .eq. 0 ) then
            write(*,*) " V7 ***  error backing up *** not converged "
            write(*,*) " new deltah2 ", deltah2,  " rel error ", rel_error
            !! " rk error x1 ", error_x(1), &
            !! " h1_x norm ", error_x(2)
         endif
      endif

      !! deltah2 = min(deltah2, delta_t_max)
      !! deltah2 = max(delta_t_min, deltah2)

      iiter = iiter + 1
      
      if ( iiter .gt. rk_MAX) then
         write(*,*) " RK V7 exceeded iteration count ", iiter
         !! should do global mpiexit
         flush(6)
         stop
      endif
   enddo

   call timer_lib_out('str')
   if(allocated(h0_old)) deallocate(h0_old) 
   call cgyro_filter

   if (var_error .gt. tol ) then
      if ( i_proc .eq. 0) then
         write(*,*) " *** HALTING local integration error *** "
         write(*,*) " total_local_error variance ", var_error
      endif
      stop
   endif
   
   delta_t_gk = delta_t_last
   total_local_error = var_error
   
!!   if ( i_proc == 0 ) then
!!        write(*,*) i_proc , " v7 converged deltah2_min, max ", &
!!             deltah2_min, deltah2_max
!!        write(*,*) i_proc , " v7 converged continuation ", delta_t_gk
!!     endif

 end subroutine cgyro_step_gk_v7

 
