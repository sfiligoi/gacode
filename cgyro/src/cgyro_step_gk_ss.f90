subroutine cgyro_step_gk_ss

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none

  double precision delta_xh2, delta_x_ht, delta_x_old1,newstep, ctheta, orig_delta_x_t
  double precision delta_xht, tolerance, Rdummy, dummy2, dummy3, total_delta_step
  double precision Rdummy1, old_delta_x, error_sum(2), delta_x_min, delta_x_max
  double precision delta_x, tau, deltah2, delta_t_last, rel_error, var_error
  double precision local_max_error, delta_t_last_step
  double precision deltah2_min, deltah2_max
  
  complex, dimension(:,:), allocatable :: h0_old_orig, h_rk4, h_rk5
  complex, dimension(:,:), allocatable :: h_rk6, rk_error_x
  complex, dimension(:,:), allocatable :: h0_old
  double precision error_x(2), tol, orig_t_current
  integer converged, conv, rkcount, rk_MAX , deadcount, iiter
  double precision scale_x
  ! double complex h0_old(nc,nv_loc)

  ! RK5(4) TS.x
  !
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

  double precision, parameter :: c2   = 16./105.
  double precision, parameter :: c3   = 8./35.
  double precision, parameter :: c4   = 9./20.
  double precision, parameter :: c5   = 2./3.
  double precision, parameter :: c6   = 7./9.
  double precision, parameter :: c7   = 1.0



  double precision, parameter :: a21  = 16./105.
  double precision, parameter :: a31  = 2./35.
  double precision, parameter :: a32  = 6./35.
  double precision, parameter :: a41  = 8793./40960.
  double precision, parameter :: a42  = -5103./8192.
  double precision, parameter :: a43 = 17577./20480.
  double precision, parameter :: a51 = 347./1458.
  double precision, parameter :: a52 = -7./20.
  double precision, parameter :: a53 = 3395./10044.
  double precision, parameter :: a54 = 49792./112995.
  double precision, parameter :: a61 = -1223224109959./9199771214400.
  double precision, parameter :: a62 = 1234787701./2523942720.
  double precision, parameter :: a63 = 568994101921./3168810084960.
  double precision, parameter :: a64  = -105209683888./891227836395.
  double precision, parameter :: a65  = 9./25.
  double precision, parameter :: a71  = 2462504862877./8306031988800.
  double precision, parameter :: a72  = -123991./287040.
  double precision, parameter :: a73  = 106522578491./408709510560.
  double precision, parameter :: a74  = 590616498832./804646848915.
  double precision, parameter :: a75  = -319138726./534081275.
  double precision, parameter :: a76  = 52758./71449.
  
  double precision, parameter :: b1 = 1093./15120.
  double precision, parameter :: b2 = 0.
  double precision, parameter :: b3 = 60025./190992.
  double precision, parameter :: b4 = 3200./20709.
  double precision, parameter :: b5 = 1611./11960.
  double precision, parameter :: b6 = 712233./2857960.
  double precision, parameter :: b7 = 3./40.

  double precision, parameter :: b1p = 84018211./991368000.
  double precision, parameter :: b2p  = 0.d0
  double precision, parameter :: b3p  = 92098979./357791680.
  double precision, parameter :: b4p  = 17606944./67891005.
  double precision, parameter :: b5p  = 3142101./235253200.
  double precision, parameter :: b6p  = 22004596809./70270091500.
  double precision, parameter :: b7p  = 9./125.

  double precision, parameter :: EPS  = 2.2d-12  !! 2.e-16 triggers ieee error
  
  !!  allocate(h_rk4(nc,nv_loc))
  !! allocate(h_rk5(nc,nv_loc))
  !! allocate(h_rk6(nc,nv_loc))
  
  allocate(h0_old(nc,nv_loc))
  
  !! allocate(h_rk6(nc,nv_loc))
  !! allocate(rk_error_x(nc,nv_loc))

  call timer_lib_in('str')
  
  !! call cgyro_step_gkssp

  orig_delta_x_t = delta_t

  iiter = 0
  total_delta_step = 0.d0
  total_local_error = 0.d0

  !! tol = .001d0
  !! if ( orig_delta_x_t .gt. 1 ) tol = tol/orig_delta_x_t
  !! if ( orig_delta_x_t .gt. 1 ) tol = tol*orig_delta_x_t

  tol = delta_t_tol

  deltah2 = delta_t_gk
  
  if ( delta_t_gk .lt. 1.d-10) then
     deltah2 = orig_delta_x_t
     delta_t_gk = deltah2
  endif
  
  delta_x_min = 1.d-10*orig_delta_x_t
  delta_x_max = orig_delta_x_t

  converged = 0
  conv = 0
  
  rk_MAX = 10000

  h0_old = h_x
  conv = 0
  delta_t_gk = 0.d0
  deltah2_min = 1.d10
  deltah2_max = 0.

  do while (total_delta_step .lt. orig_delta_x_t )

    
     if ( total_delta_step + deltah2 .gt. orig_delta_x_t ) then
        deltah2 = orig_delta_x_t - total_delta_step
        delta_t_last_step = deltah2
        if ( deltah2 .lt. 1.d-9 ) goto 1111   !! abandon
     else
        delta_t_gk = deltah2+delta_t_gk
        delta_t_last = deltah2
        deltah2_min = min(deltah2, deltah2_min)
        deltah2_max = max(deltah2, deltah2_max)
     endif
     
     if (( conv .eq. 0 ) .and. (iiter .ge. 1)) then
        !!
        !! last converged state
        !!
        h0_x = h0_old
        h_x = h0_old
        call cgyro_field_c
     else
        h0_x = h_x
        call cgyro_field_c
     endif

     if ( i_proc == 0 ) write(*,*) i_proc, " step ", iiter, " ts deltah2 ", deltah2

     !! if ( iiter .gt. 0 ) call cgyro_field_c
     
     !
     ! Stage 1
     !
     
     call cgyro_rhs(1)
     !!!$OMP PARALLEL
     h_x = h0_x + deltah2 * a21 * rhs(:,:,1)
     !!!$OMP END PARALLEL
     call cgyro_field_c
     
     ! Stage 2 ! k2

     call cgyro_rhs(2)

     !!     
     !!     h_x = h0_x + 3./32.*delta_x_t*(rhs(:,:,1) + 3.*rhs(:,:,2))
     !!
     !!!$OMP PARALLEL
     h_x = h0_x + deltah2*(a31*rhs(:,:,1) + a32*rhs(:,:,2))
     !!!$OMP END PARALLEL
     call cgyro_field_c
     
     ! Stage 3
     call cgyro_rhs(3)
     ! k4
     !!!$OMP PARALLEL
     h_x = h0_x + deltah2*(a41*rhs(:,:,1) &
          + a42*rhs(:,:,2) + a43*rhs(:,:,3))
     !!!$OMP END PARALLEL
     !
     call cgyro_field_c
     
     ! Stage 4
     call cgyro_rhs(4)
     !  k5
     !!!$OMP PARALLEL
     h_x = h0_x+deltah2*(a51*rhs(:,:,1) + a52*rhs(:,:,2) &
          + a53*rhs(:,:,3) + a54*rhs(:,:,4))
     !!!$OMP END PARALLEL
     call cgyro_field_c
     
     !  stage 5
     
     call cgyro_rhs(5)
     !!!$OMP PARALLEL
     h_x = h0_x + deltah2*(a61*rhs(:,:,1) + a62*rhs(:,:,2) &
          + a63*rhs(:,:,3) + a64*rhs(:,:,4) + a65*rhs(:,:,5))
     !!!$OMP END PARALLEL
     call cgyro_field_c

     ! stage 6
     call cgyro_rhs(6)
     !!!$OMP PARALLEL
     h_x = h0_x + deltah2*(a71*rhs(:,:,1) + a72*rhs(:,:,2) &
          + a73*rhs(:,:,3) + a74*rhs(:,:,4) &
          + a75*rhs(:,:,5) + a76*rhs(:,:,6))
     !!!$OMP END PARALLEL
     
     call cgyro_field_c

     call cgyro_rhs(7)
     
     !! soln = h_x of order 5
     !! error_x(2) = sum(abs(h0_x))
     !!!$OMP PARALLEL
     h_x = h0_x + deltah2*(b1*rhs(:,:,1) &
          +b3*rhs(:,:,3) + b4*rhs(:,:,4) &
          + b5*rhs(:,:,5) + b6*rhs(:,:,6) + b7*rhs(:,:,7))
     !!!$OMP END PARALLEL

     call cgyro_field_c
     
     !!!$OMP PARALLEL
     rhs(:,:,1) = deltah2*((b1-b1p)*rhs(:,:,1) &
          + (b3-b3p)*rhs(:,:,3) &
          + (b4-b4p)*rhs(:,:,4) &
          + (b5-b5p)*rhs(:,:,5) + (b6-b6p)*rhs(:,:,6) &
          + (b7-b7p)*rhs(:,:,7))
     !!!$OMP END PARALLEL

     !! rk_error = -1.*rk_error

     !! error_sum = 0.d0
     
     !! error_x(1) = sum(abs(rk_error_x))
     error_x(1) = sum(abs(rhs(:,:,1)))
     error_x(2) = sum(abs(h_x))     

     !! write(*,*) "me = ",  i_proc, " error before ", error, " error_sum ", error_sum

     call MPI_ALLREDUCE(error_x, error_sum, 2, MPI_DOUBLE_PRECISION,&
          MPI_SUM, MPI_COMM_WORLD, i_err)
     
!!     write(*,*) "after me = ", i_proc, " error out ", error, " error_sum ", error_sum

     error_x(1) = error_sum(1)! local truncation error
     error_x(2) = error_sum(2)! local h0

     !!
     !! delta_x = 0.84*(tol/(error+EPS)**(1./5.))  !! prevent div 0
     !!
     
     delta_x = error_x(1)
     tau = tol * max(error_x(2), 1.d0)

     !! tau = tol*error_x(2)
     
     !! scale_x = 0.9*min((tol/(error_x(1) + EPS))**(1./5.), 10.)
     !! scale_x = 0.9*min((tol/(error_x(1) + EPS))**(1./6.), delta_x_max)
     !! delta_x = ((tol*delta_x_t)/(2.*error+EPS))**(1./5.)

     !! if ( delta_x .lt. tau ) then

     rel_error = error_x(1)/(error_x(2)+1.d-12)
     var_error = sqrt(total_local_error + rel_error*rel_error)

     !! if ( error_mode .eq. 0 ) then local error
     !! if ( error_mode .eq. 1 ) then variance error
     
     ! if ( var_error .lt. tol ) then
     ! if ( i_proc == 0 ) write(*,*) " variance error mode ", var_error, total_local_error

     if ( error_x(1) .lt. tau ) then
        
        if ( i_proc == 0 ) write(*,*) " local error mode ", rel_error, " variance error", var_error

        h0_old = h0_x
        
        !! h_x = h0_x + deltah2*(b1*rhs(:,:,1)+b3*rhs(:,:,3) + b4*rhs(:,:,4) &
        !! + b5*rhs(:,:,5) + b6*rhs(:,:,6))
        
        converged = converged + 1
        conv = 1
        total_delta_step = total_delta_step + deltah2
        total_local_error = total_local_error + rel_error*rel_error

        !! t_current = t_current + delta_t
        !! scale_x = min((tol/(error_x(1) + EPS))**(1./6.), delta_x_max)
        !! scale_x = 0.9*min((tol/(error_x(1) + EPS))**(1./5.), 10.)
        
        !!delta_t = min(delta_x_max, &
        !! .8d0*deltah2*(tau/(delta_x + 1.e-14))**(1.d0/6.d0))

        !! deltah2 = 3.d0*.95d0*deltah2*(tol/(delta_x + 2.2d-12))**(.2d0)
        !! check
        
        scale_x = max(0.95d0*deltah2*(tol/(delta_x + 1.0d-12))**(.2d0), &
             0.95d0*deltah2*(tol/(delta_x + 1.0d-12))**(.25d0))
        
        deltah2 = deltah2*max(1., scale_x)
        local_max_error = max(local_max_error, rel_error)
     else
        !! h_x = h0_old
        conv = 0
        deltah2 = .5d0*deltah2          !! interpolate?
        !! deltah2 = (3./5.)*deltah2          !! interpolate?
        !! deltah2 = (2.d0/3.d0)*deltah2          !! interpolate?
        if (i_proc .eq. 0 ) then
           write(*,*) " ts ***  error backing up *** not converged " ! , &
!                " total delta_x step ", total_delta_step, " delta_x2 ", deltah2
        endif
     endif
     
     !! delta_x = 0.9*min((tol/(error + EPS))**(1./5.), 10.)
     !! deltah2 = scale_x*deltah2
     
     deltah2 = min(deltah2, delta_x_max)
     deltah2 = max(delta_x_min, deltah2)

     iiter = iiter + 1

     if (deltah2 .lt. 1.d-8) then
        if ( i_proc .eq. 0 ) &
             write(*,*) " ******* Stopping due to small substep size ", deltah2
        stop
     endif
     
     if ( iiter .gt. rk_MAX ) then
        if ( i_proc .eq. 0 ) then
           write(*,*) " rk ts max count  exceeded ", iiter
        endif
        goto 1111
     endif
  enddo
1111 continue
  call timer_lib_out('str')  
  call cgyro_field_c

!!  delta_t_gk = max(.5d0*(delta_t_gk + deltah2), delta_t/converged)

  !! delta_t_gk = max(delta_t_gk, delta_t/converged)

  !! delta_t_tol = 2.d0*min(adapt_tol, error_tol)/converged

  delta_t_gk = delta_t_last

  !! gif3 delta_t_gk = delta_t_gk/converged

  !! gif_3 delta_t_tol = 2.d0*max(total_local_error, &
       !! tol/max(converged-1, 1)
  
  !! delta_t_tol = max(delta_t_tol, min(adapt_tol, error_tol))/max(converged-1, 1)

  !! delta_t_tol = max(delta_t_tol, min(adapt_tol, error_tol))/max(converged-1, 1)

  goto 2222
  if ( i_proc .eq. 0 ) then
     write(*,*) " local error ", total_local_error
     write(*,*) " delta_t_gk ", delta_t_gk
     write(*,*) " variance local error sqrt of local error ", var_error
  endif
2222 continue
  
  if (var_error .gt. tol ) stop
  
  !! t_current = orig_t_current
  total_local_error = var_error
  if ( i_proc == 0 ) write(*,*) i_proc , " ts deltah2_min, max ", deltah2_min, deltah2_max, converged
  
  call cgyro_filter
  
  ! Filter special spectral components
  
  if(allocated(h0_old)) deallocate(h0_old)  
  
end subroutine cgyro_step_gk_ss
