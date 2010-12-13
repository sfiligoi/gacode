!---------------------------------------------------
! gyro_collision.f90 
!
! PURPOSE:
!  Collision step for real (rather than model 
!  Krook) collisions.
!---------------------------------------------------

subroutine gyro_collision

  use gyro_globals
  use gyro_pointers

  !----------------------------------------
  implicit none
  !
  integer :: ic
  integer :: n_i
  integer :: n_j
  integer :: n_k
  integer :: n_d
  integer, external :: parallel_dim
  !
  real :: CPU_temp
  !----------------------------------------

  !---------------------
  ! BEGIN Collision step
  !---------------------

  call proc_time(CPU_Ct_in)

  !-----------------------------------------------------
  ! Step 1: 
  !
  ! Advance A_parallel (A -> A+dt*A_dot) using electron 
  ! collisions only.  This is consistent with ion 
  ! collsion momentum conservation.
  !
  if (nu_s(indx_e,ir_norm) > 0.0 .and. n_field > 1) then

     ! Only electron collisions are used to update A.
     ! New value of gyro_u also computed.
     call get_ampere_coll

  endif
  !-----------------------------------------------------

  do ic=1,n_coll

     is = c_map(ic)

     !-----------------------------------------------------
     ! Step 2:
     !
     ! Calculate array f_coll to be passed to do_collision:
     !
     ! f = h-z*alpha*vp*A
     !
     ! df/dt = -z*alpha*vp*dA/dt + C(f)
     !
     ! Semi-implicit time-advance:
     !
     ! fb - f = [-z*alpha*vp*dA/dt + C(fb)]*dt
     !
     ! (1-dt*C)fb = f-dt*(z*alpha*vp*dA/dt)
     !            = h-z*alpha*vp*(A+dt*dA/dt)
     !            = h-z*alpha*vp*Ab
     !            -> f_coll
     !
     if (n_field == 2) then

        p_nek_loc = 0
        do p_nek=1+i_proc_1,n_nek_1,n_proc_1

           p_nek_loc = p_nek_loc+1

           ie = nek_e(p_nek)  
           k  = nek_k(p_nek)   

           do m=1,n_stack
              do i=1,n_x

                 f_coll(m,i,p_nek_loc) = h(m,i,p_nek_loc,is)-z(is)*&
                      alpha_s(is,i)*v_para(m,i,p_nek_loc,is)*&
                      field_tau(m,i,p_nek_loc,2)

              enddo ! i
           enddo ! m

        enddo ! p_nek

     else

        f_coll(:,:,:) =  h(:,:,:,is)

     endif
     !------------------------------------------------------

     call proc_time(CPU_Ct_out)
     CPU_temp = CPU_Ct_out-CPU_Ct_in

     !-----------------------------------------------
     ! Step 3:
     ! 
     ! Implicit collision solve:
     !
     ! Permute p_nek,i -> p_ine,k
     n_i = n_n_1*n_energy
     n_j = n_lambda
     n_k = n_x
     n_d = parallel_dim(n_k*n_i,n_proc_1)
     !
     call rTRANSP_INIT(n_i,n_j,n_k,NEW_COMM_1)
     do m=1,n_stack
        call rTRANSP_DO(f_coll(m,:,:),h_C(m,:,:))
     enddo
     call rTRANSP_CLEANUP
     !
     call proc_time(CPU_Ct_in)
     !
     ! Full collision step (pitch-angle scattering 
     ! plus specified restoring terms).
     !
     call gyro_collision_kernel(ic)
     !
     call proc_time(CPU_Ct_out)
     CPU_temp = CPU_temp+(CPU_Ct_out-CPU_Ct_in)
     !
     ! Permute p_ine,k -> p_nek,i
     n_i = n_x
     n_j = n_n_1*n_energy
     n_k = n_lambda
     n_d = parallel_dim(n_j*n_k,n_proc_1)
     !
     call fTRANSP_INIT(n_i,n_j,n_k,NEW_COMM_1)
     do m=1,n_stack
        call fTRANSP_DO(h_C(m,:,:),fb_coll(m,:,:))
     enddo
     call fTRANSP_CLEANUP
     !
     ! NOW:
     !
     ! f_coll  = h  - z*alpha*vp*Ab
     ! fb_coll = hb - z*alpha*vp*Ab 
     !------------------------------------------------

     call proc_time(CPU_Ct_in)

     !-------------------------------------------------------- 
     ! Step 4:
     !
     ! Correct for small amount of particle nonconservation.
     ! This avoids a weak long-wavelength instability for n=0.
     !
     if (is < indx_e .and. electron_method /= 3) then
        ! Total conservation
        fb_coll(:,:,:) = fb_coll(:,:,:)-f_coll(:,:,:)
        call gyro_conserve_all
        h(:,:,:,is) =  h(:,:,:,is)+fb_coll(:,:,:)
     endif
     if (is == indx_e) then 
        ! Number conservation
        h(:,:,:,is) = h(:,:,:,is)+fb_coll(:,:,:)-f_coll(:,:,:)
        call gyro_conserve_number
     endif
     !
     !--------------------------------------------------------
     call proc_time(CPU_Ct_out)
     CPU_temp = CPU_temp+(CPU_Ct_out-CPU_Ct_in)

     !---------------------
     ! END Collision step
     !---------------------

     CPU_Ct = CPU_Ct+CPU_temp
     CPU_Ct_in = CPU_Ct_out-CPU_temp

  enddo

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_collision done]'
  endif

end subroutine gyro_collision
