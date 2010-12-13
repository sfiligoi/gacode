!---------------------------------------------------
! gyro_collision_ebelli.f90 
!
! PURPOSE:
!  Collision step for real cross-species coupled collisions.
!---------------------------------------------------

subroutine gyro_collision_ebelli

  use gyro_globals
  use gyro_pointers
  use gyro_collision_private

  !----------------------------------------
  implicit none
  !
  integer :: n_i
  integer :: n_j
  integer :: n_k
  integer :: n_d
  integer, external :: parallel_dim
  integer :: test_flag = 0
  !
  real :: CPU_Ct_temp
  !----------------------------------------

  !---------------------
  ! BEGIN Collision step
  !---------------------
  call proc_time(CPU_Ct_temp)

  if (n_field > 1) then
     ! ebelli: not sure about this yet
     call catch_error('ERROR: Collision_method is not yet implemented with electromagnetic effects')
  endif

  !-----------------------------------------------------
  !
  ! Calculate array f_coll to be passed to do_collision:
  !

  ! EAB test
  if(test_flag == 1) then
     h(:,:,:,1) =  v_para(:,:,:,1)
     h(:,:,:,2) =  v_para(:,:,:,2)
  endif

  do is=1,n_kinetic

     f_coll(:,:,:) =  h(:,:,:,is)
     
     ! Permute p_nek,i -> p_ine,k
     n_i = n_n_1*n_energy
     n_j = n_lambda
     n_k = n_x
     n_d = parallel_dim(n_k*n_i,n_proc_1)
     
     call rTRANSP_INIT(n_i,n_j,n_k,NEW_COMM_1)
     do m=1,n_stack
        call rTRANSP_DO(f_coll(m,:,:),h_C_all(is,m,:,:))
     enddo
     call rTRANSP_CLEANUP

  enddo
     
  !
  ! Full collision step (pitch-angle scattering 
  ! plus specified restoring terms).
  !
  call proc_time(CPU_Ct_in)
  call gyro_collision_kernel_ebelli
  call proc_time(CPU_Ct_out)

  do is=1,n_kinetic
     ! Permute p_ine,k -> p_nek,i
     n_i = n_x
     n_j = n_n_1*n_energy
     n_k = n_lambda
     n_d = parallel_dim(n_j*n_k,n_proc_1)
     
     call fTRANSP_INIT(n_i,n_j,n_k,NEW_COMM_1)
     do m=1,n_stack
        call fTRANSP_DO(h_C_all(is,m,:,:),fb_coll(m,:,:))
     enddo
     call fTRANSP_CLEANUP
    
     !--------------------------------------------------------
     ! Correct for small amount of particle nonconservation.
     ! This avoids a weak long-wavelength instability for n=0.
     !
     ! Number conservation for all species
     h(:,:,:,is) = fb_coll(:,:,:)
     call proc_time(CPU_Ct_in)
     call gyro_conserve_number_ebelli(is)
     call proc_time(CPU_Ct_out)
     
     !--------------------------------------------------------
     
     !---------------------
     ! END Collision step
     !---------------------
     
  enddo
  
  call proc_time(CPU_Ct_out)
  
  ! EAB test
  if(test_flag ==1 .and. i_proc == 0) then
     open(unit=5,file='test.out',status='replace')
     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1
        p_nek_loc = p_nek_loc+1
        ie = nek_e(p_nek)
        k  = nek_k(p_nek)
        ck = class(k)
        do i=1,n_x
           do m=1,n_stack
              m0 = m_phys(ck,m)
              if(i==1 .and. ie == 1) then
                 write(5,*) k, m, &
                      (real(h(m,i,p_nek_loc,1))/v_para(m,i,p_nek_loc,1)-1)/dt, &
                      (real(h(m,i,p_nek_loc,2))/v_para(m,i,p_nek_loc,2)-1)/dt, &
                      nu_coll_d(i,ie,1,2)*(-1+mu(2)**2/mu(1)**2), &
                      nu_coll_d(i,ie,2,1)*(-1+mu(1)**2/mu(2)**2)
              endif
           enddo
        enddo
     enddo
     close(5)
     call catch_error('stop')
  end if

  !if(i_proc == 0) then
  !   print *, 'time_total = ', CPU_Ct_out - CPU_Ct_temp
  !endif
  
  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_collision done]'
  endif
  
end subroutine gyro_collision_ebelli
