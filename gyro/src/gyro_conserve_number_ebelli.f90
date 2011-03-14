!------------------------------------------------
! gyro_conserve_number_ebelli.f90
!
! PURPOSE:
!  Enforce number conservation by correcting h
!  for species ks 
!  (based on value in h_old) to remove the 
!  component which alters number
!  EBELLI note: This differs from gyro_conserve_number.f90
!  which only corrects the electron distribution.
!------------------------------------------------

subroutine gyro_conserve_number_ebelli(ks)

  use mpi
  use gyro_globals
  use gyro_pointers

  !---------------------------------------------------
  implicit none
  integer, intent(in) :: ks
  complex, dimension(n_blend,n_x) :: vel_sum_loc
  complex, dimension(n_blend,n_x) :: vel_sum
  !---------------------------------------------------

  !----------------------------------------------------
  ! Now, compute blending projections:
  !
  vel_sum_loc(:,:) = (0.0,0.0)
  !
  p_nek_loc = 0
  !
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     do i=1,n_x

        do m=1,n_stack

           m0 = m_phys(ck,m)

           vel_sum_loc(:,i) = vel_sum_loc(:,i)+ &
                (h_old(m,i,p_nek_loc,ks)-h(m,i,p_nek_loc,ks))* &
                cs_blend(:,m0,i,p_nek_loc)

        enddo ! m

     enddo ! i

  enddo ! p_nek_loc
  !--------------------------------------------------------------

  call MPI_ALLREDUCE(vel_sum_loc,&
       vel_sum,&
       n_blend*n_x,&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  !---------------------------------------------------
  ! Solve for blending coefficients
  !
  do i=1,n_x

     call ZGETRS('N',&
          n_blend,&
          1,&
          ff_mm(:,:,i,1),&
          n_blend,&
          ff_mm_piv(:,i,1),&
          vel_sum(:,i),&
          n_blend,&
          info)

  enddo ! i
  !---------------------------------------------------


  !---------------------------------------------------
  ! Correct distribution
  !
  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     do i=1,n_x

        do m=1,n_stack

           m0 = m_phys(ck,m)

           h(m,i,p_nek_loc,ks) = h(m,i,p_nek_loc,ks)+ &
                sum(vel_sum(:,i)*c_blend(:,m0,i,p_nek_loc))

        enddo ! m

     enddo ! i

  enddo ! p_nek_loc

end subroutine gyro_conserve_number_ebelli
