!-----------------------------------------------------
! get_ampere_coll.f90
!
! PURPOSE:
!  Solve collison-Ampere equation with periodic OR 
!  nonperiodic radial boundary conditions using 
!  exact blending-function representation.
!
! NOTES:
!  The effect of the collisional advance on the 
!  fields is quite simple.  Since the electron 
!  collision operator conserves particle number, 
!  there is no change in phi.  A_parallel gets 
!  a small increment from collisions and this 
!  we evolve explicitly.
!---------------------------------------------

subroutine get_ampere_coll

  use gyro_globals
  use gyro_collision_private
  use gyro_pointers

  !---------------------------------------------------
  implicit none
  !
  complex, dimension(n_blend,n_x) :: c0_a_temp
  complex, dimension(n_blend,n_x) :: c0_a_cdot
  complex, dimension(n_blend,n_x) :: vel_sum_loc
  !
  complex, dimension(n_stack,i1_buffer:i2_buffer) :: f1
  complex, dimension(n_stack,i1_buffer:i2_buffer) :: f2
  !
  complex, external :: BLEND_F
  !---------------------------------------------------

  include 'mpif.h'

  ix = 2

  if (n_field == 3) then
     call catch_error('ERROR: Cannot call get_ampere_coll if N_FIELD=3')
  endif

  !-----------------------------
  gyro_u(:,:,:,:)   = (0.0,0.0)
  !-----------------------------

  !----------------------------------------------------
  ! Now, compute blending projections:
  !
  ! vel_sum(j,i) -> FV[ (F_j) nu ve fe ]
  !                 
  vel_sum_loc(:,:)  = (0.0,0.0)
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

           ! collision-Ampere Equation (RHS)
           !
           ! V2[ nu v_p fe ]
           !
           ! fe = he + alpha*vp*a_parallel

           vel_sum_loc(:,i) = vel_sum_loc(:,i)+&
                cs_blend(:,m0,i,p_nek_loc)*v_para(m,i,p_nek_loc,indx_e)* &
                2.0*nu_total(i,ie,indx_e)*(h(m,i,p_nek_loc,indx_e)+ & 
                alpha_s(indx_e,i)*v_para(m,i,p_nek_loc,indx_e)*&
                field_tau(m,i,p_nek_loc,2))

        enddo ! m

     enddo ! i

  enddo ! p_nek_loc
  !--------------------------------------------------------------

  call MPI_ALLREDUCE(vel_sum_loc,&
       vel_sum_a,&
       n_x*n_blend,&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  ! NOTE:  
  !
  ! At this point, vel_sum = S_A_coll

  if (n_1(in_1) == 0 .and. boundary_method == 1) then

     ! Eliminate average in periodic case

     vel_sum_a(:,n_x) = (0.0,0.0)

  endif

  ! sparse_solve_x (below) will overwrite c0_a, so we need 
  ! to store it in a temporary variable (version 6.0.0)
  c0_a_temp(:,:) = field_blend(:,:,2)

  if (sparse_method == 1) then
     call sparse_solve_umfpack(n_ampere,n_ampere_row,2,1)
  else
     call sparse_solve_mumps(n_ampere,n_ampere_row,2,1)
  endif

  ! Collisional rate of change of c0_a
  c0_a_cdot(:,:)     = field_blend(:,:,2)
  field_blend(:,:,2) = c0_a_temp(:,:)

  ! Do explicit collisional advance of A_parallel:
  field_blend(:,:,2) = field_blend(:,:,2) + dt*c0_a_cdot(:,:)

  !---------------------------------------------------------------
  ! Interpolate a_parallel onto orbit-grid 
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

           field_tau(m,i,p_nek_loc,2) = &
                sum(field_blend(:,i,2)*c_blend(:,m0,i,p_nek_loc))

        enddo ! m

     enddo ! i

  enddo ! p_nek
  !
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! IONS:
  !
  ! gyro_u -> <phi> - v <A>
  !
  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)
     k  = nek_k(p_nek)

     ck = class(k)

     f1(:,:)= (0.0,0.0)
     f2(:,:)= (0.0,0.0)
     do i=1,n_x
        f1(:,i) = field_tau(:,i,p_nek_loc,1)
        f2(:,i) = field_tau(:,i,p_nek_loc,2)
     enddo

     do is=1,n_gk

        do m=1,n_stack

           m0 = m_phys(ck,m)

           do i_diff=-m_gyro,m_gyro-i_gyro
              do i=1,n_x

                 gyro_u(m,i,p_nek_loc,is) = gyro_u(m,i,p_nek_loc,is)+ &
                      w_gyro(m0,i_diff,i,p_nek_loc,is)*(f1(m,i_loop(i+i_diff))- &
                      v_para(m,i,p_nek_loc,is)*f2(m,i_loop(i+i_diff)))

              enddo ! i_diff

           enddo ! i
        enddo ! m

     enddo ! is

  enddo ! p_nek_loc
  !
  ! ELECTRONS:
  ! 
  ! gyro_u -> phi - v A
  !
  gyro_u(:,:,:,indx_e) = field_tau(:,:,:,1)-v_para(:,:,:,indx_e)*field_tau(:,:,:,2)
  !
  !----------------------------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[get_ampere_coll done]'
  endif

end subroutine get_ampere_coll
