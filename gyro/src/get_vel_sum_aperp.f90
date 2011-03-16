!-----------------------------------------------------
! get_vel_sum_aperp.f90
!
! PURPOSE:
!  Evaluation kernel for velocity-space perpendicular current integral.
!  Takes gyro_h_aperp and generates vel_sum_aperp.
!------------------------------------------------------------

subroutine get_vel_sum_aperp

  use mpi
  use gyro_globals
  use gyro_pointers

  !---------------------------------------------------
  implicit none
  !
  complex, dimension(n_blend,n_x) :: vel_sum_aperp_loc
  complex, dimension(n_stack,n_x) :: gz
  !---------------------------------------------------


  !----------------------------------------------------
  ! Now, compute blending projections:
  !
  vel_sum_aperp_loc(:,:) = (0.0,0.0)
  !
  p_nek_loc = 0
  !
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     gz(:,:) = (0.0,0.0)

     ! Moved energy(ie,is) coefficient to gz(:,i) loop
     ! from vel_sum_aperp_loc loop to accomodate 
     ! species variation.
     ! EMB

     do i=1,n_x
        do is=1,n_kinetic 
           gz(:,i) = gz(:,i)-gyro_h_aperp(:,i,p_nek_loc,is)*tem_s(is,i)*&
                             energy(ie,is)  ! EMB
        enddo
     enddo

     do i=1,n_x

        do m=1,n_stack

           m0 = m_phys(ck,m)

           ! Ampere Equation - Perpendicular (RHS)
           !
           ! sum_s FV[(F*_j) G_perp(hi)*(-T_s*ene*lambda)]

           do j=1,n_blend
              vel_sum_aperp_loc(j,i) = vel_sum_aperp_loc(j,i)+&
                   lambda(i,k) *&
                   gz(m,i)*cs_blend(j,m0,i,p_nek_loc) 
           enddo

        enddo ! m

     enddo ! i

  enddo ! p_nek_loc
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! Final velocity-space sum on subgroup
  !
  call MPI_ALLREDUCE(vel_sum_aperp_loc,&
       vel_sum_aperp,&
       n_x*n_blend,&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
  !--------------------------------------------------------------

  if (n_1(in_1) == 0 .and. boundary_method == 1) then

     ! Eliminate average in periodic case

     vel_sum_aperp(:,n_x) = (0.0,0.0)

  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[vel_sum_aperp done]'
  endif

end subroutine get_vel_sum_aperp
