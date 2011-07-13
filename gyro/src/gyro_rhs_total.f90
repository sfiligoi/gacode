!-----------------------------------------------------------
! gyro_rhs_total.f90
!
! PURPOSE:
!  Compute the RHS of the electron and ion GKEs for both 
!  periodic and nonperiodic boundary conditions.  This is 
!  a dramatically simplified version of the original code 
!  to compute the righthand side.
!-----------------------------------------------------------

subroutine gyro_rhs_total

  use gyro_globals
  use gyro_pointers
  use math_constants

  !-----------------------------------------------------------------------------
  implicit none
  !
  complex, dimension(:,:,:), allocatable :: z_der   
  complex, dimension(:,:,:), allocatable :: z_dis
  complex, dimension(:,:,:,:), allocatable :: cap_h
  complex, dimension(:,:,:,:), allocatable :: lit_h
  !-----------------------------------------------------------------------------

  allocate(z_der(n_stack,n_nek_loc_1,n_kinetic))
  allocate(z_dis(n_stack,n_nek_loc_1,n_kinetic))
  allocate(cap_h(n_stack,i1_buffer:i2_buffer,n_nek_loc_1,n_kinetic))
  allocate(lit_h(n_stack,i1_buffer:i2_buffer,n_nek_loc_1,n_kinetic))

  rhs(:,:,:,:)    = (0.0,0.0)
  rhs_dr(:,:,:,:) = 0.0

  !---------------------------------------------
  if (n_substep == 0) then
     call gyro_rhs_nonlinear
  endif
  !---------------------------------------------

  call proc_time(CPU_rhs_in)

  !----------------------------------------------------------------------
  ! Compute orbit-time derivative for ions ONLY, because this term is 
  ! explicit for ions.
  !
  if (kill_i_parallel_flag == 0) call gyro_tau_derivative
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! Derivative-friendly functions which have array elements outside
  ! the region 1:n_x.
  !
  cap_h(:,:,:,:) = (0.0,0.0)
  lit_h(:,:,:,:) = (0.0,0.0)

  do is=1,n_kinetic
     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   

        ck = class(k)

        do i=1,n_x
           lit_h(:,i,p_nek_loc,is) = h(:,i,p_nek_loc,is)
           cap_h(:,i,p_nek_loc,is) = h(:,i,p_nek_loc,is)+&
                z(is)*alpha_s(is,i)*gyro_u(:,i,p_nek_loc,is)
        enddo ! i

     enddo ! p_nek
  enddo ! is
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! Drift (finite-kx) plus upwind dissipation
  !
  do i=1,n_x

     z_der(:,:,:) = 0.0
     z_dis(:,:,:) = 0.0

     do i_diff=-m_dx,m_dx-i_dx

        z_der(:,:,:) = z_der(:,:,:)+omega_dr(:,i,:,:)*&
             w_d1(i_diff)*cap_h(:,i_loop(i+i_diff),:,:)

        z_dis(:,:,:) = z_dis(:,:,:)+abs(omega_dr(:,i,:,:))*&
             s_d1(i_diff)*lit_h(:,i_loop(i+i_diff),:,:)

     enddo

     ! Derivative plus dissipation
     rhs(:,i,:,:) = rhs(:,i,:,:)+z_der(:,:,:)+z_dis(:,:,:)

     ! Entropy tracking
     rhs_dr(:,i,:,:) = rhs_dr(:,i,:,:)+real(z_dis(:,:,:)*conjg(cap_h(:,i,:,:)))

  enddo ! i
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! Diagonal (ky) drift
  !
  rhs(:,:,:,:) = rhs(:,:,:,:)+i_c*omega_d1(:,:,:,:)*cap_h(:,1:n_x,:,:)
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! Diamagnetic drift
  !
  rhs(:,:,:,:) = rhs(:,:,:,:)-i_c*omega_star(:,:,:,:)*gyro_u(:,:,:,:)
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! Buffer damping
  !
  if (boundary_method == 2) then
     do i=1,n_x
        do is=1,n_kinetic
           rhs(:,i,:,is) = rhs(:,i,:,is)-explicit_damp_vec(is,i)*cap_h(:,i,:,is)
        enddo
     enddo
  endif
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! Adaptive source
  !
  if (source_flag == 1) then

     call gyro_adaptive_source

     do is=1,n_kinetic

        p_nek_loc = 0

        do p_nek=1+i_proc_1,n_nek_1,n_proc_1

           p_nek_loc = p_nek_loc+1

           ie = nek_e(p_nek)  

           do i=1,n_x

              ! In this expression, nu_source = 0 if n > 0.
              ! (see gyro_radial_operators).
              rhs(:,i,p_nek_loc,is) = rhs(:,i,p_nek_loc,is) &
                   -nu_source*h0_eq(is,ie,i)

           enddo ! i

        enddo ! p_nek

     enddo ! is

  endif
  !----------------------------------------------------------------------

  !---------------------------------------------------------
  ! Er shear
  !
  do i=1,n_x
     rhs(:,i,:,:) = rhs(:,i,:,:)-i_c*omega_eb_s(i)*h(:,i,:,:)
  enddo ! i
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Krook ion-ion collision operator:
  !
  if (krook_flag == 1) then
     call gyro_collision_krook
     rhs(:,:,:,1) = rhs(:,:,:,1)+rhs_krook(:,:,:)
  endif
  !---------------------------------------------------------

  call proc_time(CPU_rhs_out)
  CPU_RHS = CPU_RHS + (CPU_rhs_out - CPU_rhs_in)

  deallocate(z_der)
  deallocate(z_dis)
  deallocate(cap_h)
  deallocate(lit_h)

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'*[get_rhs_total done]'
  endif

end subroutine gyro_rhs_total
