!-----------------------------------------------------------------
! gyro_field_interpolation
!
! PURPOSE:
!  Compute Psi_a(R) and Chi_a(R) [Sec 3.7] of Technical Guide.
!
!  This requires first interpolating phi, A_par and B_par onto 
!  fixed theta-grids and tau-grids.  
!
!  Psi_a(R) -> gyro_uv
!
!  gyro_uv(1) = G0a [ phi ]  
!  gyro_uv(2) = G0a [ -v_par A_par ] 
!  gyro_uv(3) = G1a [ -v_perp A_perp ]
!
!  Chi_a(R) -> kyro_uv
!
!  kyro_uv(1) = G2a [ phi ]  
!  kyro_uv(2) = G2a [ -v_par A_par ] 
!  kyro_uv(3) = G3a [ -v_perp A_perp ]
!
!  This is valid for periodic OR nonperiodic boundary
!  conditions. 
!
!  NOTE:
!   This routine is expensive so some coding is inelegant 
!   for the sake of speed.
!-----------------------------------------------------------------

subroutine gyro_field_interpolation

  use gyro_globals
  use gyro_pointers

  !-----------------------------------------------------------------
  implicit none
  !
  real :: x
  !
  complex :: vtemp(n_field)
  complex :: cmplx_phase
  !
  complex, dimension(n_stack,i1_buffer:i2_buffer,n_field) :: vf
  complex, external :: BLEND_F
  !-----------------------------------------------------------------

  call gyro_timer_in('Field-interp.a')
!$acc kernels loop 
!$omp parallel do default(shared) private(cmplx_phase,j_int,x,vtemp,j)
  do i=1,n_x
     cmplx_phase = phase(in_1,i)
     do j_int=1,n_theta_int

        x = -1.0+2.0*(j_int-1)/n_theta_int

        vtemp(:) = (0.0,0.0)

        do j=1,n_blend
           vtemp(:) = vtemp(:)+field_blend(j,i,:)*BLEND_F(j,x,cmplx_phase)
        enddo
        phi(j_int,i,:) = vtemp(:) 

     enddo ! j_int
  enddo ! i
!$end parallel do
!$acc end kernels loop 

  !---------------------------------------------------------------
  ! Interpolate phi, A_par, and B_par onto orbit-grid:
  !
!$acc kernels loop
!$omp parallel do default(shared) private(p_nek_loc,p_nek,k,ck,m,m0,ix)
  do i=1,n_x
     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        k  = nek_k(p_nek)   
        ck = class(k)

        do m=1,n_stack

           m0 = m_phys(ck,m)

           do ix=1,n_field
              field_tau(m,i,p_nek_loc,ix) = &
                   sum(field_blend(:,i,ix)*c_blend(:,m0,i,p_nek_loc))
           enddo ! ix

        enddo ! m

     enddo ! p_nek
  enddo ! i
!$end parallel do
!$acc end kernels loop
  !
  !---------------------------------------------------------------

  call gyro_timer_out('Field-interp.a')
  call gyro_timer_in('Field-interp.b')

  !---------------------------------------------------------------
  ! GK SPECIES:
  !
  ! gyro_u -> <phi> - v_par <A_par> + (m*v_perp^2)/(z*e*B)*G1[B_par]
  !
  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)
     k  = nek_k(p_nek)

     ck = class(k)

     if (boundary_method == 1) then
        do i=1-m_gyro,0
           vf(:,i,:) = field_tau(:,i+n_x,p_nek_loc,:)
        enddo
        do i=n_x+1,n_x+m_gyro-i_gyro
           vf(:,i,:) = field_tau(:,i-n_x,p_nek_loc,:)
        enddo
     else
        do i=1-m_gyro,0
           vf(:,i,:) = 0.0
        enddo
        do i=n_x+1,n_x+m_gyro-i_gyro
           vf(:,i,:) = 0.0
        enddo
     endif
     do ix=1,n_field
        do i=1,n_x
           do m=1,n_stack
              vf(m,i,ix) = field_tau(m,i,p_nek_loc,ix)
           enddo
        enddo
     enddo

     do is=1,n_gk

        select case (n_field) 

        case(1)

           gyro_uv(:,:,p_nek_loc,is,1) = (0.0,0.0)
           kyro_uv(:,:,p_nek_loc,is,1) = (0.0,0.0)

!$acc kernels loop
!$omp parallel do default(shared) private(i_diff,m)
           do i=1,n_x
              do i_diff=-m_gyro,m_gyro-i_gyro
                 do m=1,n_stack

                    ! Psi_a(R) in Technical Guide
                    gyro_uv(m,i,p_nek_loc,is,1) = gyro_uv(m,i,p_nek_loc,is,1)+&
                         w_gyro0(m,i_diff,i,p_nek_loc,is)*vf(m,i+i_diff,1)

                    ! Chi_a(R) in Technical Guide
                    kyro_uv(m,i,p_nek_loc,is,1) = kyro_uv(m,i,p_nek_loc,is,1)+&
                         w_gyro2(m,i_diff,i,p_nek_loc,is)*vf(m,i+i_diff,1)

                 enddo
              enddo
           enddo
!$omp end parallel do
!$acc end kernels loop

        case (2) 

           gyro_uv(:,:,p_nek_loc,is,1:2) = (0.0,0.0)
           kyro_uv(:,:,p_nek_loc,is,1:2) = (0.0,0.0)

!$acc kernels loop
!$omp parallel do default(shared) private(i_diff,m)
           do i=1,n_x
              do i_diff=-m_gyro,m_gyro-i_gyro
                 do m=1,n_stack

                    ! Psi_a(R) in Technical Guide
                    gyro_uv(m,i,p_nek_loc,is,1) = gyro_uv(m,i,p_nek_loc,is,1)+&
                         w_gyro0(m,i_diff,i,p_nek_loc,is)*vf(m,i+i_diff,1)
                    gyro_uv(m,i,p_nek_loc,is,2) = gyro_uv(m,i,p_nek_loc,is,2)+&
                         w_gyro0(m,i_diff,i,p_nek_loc,is)* &
                         (-v_para(m,i,p_nek_loc,is)*vf(m,i+i_diff,2))

                    ! Chi_a(R) in Technical Guide
                    kyro_uv(m,i,p_nek_loc,is,1) = kyro_uv(m,i,p_nek_loc,is,1) + &
                         w_gyro2(m,i_diff,i,p_nek_loc,is)*vf(m,i+i_diff,1)
                    kyro_uv(m,i,p_nek_loc,is,2) = kyro_uv(m,i,p_nek_loc,is,2) + &
                         w_gyro2(m,i_diff,i,p_nek_loc,is)*&
                         (-v_para(m,i,p_nek_loc,is)*vf(m,i+i_diff,2))
                 enddo
              enddo
           enddo
!$omp end parallel do
!$acc end kernels loop

        case (3)

           gyro_uv(:,:,p_nek_loc,is,1:3) = (0.0,0.0)
           kyro_uv(:,:,p_nek_loc,is,1:3) = (0.0,0.0)

!$acc kernels loop
!$omp parallel do default(shared) private(i_diff,m)
           do i=1,n_x
              do i_diff=-m_gyro,m_gyro-i_gyro
                 do m=1,n_stack

                    ! Psi_a(R) in Technical Guide
                    gyro_uv(m,i,p_nek_loc,is,1) = gyro_uv(m,i,p_nek_loc,is,1)+&
                         w_gyro0(m,i_diff,i,p_nek_loc,is)*vf(m,i+i_diff,1)
                    gyro_uv(m,i,p_nek_loc,is,2) = gyro_uv(m,i,p_nek_loc,is,2)+&
                         w_gyro0(m,i_diff,i,p_nek_loc,is)* &
                         (-v_para(m,i,p_nek_loc,is)*vf(m,i+i_diff,2))
                    gyro_uv(m,i,p_nek_loc,is,3) = gyro_uv(m,i,p_nek_loc,is,3)+&
                         (w_gyro1(m,i_diff,i,p_nek_loc,is) &
                         *2.0*energy(ie)*lambda(i,k)*tem_s(is,i)/z(is))*&
                         vf(m,i+i_diff,3)

                    ! Chi_a(R) in Technical Guide
                    kyro_uv(m,i,p_nek_loc,is,1) = kyro_uv(m,i,p_nek_loc,is,1) + &
                         w_gyro2(m,i_diff,i,p_nek_loc,is)*vf(m,i+i_diff,1)
                    kyro_uv(m,i,p_nek_loc,is,2) = kyro_uv(m,i,p_nek_loc,is,2) + &
                         w_gyro2(m,i_diff,i,p_nek_loc,is)*&
                         (-v_para(m,i,p_nek_loc,is)*vf(m,i+i_diff,2))
                    kyro_uv(m,i,p_nek_loc,is,3) = kyro_uv(m,i,p_nek_loc,is,3)+&
                         w_gyro3(m,i_diff,i,p_nek_loc,is) &
                         *2.0*energy(ie)*lambda(i,k)*tem_s(is,i)/z(is)*&
                         vf(m,i+i_diff,3)
                   
                 enddo
              enddo
           enddo
!$omp end parallel do
!$acc end kernels loop

        end select

     enddo ! is

  enddo ! p_nek_loc
  !
  !----------------------------------------------------------------

  if (electron_method == 2) then

     ! Drift-kinetic ELECTRONS:
     ! 
     ! G0e   -> 1
     ! G1e   -> 1/2
     ! Psi_e -> phi - v_par A_par + (1/2) [2*ene*lambda*temp/z*B_par]
     !
     ! G2e   -> 0
     ! G3e   -> 0
     ! Chi_e -> 0

     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1
        p_nek_loc = p_nek_loc+1
        ie = nek_e(p_nek)
        k  = nek_k(p_nek)

        select case (n_field) 

        case(1)

           gyro_uv(:,:,p_nek_loc,n_spec,1) = field_tau(:,:,p_nek_loc,1)

        case(2)

!$acc kernels loop
!$omp parallel do default(shared) private(m)
           do i=1,n_x
              do m=1,n_stack
                 gyro_uv(m,i,p_nek_loc,n_spec,1) = field_tau(m,i,p_nek_loc,1)
                 gyro_uv(m,i,p_nek_loc,n_spec,2) = &
                      -v_para(m,i,p_nek_loc,n_spec)*field_tau(m,i,p_nek_loc,2)
              enddo
           enddo
!$omp end parallel do
!$acc end kernels loop

        case (3)

!$acc kernels loop
!$omp parallel do default(shared) private(m)
           do i=1,n_x
              do m=1,n_stack
                 gyro_uv(m,i,p_nek_loc,n_spec,1) = field_tau(m,i,p_nek_loc,1)
                 gyro_uv(m,i,p_nek_loc,n_spec,2) = &
                      -v_para(m,i,p_nek_loc,n_spec)*field_tau(m,i,p_nek_loc,2)
                 gyro_uv(m,i,p_nek_loc,n_spec,3) = &
                      energy(ie)*lambda(i,k)*tem_s(n_spec,i)/z(n_spec) &
                      *field_tau(m,i,p_nek_loc,3)
              enddo
           enddo
!$omp end parallel do
!$acc end kernels loop
        end select

     enddo ! p_nek

     kyro_uv(:,:,:,n_spec,:) = (0.0,0.0)

  endif


!!! added acc kernels loop by Sri V. that is not mapped to the openMP set.
!!! will have to remove the one of the (!) to get working.
!!$acc kernels loop
  do is=1,n_kinetic
     do p_nek_loc=1,n_nek_loc_1
        do i=1,n_x
           do m=1,n_stack
              gyro_u(m,i,p_nek_loc,is) = sum(gyro_uv(m,i,p_nek_loc,is,:))
           enddo
        enddo
     enddo
  enddo
!!$acc end kernels loop

  call gyro_timer_out('Field-interp.b')

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_field_interpolation done]'
  endif

end subroutine gyro_field_interpolation
