!---------------------------------------------------------
! gyro_field_interpolation
!
! PURPOSE:
!  Interpolate phi, A_par and B_par onto fixed 
!  theta-grids and tau-grids.  Also compute
!
!  gyro_uv(1) = < phi > 
!  gyro_uv(2) = < -v_par A_par> 
!  gyro_uv(3) = < -v_perp A_perp >
!
!  gyro_u = sum(gyro_uv(:))
!
!  This is valid for periodic OR nonperiodic boundary
!  conditions. 
!
!  NOTE:
!   This routine is expensive so some coding is inelegant 
!   for the sake of speed.
!---------------------------------------------------------

subroutine gyro_field_interpolation

  use gyro_globals
  use gyro_pointers

  !-----------------------------------------------------
  implicit none
  !
  real :: x
  !
  complex :: vtemp(n_field)
  complex :: cmplx_phase
  !
  complex, dimension(n_stack,i1_buffer:i2_buffer,n_field) :: vf
  complex, external :: BLEND_F
  !-----------------------------------------------------

  call gyro_timer_in('Field-interp.a')
!$acc region 
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
!$acc end region 

  !---------------------------------------------------------------
  ! Interpolate phi, A_par, and B_par onto orbit-grid:
  !
!$acc region
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
!$acc end regend 
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

!$acc region
!$omp parallel do default(shared) private(i_diff,m)
           do i=1,n_x
              do i_diff=-m_gyro,m_gyro-i_gyro
                 do m=1,n_stack
                    gyro_uv(m,i,p_nek_loc,is,1) = gyro_uv(m,i,p_nek_loc,is,1)+&
                         w_gyro(m,i_diff,i,p_nek_loc,is)*vf(m,i+i_diff,1)
                    kyro_uv(m,i,p_nek_loc,is,1) = kyro_uv(m,i,p_nek_loc,is,1) + &
                         w_gyro_rot(m,i_diff,i,p_nek_loc,is)*&
                         vf(m,i+i_diff,1)
                 enddo
              enddo
           enddo
!$omp end parallel do
!$acc end region

        case (2) 

           gyro_uv(:,:,p_nek_loc,is,:) = (0.0,0.0)
           kyro_uv(:,:,p_nek_loc,is,:) = (0.0,0.0)

!$acc region
!$omp parallel do default(shared) private(i_diff,m)
           do i=1,n_x
              do i_diff=-m_gyro,m_gyro-i_gyro
                 do m=1,n_stack
                    gyro_uv(m,i,p_nek_loc,is,1) = gyro_uv(m,i,p_nek_loc,is,1)+&
                         w_gyro(m,i_diff,i,p_nek_loc,is)*vf(m,i+i_diff,1)
                    gyro_uv(m,i,p_nek_loc,is,2) = gyro_uv(m,i,p_nek_loc,is,2)+&
                         vf(m,i+i_diff,2)*&
                         (-w_gyro(m,i_diff,i,p_nek_loc,is)*v_para(m,i,p_nek_loc,is))
                    kyro_uv(m,i,p_nek_loc,is,1) = kyro_uv(m,i,p_nek_loc,is,1) + &
                         w_gyro_rot(m,i_diff,i,p_nek_loc,is)*vf(m,i+i_diff,1)
                    kyro_uv(m,i,p_nek_loc,is,2) = kyro_uv(m,i,p_nek_loc,is,2) + &
                         w_gyro_rot(m,i_diff,i,p_nek_loc,is)*&
                         (-v_para(m,i,p_nek_loc,is)*vf(m,i+i_diff,2))
                 enddo
              enddo
           enddo
!$omp end parallel do
!$acc end region

        case (3)

           gyro_uv(:,:,p_nek_loc,is,:) = (0.0,0.0)
           kyro_uv(:,:,p_nek_loc,is,:) = (0.0,0.0)

!$acc region
!$omp parallel do default(shared) private(i_diff,m)
           do i=1,n_x
              do i_diff=-m_gyro,m_gyro-i_gyro
                 do m=1,n_stack
                    gyro_uv(m,i,p_nek_loc,is,1) = gyro_uv(m,i,p_nek_loc,is,1)+&
                         w_gyro(m,i_diff,i,p_nek_loc,is)*vf(m,i+i_diff,1)
                    gyro_uv(m,i,p_nek_loc,is,2) = gyro_uv(m,i,p_nek_loc,is,2)+&
                         (-w_gyro(m,i_diff,i,p_nek_loc,is)*&
                         v_para(m,i,p_nek_loc,is)*vf(m,i+i_diff,2))
                    gyro_uv(m,i,p_nek_loc,is,3) = gyro_uv(m,i,p_nek_loc,is,3)+&
                         (w_gyro_aperp(m,i_diff,i,p_nek_loc,is) &
                         *2.0*energy(ie,is)*lambda(i,k)*tem_s(is,i)/z(is))*&
                         vf(m,i+i_diff,3)
                    kyro_uv(m,i,p_nek_loc,is,1) = kyro_uv(m,i,p_nek_loc,is,1) + &
                         w_gyro_rot(m,i_diff,i,p_nek_loc,is)*&
                         vf(m,i+i_diff,1)
                    kyro_uv(m,i,p_nek_loc,is,2) = kyro_uv(m,i,p_nek_loc,is,2) + &
                         w_gyro_rot(m,i_diff,i,p_nek_loc,is)*&
                         (-v_para(m,i,p_nek_loc,is)*vf(m,i+i_diff,2))
                    ! Note missing third component!
                 enddo
              enddo
           enddo
!$omp end parallel do
!$acc end region

        end select

     enddo ! is

  enddo ! p_nek_loc
  !
  !----------------------------------------------------------------

  if (electron_method == 2) then

     ! DK ELECTRONS:
     ! 
     ! gyro_u -> phi - v_par A_par + ene*lambda*temp/z*B_par

     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1
        p_nek_loc = p_nek_loc+1
        ie = nek_e(p_nek)
        k  = nek_k(p_nek)

        select case (n_field) 

        case(1)

           gyro_uv(:,:,p_nek_loc,n_spec,1) = field_tau(:,:,p_nek_loc,1)

        case(2)

!$acc region
!$omp parallel do default(shared) private(m)
           do i=1,n_x
              do m=1,n_stack
                 gyro_uv(m,i,p_nek_loc,n_spec,1) = field_tau(m,i,p_nek_loc,1)
                 gyro_uv(m,i,p_nek_loc,n_spec,2) = &
                      -v_para(m,i,p_nek_loc,n_spec)*field_tau(m,i,p_nek_loc,2)
              enddo
           enddo
!$omp end parallel do
!$acc end region

        case (3)

!$acc region
!$omp parallel do default(shared) private(m)
           do i=1,n_x
              do m=1,n_stack
                 gyro_uv(m,i,p_nek_loc,n_spec,1) = field_tau(m,i,p_nek_loc,1)
                 gyro_uv(m,i,p_nek_loc,n_spec,2) = &
                      -v_para(m,i,p_nek_loc,n_spec)*field_tau(m,i,p_nek_loc,2)
                 gyro_uv(m,i,p_nek_loc,n_spec,3) = &
                      energy(ie,n_spec)*lambda(i,k)*tem_s(n_spec,i)/z(n_spec) &
                      *field_tau(m,i,p_nek_loc,3)
              enddo
           enddo
!$omp end parallel do
!$acc end region
        end select

     enddo ! p_nek

     kyro_uv(:,:,:,n_spec,:) = (0.0,0.0)

  endif


!!! added acc region by Sri V. that is not mapped to the openMP set.
!!! will have to remove the one of the (!) to get working.
!!$acc region
  do is=1,n_kinetic
     do p_nek_loc=1,n_nek_loc_1
        do i=1,n_x
           do m=1,n_stack
              gyro_u(m,i,p_nek_loc,is) = sum(gyro_uv(m,i,p_nek_loc,is,:))
           enddo
        enddo
     enddo
  enddo
!!$acc end region

  call gyro_timer_out('Field-interp.b')

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_field_interpolation done]'
  endif

end subroutine gyro_field_interpolation
