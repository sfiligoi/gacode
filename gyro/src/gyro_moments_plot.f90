!------------------------------------------------------------------------
! gyro_moments_plot.f90
!
! PURPOSE:
!  Make density, energy and V_parallel fluctuation moments (with 
!  same functional form as phi_plot) for diagnostic output file.
! 
!  This routine assumes phi_plot has been computed get_field_plot.
! 
!     density moment: moments_plot(:,:,:,1)
!      energy moment: moments_plot(:,:,:,2)
!  V_parallel moment: moments_plot(:,:,:,3)
!
! NOTES:
!  See http://fusion.gat.com/theory/Gyrousermanual for precise
!  definitions.
!--------------------------------------------------------------------

subroutine gyro_moments_plot

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  integer :: j_plot
  !
  complex, dimension(n_blend,n_x,3) :: vel_sum_loc
  complex, dimension(n_blend,n_x,3) :: vel_sum
  !
  complex, dimension(n_stack,i1_buffer:i2_buffer) :: cap_h
  complex, dimension(n_theta_plot,n_x,3) :: mom_tmp
  !---------------------------------------------------


  if (alltime_index == 0) then
     moments_plot(:,:,:,:) = (0.0,0.0)
  endif
  moments_zero_plot(:,:,:) = 0.0

  do is=1,n_kinetic

     ! First compute gyro_h = <h+z*alpha*<U>>
     !
     ! Note that gyro_h is really a temporary variable 
     ! here; only gyro_h(:,:,:,1) is used.

     gyro_h(:,:,:,1) = (0.0,0.0)

     p_nek_loc = 0

     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)
        k  = nek_k(p_nek)

        ck = class(k)

        cap_h(:,:) = (0.0,0.0)
        do i=1,n_x
           cap_h(:,i) = h(:,i,p_nek_loc,is)+&
                z(is)*alpha_s(is,i)*gyro_u(:,i,p_nek_loc,is)
        enddo

        do m=1,n_stack

           m0 = m_phys(ck,m)

           if (is <= n_gk) then
              do i=1,n_x
                 do i_diff=-m_gyro,m_gyro-i_gyro

                    ip = i+i_diff

                    gyro_h(m,i,p_nek_loc,1) = gyro_h(m,i,p_nek_loc,1)+ &
                         w_gyro(m0,i_diff,i,p_nek_loc,is)*cap_h(m,i_loop(ip))

                 enddo ! i_diff
              enddo ! i
           else
              do i=1,n_x
                 gyro_h(m,i,p_nek_loc,1) = cap_h(m,i)
              enddo ! i
           endif

        enddo ! m

     enddo ! p_nek

     do i=1,n_x
        gyro_h(:,i,:,1) = gyro_h(:,i,:,1) &
             -z(is)*alpha_s(is,i)*field_tau(:,i,:,1)
     enddo

     !----------------------------------------------------
     ! Now, compute blending projections:
     !
     ! vel_sum(j,i) -> FV[ (F*_j) gyro_h ]
     !
     vel_sum_loc = (0.0,0.0)
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

              ! 1: density moment
              vel_sum_loc(:,i,1) = vel_sum_loc(:,i,1)+&
                   gyro_h(m,i,p_nek_loc,1)*cs_blend(:,m0,i,p_nek_loc)

              ! 2: energy moment
              vel_sum_loc(:,i,2) = vel_sum_loc(:,i,2)+&
                   energy(ie,is)*tem_s(is,i)*&
                   gyro_h(m,i,p_nek_loc,1)*cs_blend(:,m0,i,p_nek_loc)

              ! 3: v_parallel moment
              vel_sum_loc(:,i,3) = vel_sum_loc(:,i,3)+&
                   v_para(m,i,p_nek_loc,is)*&
                   gyro_h(m,i,p_nek_loc,1)*cs_blend(:,m0,i,p_nek_loc)

           enddo ! m
        enddo ! i

     enddo
     !--------------------------------------------------------------

     call MPI_ALLREDUCE(vel_sum_loc,&
          vel_sum,&
          size(vel_sum),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     !---------------------------------------------------
     ! Solve for blending coefficients:
     !
     do i=1,n_x

        call ZGETRS('N',&
             n_blend,&
             3,&
             ff_mm(:,:,i,1),&
             n_blend,&
             ff_mm_piv(:,i,1),&
             vel_sum(:,i,:),&
             n_blend,&
             info)

     enddo ! i
     !---------------------------------------------------

     !---------------------------------------------------------------
     ! Theta-interpolate moments_plot from vel_sum projections: 
     !
     do ix=1,3
        do i=1,n_x
           do j_plot=1,n_theta_plot
              mom_tmp(j_plot,i,ix) = sum(vel_sum(:,i,ix)*blend_plot(:,j_plot,i))
           enddo ! j_plot
        enddo ! i
     enddo ! ix
     !----------------------------------------------------------------

     !---------------------------------------------------------------
     ! Integrate to obtain flux-surface averages
     !
     do ix=1,2
        do i=1,n_x
           do j=1,n_blend
              moments_zero_plot(i,is,ix) = &
                   moments_zero_plot(i,is,ix)+c_fluxave(j,i)*real(vel_sum(j,i,ix))
           enddo ! j
        enddo ! i
     enddo ! ix
     !----------------------------------------------------------------

     !----------------------------------------------------------------
     ! Compute moments_plot with averaging:
     !
     moments_plot(:,:,is,:) = moments_plot(:,:,is,:)+&
          w_time(alltime_index+1)*mom_tmp(:,:,:)
     !----------------------------------------------------------------

  enddo ! is

end subroutine gyro_moments_plot
