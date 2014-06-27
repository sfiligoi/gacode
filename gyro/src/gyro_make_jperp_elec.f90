!------------------------------------------------------------------
! gyro_make_jperp_elec.f90 
!
! PURPOSE:
!  Generate matrix of velocity-space averaged blending projection 
!  integrals for jperp.
!------------------------------------------------------------------

subroutine gyro_make_jperp_elec

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  complex, dimension(n_blend,n_x) :: vel_sum_loc1
  complex, dimension(n_blend,n_x) :: vel_sum_loc2
  complex, dimension(n_blend,n_x) :: vel_sum_glob
  !---------------------------------------------------


  if (electron_method > 2) then     
     ! Separate treatment of GK electrons not required
     coll_vel_perp1(:,:,:) = (0.0,0.0)
     coll_vel_perp2(:,:,:) = (0.0,0.0)
     return
  endif

  !----------------------------------------------------
  ! Now, compute matrix of blending projections:
  !
  do jp=1,n_blend

     vel_sum_loc1 = (0.0,0.0)
     vel_sum_loc2 = (0.0,0.0)

     p_nek_loc = 0

     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)
        k  = nek_k(p_nek)

        ck = class(k)

        do i=1,n_x

           do m=1,n_stack

              m0 = m_phys(ck,m)

              ! ** Includes electron current only: 
              !
              !  FV[ (F*_j) (2*n_i*T_i*V[(ene*lambda*G_perp)^2] (F_jp) ]
              ! 
              ! because coll_vel is always added to aa_mm,
              ! which already includes the ion current. 

              vel_sum_loc1(:,i) = vel_sum_loc1(:,i) +&
                   2.0*den_s(n_spec,i)*tem_s(n_spec,i) * 0.5**2 &
                   * energy(ie)**2 * lambda(i,k)**2 &
                   * cs_blend(:,m0,i,p_nek_loc) &
                   * c_blend(jp,m0,i,p_nek_loc)

              !
              ! FV[ (F*_j) (z_i*n_i*V[(ene*lambda*G*G_perp)]) (F_jp) ]
              !

              vel_sum_loc2(:,i) = vel_sum_loc2(:,i)+&
                   z(n_spec)*den_s(n_spec,i) * 0.5 &
                   * energy(ie) * lambda(i,k) &
                   * cs_blend(:,m0,i,p_nek_loc) &
                   * c_blend(jp,m0,i,p_nek_loc)

           enddo ! m

        enddo ! i

     enddo ! p_nek_loc
     !----------------------------------------------------

     call MPI_ALLREDUCE(vel_sum_loc1,&
          vel_sum_glob,&
          n_x*n_blend,&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     do i=1,n_x
        coll_vel_perp1(i,:,jp) = vel_sum_glob(:,i)
     enddo

     call MPI_ALLREDUCE(vel_sum_loc2,&
          vel_sum_glob,&
          n_x*n_blend,&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     do i=1,n_x
        coll_vel_perp2(i,:,jp) = vel_sum_glob(:,i)
     enddo

  enddo ! jp

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_make_jperp_elec done]'
  endif

end subroutine gyro_make_jperp_elec
