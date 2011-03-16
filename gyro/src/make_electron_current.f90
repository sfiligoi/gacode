!-----------------------------------------------------
! make_electron_current.f90 [caller BigScience]
!
! PURPOSE:
!  Generate matrix of velocity-space averaged 
!  blending projection integrals.
!---------------------------------------------

subroutine make_electron_current(i_print)

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: i_print
  !
  complex, dimension(n_blend,n_x) :: vel_sum_loc
  complex, dimension(n_blend,n_x) :: vel_sum_glob
  !---------------------------------------------------


  if (electron_method > 2) then     
     ! Separate treatment of GK electrons not required
     coll_vel(:,:,:) = (0.0,0.0)
     return
  endif

  !----------------------------------------------------
  ! Now, compute matrix of blending projections:
  !
  do jp=1,n_blend

     vel_sum_loc = (0.0,0.0)

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
              !  alpha_e FV[(F*_j) ve^2 (F_jp)]
              ! 
              ! because coll_vel is always added to aa_mm,
              ! which already includes the ion current. 

              vel_sum_loc(:,i) = vel_sum_loc(:,i)+ &
                   v_para(m,i,p_nek_loc,n_spec)**2*alpha_s(n_spec,i)* &
                   c_blend(jp,m0,i,p_nek_loc)* &
                   cs_blend(:,m0,i,p_nek_loc)

           enddo ! m

        enddo ! i

     enddo ! p_nek_loc
     !----------------------------------------------------

     call MPI_ALLREDUCE(vel_sum_loc,&
          vel_sum_glob,&
          n_x*n_blend,&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     do i=1,n_x
        coll_vel(i,:,jp) = vel_sum_glob(:,i)
     enddo

  enddo ! jp

  if (i_print == 1) then

     d_theta = pi_2/n_blend

     do i=1,n_x
        print *,i
        do j=1,n_blend
           print '(32(f10.6,1x))',real(coll_vel(i,j,:)) &
                /(mu(2)**2+1.0)/d_theta
        enddo
        do j=1,n_blend
           print '(32(f10.6,1x))',aimag(coll_vel(i,j,:)) &
                /(mu(2)**2+1.0)/d_theta
        enddo
     enddo
  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_electron_current done]'
  endif

end subroutine make_electron_current
