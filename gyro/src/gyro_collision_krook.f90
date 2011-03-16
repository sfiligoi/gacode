!------------------------------------------------
! gyro_collision_krook.f90
!
! PURPOSE:
!  Evaluate the Krook ion collision operator.
!------------------------------------------------

subroutine gyro_collision_krook

  use mpi
  use gyro_globals
  use gyro_pointers

  !---------------------------------------------------
  implicit none
  !
  complex, dimension(n_blend,n_x) :: sum_loc
  complex, dimension(n_blend,n_x) :: sum_glob
  complex, dimension(2*n_blend,n_x) :: sum2_loc
  complex, dimension(2*n_blend,n_x) :: sum2_glob
  complex, dimension(n_stack,n_x) :: h_tmp
  !---------------------------------------------------


  !----------------------------------------------------
  ! Now, compute blending projections:
  !
  sum_loc(:,:)  = (0.0,0.0)
  sum2_loc(:,:) = (0.0,0.0)
  !
  p_nek_loc = 0
  !
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     h_tmp(:,:) = h(:,:,p_nek_loc,1)

!DIR$ PREFERVECTOR
     do i=1,n_x

        do m=1,n_stack

           m0 = m_phys(ck,m)

           do j=1,n_blend

              ! Poisson Equation (RHS)
              !
              ! FV[(F*_j) he]

              ! Moment: 1

              sum2_loc(j,i) = sum2_loc(j,i)+ &
                   h_tmp(m,i)*&
                   cs_blend(j,m0,i,p_nek_loc)

              ! Moment: v_parallel

              sum_loc(j,i) = sum_loc(j,i)+ &
                   h_tmp(m,i)*v_para(m,i,p_nek_loc,1)*&
                   cs_blend(j,m0,i,p_nek_loc)

              ! Moment: energy

              sum2_loc(j+n_blend,i) = sum2_loc(j+n_blend,i)+ &
                   h_tmp(m,i)*energy(ie,1)*&
                   cs_blend(j,m0,i,p_nek_loc)

           enddo ! j

        enddo ! m

     enddo ! i

  enddo ! p_nek_loc
  !--------------------------------------------------------------

  call MPI_ALLREDUCE(sum_loc,&
       sum_glob,&
       size(sum_glob),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  call MPI_ALLREDUCE(sum2_loc,&
       sum2_glob,&
       size(sum2_glob),&
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
          ff_mm(:,:,i,2),&
          n_blend,&
          ff_mm_piv(:,i,2),&
          sum_glob(:,i),&
          n_blend,&
          info)

     call ZGETRS('N',&
          2*n_blend,&
          1,&
          ff2_mm(:,:,i),&
          2*n_blend,&
          ff2_mm_piv(:,i),&
          sum2_glob(:,i),&
          2*n_blend,&
          info)

  enddo ! i
  !---------------------------------------------------


  !---------------------------------------------------
  ! Get collision term:
  !
  RHS_krook(:,:,:) = (0.0,0.0)
  !
  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

!DIR$ PREFERVECTOR
     do i=1,n_x

        do m=1,n_stack

           m0 = m_phys(ck,m)

           do j=1,n_blend

              RHS_krook(m,i,p_nek_loc) = &
                   RHS_krook(m,i,p_nek_loc)+ &
                   (sum2_glob(j,i)+ &
                   sum_glob(j,i)*v_para(m,i,p_nek_loc,1)+ &
                   sum2_glob(j+n_blend,i)*energy(ie,1))*c_blend(j,m0,i,p_nek_loc)

           enddo ! j

        enddo ! m

     enddo ! i

  enddo ! p_nek_loc

  do i=1,n_x
     RHS_krook(:,i,:) = -nu_i_krook*(h(:,i,:,1)-RHS_krook(:,i,:))
  enddo

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[do_collision_krook done]'
  endif

end subroutine gyro_collision_krook
