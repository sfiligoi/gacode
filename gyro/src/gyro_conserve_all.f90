!------------------------------------------------
! gyro_conserve_all.f90
!
! PURPOSE:
!  Correct a collision operator for number,
!  momentum and energy conservation.
!------------------------------------------------

subroutine gyro_conserve_all

  use gyro_globals
  use gyro_pointers

  !---------------------------------------------------
  implicit none
  !
  complex, dimension(n_blend,n_x) :: sum_loc
  complex, dimension(n_blend,n_x) :: sum_glob
  complex, dimension(2*n_blend,n_x) :: sum2_loc
  complex, dimension(2*n_blend,n_x) :: sum2_glob
  !
  complex, dimension(n_stack,n_x,n_nek_loc_1) :: cf_new
  !---------------------------------------------------

  include 'mpif.h'
 
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

     do i=1,n_x

        do m=1,n_stack

           m0 = m_phys(ck,m)

           do j=1,n_blend

              ! Poisson Equation (RHS)
              !
              ! FV[(F*_j) he]

              ! Moment: 1

              sum2_loc(j,i) = sum2_loc(j,i)+ &
                   fb_coll(m,i,p_nek_loc)*&
                   cs_blend(j,m0,i,p_nek_loc)

              ! Moment: v_parallel

              sum_loc(j,i) = sum_loc(j,i)+ &
                   fb_coll(m,i,p_nek_loc)*v_para(m,i,p_nek_loc,is)*&
                   cs_blend(j,m0,i,p_nek_loc)

              ! Moment: energy

              sum2_loc(j+n_blend,i) = sum2_loc(j+n_blend,i)+ &
                   fb_coll(m,i,p_nek_loc)*energy(ie,1)*&
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
  cf_new(:,:,:) = fb_coll(:,:,:)
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

           do j=1,n_blend

              cf_new(m,i,p_nek_loc) = &
                   cf_new(m,i,p_nek_loc)- &
                   (sum2_glob(j,i)+ &
                   sum_glob(j,i)*v_para(m,i,p_nek_loc,is)+ &
                   sum2_glob(j+n_blend,i)*energy(ie,1))*c_blend(j,m0,i,p_nek_loc)

           enddo ! j

        enddo ! m

     enddo ! i

  enddo ! p_nek_loc

  fb_coll(:,:,:) = cf_new(:,:,:)

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_conserve_all done]'
  endif

end subroutine gyro_conserve_all
