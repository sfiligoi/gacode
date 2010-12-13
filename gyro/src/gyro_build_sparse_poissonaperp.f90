!------------------------------------------------------------
! gyro_build_sparse_poissonaperp.f90 [caller: sparse_solve_*]
!
! PURPOSE:
!  Fill elements of sparse Poisson-aperp submatrix.
!------------------------------------------------------------

subroutine gyro_build_sparse_poissonaperp

  use gyro_globals
  use math_constants
  use poissonaperp_private

  !------------------------
  implicit none
  !------------------------

  ! Zero only the in_1 components:
  m_poissonaperp(:)    = (0.0,0.0)
  indx_poissonaperp(:) = 0

  k_counter = 0

  !-----------------------------------------------------------
  ! Compute block M(1,1) 
  !
  do i=1,n_x_max
     do i_diff=-m_gyro,m_gyro-i_gyro

        ip = i_loop(i+i_diff)

        if (ip >= 1 .and. ip <= n_x) then

           do j=1,n_blend
              do jp=1,n_blend

                 val = ap_mm(i,i_diff,j,jp)

                 ij  = i + (j-1)*n_x
                 ijp = ip + (jp-1)*n_x

                 k_counter = k_counter+1

                 indx_poissonaperp(k_counter)           = ij
                 indx_poissonaperp(k_counter+n_poissonaperp) = ijp

                 m_poissonaperp(k_counter) = val

              enddo ! jp
           enddo ! j

        endif

     enddo ! ip
  enddo ! i

  !-----------------------------------------------------------
  ! Compute block M(2,2)
  !
  do i=1,n_x_max
     do i_diff=-m_gyro,m_gyro-i_gyro

        ip = i_loop(i+i_diff)

        if (ip >= 1 .and. ip <= n_x) then

           do j=1,n_blend
              do jp=1,n_blend

                 val = ab_mm(i,i_diff,j,jp) &
                      + coll_vel_perp1(i,j,jp)*w_g0(i_diff)

                 ij  = i + (j-1)*n_x   + n_x*n_blend
                 ijp = ip + (jp-1)*n_x + n_x*n_blend

                 k_counter = k_counter+1

                 indx_poissonaperp(k_counter)           = ij
                 indx_poissonaperp(k_counter+n_poissonaperp) = ijp

                 m_poissonaperp(k_counter) = val

              enddo ! jp
           enddo ! j

        endif

     enddo ! ip
  enddo ! i

  !-----------------------------------------------------------
  ! Compute block M(1,2)
  !
  do i=1,n_x_max
     do i_diff=-m_gyro,m_gyro-i_gyro

        ip = i_loop(i+i_diff)

        if (ip >= 1 .and. ip <= n_x) then

           do j=1,n_blend
              do jp=1,n_blend

                 val = -2.0*abp_mm(i,i_diff,j,jp) &
                      -2.0*coll_vel_perp2(i,j,jp)*w_g0(i_diff)

                 ij  = i + (j-1)*n_x
                 ijp = ip + (jp-1)*n_x + n_x*n_blend

                 k_counter = k_counter+1

                 indx_poissonaperp(k_counter)           = ij
                 indx_poissonaperp(k_counter+n_poissonaperp) = ijp

                 m_poissonaperp(k_counter) = val

              enddo ! jp
           enddo ! j

        endif

     enddo ! ip
  enddo ! i

  !-----------------------------------------------------------
  ! Compute block M(2,1)
  !
  do i=1,n_x_max
     do i_diff=-m_gyro,m_gyro-i_gyro

        ip = i_loop(i+i_diff)

        if (ip >= 1 .and. ip <= n_x) then

           do j=1,n_blend
              do jp=1,n_blend

                 val = abp_mm(i,i_diff,j,jp)  &
                      + coll_vel_perp2(i,j,jp)*w_g0(i_diff)

                 ij  = i + (j-1)*n_x  + n_x*n_blend
                 ijp = ip + (jp-1)*n_x

                 k_counter = k_counter+1

                 indx_poissonaperp(k_counter)           = ij
                 indx_poissonaperp(k_counter+n_poissonaperp) = ijp

                 m_poissonaperp(k_counter) = val

              enddo ! jp
           enddo ! j

        endif

     enddo ! ip
  enddo ! i

  if (boundary_method == 1 .and. n_1(in_1) == 0) then

     !------------------------------------------------
     ! Here, we need to write a line of 1's in last 
     ! row (i=n_x) of blocks M(1,1), M(2,2).  A 
     ! line of zeros in the last rows of M(x,y) for y!=x
     ! will of course be automatic.
     !     
     ! We do this to ensure that all n=0 modes have 
     ! no radial average.
     !------------------------------------------------

     !------------------------------------------------
     ! Line of 1's in M(1,1)
     ! 
     i = n_x 
     do j=1,n_blend
        do ip=1,n_x

           ! Diagonal in (j,jp)
           jp = j

           ij  = i + (j-1)*n_x
           ijp = ip + (jp-1)*n_x

           k_counter = k_counter+1

           indx_poissonaperp(k_counter)                = ij
           indx_poissonaperp(k_counter+n_poissonaperp) = ijp

           m_poissonaperp(k_counter) = (1.0,0.0)

        enddo ! ip
     enddo ! j

     !------------------------------------------------
     ! Line of 1's in M(2,2)
     ! 
     i = n_x
     do j=1,n_blend
        do ip=1,n_x

           ! Diagonal in (j,jp)
           jp = j

           ij  = i + (j-1)*n_x + n_x*n_blend
           ijp = ip + (jp-1)*n_x + n_x*n_blend

           k_counter = k_counter+1

           indx_poissonaperp(k_counter)                = ij
           indx_poissonaperp(k_counter+n_poissonaperp) = ijp

           m_poissonaperp(k_counter) = (1.0,0.0)

        enddo ! ip
     enddo ! j

  endif

  if (k_counter /= n_poissonaperp) then
     print *,k_counter,n_poissonaperp
     call catch_error("Element count mismatch in gyro_build_sparse_poissonaperp")
  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_build_sparse_poissonaperp done]'
  endif

end subroutine gyro_build_sparse_poissonaperp

