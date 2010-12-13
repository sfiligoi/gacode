!-------------------------------------------------------
! gyro_build_sparse_maxwell.f90 [caller: sparse_solve_*]
!
! PURPOSE:
!  Fill elements of sparse Maxwell matrix.
! 
! NOTES:
!  The Maxwell matrix has the block structure:
!  
!                 / M(1,1)  M(1,2) M(1,3) \
!     m_maxwell = |                       |
!                 | M(2,1)  M(2,2) M(2,3) |    
!                 |                       |
!                 \ M(3,1)  M(3,2) M(3,3) /
!
!                                / phi   \
!  This acts on a field vector = |       |
!                                | A_par |
!                                |       |
!                                \ B_par /
! 
!-------------------------------------------------------

subroutine gyro_build_sparse_maxwell

  use gyro_globals
  use math_constants
  use maxwell_private

  !------------------------
  implicit none
  !------------------------


  ! Zero only the in_1 components:
  m_maxwell(:)    = (0.0,0.0)
  indx_maxwell(:) = 0

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

                 if (i_diff == 0 .and. electron_method == 2) then
                    val = val-imp(i,j,jp,1) 
                 endif

                 ij  = i + (j-1)*n_x
                 ijp = ip + (jp-1)*n_x

                 k_counter = k_counter+1

                 indx_maxwell(k_counter)           = ij
                 indx_maxwell(k_counter+n_maxwell) = ijp

                 m_maxwell(k_counter) = val

              enddo ! j
           enddo ! jp

        endif

     enddo ! i_diff
  enddo ! i

  !--------------------------------------
  ! Continue past here for EM simulations.
  !--------------------------------------

  ! First A_parallel

  if (n_field > 1) then 

     !-----------------------------------------------------------
     ! Compute block M(2,2)
     !
     do i=1,n_x_max
        do i_diff=-mg_dx,mg_dx-ig_dx

           ip = i_loop(i+i_diff)

           if (ip >= 1 .and. ip <= n_x) then

              do j=1,n_blend
                 do jp=1,n_blend

                    val = aa_mm(i,i_diff,j,jp) 

                    if (i_diff == 0) then
                       val = val-imp(i,j,jp,4)
                    endif

                    ij  = i + (j-1)*n_x + n_x*n_blend
                    ijp = ip + (jp-1)*n_x + n_x*n_blend

                    k_counter = k_counter+1

                    indx_maxwell(k_counter)           = ij
                    indx_maxwell(k_counter+n_maxwell) = ijp

                    m_maxwell(k_counter) = val

                 enddo ! jp
              enddo ! j

           endif

        enddo ! ip
     enddo ! i
     !-----------------------------------------------------------


     !-----------------------------------------------------------
     ! Compute block M(1,2)
     !
     do i=1,n_x_max
        do j=1,n_blend
           do jp=1,n_blend

              ij  = i + (j-1)*n_x
              ijp = i + (jp-1)*n_x + n_x*n_blend

              k_counter = k_counter+1

              indx_maxwell(k_counter)           = ij
              indx_maxwell(k_counter+n_maxwell) = ijp

              m_maxwell(k_counter) = -imp(i,j,jp,2)

           enddo ! jp
        enddo ! j
     enddo ! i
     !-----------------------------------------------------------

     !-----------------------------------------------------------
     ! Compute block M(2,1)
     !
     do i=1,n_x_max
        do j=1,n_blend
           do jp=1,n_blend

              ij  = i + (j-1)*n_x + n_x*n_blend
              ijp = i + (jp-1)*n_x 

              k_counter = k_counter+1

              indx_maxwell(k_counter)           = ij
              indx_maxwell(k_counter+n_maxwell) = ijp

              m_maxwell(k_counter) = -imp(i,j,jp,3)

           enddo ! jp
        enddo ! j
     enddo ! i
     !-----------------------------------------------------------

  endif

  ! Now B_par

  if (n_field > 2) then

     !-----------------------------------------------------------
     ! Compute block M(3,3)
     !
     do i=1,n_x_max
        do i_diff=-m_gyro,m_gyro-i_gyro

           ip = i_loop(i+i_diff)

           if (ip >= 1 .and. ip <= n_x) then

              do j=1,n_blend
                 do jp=1,n_blend

                    val = ab_mm(i,i_diff,j,jp)

                    if (i_diff == 0) then
                       val = val-imp(i,j,jp,8) 
                    endif

                    ij  = i + (j-1)*n_x   + 2*n_x*n_blend
                    ijp = ip + (jp-1)*n_x + 2*n_x*n_blend

                    k_counter = k_counter+1

                    indx_maxwell(k_counter)           = ij
                    indx_maxwell(k_counter+n_maxwell) = ijp

                    m_maxwell(k_counter) = val

                 enddo ! jp
              enddo ! j

           endif

        enddo ! i_diff
     enddo ! i
     !-----------------------------------------------------------

     !-----------------------------------------------------------
     ! Compute block M(1,3)
     !
     do i=1,n_x_max
        do i_diff=-m_gyro,m_gyro-i_gyro

           ip = i_loop(i+i_diff)

           if (ip >= 1 .and. ip <= n_x) then

              do j=1,n_blend
                 do jp=1,n_blend

                    val = -2.0*abp_mm(i,i_diff,j,jp)

                    if (i_diff == 0) then
                       val = val+2.0*imp(i,j,jp,6) 
                    endif

                    ij  = i + (j-1)*n_x
                    ijp = ip + (jp-1)*n_x + 2*n_x*n_blend

                    k_counter = k_counter+1

                    indx_maxwell(k_counter)           = ij
                    indx_maxwell(k_counter+n_maxwell) = ijp

                    m_maxwell(k_counter) = val

                 enddo ! jp
              enddo ! j

           endif

        enddo ! ip
     enddo ! i
     !-----------------------------------------------------------

     !-----------------------------------------------------------
     ! Compute block M(3,1)
     !
     do i=1,n_x_max
        do i_diff=-m_gyro,m_gyro-i_gyro

           ip = i_loop(i+i_diff)

           if (ip >= 1 .and. ip <= n_x) then

              do j=1,n_blend
                 do jp=1,n_blend

                    val = abp_mm(i,i_diff,j,jp)

                    if (i_diff == 0) then
                       val = val-imp(i,j,jp,6) 
                    endif

                    ij  = i + (j-1)*n_x  + 2*n_x*n_blend
                    ijp = ip + (jp-1)*n_x

                    k_counter = k_counter+1

                    indx_maxwell(k_counter)           = ij
                    indx_maxwell(k_counter+n_maxwell) = ijp

                    m_maxwell(k_counter) = val

                 enddo ! jp
              enddo ! j

           endif

        enddo ! ip
     enddo ! i
     !-----------------------------------------------------------

     !-----------------------------------------------------------
     ! Compute block M(2,3)
     !
     do i=1,n_x_max
        do j=1,n_blend
           do jp=1,n_blend

              ij  = i + (j-1)*n_x  + n_x*n_blend 
              ijp = i + (jp-1)*n_x + 2*n_x*n_blend

              k_counter = k_counter+1

              indx_maxwell(k_counter)           = ij
              indx_maxwell(k_counter+n_maxwell) = ijp

              m_maxwell(k_counter) = -imp(i,j,jp,5)

           enddo ! jp
        enddo ! j
     enddo ! i
     !-----------------------------------------------------------

     !-----------------------------------------------------------
     ! Compute block M(3,2)
     !
     do i=1,n_x_max
        do j=1,n_blend
           do jp=1,n_blend

              ij  = i + (j-1)*n_x  + 2*n_x*n_blend
              ijp = i + (jp-1)*n_x + n_x*n_blend

              k_counter = k_counter+1

              indx_maxwell(k_counter)           = ij
              indx_maxwell(k_counter+n_maxwell) = ijp

              m_maxwell(k_counter) = -imp(i,j,jp,7)

           enddo ! jp
        enddo ! j
     enddo ! i
     !-----------------------------------------------------------

  endif

  if (boundary_method ==1 .and. n_1(in_1) == 0) then

     !------------------------------------------------
     ! Here, we need to write a line of 1's in last 
     ! row (i=n_x) of blocks M(1,1), M(2,2), M(3,3).  A 
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

           indx_maxwell(k_counter)           = ij
           indx_maxwell(k_counter+n_maxwell) = ijp

           m_maxwell(k_counter) = (1.0,0.0)

        enddo ! ip
     enddo ! j
     !------------------------------------------------

     if (n_field > 1) then

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

              indx_maxwell(k_counter)           = ij
              indx_maxwell(k_counter+n_maxwell) = ijp

              m_maxwell(k_counter) = (1.0,0.0)

           enddo ! ip
        enddo ! j
        !------------------------------------------------

     endif

     if (n_field > 2) then

        !------------------------------------------------
        ! Line of 1's in M(3,3)
        ! 
        i = n_x
        do j=1,n_blend
           do ip=1,n_x

              ! Diagonal in (j,jp)
              jp = j

              ij  = i + (j-1)*n_x   + 2*n_x*n_blend
              ijp = ip + (jp-1)*n_x + 2*n_x*n_blend

              k_counter = k_counter+1

              indx_maxwell(k_counter)           = ij
              indx_maxwell(k_counter+n_maxwell) = ijp

              m_maxwell(k_counter) = (1.0,0.0)

           enddo ! ip
        enddo ! j
        !------------------------------------------------

     endif

  endif

  if (k_counter /= n_maxwell) then
     print *,k_counter,n_maxwell
     call catch_error("Element count mismatch in gyro_build_sparse_maxwell")
  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_build_sparse_maxwell done]'
  endif

end subroutine gyro_build_sparse_maxwell

