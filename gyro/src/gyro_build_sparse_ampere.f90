!------------------------------------------------------
! gyro_build_sparse_ampere.f90 [caller: sparse_solve_*]
!
! PURPOSE:
!  Fill elements of sparse Ampere submatrix.
!------------------------------------------------------

subroutine gyro_build_sparse_ampere

  use gyro_globals
  use gyro_collision_private
  use math_constants

  !------------------------
  implicit none
  !------------------------


  ! Zero only the in_1 components:
  m_ampere(:) = (0.0,0.0)
  indx_ampere = 0

  k_amp = 0

  do i=1,n_x_max
     do i_diff=-mg_dx,mg_dx-ig_dx

        ip = i_loop(i+i_diff)

        if (ip >= 1 .and. ip <= n_x) then

           do j=1,n_blend
              do jp=1,n_blend

                 val = aa_mm(i,i_diff,j,jp) &
                      +coll_vel(i,j,jp)*w_gd0(i_diff)

                 ij  = i + (j-1)*n_x 
                 ijp = ip + (jp-1)*n_x

                 k_amp = k_amp+1

                 indx_ampere(k_amp)          = ij
                 indx_ampere(k_amp+n_ampere) = ijp

                 m_ampere(k_amp) = val

              enddo ! jp
           enddo ! j

        endif

     enddo ! ip
  enddo ! i

  if (boundary_method == 1 .and. n_1(in_1) == 0) then

     ! need line of 1's in last row (i=n_x), since 
     ! all n=0 modes have no average.

     i = n_x
     do j=1,n_blend
        do ip=1,n_x

           ! Diagonal in (j,jp)
           jp = j

           ij  = i + (j-1)*n_x
           ijp = ip + (jp-1)*n_x

           k_amp = k_amp+1

           indx_ampere(k_amp) = ij
           indx_ampere(k_amp+n_ampere) = ijp

           m_ampere(k_amp) = (1.0,0.0)

        enddo ! ip
     enddo ! j

  endif

  if (k_amp /= n_ampere) then
     print *,k_amp,n_ampere
     call catch_error("Element count mismatch in gyro_build_sparse_ampere")
  endif

end subroutine gyro_build_sparse_ampere

