!-------------------------------------------------------
! gyro_g_squared.f90 [caller: get_diffusivity.f90]
!
! PURPOSE:
!  Calculates the flux-surface average of g_star*g
!  summed over velocity space and species.
!  
! NOTES:
!  Used to get analog of phi_QL_squared_n (phi_norm) 
!  similar to Staebler v-norm for lindiff_method = 6  
!  (quasilinear (QL) norm).
!--------------------------------------------------------------------------

subroutine gyro_g_squared

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants
  use ompdata

  !---------------------------------------------------
  implicit none
  !
  real, dimension(2,n_x) :: moment
  !
  complex, dimension(i1_buffer:i2_buffer,n_kinetic,n_stack) :: g_loc 
  complex, dimension(n_x) :: phi_loc
  complex, dimension(n_x,n_kinetic) :: ave_g
  !
  !--------------------------------------------------  
!$omp parallel private(p_nek_loc,ie,k,ck,m0)
  !--------------------------------------------------------
  moment(:,ibeg:iend) = 0.0
  !--------------------------------------------------------

  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek) 
     k  = nek_k(p_nek)

     ck = class(k)

!$omp barrier
!$omp single
     g_loc(i1_buffer:0,:,:) = (0.0,0.0)
     g_loc(n_x+1:i2_buffer,:,:) = (0.0,0.0)
!$omp end single
     do m=1,n_stack
        !---------------------------------------------------------
        ! Nonadiabatic distribution, g.
        !
        do is=1,n_kinetic
           do i = ibeg, iend
              g_loc(i,is,m) = h(m,i,p_nek_loc,is)+&
                   z(is)*alpha_s(is,i)*gyro_u(m,i,p_nek_loc,is)
           enddo
        enddo
     end do
!$omp barrier  ! ensure all g_loc values are available

     do m=1,n_stack

        m0 = m_phys(ck,m)

        !
        ! Potential (not gyroaveraged).
        do i = ibeg, iend
           phi_loc(i) = field_tau(m,i,p_nek_loc,1)
        enddo
        !---------------------------------------------------------

        !---------------------------------------------------------
        ! Compute gyroaverages
        !
        do i = ibeg, iend
           ave_g(i,:) = 0.0
           do i_diff=-m_gyro,m_gyro-i_gyro

              ip = i+i_diff

              ave_g(i,1:n_gk) = ave_g(i,1:n_gk)+ &
                   w_gyro(m0,i_diff,i,p_nek_loc,1:n_gk)*&
                   g_loc(i_loop(ip),1:n_gk,m)

           enddo ! i_diff
        enddo ! i

        if (electron_method == 2) then
           ave_g(ibeg:iend,indx_e) = g_loc(ibeg:iend,indx_e,m)
        endif
        !---------------------------------------------------------

        !---------------------------------------------------------
        ! Velocity-space integration
        ! 
        do is=1,n_kinetic
           do i = ibeg, iend
              moment(1,i) = moment(1,i)+ &
                   real(ave_g(i,is)*conjg(ave_g(i,is)))*w_p(ie,i,k,is)
           enddo ! i
        enddo ! is

        moment(1,ibeg:iend) = moment(1,ibeg:iend)/2.0

        ! 1/2 -> conforms to definition of entropy <|h|**2/2>

        do i = ibeg, iend
           moment(2,i) = moment(2,i)+ &
                real(phi_loc(i)*conjg(phi_loc(i)))*w_p(ie,i,k,1)
        enddo ! i
        !---------------------------------------------------------

     enddo ! m

  enddo ! p_nek
!$omp end parallel

  !----------------------------------------------------------------------
  ! ALLREDUCE to obtain full velocity-space sums.  The nonlinear
  ! fluxes will be correct on all processors.
  !
  call MPI_ALLREDUCE(moment, &
       g_squared, &
       size(g_squared), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)
  !----------------------------------------------------------------------

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_g_squared called]'
  endif

end subroutine gyro_g_squared
