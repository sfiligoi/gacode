!-------------------------------------------------------
! get_g_squared.f90 [caller: get_diffusivity.f90]
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

subroutine get_g_squared

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  real, dimension(2,n_x) :: moment
  !
  complex, dimension(i1_buffer:i2_buffer,n_kinetic) :: g_loc 
  complex, dimension(n_x) :: phi_loc
  complex, dimension(n_x,n_kinetic) :: ave_g
  !
  !--------------------------------------------------  

  !--------------------------------------------------------
  moment(:,:) = 0.0
  !--------------------------------------------------------

  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek) 
     k  = nek_k(p_nek)

     ck = class(k)

     do m=1,n_stack

        m0 = m_phys(ck,m)

        g_loc(:,:) = (0.0,0.0)

        !---------------------------------------------------------
        ! Nonadiabatic distribution, g.
        !
        do is=1,n_kinetic
           do i=1,n_x
              g_loc(i,is) = h(m,i,p_nek_loc,is)+&
                   z(is)*alpha_s(is,i)*gyro_u(m,i,p_nek_loc,is)
           enddo
        enddo
        !
        ! Potential (not gyroaveraged).
        do i=1,n_x
           phi_loc(i) = field_tau(m,i,p_nek_loc,1)
        enddo
        !---------------------------------------------------------

        !---------------------------------------------------------
        ! Compute gyroaverages
        !
        ave_g(:,:) = 0.0

        do i=1,n_x
           do i_diff=-m_gyro,m_gyro-i_gyro

              ip = i+i_diff

              ave_g(i,1:n_gk) = ave_g(i,1:n_gk)+ &
                   w_gyro(m0,i_diff,i,p_nek_loc,1:n_gk)*&
                   g_loc(i_loop(ip),1:n_gk)

           enddo ! i_diff
        enddo ! i

        if (electron_method == 2) then
           ave_g(1:n_x,indx_e) = g_loc(1:n_x,indx_e)
        endif
        !---------------------------------------------------------

        !---------------------------------------------------------
        ! Velocity-space integration
        ! 
        do is=1,n_kinetic
           do i=1,n_x
              moment(1,i) = moment(1,i)+ &
                   real(ave_g(i,is)*conjg(ave_g(i,is)))*w_p(ie,i,k,is)
           enddo ! i
        enddo ! is

        moment(1,:) = moment(1,:)/2.0

        ! 1/2 -> conforms to definition of entropy <|h|**2/2>

        do i=1,n_x
           moment(2,i) = moment(2,i)+ &
                real(phi_loc(i)*conjg(phi_loc(i)))*w_p(ie,i,k,1)
        enddo ! i
        !---------------------------------------------------------

     enddo ! m

  enddo ! p_nek

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
     print *,'[get_g_squared called]'
  endif

end subroutine get_g_squared
