!----------------------------------------------------------
! gyro_nonlinear_flux_velocity.f90
!
! PURPOSE:
!  Calculate the velocity dependence of the flux integrals.
!----------------------------------------------------------

subroutine gyro_nonlinear_flux_velocity

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  real, dimension(n_energy,n_lambda,n_kinetic,n_field,n_moment) :: moment
  !
  complex, dimension(i1_buffer:i2_buffer,n_kinetic) :: hh
  complex, dimension(n_x,n_energy,n_lambda,n_kinetic,n_field) :: ave_h
  !
  complex, dimension(n_x,n_energy,n_lambda,n_field) :: f_w
  complex, dimension(n_x) :: ikrho
  !--------------------------------------------------  


  if (n_field == 3) then
     call catch_error(&
          'ERROR: Velocity-space diagnostics not available for N_FIELD=3.')
  endif

  !--------------------------------------------------------
  moment(:,:,:,:,:) = 0.0
  !--------------------------------------------------------

  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek) 
     k  = nek_k(p_nek)

     ck = class(k)

     do m=1,n_stack

        !---------------------------------------------------------
        !
        hh(:,:) = (0.0,0.0)
        do i=1,n_x
           hh(i,:) = h(m,i,p_nek_loc,:)
        enddo

        m0 = m_phys(ck,m)

        !-----------------------------
        ! Index reminder:
        !
        ! ave_h(n_x,n_kinetic,n_field)

        ave_h(:,:,:,:,1) = (0.0,0.0)

        do i=1,n_x
           do i_diff=-m_gyro,m_gyro-i_gyro

              ip = i+i_diff

              ave_h(i,ie,k,1:n_gk,1) = ave_h(i,ie,k,1:n_gk,1)+ &
                   w_gyro(m0,i_diff,i,p_nek_loc,1:n_gk)*&
                   hh(i_loop(ip),1:n_gk)

           enddo ! i_diff
        enddo ! i

        if (electron_method == 2) then
           ave_h(:,ie,k,n_spec,1) = hh(1:n_x,n_spec)
        endif

        if (n_field == 2) then
           do is=1,n_kinetic
              ave_h(:,ie,k,is,2) = &
                   ave_h(:,ie,k,is,1)*(-v_para(m,:,p_nek_loc,is))
           enddo
        endif

        f_w(:,:,:,:) = (0.0,0.0)
        if (n_field == 1) then
           do j=1,n_blend
              f_w(:,ie,k,1) = f_w(:,ie,k,1)+&
                   conjg(field_blend(j,:,1))*cs_blend(j,m0,:,p_nek_loc)
           enddo ! j
        else
           do j=1,n_blend
              f_w(:,ie,k,1) = f_w(:,ie,k,1) &
                   +conjg(field_blend(j,:,1))*cs_blend(j,m0,:,p_nek_loc)
              f_w(:,ie,k,2) = f_w(:,ie,k,2) &
                   +conjg(field_blend(j,:,2))*cs_blend(j,m0,:,p_nek_loc)
           enddo ! j
        endif

        !-------------------------------------------------
        ! Keep running total of moment(:,:) over p_nek_loc 
        ! and m:
        !-------------------------------------------------

        ikrho(:) = i_c*krho_i(in_1,:)

        !-----------------------------------------------------
        ! Main diffusivities:
        !
        do ix=1,n_field
           do is=1,n_kinetic

              ! Moment 1: density
              moment(ie,k,is,ix,1) = moment(ie,k,is,ix,1)+&
                   sum(real(ikrho(:)*f_w(:,ie,k,ix)*&
                   ave_h(:,ie,k,is,ix)))

              ! Moment 2: energy
              moment(ie,k,is,ix,2) = moment(ie,k,is,ix,2)+&
                   sum(real(ikrho(:)*f_w(:,ie,k,ix)*&
                   ave_h(:,ie,k,is,ix))*energy(ie,is)*tem_s(is,:))

           enddo ! is
        enddo ! ix
        !-----------------------------------------------------

     enddo ! m

  enddo ! p_nek

  !-----------------------------------------------------------
  ! Correct for Jacobian
  !
  do ie=1,n_energy
     do k=1,n_lambda
        do is=1,n_kinetic
           moment(ie,k,is,:,:) = moment(ie,k,is,:,:)*&
                sqrt(energy(ie,is))*exp(-energy(ie,is))/w_energy(ie,is)*&
                (1.0/omega(ir_norm,k))/w_lambda(ir_norm,k)
        enddo
     enddo
  enddo
  !-----------------------------------------------------------

  !----------------------------------------------------------------------
  ! ALLREDUCE to obtain full velocity-space sums.  The nonlinear
  ! fluxes will be correct on all processors.
  !
  call MPI_ALLREDUCE(moment, &
       nonlinear_flux_velocity, &
       size(nonlinear_flux_velocity), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[get_nonlinear_flux_energy called]'
  endif

end subroutine gyro_nonlinear_flux_velocity
