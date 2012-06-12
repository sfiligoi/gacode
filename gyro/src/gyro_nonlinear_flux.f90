!--------------------------------------------------------------------------
! gyro_nonlinear_flux.f90 [caller: gyro_write_master]
!
! PURPOSE:
!  Evaluate primitive nonlinear fluxes and turbulent exchange using
!  Sugama 1998 forms, where 
!
!  h_a^(Sugama)(R) = H_a(R) -> cap_h
!
! NOTES:
!  Index 1: Particle flux (gyroBohm units)
!  Index 2: Energy flux (gyroBohm units)
!  Index 3: Total momentum flux (gyroBohm units)
!  Index 4: Turbulent exchange flux/length (gyroBohm units)
!--------------------------------------------------------------------------

subroutine gyro_nonlinear_flux

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  real, dimension(n_x,n_kinetic,n_field,p_moment,2) :: moment
  real, dimension(n_kinetic,4) :: excparts
  real, dimension(n_kinetic,3) :: momparts
  real, dimension(4) :: exctemp
  real, dimension(3) :: momtemp

  complex, dimension(i1_buffer:i2_buffer,n_kinetic) :: cap_h
  complex, dimension(i1_buffer:i2_buffer,n_kinetic) :: h1
  complex, dimension(i1_buffer:i2_buffer,n_kinetic) :: h2
  complex, dimension(i1_buffer:i2_buffer,n_kinetic) :: gyro_h_cap
  complex, dimension(i1_buffer:i2_buffer,n_kinetic) :: gyro_h_cap_dot
  complex, dimension(n_x) :: ikrho
  !--------------------------------------------------  


  moment(:,:,:,:,:) = 0.0
  excparts(:,:)     = 0.0
  momparts(:,:)     = 0.0
  exctemp(:)        = 0.0
  momtemp(:)        = 0.0

  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek) 
     k  = nek_k(p_nek)

     ck = class(k)

     do m=1,n_stack

        cap_h(:,:) = (0.0,0.0)
        cap_h(1:n_x,:) = h_cap(m,:,p_nek_loc,:)

        do i=1,n_x
           h1(i,:) = h_cap(m,i,p_nek_loc,:)
           h2(i,:) = h_cap_dot(m,i,p_nek_loc,:)
        enddo
        if (boundary_method == 1) then
           do i=1-m_gyro,0
              h1(i,:) = h_cap(m,i+n_x,p_nek_loc,:)
              h2(i,:) = h_cap_dot(m,i+n_x,p_nek_loc,:)
           enddo
           do i=n_x+1,n_x+m_gyro-i_gyro
              h1(i,:) = h_cap(m,i-n_x,p_nek_loc,:)
              h2(i,:) = h_cap_dot(m,i-n_x,p_nek_loc,:)
           enddo
        else
           do i=1-m_gyro,0
              h1(i,:) = 0.0
              h2(i,:) = 0.0
           enddo
           do i=n_x+1,n_x+m_gyro
              h1(i,:) = 0.0
              h2(i,:) = 0.0
           enddo
        endif

        gyro_h_cap(:,:) = (0.0,0.0)
        do is=1,n_kinetic
           do i=1,n_x
              do i_diff=-m_gyro,m_gyro-i_gyro
                 gyro_h_cap(i,is) = gyro_h_cap(i,is)+&
                      w_gyro(m,i_diff,i,p_nek_loc,is)*h1(i+i_diff,is)
              enddo ! i_diff
           enddo
        enddo

        gyro_h_cap_dot(:,:) = (0.0,0.0)
        do is=1,n_kinetic
           do i=1,n_x
              do i_diff=-m_gyro,m_gyro-i_gyro
                 gyro_h_cap_dot(i,is) = gyro_h_cap_dot(i,is)+&
                      w_gyro(m,i_diff,i,p_nek_loc,is)*h2(i+i_diff,is)
              enddo ! i_diff
           enddo
        enddo

        !-----------------------------------------------------
        ! Compute basic fluxes
        !
        ikrho(:) = i_c*krho_i(in_1,:)
        !
        do ix=1,n_field
           do is=1,n_kinetic
              do i=1,n_x

                 ! Keep running total of moment over p_nek_loc and m:

                 ! Moment 1: density
                 moment(i,is,ix,1,ck) = moment(i,is,ix,1,ck)+real(-ikrho(i)* &
                      conjg(cap_h(i,is))*gyro_uv(m,i,p_nek_loc,is,ix)*w_p(ie,i,k,is))

                 ! Moment 2: energy
                 moment(i,is,ix,2,ck) = moment(i,is,ix,2,ck)+real(-ikrho(i)* &
                      conjg(cap_h(i,is))*gyro_uv(m,i,p_nek_loc,is,ix)*w_p(ie,i,k,is))* &
                      energy(ie,is)*tem_s(is,i)

                 ! Moment 3: momentum

                 ! - 1. parallel part

                 momtemp(1) = real(-ikrho(i)* &
                      conjg(cap_h(i,is))*gyro_uv(m,i,p_nek_loc,is,ix)*w_p(ie,i,k,is))* &
                      v_para(m,i,p_nek_loc,is)/mu(is)**2*&
                      bigr_t(i,k,m)*bt_t(i,k,m)/b0_t(i,k,m)

                 ! - 2. perpendicular part

                 momtemp(2) = -real(-ikrho(i)* &
                      conjg(cap_h(i,is))*kyro_uv(m,i,p_nek_loc,is,ix)*w_p(ie,i,k,is))* &
                      v_perp(m,i,p_nek_loc,is)/mu(is)**2*&
                      bigr_t(i,k,m)*bp_t(i,k,m)/b0_t(i,k,m)

                 ! - 3. convective_part

                 momtemp(3) = real(-ikrho(i)* &
                      conjg(cap_h(i,is))*gyro_uv(m,i,p_nek_loc,is,ix)*w_p(ie,i,k,is))* &
                      mach_s(i)*bigr_t(i,k,m)/rmaj_s(i)/mu(is)**2*&
                      sqrt(tem_s(indx_e,i))

                 moment(i,is,ix,3,ck) = moment(i,is,ix,3,ck)+sum(momtemp(:))

                 ! Momentum breakdown (sum over i, ix, ck)

                 momparts(is,:) = momparts(is,:)+momtemp(:)/n_x 

                 ! Moment 4: Exchange
                 moment(i,is,ix,4,ck) = moment(i,is,ix,4,ck)+z(is)*real( &
                      conjg(cap_h(i,is))*gyro_uv_dot(m,i,p_nek_loc,is,ix)*w_p(ie,i,k,is))

                 ! Exchange breakdown (sum over i, ix, ck): ** valid for ES only **

                 exctemp(1) = z(is)*real( &
                      conjg(cap_h(i,is))*gyro_uv_dot(m,i,p_nek_loc,is,ix)*w_p(ie,i,k,is))
                 exctemp(2) = z(is)*real( &
                      conjg(h_cap_dot(m,i,p_nek_loc,is))*gyro_uv(m,i,p_nek_loc,is,ix)*w_p(ie,i,k,is))
                 exctemp(3) = z(is)*real( &
                      conjg(gyro_h_cap(i,is))*field_tau_dot(m,i,p_nek_loc,ix)*w_p(ie,i,k,is))
                 exctemp(4) = z(is)*real( &
                      conjg(gyro_h_cap_dot(i,is))*field_tau(m,i,p_nek_loc,ix)*w_p(ie,i,k,is))

                 excparts(is,:) = excparts(is,:)+exctemp(:)/n_x 

              enddo ! i
           enddo ! is
        enddo ! ix
        !-----------------------------------------------------

     enddo ! m

  enddo ! p_nek

  !----------------------------------------------------------------------
  ! ALLREDUCE to obtain full velocity-space sums.  The nonlinear
  ! fluxes will be correct on all processors.
  !
  call MPI_ALLREDUCE(moment(:,:,:,:,1), &
       nonlinear_flux_passing, &
       size(nonlinear_flux_passing), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  call MPI_ALLREDUCE(moment(:,:,:,:,2), &
       nonlinear_flux_trapped, &
       size(nonlinear_flux_trapped), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  call MPI_ALLREDUCE(momparts(:,:), &
       nonlinear_flux_momparts, &
       size(nonlinear_flux_momparts), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  call MPI_ALLREDUCE(excparts(:,:), &
       nonlinear_flux_excparts, &
       size(nonlinear_flux_excparts), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)
  !---------------------------------------------------------------------- 

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_nonlinear_flux called]'
  endif

end subroutine gyro_nonlinear_flux
