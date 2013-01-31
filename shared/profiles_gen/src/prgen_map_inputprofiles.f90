!------------------------------------------------------------
! prgen_map_inputprofiles.f90
!
! PURPOSE:
!  Map input profiles plus gfile data.
!------------------------------------------------------------

subroutine prgen_map_inputprofiles

  use prgen_globals
  use EXPRO_interface

  implicit none

  ! Compute rho, bref and arho:
  call prgen_get_chi(nx,q_gato,kappa,rmin,dpsi,rho,EXPRO_b_ref,EXPRO_arho)

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  allocate(vec(n_indx,nx))
  vec(1,:)  = rho(:)
  vec(2,:)  = rmin(:)
  vec(3,:)  = rmaj(:)
  vec(4,:)  = q(:)
  vec(5,:)  = kappa(:)
  vec(6,:)  = delta(:)
  vec(7,:)  = EXPRO_te(:)
  vec(8,:)  = EXPRO_ne(:)
  vec(9,:)  = EXPRO_z_eff(:)
  vec(10,:) = EXPRO_w0(:)
  vec(11,:) = EXPRO_flow_mom(:)
  vec(12,:) = EXPRO_pow_e(:)
  vec(13,:) = EXPRO_pow_i(:)
  vec(14,:) = EXPRO_pow_ei(:)
  vec(15,:) = zeta(:)
  vec(16,:) = EXPRO_flow_beam(:)
  vec(17,:) = EXPRO_flow_wall(:)
  vec(18,:) = zmag(:)
  vec(19,:) = EXPRO_ptot(:)
  vec(20,:) = dpsi(:)

  ! Ion temperatures and densities
  vec(21:25,:) = EXPRO_ni(1:5,:)
  vec(26:30,:) = EXPRO_ti(1:5,:)

  ! Ion velocities
  vec(31:35,:) = EXPRO_vtor(1:5,:)
  vec(36:40,:) = EXPRO_vpol(1:5,:)

end subroutine prgen_map_inputprofiles

