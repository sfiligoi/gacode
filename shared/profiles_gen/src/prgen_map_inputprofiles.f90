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

  integer :: i
  real :: xnew(nx)

  ! Compute rho, bref and arho based on GATO dpsi and q.
  call prgen_get_chi(nx,q,kappa,rmin,dpsi,rho,EXPRO_b_ref,EXPRO_arho)

  ! Align with new (GATO) rho grid, not old rho grid, EXPRO_rho:
  call cub_spline(EXPRO_rho,EXPRO_te,nx,rho,xnew,nx)
  EXPRO_te = xnew
  call cub_spline(EXPRO_rho,EXPRO_ne,nx,rho,xnew,nx)
  EXPRO_ne = xnew
  call cub_spline(EXPRO_rho,EXPRO_z_eff,nx,rho,xnew,nx)
  EXPRO_z_eff = xnew
  call cub_spline(EXPRO_rho,EXPRO_w0,nx,rho,xnew,nx)
  EXPRO_w0 = xnew
  call cub_spline(EXPRO_rho,EXPRO_flow_mom,nx,rho,xnew,nx)
  EXPRO_flow_mom = xnew
  call cub_spline(EXPRO_rho,EXPRO_pow_e,nx,rho,xnew,nx)
  EXPRO_pow_e = xnew
  call cub_spline(EXPRO_rho,EXPRO_pow_i,nx,rho,xnew,nx)
  EXPRO_pow_i = xnew
  call cub_spline(EXPRO_rho,EXPRO_pow_ei,nx,rho,xnew,nx)
  EXPRO_pow_ei = xnew
  call cub_spline(EXPRO_rho,EXPRO_flow_beam,nx,rho,xnew,nx)
  EXPRO_flow_beam = xnew
  call cub_spline(EXPRO_rho,EXPRO_flow_wall,nx,rho,xnew,nx)
  EXPRO_flow_wall = xnew
  call cub_spline(EXPRO_rho,EXPRO_ptot,nx,rho,xnew,nx)
  EXPRO_ptot = xnew
  do i=1,5
     call cub_spline(EXPRO_rho,EXPRO_ni(i,:),nx,rho,xnew,nx)
     EXPRO_ni(i,:) = xnew(:)
     call cub_spline(EXPRO_rho,EXPRO_ti(i,:),nx,rho,xnew,nx)
     EXPRO_ti(i,:) = xnew(:)
     call cub_spline(EXPRO_rho,EXPRO_vtor(i,:),nx,rho,xnew,nx)
     EXPRO_vtor(i,:) = xnew(:)
     call cub_spline(EXPRO_rho,EXPRO_vpol(i,:),nx,rho,xnew,nx)
     EXPRO_vpol(i,:) = xnew(:)
  enddo

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

  ! Deallocate arrays
  call EXPRO_alloc('./',0)

end subroutine prgen_map_inputprofiles

