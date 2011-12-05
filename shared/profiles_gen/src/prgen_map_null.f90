!------------------------------------------------------------
! prgen_map_null.f90
!
! PURPOSE:
!  Map shape data only.  
!------------------------------------------------------------

subroutine prgen_map_null

  use prgen_read_globals

  implicit none

  ! Compute rho, bref and arho:
  call prgen_get_chi(nx,q_gato,kappa,rmin,dpsi,rho,null_bref,null_arho)

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
  vec(7,:)  = 1.0
  vec(8,:)  = 1.0
  vec(9,:)  = 0.0
  vec(10,:) = 0.0
  vec(11,:) = 0.0
  vec(12,:) = 0.0
  vec(13,:) = 0.0
  vec(14,:) = 0.0
  vec(15,:) = zeta(:)
  vec(16,:) = 0.0
  vec(17,:) = 0.0
  vec(18,:) = zmag(:)
  vec(19,:) = 0.0
  vec(20,:) = dpsi(:)

  ! Ion temperatures and densities
  vec(21:25,:) = 1.0
  vec(26:30,:) = 1.0

  ! Ion velocities
  vec(31:40,:) = 0.0

 end subroutine prgen_map_null

