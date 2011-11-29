!------------------------------------------------------------
! prgen_map_null.f90
!
! PURPOSE:
!  Map minimal data.  
!------------------------------------------------------------

subroutine prgen_map_null

  use prgen_read_globals

  implicit none

  integer :: i 

  !---------------------------------------------------------
  ! Compute rho by integrating d(chi_t) = q d(psi_p)
  ! using the trapezoidal rule
  !
  ! NOTE: rho = sqrt(2*chi_t/bref)
  !
  rho(1) = 0.0
  do i=2,nx
     rho(i) = rho(i-1) + 0.5*(q_gato(i)+q_gato(i-1))*(dpsi(i)-dpsi(i-1))
  enddo
  !
  ! Choose b_ref = 1:
  !
  null_bref = 1.0
  null_arho = sqrt(2.0*rho(nx)/null_bref)
  !
  ! Replace rho with normalized rho:
  !
  rho(:) = sqrt(rho(:)/rho(nx))
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  allocate(vec(n_indx,nx))
  vec(:,:) = 0.0
  !
  vec(1,:)  = rho(:)
  vec(2,:)  = rmin(:)
  vec(3,:)  = rmaj(:)
  vec(4,:)  = q(:)
  vec(5,:)  = kappa(:)
  vec(6,:)  = delta(:)
  vec(7,:)  = 0.0
  vec(8,:)  = 0.0
  vec(9,:)  = 0.0
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

 end subroutine prgen_map_null

