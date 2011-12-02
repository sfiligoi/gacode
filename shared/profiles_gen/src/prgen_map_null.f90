!------------------------------------------------------------
! prgen_map_null.f90
!
! PURPOSE:
!  Map shape data only.  
!------------------------------------------------------------

subroutine prgen_map_null

  use prgen_read_globals

  implicit none

  real, dimension(:), allocatable :: chi_t
  integer :: i 

  !---------------------------------------------------------
  ! Compute rho by integrating d(chi_t) = q d(psi_p)
  ! using the trapezoidal rule
  !
  ! NOTE: chi_t/bref)
  !
  allocate(chi_t(nx))  
  !
  chi_t(1) = 0.0
  do i=2,nx
     chi_t(i) = chi_t(i-1) + 0.5*(q_gato(i)+q_gato(i-1))*(dpsi(i)-dpsi(i-1))
  enddo
  !
  ! Even though bref is arbitrary, choose physically sensible values 
  ! rather than, say, bref=1.  This is exact in the limit B is constant 
  ! and the flux surfaces are elliptical with fixed kappa.
  ! 
  null_bref = 2*chi_t(nx)/(kappa(nx)*rmin(nx)**2)
  null_arho = sqrt(kappa(nx))*rmin(nx)
  !
  ! Normalized root of chi_t:
  !
  rho(:) = sqrt(chi_t(:)/chi_t(nx))
  deallocate(chi_t)
  !---------------------------------------------------------

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

