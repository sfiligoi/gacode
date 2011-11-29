!------------------------------------------------------------
! prgen_map_peqdsk.f90
!
! PURPOSE:
!  Map native peqdsk data onto input.profiles standard.  
!------------------------------------------------------------

subroutine prgen_map_peqdsk

  use prgen_read_globals

  implicit none

  integer :: i

  !---------------------------------------------------------
  ! Compute rho from: d chi_t = q d psi_p
  ! use the trapezoidal rule
  ! rho = sqrt(2 * chi_t / bref)
  rho(1) = 0.0
  do i=2,peqdsk_nj
     rho(i) = rho(i-1) + 0.5*(q_gato(i)+q_gato(i-1))*(dpsi(i)-dpsi(i-1))
  enddo
  ! choose b_ref = 1
  peqdsk_bref = 1.0
  peqdsk_arho = sqrt(2.0 * rho(peqdsk_nj) / peqdsk_bref)
  ! normalized rho
  rho(:) = sqrt(rho(:) / rho(peqdsk_nj))
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  allocate(vec(n_indx,peqdsk_nj))
  vec(:,:) = 0.0
  !
  vec(1,:)  = rho(:)
  vec(2,:)  = rmin(:)
  vec(3,:)  = rmaj(:)
  vec(4,:)  = q(:)
  vec(5,:)  = kappa(:)
  vec(6,:)  = delta(:)
  vec(7,:)  = peqdsk_te(:)
  vec(8,:)  = peqdsk_ne(:)*10
  vec(9,:)  = 0.0      ! zeff
  vec(10,:) = -1e3*peqdsk_omegat(:) ! Omega_tor
  vec(11,:) = 0.0      ! flow_mom
  vec(12,:) = 0.0      ! pow_e
  vec(13,:) = 0.0      ! pow_i 
  vec(14,:) = 0.0      ! pow_ei_exp
  vec(15,:) = zeta(:)
  vec(16,:) = 0.0      ! flow_beam
  vec(17,:) = 0.0      ! flow_wall_exp
  vec(18,:) = zmag(:)  
  vec(19,:) = 0.0      
  vec(20,:) = dpsi(:)

  ! ni
  i=1
  vec(21+i-1,:) = peqdsk_ni(:)*10

  ! ti
  vec(26,:) = peqdsk_ti(:)

  ! vphi
  vec(31,:) = 0.0
  vec(32,:) = peqdsk_omegat(:)*1000*(rmaj(:)+rmin(:))
  vec(33,:) = 0.0
  vec(34,:) = 0.0
  vec(35,:) = 0.0

  ! vpol
  vec(36,:) = 0.0
  vec(37,:) = 0.0
  vec(38,:) = 0.0
  vec(39,:) = 0.0
  vec(40,:) = 0.0

  signpsi = abs(dpsi(peqdsk_nj)-dpsi(1))/&
       (dpsi(peqdsk_nj)-dpsi(1))
  
end subroutine prgen_map_peqdsk
