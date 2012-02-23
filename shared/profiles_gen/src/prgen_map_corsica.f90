!------------------------------------------------------------
! prgen_map_corsica.f90
!
! PURPOSE:
!  Map native corsica data onto input.profiles standard.  
!------------------------------------------------------------

subroutine prgen_map_corsica

  use prgen_read_globals

  implicit none

  ! Compute rho, bref and arho:
  call prgen_get_chi(nx,q_gato,kappa,rmin,dpsi,rho,peqdsk_bref,peqdsk_arho)

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  allocate(vec(n_indx,corsica_nvals))
  vec(:,:) = 0.0
  !
  vec(1,:)  = corsica_rho(:)
  vec(2,:)  = rmin(:)
  vec(3,:)  = rmaj(:)
  vec(4,:)  = corsica_q(:)
  vec(5,:)  = kappa(:)
  vec(6,:)  = delta(:)
  vec(7,:)  = corsica_te(:)
  vec(8,:)  = corsica_ne(:)*10.
  vec(9,:)  = corsica_zeff(:)
  vec(10,:) = 0.0      ! omega
  vec(11,:) = 0.0      ! flow_mom
  vec(12,:) = 0.0      ! pow_e
  vec(13,:) = 0.0      ! pow_i 
  vec(14,:) = 0.0      ! pow_ei_exp
  vec(15,:) = zeta(:)
  vec(16,:) = 0.0      ! flow_beam
  vec(17,:) = 0.0      ! flow_wall_exp
  vec(18,:) = zmag(:)  
  vec(19,:) = 0.0      ! ptot
  vec(20,:) = dpsi(:)

  ! ni
  vec(21,:) = corsica_ndt(:)*10.
  
  ! ti
  vec(26,:) = corsica_ti(:)

  ! vphi
  vec(31,:) = 0.0
  vec(32,:) = 0.0
  vec(33,:) = 0.0
  vec(34,:) = 0.0
  vec(35,:) = 0.0

  ! vpol
  vec(36,:) = 0.0
  vec(37,:) = 0.0
  vec(38,:) = 0.0
  vec(39,:) = 0.0
  vec(40,:) = 0.0

  signpsi = abs(dpsi(nx)-dpsi(1))/&
       (dpsi(nx)-dpsi(1))
  
end subroutine prgen_map_corsica
