!--------------------------------------------------------------
! prgen_map_plasmastate.f90
!
! PURPOSE:
!  Map native plasmastate data onto input.profiles standard.
!
! NOTES:
!  See prgen_map_iterdb.f90 for analogous routine for 
!  iterbd data.
!--------------------------------------------------------------

subroutine prgen_map_plasmastate

  use prgen_read_globals

  implicit none

  integer :: i
  integer :: ip
  real, dimension(nx) :: dphidpsi
  real :: sign_omega0
  real :: sign_plst_omega

  !--------------------------------------------------------------------
  ! Calculate transport sources:
  !
  pow_e(1)     = 0.0
  pow_i(1)     = 0.0
  pow_ei_exp(1)= 0.0
  flow_mom(1)  = 0.0
  flow_beam(1) = 0.0

  do i=2,nx

     ! pow_e(i) = pow_e(i-1)+1e-6*plst_pbe(i-1)
     ! pow_i(i) = pow_i(i-1)+1e-6*plst_pbi(i-1)

     pow_e(i) = pow_e(i-1)+1e-6*plst_pe_trans(i-1)
     pow_i(i) = pow_i(i-1)+1e-6*plst_pi_trans(i-1)

     ! plst_pei_trans ~ pow_ie NOT pow_ei, so minus
     pow_ei_exp(i) = pow_ei_exp(i-1)-1e-6*plst_pei_trans(i-1)

     ! tq_trans already in Nm.

     flow_mom(i) = flow_mom(i-1)+plst_tq_trans(i-1)

     ! MW/keV = 0.624e22/s

     flow_beam(i) = flow_beam(i-1)+plst_sn_trans(i-1)/0.624e22

  enddo
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! Convert potential to Omega
  !
  ! omega0 = -c d(Phi)/dpsi
  !
  ! phi [statvolt] = (10/3) phi (kV)
  ! psi [Maxwell]  = 10^8 psi [Weber]
  ! c [cm/s] = 2.9979e10
  !
  call bound_deriv(dphidpsi,plst_epot,plst_psipol,nx)
  !
  signpsi = abs(plst_psipol(nx)-plst_psipol(1))/&
       (plst_psipol(nx)-plst_psipol(1))
  !
  ! As computed, its not clear if this has the correct sign:
  omega0 = -2.9979e10*dphidpsi*(10.0/3.0)/1e8
  !
  ! Get sign of desired omega computed above
  if (omega0(1) < 0.0) then
     sign_omega0 = -1.0
  else
     sign_omega0 = 1.0
  endif
  !
  ! Get true sign of approximate omega (plst_omegat)
  if (plst_omegat(1) < 0.0) then
     sign_plst_omega = -1.0
  else
     sign_plst_omega = 1.0
  endif
  !
  ! Ensure sign of omega is opposite plst_omegat (since plst 
  ! definition of toroidal angle is opposite GYRO/NEO):
  omega0 = omega0*sign_omega0*(-sign_plst_omega)
  !--------------------------------------------------------------------

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  allocate(vec(n_indx,nx))
  vec(:,:) = 0.0
  !
  vec(1,:)  = plst_rho(:)
  vec(2,:)  = rmin(:)
  vec(3,:)  = rmaj(:)
  vec(4,:)  = q(:)
  vec(5,:)  = kappa(:)
  vec(6,:)  = delta(:)
  vec(7,:)  = plst_ts(:,1)
  vec(8,:)  = plst_ns(:,1)*1e-19
  vec(9,:)  = plst_zeff_th(:)
  vec(10,:) = omega0(:) 
  vec(11,:) = flow_mom(:)
  vec(12,:) = pow_e(:)
  vec(13,:) = pow_i(:)
  vec(14,:) = pow_ei_exp(:)
  vec(15,:) = zeta(:)
  vec(16,:) = flow_beam(:)
  vec(17,:) = 0.0 ! flow_wall
  vec(18,:) = zmag(:)
  vec(19,:) = plst_ptowb ! total pressure, thermal + fast ion
  vec(20,:) = dpsi(:)

  ! ni
  do i=2,min(plst_dp1_nspec_th+1,6)
     ip = reorder_vec(i-1)+1
     vec(21+i-2,:) = plst_ns(:,ip)*1e-19
  enddo

  ! ti
  do i=2,min(plst_dp1_nspec_th+1,6)
     ip = reorder_vec(i-1)+1
     vec(26+i-2,:) = plst_ts(:,ip)
  enddo

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
  !---------------------------------------------------------

  ! Ion reordering diagnostics

  print '(a)','INFO: (prgen) Found these ion species'
  do i=2,plst_dp1_nspec_alla
     ip = reorder_vec(i-1)+1
     if (i <= 6) then
        print '(t6,i2,1x,3(a))',&
             i-1,trim(plst_alla_name(i)),' -> ',trim(plst_alla_name(ip))
     else
        print '(t6,i2,1x,3(a))',&
             i-1,trim(plst_alla_name(i)),' [unmapped]'
     endif
  enddo

end subroutine prgen_map_plasmastate
