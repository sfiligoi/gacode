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

  use prgen_globals

  implicit none

  integer :: i
  integer :: ip
  real, dimension(nx) :: dphidpsi

  !--------------------------------------------------------------------
  ! Calculate integrated powers from input sources
  !
  pow_e(1)     = 0.0
  pow_i(1)     = 0.0
  pow_ei(1)    = 0.0
  flow_mom(1)  = 0.0
  flow_beam(1) = 0.0

  pow_e_fus(1) = 0.0 
  pow_i_fus(1) = 0.0

  do i=2,nx

     ! Total powers to electrons and ions "per zone"
     ! Integrated power is thus a partial sum.
     ! Factor of 1e-6 converts plasmastate (W) to input.profiles (MW).

     pow_e(i) = pow_e(i-1)+1e-6*plst_pe_trans(i-1)
     pow_i(i) = pow_i(i-1)+1e-6*plst_pi_trans(i-1)

     pow_i_fus(i) = pow_i_fus(i-1)+1e-6*plst_pfusi(i-1)
     pow_e_fus(i) = pow_e_fus(i-1)+1e-6*plst_pfuse(i-1)

     ! Collisional exchange
     !
     ! plst_qie : power from ions to electrons
     ! pow_ei   : power from electrons to ions
     ! 
     ! Thus, we need negative sign here:
     pow_ei(i) = pow_ei(i-1)-1e-6*plst_qie(i-1)

     ! Momentum source
     !
     ! tq_trans already in Nm.
     ! COORDINATES: -ipccw accounts for plasmastate toroidal angle convention
     flow_mom(i) = flow_mom(i-1)+plst_tq_trans(i-1)*(-ipccw)

     ! Particle source
     !
     ! MW/keV = 0.624e22/s

     flow_beam(i) = flow_beam(i-1)+plst_sn_trans(i-1)/0.624e22

  enddo
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! COORDINATES: set sign of poloidal flux 
  dpsi(:) = abs(dpsi(:))*(-ipccw)
  ! 
  ! Convert potential to Omega
  !
  ! omega0 = -c d(Phi)/dpsi
  !
  ! phi [statvolt] = (10/3) phi (kV)
  ! psi [Maxwell]  = 10^8 psi [Weber]
  ! c [cm/s] = 2.9979e10
  !
  ! NOTE: dpsi = plst_psipol-plst_psipol(1)
  !
  call bound_deriv(dphidpsi,plst_epot,dpsi,nx)
  !
  omega0 = -2.9979e10*dphidpsi*(10.0/3.0)/1e8
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
  ! COORDINATES: set sign of q
  vec(4,:)  = abs(q(:))*ipccw*btccw
  vec(5,:)  = kappa(:)
  vec(6,:)  = delta(:)
  vec(7,:)  = plst_ts(:,1)
  vec(8,:)  = plst_ns(:,1)*1e-19
  vec(9,:)  = plst_zeff(:)
  vec(10,:) = omega0(:) 
  vec(11,:) = flow_mom(:)
  vec(12,:) = pow_e(:)
  vec(13,:) = pow_i(:)
  vec(14,:) = pow_ei(:)
  vec(15,:) = zeta(:)
  vec(16,:) = flow_beam(:)
  vec(17,:) = 0.0 ! flow_wall
  vec(18,:) = zmag(:)
  vec(19,:) = plst_ptowb ! total pressure, thermal + fast ion
  ! COORDINATES: This poloidal flux has correct sign (see above).
  vec(20,:) = dpsi(:)

  print *,plst_dp1_nspec_th+1

  ! ni,ti
  do i=1,5
     ip = reorder_vec(i)
     if (ip < plst_dp1_nspec_th+1) then
        vec(21+i-1,:) = plst_ns(:,ip+1)*1e-19
        vec(26+i-1,:) = plst_ts(:,ip+1)
     endif
  enddo

  ! vphi
  do i=1,5
     ip = reorder_vec(i)+1
     if (ip <= plst_dp1_nspec_th) then
        if (trim(plst_all_name(ip)) == 'C') then
           ! COORDINATES: -ipccw accounts for plasmastate toroidal angle convention
           vec(31+i-1,:) = -ipccw*plst_omegat(:)*(rmaj(:)+rmin(:))
        endif
     endif
  enddo

  ! vpol
  vec(36:40,:) = 0.0

  ! Extra powers
  vec(41,:) = pow_e_fus(:)
  vec(42,:) = pow_i_fus(:)
  !---------------------------------------------------------

  ! Ion reordering diagnostics

  print '(a)','INFO: (prgen) Found these ion species:'
  do i=2,plst_dp1_nspec_all
     print '(t6,i2,1x,3(a))', i-1,trim(plst_all_name(i))
  enddo
  print '(a)','INFO: (prgen) Created these species:'
  do i=1,5
     ip = reorder_vec(i)
     if (ip >= plst_dp1_nspec_all) then
        print '(t6,i2,1x,3(a))',i,'[null]'
     else
        print '(t6,i2,1x,3(a))',i,trim(plst_all_name(ip+1))
     endif
  enddo

end subroutine prgen_map_plasmastate
