!------------------------------------------------------------
! prgen_map_iterdb.f90
!
! PURPOSE:
!  Map native iterdb data onto input.profiles standard.  
!  
! NOTES:
!  - This routine is common to both text and NetCDF formats.
!  - See map_plasmastate.f90 for analogous routine for 
!    plasmastate data.
!------------------------------------------------------------

subroutine prgen_map_iterdb

  use prgen_read_globals

  implicit none

  integer :: i
  integer :: ip
  integer :: n0
  real :: xoh_exp
  real :: xrad_exp
  real :: xfus_exp
  real :: xwdot

  do i=1,onetwo_nj
     rho(i) = (i-1)/(onetwo_nj-1.0)
  enddo

  ! Volume integrations to obtain integrated 
  ! powers and flows

  xoh_exp  = 1.0 !counts powe_oh_exp
  xrad_exp = 1.0 !counts powe_rad_exp
  xfus_exp = 1.0 !counts any non-zero powe_fus_exp and powi_fus_exp
  xwdot    = 1.0 !discounts any non-zero wdot_exp

  call volint(1e-6*onetwo_qbeame,powe_beam_exp)
  call volint(1e-6*onetwo_qrfe,powe_rf_exp)
  call volint(1e-6*onetwo_qohm,powe_oh_exp)
  call volint(1e-6*onetwo_qrad,powe_rad_exp)
  call volint(1e-6*onetwo_qione,powe_ion_exp)
  call volint(1e-6*onetwo_dpedt,powe_wdot_exp)
  call volint(1e-6*onetwo_qfuse,powe_fus_exp)
  call volint(1e-6*onetwo_qbeami,powi_beam_exp)
  call volint(1e-6*onetwo_qrfi,powi_rf_exp)
  call volint(1e-6*onetwo_qioni,powi_ion_exp)
  call volint(1e-6*onetwo_qcx,powi_cx_exp)
  call volint(1e-6*onetwo_dpidt,powi_wdot_exp)
  call volint(1e-6*onetwo_qfusi,powi_fus_exp)
  ! Note sign change to keep pow_ei_exp positive for te>ti
  call volint(-1e-6*onetwo_qdelt,pow_ei_exp)

  ! Wall ion flow (scale is undefined) (MW/KeV = Kamp)
  call volint(kevdsecpmw*sion_d,flow_wall_exp)

  ! Beam ion flow
  call volint(kevdsecpmw*(onetwo_sbeam+sbcx_d),flow_beam)

  ! Torque (TAM flow) (nt-m)
  call volint(onetwo_storqueb,flow_mom)

  ! Total transport power (MW) to electrons: pow_e

  ! powe_beam_exp (MW): NBI power added to electrons
  ! powe_rf_exp (MW)  : RF power added to electrons
  ! powe_oh_exp (MW)  : ohmic power added to electrons
  ! powe_rad-exp (MW) : [negative] radiation loss added to electrons
  ! powe_ion_exp (MW) : [negative] neutral ionization loss added
  ! powe_wdot_exp (MW): [zero] W-dot power subtracted
  ! pow_ei_exp (MW)   : [positive for Te > Ti] power *subtracted*
  ! powe_fus_exp (MW) : [zero] any fusion power added

  pow_e(:) = powe_beam_exp(:) &     
       +powe_rf_exp(:) &                
       +xoh_exp*powe_oh_exp(:) &      
       +xrad_exp*powe_rad_exp(:) &    
       +powe_ion_exp(:) &
       -(1.0-xwdot)*powe_wdot_exp(:) &
       -pow_ei_exp(:) &
       +xfus_exp*powe_fus_exp(:)

  ! Total transport power (MW) to ions: pow_i

  ! powi_beam_exp (MW): NBI power added to electrons
  ! powi_rf_exp (MW)  : RF power added to electrons
  ! powi_ion_exp (MW) : [positive] neutral ionization loss added
  ! pow_i (MW)        : [negative] neutral charge exch loss added
  ! powi_wdot_exp (MW): [zero] W-dot power subtracted
  ! pow_ei_exp (MW)   : [positive for Te > Ti] power *added*
  ! powi_fus_exp (MW) : [zero] any fusion power added


  pow_i(:) = powi_beam_exp(:) &
       +powi_rf_exp(:) &
       +powi_ion_exp(:) &
       +powi_cx_exp(:) &
       -(1.0-xwdot)*powi_wdot_exp(:) &
       +pow_ei_exp(:) &
       +xfus_exp*powi_fus_exp(:)


  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  allocate(vec(n_indx,onetwo_nj))
  vec(:,:) = 0.0
  !
  vec(1,:)  = rho(:)
  vec(2,:)  = rmin(:)
  vec(3,:)  = rmaj(:)
  vec(4,:)  = q(:)
  vec(5,:)  = kappa(:)
  vec(6,:)  = delta(:)
  vec(7,:)  = onetwo_te(:)
  vec(8,:)  = onetwo_ene(:)*1e-19
  vec(9,:)  = onetwo_zeff(:)
  vec(11,:) = flow_mom(:)
  vec(12,:) = pow_e(:)
  vec(13,:) = pow_i(:)
  vec(14,:) = pow_ei_exp(:)
  vec(15,:) = zeta(:)
  vec(16,:) = flow_beam(:)
  vec(17,:) = flow_wall_exp(:)
  vec(18,:) = zmag(:)
  vec(19,:) = 0.0
  vec(20,:) = dpsi(:)

  !-----------------------------------------------------------------
  ! Construct ion densities and temperatures with reordering
  ! in general case.  Use vphi and vpol as temporary arrays.
  !
  do i=1,onetwo_nion
     ! ni
     vec(31+i-1,:) = onetwo_enion(:,i)*1e-19
     ! Ti
     vec(36+i-1,:) = onetwo_ti(:)
  enddo
  do i=1,onetwo_nbion
     ! ni
     vec(31+i+onetwo_nion-1,:) = onetwo_enbeam(:,i)*1e-19
     ! Ti: T[keV] = (p/n)[J]/1.6022e-16[J/eV]
     vec(36+i+onetwo_nion-1,:) = onetwo_pressb(:,i)/onetwo_enbeam(:,i)/&
          1.6022e-16
  enddo

  ! reorder
  do i=1,5 
     vec(20+i,:) = vec(30+reorder_vec(i),:)
     vec(25+i,:) = vec(35+reorder_vec(i),:)
  enddo

  ! vphi
  vec(31:35,:) = 0.0

  ! Carbon velocity at theta=0 ONLY in typical case
  ! of one main ion with Carbon impurity

   allocate(vphi_carbon(nx))

  if (trim(onetwo_namei(1)) == 'c') then
     ! Negative sign to account for DIII-D convention
     vphi_carbon(:) = -onetwo_angrot(:)*(rmaj(:)+rmin(:))
  else
     vphi_carbon(:) = 0.0
  endif

  ! Insert carbon toroidal velocity
  do i=1,5
     if (reorder_vec(i) == onetwo_nprim+1) then
        vec(30+i,:) = vphi_carbon(:)
     endif
  enddo

  ! vpol
  vec(36,:) = 0.0
  vec(37,:) = 0.0
  vec(38,:) = 0.0
  vec(39,:) = 0.0
  vec(40,:) = 0.0

  !---------------------------------------------------
  ! Read the cer file and overlay
  !
  if (cer_file /= "null") then
     allocate(vpolc_exp(nx))
     allocate(vtorc_exp(nx))
     call prgen_read_cer
     vec(10,:) = omega0(:)
     vec(31+onetwo_nprim,:) = vtorc_exp(:)
     vec(36+onetwo_nprim,:) = vpolc_exp(:)
  endif
  !---------------------------------------------------

  ! Additional transport power breakdown
  if (n_indx2 == 15) then

     allocate(vec2(n_indx2,onetwo_nj))
     vec2(:,:) = 0.0

     vec2(1,:) = powe_beam_exp(:)
     vec2(2,:) = powe_RF_exp(:)
     vec2(3,:) = powe_oh_exp(:)
     vec2(4,:) = powe_rad_exp(:)
     vec2(5,:) = powe_ion_exp(:)
     vec2(6,:) = powe_wdot_exp(:)
     vec2(7,:) = powe_fus_exp(:)
     vec2(8,:) = pow_e(:)      !net transport power in electrons
     vec2(9,:) = pow_ei_exp(:)    !electron to ion energy exchange
     vec2(10,:) = pow_i(:)     !net transport power in ions
     vec2(11,:) = powi_beam_exp(:)
     vec2(12,:) = powi_ion_exp(:)
     vec2(13,:) = powi_wdot_exp(:)
     vec2(14,:) = powi_fus_exp(:)
     vec2(15,:) = powi_cx_exp(:)

  endif
  !---------------------------------------------------------

  ion_name(1:onetwo_nprim) = onetwo_namep(1:onetwo_nprim)

  n0 = onetwo_nprim

  ion_name(1+n0:onetwo_nimp+n0) = &
       onetwo_namei(1:onetwo_nimp)

  n0 = n0+onetwo_nimp

  ion_name(1+n0:onetwo_nbion+n0) = &
       onetwo_nameb(1:onetwo_nbion)

  ! Sign of poloidal flux

  signpsi = abs(onetwo_psi(nx)-onetwo_psi(1))/&
       (onetwo_psi(nx)-onetwo_psi(1))

end subroutine prgen_map_iterdb

!------------------------------------------------------------------
! Routine to perform Volume integration
!------------------------------------------------------------------

subroutine volint(f,fdv)

  use prgen_read_globals, &
       only : onetwo_rho_grid,pi,onetwo_R0,onetwo_hcap,onetwo_nj

  implicit none 

  integer :: i
  real, intent(in) :: f(onetwo_nj)
  real, intent(out) :: fdv(onetwo_nj) 
  real :: drho
  real :: dvoldr_p,dvoldr_m

  fdv(1) = 0.0

  do i=2,onetwo_nj

     drho = onetwo_rho_grid(i)-onetwo_rho_grid(i-1)

     dvoldr_p = 2.0*pi*onetwo_rho_grid(i)*2.0*pi*onetwo_R0*onetwo_hcap(i)
     dvoldr_m = 2.0*pi*onetwo_rho_grid(i-1)*2.0*pi*onetwo_R0*onetwo_hcap(i-1)

     fdv(i) = fdv(i-1)+0.5*(dvoldr_p*f(i)+dvoldr_m*f(i-1))*drho

  enddo

end subroutine volint

