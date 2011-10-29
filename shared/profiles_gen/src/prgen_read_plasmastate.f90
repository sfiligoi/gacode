!--------------------------------------------------------------
! prgen_read_plasmastate.f90
!
! PURPOSE:
!  Extract plasmastate data from native NETCDF file.
!--------------------------------------------------------------

subroutine prgen_read_plasmastate

  use prgen_read_globals
  use netcdf

  implicit none

  ! NetCDF variables
  integer :: i
  integer :: ncid
  integer :: varid
  integer :: err
  real :: dummy

  ! Open the file (NF90_NOWRITE means read-only)
  err = nf90_open(raw_data_file,NF90_NOWRITE,ncid)

  ! Shot
  plst_tag = 'shot_number'
  err = nf90_inq_varid(ncid,trim(plst_tag),varid)
  err = nf90_get_var(ncid,varid,plst_shot_number)
  if (verbose_flag == 1) print *,err,plst_tag,plst_shot_number

  ! Tokamak
  plst_tag = 'tokamak_id'
  err = nf90_inq_varid(ncid,trim(plst_tag),varid)
  err = nf90_get_var(ncid,varid,plst_tokamak_id)
  if (verbose_flag == 1) print *,err,plst_tag,plst_tokamak_id

  ! Number of thermal species
  plst_tag = 'dp1_nspec_th'
  err = nf90_inq_dimid(ncid,trim(plst_tag),varid)
  err = nf90_Inquire_Dimension(ncid,varid,len=plst_dp1_nspec_th)
  if (verbose_flag == 1)  print *,err,plst_tag,plst_dp1_nspec_th

  ! Number of thermal species (abridged)
  plst_tag = 'dp1_nspec_tha'
  err = nf90_inq_dimid(ncid,trim(plst_tag),varid)
  err = nf90_Inquire_Dimension(ncid,varid,len=plst_dp1_nspec_tha)
  if (verbose_flag == 1)  print *,err,plst_tag,plst_dp1_nspec_tha

  ! Number of species (thermal + fast)
  plst_tag = 'dp1_nspec_all'
  err = nf90_inq_dimid(ncid,trim(plst_tag),varid)
  err = nf90_Inquire_Dimension(ncid,varid,len=plst_dp1_nspec_all)
  if (verbose_flag == 1) print *,err,plst_tag,plst_dp1_nspec_all

  ! Number of species (thermal + fast, abridged)
  plst_tag = 'dp1_nspec_alla'
  err = nf90_inq_dimid(ncid,trim(plst_tag),varid)
  err = nf90_Inquire_Dimension(ncid,varid,len=plst_dp1_nspec_alla)
  if (verbose_flag == 1) print *,err,plst_tag,plst_dp1_nspec_alla

  ! Gridpoints
  plst_tag = 'dim_nrho'
  err = nf90_inq_dimid(ncid,trim(plst_tag),varid)
  err = nf90_Inquire_Dimension(ncid,varid,len=plst_dim_nrho)
  if (verbose_flag == 1) print *,err,plst_tag,plst_dim_nrho

  ! Gridpoints (geometry)
  plst_tag = 'dim_nrho_eq_geo'
  err = nf90_inq_dimid(ncid,trim(plst_tag),varid)
  err = nf90_Inquire_Dimension(ncid,varid,len=plst_dim_nrho_eq_geo)
  if (verbose_flag == 1) print *,err,plst_tag,plst_dim_nrho_eq_geo
  if (plst_dim_nrho_eq_geo /= plst_dim_nrho) then
     print *,'ERROR: plst_dim_nrho_eq_geo /= plst_dim_nrho'
     stop
  endif

  nx = plst_dim_nrho

  call allocate_internals
  call allocate_plasmastate_vars

  ! Species names (abridged)
  plst_tag = 'ALLA_name'
  err = nf90_inq_varid(ncid,trim(plst_tag),varid)
  err = nf90_get_var(ncid,varid,plst_alla_name(1:plst_dp1_nspec_alla))

  ! Species names
  plst_tag = 'ALL_name'
  err = nf90_inq_varid(ncid,trim(plst_tag),varid)
  err = nf90_get_var(ncid,varid,plst_all_name(1:plst_dp1_nspec_all))

  ! Flux-surface volume
  err = nf90_inq_varid(ncid,trim('vol'),varid)
  err = nf90_get_var(ncid,varid,plst_vol)

  ! Normalizing B for toroidal flux
  err = nf90_inq_varid(ncid,trim('B_axis_vac'),varid)
  err = nf90_get_var(ncid,varid,plst_b_axis_vac)

  ! B_phi orientation
  err = nf90_inq_varid(ncid,trim('kccw_Bphi'),varid)
  err = nf90_get_var(ncid,varid,plst_kccw_bphi)

  ! J_phi orientation
  err = nf90_inq_varid(ncid,trim('kccw_Jphi'),varid)
  err = nf90_get_var(ncid,varid,plst_kccw_jphi)

  ! Root of normalized toroidal flux (rho)
  err = nf90_inq_varid(ncid,trim('rho'),varid)
  err = nf90_get_var(ncid,varid,plst_rho(:))

  ! Grad(rho)
  err = nf90_inq_varid(ncid,trim('grho1'),varid)
  err = nf90_get_var(ncid,varid,plst_grho1(:))

  ! Toroidal flux
  err = nf90_inq_varid(ncid,trim('phit'),varid)
  err = nf90_get_var(ncid,varid,plst_phit(:))

  ! Poloidal flux/(2pi)
  err = nf90_inq_varid(ncid,trim('psipol'),varid)
  err = nf90_get_var(ncid,varid,plst_psipol(:))

  ! Elongation
  err = nf90_inq_varid(ncid,trim('elong'),varid)
  err = nf90_get_var(ncid,varid,plst_elong(:))

  ! Triangularity
  err = nf90_inq_varid(ncid,trim('triang'),varid)
  err = nf90_get_var(ncid,varid,plst_triang(:))

  ! 1/q
  err = nf90_inq_varid(ncid,trim('iota'),varid)
  err = nf90_get_var(ncid,varid,plst_iota(:))
  q(:) = 1.0/plst_iota(:)

  ! R_major (-pi)
  err = nf90_inq_varid(ncid,trim('R_midp_in'),varid)
  err = nf90_get_var(ncid,varid,plst_r_midp_in(:))

  ! R_major (0)
  err = nf90_inq_varid(ncid,trim('R_midp_out'),varid)
  err = nf90_get_var(ncid,varid,plst_r_midp_out(:))

  ! Elevation
  err = nf90_inq_varid(ncid,trim('Z_midp'),varid)
  err = nf90_get_var(ncid,varid,plst_z_midp(:))

  ! Z_eff
  err = nf90_inq_varid(ncid,trim('Zeff'),varid)
  err = nf90_get_var(ncid,varid,plst_zeff(1:nx-1))
  plst_zeff(nx) = 0.0

  ! Temperatures
  err = nf90_inq_varid(ncid,trim('Ts'),varid)
  err = nf90_get_var(ncid,varid,plst_ts(1:nx-1,1:plst_dp1_nspec_th)) 

  ! Temperature (Te at r/a=1)
  err = nf90_inq_varid(ncid,trim('Te_bdy'),varid)
  err = nf90_get_var(ncid,varid,dummy)

  plst_ts(nx,1) = dummy

  ! Temperature (Ti at r/a=1)
  err = nf90_inq_varid(ncid,trim('Ti_bdy'),varid)
  err = nf90_get_var(ncid,varid,dummy)

  plst_ts(nx,2:plst_dp1_nspec_th) = dummy

  ! Densities
  err = nf90_inq_varid(ncid,trim('ns'),varid)
  err = nf90_get_var(ncid,varid,plst_ns(1:nx-1,1:plst_dp1_nspec_th)) 

  ! Densities (n at r/a=1)
  err = nf90_inq_varid(ncid,trim('ns_bdy'),varid)
  err = nf90_get_var(ncid,varid,plst_ns(nx,1:plst_dp1_nspec_th)) 
  
  ! Bean density
  err = nf90_inq_varid(ncid,trim('nbeami'),varid)
  err = nf90_get_var(ncid,varid,plst_nb(1:nx-1))

  ! Beam perpendicular energy
  err = nf90_inq_varid(ncid,trim('eperp_beami'),varid)
  err = nf90_get_var(ncid,varid,plst_eperpb(1:nx-1))

  ! Beam parallel energy
  err = nf90_inq_varid(ncid,trim('epll_beami'),varid)
  err = nf90_get_var(ncid,varid,plst_eparb(1:nx-1))

  plst_ns(1:nx-1,plst_dp1_nspec_th+1) = plst_nb
  plst_ts(1:nx-1,plst_dp1_nspec_th+1) = 2.0/3.0*(plst_eparb + plst_eperpb)

  plst_ns(nx,plst_dp1_nspec_th+1) = plst_ns(nx-1,plst_dp1_nspec_th+1)
  plst_ts(nx,plst_dp1_nspec_th+1) = plst_ts(nx-1,plst_dp1_nspec_th+1)

  ! Total plasma pressure, thermal + fast ions
  err = nf90_inq_varid(ncid,trim('P_eq'),varid)
  err = nf90_get_var(ncid,varid,plst_ptowb(:))

  ! Radial electrostatic (equilibrium) potential (kV not keV!) 
  err = nf90_inq_varid(ncid,trim('Epot'),varid)
  err = nf90_get_var(ncid,varid,plst_epot(:))

  ! Rotation frequency (1/s) 
  err = nf90_inq_varid(ncid,trim('omegat'),varid)
  err = nf90_get_var(ncid,varid,plst_omegat(1:nx-1))

  ! Rotation frequency at r/a=1 (1/s) 
  err = nf90_inq_varid(ncid,trim('omegat_bdy'),varid)
  err = nf90_get_var(ncid,varid,dummy)

  plst_omegat(nx) = dummy

  ! Electron power
  err = nf90_inq_varid(ncid,trim('pbe'),varid)
  err = nf90_get_var(ncid,varid,plst_pbe(1:nx-1))
  plst_pbe(nx) = 0.0

  ! Ion power
  err = nf90_inq_varid(ncid,trim('pbi'),varid)
  err = nf90_get_var(ncid,varid,plst_pbi(1:nx-1))
  plst_pbi(nx) = 0.0

  ! Electron source power
  err = nf90_inq_varid(ncid,trim('pe_trans'),varid)
  err = nf90_get_var(ncid,varid,plst_pe_trans(1:nx-1))
  plst_pe_trans(nx) = 0.0

  ! Ion source power
  err = nf90_inq_varid(ncid,trim('pi_trans'),varid)
  err = nf90_get_var(ncid,varid,plst_pi_trans(1:nx-1))
  plst_pi_trans(nx) = 0.0

  ! Electron-ion Coupling
  err = nf90_inq_varid(ncid,trim('qie'),varid)
  err = nf90_get_var(ncid,varid,plst_pei_trans(1:nx-1))
  plst_pei_trans(nx) = 0.0

  ! Angular momentum source torque 
  err = nf90_inq_varid(ncid,trim('tq_trans'),varid)
  err = nf90_get_var(ncid,varid,plst_tq_trans(1:nx-1))
  plst_tq_trans(nx) = 0.0

  ! Particle source 
  err = nf90_inq_varid(ncid,trim('sn_trans'),varid)
  err = nf90_get_var(ncid,varid,plst_sn_trans(1:nx-1))
  plst_sn_trans(nx) = 0.0

  err = nf90_close(ncid)

  ! Ensure zero of flux and correct sign
  dpsi(:) = abs(plst_psipol(:)-plst_psipol(1))*(-plst_kccw_jphi)

  ! Compute rmin and rmaj based on outer and 
  ! inner major radii at midplane (of course, 
  ! this is not correct for elevated plasmas).
  !
  rmin(:)  = 0.5*(plst_r_midp_out-plst_r_midp_in)
  rmaj(:)  = 0.5*(plst_r_midp_out+plst_r_midp_in)

  zmag(:)  = plst_z_midp(:)
  kappa(:) = plst_elong(:)
  delta(:) = plst_triang(:)

  !----------------------------------------------
  ! Error check for missing/zero boundary (n,T)
  !
  do i=1,plst_dp1_nspec_th  
     call boundary_fix(plst_rho,plst_ts(:,i),nx)
     call boundary_fix(plst_rho,plst_ns(:,i),nx)
  enddo
  !----------------------------------------------

  ! No squareness 
  zeta(:) = 0.0

end subroutine prgen_read_plasmastate

subroutine boundary_fix(x,f,n)

  implicit none

  integer, intent(in) :: n
  real, intent(in), dimension(n) :: x
  real, intent(inout), dimension(n) :: f
  real, parameter :: edge_tol = 1e-6
  real :: dx
  real :: xlen

  if (f(n)/f(1) > edge_tol) return

  dx   = x(n-1)-x(n-2)
  xlen = log(f(n-2)/f(n-1))
  f(n) = f(n-1)*exp(-dx/xlen) 

end subroutine boundary_fix

