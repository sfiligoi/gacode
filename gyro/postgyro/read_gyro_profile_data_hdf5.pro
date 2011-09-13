FUNCTION read_gyro_profile_data_HDF5, simdir, profile_data
;
; C. Holland, UCSD
;
; v1.0: 9.13.2011: load in profile, units data from out.gyro.initdata.h5

  dirpath = GETENV('GYRO_DIR') + '/sim/' + simdir
  filepath =  dirpath + '/out.gyro.initdata.h5'  
  data = H5_PARSE(filepath,/READ_DATA)

  n_r = data.n_x._data
  n_theta_section = data.n_theta_section._data
  n_pass = data.n_pass._data
  n_trap  = data.n_trap._data
  n_energy = data.n_energy._data
  n_theta_plot = data.n_theta_plot._data
  n_0 = data.n0._data
  n_n = data.n_n._data
  n_dn = data.d_n._data
   n_field = data.n_field._data
  n_ion = data.n_ion._data
  n_kinetic = data.n_kinetic._data
  n_spec = data.n_spec._data

  ;TO BE FIXED
  n_bnd = 0  ;data.n_explicit_damp._data
  n_fine = n_theta_plot ;data.n_fine._data
  theta_mult = 1

  r = data.r._data
  q = data.q._data
  r_s = data.r_s._data
  q_s = data.q_s._data
  dlntdr_s = data.dlntdr_s._data
  dlnndr_s = data.dlnndr_s._data
  tem_s = data.tem_s._data
  den_s = data.den_s._data
  aspect_s = data.aspect_s._data
  delta_s = data.delta_s._data
  zeta_s = data.zeta_s._data
  kappa_s = data.kappa_s._data
  shift_s = data.drmaj_s._data
  shat_s = data.shat_s._data
  s_delta_s = data.s_delta_s._data
  s_zeta_s = data.s_zeta_s._data
  s_kappa_s = data.s_kappa_s._data
  zmag_s = data.zmag_s._data
  dzmag_s = data.dzmag_s._data
  beta_unit_s = data.beta_unit_s._data
  gamma_e_s = data.gamma_e_s._data
  gamma_p_s = data.gamma_p_s._data
  mach_s = data.mach_s._data
  b_unit_s = data.b_unit_s._data
  dr_eodr = data.dr_eodr._data
  z_eff_s = data.z_eff_s._data
  nu_s = data.nu_s._data
  w0_s = data.w0_s._data

  lambda = data.lambda._data
  energy = data.energy._data
  lambda_tp = data.lambda_tp._data
  kt_rho = data.krho_collect._data
  rho_s = data.rho_sd._data
  zcharge = data.zcharge._data

  ;;units
  Aphys = data.a._data
  csda = data.csda_norm._data  ;csda_norm_d
  Gamma_gB = data.gamma_gbd._data
  Q_gB = data.Q_gBd._data
  Pi_gB = data.Pi_gBd._data
  S_gB = data.S_gBd._data

  exists_nu_geo = 1
  nu_geo = data.nu._data

  ;return useful data in structure, easy to add to as needed
  profile_data = {n_r:n_r, $    ;# of radial grid points
                  n_theta_plot:n_theta_plot, $ ;# of _saved_ theta points
                  n_n:n_n, $    ;n = n0 + n_dn*INDGEN(n_n)
                  n_0:n_0, $
                  n_dn:n_dn, $
                  n_field:n_field, $ ;1 for phi only, 2 for phi and A_||
                  n_ion:n_ion, $ ;# ions evolved
                  n_spec:n_spec, $ ;total # ions and e-
                  n_kinetic:n_kinetic, $ ;# of evolved ions and e-
                  ktheta: kt_rho, $ ;k_theta rho_s/i as function of n
                  rho_s: rho_s, $     ;rho_star

                  ;Geometry variables
                  R0:aspect_s*r_s, $ ;R_0(r)/a
                  r:r, $        ;r/a
                  q:q, $        ;q(r) (flux-surface avg?)
                  shear:shat_s, $ ;r/q dq/dr (fs avg or midplane?)
                  kappa:kappa_s, $ ;kappa(r)
                  s_kappa:s_kappa_s, $ ;s_kappa = r dln(kappa)/dr
                  delta:delta_s, $
                  s_delta:s_delta_s, $
                  zeta:zeta_s, $ ;
                  s_zeta:s_zeta_s, $ ;
                  zmag:zmag_s, $ ;
                  dzmag:dzmag_s, $ ;
                  theta_mult:theta_mult, $
                  nu_geo:nu_geo, $
                  exists_nu_geo:exists_nu_geo, $

                  ;mean profiles
                  n:den_s, $    ;density profiles
                  dlnndr: dlnndr_s, $ ;a/Ln
                  T:tem_s, $    ;temperature profiles
                  dlnTdr: dlnTdr_s, $ ;a/LT
                  Z:zcharge, $  ;Z(species)
                  B_unit:B_unit_s, $ ;effective_B
                  gamma_e: gamma_e_s, $ ;ExB shear
                  gamma_p: gamma_p_s, $ ;parallel shear
                  mach: mach_s, $ ;Mach #
                  w0:w0_s, $    ;rotation frequency

                  Aphys: Aphys*100,  $ ;cm
                  csda: csda*1e-3,      $ ;kHz
                  Gamma_gB: Gamma_gB*624., $ ;10**19/m**2/s
                  Q_gB: Q_gB*100., $ ;W/cm**2
                  Pi_gB: Pi_gB, $ ;Nm/m**2
                  S_gB: S_gB, $ ;W/cm**3
                  
                  n_bnd: n_bnd, $ ;# pts in boundary layer
                  nonlinear_flag: FIX(data.nonlinear_flag._data) $
                 }

  PRINT, 'Read ' + filepath

  RETURN, 1  ;if we got here then loaded ok

END ;read_gyro_profile_data_HDF5
