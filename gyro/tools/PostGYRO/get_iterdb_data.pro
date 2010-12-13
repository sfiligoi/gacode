FUNCTION get_iterdb_data, filename, DIRNAME = dirname, KINEFIT = kinefit
;
; C. Holland, UCSD
; v1.0: 5/11/2007
; v2.0: 8/30/2007: added dirname option
;
; opens a ascii iterdb file (from autoonetwo), reads in data, and
; returns a structure with some of the more useful fields
;

;  DEFAULT, dirname,  '~/research/BES/MachNumber06/iterdb_files/'
  DEFAULT, dirname,  '~/research/CECE/'
  DEFAULT, filename, 'iterdb.126285_01500'
  DEFAULT, kinefit, 0

  tmpstr = ' '
  OPENR, lun, dirname+'/' + filename, /GET_LUN
  READF, lun, tmpstr
  READF, lun, tmpstr
  READF, lun, shot_number
  READF, lun, tmpstr
  READF, lun, NR
  READF, lun, tmpstr
  READF, lun, n_ion
  READF, lun, tmpstr
  READF, lun, n_prime
  READF, lun, tmpstr
  READF, lun, n_impurity
  READF, lun, tmpstr
  READF, lun, n_neutral
  READF, lun, tmpstr
  READF, lun, beam_ion_idx

  FOR ii = 0, n_prime-1 DO BEGIN
      READF, lun, tmpstr        
      READF, lun, tmpstr
  ENDFOR

  FOR ii = 0, n_impurity-1 DO BEGIN
      READF, lun, tmpstr        
      READF, lun, tmpstr
  ENDFOR

  FOR ii = 0, n_neutral-1 DO BEGIN
      READF, lun, tmpstr        
      READF, lun, tmpstr
  ENDFOR

  READF, lun, tmpstr
  READF, lun, time
  READF, lun, tmpstr
  READF, lun, R_geo
  READF, lun, tmpstr
  READF, lun, R_mag
  READF, lun, tmpstr
  READF, lun, R0
  READF, lun, tmpstr
  READF, lun, kappa
  READF, lun, tmpstr
  READF, lun, delta
  READF, lun, tmpstr
  READF, lun, pindent
  READF, lun, tmpstr
  READF, lun, volo
  READF, lun, tmpstr
  READF, lun, xcareao
  READF, lun, tmpstr
  READF, lun, Btor     ;vaccum Btor at R_mag(?)
  READF, lun, tmpstr
  READF, lun, tot_j, tot_j_ohmic, tot_j_bootstrap, tot_j_beam, tot_j_RF
  READF, lun, tmpstr
  READF, lun, betap  ;polodial beta
  READF, lun, tmpstr
  READF, lun, beta   ;toroidal beta
  READF, lun, tmpstr
  READF, lun, ali    ;inductance
  READF, lun, tmpstr
  READF, lun, te0
  READF, lun, tmpstr
  READF, lun, ti0

  psi_on_rho = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, psi_on_rho

  rho = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, rho

  fcap = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, fcap

  gcap = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, gcap

  hcap = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, hcap

  Te = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, Te

  Ti = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, Ti

  q = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, q

  nelec = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, nelec

  ni_prim = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, ni_prim

  ni_imp = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, ni_imp

  ;ionization particle source
  sion = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, sion

  ;recombination particle source
  srecom = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, srecom

  ;cx thermal neutral particle source
  scx = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, scx

  ;cx beam neutral particle sink
  sbcx = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, sbcx

  ;total particle source
  source = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, source

  dudt = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, dudt

  ni_fast = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, ni_fast

  ;total neutral density
  n_neutral = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, n_neutral

  ;wall neutral density
  n_neutral_wall = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, n_neutral_wall

  ;volume neutral density
  n_neutral_volume = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, n_neutral_volume

  ;volume neutral source
  volume_neutral_src = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, volume_neutral_src

  ;beam electron source
  sbelec = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, sbelec

  ;beam thermal ion source
  sbion = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, sbion

  j_tot = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, j_tot

  j_ohmic = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, j_ohmic

  j_bootstrap = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, j_bootstrap

  j_beam = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, j_beam

  j_RF = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, j_RF

  Er_tor = FLTARR(NR)
  IF (kinefit EQ 0) THEN BEGIN
      READF, lun, tmpstr
      READF, lun, Er_tor
  ENDIF

  ugly_geo_factor = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, ugly_geo_factor

  Zeff = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, Zeff

  ang_rot = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, ang_rot

  chi_e = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, chi_e

  chi_i = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, chi_i

  kappa_i_neo = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, kappa_i_neo

  wdot_e = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, wdot_e

  wdot_i = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, wdot_i

  cond_elec = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, cond_elec

  cond_ion = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, cond_ion

  conv_elec = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, conv_elec

  conv_ion = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, conv_ion

  beam_pow_to_elec = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, beam_pow_to_elec

  ;elec-ion equilibration watts/meter**3
  qdelt = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, qdelt

  ;beam power to ions, watts/meter**3
  beam_pow_to_ions = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, beam_pow_to_ions

  ;RF electron heating, watts/meter**3
  qrfe = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, qrfe

  ;-RF ion heating, watts/meter**3
  qrfi = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, qrfi

  ;electron heating from recomb/ion
  qione = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, qione

  qioni = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, qioni

  qcx = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, qcx

  elec_heat_2d = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, elec_heat_2d

  ion_heat_2d = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, ion_heat_2d

  fusion_elec_heat = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, fusion_elec_heat

  fusion_ion_heat = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, fusion_ion_heat

  beam_fusion_elec_heat = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, beam_fusion_elec_heat

  beam_fusion_ion_heat = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, beam_fusion_ion_heat

  qmag_elec_heat = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, qmag_elec_heat

  sawtooth_elec_heat = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, sawtooth_elec_heat

  sawtooth_ion_heat = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, sawtooth_ion_heat

  rad_pow_dens = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, rad_pow_dens

  elec_ohm_pow_dens = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, elec_ohm_pow_dens

  fs_Rmaj = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, fs_Rmaj

  fs_rmin = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, fs_rmin

  fs_vol = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, fs_vol

  fs_kappa = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, fs_kappa

  fs_delta_plus = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, fs_delta_plus

  fs_delta_minus = FLTARR(NR)
  IF (kinefit EQ 1) THEN fs_delta_minus = fs_delta_plus ELSE BEGIN
      READF, lun, tmpstr
      READF, lun, fs_delta_minus
  ENDELSE

  fs_indentation = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, fs_indentation

  fs_surfarea = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, tmpstr
  READF, lun, fs_surfarea

  fs_crosssecarea = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, fs_crosssecarea

  fs_absgradrho = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, fs_absgradrho

  fs_gradrhosq = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, fs_gradrhosq

  READF, lun, tmpstr
  READF, lun, nbdry

  r_bdry = FLTARR(nbdry)
  READF, lun, tmpstr
  READF, lun, r_bdry

  z_bdry = FLTARR(nbdry)
  READF, lun, tmpstr
  READF, lun, z_bdry

  beam_torque_bdry = FLTARR(NR)
  READF, lun, tmpstr
  READF, lun, beam_torque_bdry

  FREE_LUN, lun

  delta0 = 0 ;(fs_delta_plus + fs_delta_minus)/2
  delta1 = 0 ;(fs_delta_plus - fs_delta_minus)/2
  results = {shot:shot_number, NR:NR, psi:psi_on_rho, rho:rho, rmin:fs_rmin, $
             Rmaj:fs_Rmaj, kappa:fs_kappa, delta0:delta0, delta1:delta1, q:q, $
             Te:Te, Ti:Ti, n_e:nelec, n_ion:ni_prim, n_fastion:ni_fast, $
             n_imp:ni_imp, Zeff:Zeff, chi_e:chi_e, chi_i:chi_i, $
             kappa_i_neo:kappa_i_neo, j_tot:j_tot, $
             ang_rot:ang_rot, Er_tor:Er_tor, fs_vol: fs_vol}

  RETURN, results
END ;get_iterdb_data
