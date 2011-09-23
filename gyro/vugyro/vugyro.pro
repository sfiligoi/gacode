;---------------------------------------------------------
; vugyro.pro
;
; PURPOSE:
;  Toplevel IDL file for vugyro suite: an IDL widget-based 
;  visualization utility for the gyro code.
;--------------------------------------------------------

pro vugyro, _EXTRA=extra

  common GLOBAL,$
    ave_poloidal, $
    ax_c, $
    az_c, $
    boundary_method, $
    box_multiplier, $
    chi_e_exp, $
    chi_i_exp, $
    chi_string, $
    collision_wid, $
    csa_string,$
    c_size, $
    CTab, $
    c_table_max,$
    c_table_min,$
    den_s, $
    diff, $
    diff_e, $
    diff_i, $
    diff_i_trapped, $
    diff_n, $
    diff_e_n, $
    diff_QL_n, $
    phi_squared_QL_n, $
    g_squared_QL_n, $
    diff_ne_exp, $
    diff_to_flow_e1, $
    diff_to_flow_e2, $
    diff_to_flow_heating, $
    diff_to_flow_mi, $ 
    diff_to_flow_ne, $ 
    diff_trapped, $
    dir, $
    dist_contour_wid, $
    dn_ss, $
    electron_method, $
    energy,$
    entropy, $
    eta_i_diff_exp, $
    eta_i_tot_exp, $
    exists_balloon, $
    exists_coll, $
    exists_diff, $
    exists_diff_i, $
    exists_diff_i_t, $
    exists_diff_n, $
    exists_diff_QL_n, $
    exists_diff_t, $
    exists_entropy,$
    exists_exp_profile, $
    exists_exp_derived, $
    exists_field_r0, $
    exists_gbflux, $
    exists_gbflux_n, $
    exists_gbflux_i, $
    exists_gbflux_trapped, $
    exists_gbflux_i_trapped, $
    exists_g_squared_QL_n, $
    exists_geometry, $
    exists_h, $
    exists_he, $
    exists_kxkyspec, $
    exists_mom_e, $
    exists_mom_n, $
    exists_mom_v, $
    exists_nl_transfer, $
    exists_omega, $
    exists_phi_squared_QL_n, $
    exists_profile, $
    exists_source, $
    exists_time, $
    exists_transport_ne_te_ti, $
    exists_u, $
    exists_units, $
    exists_velocity, $
    exists_zerobar, $
    exp_profile, $
    exp_derived, $
    exp_profile_indx, $
    exp_profile_label, $ 
    field_r0, $
    field_r0_grid, $
    field_r0_wid, $
    field_r_wid, $
    f_nek, $
    gbflux, $
    gbflux_i, $
    gbflux_trapped, $
    gbflux_i_trapped, $
    gbflux_n, $
    geometry, $
    geometry_label, $
    harmonics_wid, $
    hp, $
    hpe, $
    ht, $
    hte, $
    i_abs, $
    i_all, $ 
    i_c, $
    i_f, $
    ik_c, $
    i_log, $
    i_loglog, $
    i_moment, $               
    in_c, $
    i_p, $
    i_ptype, $
    i_r_dr,$
    i_spec, $             
    it1, $
    it2, $
    i_units,$
    i_zero, $
    j_c, $
    kr_rho, $
    kr_rho_string,$
    kt_rho, $
    kt_rho_string,$
    kxkyspec, $
    lambda,$
    lambda_tp, $
    m_c, $
    mom_e,$
    mom_n,$
    mom_v,$
    n0, $
    nbin, $
    n_bnd, $
    n_CTab, $
    ndir, $
    n_dn, $
    n_energy, $
    n_exp_profile, $
    n_exp_derived, $
    n_field, $
    n_fine, $
    n_geometry, $
    n_ion, $
    n_kinetic, $
    n_lambda, $
    nLevels, $
    n_moment, $
    n_n, $
    n_nLevels, $
    nonlinear_flag, $
    nonlinear_transfer_wid, $
    nl_transfer, $
    n_pass, $
    n_profile_label, $
    n_r, $
    n_rho, $
    n_spec, $
    n_ss, $
    n_stack,$
    n_theta_p,$
    n_theta_plot, $
    n_theta_section,$
    n_theta_t,$
    n_time, $
    n_time1, $
    n_tor, $ 
    n_trap, $
    n_trho, $
    p_moment, $
    phi, $
    phi_doppler, $
    phi_star, $
    plot_export, $
    plot_mode, $
    plot_units, $
    poloidal_wid, $
    profile_label, $
    ps_color, $
    ps_label, $
    p_squ, $
    ps_thick, $
    ps_val, $
    q_i, $
    quiet, $
    r, $
    remotedir,$
    remotedir_flag,$
    r_from_rho, $
    rho_hat, $
    rho_s, $
    s_CTab, $
    simdir, $
    simroot, $
    skip_large, $
    s_nLevels, $
    source, $
    spec_label,$
    sx, $       
    sy, $
    t, $
    t_c, $
    tem_s, $
    theta, $
    thetad_plot, $
    thetad_r0_plot, $
    theta_plot, $
    theta_t_p,$
    theta_t_t,$
    t_max, $
    t_min, $
    transport_ne_te_ti, $
    t_string, $
    u, $
    v_CTab, $
    version, $
    version_date, $
    version_tag, $
    v_nLevels, $
    w,$
    xunits, $
    xi_p,$
    xi_t,$
    xsize, $
    ysize, $
    zcharge, $
    zerobar

  
;  Reduce name conflicts for those who use IDL in non-GYRO contexts:
;  need to have GYRO come before <IDL_DEFAULT> on LOKI
!PATH=expand_path("$GACODE_ROOT/gyro/vugyro:$GACODE_ROOT/tgyro/tools/idl:<IDL_DEFAULT>")
;  while vugyro is run, this over-rides the IDL_PATH set in the user's
;  dot-files, which may include many directories with experimental
;  ananlysis routines. <IDL_DEFAULT> is recommended for use with 
;  expand_path, and corresponds to the envrionment variable IDL_DIR.

;;----------------------------------;
  ;; Set display characteristics
  ;;
  device,get_visual_depth=deep
  ;;  print,deep
  device,true_color=24,decomposed=0
  ;;-----------------------------------;

  ;;---------------------------------------;
  ;; Filled circles for node points,;
  ;;
  a = findgen(8)*(!Pi*2/8.0)
  usersym,0.8*cos(a),0.8*sin(a),/fill
  ;;---------------------------------------;

  spawn,getenv("GACODE_ROOT")+'/shared/bin/'+'gacode_getversion'

  ;;---------------------------------------;
  ;; Some control parameters; 
  ;;
  user_control, _EXTRA=extra
  ;;
  home = getenv('GACODE_ROOT')+'/gyro/vugyro/'

  if (remotedir_flag eq 0) then begin
     simroot = getenv('PWD')
  endif else begin
     simroot = remotedir
  endelse
  ;;
  dirtext = simroot
  ;;
  version = 'null'
  openr,1,getenv("GACODE_ROOT")+'/.VERSION',error=i_err
  if (i_err eq 0) then begin
     readf,1,version
  endif else begin
     print,'No version information. Setting to 1.0.0'
     version = '1.0.0'
  endelse  
  close,1
  ;;---------------------------------------;

  term_type = getenv('TERM')
  print,'Terminal type: ',term_type

  ;;--------------------------------------------
  ;; COMMON BLOCK INITIALIZATIONS
  ;;
  balloon_common
  midplane_common
  tag_common
  plot_variables_common
  poloidal_common
  zmoment_common
  t_error_common
  collision_common
  profile_sim_common
  color_table_setup
  ;;--------------------------------------------

  initialize

  ;;-----------------------------------;
  ;; Get list of simulation directories
  ;; 
  get_dir_list
  ;;-----------------------------------;

  ;;-------------------------------------------------------;
  ;; Build toplevel widget:
  ;;   
  mainwindow = widget_base(title='vugyro '+version+dirtext,$
                           mbar=menubar,$
                           tlb_frame_attr=1)

  base = widget_base(mainwindow,$
                     /column)

  ;;----------------------------------------------------------
  ;; OPTIONS
  ;;
  top = widget_button(menubar,$
                      value='Options',$
                      /menu)

  x = widget_button(top,$
                    value='MASTER options',$
                    uvalue='master_see')

  x = widget_button(top,$
                    value='Postscript options',$
                    uvalue='ps_options_see')

  x = widget_button(top,$
                    value='Contour options',$
                    uvalue='contour_options_see')

  x = widget_button(top,$
                    value='Exit vuGyro',$
                    uvalue='exit')
  ;;----------------------------------------------------------
  

  ;;----------------------------------------------------------
  ;; SIMULATION selector:
  ;;
  top = widget_button(menubar,$
                      value='Data',$
                      /menu)

  x = widget_button(top,$
                    value='Choose via dialog...',$
                    uvalue='sim_pick_dialog')
  
  for i=0,ndir-1 do begin
     x = widget_button(top,$
                       value=dir(i),$
                       uvalue=strtrim(string(i),2))  
  endfor
  ;;----------------------------------------------------------
  

  ;;----------------------------------------------------------
  ;; EQUILIBRIUM
  ;;
  top = widget_button(menubar,$
                      value='Equilibrium',$
                      /menu)

  x = widget_button(top,$
                    value='Raw profiles',$
                    uvalue='exp_profiles_see')

  x = widget_button(top,$
                    value='Simulation profiles',$
                    uvalue='all_profiles_see')

  x = widget_button(top,$
                    value='Geometry arrays',$
                    uvalue='geometry_see')

  x = widget_button(top,value='Collision stencil',$
                    uvalue='collision_see')

  ;;----------------------------------------------------------


  ;;------------------------------------------------------------
  ;; FLUCTUATIONS
  ;;
  top = widget_button(menubar,$
                      value='Fluctuations',$
                      /menu)

  x = widget_button(top,$
                    value='Linear frequency',$
                    /menu)

  x1 = widget_button(x,$
                     value='frequency vs. time',$
                     uvalue='linear_freq_see')

  x1 = widget_button(x,$
                     value='frequency vs. ky',$
                     uvalue='linear_spec_see')
  
  x = widget_button(top,$
                    value='Fields',$
                    /menu)

  x1 = widget_button(x,$
                     value='ballooning space',$
                     uvalue='balloon_see')

  x1 = widget_button(x,$
                     value='theta dependence',$
                     uvalue='field_theta_see')

  x1 = widget_button(x,$
                     value='radial dependence',$
                     uvalue='field_r_see')

  x1 = widget_button(x,$
                     value='poloidal harmonics',$
                     uvalue='harmonics_see')

  x1 = widget_button(x,$
                     value='radial shear',$
                     uvalue='field_shear_see')

  x1 = widget_button(x,$
                     value='theta dependence at r=r0',$
                     uvalue='field_r0_theta_see')

  x1 = widget_button(x,$
                     value='time trace',$
                     uvalue='time_trace_see')

  x = widget_button(top,$
                    value='Fluxes',$
                    /menu)

  x1 = widget_button(x,$
                     value='flux vs. time',$
                     uvalue='gbflux_see')

  x1 = widget_button(x,$
                     value='flux vs. radius',$
                     uvalue='gbflux_i_see')

  x1 = widget_button(x,$
                     value='flux vs. ky',$
                     uvalue='gbflux_n_see')

  x1 = widget_button(x,$
                     value='velocity dependence',$
                     uvalue='velocity_see')

  x = widget_button(top,$
                    value='Diffusion - box average',$
                    /menu)

  x1 = widget_button(x,$
                     value='particle/energy + exp (box-car ave)',$
                     uvalue='diffusion_bca_see')

  x1 = widget_button(x,$
                     value='particle/energy + exp (stat ave)',$
                     uvalue='diffusion_ave_see')

  x1 = widget_button(x,$
                     value='n-dependence and QL test',$
                     uvalue='each_diffusion_see')

  x = widget_button(top,$
                    value='Diffusion - radial profile',$
                    /menu)

  x1 = widget_button(x,$
                     value='particle/energy',$
                     uvalue='diffusion_i_ave_see')

  x1 = widget_button(x,$
                     value='effective energy',$
                     uvalue='diffusion_i_eff_see')

  x1 = widget_button(x,$
                     value='particle/energy in (r,t)',$
                     uvalue='diffusion_i_rt_see')

  x = widget_button(top,$
                    value='Entropy production',$
                    uvalue='entropy_balance_see')

  x = widget_button(top,value='Spectral Analysis',/menu)

  x1 = widget_button(x,$
                     value='2D (n,p)(t) spectrum',$
                     uvalue='spectrum_np_see')

  x1 = widget_button(x,$
                     value='2D <(n,p)> spectrum',$
                     uvalue='spectrum_np_ave_see')

  x1 = widget_button(x,$
                     value='1D f(n)(t)',$
                     uvalue='spectrum_n_see')

  x1 = widget_button(x,$
                     value='1D f(n,p)(t)',$
                     uvalue='each_np_see')

  x1 = widget_button(x,$
                     value='1D <f(n)>',$
                     uvalue='spectrum_n_ave_see')

  x = widget_button(top,$
                    value='Midplane data',$
                    /menu)

  x1 = widget_button(x,$
                     value='1D Total power spectra',$
                     uvalue='midplane_power_see')

  x1 = widget_button(x,$
                     value='radial power spectra per n',$
                     uvalue='midplane_n_power_see')

  x1 = widget_button(x,$
                     value='averaged radial shearing',$
                     uvalue='midplane_deriv_see')

  x1 = widget_button(x,$
                     value='temporal spectrum per n',$
                     uvalue='midplane_freq_see')

  x1 = widget_button(x,$
                     value='2D Power spectrum (new)',$
                     uvalue='midplane_avgkspect_see')

  x1 = widget_button(x,$
                     value='2D Power spectrum (old)',$
                     uvalue='midplane_xypower_see')

  x1 = widget_button(x,$
                     value='1D Power spectrum',$
                     uvalue='midplane_each_power_see')

  x1 = widget_button(x,$
                     value='2D fluctuations',$
                     uvalue='midplane_fluc_see')

  x1 = widget_button(x,$
                     value='RMS fluctuation levels',$
                     uvalue='midplane_rms_flucamp_see')

  x1 = widget_button(x,$
                     value='correlation functions',$
                     uvalue='midplane_corr_see')

  x1 = widget_button(x,$
                     value='n=0 (kr,omega) spectrum',$
                     uvalue='midplane_zfspect_see')

  x = widget_button(top,$
                    value='Profile evolution',$
                    /menu)

  x1 = widget_button(x,$
                     value='profile fluctuations',$
                     uvalue='osc_see')

  x1 = widget_button(x,$
                     value='profile fluctuations (from source)',$
                     uvalue='osc_h_see')

  x1 = widget_button(x,$
                     value='potential fluctuations',$
                     uvalue='zerobar_see')

  x1 = widget_button(x,$
                     value='source shape',$
                     uvalue='source_see')

  x1 = widget_button(x,$
                     value='source t-dependence',$
                     uvalue='source_t_see')

  x1 = widget_button(x,$
                     value='fluctuations in r-t plane',$
                     uvalue='zmoment_see')

  x = widget_button(top,$
                    value='Poloidal contours',$
                    uvalue='poloidal_see')

  x = widget_button(top,$
                    value='Time integration error',$
                    uvalue='t_error_see')
  ;;------------------------------------------------------------
  
  ;;------------------------------------------------------------
  ;; Special
  ;;
  top = widget_button(menubar,$
                      value='Extras',$
                      /menu)

  x = widget_button(top,$
                    value='Distribution function',$
                    /menu)

  x1 = widget_button(x,$
                     value='theta-dependence',$
                     uvalue='dist_see')

  x1 = widget_button(x,$
                     value='(k,theta)-dependence',$
                     uvalue='dist_k_see')

  x = widget_button(top,$
                    value='Nonlinear transfer',$
                    /menu)

  x1 =  widget_button(x,$
                     value='kx-ky transfer spectrum',$
                     uvalue='nonlinear_transfer_see')

  x = widget_button(top,$
                    value='Transport profiles',$
                    /menu)

  x1 = widget_button(x,$
                     value='transport_ne_te_ti gradients',$
                     uvalue='transport_ne_te_ti_gradients_see')

  x1 = widget_button(x,$
                     value='transport_ne_te_ti flows',$
                     uvalue='transport_ne_te_ti_flows_see')

  x1 = widget_button(x,$
                     value='transport_ne_te_ti profiles',$
                     uvalue='transport_ne_te_ti_profiles_see')

  x1 = widget_button(x,$
                     value='transport_ne_te_ti flows_vs_grads',$
                     uvalue='transport_ne_te_ti_flows_vs_grads_see')
  
  ;;------------------------------------------------------------

  ;;------------------------------------------------------------
  ;; HELP
  ;;
  top = widget_button(menubar,$
                      value='Browse',$
                      /help, $
                      /menu)

  x = widget_button(top,$
                    value='Contents of out.gyro.run',$
                    uvalue='view_run')

  x = widget_button(top,$
                    value='Contents of efficiency.out',$
                    uvalue='view_efficiency')
  ;;------------------------------------------------------------


  ;; Set a fixed-size draw widget
  drawlogo = widget_draw(base,xsize=500,ysize=1,retain=2)

  ;; Help messages
  info = widget_text(base,ysize=4,font='8x13',/wrap,value=$
;   '123456789012345678901234567890123456789012345678901234567890'$
    'NOTES AND REMINDERS:                                        '$
   +'                                                            '$
   +'If you experience an error reading data, run                '$
   +'gyro -t <simdir> to update your simulation directory.      ')

  ;; Render the top-level widget
  widget_control, mainwindow, /hourglass, /realize

  ;; Set line-plotting colors:

  set_line_colors

  ;; Use X-windows cursor:

  device,/cursor_original

  xmanager,'vugyro',mainwindow,event_handler='gyro_choices'

end
