pro initialize

  common GLOBAL
  common PLOT_VARIABLES

  ;;---------------------------------------------
  ;; Color table defaults (see color_table_setup)
  ;; for details.
  ;; 
  nlevels = 24
  ctab    = 5
  ;;---------------------------------------------

  ;;---------------------------------------------
  ;; Default rotation angles
  ;;
  ax_c = 40
  az_c = 30
  ;;---------------------------------------------
  
  ;;---------------------------------------------
  ;;---------------------------------------------

  ;;---------------------------------------------
  ;; character size
  ;;
  c_size      = 1.6
  !p.charsize = c_size
  dotsize     = 1.5
  annosize    = 0.8
  c_scale     = 1.0
  ;;
  ;;---------------------------------------------

  ;;---------------------------------------------
  ;; Set initial simulation directory 
  ;; to 'null' (traps error).
  
  simdir = 'null'  

  t_c   = 0
  in_c  = 0
  m_c   = 0
  p_squ = 0
  ik_c  = 0
  i_c   = 0 
  j_c   = 0 
  i_p   = 0
  ;;---------------------------------------------
  
  ;;---------------------------------------------
  ;; Miscellaneous initializations:
  ;;
  ;; plot_mode: (1=screen,2=PS)
  plot_mode = 1
  ;;
  plot_export = 0
  
  ;; i_units: (0=gB, 1=m^2/s)
  i_units = 0

  i_zero   = 0
  i_abs    = 0
  i_r_dr   = 0
  i_log    = 0
  i_loglog = 0
  i_zero   = 0
  i_spec   = 0 
  i_ptype  = 0
  i_moment = 1
  i_f      = 0
  i_all    = 0

  nbin     = 10

  c_table_min = -1.0
  c_table_max = -1.0
  ;;---------------------------------------------

end
