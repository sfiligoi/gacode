pro map_exp_profiles

  common GLOBAL

  exp_profile = transpose(exp_profile)

  exp_profile_label = strarr(n_exp_profile)
  exp_profile_indx  = intarr(n_exp_profile)
  
  for i=0,17 do begin
     exp_profile_indx[i] = i
  endfor

  exp_profile_label[0] = 'rhogrid_exp'
  exp_profile_label[1] = 'rmin_exp'
  exp_profile_label[2] = 'rmaj_exp'
  exp_profile_label[3] = 'q_exp'
  exp_profile_label[4] = 'kappa_exp'
  exp_profile_label[5] = 'delta_exp'
  exp_profile_label[6] = 'te_exp'
  exp_profile_label[7] = 'ne_exp'
  exp_profile_label[8] = 'z_eff_exp'
  exp_profile_label[9] = 'w0_exp'
  exp_profile_label[10] = 'flow_mom_exp'
  exp_profile_label[11] = 'pow_e_exp'
  exp_profile_label[12] = 'pow_i_exp'
  exp_profile_label[13] = 'pow_ei_exp'
  exp_profile_label[14] = 'zeta_exp'
  exp_profile_label[15] = 'flow_beam_exp'
  exp_profile_label[16] = 'flow_wall_exp'
  exp_profile_label[17] = 'zmag_exp'

  ;; Restack ion density and temperature:
  
  for i=1,n_ion do begin

     j = 17+4*i-3
     
     exp_profile_label[j] = 'ni'+strtrim(string(i),1)+'_exp'
     exp_profile_indx[j] = 20+i-1

     exp_profile_label[j+1] = 'ti'+strtrim(string(i),1)+'_exp'
     exp_profile_indx[j+1] = 25+i-1

     exp_profile_label[j+2] = 'vphi_i'+strtrim(string(i),1)+'_exp'
     exp_profile_indx[j+2] = 30+i-1

     exp_profile_label[j+3] = 'vpol_i'+strtrim(string(i),1)+'_exp'
     exp_profile_indx[j+3] = 35+i-1

  endfor
  
  n_exp_profile = 18+4*n_ion

  ;; Generate r over entire range of rho.

  r_from_rho = exp_profile[1,*]/exp_profile[1,n_rho-1]
  rho_hat = exp_profile[0,*]/exp_profile[0,n_rho-1]

  ;;

  return

end
