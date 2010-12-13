pro all_profiles_event, all_profiles_obj

  common GLOBAL
  common PLOT_VARIABLES
  common PROFILE_SIM_DATA
  common PRIVATE_ALL_PROFILES

  widget_control, all_profiles_obj.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  if uvalue ge 40 then begin 

     case (uvalue) of

        40: begin
           i_spec = i_spec+1
           if i_spec gt n_spec-1 then i_spec = 0
           goto, plot_it
        end

        41: begin
           i_r_dr = 1-i_r_dr
           goto, plot_it
        end

        42: begin
           plot_export = 1
           goto, plot_it
        end

        43: begin
           plot_mode = 2
           goto, plot_it
        end

        44: widget_control, all_profiles_obj.top, /destroy

     endcase

     return

  endif else begin

     i_all_profile = uvalue

  endelse


  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  pname = profile_label[i_all_profile-1]
  if (i_all_profile eq 5) or $
    (i_all_profile eq 6) or $
    (i_all_profile eq 7) or $
    (i_all_profile eq 8) or $
    (i_all_profile eq 22) then begin
     pname = spec_label[i_spec]+pname
  endif

  plot_def_new,pname

  y = fltarr(n_r)

  case (i_all_profile) of 

     1: y[*] = r[*]
     2: y[*] = q_i[*]
     3: y[*] = r_m[*]
     4: y[*] = q_m[*]
     5: y[*] = dlntdr_s[i_spec,*]
     6: y[*] = dlnndr_s[i_spec,*]
     7: y[*] = tem_s[i_spec,*]
     8: y[*] = den_s[i_spec,*]
     9: y[*] = phi_doppler[*]
     10: y[*] = aspect_s[*]
     11: y[*] = delta_s[*]
     12: y[*] = zeta_s[*]
     13: y[*] = kappa_s[*]
     14: y[*] = shift_s[*]
     15: y[*] = shat_s[*]
     16: y[*] = s_delta_s[*]
     17: y[*] = s_zeta_s[*]
     18: y[*] = s_kappa_s[*]
     19: y[*] = zmag_s[*]
     20: y[*] = dzmag_s[*]
     21: y[*] = beta_unit_s[*]
     22: y[*] = pgamma_s[i_spec,*]
     23: y[*] = b_unit_s[*]
     24: y[*] = dr_eodr[*]
     25: y[*] = z_eff_s[*]
     26: y[*] = nu_s[*]
     27: y[*] = gamma_eb_s[*]
     28: y[*] = w0_s[*]

     31: begin
        y[*] = den_s[n_spec-1,*]
        for i=0,n_ion-1 do begin
           y[*] = y[*]-zcharge[i]*den_s[i,*]
        endfor
     end

  endcase

  if (i_r_dr eq 1) then begin
     y_prim = fltarr(n_r)
     scalar_deriv,r,y,y_prim,boundary_method
     y = -y_prim
     pname = pname+'-prime'
  endif

  dy = 0.1*(max(y)-min(y))

  if dy lt 0.01 then dy = 0.1

  plot,[0],[0],$
    /nodata,$
    ytitle=pname,$
    xstyle=1,$
    xminor=0,$
    xrange=[min(r),max(r)],$
    xtitle='!3r/a',$
    ystyle=1,$
    yminor=0,$
    yrange=[min(y)-dy,max(y)+dy],$
    color=line,$
    charsize=csize

  oplot,r,y,color=color_vec[1]
  oplot,r,y,psym=8,symsize=dotsize

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for i=0,n_r-1 do begin
        printf,1,r[i],y[i]
     endfor
     print,'Exported data to ',pname+'.idlout'
     close,1
     plot_export = 0
  endif
  ;;---------------------------------------------------

  plot_finish

  return

end
