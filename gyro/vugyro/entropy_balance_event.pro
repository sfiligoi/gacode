pro entropy_balance_event, entropy_balance_obj

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_ENTROPY_BALANCE
  
  widget_control, entropy_balance_obj.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: goto, plot_it

     2: begin
        i_spec = i_spec+1
        if (i_spec ge n_kinetic) then i_spec = 0
        goto, plot_it
     end   

     3: begin
        entropy_balance_max = entropy_balance_max*2
        goto, plot_it
     end

     4: begin
        entropy_balance_max = entropy_balance_max/2.0
        goto, plot_it
     end

     5: begin
        plot_mode = 2
        goto, plot_it
     end

     6: begin
        plot_export = 1
        goto, plot_it
     end

     7: widget_control, entropy_balance_obj.top, /destroy
     
  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  if (i_spec lt n_ion) and (electron_method ne 3) then begin
     tag_s  = 'ion'+strtrim(string(i_spec+1),2)
  endif else begin
     tag_s  = 'elec'
  endelse

  pname = tag_s+'_entropy'
  title = tag_s+'!3 Entropy Balance'

  plot_def_new,pname

  ytitle = '!3 sigma - H!u*!n dH/dt + D = 0'

  plot,[0],[0],$
    /nodata,$
    title=title,$
    ytitle=ytitle,$
    xstyle=1,$
    xrange=[0.0,max(t)],$
    xtitle=csa_string,$
    ystyle=1,$
    yrange=entropy_balance_max*[-1,1],$
    color=line


  oplot,t,entropy[i_spec,0,*],color=color_vec[0]
  oplot,t,entropy[i_spec,1,*],color=color_vec[1]
  oplot,t,entropy[i_spec,2,*],color=color_vec[2]
  oplot,t,entropy[i_spec,3,*],color=color_vec[3]
  oplot,t,entropy[i_spec,4,*],color=color_vec[4]

  entropy_tot = fltarr(n_time)

  entropy_tot[*] = entropy[i_spec,0,*]+$
                   entropy[i_spec,1,*]+$
                   entropy[i_spec,2,*]+$
                   entropy[i_spec,3,*]+$
                   entropy[i_spec,4,*]

  oplot,t,entropy_tot

  ;; sigma

  diff_stat_fast,entropy[i_spec,0,*],it1,it2,ave_y

  xyouts,0.55*max(t),$
    0.9*entropy_balance_max,'!3sigma: '+strtrim(string(ave_y),2),$
    color=color_vec[0]
 
  ;; - H* dH/dt

  diff_stat_fast,entropy[i_spec,1,*],it1,it2,ave_y

  xyouts,0.55*max(t),$
    0.8*entropy_balance_max,'!3-H!u*!n dH/dt: '+strtrim(string(ave_y),2),$
    color=color_vec[1]

  ;; D_coll

  diff_stat_fast,entropy[i_spec,2,*],it1,it2,ave_y

  xyouts,0.55*max(t),$
    0.7*entropy_balance_max,'!3D!dcoll!n: '+strtrim(string(ave_y),2),$
    color=color_vec[2]

  ;; D_tau

  diff_stat_fast,entropy[i_spec,3,*],it1,it2,ave_y

  xyouts,0.55*max(t),$
    0.6*entropy_balance_max,'!3D!dtau!n: '+strtrim(string(ave_y),2),$
    color=color_vec[3]

  ;; D_r

  diff_stat_fast,entropy[i_spec,4,*],it1,it2,ave_y

  xyouts,0.55*max(t),$
    0.5*entropy_balance_max,'!3D!dr!n: '+strtrim(string(ave_y),2),$
    color=color_vec[4]

  ;; Total (should vanish)

  diff_stat_fast,entropy_tot[*],it1,it2,ave_y

  xyouts,0.55*max(t),$
    0.4*entropy_balance_max,'!3Total: '+strtrim(string(ave_y),2),$
    color=line

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     printf,1,t,entropy_tot
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  endif
  ;;---------------------------------------------------

  plot_finish

  return

end
