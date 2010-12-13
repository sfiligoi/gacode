pro exp_profiles_event, exp_profiles_obj

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_EXP_PROFILES
  
  widget_control, exp_profiles_obj.id, $
                  get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  if uvalue ge 40 then begin 

     case (uvalue) of

        40: begin
           i_r_dr = 1-i_r_dr
           goto, plot_it
        end

        41: begin
           i_r_rho = 1-i_r_rho
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

        44: widget_control, exp_profiles_obj.top, /destroy

     endcase

     return

  endif else begin

     i_profile = uvalue

  endelse

  plot_it:
  
  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  pname = exp_profile_label[i_profile]

  plot_def_new,pname

  x = fltarr(n_rho)
  y = fltarr(n_rho)
  
  bound = fltarr(2)
  bound[0] = min(r)
  bound[1] = max(r)

  if (i_r_rho eq 0) then begin
     x[*] = r_from_rho[*]
     xtitle='!3r'
  endif else begin
     x[*] = rho_hat[*]
     xtitle = '!4q!3'
     bound = interpol(rho_hat,r_from_rho,bound)
  endelse

  y[*] = exp_profile[exp_profile_indx[i_profile],*]

  if (i_r_dr eq 1) then begin
     y_prim = fltarr(n_rho)
     scalar_deriv,x,y,y_prim,2
     y = -y_prim
     pname = pname+'-prime'
  endif

  plot,[0],[0],$
       /nodata,$
       ytitle=pname,$
       xstyle=1,$
       xminor=0,$
       xrange=[0,1],$
       xtitle=xtitle,$
       ystyle=1,$
       yminor=0,$
       yrange=[min(y)-0.1,max(y)+0.1],$
       color=line,$
       charsize=csize

  oplot,x,y,color=color_vec[0]
  ;;oplot,x,y,psym=8,symsize=dotsize

  oplot,bound[0]*[1,1],100*[-1,1],linestyle=1
  oplot,bound[1]*[1,1],100*[-1,1],linestyle=1

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for i=0,n_rho-1 do begin
        printf,1,x[i],y[i]
     endfor
     print,'Exported data to ',pname+'.idlout'
     close,1
     plot_export = 0
  endif
  ;;---------------------------------------------------

  plot_finish
  
  return

end
