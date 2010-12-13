pro osc_h_event, osc_h_obj

  common GLOBAL
  common PROFILE_SIM_DATA
  common PLOT_VARIABLES
  common PRIVATE_OSC_H
  
  widget_control, osc_h_obj.id, $
    get_uvalue=uvalue

  wset,widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------
  
  case (uvalue) of 
     
     1: goto, plot_it

     2: begin
        gradient_flag = 1-gradient_flag
        goto, plot_it
     end   

     3: begin
        i_osc_mom = 1-i_osc_mom
        goto, plot_it
     end   

     5: begin
        zoom = zoom/2
        goto, plot_it
     end

     6: begin
        zoom = zoom*2
        goto, plot_it
     end

     7: begin
        equil_flag = 1-equil_flag
        goto, plot_it
     end

     8: begin
        plot_mode = 2
        goto, plot_it
     end

     9: widget_control, osc_h_obj.top, /destroy
     
  endcase
  
  ;;----------------------------------------------------------;

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  indx = gradient_flag+2*i_osc_mom
  ;; 0=density
  ;; 1=density'
  ;; 2=energy
  ;; 3=energy'

  y    = fltarr(n_r,n_time)
  y_eq = fltarr(n_r)

  h0_n = fltarr(n_r,n_time)
  h0_e = fltarr(n_r,n_time)

  ; density moment
  h0_n[*,*] = source[0,0,*,*]

  ; energy moment
  h0_e[*,*] = source[0,1,*,*]

  case (indx) of 

     0: begin
        
        title  = '!3Ion gyrocenter density'
        ytitle = '!3<n!di!n>'
        pname  = 'h-ni'

        y_eq[*] = equil_flag*den_s[0,*]
        for i=0,n_r-1 do begin
           y[i,*] = y_eq[i]+h0_n[i,*]
        endfor

     end

     1: begin

        title  = '!3Ion gyrocenter density log-grad'
        ytitle = '!3<dln(n!di!n)/dr>'
        pname  = 'h-dlnnidr'

        y_eq[*] = dlnndr_s[0,*]

        yx = fltarr(n_r,n_time)
        yx[*,*] = h0_n[*,*]

        vector_deriv,r,yx,y,boundary_method

        for i=0,n_r-1 do begin
           y[i,*] = -y[i,*]+den_s[0,i]*y_eq[i]
           y[i,*] = y[i,*]/(den_s[0,i]+pert_flag*yx[i])
        endfor

     end

     2: begin

        title = '!3Ion gyrocenter energy'
        ytitle = '<E!di!n>'
        pname = 'h-ti'

        y_eq[*] = equil_flag*tem_s[0,*]

        for i=0,n_r-1 do begin
           y[i,*] = (h0_e[i,*]-1.5*h0_n[i,*]*tem_s[0,i])/ $
             (1.5*den_s[0,i])
           y[i,*] = y_eq[i]+y[i,*]
        endfor

     end

     3: begin

        title = '!3Ion gyrocenter energy log-grad'
        ytitle = '<dln(E!di!n)/dr>'
        pname = 'h-dlnEidr'
        
        y_eq[*] = dlntdr_s[0,*]

        yx = fltarr(n_r,n_time)
        for i=0,n_r-1 do begin
           yx[i,*] = (h0_e[i,*]-1.5*h0_n[i,*]*tem_s[0,i])/ $
             (1.5*den_s[0,i])
        endfor

        vector_deriv,r,yx,y,boundary_method

        print,pert_flag

        for i=0,n_r-1 do begin
           y[i,*] = -y[i,*]+tem_s[0,i]*y_eq[i]
           y[i,*] = y[i,*]/(tem_s[0,i]+pert_flag*yx[i])
        endfor

     end

  endcase

  plot_def_new,pname

  y_mid = y_eq[n_r/2]

  y0    = fltarr(n_time)
  ave_y = fltarr(n_r)

  for ii=0,n_r-1 do begin
     y0[*] = y[ii,*]
     diff_stat_fast,y0,it1,it2,ave
     ave_y[ii] = ave
  endfor
  
  xtitle = '!3r/a with '+t_string

  plot,[0],[0],$
    /nodata,$
    title=title,$
    xstyle=1,$
    xminor=0,$
    xrange=[min(r),max(r)],$
    xtitle=xtitle,$
    ystyle=1,$
    yminor=0,$
    yrange=y_mid+[-1,1]*zoom,$
    ytitle=ytitle,$
    color=line

  ;; Equilibrium

  oplot,r,y_eq,color=color_vec[2],psym=8,symsize=0.5 

  ;; Equilibrium plus fluctuations

  oplot,r,ave_y,color=line
  
  plot_bnd 
  plot_finish
  
  return

end
