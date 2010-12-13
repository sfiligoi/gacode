pro zerobar_event, zerobar_obj

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_ZEROBAR

  widget_control, zerobar_obj.id, $
                  get_uvalue=uvalue

  wset, widget

  ;;------------------------------------------------------------------
  ;; MENU FUNCTION EVENT HANDLERS
  ;;
  case (uvalue) of 
     
     1: begin
        active_zerobar = 1
        title  = '!3Flux-surface, t-average of !4u!3'
        ytitle = '!3<!4u!3>'
        pname  = 'phi_fluxave'
        goto, plot_it
     end   
     
     2: begin
        active_zerobar = 2
        title  = 'Flux-surface, t-average of A'
        ytitle = '!3<A>'
        pname = 'a_fluxave'
        goto, plot_it
     end

     3: begin
        active_zerobar = 3
        title  = 'Flux-surface, t-average of d(!4u!3)/dr'
        ytitle = '!3<d(!4u!3)/dr>'
        pname  = 'dphi_fluxave'
        goto, plot_it
     end   

     4: begin
        active_zerobar = 4
        title  = 'Flux-surface, t-average of d!u2!n(!4u!3)/dr!u2!n'
        ytitle = '!3<d!u2!n!4u!3/dr!u2!n>'
        pname  = 'ddphi_fluxave'
        goto, plot_it
     end   

     5: begin
        zoom = zoom*2
        goto, plot_it
     end

     6: begin
        zoom = zoom/2
        goto, plot_it
     end

     7: begin
        plot_mode = 2
        goto, plot_it
     end

     8: widget_control, zerobar_obj.top, /destroy
     
  endcase
  ;;
  ;;-------------------------------------------------
  
  return

  plot_it:

  plot_def_new,pname

  y  = fltarr(n_r,n_time)
  dy = fltarr(n_r)

  phi1 = fltarr(n_r)
  phi2 = fltarr(n_r)

  scalar_deriv,r,phi_doppler,phi1,boundary_method
  scalar_deriv,r,phi1,phi2,boundary_method

  case (active_zerobar) of

     1: begin
        y[*,*] = zerobar[0,*,*]
        dy[*]  = phi_doppler[*]
     end

     2: begin
        y[*,*] = zerobar[1,*,*]
        dy[*]  = 0.0
     end

     3: begin
        fr = fltarr(n_r,n_time)
        fr[*,*] = zerobar[0,*,*]
        vector_deriv,r,fr,y,boundary_method
        dy[*] = phi1[*]
     end
     4: begin
        fr = fltarr(n_r,n_time)
        fr[*,*] = zerobar[0,*,*]
        vector_deriv,r,fr,y,boundary_method
        vector_deriv,r,y,fr,boundary_method
        y[*,*] = rho_s*fr[*,*]
        dy[*]  = rho_s*phi2[*]
     end

  endcase

  y0 = fltarr(n_time)
  y1 = fltarr(n_r)

  for ii=0,n_r-1 do begin
     y0[*] = y[ii,*]
     diff_stat_fast,y0,it1,it2,A_f
     y1[ii] = A_f
  endfor

  print,sqrt(total(y1[*]*y1[*])/n_r)

  plot,[0],[0],$
       /nodata,$
       title=title,$
       xstyle=1,$
       xminor=0,$
       xrange=[min(r),max(r)],$
       xtitle='!3r/a',$
       ystyle=1,$
       yminor=0,$
       yrange=[-1,1]*max(abs(y1+dy))/zoom,$
       color=line

  ;; Equilibrium

  oplot,r,dy,color=color_vec[2],psym=8,symsize=0.5

  ;; Equilibrium plus fluctuations

  oplot,r,dy+y1,color=line

  plot_bnd
  plot_finish

  return

end
