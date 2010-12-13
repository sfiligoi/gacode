;---------------------------------------------------------
; harmonics_event.pro
;
; PURPOSE:
;  Extract poloidal harmonics from theta-dependent 
;  eigenfunction
; 
; REVISIONS
; 15 May 03: jc
;  Comments added, fixed normalization.
;--------------------------------------------------------

pro harmonics_event, harmonics

  common GLOBAL
  common PLOT_VARIABLES
  
  widget_control, harmonics.id, $
      get_uvalue=uvalue

  wset, harmonics_wid

  ;;-------------------------------------------------------
  ;; Now, begin menu selection
  ;;
  top:
  
  case (uvalue) of 
    
    1: begin

      str_m = strtrim(string(fix(m_c)),2)
      str_n = strtrim(string(fix(n_tor[in_c])),2)

      ;;-------------------------------------------------
      ;; widget control for label functions
      ;;
      widget_control, harmonics.top, $
          get_uvalue=state, $
          /no_copy
      
      t_str = 't = '+strtrim(string(t[t_c]),2)
      m_str = 'm = '+str_m
      n_str = 'n = '+str_n

      widget_control, state.t_label, set_value=t_str
      widget_control, state.m_label, set_value=m_str
      widget_control, state.n_label, set_value=n_str      

      widget_control, harmonics.top, $
          set_uvalue=state, $
          /no_copy
      ;;-------------------------------------------------

      plot_title = 'harmonic_'+str_m+'_'+str_n

      plot_def_new,plot_title

      f_m    = complexarr(n_r)
      r_surf = fltarr(2)

      get_harmonic,f_m,r_surf

      xtitle = '!3r/a'
      if (plot_mode eq 2) then begin
        title  = '!3m='+str_m+', n='+str_n
        ytitle = '!4u!3!dmn!n'
      endif else begin
        title  = '!4u!3(m='+str_m+',n='+str_n+')!3'
        ytitle = '!3aqua=Re ; yellow=Im'
      endelse

      plot,[0],[0],$
          /nodata,$
          title=title,$
          xstyle=1,$
          xminor=0,$
          xrange=[min(r),max(r)],$
          xtitle=xtitle,$
          ystyle=1,$
          yminor=0,$
          yrange=max(abs(f_m))*[-1,1],$
          ytitle=ytitle,$
          color=line

      ;; Plot harmonic

      oplot,r,float(f_m),color=color_vec[0]
      oplot,r,imaginary(f_m),color=color_vec[1]

      ;; Plot singular surface(s)

      oplot,r_surf[0]*[1,1],[-10^4,10^4],linestyle=1,color=line
      oplot,r_surf[1]*[1,1],[-10^4,10^4],linestyle=1,color=line

      if (plot_export eq 1) then begin
        openw,1,plot_title+'.out'
        printf,1,r,float(f_m),imaginary(f_m)
        close,1
        print,'Exported data to ',plot_title+'.out'
        plot_export = 0
      endif

      plot_finish

    end   
    
    2: begin
      counter_up,t_c,n_time1,1
      goto, plot_it
    end

    3: begin
      counter_dn,t_c,0,1
      goto, plot_it
    end

    4: begin
      counter_up,t_c,n_time1,10 
      goto, plot_it
    end

    5: begin 
      counter_dn,t_c,0,10
      goto, plot_it
    end

    6: begin
      counter_up,in_c,n_n-1,1 
      goto, plot_it
    end

    7: begin
      counter_dn,in_c,0,1  
      goto, plot_it
    end

    8: begin
      counter_up,m_c,1000,1 
      goto, plot_it
    end

    9: begin
      counter_dn,m_c,-1000,1
      goto, plot_it
    end

    10: begin
      if (n_field gt 1) then i_f = 1-i_f
      goto, plot_it
    end

    11: begin
      plot_mode = 2
      goto, plot_it
    end

    12: begin
      plot_export = 1
      goto, plot_it
    end

    13: widget_control, harmonics.top, /destroy

  endcase
  
  return

  plot_it:
  
  uvalue=1
  goto,top

end
