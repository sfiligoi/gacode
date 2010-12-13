pro transport_ne_te_ti_flows_event, transport_ne_te_ti_flows 
;
; CH 6-19-08- minor bug fixes, added time-avg plot option
;

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_TR_FLOWS
  
  widget_control, transport_ne_te_ti_flows.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; Now, begin menu selection
  
  case (uvalue) of 
     
     1: begin
        i_transport_ne_te_ti = 1
        title  = 'Electron particle flow'
        pname  = 'transport_ne_flows'
        ytitle = '!3Electron particle flow!d!n'
        goto, plot_it
     end   

     2: begin
        i_transport_ne_te_ti =2 
        title = 'Electron power flow'
        ytitle = '!3Electron Power!d!n'
        pname = 'transport_te_flows'
        goto, plot_it
     end   

     3: begin
        i_transport_ne_te_ti =3
        title = 'Ion power flow'
        ytitle = '!3Ion Power!d!n'
        pname = 'transport_ti_flows'
        goto, plot_it
     end

     4: begin
        counter_up,t_c,n_time1,1
        goto, plot_it
     end

     5: begin
        counter_dn,t_c,0,1
        goto, plot_it
     end

     6: begin
        counter_up,t_c,n_time1,10
        goto, plot_it
     end

     7: begin
        counter_dn,t_c,0,10
        goto, plot_it
     end

     8: begin
        i_flow_avg_flag = 0
        i_all = 1-i_all 
        goto, plot_it
     end

     9: begin
        i_flow_avg_flag = 1 - i_flow_avg_flag
        i_all = 0
        goto, plot_it
     end

     14: begin
        zoom = zoom/2
        goto, plot_it
     end
     15: begin
        zoom = zoom*2
        goto, plot_it
     end

     11: begin
        plot_mode = 2
        goto, plot_it
     end

     12: begin
        i_ratio_value = 1-i_ratio_value
        goto, plot_it
     end

     13: widget_control, transport_ne_te_ti_flows.top, /destroy
     
  endcase
  
  ;;----------------------------------------------------------;

  return

  plot_it:

  bound = fltarr(2)
  bound[0] = min(r)
  bound[1] = max(r)

  plot_def_new,pname

  exp_value = fltarr(n_trho)

  if (i_transport_ne_te_ti eq 1) then begin
     for ii = 0,n_trho-1 do begin
 ;;        exp_value[ii] = exp_profile[19,ii+1]
        exp_value[ii] = exp_profile[15,ii+1]
     endfor
  endif
  if (i_transport_ne_te_ti eq 2) then begin
     for ii = 0,n_trho-1 do begin
 ;;        exp_value[ii] = exp_profile[15,ii+1]
        exp_value[ii] = exp_profile[11,ii+1]
     endfor
  endif
  if (i_transport_ne_te_ti eq 3) then begin
     for ii = 0,n_trho-1 do begin
 ;;       exp_value[ii] = exp_profile[16,ii+1]
        exp_value[ii] = exp_profile[12,ii+1]
     endfor
  endif


  
  ratio_value = fltarr(n_trho)
  
  if (i_ratio_value eq 0) then begin
     for ii = 0,n_trho-1 do begin
        ratio_value[ii] = 1.
     endfor
     units = ' [units of MW]'
     if (i_transport_ne_te_ti eq 1) then begin
        units = ' [units of MW/Kev]'
     endif
  endif else begin
     units = ' [ratio tr/exp]'
     if (i_transport_ne_te_ti eq 1) then begin
        for ii = 0,n_trho-1 do begin
           ratio_value[ii] = exp_profile[19,ii+1]
        endfor
     endif
     if (i_transport_ne_te_ti eq 2) then begin
        for ii = 0,n_trho-1 do begin
           ratio_value[ii] = exp_profile[15,ii+1]
        endfor
     endif
     if (i_transport_ne_te_ti eq 3) then begin
        for ii = 0,n_trho-1 do begin
           ratio_value[ii] = exp_profile[16,ii+1]
        endfor
     endif
  endelse

  x_plot = fltarr(n_trho)
  x_plot[*] = transport_ne_te_ti[0,*,t_c]
  y_plot = fltarr(n_trho)
  y_exp_plot = fltarr(n_trho)

  if (i_flow_avg_flag eq 1) then begin
      y_plot = TOTAL(transport_ne_te_ti[i_transport_ne_te_ti,*,it1:it2],3)/$
               (it2-it1+1)
      y_plot = REFORM(y_plot)/ratio_value
      y_exp_plot = exp_value/ratio_value
      xtitle = '!3r/a over' + t_string

     plot,[0],[0],$
          /nodata,$
          title=title,$
          xstyle=1,$
          xrange=[min(x_plot),max(x_plot)],$
          xtitle=xtitle,$
          ystyle=1,$
          yrange=[0,zoom],$
          ytitle=ytitle+units,$
          color=line

     oplot,x_plot,y_exp_plot,color=color_vec[0]
     oplot,x_plot,y_plot,color=color_vec[1]

  endif else if (i_all eq 0) then begin
     y_plot[*] = transport_ne_te_ti[i_transport_ne_te_ti,*,t_c]


     for ii = 0,n_trho-1 do begin
        y_plot[ii] = y_plot[ii]/ratio_value[ii]
        y_exp_plot[ii] = exp_value[ii]/ratio_value[ii]
     endfor


     t_str = ' t = '+strtrim(string(t[t_c]),2)
     xtitle = '!3r/a at '+t_str

     zoom0 = 0.0
     if (i_transport_ne_te_ti eq 1 ) then begin
        zoom0 = -zoom
     endif

     plot,[0],[0],$
       /nodata,$
       title=title,$
       xstyle=1,$
       xrange=[min(x_plot),max(x_plot)],$
       xtitle=xtitle,$
       ystyle=1,$
       yrange=[zoom0,zoom],$
       ytitle=ytitle+units,$
       color=line

     oplot,x_plot,y_exp_plot,color=color_vec[0]
     oplot,x_plot,y_plot,color=color_vec[1]

  endif

  if ((i_flow_avg_flag eq 0) and ( i_all eq 1)) then begin
     
     it_max = 0
     for it = 0,n_time-1,100 do begin
        it_max = it_max + 1
     endfor 

     y_plot_all = fltarr(n_trho,it_max)
     iit = 0
     for it = 0,n_time-1,100 do begin
        y_plot_all[*,iit] = transport_ne_te_ti[i_transport_ne_te_ti,*,it]
        iit = iit +1
     endfor


     for ii = 0,n_trho-1 do begin
        y_plot_all[ii,*] = y_plot_all[ii,*]/ratio_value[ii]
        y_exp_plot[ii] = exp_value[ii]/ratio_value[ii]
     endfor

     zoom0 = 0.0
     if (i_transport_ne_te_ti eq 1) then begin
        zoom0 = -zoom
     endif

     xtitle = '!3r/a '

     plot,[0],[0],$
       /nodata,$
       title=title,$
       xstyle=1,$
       xrange=[min(x_plot),max(x_plot)],$
       xtitle=xtitle,$
       ystyle=1,$
       yrange=[zoom0,zoom],$
       ytitle=ytitle+units,$
       color=line

     oplot,x_plot,y_exp_plot,thick=8,color=color_vec[0]
     
     iit = 0
     for it = 0,n_time-1,100 do begin
        y_plot[*]= y_plot_all[*,iit] 
        if(it ne 0 ) then begin
           oplot,x_plot,y_plot,thick=1,color=color_vec[1]
        endif
        iit = iit +1
     endfor

     oplot,x_plot,y_plot,thick=8,color=color_vec[2]

  endif

  oplot,bound[0]*[1,1],100*[-1,1],linestyle=1
  oplot,bound[1]*[1,1],100*[-1,1],linestyle=1


  plot_finish

  ;;----------------------------------------------------------------

  return

end
