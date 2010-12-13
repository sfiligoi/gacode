pro transport_ne_te_ti_gradients_event, transport_ne_te_ti_gradients 
;
; CH 6-19-08- minor bug fixes, added time-avg plot option
;

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_GRADS
  
  widget_control, transport_ne_te_ti_gradients.id, $
                  get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; Now, begin menu selection

  
  case (uvalue) of 
     
     1: begin
        i_transport_ne_te_ti = 4 
        title  = 'Electron density gradient'
        pname  = 'transport_ne_gradients'
        goto, plot_it
     end   

     2: begin
        i_transport_ne_te_ti = 5
        title = 'Electron temperature gradient'
        pname = 'transport_te_gradients'
        goto, plot_it
     end   
     
     3: begin
        i_transport_ne_te_ti= 6
        title = 'Ion temperature gradient'
        pname = 'transport_ti_gradients'
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

     9: begin
         i_grad_avg_flag = 1 - i_grad_avg_flag
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

     13: widget_control, transport_ne_te_ti_gradients.top, /destroy
     
  endcase
  
  ;;----------------------------------------------------------;

  return

  plot_it:

  bound = fltarr(2)
  bound[0] = min(r)
  bound[1] = max(r)

  plot_def_new,pname

  exp_value = fltarr(n_trho)
  exp_value0 = fltarr(n_trho)
  z_tr = fltarr(n_trho)

  if (i_transport_ne_te_ti eq 4) then begin
     z_tr[*] = transport_ne_te_ti[i_transport_ne_te_ti,*,0]
     ;;for ii = 0,n_trho-1 do begin
     for ii = 1,n_trho-1 do begin

        exp_value[ii] = (exp_profile[7,ii]-exp_profile[7,ii+1])  $
                        *2./(exp_profile[7,ii]+exp_profile[7,ii+1])  $
                        /(exp_profile[1,ii+1]-exp_profile[1,ii]) $
                        *exp_profile[1,n_trho-1]

        
        exp_value0[ii] = (z_tr[ii-1]-z_tr[ii])  $
                  *2./(z_tr[ii-1]+z_tr[ii])  $
                  /(exp_profile[1,ii+1]-exp_profile[1,ii]) $
                  *exp_profile[1,n_trho-1]

     endfor
  endif
  if (i_transport_ne_te_ti eq 5) then begin
     z_tr[*] = transport_ne_te_ti[i_transport_ne_te_ti,*,0]

     ;;for ii = 0,n_trho-1 do begin
       for ii = 1,n_trho-1 do begin
        exp_value[ii] = (exp_profile[6,ii]-exp_profile[6,ii+1])  $
                        *2./(exp_profile[6,ii]+exp_profile[6,ii+1])  $
                        /(exp_profile[1,ii+1]-exp_profile[1,ii]) $
                        *exp_profile[1,n_trho-1]

        exp_value0[ii] = (z_tr[ii-1]-z_tr[ii])  $
                  *2./(z_tr[ii-1]+z_tr[ii])  $
                  /(exp_profile[1,ii+1]-exp_profile[1,ii]) $
                  *exp_profile[1,n_trho-1]

     endfor
  endif
  if (i_transport_ne_te_ti eq 6) then begin
      z_tr[*] = transport_ne_te_ti[i_transport_ne_te_ti,*,0]

     ;;for ii = 0,n_trho-1 do begin
       for ii = 1,n_trho-1 do begin
 
        exp_value[ii] = (exp_profile[25,ii]-exp_profile[25,ii+1])  $
                        *2./(exp_profile[25,ii]+exp_profile[25,ii+1])  $
                        /(exp_profile[1,ii+1]-exp_profile[1,ii]) $
                        *exp_profile[1,n_trho-1]

        exp_value0[ii] = (z_tr[ii-1]-z_tr[ii])  $
                  *2./(z_tr[ii-1]+z_tr[ii])  $
                  /(exp_profile[1,ii+1]-exp_profile[1,ii]) $
                  *exp_profile[1,n_trho-1]



     endfor
  endif
  
  ratio_value = fltarr(n_trho)
  
  if (i_ratio_value eq 0) then begin
     for ii = 0,n_trho-1 do begin
        ratio_value[ii] = 1.
        units = 'tr/exp'
     endfor
     if (i_transport_ne_te_ti eq 4) then begin
        ytitle = '!3dlnNedr!d!n'
     endif
     if (i_transport_ne_te_ti eq 5) then begin
        ytitle = '!3dlnTedr!d!n'
     endif
     if (i_transport_ne_te_ti eq 6) then begin
        ytitle = '!3dlnTidr!d!n'
     endif
  endif
  if (i_ratio_value eq 1) then begin
     ratio_value[*] = exp_profile0[*]
     if (i_transport_ne_te_ti eq 4) then begin
        ytitle = '!3dlnNedr!d!n'
     endif
     if (i_transport_ne_te_ti eq 5) then begin
        ytitle = '!3dlnTedr!d!n'
     endif
     if (i_transport_ne_te_ti eq 6) then begin
        ytitle = '!3dlnTidr!d!n'
     endif
     units = ' [units 1 ]'
  endif
  
  print, 't=',t[t_c]

  x_plot = fltarr(n_trho)
  x_plot[*] = transport_ne_te_ti[0,*,t_c]
  y_plot = fltarr(n_trho)
  y_exp_plot = fltarr(n_trho)
  y_tr = fltarr(n_trho)

;  y_tr[*] = transport_ne_te_ti[i_transport_ne_te_ti,*,t_c]
;  for ii = 1,n_trho-1 do begin
;     y_plot[ii] = (y_tr[ii-1]-y_tr[ii])  $
;                  *2./(y_tr[ii-1]+y_tr[ii])  $
;                  /(exp_profile[1,ii+1]-exp_profile[1,ii]) $
;                  *exp_profile[1,n_trho-1]
;  endfor
;  y_plot[0] = 0.

  if (i_grad_avg_flag eq 0) then begin
      itmin = t_c
      itmax = t_c
      n_avg = 1
      t_str = ' t = '+strtrim(string(t[t_c]),2)
      xtitle = '!3r/a at '+t_str
  endif else begin
      itmin = it1
      itmax = it2
      n_avg = it2-it1+1
      xtitle = '!3r/a over' + t_string
  endelse

  for iit = itmin, itmax do begin
      y_tr[*] = transport_ne_te_ti[i_transport_ne_te_ti,*,iit]
      for ii = 1,n_trho-1 do begin
          y_plot[ii] += (y_tr[ii-1]-y_tr[ii])  $
                       *2./(y_tr[ii-1]+y_tr[ii])  $
                       /(exp_profile[1,ii+1]-exp_profile[1,ii]) $
                       *exp_profile[1,n_trho-1]
      endfor
  endfor
  y_plot /= n_avg
  y_plot[0] = 0.


  for ii = 0,n_trho-1 do begin
     y_plot[ii] = y_plot[ii]/ratio_value[ii]
     y_exp_plot[ii] = exp_value0[ii]/ratio_value[ii]
  endfor



  plot,[0],[0],$
       /nodata,$
       title=title,$
       xstyle=1,$
       xrange=[min(x_plot),max(x_plot)],$
       xtitle=xtitle,$
       ystyle=1,$
       yrange=[0,zoom],$
       ytitle=ytitle,$
       color=line

  oplot,x_plot,y_exp_plot,color=color_vec[0]
  oplot,x_plot,y_plot,color=color_vec[1]

  oplot,bound[0]*[1,1],100*[-1,1],linestyle=1
  oplot,bound[1]*[1,1],100*[-1,1],linestyle=1


  plot_finish

  ;;----------------------------------------------------------------

  return

end
