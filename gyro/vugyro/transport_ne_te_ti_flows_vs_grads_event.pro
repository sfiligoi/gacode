pro transport_ne_te_ti_flows_vs_grads_event, transport_ne_te_ti_flows_vs_grads 

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_TR_FLOWS_VS_GRADS
  
  widget_control, transport_ne_te_ti_flows_vs_grads.id, $
                  get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; Now, begin menu selection

  case (uvalue) of 
     
     1: begin
        i_transport_ne_te_ti = 1
        j_transport_ne_te_ti = 4
        title  = 'Electron particle flow'
        pname  = 'transport_ne_flows_vs_grads'
        ytitle = '!3Electron particle flow!d!n'
        goto, plot_it
     end   

     2: begin
        i_transport_ne_te_ti =2 
        j_transport_ne_te_ti =5
        title = 'Electron power flow'
        ytitle = '!3Electron Power!d!n'
        pname = 'transport_te_flows_vs_grads'
        goto, plot_it
     end   

     3: begin
        i_transport_ne_te_ti =3
        j_transport_ne_te_ti =6
        title = 'Ion power flow'
        ytitle = '!3Ion Power!d!n'
        pname = 'transport_ti_flows_vs_grads'
        goto, plot_it
     end

     4: begin
        i_rho_offset = i_rho_offset-1 
        goto, plot_it
     end

     5: begin
        i_rho_offset = i_rho_offset+1
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

     13: widget_control, transport_ne_te_ti_flows_vs_grads.top, /destroy
     
  endcase
  
  ;;----------------------------------------------------------;

  return

  plot_it:

  plot_def_new,pname

  exp_flows = fltarr(n_trho)
  exp_grads = fltarr(n_trho)
  exp_grads0 = fltarr(n_trho)
  z_tr = fltarr(n_trho)
  
  i_rho_pivot = 28
  i_rho_left = 13
  i_rho_right = 37
  ii = i_rho_pivot+i_rho_offset


  if(ii gt i_rho_right) then begin
     ii = i_rho_right
  endif
  if(ii lt i_rho_left) then begin
     ii = i_rho_left
  endif

  print, 'ii=',ii,'i_rho_offset=',i_rho_offset,'r=',exp_profile[1,ii]/exp_profile[1,n_trho]

  if (i_transport_ne_te_ti eq 1) then begin
     z_tr[*] = transport_ne_te_ti[4,*,0]
     exp_flows[ii] = exp_profile[15,ii+1]

     exp_grads[ii] = (exp_profile[7,ii]-exp_profile[7,ii+1])  $
                     *2./(exp_profile[7,ii]+exp_profile[7,ii+1])  $
                     /(exp_profile[1,ii+1]-exp_profile[1,ii])

     exp_grads0[ii] = (z_tr[ii-1]-z_tr[ii])  $
                  *2./(z_tr[ii-1]+z_tr[ii])  $
                  /(exp_profile[1,ii+1]-exp_profile[1,ii]) 


  endif
  if (i_transport_ne_te_ti eq 2) then begin
     z_tr[*] = transport_ne_te_ti[5,*,0]
     exp_flows[ii] = exp_profile[11,ii+1]

     exp_grads[ii] = (exp_profile[6,ii]-exp_profile[6,ii+1])  $
                     *2./(exp_profile[6,ii]+exp_profile[6,ii+1])  $
                     /(exp_profile[1,ii+1]-exp_profile[1,ii])
 
     exp_grads0[ii] = (z_tr[ii-1]-z_tr[ii])  $
                  *2./(z_tr[ii-1]+z_tr[ii])  $
                  /(exp_profile[1,ii+1]-exp_profile[1,ii])

  endif
  if (i_transport_ne_te_ti eq 3) then begin
     z_tr[*] = transport_ne_te_ti[6,*,0]
     exp_flows[ii] = exp_profile[12,ii+1]

     exp_grads[ii] = (exp_profile[25,ii]-exp_profile[25,ii+1])  $
                     *2./(exp_profile[25,ii]+exp_profile[25,ii+1])  $
                     /(exp_profile[1,ii+1]-exp_profile[1,ii])

     exp_grads0[ii] = (z_tr[ii-1]-z_tr[ii])  $
                  *2./(z_tr[ii-1]+z_tr[ii])  $
                  /(exp_profile[1,ii+1]-exp_profile[1,ii])


  endif

  x_plot = fltarr(n_time)
  x_plot[*] =t[*]

  yf_ave = fltarr(n_time)
  yf_ave[*] = 0.0
  for iia = i_rho_left,i_rho_right do begin
     yf_ave[*] =  yf_ave[*]+transport_ne_te_ti[i_transport_ne_te_ti,iia,*]
  endfor
  yf_ave[*] =  yf_ave[*]/(i_rho_right-i_rho_left+1)
  
  yf_plot = fltarr(n_time)
  yf_plot[*] = transport_ne_te_ti[i_transport_ne_te_ti,ii,*]/exp_flows[ii]

  yg_plot = fltarr(n_time)
  yg_plot[*] = (transport_ne_te_ti[j_transport_ne_te_ti,ii-1,*] $
                -transport_ne_te_ti[j_transport_ne_te_ti,ii,*]) $
               *2./(transport_ne_te_ti[j_transport_ne_te_ti,ii-1,*] $
                    +transport_ne_te_ti[j_transport_ne_te_ti,ii,*]) $
               /(exp_profile[1,ii+1]-exp_profile[1,ii])/exp_grads0[ii]
  u_plot = fltarr(n_time)
  u_plot[*] = 1.

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
       ytitle=ytitle,$
       color=line

  oplot,x_plot,yf_plot,color=color_vec[0]
  oplot,x_plot,yg_plot,color=color_vec[1]
  oplot,x_plot,u_plot,color=color_vec[2]


  plot_finish

  ;;----------------------------------------------------------------

  return

end
