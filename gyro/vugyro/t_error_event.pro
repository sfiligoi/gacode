pro t_error_event, t_error_obj

  common GLOBAL
  common PLOT_VARIABLES
  common T_ERROR_DATA
  
; In auto mode  t_error_obj won't be defined
if(n_elements(t_error_obj) le 0) then goto,plot_it

 widget_control, t_error_obj.id, $
                  get_uvalue=uvalue

  wset, t_error_wid

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: goto, plot_it

     2: begin
        t_error_min_fac = t_error_min_fac*10.0
        goto, plot_it
     end

     3: begin
        t_error_min_fac = t_error_min_fac/10.0
        goto, plot_it
     end

     4: begin
        plot_mode = 2
        goto, plot_it
     end

     5: widget_control, t_error_obj.top, /destroy
     
  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  title = '!3Time-integration Error'
  ytitle = '!4d!3h/h!3'

  plot_def_new,'t_error'

temp=min(t_error[0:n_kinetic-1,1:n_time1])
err_min_arr=[1.e-9,1.e-8,1.e-7,1.e-6,1.e-5,1.e-4,1.e-3]
ierr=where(err_min_arr le temp, icount)
if(icount le 0) then imin=0 else imin=ierr[icount-1]
t_error_min=t_error_min_fac*err_min_arr[ierr[imin]]

  plot_io,[0],[0],$
          /nodata,$
          title=title,$
          ytitle=ytitle,$
          xstyle=1,$
          xminor=0,$
          xrange=[min(t),max(t)],$
          xtitle=csa_string,$
          ystyle=1,$
          yminor=0,$
          yrange=[t_error_min,2*max(t_error)],$
          color=line

  for i=0,n_kinetic-1 do begin
     oplot,t[1:n_time1],t_error[i,1:n_time1],$
           color=color_vec[i],$
           linestyle=line_vec[i]
  endfor

  plot_finish

  return

end
