pro source_t_event, source_t_obj

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_SOURCE_T

  widget_control, source_t_obj.id, $
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
        i_moment = 1-i_moment
        goto, plot_it
     end   
     
     4: begin
        zoom = 2*zoom
        goto, plot_it
     end

     5: begin
        zoom = zoom/2
        goto, plot_it
     end

     6: begin
        plot_mode = 2
        goto, plot_it
     end

     7: begin
        plot_export = 1
        goto, plot_it
     end

     8: widget_control, source_t_obj.top, /destroy
     
  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  yrt = fltarr(n_r,n_time)

  if (i_spec lt n_ion) then begin
     tag_s = 'Ion '+strtrim(string(i_spec+1),2)
     ptag_s = 'ion'+strtrim(string(i_spec+1),2)
  endif else begin
     tag_s = 'Electron'
     ptag_s = 'elec'
  endelse
  if (i_moment eq 0) then begin
     tag_m = 'Density'
     ptag_m = 'density'
  endif else begin
     tag_m = 'Energy'
     ptag_m = 'energy'
  endelse
 
  title = tag_s+' '+tag_m+' Source'
  pname = 'source_t-'+ptag_s+ptag_m

  yrt[*,*] = source[i_spec,2+i_moment,*,*]

  plot_def_new,pname

  y = fltarr(n_time)
  y[*] = yrt[n_r/2,*]

  plot,[0],[0],$
       /nodata,$
       title=title,$
       xstyle=1,$
       xminor=0,$
       xrange=[min(t),max(t)],$
       xtitle=csa_string,$
       ystyle=1,$
       yminor=0,$
       yrange=[-1,1]*max(y[*])/zoom,$
       ytitle=ytitle,$
       color=line

  oplot,t,y,color=color_vec[0]   

  plot_finish

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for i=0,n_r-1 do begin
        printf,1,t[i],y[i]
     endfor
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  endif
  ;;---------------------------------------------------

  return

end
