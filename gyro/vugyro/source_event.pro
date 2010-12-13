pro source_event, source_obj

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_SOURCE

  widget_control, source_obj.id, $
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

     8: widget_control, source_obj.top, /destroy
     
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
  pname = 'source-'+ptag_s+ptag_m

  yrt[*,*] = source[i_spec,2+i_moment,*,*]

  plot_def_new,pname

  yt = fltarr(n_time)
  y  = fltarr(n_r)

  for ii=0,n_r-1 do begin
     yt[*] = yrt[ii,*]
     diff_stat_fast,yt,it1,it2,ave_y
     y[ii] = ave_y
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
       yrange=[-1,1]*max(y[*])/zoom,$
       ytitle=ytitle,$
       color=line

  oplot,r,y,color=color_vec[0]   
  oplot,r,0*r,linestyle=1

  plot_finish

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for i=0,n_r-1 do begin
        printf,1,r[i],y[i]
     endfor
     printf,1,'Average:',t_min,t_max
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  endif
  ;;---------------------------------------------------

  return

end
