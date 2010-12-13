pro zmoment_event, zmoment_obj

  common GLOBAL
  common PLOT_VARIABLES
  common ZMOMENT_DATA
  
  widget_control, zmoment_obj.id, $
      get_uvalue=uvalue

  wset, zmoment_wid

  if (i_f eq n_field) then i_f = n_field-1

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 

    1: begin
      i_spec = i_spec+1
      if (i_spec eq n_kinetic) then i_spec = 0
      goto, plot_it
    end   

    2: begin
      i_moment = 1-i_moment
      goto, plot_it
    end   

    3: begin
      i_f = i_f+1
      if (i_f eq n_field) then i_f = 0
      goto, plot_it
      goto, plot_it
    end

    5: begin
      z_max = z_max*0.5
      goto, plot_it
    end
    6: begin
      z_max = z_max*2.0
      goto, plot_it
    end

    7: begin
      i_ptype = 1-i_ptype
      goto, plot_it
    end

    8: begin
      plot_mode = 2
      goto, plot_it
    end

    9: widget_control, zmoment_obj.top, /destroy

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  y = fltarr(n_r,n_time)

  smf_tag,title,pname,ytitle
  title = title+' Moment'
  pname = 'moment-'+pname

  y[*,*] = zmoment[*,i_spec,i_moment,*]

  if (i_ptype eq 0) then begin

    plot_def_new,pname+'_CMAP'

    loadct,CTab

    !p.region = [0,0.15,1.0,1.0]

    contour,transpose(y),t,r, $
        nlevels=nLevels,$
        /fill, $
        ytitle='!3r/a', $
        ystyle=1,yrange=[max(r),min(r)],$
        xtitle=csa_string, $
        xstyle=1,xrange=[t_min,t_max],$
        zstyle=1,zrange=[z_min*min(y),z_max*max(y)],$
        title=title, $
        charsize=c_size*c_scale

    !p.region = [0,0,1.0,0.15]

    color_map,z_min*min(y),z_max*max(y),nLevels

    !p.region = 0

    set_line_colors

  endif else begin

    plot_def_new,pname+'_LINE'

    plot,[0],[0],$
        /nodata, $
        ytitle='!3r/a', $
        ystyle=1,yrange=[max(r),min(r)],$
        xtitle=csa_string, $
        xstyle=1,xrange=[t_min,t_max],$
        title=title

    dr = max(r)-min(r) 

    y_mult = dr/max(abs(y))/z_max

    for i=0,n_r-1 do begin
      oplot,t,y[i,*]*y_mult+r[i]
    endfor

  endelse

  plot_finish

  return

end
