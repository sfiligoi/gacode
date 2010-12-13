;
; C. Holland, UCSD
; v1.0: July 19, 2007: plots flux-surface avg chi/D as function of
; (r,t)
;
PRO DIFFUSION_I_RT_PLOT, SAVEFILE = savefile, ASCIIFILE = asciifile

  common GLOBAL
  common PLOT_VARIABLES
  common PROFILE_SIM_DATA
  

  ;;-------------------------------------------------------
  ;; Text labels
  ;;
  smf_tag,title,pname,ytitle

  title = title+' Diffusion'
  pname = 'diff_i-'+pname+'-rt'

  if (i_units eq 0) then begin
     plot_units = 1.0
     units = ' [units of (c!ds!n/a)!4q!3!ds!n!u2!n]'
  endif else begin
     plot_units = xunits[8]
     units = ' [units of (m!u2!n/sec)]'
  endelse

  title = title+' ['+strtrim(string(n_ss),2)+']'
  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; Pick out correct diffusivity
  ;;

  y_i = fltarr(n_r,n_time)

  if i_f ge 0 then begin
     y_i[*,*] = diff_i[i_spec,i_f,i_moment,*,*]
  endif else begin
     y_i[*,*] = diff_i[i_spec,0,i_moment,*,*]+$
       diff_i[i_spec,1,i_moment,*,*]
  endelse
  y_i *= plot_units
  ymin = MIN(y_i, MAX=ymax)

  plot_def_new,pname

  if (c_table_max > 0.0) then begin
     clevels = c_table_min+$
       findgen(nlevels)*(c_table_max-c_table_min)/(nlevels-1.0)
  endif
  loadct, CTab, /SILENT

  plot_position = [0.15,0.15,0.85,0.85]
  colorbar_position = FLTARR(4)
  colorbar_position[0] = plot_position[2] + $
    0.02 * (plot_position[2] - plot_position[0])
  colorbar_position[1] = plot_position[1]
  colorbar_position[2] = plot_position[2] + $
    0.05 * (plot_position[2] - plot_position[0])
  colorbar_position[3] = plot_position[3]


  contour,y_i,r,t,$
    nlevels=nlevels,$
    levels=clevels,$
    xtitle='!3r/a', $
    xstyle=1,$
    xrange=[min(r),max(r)],$
    ytitle=ytitle+units, $
    ystyle=1,$
    title=title,$
    /fill, $
    position=plot_position

  ;plot ss locations
  set_line_colors
  get_singsurf_vec,n_ss,r_surf  
  n_surf = n_elements(r_surf)
  for i=0,n_surf-1 do begin
     oplot,r_surf[i]*[1,1],[min(t),max(t)], color=color_vec[3]
  endfor

  ;plot lines denoting damping boundaries
  plot_bnd

  ;plot colorbar
  loadct, CTab, /SILENT
; 11/20/97- out until IDL bug fixed
;  COLORBAR, NCOLORS = !D.TABLE_SIZE, MINRANGE = ymin, MAXRANGE = ymax, $
;	/VERTICAL, /RIGHT, POSITION = colorbar_position, $
;	FORMAT = '(G10.2)'

  plot_finish
  set_line_colors

  ;write .sav file
  IF N_ELEMENTS(savefile) NE 0 THEN BEGIN
      SAVE, y_i, FILENAME = savefile
      MESSAGE, 'Wrote to ' + savefile, /INFO
  ENDIF

  ;write text file
  IF N_ELEMENTS(asciifile) NE 0 THEN BEGIN
      OPENW, lun, asciifile, /GET_LUN
      PRINTF, lun, y_i
      FREE_LUN, lun
      MESSAGE, 'Wrote to ' + asciifile, /INFO
  ENDIF

END ;DIFFSUION_I_RT_PLOT
