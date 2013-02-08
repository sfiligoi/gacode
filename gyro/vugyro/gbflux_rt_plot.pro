;
; C. Holland, UCSD
; v1.0: May 10, 2012: plots flux-surface avg fluxes as function of
; (r,t)
; v2.0: Dec 20, 2012: updated for new input.profiles.extra structure
;
PRO GBFLUX_RT_PLOT, SAVEFILE = savefile

  common GLOBAL
  common PLOT_VARIABLES
  common PROFILE_SIM_DATA

  ;;-------------------------------------------------------
  ;; Text labels
  ;;
  gbflux_tag,title,pname,fluxtitle,units1,units2,units3
  pname = 'gbflux_rt-'+pname

   case (i_units) of

     0: begin
        ;; flux in gyrobohm units
        unit_norm = FLTARR(n_r) + 1.0
        units = units1
     end

     1: begin
        ;; flux in W/m^2, etc.
        unit_norm = FLTARR(n_r) + xunits[9+i_moment]
        units = units2
     end 

     2: begin 
        ;; flux in W, etc.
        unit_norm = xunits[9+i_moment]
        if (exists_exp_derived) then begin
            ;; Multiply by V'
;            unit_norm *= INTERPOL(exp_derived[22,*],r_from_rho,r)
; CH fix 12.20.2012: update for new input.profiles.extra structure
            unit_norm *= INTERPOL(exp_derived[23,*],r_from_rho,r)
        endif
        units = units3
     end

  endcase

  xtitle = '!3r/a'

  plot_def_new,pname
 ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; Pick out correct flux
  y_rt = fltarr(n_r,n_time)
  if i_f ge 0 then begin
     y_rt[*,*] = gbflux_i[i_spec,i_f,i_moment,*,*]
  endif else begin
     y_rt[*,*] = 0.0
     for ix=0,n_field-1 do begin
        y_rt[*,*] += gbflux_i[i_spec,ix,i_moment,*,*]
     endfor
  endelse

  ;;apply units correction
  y_rt *= (unit_norm # (FLTARR(n_time) + 1))

  ;;-------------------------------------------------------


  ;;-------------------------------------------------------
  ;; set up plot
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

  xtitle = '!3r/a'
  ytitle = csa_string
  title = fluxtitle + units + ' ' + title
  contour,y_rt,r,t,$
          nlevels=nlevels,$
          levels=clevels,$
          xtitle=xtitle, $
          xstyle=1,$
          xrange=[min(r),max(r)],$
          ytitle=ytitle, $
          ystyle=1,$
          title = title+' ['+strtrim(string(n_ss),2)+']',$
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
  ymin = MIN(y_rt, MAX=ymax)
; 11/20/97- out until IDL bug fixed
  COLORBAR, NCOLORS = !D.TABLE_SIZE, MINRANGE = ymin, MAXRANGE = ymax, $
	/VERTICAL, /RIGHT, POSITION = colorbar_position, $
	FORMAT = '(G10.2)'

  plot_finish
  set_line_colors

  ;write .sav file
  IF N_ELEMENTS(savefile) NE 0 THEN BEGIN
      SAVE, y_rt, r, t, xtitle, ytitle, title, FILENAME = savefile
      MESSAGE, 'Wrote to ' + savefile, /INFO
  ENDIF

END ;GBFLUX_RT_PLOT
