PRO MIDPLANE_RMS_FLUCAMP_PLOT, LOCAL_NORM = local_norm, FINITE_N=finite_n, $
  SAVEFILE = savefile, ASCIIFILE = asciifile
;
; 10/2/2007: updated calc routine
; 11/9/2007: added finite-n option
; 11/20/2007: added toroidal averaging
; 9/26/2008: updated to display magnetic fluctuations
; 5/11/2012: updated to plot in (%) rathern than absolute levels

  common GLOBAL
  common PROFILE_SIM_DATA
  common MIDPLANE_DATA
  common PLOT_VARIABLES

  pname = 'rms_flucamp_'+tag[i_pwr]
  IF KEYWORD_SET(finite_n) THEN  title = 'RMS finite-n ' $
  ELSE title = 'RMS ' 
  title += tag[i_pwr] + ' fluctuations '+$
          ' ['+strtrim(string(n_ss),2)+']'
  plot_def_new, pname

  ampfield = REFORM(pwr[i_pwr,*,*,it1:it2])

  flucpwr = 2*TOTAL(ABS(ampfield[*,1:*,*])^2,2)
  IF NOT(KEYWORD_SET(finite_n)) THEN flucpwr += $
	REFORM(ABS(ampfield[*,0,*])^2)
  rmsamp = SQRT(TOTAL(flucpwr,2)/(it2-it1+1))

  ;figure out ytitle and appropriate local norm
  IF (STRPOS(tag[i_pwr],'tem') NE -1) THEN BEGIN
      IF (STRPOS(tag[i_pwr],'elec') NE -1) THEN BEGIN
          ytitle = '!4d!XT!De!N(r)/T!De0!N'
          norm = REFORM(tem_s[n_spec-1,*])
      ENDIF ELSE BEGIN
          ytitle = '!4d!XT!Di!N(r)/T!Di0!N'
          norm = REFORM(tem_s[0,*])
      ENDELSE
  ENDIF ELSE IF (STRPOS(tag[i_pwr],'energy') NE -1) THEN BEGIN
      IF (STRPOS(tag[i_pwr],'elec') NE -1) THEN BEGIN
          ytitle = '!4d!Xp!De!N(r)/p!De0!N'
          norm = REFORM(tem_s[n_spec-1,*]*den_s[n_spec-1,*])
      ENDIF ELSE BEGIN
          ytitle = '!4d!Xp!Di!N(r)/p!Di0!N'
          norm = REFORM(tem_s[0,*]*den_s[0,*])
      ENDELSE
  ENDIF ELSE IF (STRPOS(tag[i_pwr],'density') NE -1) THEN BEGIN
      IF (STRPOS(tag[i_pwr],'elec') NE -1) THEN BEGIN
          ytitle = '!4d!Xn1De!N(r)/n!De0!N'
          norm = REFORM(den_s[n_spec-1,*])
      ENDIF ELSE BEGIN
          ytitle = '!4d!Xn!Di!N(r)/n!Di0!N'
          norm = REFORM(den_s[0,*])
      ENDELSE
  ENDIF ELSE IF (STRPOS(tag[i_pwr],'A_par') NE -1) THEN BEGIN
      ytitle = '(!Sc!R!A-!Ds!N/c)e!4d!XA(r)!D!9||!X!N/T!De0!N'
      norm = REFORM(tem_s[n_spec-1,*]) ;Te
  ENDIF ELSE IF (STRPOS(tag[i_pwr],'potential') NE -1) THEN BEGIN
      ytitle = 'e!4du!X(r)/T!De0!N'
      norm = REFORM(tem_s[n_spec-1,*]) ;Te
  ENDIF ELSE ytitle = "**Didn't find mean profile!**"
  IF KEYWORD_SET(local_norm) THEN BEGIN
      rmsamp /= norm
      ytitle = ytitle + '(r)'
  ENDIF ELSE ytitle = ytitle + '(r!D0!N)'

  ;convert to % rather than abs. value
  rmsamp *= 100.
  ytitle += ' (%)'

  mean_amp = MEAN(rmsamp[n_bnd:n_r-1-n_bnd])
  mean_sdom = SDOM(rmsamp[n_bnd:n_r-1-n_bnd])
  title += ' !3['+NUMTOSTRING(mean_amp)+'!9+!3'+NUMTOSTRING(mean_sdom)+']'

  PLOT, r, rmsamp, /XS, XTITLE = '!3r/a with' + t_string, TITLE=title, $
        YTITLE = ytitle

  OPLOT, [r[n_bnd],r[n_r-1-n_bnd]], [1,1]*mean_amp, color=color_vec[1]
  OPLOT, [r[n_bnd],r[n_r-1-n_bnd]], [1,1]*(mean_amp-mean_sdom), $
         color=color_vec[2]
  OPLOT, [r[n_bnd],r[n_r-1-n_bnd]], [1,1]*(mean_amp+mean_sdom), $
         color=color_vec[2]

  ;plot ss locations
  set_line_colors
  get_singsurf_vec,n_ss,r_surf  
  n_surf = n_elements(r_surf)
  for i=0,n_surf-1 do begin
     oplot,r_surf[i]*[1,1],[min(t),max(t)], color=color_vec[3]
  endfor

  ;plot lines denoting damping boundaries
  plot_bnd

  plot_finish

  ;write .sav file
  IF N_ELEMENTS(savefile) NE 0 THEN BEGIN
      SAVE, r, rmsamp, FILENAME = savefile
      MESSAGE, 'Wrote to ' + savefile, /INFO
  ENDIF

  ;write text file
  IF N_ELEMENTS(asciifile) NE 0 THEN BEGIN
      OPENW, lun, asciifile, /GET_LUN
      PRINTF, lun, r, rmsamp
      FREE_LUN, lun
      MESSAGE, 'Wrote to ' + asciifile, /INFO
  ENDIF

END ;RMS_FLUCAMP_PLOT
