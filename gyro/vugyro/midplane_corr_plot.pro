;
; C.Holland, UCSD
; v1.0: July 19, 2007
;
; Plots midplane (r,alpha) correlation functions
;
; 11/9/2007: added finite-n option
; 12/12/2007: CH- change radial correlation scale length to use rho_s rather
;                 than minor radius a as norm.  Corrected box-averaged r_corr
;                 to use only points inside n_bnd.  Corrected dalpha axis to 
;                 go only from -pi/n_dn to pi/n_dn (did -pi to pi before).
;                 Normalized l_alpha corr scale to pi.

PRO midplane_corr_plot, i_r0, which_plot, FINITE_N=finite_n, $
  LOCAL_CORR = local_corr, SAVEFILE = savefile, ASCIIFILE = asciifile

  common GLOBAL
  common PROFILE_SIM_DATA
  common PLOT_VARIABLES
  common MIDPLANE_DATA

  ;create plot data
  field = REFORM(pwr[i_pwr,*,*,it1:it2])
  IF KEYWORD_SET(finite_n) THEN BEGIN
      field[*,0,*] = 0.         ;remove ZF component
      title = 'finite-n '
  ENDIF ELSE BEGIN
      FOR i_r = 0, n_r-1 DO field[i_r,0,*] -= MEAN(field[i_r,0,*])
      title = ''
  ENDELSE


  IF (which_plot EQ 0) THEN BEGIN ;C(dr,alpha=0)
      xtitle = '!4D!X!3r/a'
      IF KEYWORD_SET(local_corr) THEN BEGIN ; C(r0,dr)
          ;least error-prone way
          n_r2 = n_r - 2*n_bnd
          r0 = r[i_r0]
          xaxis = FLTARR(n_r2)
          corr = FLTARR(n_r2)
          FOR i_r = n_bnd, n_r-1-n_bnd DO BEGIN
              xaxis[i_r-n_bnd] = r[i_r] - r0

              corr[i_r-n_bnd] += FLOAT(TOTAL(CONJ(field[i_r0,0,*])*$
                                       field[i_r,0,*]))
              corr[i_r-n_bnd] += 2*FLOAT(TOTAL(TOTAL(CONJ(field[i_r0,1:*,*])*$
                                                     field[i_r,1:*,*],3),2))
          ENDFOR
          corr /= MAX(corr)

          title += tag[i_pwr] + ' !12<!XC(!4D!Xr, r0 = ' + $
                  NUMTOSTRING(r[i_r0],LEN=5) + ')!12>!X!D!4a!X,t!N with' + $
                  t_string
          lcorr = EXP_FIT(xaxis, corr,0.05,/ABSX)
          corr_fit = EXP(-ABS(xaxis)/lcorr)
          corr_label = '!8l!X!Dr!N/!4q!X!Ds!N = ' + NUMTOSTRING(lcorr/rho_s,LEN=6)
      ENDIF ELSE BEGIN ; <C(dr)>_r
          field_kn = FFT(field[n_bnd:n_r-1-n_bnd,*,*],DIM=1)
          field_kn[*,1:*,*] *= SQRT(2) ;account for power diff in finite-n
          avg_kpow = TOTAL(TOTAL(ABS(field_kn[*,*,*])^2,3),2)
          corr = FLOAT(FFT(avg_kpow,/INVERSE))
          corr = SHIFT(corr, N_R/2-n_bnd)/MAX(corr)
          xaxis = r[n_bnd:n_r-1-n_bnd] - r[N_R/2]
          title += tag[i_pwr] + $
                  ' !12<!XC(!4D!Xr)!12>!X!Dr,!4a!X,t!N'+$
                  ' with' + t_string
          pname = 'midplane_avg_rcorr'

          lcorr = EXP_FIT(xaxis, corr,0.05,/ABSX)
          corr_fit = EXP(-ABS(xaxis)/lcorr)
          corr_label = '!8l!X!Dr!N/!4q!X!Ds!N = ' + NUMTOSTRING(lcorr/rho_s,LEN=6)
      ENDELSE
  ENDIF ELSE IF (which_plot EQ 1) THEN BEGIN ;C(dalpha)
      IF KEYWORD_SET(local_corr) THEN BEGIN  ;C(r0, dalpha)
          powspect = FLTARR(3*n_n)
          powspect[0:n_n-1] = TOTAL(ABS(field[i_r0,*,*])^2,3)/(it2-it1+1)
          powspect[2*n_n:3*n_n-2] = REVERSE(powspect[1:n_n-1])

          corr = SHIFT(FLOAT(FFT(powspect)),3*n_n/2)
          corr /= MAX(corr)     ;normalize autocorrelation to 1
          xaxis = (2*!PI*FINDGEN(3*n_n)/(3*n_n) - !PI)/n_dn

          xtitle = '!4Da!X'
          title += tag[i_pwr] + ' !12<!XC(!4Da!X, r0 = ' + $
                  NUMTOSTRING(r[i_r0],LEN=5) + ')!12>!X!D!4a!X,t!N with' + $
                  t_string
          pname = 'midplane_local_alphacorr'

          corr_env = ENVELOPE(corr)
          lcorr = EXP_FIT(xaxis, corr_env,0.05,/ABSX)
          corr_fit = EXP(-ABS(xaxis)/lcorr)
          corr_label = '!8l!X!D!4a!X!N/!4p!X = ' + NUMTOSTRING(lcorr/!PI,LEN=6)
      ENDIF ELSE BEGIN ;<C(dalpha)>_r
          corr = FLTARR(3*n_n)
          FOR i_r = n_bnd, N_R-1-n_bnd DO BEGIN
              powspect = FLTARR(3*n_n)
              powspect[0:n_n-1] += TOTAL(ABS(field[i_r,*,*])^2,3)/(it2-it1+1)
              powspect[2*n_n:3*n_n-2] = REVERSE(powspect[1:n_n-1])
              i_rcorr = SHIFT(FLOAT(FFT(powspect)),3*n_n/2)
              i_rcorr /= MAX(i_rcorr) ;normalize autocorrelation to 1
              corr += i_rcorr
          ENDFOR
          corr /= N_R-2*n_bnd
          xaxis = (2*!PI*FINDGEN(3*n_n)/(3*n_n) - !PI)/n_dn


          xtitle = '!4Da!X'
          title += tag[i_pwr] + $
                  ' !12<!XC(!4Da!X)!12>!X!Dr,!4a!X,t!N with' + t_string
          pname = 'midplane_avg_alphacorr'
      
          corr_env = ENVELOPE(corr)
          lcorr = EXP_FIT(xaxis, corr_env,0.05,/ABSX)
          corr_fit = EXP(-ABS(xaxis)/lcorr)
          corr_label = '!8l!X!D!4a!X!N/!4p!X = ' + NUMTOSTRING(lcorr/!PI,LEN=6)
      ENDELSE
  ENDIF ELSE IF (which_plot EQ 2) THEN BEGIN ;C(dr,dalpha)
      IF KEYWORD_SET(local_corr) THEN BEGIN
          xspect = COMPLEXARR(n_r,3*n_n)
          FOR it=0,it2-it1 DO FOR i_r = n_bnd, n_r-1-n_bnd DO $
            xspect[i_r,0:n_n-1] += $
              CONJ(field[i_r0,*,it])*field[i_r,*,it]
          xspect[*,2*n_n:3*n_n-2] = REVERSE(CONJ(xspect[*,1:n_n-1]),2)
          corr = SHIFT(FLOAT(FFT(xspect,DIM=2)),0,3*n_n/2)
          corr /= MAX(corr)

          xaxis = r - r[i_r0]
          xtitle = '!4D!Xr/a'
          yaxis = (2*!PI*FINDGEN(3*n_n)/(3*n_n) - !PI)/n_dn
          ytitle = '!4Da!X'
          title += tag[i_pwr] + ' !12<!XC(!4D!Xr,!4Da!X, r0 = ' + $
                  NUMTOSTRING(r[i_r0],LEN=5) + ')!12>!X!D!4a!X,t!N with' + $
                  t_string

      ENDIF ELSE BEGIN ;r-avg C(dr,dalpha)
          n_r2 = n_r - 2*n_bnd
          xspect = FLTARR(n_r2,3*n_n)
          FOR it=0,it2-it1 DO BEGIN
              field_kn = FFT(field[n_bnd:n_r-1-n_bnd,*,it], DIM=1)
              xspect[*,0:n_n-1] += ABS(field_kn)^2
          ENDFOR
          xspect[*,2*n_n:3*n_n-2] = REVERSE(CONJ(xspect[*,1:n_n-1]),2)
          corr = SHIFT(FLOAT(FFT(xspect)),n_r2/2,3*n_n/2)
          corr /= MAX(corr)

          xaxis = r[n_bnd:n_r-1-n_bnd] - r[i_r0]
          xtitle = '!4D!Xr/a'
          yaxis = (2*!PI*FINDGEN(3*n_n)/(3*n_n) - !PI)/n_dn
          ytitle = '!4Da!X'
          title += tag[i_pwr] + ' !12<!XC(!4D!Xr,!4Da!X)!12>!X!Dr,!4a!X,t!N'+$
                   ' with ' + t_string
      ENDELSE
  ENDIF

  ;do plot
  plot_def_new, pname+'_'+tag[i_pwr]

  IF (which_plot EQ 2) THEN BEGIN ;2d plot
      loadct, CTab, /SILENT
      if (c_table_max > 0.0) then begin
          clevels = c_table_min+$
                    findgen(nlevels)*(c_table_max-c_table_min)/(nlevels-1.0)
      endif

      plot_position = [0.15,0.15,0.85,0.85]
      colorbar_position = FLTARR(4)
      colorbar_position[0] = plot_position[2] + $
        0.02 * (plot_position[2] - plot_position[0])
      colorbar_position[1] = plot_position[1]
      colorbar_position[2] = plot_position[2] + $
        0.05 * (plot_position[2] - plot_position[0])
      colorbar_position[3] = plot_position[3]

      contour,corr,xaxis,yaxis,$
              nlevels=nlevels,$
              levels = -1. + 2*FINDGEN(nlevels)/(nlevels-1), $
              c_colors=clevels,$
              xtitle=xtitle, $
              xstyle=1,$
              ytitle=ytitle, $
              ystyle=1,$
              title=title,$
              /fill, $
              position=plot_position

      contour, corr, xaxis, yaxis, nlevels=1, levels=0., /overplot, $
               position=plot_position, c_linestyle=2

     ;plot colorbar
;11/20/07- out until IDL bug fixed
      COLORBAR, NCOLORS = !D.TABLE_SIZE, MINRANGE = -1., $ ;MIN(corr), $
	MAXRANGE = 1, /VERTICAL, /RIGHT, POSITION = colorbar_position, $
       FORMAT = '(G10.2)'

  ENDIF ELSE BEGIN ;1d plot
      PLOT, xaxis, corr, /XS, XTITLE=xtitle, TITLE = title
      
      ;correlation envelope and fit
      OPLOT, xaxis, corr_fit, LINESTYLE=2, COLOR=100
      XYOUTS, MIN(xaxis) + 0.1*(MAX(xaxis)-MIN(xaxis)), 0.9, corr_label

      ;line to denote C=0
      OPLOT, xaxis, FLTARR(N_ELEMENTS(xaxis)), LINESTYLE=1
  ENDELSE
  plot_finish
  set_line_colors

  ;write .sav file
  IF N_ELEMENTS(savefile) NE 0 THEN BEGIN
      IF N_ELEMENTS(yaxis) NE 0 THEN SAVE, xaxis, yaxis, corr, $
        FILENAME=savefile $
      ELSE SAVE, xaxis, corr, FILENAME = savefile
      MESSAGE, 'Wrote to ' + savefile, /INFO
  ENDIF

  ;write text file
  IF N_ELEMENTS(asciifile) NE 0 THEN BEGIN
      OPENW, lun, asciifile, /GET_LUN
      PRINTF, lun, xaxis
      IF N_ELEMENTS(yaxis) NE 0 THEN PRINTF, lun, yaxis
      PRINTF, lun, corr
      FREE_LUN, lun
      MESSAGE, 'Wrote to ' + asciifile, /INFO
  ENDIF

END ;midplane_corr_plot
