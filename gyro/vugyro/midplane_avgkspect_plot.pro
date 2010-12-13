PRO midplane_avgkspect_plot, pwr_idx, krmax, USE_N = use_n, $
  NO_KR = no_kr, PLOTMIN = plotmin, PLOTMAX = plotmax,$
  SAVEFILE = savefile, ASCIIFILE = asciifile
;
; C. Holland, UCSD
; Plotting routine for improved midplane <S(kx,ky)>
;
; v1.0: Jan 26, 2007- debugged (I think), ready for release into the
;                     wild
;
; v1.1: Jul 19, 2007- updated to use vuGyro postscript plotting
;                     options, no_kr option
; v1.2: Aug 16, 2007- added use_n feature to allow n or r-dep k_theta plots

  common GLOBAL
  common PROFILE_SIM_DATA
  common PLOT_VARIABLES
  common MIDPLANE_DATA


  
  ;create spectrum and axes
  field = REFORM(pwr[pwr_idx,*,*,*])
  title = '!3log!D10!N!12<!X|' + tag[pwr_idx]
  IF KEYWORD_SET(no_kr) THEN BEGIN
      spect = ALOG10(TOTAL(ABS(field[*,*,it1:it2])^2,3)/(it2-it1+1))

      xaxis = r # (FLTARR(n_n) + 1)
      xtitle = '!3r/a'
      xrange = [MIN(r),MAX(r)]

      IF KEYWORD_SET(use_n) THEN BEGIN
          yaxis = (FLTARR(n_r) + 1) # n_tor
          ytitle = '!3n'
          title = title + '(r,n)|!U2!N!12>!X'
      ENDIF ELSE BEGIN
          yaxis = ((q_m/r_m)*rho_s*(sqrt(tem_s(n_spec-1,*))/b_unit_s))#n_tor
          ytitle = kt_rho_string
          title = title + '(r,k!D!4h!X!N)|!U2!N!12>!X'
      ENDELSE
  ENDIF ELSE BEGIN
      spect = FLTARR(n_r, n_n)
      FOR i_time = it1, it2 DO FOR i_n = 0, n_n-1 DO $
        spect[*,i_n] += ABS(FFT(field[*,i_n,i_time]))^2
      spect = ALOG10(SHIFT(spect, n_r/2-1,0)/(it2 - it1 + 1))

      xaxis = kr_rho
      xtitle = '!3' + kr_rho_string
      xrange = [-krmax, krmax]
      IF KEYWORD_SET(use_n) THEN BEGIN
          yaxis = n_tor
          ytitle = '!3n'
          title = title + '(k!Dr!N,n)|!U2!N!12>!X'
      ENDIF ELSE BEGIN
          yaxis = kt_rho
          ytitle = '!3' + kt_rho_string
          title = title + '(k!Dr!N,k!D!4h!X!N)|!U2!N!12>!X'
      ENDELSE
  ENDELSE
  title = title + ' with' + t_string

  pname = 'midplane_' + tag[pwr_idx] + '_avgkspect'
  plot_position = [0.15,0.15,0.85,0.85]

  ;color table and bar info
  LOADCT, CTab, /SILENT
  DEFAULT, plotmin, MIN(spect)
  DEFAULT, plotmax, MAX(spect)
  spect = spect > plotmin
  spect = spect < plotmax
  cmin = !D.TABLE_SIZE*(MIN(spect) - plotmin)/(plotmax - plotmin)
  cmax = !D.TABLE_SIZE*(MAX(spect) - plotmin)/(plotmax - plotmin) - 1
  color_array = ROUND((cmax - cmin)*(FINDGEN(nlevels)/(nlevels-1)) + cmin)
  colorbar_position = FLTARR(4)
  colorbar_position[0] = plot_position[2] + $
    0.02 * (plot_position[2] - plot_position[0])
  colorbar_position[1] = plot_position[1]
  colorbar_position[2] = plot_position[2] + $
    0.05 * (plot_position[2] - plot_position[0])
  colorbar_position[3] = plot_position[3]
 
  plot_def_new, pname
 
  ;do plot
  CONTOUR, spect, xaxis, yaxis, /XS, /YS, NLEVELS = nlevels, /FILL, $
           XTITLE = xtitle, YTITLE = ytitle, TITLE = title, $
           C_COLORS = color_array, XRANGE = xrange, POSITION = plot_position
  IF KEYWORD_SET(no_kr) THEN BEGIN
      set_line_colors
      plot_bnd
  ENDIF

  ;colorbar
  LOADCT, CTab, /SILENT
; 11/20/07- out until IDL bug fixed
;  COLORBAR, NCOLORS = !D.TABLE_SIZE, MINRANGE = plotmin, $
;	MAXRANGE = plotmax, /VERTICAL, /RIGHT, POSITION = colorbar_position, $
;       FORMAT = '(G10.2)'

  plot_finish
  set_line_colors

  ;write .sav file
  IF N_ELEMENTS(savefile) NE 0 THEN BEGIN
      SAVE, kr_rho, kt_rho, spect, kr_rho_string, kt_rho_String, $
            FILENAME = savefile
      MESSAGE, 'Wrote to ' + savefile, /INFO
  ENDIF

  ;write text file
  IF N_ELEMENTS(asciifile) NE 0 THEN BEGIN
      OPENW, lun, asciifile, /GET_LUN
      PRINTF, lun, kr_rho
      PRINTF, lun, kt_rho
      PRINTF, lun, spect
      PRINTF, lun, kr_rho_string
      PRINTF, lun, kt_rho_string
      FREE_LUN, lun
      MESSAGE, 'Wrote to ' + asciifile, /INFO
  ENDIF

END ;midplane_avgkspect_plot
