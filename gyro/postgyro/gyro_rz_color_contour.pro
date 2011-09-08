PRO GYRO_RZ_COLOR_CONTOUR, field, R, Z, XRANGE=xrange, YRANGE=yrange, $
  PLOTRANGE = plotrange, Aphys = Aphys, THICK = thick, CHARSIZE = cs, $
  PLOT_POSITION = full_plot_position, _EXTRA = extra
;
; C. Holland, UCSD
;
; makes contour (R,Z) plot of field, with colors evenly distributed
; between (-plotrange,plotrange).  Default is (R-R0)/a, Z/a, but
; setting Aphys = XX cm changes axes and labels to R-R0 and Z in cm.
;
; v1.0 11/22/06
; v2.0 5/22/07: added support for plotrange = 2-element vector to
; explicitly specify plotmin/max, 1 element for equal&opposite values
; v3.0 9/7/11: updated default xrange, yrange calculations to give 
; accurate aspect ratio visulaization i.e. LX=LY by default now
;

  IF (N_ELEMENTS(Aphys) NE 0) THEN BEGIN
      Rplot = R*Aphys
      Zplot = Z*Aphys
      xtitle = 'R (cm)'
      ytitle = 'Z (cm)'
  ENDIF ELSE BEGIN
      Rplot = R
      Zplot = Z
      xtitle = 'R/a'
      ytitle = 'Z/a'
  ENDELSE
  DEFAULT, yrange, 1.1*MAX(ABS(Zplot))*[-1,1]
  DEFAULT, xrange, yrange + (MIN(Rplot)+MAX(Rplot))/2

  IF (N_ELEMENTS(plotrange) EQ 2) THEN BEGIN
      plotmin = plotrange[0]
      plotmax = plotrange[1]
  ENDIF ELSE BEGIN
      DEFAULT, plotrange, MAX(ABS(field))
      plotmax = plotrange
      plotmin = -plotmax
  ENDELSE
  plot_Data = field
  plot_Data = (field > plotmin) < plotmax

  nlevels = 29
  cmin = !D.TABLE_SIZE*(MIN(plot_data) - plotmin)/(plotmax - plotmin)
  cmax = !D.TABLE_SIZE*(MAX(plot_data) - plotmin)/(plotmax - plotmin) - 1
  color_array = ROUND((cmax - cmin)*(FINDGEN(nlevels)/(nlevels-1)) + cmin)

   IF KEYWORD_SET(full_plot_position) THEN BEGIN
	plot_position = FLTARR(4);
	plot_position[0] = full_plot_position[0] + $
		0.10 * (full_plot_position[2] - full_plot_position[0])
	plot_position[1] = full_plot_position[1] + $
		0.15 * (full_plot_position[3] - full_plot_position[1])
	plot_position[2] = full_plot_position[0] + $
		0.80 * (full_plot_position[2] - full_plot_position[0])
	plot_position[3] = full_plot_position[1] + $
		0.85 * (full_plot_position[3] - full_plot_position[1])
  ENDIF ELSE plot_position = [0.15,0.2,0.8,0.85]
  colorbar_position = FLTARR(4)
  colorbar_position[0] = plot_position[2] + $
	0.02 * (plot_position[2] - plot_position[0])
  colorbar_position[1] = plot_position[1]
  colorbar_position[2] = plot_position[2] + $
	0.05 * (plot_position[2] - plot_position[0])
  colorbar_position[3] = plot_position[3]

  DEFAULT, thick, 1
  DEFAULT, cs, 1.5
  CONTOUR, TRANSPOSE(plot_data), Rplot, Zplot, /FILL, $
           NLEVELS = nlevels, C_COLORS=color_array, $
           POSITION = plot_position, $
           XRANGE = xrange, /XS, $
           YRANGE = yrange, /YS, $
           XTITLE = xtitle, YTITLE = ytitle, $
           THICK=thick, XTHICK=thick, YTHICK=thick, CHARTHICK=thick,$
           CHARSIZE=cs, _EXTRA = extra

  COLORBAR, NCOLORS = !D.TABLE_SIZE, MINRANGE = plotmin, MAXRANGE = plotmax, $
	/VERTICAL, /RIGHT, CHARSIZE = cs, POSITION = colorbar_position, $
	FORMAT = '(G10.2)', XTHICK = thick, YTHICK = thick, CHARTHICK = thick

END ;GYRO_RZ_COLOR_CONTOUR
