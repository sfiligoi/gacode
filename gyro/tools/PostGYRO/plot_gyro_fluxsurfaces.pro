PRO plot_gyro_fluxsurfaces, data, SF = sf, N_LEVELS = n_levels, APHYS=aphys, $
  OPLOT=overplot,THICK=thick, _EXTRA=extra
;
; C. Holland, UCSD
; v1.0 2/21/2007
;
; Given a data set from GET_GYRO_DATA(), draws contours of flux
; surfaces
;
; SF determines interplotation used to generate theta gridpoints
; N_LEVELS is the number of contours to plot; default is inner and
; outer surfaces
;
; v2.0: 5/23/2007
; updated to use Aphys, OPLOT, _EXTRA keywords
;

  DEFAULT, sf, 1
  DEFAULT, n_levels, 2
  DEFAULT, thick, 1
  n_y = sf*data.n_theta_plot
  n_r = data.profile_data.n_r
  GYRO_GENERATE_RZCOORDS, data.profile_data, n_y, R, Z

  idx = ((data.N_R-1)/(n_levels-1))*INDGEN(n_levels)
  i_r = n_r - 1
  i_r = idx[n_levels-1]

  R0 = data.profile_data.R0
  rr = MAX(ABS(R-R0[i_r]))
  zr = MAX(ABS(Z))
  pr = zr > rr

  IF (N_ELEMENTS(Aphys) EQ 0) THEN BEGIN
      xtitle = 'R/a'
      ytitle = 'Z/a'
  ENDIF ELSE BEGIN
      R *= Aphys
      R0 *= Aphys
      Z *= Aphys
      pr *= Aphys
      xtitle = 'R (cm)'
      ytitle = 'Z (cm)'
  ENDELSE

  IF KEYWORD_SET(overplot) THEN OPLOT, R[i_r,*], Z[i_r,*], $
    THICK=thick, _EXTRA = extra $
  ELSE PLOT, R[i_r,*], Z[i_r,*], /XS, /YS, XTITLE = xtitle, YTITLE = ytitle, $
             XRANGE=[R0[i_r]-pr,R0[i_r]+pr], YRANGE=[-pr,pr], $
             THICK=thick, XTHICK=thick,YTHICK=thick, CHARTHICK=thick,$
             _EXTRA = extra

  FOR ii = 0, n_levels-2 DO BEGIN
      OPLOT, R[idx[ii],*], Z[idx[ii],*], THICk=thick, _EXTRA = extra
  ENDFOR

END ;plot_gyro_fluxsurfaces

