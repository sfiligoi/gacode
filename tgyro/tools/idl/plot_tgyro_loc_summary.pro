PRO plot_tgyro_loc_summary, simdir, N_it = N_it, RHO=rho, MKS=mks,$
  PS = ps, DATA2 = data2, Nit2=N_it2, SUMMARY=summary, $
  OLD=old, _EXTRA=extra
;
; C. Holland, UCSD
;
; v1.0: Oct. 29, 2008
;
; basic summary plotting routine.  Plots Ti, Te, Qi, and Qe of a
; given iteration vs. their iteration 0 values, and the corresponding
; target fluxes (as well as the neoclassical fluxes).  Recommend making
; IDL window at least 768x768 in size to plot clearly.
;
; Basic Usage
; --------------
; To plot the sample TGYRO simulation wtest_1 (from tgyro -g):
; IDL> plot_tgyro_loc_summary, 'wtest_1'
;
; This will plot temperature/density/rotation profiles on the left, and fluxes
; in gyroBohm units on the right, all against normalized midplane minor
; radius.  The results will correspond to the last iteration of the simulation.
;
; To plot the results of the 1st iteration, against normalized torodial flux,
; with fluxes in MKS units, enter
;
; IDL> plot_tgyro_loc_summary, 'wtest_1', N_it=1, /rho, /mks
;
; note that N_it=0 corresponds to the inital calculation before TGYRO has
; started iterating.
;
; To compare the results of two different times (or different simulations),
; use the following commands:
; IDL> data2 = get_tgyro_loc_data('wtest_1')
; IDL> plot_tgyro_loc_summary, 'wtest_1', data2=data2, n_it2 = 0
;
; The example above would plot the results of the last iteration of wtest_1 in
; red, and the inital pre-iteration results in yellow (assuming IDL color_table 5
; is loaded).
;
;
; KEYWORDS
; simdir: string containing name of valid directory in $TGYRO_DIR/sim
; N_it: iteration # to plot, defaults to last iteration
; RHO: set=1 to plot quantities vs rho, 0 (default) vs rmin/a
; MKS: set=1 to plot fluxes in MW/m**2, 0 (default) in gB units
; PS: set=1 to make characters slightly larger and lines thicker for
;     sending to Postscript files.
; DATA2: data structure from get_tgyro_loc_data to overplot (at Nit2)
; in yellow
;
; v1.0.1: Nov. 11, 2008
; Added SUMMARY flag.  Set equal to 1 to make a file called
; summary.eps in the simulation directory with the plotted figures. 
;
; v1.1: Apr. 15, 2009
; Added _EXTRA flag to pass to additional keywords to GET_TGYRO_LOC_DATA
;
; v1.2: May 11, 2009
; Added plot_n window to also plot n_e, gamma_e profiles
;
; v1.2.1: July 20, 2009
; Made gB unit plots on log-lin scale, rather than lin-lin
;
; v2.0: September 3, 2009
; Updated to plot rotation/momentum fluxes, compatibility with latest
; get_tgyro_loc_data
;

  IF KEYWORD_SET(ps) THEN BEGIN
     thick = 6
     cs = 2
  ENDIF ELSE BEGIN
     thick = 1
     cs = 2
  ENDELSE
  IF KEYWORD_SET(plot_n) THEN cs *= 2

  IF KEYWORD_SET(old) THEN data = GET_TGYRO_LOC_DATA_OLD(simdir) $
  ELSE data = GET_TGYRO_LOC_DATA(simdir)
  IF (N_ELEMENTS(data2) GT 0) THEN d2flag = 1 ELSE d2flag = 0
  DEFAULT, N_it, data.N_it-1
  DEFAULT, N_it2, N_it

  IF KEYWORD_SET(rho) THEN BEGIN
      x_exp = data.exp_rho
      x = data.rho
      xtitle = 'norm. toroidal flux !4q!X'
      IF (d2flag) THEN x2 = data2.rho
  ENDIF ELSE BEGIN
      x_exp = data.exp_rmin
      x = data.rmin
      xtitle = 'r!Dmin!N/a'
      IF (d2flag) THEN x2 = data2.rmin
  ENDELSE

  Qi_target = data.Qi_target[*,N_it]
  Qi_tot = data.Qi_tot[*,N_it]
  Qi_neo = data.Qi_neo[*,N_it]
  Qi0 = data.Qi_target[*,0]
  Qe_target = data.Qe_target[*,N_it]
  Qe_tot = data.Qe_tot[*,N_it]
  Qe_neo = data.Qe_neo[*,N_it]
  Qe0 = data.Qe_target[*,0]
  Ge_target = data.Ge_target[*,N_it]
  Ge_tot = data.Ge_tot[*,N_it]
  Ge_neo = data.Ge_neo[*,N_it]
  Mflux_target = data.Mflux_target[*,N_it]
  Mflux_tot = data.Mflux_tot[*,N_it]
  Mflux_neo = data.Mflux_i_neo[*,N_it] + data.Mflux_e_neo[*,N_it]
  Mflux0 = data.Mflux_target[*,N_it]
  units = '(gB)'
  ge_units = ' (gB)'
  mflux_units = ' (gB)'
  x0 = 0.2
  d2norm = 1.
  d2gnorm = 1.
  d2mnorm = 1.
  IF KEYWORD_SET(MKS) THEN BEGIN
      Qi_target *= data.Q_gB[*,N_it]
      Qi_tot *= data.Q_gB[*,N_it]
      Qi0 *= data.Q_gB[*,0]
      Qi_neo *= data.Q_gB[*,N_it]
      Qe_target *= data.Q_gB[*,N_it]
      Qe_tot *= data.Q_gB[*,N_it]
      Qe0 *= data.Q_gB[*,0]
      Qe_neo *= data.Q_gB[*,N_it]
      Ge_target *= data.G_gB[*,N_it]
      Ge_tot *= data.G_gB[*,N_it]
      Ge_neo *= data.G_gB[*,N_it]
      Mflux_target *= data.Pi_gB[*,N_it]
      Mflux_tot *= data.Pi_gB[*,N_it]
      Mflux_neo *= data.Pi_gB[*,N_it]
      Mflux0 *= data.Pi_gB[*,0]
      units = '(MW/m!U2!N)'
      ge_units = ' (10!U19!N/m!U2!N s)'
      Mflux_units = ' (J/m!U2!N)'
      IF (d2flag) THEN d2norm = data2.Q_gB[*,N_it2]
      IF (d2flag) THEN d2gnorm = data2.G_gB[*,N_it2]
      IF (d2flag) THEN d2mnorm = data2.Pi_gB[*,N_it2]
      x0 = 0.5
  ENDIF

  !P.MULTI = [0,2,4]

  ymax = MAX(data.exp_ti) > MAX(data.ti[*,N_it])
  PLOT, x_exp, data.exp_ti, XRANGE = [0,1], TITLE = data.simdir, $
        THICK=thick,XTHICK=thick,YTHICK=thick,CHARTHICK=thick, $
        CHARSIZE=cs, YRANGE=[0,ymax]
  IF (d2flag) THEN OPLOT, x2, data2.ti[*,N_it2], COLOR=150, PSYM=-2,THICK=thick
  OPLOT, x, data.ti[*,N_it], COLOR=100,PSYM=-4, THICK=thick
  XYOUTS, 0.6, 0.8*ymax, 'T!Di!N (keV)',$
          CHARSIZE=2,CHARTHICK=thick
  IF (d2flag) THEN XYOUTS, 0.5, ymax, data2.simdir, ALIGN=0.5, $
    COLOR=150, CHARTHICK=thick, CHARSIZE=cs

  ymax = MAX(Qi_target) > MAX(Qi_tot)
  IF KEYWORD_SET(mks) THEN BEGIN
     ymin = 0
     ylog = 0
  ENDIF ELSE BEGIN
     ymin = 0.1*(MIN(Qi_target[1:*]) < MIN(Qi_tot[1:*]))
     ylog = 1
  ENDELSE
  PLOT, x, Qi_target, XRANGE=[0,1], LINESTYLE=2, YLOG = ylog, $
        TITLE = 'B/W: iteration 0, red: TGYRO', YRANGE=[ymin,ymax], $
        CHARSIZE=cs, $
        THICK=thick,XTHICK=thick,YTHICK=thick,CHARTHICK=thick
  IF (d2flag) THEN OPLOT, x2, data2.Qi_tot[*,N_it2]*d2norm, COLOR=150, $
    PSYM=-2,THICK=thick
  OPLOT, x, Qi0, THICK=thick
  OPLOT, x, Qi_tot, COLOR=100,PSYM=-4, THICK=thick
  OPLOT, x, Qi_neo, COLOR=50, PSYM=-4, THICK=thick
  XYOUTS, x0, 0.8*ymax, 'Q!Di!N ' + units,$
          CHARSIZE=2,CHARTHICK=thick

  ymax = MAX(data.exp_te) > MAX(data.te[*,N_it])
  PLOT, x_exp, data.exp_te, XRANGE = [0,1], CHARSIZE=cs, $
        THICK=thick,XTHICK=thick,YTHICK=thick,CHARTHICK=thick, $
        TITLE = 'iteration # ' + NUMTOSTRING(N_it), YRANGE=[0,ymax]
  IF (d2flag) THEN OPLOT, x2, data2.te[*,N_it2], COLOR=150, PSYM=-2,THICK=thick
  OPLOT, x, data.te[*,N_it], COLOR=100,PSYM=-4, THICK=thick
  XYOUTS, 0.6, 0.8*ymax, 'T!De!N (keV)',$
          CHARSIZE=2,CHARTHICK=thick
  IF (d2flag) THEN XYOUTS, 0.5, ymax, 'iteration # ' + NUMTOSTRING(N_it2), $
    ALIGN=0.5, COLOR=150, CHARTHICK=thick, CHARSIZE=cs

  ymax = MAX(Qe_target) > MAX(Qe_tot)
  IF KEYWORD_SET(mks) THEN BEGIN
     ymin = 0
     ylog = 0
  ENDIF ELSE BEGIN
     ymin = 0.1*(MIN(Qe_target[1:*]) < MIN(Qe_tot[1:*]))
     ylog = 1
  ENDELSE
  PLOT, x, Qe_target, XRANGE=[0,1], XTITLE=xtitle, LINESTYLE=2, YLOG=ylog, $
        TITLE = 'dashed: target, blue: neo', YRANGE=[ymin,ymax], $
        CHARSIZE=cs, $
        THICK=thick,XTHICK=thick,YTHICK=thick,CHARTHICK=thick
  IF (d2flag) THEN OPLOT, x2, data2.Qe_tot[*,N_it2]*d2norm, COLOR=150, $
    PSYM=-2,THICK=thick
  OPLOT, x, Qe0, THICK=thick
  OPLOT, x, Qe_tot, COLOR=100,PSYM=-4, THICK=thick
  OPLOT, x, Qe_neo, COLOR=50, PSYM=-4, THICK=thick
  XYOUTS, x0, 0.8*ymax, 'Q!De!N ' + units,$
          CHARSIZE=2,CHARTHICK=thick

  ymin = MIN(data.exp_ne) < MIN(data.n_e[*,N_it])
  ymax = MAX(data.exp_ne) > MAX(data.n_e[*,N_it])
  PLOT, x_exp, data.exp_ne, XRANGE = [0,1], $
        THICK=thick,XTHICK=thick,YTHICK=thick,CHARTHICK=thick, $
        CHARSIZE=cs, YRANGE=[ymin,ymax]
  IF (d2flag) THEN OPLOT, x2, data2.n_e[*,N_it2], COLOR=150, PSYM=-2,$
                          THICK=thick
  OPLOT, x, data.n_e[*,N_it], COLOR=100,PSYM=-4, THICK=thick
  XYOUTS, 0.4, 0.8*ymax, 'n!De!N (10!U19!N/m!U3!N)',$
          CHARSIZE=2,CHARTHICK=thick
  
  ymax = MAX(ABS(Ge_tot)) > MAX(Ge_target)
  PLOT, x, Ge_target, XRANGE=[0,1], LINESTYLE=2, $
        YRANGE=[-ymax,ymax], CHARSIZE=cs, $
        THICK=thick,XTHICK=thick,YTHICK=thick,CHARTHICK=thick
  OPLOT, x, Ge_tot, COLOR=100,PSYM=-4, THICK=thick
  OPLOT, x, Ge_neo, COLOR=50, PSYM=-4, THICK=thick
  IF (d2flag) THEN OPLOT, x2, data2.ge_tot[*,N_it2]*d2gnorm, $
                          COLOR=150, PSYM=-4, THICK=thick
  XYOUTS, 0.4, 0.7*ymax, '!4C!X!De!N ' + ge_units,$
          CHARSIZE=2,CHARTHICK=thick

  rot0 = data.M[*,0]*data.c_s[*,0]/1e5
  rot = data.M[*,N_it]*data.c_s[*,N_it]/1e5
  ymin = MIN(rot0) < MIN(rot)
  ymax = MAX(rot0) > MAX(rot)
  PLOT, x, rot0, XRANGE = [0,1], $
        THICK=thick,XTHICK=thick,YTHICK=thick,CHARTHICK=thick, $
        CHARSIZE=cs, XTITLE=xtitle, YRANGE=[ymin,ymax]
  OPLOT, x, rot, COLOR=100, PSYM=-4, THICK=thick
  IF (d2flag) THEN OPLOT, x2, data2.M[*,N_it2]*data2.c_s[*,N_it2]/1e5, COLOR=150, $
                          PSYM=-2, THICK=thick
  XYOUTS, 0.4, 0.8*ymax, 'R!4x!X (10!U5!N m/s)', CHARSIZE=2,CHARTHICK=thick

  ymax = MAX(Mflux_target) > MAX(Mflux_tot)
  IF KEYWORD_SET(mks) THEN ymin = 0 ELSE $
     ymin = 0.1*(MIN(Mflux_target[1:*]) < MIN(Mflux_tot[1:*]))

  PLOT, x, Mflux_target, XRANGE=[0,1], XTITLE=xtitle, LINESTYLE=2, $
        YRANGE=[ymin,ymax], CHARSIZE=cs, $
        THICK=thick,XTHICK=thick,YTHICK=thick,CHARTHICK=thick
  OPLOT, x, Mflux0, THICK=thick
  OPLOT, x, Mflux_tot, COLOR=100,PSYM=-4, THICK=thick
  OPLOT, x, Mflux_neo, COLOR=50, PSYM=-4, THICK=thick
  IF (d2flag) THEN OPLOT, x2, data2.mflux_tot[*,N_it2]*d2mnorm, $
                          COLOR=150, PSYM=-4, THICK=thick
  XYOUTS, x0, 0.8*ymax, '!4P!X' + mflux_units,$
          CHARSIZE=2,CHARTHICK=thick

  !P.MULTI=0

  IF KEYWORD_SET(summary) THEN BEGIN
      SET_PLOT, 'PS'
      DEVICE, XS=30,YS=30, /ENCAPS, /COLOR, BITS=8, $
              FILE = GETENV('TGYRO_DIR') + '/sim/' + simdir + '/summary.eps'
      PLOT_TGYRO_LOC_SUMMARY, simdir, N_it = N_it, RHO=rho, MKS=mks,$
        PS = 1, DATA2 = data2, Nit2=N_it2, PLOT_N=plot_n, _EXTRA=extra
      DEVICE, /CLOSE
      SET_PLOT, 'X'
  ENDIF
END ;plot_tgyro_loc_summary
