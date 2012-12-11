FUNCTION CALCULATE_GYRO_TIME_ENS, y, itmin, itmax, itavg, T_ENS = t_ens
;
; Given real timeseires y(t), break up into subwindows itavg long between itmin and itmax
; return y_ens = ensemble of averages
; T_ENS is starting points of each averaging window used for y_ens
  z = y[itmin:itmax]
  N_T = itmax-itmin+1
  NM = FIX(N_T/itavg)

  y_ens = FLTARR(NM)
  t_ens = itmin + itavg*FINDGEN(NM)
  FOR im=0,NM-1 DO y_ens[im] = TOTAL(z[im*itavg:(im+1)*itavg-1])/itavg
  RETURN, y_ens

END ;calculate_gyro_time_ens

PRO PLOT_GYRO_SUMMARY, simdir, ITMIN=itmin, ITMAX=itmax, ITAVG=itavg, USE_LX=use_lx,$
	THICK=th, CHARSIZE=cs
; USE_LX keyword plots radial profiles vs Dx/rhos instead of rmin/a
; good EPS versions via
; IDL> SET_PLOT, 'PS'
; IDL> DEVICE, XS=30, YS=30, /encaps, /color, bits=8,file='<filename>.eps'
; IDL> PLOT_GYRO_SUMMARY, <simdir>, THICK=6, CHARSIZE=3, ...
; IDL> DEVICE, /CLOSE
; IDL> SET_PLOT, 'X'

  data = GET_GYRO_DATA(simdir)
  Q_t = TOTAL(data.q_t,2)
  Qi_t = TOTAL(Q_t[0:data.n_ion-1,*],1)
  IF (data.n_kinetic EQ 1) THEN BEGIN
	Qe_t = data.t*0. 
	Ge_t = data.t*0. 
  ENDIF ELSE BEGIN
	Qe_t = REFORM(q_t[data.n_kinetic-1,*])
	Ge_t = REFORM(TOTAL(data.gamma_t[data.n_kinetic-1,*,*],2))
  ENDELSE
  Pi_t = TOTAL(TOTAL(data.pi_t[0:data.n_ion-1,*,*],2),1)
  ymin = MIN(Q_t) < 1.5*MIN(Ge_t) < MIN(Pi_t) < 0
  ymax = MAX(Q_t) > 1.5*MAX(Ge_t) > MAX(Pi_t) > 0

  DEFAULT, itmin, data.n_time/2
  DEFAULT, itmax, data.n_time-1
  DEFAULT, itavg, 1
  
  dcolor = !D.TABLE_SIZE/5
  DEFAULT, th, 1
  DEFAULT, cs, 2
  !P.MULTI=[0,0,3]

  PLOT, data.t, Qi_t, YRANGE=[ymin, ymax], /XS, $
	XTITLE = 'time (a/c!Ds!N)', CHARSIZE=cs, TITLE=simdir, $
	THICK=th, XTHICK=th, YTHICK=th, CHARTHICK=th
  Qi_t_ens = CALCULATE_GYRO_TIME_ENS(Qi_t,itmin,itmax,itavg,T_ENS=t_ens)
  NM = N_ELEMENTS(t_ens)
  OPLOT, [t_ens[0], t_ens[NM-1]+itavg-1], [1,1]*MEAN(Qi_t_ens), THICK=th
  XYOUTS, t_ens[0], 1.1*MEAN(Qi_t_ens), 'Q!Di!N/Q!DgB!N', $
	CHARSIZE=2, CHARTHICK=th
  OPLOT, data.t, data.t*0, LINESTYLE=1, THICK=th

  OPLOT, data.t, Pi_t, COLOR=dcolor, THICK=th
  Pi_t_ens = CALCULATE_GYRO_TIME_ENS(Pi_t,itmin,itmax,itavg)
  OPLOT, [t_ens[0], t_ens[NM-1]+itavg-1], [1,1]*MEAN(Pi_t_ens), $
	COLOR=dcolor, THICK=th
  XYOUTS, t_ens[0], 1.1*MEAN(Pi_t_ens), '!4P!X!Di!N/!4P!X!DgB!N', $
	COLOR=dcolor, CHARSIZE=2, CHARTHICK=th

  IF (data.n_kinetic GT data.n_ion) THEN BEGIN
	OPLOT, data.t, Qe_t, COLOR=2*dcolor, THICK=th
  	Qe_t_ens = CALCULATE_GYRO_TIME_ENS(Qe_t,itmin,itmax,itavg)
  	OPLOT, [t_ens[0], t_ens[NM-1]+itavg-1], [1,1]*MEAN(Qe_t_ens), $
		COLOR=2*dcolor, THICK=th
	XYOUTS, t_ens[0], 1.1*MEAN(Qe_t_ens), 'Q!De!N/Q!DgB!N', $
		COLOR=2*dcolor, CHARSIZE=2, CHARTHICK=th

	OPLOT, data.t, 1.5*Ge_t, COLOR=3*dcolor, THICK=th
  	Ge_t_ens = CALCULATE_GYRO_TIME_ENS(Ge_t,itmin,itmax,itavg)
  	OPLOT, [t_ens[0], t_ens[NM-1]+itavg-1], [1,1]*1.5*MEAN(Ge_t_ens), $
		COLOR=3*dcolor, THICK=th
	XYOUTS, t_ens[0], 1.1*1.5*MEAN(Ge_t_ens), '1.5!4C!X!De!N/!4C!X!DgB!N', $
		COLOR=3*dcolor, CHARSIZE=2, CHARTHICK=th
  ENDIF ELSE BEGIN
	Qe_t_ens = FLTARR(NM)
	Ge_t_ens = FLTARR(NM)
  ENDELSE

  PRINT, 'NM: ', NM
  PRINT, 'it=['+NUMTOSTRING(itmin)+':'+NUMTOSTRING(t_ens[NM-1]+itavg-1)+']'
  PRINT, 'Q_ion: ', MEAN(Qi_t_ens), '+/-', SDOM(Qi_t_ens)
  PRINT, 'Q_e: ', MEAN(Qe_t_ens), '+/-', SDOM(Qe_t_ens)
  PRINT, 'Gamma_e: ', MEAN(Ge_t_ens), '+/-', SDOM(Ge_t_ens)
  PRINT, 'Pi_ion: ', MEAN(Pi_t_ens), '+/-', SDOM(Pi_t_ens)

  ;spectra
  Qi_kt = TOTAL(TOTAL(data.q_nt[0:data.n_ion-1,*,*,*],2),1)
  Pi_kt = TOTAL(TOTAL(data.pi_nt[0:data.n_ion-1,*,*,*],2),1)
  Qe_kt = REFORM(TOTAL(data.q_nt[data.n_kinetic-1,*,*,*],2))
  Ge_kt = REFORM(TOTAL(data.gamma_nt[data.n_kinetic-1,*,*,*],2))
  Qi_kt_ens = FLTARR(data.n_n,NM)
  Pi_kt_ens = FLTARR(data.n_n,NM)
  Qe_kt_ens = FLTARR(data.n_n,NM)
  Ge_kt_ens = FLTARR(data.n_n,NM)
  Qi_kt_mean = FLTARR(data.n_n)
  Pi_kt_mean = FLTARR(data.n_n)
  Qe_kt_mean = FLTARR(data.n_n)
  Ge_kt_mean = FLTARR(data.n_n)
  Qi_kt_sdom = FLTARR(data.n_n)
  Pi_kt_sdom = FLTARR(data.n_n)
  Qe_kt_sdom = FLTARR(data.n_n)
  Ge_kt_sdom = FLTARR(data.n_n)
  IF (data.n_kinetic EQ data.n_ion) THEN BEGIN
	Qe_kt[*,*] = 0.
	Ge_kt[*,*] = 0.
  ENDIF

  FOR i_n = 0, data.n_n-1 DO BEGIN
	Qi_kt_ens[i_n,*] = CALCULATE_GYRO_TIME_ENS(REFORM(Qi_kt[i_n,*]),itmin,itmax,itavg)
	Qi_kt_mean[i_n] = MEAN(Qi_kt_ens[i_n,*])
	Qi_kt_sdom[i_n] = SDOM(Qi_kt_ens[i_n,*])
	Pi_kt_ens[i_n,*] = CALCULATE_GYRO_TIME_ENS(REFORM(Pi_kt[i_n,*]),itmin,itmax,itavg)
	Pi_kt_mean[i_n] = MEAN(Pi_kt_ens[i_n,*])
	Pi_kt_sdom[i_n] = SDOM(Pi_kt_ens[i_n,*])
	Qe_kt_ens[i_n,*] = CALCULATE_GYRO_TIME_ENS(REFORM(Qe_kt[i_n,*]),itmin,itmax,itavg)
	Qe_kt_mean[i_n] = MEAN(Qe_kt_ens[i_n,*])
	Qe_kt_sdom[i_n] = SDOM(Qe_kt_ens[i_n,*])
	Ge_kt_ens[i_n,*] = CALCULATE_GYRO_TIME_ENS(REFORM(Ge_kt[i_n,*]),itmin,itmax,itavg)
	Ge_kt_mean[i_n] = MEAN(Ge_kt_ens[i_n,*])
	Ge_kt_sdom[i_n] = SDOM(Ge_kt_ens[i_n,*])
  ENDFOR

  ymin = MIN(Qi_kt_mean) < MIN(Qe_kt_mean) < MIN(Pi_kt_mean) < 1.5*MIN(Ge_kt_mean) < 0
  ymax = MAX(Qi_kt_mean) > MAX(Qe_kt_mean) > MAX(Pi_kt_mean) > 1.5*MAX(Ge_kt_mean) > 0

  PLOT, data.ktheta, Qi_kt_mean, PSYM=-4, XTITLE = 'k!D!4h!Nq!X!Ds!N', CHARSIZE=cs, $
	YRANGE=[ymin, ymax], /XS, $
	TITLE='t=['+NUMTOSTRING(data.t[itmin])+':'+NUMTOSTRING(data.t[t_ens[NM-1]+itavg-1])+'], itavg ='+$
	NUMTOSTRING(itavg) + ', NM = ' +NUMTOSTRING(NM), $
        THICK=th, XTHICK=th, YTHICK=th, CHARTHICK=th
  ERRPLOT, data.ktheta, Qi_kt_mean-Qi_kt_sdom, Qi_kt_mean+Qi_kt_sdom, THICK=th
  OPLOT, data.t, data.t*0, LINESTYLE=1, THICK=th

  OPLOT, data.ktheta, Pi_kt_mean, PSYM=-2, COLOR=dcolor, THICK=th
  ERRPLOT, data.ktheta, Pi_kt_mean-Pi_kt_sdom, Pi_kt_mean+Pi_kt_sdom, COLOR=dcolor, THICK=th
  OPLOT, data.ktheta, Qe_kt_mean, PSYM=-1, COLOR=2*dcolor, THICK=th
  ERRPLOT, data.ktheta, Qe_kt_mean-Qe_kt_sdom, Qe_kt_mean+Qe_kt_sdom, COLOR=2*dcolor, THICK=th
  OPLOT, data.ktheta, 1.5*Ge_kt_mean, PSYM=-1, COLOR=3*dcolor, THICK=th
  ERRPLOT, data.ktheta, 1.5*(Ge_kt_mean-Ge_kt_sdom), 1.5*(ge_kt_mean+Ge_kt_sdom), COLOR=3*dcolor, THICK=th

  XYOUTS, data.ktheta[data.n_n/2], 0.95*ymax, 'Q!Di!N/Q!DgB!N: ' + NUMTOSTRING(MEAN(Qi_t_ens)) + '!9 + !X' + $
	NUMTOSTRING(SDOM(Qi_t_ens)), CHARTHICK=th, CHARSIZE=2
  XYOUTS, data.ktheta[data.n_n/2], 0.75*ymax, 'Q!De!N/Q!DgB!N: ' + NUMTOSTRING(MEAN(Qe_t_ens)) + '!9 + !X' + $
	NUMTOSTRING(SDOM(Qe_t_ens)), CHARTHICK=th, CHARSIZE=2, COLOR=2*dcolor
  XYOUTS, data.ktheta[data.n_n/2], 0.55*ymax, '!4C!X!De!N/!4C!X!DgB!N: ' + NUMTOSTRING(MEAN(Ge_t_ens)) + '!9 + !X' + $
	NUMTOSTRING(SDOM(Ge_t_ens)), CHARTHICK=th, CHARSIZE=2, COLOR=3*dcolor
  XYOUTS, data.ktheta[data.n_n/2], 0.35*ymax, '!4P!X!Di!N/!4P!X!DgB!N: ' + NUMTOSTRING(MEAN(Pi_t_ens)) + '!9 + !X' + $
	NUMTOSTRING(SDOM(Pi_t_ens)), CHARTHICK=th, CHARSIZE=2, COLOR=dcolor

  ;profiles
  Qi_rt = TOTAL(TOTAL(data.q_rt[0:data.n_ion-1,*,*,*],2),1)
  Pi_rt = TOTAL(TOTAL(data.pi_rt[0:data.n_ion-1,*,*,*],2),1)
  Qe_rt = REFORM(TOTAL(data.q_rt[data.n_kinetic-1,*,*,*],2))
  Ge_rt = REFORM(TOTAL(data.gamma_rt[data.n_kinetic-1,*,*,*],2))
  Qi_rt_ens = FLTARR(data.n_r,NM)
  Pi_rt_ens = FLTARR(data.n_r,NM)
  Qe_rt_ens = FLTARR(data.n_r,NM)
  Ge_rt_ens = FLTARR(data.n_r,NM)
  Qi_rt_mean = FLTARR(data.n_r)
  Pi_rt_mean = FLTARR(data.n_r)
  Qe_rt_mean = FLTARR(data.n_r)
  Ge_rt_mean = FLTARR(data.n_r)
  Qi_rt_sdom = FLTARR(data.n_r)
  Pi_rt_sdom = FLTARR(data.n_r)
  Qe_rt_sdom = FLTARR(data.n_r)
  Ge_rt_sdom = FLTARR(data.n_r)
  IF (data.n_kinetic EQ data.n_ion) THEN BEGIN
	Qe_rt[*,*] = 0.
	Ge_rt[*,*] = 0.
  ENDIF

  FOR i_r = 0, data.n_r-1 DO BEGIN
	Qi_rt_ens[i_r,*] = CALCULATE_GYRO_TIME_ENS(REFORM(Qi_rt[i_r,*]),itmin,itmax,itavg)
	Qi_rt_mean[i_r] = MEAN(Qi_rt_ens[i_r,*])
	Qi_rt_sdom[i_r] = SDOM(Qi_rt_ens[i_r,*])
	Pi_rt_ens[i_r,*] = CALCULATE_GYRO_TIME_ENS(REFORM(Pi_rt[i_r,*]),itmin,itmax,itavg)
	Pi_rt_mean[i_r] = MEAN(Pi_rt_ens[i_r,*])
	Pi_rt_sdom[i_r] = SDOM(Pi_rt_ens[i_r,*])
	Qe_rt_ens[i_r,*] = CALCULATE_GYRO_TIME_ENS(REFORM(Qe_rt[i_r,*]),itmin,itmax,itavg)
	Qe_rt_mean[i_r] = MEAN(Qe_rt_ens[i_r,*])
	Qe_rt_sdom[i_r] = SDOM(Qe_rt_ens[i_r,*])
	Ge_rt_ens[i_r,*] = CALCULATE_GYRO_TIME_ENS(REFORM(Ge_rt[i_r,*]),itmin,itmax,itavg)
	Ge_rt_mean[i_r] = MEAN(Ge_rt_ens[i_r,*])
	Ge_rt_sdom[i_r] = SDOM(Ge_rt_ens[i_r,*])
  ENDFOR

  ymin = MIN(Qi_rt_mean) < MIN(Qe_rt_mean) < MIN(Pi_rt_mean) < 1.5*MIN(Ge_rt_mean) < 0
  ymax = MAX(Qi_rt_mean) > MAX(Qe_rt_mean) > MAX(Pi_rt_mean) > 1.5*MAX(Ge_rt_mean) > 0

  x = data.r
  xtitle = 'r!Dmin!N/a'
  IF KEYWORD_SET(USE_LX) THEN BEGIN
	x = (data.r - data.r[data.n_r/2])/data.rho_s
	xtitle = '!4D!Xx/!4q!X!Ds!N'
  ENDIF

  PLOT, x, Qi_rt_mean, CHARSIZE=cs, YRANGE=[ymin,ymax], /XS, XTITLE=xtitle, $
        THICK=th, XTHICK=th, YTHICK=th, CHARTHICK=th
  OPLOT, x, Qi_rt_mean + Qi_rt_sdom, LINE=2, THICK=th
  OPLOT, x, Qi_rt_mean - Qi_rt_sdom, LINE=2, THICK=th
  OPLOT, [x[data.n_bnd],x[data.n_r-1-data.n_bnd]],[1,1]*MEAN(Qi_t_ens), THICK=th

  OPLOT, x, x*0, linestyle=1, THICK=th
  OPLOT, [1,1]*x[data.n_bnd], 10.*[ymin,ymax], linestyle=2, THICK=th
  OPLOT, [1,1]*x[data.n_r-1-data.n_bnd], 10.*[ymin,ymax], linestyle=2, THICK=th

  OPLOT, x, Pi_rt_mean, COLOR=dcolor, THICK=th
  OPLOT, x, Pi_rt_mean + Pi_rt_sdom, LINE=2, COLOR=dcolor, THICK=th
  OPLOT, x, pi_rt_mean - Pi_rt_sdom, LINE=2, COLOR=dcolor, THICK=th
  OPLOT, [x[data.n_bnd],x[data.n_r-1-data.n_bnd]],[1,1]*MEAN(Pi_t_ens), $
	COLOR=dcolor, THICK=th

  OPLOT, x, Qe_rt_mean, COLOR=2*dcolor, THICK=th
  OPLOT, x, Qe_rt_mean + Qe_rt_sdom, LINE=2, COLOR=2*dcolor, THICK=th
  OPLOT, x, Qe_rt_mean - Qe_rt_sdom, LINE=2, COLOR=2*dcolor, THICK=th
  OPLOT, [x[data.n_bnd],x[data.n_r-1-data.n_bnd]],[1,1]*MEAN(Qe_t_ens), $
	COLOR=2*dcolor, THICK=th

  OPLOT, x, 1.5*Ge_rt_mean, COLOR=3*dcolor, THICK=th
  OPLOT, x, 1.5*(Ge_rt_mean + Ge_rt_sdom), LINE=2, COLOR=3*dcolor, THICK=th
  OPLOT, x, 1.5*(Ge_rt_mean - Ge_rt_sdom), LINE=2, COLOR=3*dcolor, THICK=th
  OPLOT, [x[data.n_bnd],x[data.n_r-1-data.n_bnd]],[1,1]*1.5*MEAN(Ge_t_ens), $
	COLOR=3*dcolor, THICK=th
  !P.MULTI=0
END ;PLOT_GYRO_SUMMARY
