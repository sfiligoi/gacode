PRO gyro_quicksum, simdir, ITMIN=itmin, ITMAX=itmax, $
                   PLOT_PHIRMS = plot_phirms, $
                   PRINT_FRACDIFF=print_fracdiff, PLOT_FRACDIFF=plot_fracdiff
;C. Holland, UCSD
;v1.0: 8.5.09: quick summary tool plotting+printing avg. diff and diff-n.  

  data = GET_GYRO_DATA(simdir, READ_LARGE=5)
  DEFAULT, itmin, data.n_time/2
  DEFAULT, itmax, data.n_time-1

  n_kinetic = data.profile_data.n_kinetic
  avgchi = FLTARR(data.profile_data.n_kinetic)
  FOR ii=0, n_kinetic-1 DO $
         avgchi[ii] = MEAN(TOTAL(data.chi[ii,*,itmin:itmax],2))
  avgDe = MEAN(TOTAL(data.D[n_kinetic-1,*,itmin:itmax],2))

  PRINT, 'Summary for ' + data.simdir
  PRINT, 'itmin: ', itmin
  PRINT, 'itmax: ', itmax
  PRINT, 'chi_i (gB) (sim, expt): ', avgchi[0], $
         data.profile_Data.chi_i_exp[data.n_r/2]
  PRINT, 'chi_e (gB) (sim, expt): ', avgchi[n_kinetic-1], $
         data.profile_Data.chi_e_exp[data.n_r/2]
  PRINT, 'Dn_e (gB)  (sim, expt): ', avgDe, $
         data.profile_Data.diff_ne_exp[data.n_r/2]

  ;diff-n
  n_n = data.n_n
  avgchi_kt = FLTARR(3,n_n)
  ie = data.profile_data.n_kinetic-1
  FOR i_n = 0, n_n-1 DO BEGIN
      avgchi_kt[0,i_n] = MEAN(TOTAL(data.profile_data.chi_kt[0,*,i_n,$
                                                               itmin:itmax],1))
      avgchi_kt[1,i_n] = MEAN(TOTAL(data.profile_data.chi_kt[ie,*,i_n,$
                                                               itmin:itmax],1))
      avgchi_kt[2,i_n] = MEAN(TOTAL(data.profile_data.D_kt[ie,*,i_n,$
                                                               itmin:itmax],1))
  ENDFOR
  avgchi_kt[0,*] /= avgchi[0]
  avgchi_kt[1,*] /= avgchi[1]
  avgchi_kt[2,*] /= avgDe
  PRINT, 'fractional diffusivities vs ky rho_s'
  IF KEYWORD_SET(print_fracdiff) THEN BEGIN
      PRINT, '      ky rho_s     chi_i        chi_e        Dne' 
      FOR ii=0,n_n-1 DO PRINT, data.profile_data.ktheta[ii], avgchi_kt[0,ii], $
        avgchi_kt[1,ii], avgchi_kt[2,ii]
  ENDIF ELSE PRINT, '                   chi_i        chi_e        Dne'
  chimax = MAX(avgchi_kt[0,*],idx0)
  chimax = MAX(avgchi_kt[1,*],idx1)
  chimax = MAX(ABS(avgchi_kt[2,*]),idx2)
  PRINT, 'peak ky:  ', data.profile_Data.ktheta[idx0], $
         data.profile_Data.ktheta[idx1], data.profile_Data.ktheta[idx2]
  PRINT, 'max ky: ', data.profile_data.ktheta[n_n-1]
  IF KEYWORD_SET(plot_fracdiff) THEN BEGIN
      yrange = [MIN(avgchi_kt) < 0, MAX(avgchi_kt)]
      PLOT, data.profile_data.ktheta, avgchi_kt[0,*], PSYM=-4,$
            XTITLE='ky r_s', YTITLE='frac. diffusivities', /XS, $
            TITLE='chi_i: B/W, chi_e: red, Dne: blue', YRANGE=yrange
      OPLOT, data.profile_data.ktheta, FLTARR(n_n), LINESTYLE=2
      OPLOT, data.profile_data.ktheta, avgchi_kt[1,*], PSYM=-2, COLOR=100
      OPLOT, data.profile_data.ktheta, avgchi_kt[2,*], PSYM=-1, COLOR=50
  ENDIF

  IF KEYWORD_SET(plot_phirms) THEN BEGIN
      kxkyspec = FLTARR(data.n_r, n_n, data.n_time)
      OPENR, lun, GETENV('GYRO_DIR') + '/sim/' + data.simdir+'/kxkyspec.out', $
             /GET_LUN
      READF, lun, kxkyspec
      FREE_LUN, lun
  
      phisq = FLTARR(n_n, data.n_time)
      FOR ii=0,n_n-1 DO phisq[ii,*] = TOTAL(kxkyspec[*,ii,*],1)
      y = FLTARR(2,data.n_time)
      y[0,*] = SQRT(phisq[0,*])
      y[1,*] = SQRT(2*TOTAL(phisq[1:*,*],1))
      zonalmean = MEAN(y[0,itmin:itmax])
      turbmean = MEAN(y[1,itmin:itmax])

      tfit = data.t[itmin:itmax]
      PLOT, data.t, y[0,*], /XS, XTITLE='time', $
            TITLE='RMS n=0 (B/W) and n>0 phi (red)', $
            YRANGE=[0, MAX(y)]
      OPLOT, data.t, y[1,*], COLOR=100
      OPLOT, tfit, FLTARR(N_ELEMENTS(tfit)) + zonalmean, COLOR=50
      OPLOT, tfit, FLTARR(N_ELEMENTS(tfit)) + turbmean, COLOR=150

      res = LINFIT(tfit,y[0,itmin:itmax],YFIT=zonalfit)
      PRINT, 'zonal % change: ', (MAX(zonalfit)-MIN(zonalfit))/zonalmean
      res = LINFIT(tfit,y[1,itmin:itmax],YFIT=turbfit)
      PRINT, 'turb % change: ', (MAX(turbfit)-MIN(turbfit))/turbmean
      OPLOT, tfit, zonalfit, COLOR=50, LINESTYLE=2
      OPLOT, tfit, turbfit, COLOR=150, LINESTYLE=2
  ENDIF
END ;gyro_quicksum
