PRO plot_gyro_linfreq, simdir, ITMIN=itmin, ITMAX=itmax, FILENAME=filename

  data = GET_GYRO_DATA(simdir)
  NN = data.n_n
  N_TIME = data.n_time

  arr = FLTARR(4,NN,N_TIME)
  OPENR, 1, GETENV('GYRO_DIR') + '/sim/' + simdir + '/out.gyro.freq'
  READF, 1, arr
  CLOSE, 1
  kt = data.ktheta
  freq = FLTARR(NN)
  freq_stddev = FLTARR(NN)
  gamma = FLTARR(NN)
  gamma_stddev = FLTARR(NN)

  DEFAULT, itmin, data.n_time/2
  DEFAULT, itmax, data.n_time-1
  FOR ii=0,NN-1 DO BEGIN
	freq[ii] = MEAN(arr[0,ii,itmin:itmax])
	freq_stddev[ii] = STDDEV(arr[0,ii,itmin:itmax])
	gamma[ii] = MEAN(arr[1,ii,itmin:itmax])
	gamma_stddev[ii] = STDDEV(arr[1,ii,itmin:itmax])
  ENDFOR

  !P.MULTI=[0,0,2]
  PLOT, kt, freq, psym=4, YTITLE = 'omega (a/c_s)', $
	TITLE='it=['+NUMTOSTRING(itmin)+':'+NUMTOSTRING(itmax)+']'
  ERRPLOT, kt, freq-freq_stddev,freq+freq_stddev
  OPLOT, kt, 0.*kt, LINESTYLE=2
  PLOT, kt, gamma, psym=4, YTITLE = 'gamma (a/c_s)'
  ERRPLOT, kt, gamma-gamma_stddev,gamma+gamma_stddev
  OPLOT, kt, 0.*kt, LINESTYLE=2
  !P.MULTI=0

  IF N_ELEMENTS(filename) NE 0 THEN BEGIN
	OPENW, 1, GETENV('GYRO_DIR') + '/sim/' + simdir + '/' + filename
	PRINTF, 1, '#ktrho_s   omega   omega_stddev    gamma   gamma_stddev'
	FOR ii=0,nn-1 DO PRINTF, 1, kt[ii], freq[ii], freq_stddev[ii], $
		gamma[ii], gamma_stddev[ii]
	CLOSE, 1
  ENDIF
END ;plot_gyro_linfreq
