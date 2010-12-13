PRO write_tgyro2gyro_input, data, ii, IT=it
;
; C. Holland, UCSD
; v1.0 10.18.09
; 
; writes an GYRO INPUT file using parameters for a local TGYRO iteration.  
; To use, 1st call in IDL
;
; IDL> data = GET_TGYRO_LOC_DATA('simdir')
;
; where simdir is the name of the tgyro simulation you want.  You then call
;
; IDL> WRITE_TGYRO2GYRO_INPUT, data, ii
;
; where ii is the radial index of the location you want.  Routine defaults to
; using the last iteration; set IT to specify a different one
;
; a file named INPUT_ii will be written in the the TGYRO simdir directory.
;

  DEFAULT, it, data.n_it-1

  OPENW, 1, GETENV('TGYRO_DIR') + '/sim/' + data.simdir + $
	'/INPUT_'+NUMTOSTRING(ii)
  PRINTF, 1, '###PARAMETERS FROM ' + data.simdir + ' iteration ' + $
          NUMTOSTRING(it) + '###'
  PRINTF, 1, 'RADIUS='+NUMTOSTRING(data.rmin[ii])
  PRINTF, 1, 'ASPECT_RATIO='+NUMTOSTRING(data.rmaj[ii])
  PRINTF, 1, 'SHIFT='+NUMTOSTRING(data.shift[ii])
  PRINTF, 1, 'ZMAG='+NUMTOSTRING(data.z_mag[ii])
  PRINTF, 1, 'DZMAG='+NUMTOSTRING(data.dzmag[ii])
  PRINTF, 1, 'KAPPA='+NUMTOSTRING(data.kappa[ii])
  PRINTF, 1, 'S_KAPPA='+NUMTOSTRING(data.s_kappa[ii])
  PRINTF, 1, 'DELTA='+NUMTOSTRING(data.delta[ii])
  PRINTF, 1, 'S_DELTA='+NUMTOSTRING(data.s_delta[ii])
  PRINTF, 1, 'ZETA='+NUMTOSTRING(data.zeta[ii])
  PRINTF, 1, 'S_ZETA='+NUMTOSTRING(data.s_zeta[ii])
  PRINTF, 1, 'SAFETY_FACTOR='+NUMTOSTRING(data.q[ii])
  PRINTF, 1, 'SHEAR='+NUMTOSTRING(data.s_hat[ii])
  PRINTF, 1, 'RHO_STAR='+NUMTOSTRING(data.rho_star[ii,it])
  PRINTF, 1, 'Z_EFF='+NUMTOSTRING(INTERPOL(data.exp_zeff, data.exp_rmin, $
                                           data.rmin[ii]))
  PRINTF, 1, 'MACH='+NUMTOSTRING(data.M[ii,it])
  PRINTF, 1, 'PGAMMA='+NUMTOSTRING(data.gamma_p[ii,it])
  PRINTF, 1, 'GAMMA_E='+NUMTOSTRING(data.gamma_e[ii,it])

  PRINTF, 1, '# Ion  1'
  PRINTF, 1, 'NI_OVER_NE=1.0'
  PRINTF, 1, 'TI_OVER_TE='+NUMTOSTRING(data.ti[ii,it]/data.te[ii,it])
  PRINTF, 1, 'DLNNDR='+NUMTOSTRING(data.a_over_Lni[ii,it])
  PRINTF, 1, 'DLNTDR='+NUMTOSTRING(data.a_over_Lti[ii,it])

  PRINTF, 1, '# Electrons'
  PRINTF, 1, 'DLNNDR_ELECTRON='+NUMTOSTRING(data.a_over_Lne[ii,it])
  PRINTF, 1, 'DLNTDR_ELECTRON='+NUMTOSTRING(data.a_over_Lte[ii,it])
  PRINTF, 1, 'BETAE_UNIT='+NUMTOSTRING(data.betae_unit[ii,it])
  PRINTF, 1, 'NU_EI='+NUMTOSTRING(data.nu_ee[ii,it])

  PRINTF, 1, ' '
  OPENR, 2, GETENV('TGYRO_DIR') + '/tools/idl/INPUT.template'
  s = ' '
  WHILE (~EOF(2)) DO BEGIN
     READF, 2, s
     PRINTF, 1, s
  ENDWHILE
  CLOSE, 2
  CLOSE, 1

END ;write_tgyro2gyro_input
