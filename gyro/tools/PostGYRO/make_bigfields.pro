;
; data is from get_gyro_data(XXX, read_large=5)
; 1-7-2009: generalized
;

PRO make_bigne, data, BIGDIR = bigdir
;
; bigdir is where moment_n.out is saved, defaults to data.simdir
;
  DEFAULT, bigdir, getenv('GYRO_DIR') +'/sim/'+ data.simdir  
  PRINT, 'simdir: ', data.simdir
  PRINT, 'bigdir: ', bigdir

  OPENR, 1, bigdir + '/moment_n.out'
  ntp = data.n_theta_plot
  n_kin = data.profile_data.n_kinetic
  array = FLTARR(2,ntp,data.n_r,n_kin, data.n_n)
  n_e = COMPLEXARR(ntp, data.n_r, data.n_n, data.n_time)
  help, array, n_e

  FOR it = 0, data.n_time-1 DO BEGIN  
      print, it, data.n_time-1
      READF, 1, array
      n_e[*,*,*,it] = COMPLEX(REFORM(array[0,*,*,n_kin-1,*]), $
                              REFORM(array[1,*,*,n_kin-1,*]))
  ENDFOR
  CLOSE, 1

  SAVE, n_e, file=getenv('GYRO_DIR')+'/sim/' + data.simdir + '/n_e.sav'
END ;make_bigne

PRO make_bigte, data, BIGDIR=bigdir
;
; bigdir is where moment_n.out is saved, defaults to data.simdir
; assumes n_e.sav exists in data.simdir
;
  DEFAULT, bigdir, getenv('GYRO_DIR')+'/sim/' + data.simdir  
  PRINT, 'simdir: ', data.simdir
  PRINT, 'bigdir: ', bigdir

  ;Te = (2/3)mom_e - mom_n
  ;first set Te = (2/3)mom_e from big file
  OPENR, 1, bigdir + '/moment_e.out'
  ntp = data.n_theta_plot
  n_kin = data.profile_data.n_kinetic
  array = FLTARR(2,ntp,data.n_r,n_kin, data.n_n)
  Te = COMPLEXARR(ntp, data.n_r, data.n_n, data.n_time)
  help, array, Te

  ;Te = (2/3)mom_e - mom_n
  ;first set Te = (2/3)mom_e from big file
  FOR it = 0, data.n_time-1 DO BEGIN
      print, it, data.n_time-1
      READF, 1, array
      Te[*,*,*,it] = (2./3)*COMPLEX(REFORM(array[0,*,*,n_kin-1,*]), $
                                    REFORM(array[1,*,*,n_kin-1,*]))
  ENDFOR
  CLOSE, 1
  RESTORE, getenv('GYRO_DIR')+'/sim/' + data.simdir + '/n_e.sav'
  Te -= n_e

  SAVE, Te, file=getenv('GYRO_DIR')+'/sim/' + data.simdir + '/Te.sav'
END ;make_bigte

PRO make_bigfields, data, BIGDIR=bigdir

  MAKE_BIGNE, data, BIGDIR=bigdir
  MAKE_BIGTE, data, BIGDIR=bigdir

END ;make_bigfields
