;
; C. Holland, UCSD
; 8.25.2011- updated to be compatible with gacode file naming conventions
;
; simple script to extract electron density, temperature fluctuation
; fields and save them to IDL .sav files for synthetic diagnostic
; calculations.
;
;  to use, do
;  IDL> data = GET_GYRO_DATA('<simdir>', READ_LARGE=0)
;  IDL> MAKE_BIGFIELDS, data
;

PRO make_bigne, data, BIGDIR = bigdir
;
; bigdir is where moment_n.out is saved, defaults to data.simdir
;
  DEFAULT, bigdir, GETENV('GYRO_DIR') + '/sim/' + data.simdir
  PRINT, 'simdir: ', data.simdir
  PRINT, 'bigdir: ', bigdir

  OPENR, 1, bigdir + '/out.gyro.moment_n'
  ntp = data.n_theta_plot
  n_kin = data.n_kinetic
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

  SAVE, n_e, FILE=GETENV('GYRO_DIR') + '/sim/' + data.simdir + '/n_e.sav'
END ;make_bigne

PRO make_bigte, data, BIGDIR=bigdir
;
; bigdir is where moment_n.out is saved, defaults to data.simdir
; assumes n_e.sav exists in data.simdir
;
  DEFAULT, bigdir, GETENV('GYRO_DIR') + '/sim/' + data.simdir
  PRINT, 'simdir: ', data.simdir
  PRINT, 'bigdir: ', bigdir

  ;Te = (2/3)mom_e - mom_n
  ;first set Te = (2/3)mom_e from big file
  OPENR, 1, bigdir + '/out.gyro.moment_e'
  ntp = data.n_theta_plot
  n_kin = data.n_kinetic
  array = FLTARR(2,ntp,data.n_r,n_kin, data.n_n)
  Te = COMPLEXARR(ntp, data.n_r, data.n_n, data.n_time)
  help, array, Te

  ;Te = (2/3)mom_e - mom_n, 
  ;don't need any equilbirum n or T corrections because only e- fields
  ;first set Te = (2/3)mom_e from big file
  FOR it = 0, data.n_time-1 DO BEGIN
      print, it, data.n_time-1
      READF, 1, array
      Te[*,*,*,it] = (2./3)*COMPLEX(REFORM(array[0,*,*,n_kin-1,*]), $
                                    REFORM(array[1,*,*,n_kin-1,*]))
  ENDFOR
  CLOSE, 1
  RESTORE, GETENV('GYRO_DIR') + '/sim/' + data.simdir + '/n_e.sav'
  Te -= n_e

  SAVE, Te, FILE=GETENV('GYRO_DIR') + '/sim/' + data.simdir + '/Te.sav'
END ;make_bigte

PRO make_bigfields, data, BIGDIR=bigdir

  MAKE_BIGNE, data, BIGDIR=bigdir
  MAKE_BIGTE, data, BIGDIR=bigdir

END ;make_bigfields
