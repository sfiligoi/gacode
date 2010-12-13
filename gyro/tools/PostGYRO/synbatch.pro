PRO synbatch, simdir, DIRPATH=dirpath, MAKE_BIG = make_big, $
              ITMIN=itmin, ITMAX=itmax, OMEGA0=omega0, $
              N_TOR_FRAC=n_tor_frac, _EXTRA=extra
;
; simple batch script for making synthetic BES/CECE/nT_Phase results.  Key
; keywords:
;
; ITMIN: index of time array to start at, defaults to n_time/2
; ITMAX: indx of time array to stop at, defaults to itmin
; OMEGA0: equilibrium rotation rate (in c_s/a units), from run.out,
; defaults to 0
; N_TOR_FRAC: # of toroidal locations to calculate at, defaults to 1
;

 
  data = GET_GYRO_DATA(simdir, READ_LARGE=5)

  IF KEYWORD_SET(make_big) THEN MAKE_BIGFIELDS, data

  dirpath = GETENV('GYRO_DIR') + '/sim/' + data.simdir + '/'
  RESTORE, dirpath + 'n_e.sav'
  RESTORE, dirpath + 'Te.sav'
  ;print,findfile(dirpath+'syncece.sav')
  IF (findfile(dirpath+'syncece.sav'))[0] eq '' then begin
    results = MAKE_SYNCECE_ARRAY(data, Te, ITMIN=itmin, ITMAX=itmax, $
                                OMEGA0=omega0, N_TOR_FRAC=n_tor_frac, $
                                _EXTRA = extra)
    SAVE, results, FILENAME = dirpath + 'syncece.sav'
  ENDIF
  if (findfile(dirpath+'psf.sav'))[0] eq '' then begin
    print,'You need to run auto_psfu, then save that structure to'
    print,dirpath+'psf.sav'
  endif else if (findfile(dirpath+'synbes.sav'))[0] eq '' then begin
    restore,dirpath+'psf.sav'
    results = make_synbes_array(data,n_e,psf,n_tor_frac=n_tor_frac,itmax=itmax,$
                                 itmin=itmin,omega0=omega0,_EXTRA=extra)
    save,results,filename = dirpath + 'synbes.sav'
  endif
  if (findfile(dirpath+'synrefl.sav'))[0] eq '' then begin
    results = make_synrefl_array(data, n_e, ITMIN = itmin, ITMAX = itmax, $
      OMEGA0 = omega0, N_tor_frac = n_tor_frac, _EXTRA=extra)
    SAVE, results, FILE = dirpath + 'synrefl.sav'
  endif
END ;synbatch
