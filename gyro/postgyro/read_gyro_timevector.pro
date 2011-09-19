FUNCTION read_gyro_timevector, simdir, SILENT = silent
;
; read t.out and return time vector
; cut and past of read_time_vector.pro from vugyro
;
  filepath = simdir + '/out.gyro.t'

  OPENR, 1, filepath, ERROR=i_err
  IF (i_err eq 0) THEN BEGIN
    
    n_time = 0
    WHILE (NOT EOF(1)) DO BEGIN

      READF,1,dummy
      n_time = n_time+1

    END 
    CLOSE, 1
    
    n_time1 = n_time-1

    ;;-----------------------------------

    temp = fltarr(2,n_time)

    ;;-----------------------------------
    ;; Next, get time vector:
    ;;
    OPENR,1,filepath
    READF,1,temp
    CLOSE,1
    ;;

    t = temp[1,*]
    IF NOT KEYWORD_SET(silent) THEN PRINT, 'Read ' + filepath
  ENDIF ELSE BEGIN
      PRINT, "Couldn't load " + filepath + ' !!!'
      t = -1
  ENDELSE

  RETURN, t
END ;read_gyro_timevector

