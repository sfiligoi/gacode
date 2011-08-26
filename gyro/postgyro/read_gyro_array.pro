FUNCTION read_gyro_array, array, filename
;
; C. Holland, UCSD
;
; v1.0 1-30-2007
;
; Reads in a GYRO array file; straight copy of read_array.pro from
; vugyro.  Returns 1 on success, 0 on failure.
;
; Arguments:
; Array: data; needs to be allocated before read_gyro_array called.
; Filename: full pathname to file (i.e. not just u.out)
;
; v1.1 4-4-2007: changed to allocate lun rather than hard-coded value
; of 1

  array_exsists = 0

  OPENR, lun, filename, ERR = ierr, /GET_LUN
  IF (ierr EQ 0) THEN BEGIN
      READF, lun, array
      array_exsists = 1
      PRINT, 'Read ' + filename
      FREE_LUN, lun
 ENDIF ELSE BEGIN
      array_exsists = 0
      PRINT, "Couldn't read " + filename + ' !!!'
  ENDELSE
 
  RETURN, array_exsists
END ;read_gyro_array

