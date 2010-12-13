;
; sim_version_test
;
; vugyro version test, adapted to run in for a specific simulation
; specified by local directory name simdir.  version_min should be a
; float or string, not integer
;
;; Return:  flag = 0 if sim version is before version_min
;;               = 1 if sim_version is not before version_min

FUNCTION sim_version_test, simdir, version_min

  ;; read in version tag, default to local gyro version
  version_tag = ' '
  dirpath = GETENV('GYRO_DIR') + '/sim/' + simdir
  OPENR, lun, dirpath + '/VERSION_tag', ERR=ierr, /GET_LUN
  IF (ierr EQ 0) THEN BEGIN
      READF, lun, version_tag
      FREE_LUN, lun
  ENDIF ELSE BEGIN
      MESSAGE, 'VERSION_tag not found, defaulting to local GYRO version', $
               /INFO
      OPENR, lun, GETENV('GYRO_DIR') + '/VERSION', /GET_LUN
      READF, lun, version_tag
      FREE_LUN, lun
  ENDELSE

  ;; Strip simulation version

  a = version_tag

  d1 = strpos(a,'.')
  d2 = strpos(a,'.',d1+1)
  d3 = strpos(a,'-pre',d1+1)

  if (d3 eq -1) then begin
     d3 = strlen(a)
  endif

  v1_sim = float(strmid(a,0,d1))
  v2_sim = float(strmid(a,d1+1,d2-d1-1))
  v3_sim = float(strmid(a,d2+1,d3-d2-1))

  ;; Strip minimum version

  a = version_min

  d1 = strpos(a,'.')
  d2 = strpos(a,'.',d1+1)
  d3 = strpos(a,'-pre',d1+1)

  if (d3 eq -1) then begin
     d3 = strlen(a)
  endif

  v1_min = float(strmid(a,0,d1))
  v2_min = float(strmid(a,d1+1,d2-d1-1))
  v3_min = float(strmid(a,d2+1,d3-d2-1))

  if (v1_sim lt v1_min) then begin
     flag = 0
     return, flag
  endif

  if (v1_sim gt v1_min) then begin
     flag = 1
     return, flag
  endif

  ;; v1s are equal

  if (v2_sim lt v2_min) then begin
     flag = 0
     return, flag
  endif

  if (v2_sim gt v2_min) then begin
     flag = 1
     return, flag
  endif

  ;; v2s are equal

  if (v3_sim lt v3_min) then begin
     flag = 0
     return, flag
  endif

  flag = 1

  return, flag

end ;sim_version_test
