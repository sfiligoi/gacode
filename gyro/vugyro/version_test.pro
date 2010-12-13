;; Return:  flag = 0 if sim version is before version_min
;;               = 1 if sim_version is not before version_min

pro version_test,version_min,flag

  common GLOBAL

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
     return
  endif

  if (v1_sim gt v1_min) then begin
     flag = 1
     return
  endif

  ;; v1s are equal

  if (v2_sim lt v2_min) then begin
     flag = 0
     return
  endif

  if (v2_sim gt v2_min) then begin
     flag = 1
     return
  endif

  ;; v2s are equal

  if (v3_sim lt v3_min) then begin
     flag = 0
     return
  endif

  flag = 1

  return

end
