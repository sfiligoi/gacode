FUNCTION numtostring, X, LENGTH = length
; takes number X, returns string of length characters
; updated 8/2/03 to include length

  IF N_ELEMENTS(x) EQ 0 THEN RETURN, '' ELSE BEGIN
	xstring = STRCOMPRESS(STRING(X), /REMOVE_ALL)

  	IF (N_ELEMENTS(length) EQ 0) THEN $
	xlen = STRLEN(xstring) ELSE xlen = length

  	RETURN, STRMID(xstring,0,xlen)
  ENDELSE
END
