;
; PLOTLABEL
;
; C. Holland, UCSD, 6/19/02
;
; simple routine for adding label to 2-D plot
; draws line from x0,y0 to x1,y1 (data coordinates),
; then puts label string at x1,y1
;
; Upadated 9/11/02 to use _EXTRA keyword
;

PRO plotlabel, x0, y0, x1, y1, label, _EXTRA = extra

  PLOTS, x0, y0, /DATA
  PLOTS, x1, y1, /DATA, /CONTINUE, _EXTRA = extra
  XYOUTS, x1, y1, label, _EXTRA = extra

END