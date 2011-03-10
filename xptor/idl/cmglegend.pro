; CMGLEGEND                   C.M. Greenfield             June 23, 1999
;                             Last modification          March 10, 2000
;
; NAME:
;       CMGLEGEND
;
; PURPOSE:
;       Add a legend to a plot
;
; CATEGORY:
;       Data display
;
; CALLING SEQUENCE:
;       cmglegend,labels[,COLORS=COLORS][,LINE=LINE][,CSIZE=CSIZE] $
;               [,RIGHT=RIGHT][,BOTTOM=BOTTOM][,/NOLINES]
;               [,SYMBOL=SYMBOL][,SYMSIZE=SYMSIZE]
;
; INPUTS:
;       LABELS  Array of labels
;
; KEYWORD PARAMETERS:
;	COLORS  Array of color indices. Default is all drawn in color
;               !p.color
;       LINE    Array of line style indices. Default is all drawn in
;               style !p.linestyle
;       CSIZE   Character size. Default is 1.
;       RIGHT   Draw legend on right side of plot. Default is left
;               side.
;       BOTTOM  Draw legend at bottom of plot. Default is top.
;       NOLINES Only write text, skip drawing lines
;       SYMBOL  Array of plot symbols. If SYMBOL is a string array,
;               SYMSET will be called to define a user symbol,
;               otherwise, the symbol index will be used. 
;       SYMSIZE Symbol size. Default is 1.0
;       THICK   Line thickness
;
; MODIFICATIONS
;       9/10/99 CMG - Change name to CMGLEGEND to avoid conflict with other 
;                     routines (I should have known better!).
;       8/26/99 CMG - Allow symbols
;       8/25/99 CMG - Increased thickness of characters on X-windows device.
;       6/23/99 CMG - Changes to allow legends in log plots.
;
;-----------------------------------------------------------------------------
pro cmglegend,labels,colors=colors,line=line,sym=sym,csize=csize, $
           right=right,bottom=bottom,nolines=nolines,symbol=symbol, $
           symsize=symsize,thick=thick
psave=!p
n=n_elements(labels)
if keyword_set(line) then ls=line else ls=replicate(!p.linestyle,n)
if (not keyword_set(csize)) then csize=1.
;Determine size of characters
cs=convert_coord([!d.x_ch_size,0],[!d.y_ch_size,0],/device,/to_data)
if (!x.type eq 1) then cs[0,*]=alog10(cs[0,*])
if (!y.type eq 1) then cs[1,*]=alog10(cs[1,*])
if (not keyword_set(csize)) then csize=1.
cs=cs*csize
cx=cs[0,0]-cs[0,1]
;Set up line spacing
case !d.name of
    'X':cy=0.7*(cs[1,0]-cs[1,1])
    'PS':cy=cs[1,0]-cs[1,1]
endcase
dy=1.25*cy
ddy=cy/2
if keyword_set(bottom) then y=!y.crange[0]+(float(n)+0.5)*dy  $
else y=!y.crange[1]
if (keyword_set(nolines) and (n_elements(symbol) gt 0)) then  $
  shx=[0.01,0.06,0.06] $
  else shx=0
if keyword_set(right) then begin
    x=!x.crange[1]-([0.03,0.11,0.13]-shx)*(!x.crange[1]-!x.crange[0])
    align=1.0
endif else begin
    x=!x.crange[0]+([0.03,0.11,0.13]-shx)*(!x.crange[1]-!x.crange[0])
    align=0.0
endelse
if (!x.type eq 1) then x=10.^x
lines=not keyword_set(nolines)
symbols=n_elements(SYMBOL) gt 0
if symbols then begin
    string=(size(SYMBOL))[2] eq 7
    xsym=(x[0]+x[1])/2.
    if (n_elements(symsize) eq 0) then symsize=1
endif
for i=0,n-1 do begin
    y=y-dy
    yl=y
    yp=y+ddy
    if (!y.type eq 1) then begin
        yl=10.^yl
        yp=10.^yp
    endif
    if keyword_set(colors) then !p.color=colors(i)
    case !d.name of
        'PS':ct=1.0
        'X':ct=1.5
        ELSE:ct=1.0
    endcase
    if (not (lines or symbols)) then  $
      xyouts,x[0],yl,labels[i],charsize=csize,align=align,charthick=ct $
    else begin
        if (lines) then oplot,x[0:1],[yp,yp],line=ls[i],thick=thick
        if (symbols) then begin
            if (string) then begin
                symset,symbol[i],size=symsize
                oplot,[xsym],[yp],psym=8
                endif else oplot,[xsym],[yp],psym=symbol[i],symsize=symsize
        endif
        xyouts,x[2],yl,labels[i],charsize=csize,align=align,charthick=ct
    endelse
endfor
!p=psave

end
