pro color_map,f_min,f_max,nc

  common PLOT_VARIABLES

  ;; levels must be in ascending order

  if f_max < f_min then begin
    temp  = f_min
    f_min = f_max
    f_max = temp
  endif

  x_scale = fltarr(nc)
  for i=0,nc-1 do begin
    x_scale[i] = f_min+i*(f_max-f_min)/(nc-1.0)
  endfor

  f_scale = fltarr(nc,2)
  for i=0,nc-1 do begin
    f_scale[i,*] = x_scale[i]
  endfor

  contour,f_scale,x_scale,[1,2], $
      nlevels=nc,$
      /fill, $
      ystyle=1,$
      yticks=1,$
      ytickname=[' ',' '],$
      xstyle=1,$
      xticks=3,$
      xrange=[min(x_scale),max(x_scale)],$
      charsize=c_scale,$
      xmargin=[12,8],$
      /noerase

  oplot,[0,0],[1,2]

end
