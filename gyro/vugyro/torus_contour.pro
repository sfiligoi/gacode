pro torus_contour,tt

  common GLOBAL
  common PLOT_VARIABLES
  common POLOIDAL_DATA

  x_cut = 0.86

  ny = mUt*n_theta_plot

  ;; Radii

  g = fltarr(ny,n_r)

  c = findgen(n_r)*r[0]/(n_r-1)

  g[*,*] = 0.0

  !p.region = [0,0,1.0,x_cut]

  ;;draw contour with messy inner circle

  loadct,CTab
  polar_contour,a,th,r, $
    nlevels=nlevels,$
    levels=clevels,$
    /fill, $
    xtitle='!3x/a', $
    ytitle='!3Z/a', $
    title=title+'   t = '+strtrim(string(t[tt]),2)

  ;;load fixed color table, and fill hole

  loadct,1

  if plot_mode eq 1 then begin
     hole_color=0
  endif else begin
     hole_color=255
  endelse

  polar_contour,g,th,c, $
    levels=[0],$
    nlevels=1,$
    /fill, $ 
    /overplot, $
    color=hole_color

  ;;reset colortable

  loadct,CTab

  !p.region = [0,x_cut,1.0,1.0]

  color_map,min(clevels),max(clevels),nlevels

  !p.region = 0

  plot_finish
  set_line_colors

end 
