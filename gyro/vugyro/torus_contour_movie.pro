pro torus_contour_movie,tt

  common GLOBAL
  common PLOT_VARIABLES
  common POLOIDAL_DATA

  x_cut = 0.86

  ny = mUt*n_theta_plot

  ;; Radii

  g = fltarr(ny,n_r)

  c = findgen(n_r)*r[0]/(n_r-1)

  g[*,*] = 0.0

  window,4,xsize=500,ysize=530

  r_0 = 1.05*max(r)

  loadct,CTab

  polar_contour,a,th,r, $
    nlevels=nlevels,$
    levels=clevels,$
    /fill,xstyle=1,ystyle=1,$
    xrange=[-r_0,r_0],yrange=[-r_0,r_0],$
    xmargin=[0,0],ymargin=[0,0]

  polar_contour,g,th,c, $
    levels=[0],$
    nlevels=1,$
    /fill, $ 
    /overplot, $
    color=(1-line)*255

  xyouts,0.3,-0.88,'t = '+strtrim(string(t[tt]),2),size=2,color=line

  plot_finish

  set_line_colors

  return

end 
