 ;; Draw horizontal lines bracketing undamped region:

pro plot_bnd

  common GLOBAL
  common PLOT_VARIABLES

  oplot,r[n_bnd]*[1,1],1e3*[-1,1],color=line,linestyle=1
  oplot,r[n_r-1-n_bnd]*[1,1],1e3*[-1,1],color=line,linestyle=1

end
