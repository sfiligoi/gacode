pro view_efficiency

  common GLOBAL

  ;;---------------------------------------------------------
  ;; Return conditions
  ;;
  if simdir eq 'null' then return
  ;;---------------------------------------------------------
 
  spawn,'cat efficiency.out > vugyro.tmp'
  xdisplayfile,'vugyro.tmp'

end
