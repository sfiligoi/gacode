pro helpme

  common GLOBAL

  ;;---------------------------------------------------------
  ;; Return conditions
  ;;
  if simdir eq 'null' then return
  ;;---------------------------------------------------------
 
  spawn,'cat run.out > vugyro.tmp'
  xdisplayfile,'vugyro.tmp'

end
