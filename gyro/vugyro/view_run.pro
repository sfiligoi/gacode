pro view_run

  common GLOBAL

  ;;---------------------------------------------------------
  ;; Return conditions
  ;;
  if simdir eq 'null' then return
  ;;---------------------------------------------------------
 
  spawn,'cat out.gyro.run > vugyro.tmp'
  xdisplayfile,'vugyro.tmp'

end
