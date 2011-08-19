pro read_collision

  common GLOBAL
  common COLLISION_DATA

  file = 'out.gyro.collision_grid'  
  
  openr,1,file,error=i_err
  if (i_err eq 0) then begin

     readf,1,dummy
     readf,1,dummy
     readf,1,dummy
     readf,1,dummy
     readf,1,n_r_write
     i_r_write = intarr(n_r_write)
     readf,1,i_r_write

     yy = fltarr(n_stack,n_lambda,n_r_write)
     readf,1,yy
     xx = fltarr(n_stack,n_lambda,n_r_write)
     readf,1,xx

     print_found,file,1     
     exists_coll = 1

  endif else begin

     print_found,file,0
     exists_coll = 0

  endelse
  close,1

  return

end

