;;------------------------------------------
;; Read array (array) from input file (file)
;; and mark existence in tag (exist_flag)
;;------------------------------------------

pro read_array,array,file,exist_flag

  openr,1,file,error=i_err
  if (i_err eq 0) then begin

     readf,1,array
     print_found,file,1
     exist_flag = 1

  endif else begin
     
     array = fltarr(1)
     print_found,file,0
     exist_flag = 0

  endelse

  close,1

  return

end

  
