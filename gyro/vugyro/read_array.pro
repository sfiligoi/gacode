;;------------------------------------------
;; Read array (array) from input file (file)
;; and mark existence in tag (exist_flag)
;;------------------------------------------
;;
;; CH fix 3.14.13- skip any lines at beginning starting with #
;;
pro read_array,array,file,exist_flag

  openr,1,file,error=i_err
  if (i_err eq 0) then begin

     s = ' '
     n_cmt = 0
     readf, 1, s
     while (strpos(s,'#') EQ 0) do begin
	n_cmt += 1
	readf, 1, s
     endwhile 
     close, 1
     openr, 1, file
     for ii = 0, n_cmt-1 do readf, 1, s

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

  
