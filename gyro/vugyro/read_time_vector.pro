pro read_time_vector

  common GLOBAL

  file = 't.out'

  openr,1,file,error=i_err
  if (i_err eq 0) then begin
    
    n_time = 0
    while (not eof(1)) do begin

      readf,1,dummy
      n_time = n_time+1

    end 
    close,1
    
    n_time1 = n_time-1

    ;;-----------------------------------

    temp = fltarr(2,n_time)

    ;;-----------------------------------
    ;; Next, get time vector:
    ;;
    openr,1,file
    readf,1,temp
    close,1
    ;;
    print_found,file,1
    ;;-----------------------------------

    t = temp[1,*]

    t_min = max(t)/2.0
    t_max = max(t)

    ;; initialize title string

    get_t_string
    t_indices,t_min,t_max,it1,it2
    exists_time = 1

  endif else begin

    print_found,file,0
    exists_time = 0

    n_time = 0

  endelse

  return

end
