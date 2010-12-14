pro get_dir_list

  common GLOBAL

  ;;-------------------------------------------------------------
  ;; Get directories  
  ;;
  print,simroot
  spawn,getenv("GACODE_ROOT")+'/gyro/bin/'+'vugyro_dir_init '+simroot
  openr,1,'~/.vugyrorc'

  dir = strarr(200)
  x = " "
  ndir = 0

  while not eof(1) do begin
    readf,1,x
    dir(ndir) = x
    ndir = ndir+1
  endwhile
  close,1
  ;;
  ;;-------------------------------------------------------------

end
