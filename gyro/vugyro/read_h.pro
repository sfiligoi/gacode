pro read_h

  common GLOBAL

  file='theta_m.out'

  openr,1,file,error=i_err
  if (i_err eq 0) then begin

     if (n_pass gt 0) then begin
        theta_t_p = fltarr(n_theta_p,n_pass)
        dum_p = fltarr(n_theta_p)
        for k=1,n_pass do begin
           readf,1,dum_p
           theta_t_p[*,k-1] = dum_p[*]
        endfor
     endif

     if (n_trap gt 0) then begin
        theta_t_t = fltarr(n_theta_t,n_trap)    
        dum_t = fltarr(n_theta_t)
        for k=1,n_trap do begin
           readf,1,dum_t
           theta_t_t[*,k-1] = dum_t[*]
        endfor
     endif

     if (n_pass gt 0) then begin
        xi_p = fltarr(n_theta_p,n_pass)
        dum_p = fltarr(n_theta_p)
        for k=1,n_pass do begin
           readf,1,dum_p
           xi_p[*,k-1] = dum_p[*]
        endfor
     endif

     if (n_trap gt 0) then begin
        xi_t = fltarr(n_theta_t,n_trap)    
        dum_t = fltarr(n_theta_t)
        for k=1,n_trap do begin
           readf,1,dum_t
           xi_t[*,k-1] = dum_t[*]
        endfor
     endif

     print_found,file,1

  endif else begin

     print_found,file,0
     
  endelse

  close,1

  file = 'hp.out'

  if (n_pass gt 0) then begin
     openr,1,file,error=i_err
     if (i_err eq 0) then begin

        hp = fltarr(2,n_theta_p,n_pass,n_kinetic,n_time)
        readf,1,hp

        print_found,file,1
        exists_h = 1

     endif else begin

        print_found,file,0
        exists_h = 0

     endelse
     close,1
  endif

  file = 'ht.out'

  if (n_trap gt 0) then begin
     openr,1,file,error=i_err
     if (i_err eq 0) then begin
        
        ht = fltarr(2,n_theta_t,n_trap,n_kinetic,n_time)
        readf,1,ht

        print_found,file,1
        exists_h = 1

     endif else begin

        print_found,file,0
        exists_h = 0

     endelse
     close,1
  endif

end
