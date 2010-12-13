pro read_ballooning_mode

  common GLOBAL
  common PRIVATE_BALLOON

  n_balloon = 0

  balloon_tag = strarr(9)

  openr,1,'balloon_phi.out',error=i_err
  if (i_err eq 0) then begin
     balloon_tag[n_balloon] = 'balloon_phi'
     n_balloon = n_balloon+1
     close,1
  endif

  openr,1,'balloon_a.out',error=i_err
  if (i_err eq 0) then begin
     balloon_tag[n_balloon] = 'balloon_a'
     n_balloon = n_balloon+1
     close,1
  endif

  openr,1,'balloon_aperp.out',error=i_err
  if (i_err eq 0) then begin
     balloon_tag[n_balloon] = 'balloon_aperp'
     n_balloon = n_balloon+1
     close,1
  endif

  openr,1,'balloon_E_par.out',error=i_err
  if (i_err eq 0) then begin
     balloon_tag[n_balloon] = 'balloon_epar'
     n_balloon = n_balloon+1
     close,1
  endif

  openr,1,'balloon_n_ion.out',error=i_err
  if (i_err eq 0) then begin
     balloon_tag[n_balloon] = 'balloon_n_ion'
     n_balloon = n_balloon+1
     close,1
  endif
  openr,1,'balloon_n_elec.out',error=i_err
  if (i_err eq 0) then begin
     balloon_tag[n_balloon] = 'balloon_n_elec'
     n_balloon = n_balloon+1
     close,1
  endif

  openr,1,'balloon_e_ion.out',error=i_err
  if (i_err eq 0) then begin
     balloon_tag[n_balloon] = 'balloon_e_ion'
     n_balloon = n_balloon+1
     close,1
  endif
  openr,1,'balloon_e_elec.out',error=i_err
  if (i_err eq 0) then begin
     balloon_tag[n_balloon] = 'balloon_e_elec'
     n_balloon = n_balloon+1
     close,1
  endif

  openr,1,'balloon_v_ion.out',error=i_err
  if (i_err eq 0) then begin
     balloon_tag[n_balloon] = 'balloon_v_ion'
     n_balloon = n_balloon+1
     close,1
  endif
  openr,1,'balloon_v_elec.out',error=i_err
  if (i_err eq 0) then begin
     balloon_tag[n_balloon] = 'balloon_v_elec'
     n_balloon = n_balloon+1
     close,1
  endif

  ;; Time-independent bits of the above

  if n_balloon gt 0 then begin

     balloon_m    = fix(box_multiplier)
     balloon_np   = n_r/balloon_m
     balloon_next = n_theta_plot*balloon_np

     x1 = fltarr(2,balloon_next,balloon_m)
     balloon_plot = fltarr(2,balloon_next,balloon_m,n_balloon,n_time)

     for i=0,n_balloon-1 do begin
        openr,1,balloon_tag[i]+'.out'
        for tt=0,n_time1 do begin
           readf,1,x1
           balloon_plot[*,*,*,i,tt] = x1[*,*,*]
        endfor
        print_found,balloon_tag[i]+'.out',1
        close,1
     endfor

     exists_balloon = 1

  endif else begin

     print_found,'balloon*.out',0
     exists_balloon = 0

  endelse

  return

end

  
