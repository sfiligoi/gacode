pro nonlinear_transfer_event, nonlinear_transfer

  common GLOBAL
  common PRIVATE_NONLINEAR_TRANSFER

  widget_control, nonlinear_transfer.id, $
    get_uvalue=uvalue

  wset,  nonlinear_transfer_wid

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 

     1: begin
        goto, plot_it
     end   

     11: begin
        print, '************************************************************'
        print, '************************************************************'
        print, '* Hit NLTrEnGam 4 times to get first entropy "S" nonlinear *'
        print, '*  transfer Tr_k = -2*gamma_S/S_k followed by              *'
        print, '*  S_k,  rate 2*gamma_S, then ExB shear rate               *' 
        print, '*  kx-ky spectrum                                          *'
        print, '************************************************************'
        print, '************************************************************'

        goto, plot_it
     end

     2: begin
        zoom_TrEngGam = zoom_TrEngGam*2
        goto, plot_it
     end

     3: begin
        zoom_TrEngGam = zoom_TrEngGam/2
        goto, plot_it
     end


     4: begin
        TrEngGam_flag = 1+TrEngGam_flag
        if TrEngGam_flag eq 4 then begin
           TrEngGam_flag=0
        endif
        goto, plot_it
     end

     5: begin
        x2y_flag = 1+x2y_flag
        if x2y_flag eq 3 then begin
           x2y_flag=0
        endif
        goto, plot_it
     end



     6: begin
        plot_mode = 2
        goto, plot_it
     end

     8: widget_control, nonlinear_transfer.top, /destroy

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  ;; declare working variables for nonlinear transfer plots
  Tr_np_ave = fltarr(n_r,n_n)
  Eng_np_ave = fltarr(n_r,n_n)
  Tr_np_ave[*,*] = 0.0
  Eng_np_ave[*,*] = 0.0

  nlevels = 48

  z = fltarr(n_r,n_n)     
  x = fltarr(n_r)
  y = fltarr(n_n)

  xkr_rho = fltarr(n_r)
  for i=1,n_r-1 do begin
     xkr_rho[i]=kr_rho[i-1]
  endfor
  xkr_rho[0]=-kr_rho[n_r-1]


  ZZ = fltarr(n_r,n_n,n_time)

  if (TrEngGam_flag eq 0) then begin
     TrEngGamstr = '_entropy_nonlinear_transfer'
     ZZ[*,*,*] = nl_transfer[*,0,*,*]
  endif 

  if (TrEngGam_flag eq 1) then begin
     TrEngGamstr = '_turbulent_entropy'
     ZZ[*,*,*] = nl_transfer[*,1,*,*]
  endif 

  if (TrEngGam_flag eq 2) then begin
     TrEngGamstr = '_entropy_rate_2x_gamma_S'
     if (Tr_np_ave[n_r/2,n_n-1] ne 0.0) and (Tr_np_ave[n_r/2,n_n-1] ne 0.0) then begin
        for i_n = 1,n_n-1 do begin
           for i = 0,n_r-1 do begin
              z[i,i_n] = -Tr_np_ave[i,i_n]/Eng_np_ave[i,i_n]
           endfor
        endfor
        for i = 1,n_r-1 do begin
           z[i,0]= -Tr_np_ave[i,0]/Eng_np_ave[i,0]
        endfor
        for i = 0,-1 do begin
           z[i,0]= -Tr_np_ave[i,0]/Eng_np_ave[i,0]
        endfor
        z[n_r/2,0] = 0.0
     endif
  endif 

  if (TrEngGam_flag eq 3) then begin
     TrEngGamstr = '_ExB_shear_rate'

     print, 'x2y_flag=',x2y_flag

     if x2y_flag eq 0 then begin
        for i_n = 0,n_n-1 do begin
           for i=0,n_r-1 do begin
              ZZ[i,i_n,*]= (xkr_rho[i])^2*(kxkyspec[i,i_n,*])^0.5/rho_s
           endfor
        endfor
     endif
     if x2y_flag eq 1 then begin
        for i_n = 0,n_n-1 do begin
           for i=0,n_r-1 do begin
              ZZ[i,i_n,*]= (kt_rho[i_n]^2)*(kxkyspec[i,i_n,*])^0.5/rho_s
           endfor
        endfor
     endif
     if x2y_flag eq 2 then begin
        for i_n = 0,n_n-1 do begin
           for i=0,n_r-1 do begin
              ZZ[i,i_n,*]= (xkr_rho[i]^2+kt_rho[i_n]^2)*(kxkyspec[i,i_n,*])^0.5/rho_s
           endfor
        endfor
     endif

  endif


  title = TrEngGamstr+'_kx-ky_spectrum'
  pname = 'NLtr'+TrEngGamstr

  if (TrEngGam_flag ne 2) then begin
     z[*,*] = 0.0
     for i_n=0,n_n-1 do begin 
        for i_time=it1,it2 do begin       
           z[*,i_n] = z[*,i_n]+ZZ[*,i_n,i_time]
        endfor
        z[*,i_n] = z[*,i_n]/(it2-it1+1)
        ;;z[*,i_n] = shift(z[*,i_n],n_r/2-1)
        ;;  z[*,i_n] = shift(z[*,i_n],-1)
     endfor
  endif

  if (TrEngGam_flag eq 0) then begin
     Tr_np_ave[*,*] = z[*,*]
  endif
  if (TrEngGam_flag eq 1) then begin
     Eng_np_ave[*,*] = z[*,*]
  endif



  y[*] = kt_rho[*]
  x[*] = xkr_rho[*]
  

  plot_def_new,pname
  
  minz=min(z)
  maxz=max(z)

  print, 'minz=',minz,'maxz=',maxz,'nlevels=',nlevels

  print, 'kr_rho_max=',xkr_rho[n_r-1]
  print, 'kr_rho_min=',xkr_rho[0]

  ;; print data
  if (TrEngGam_flag eq 0) then begin

     Tr_tot = 0.0
     Tr_neg_tot = 0.0
     for i_n = 0, n_n-1 do begin
        for i=0,n_r-1 do begin
           Tr_tot = Tr_tot + z[i,i_n]
           if z[i,i_n] lt 0 then begin
              Tr_neg_tot = Tr_neg_tot + z[i,i_n]
           endif
        endfor
     endfor

     print, 'Tr_tot=',Tr_tot
     print, 'Tr_neg_tot=',Tr_neg_tot
     print, 'Max_Tr=',max(z) 
     print, 'Min_Tr=',min(z)

     openw,1,'entropy_nonlinear_transfer_data'

     printf,1, 'Tr_tot=',Tr_tot
     printf,1, 'Tr_neg_tot=',Tr_neg_tot
     printf,1, 'Max_Tr=',max(z)
     printf,1, 'Min_Tr=',min(z)


     printf,1, 0
     for i_n=0,n_n-1 do begin
        printf,1, kt_rho[i_n]
     endfor
     printf,1, 1
     for i_n=0,n_n-1 do begin
        printf,1, z[n_r/2,i_n]
     endfor
     printf,1, 2
     for i=0, n_r-1 do begin
        printf,1, xkr_rho[i]
     endfor
     printf,1, 3
     for i=n_r/2, n_r-1 do begin
        printf,1, z[i,0]
     endfor

     printf,1, '-------spectrum ---------------------------------------------'
     printf,1 , '--------------------------------------------'
     printf,1 ,  TrEngGamstr
     printf,1 , '--------------------------------------------'
     for i_n=0,n_n-1 do begin
        printf,1 , ' n=',i_n
        printf,1 , '--------------------------------------------'
        for i=0,n_r-1 do begin
           printf,1 , 'p=',i,' z=', z[i,i_n] , ' kr_rho=',xkr_rho[i]
        endfor
     endfor
     print , '--------------------------------------------'
     printf,1, '--------------------------------------------------------------'


     close,1

  endif


  if (TrEngGam_flag eq 1) then begin

     Eng_tot = 0.0
     for i_n = 0, n_n-1 do begin
        for i=0,n_r-1 do begin
           Eng_tot = Eng_tot + z[i,i_n]
        endfor
     endfor
     print, 'ave_S=',Eng_tot/(n_r*n_n)
     print, 'min_S=',min(z)
     print, 'max_S=',max(z)

     
     openw,2,'turbulent_entropy_data'

     printf,2, 'ave_S=',Eng_tot/(n_r*n_n)
     printf,2, 'min_S=',min(z)
     printf,2, 'max_S=',max(z)

     
     printf,2, 0
     for i_n=0,n_n-1 do begin
        printf,2, kt_rho[i_n]
     endfor
     printf,2, 1
     for i_n=0,n_n-1 do begin
        printf,2, z[n_r/2,i_n]
     endfor
     printf,2, 2
     for i=0, n_r-1 do begin
        printf,2, xkr_rho[i]
     endfor
     printf,2, 3
     for i=n_r/2, n_r-1 do begin
        printf,2, z[i,0]
     endfor

     printf,2, '-------spectrum ---------------------------------------------'
     printf,2 , '--------------------------------------------'
     printf,2 ,  TrEngGamstr
     printf,2 , '--------------------------------------------'
     for i_n=0,n_n-1 do begin
        printf,2 , ' n=',i_n
        printf,2 , '--------------------------------------------'
        for i=0,n_r-1 do begin
           printf,2 , 'p=',i,' z=', z[i,i_n] , ' kr_rho=',xkr_rho[i]
        endfor
     endfor
     print , '--------------------------------------------'
     printf,2, '--------------------------------------------------------------'


     close,2
     
  endif

  if (TrEngGam_flag eq 2) then begin

     print, 'max (2*gamma_S)=',max(z)
     print, 'min (2*gamma_S)=',min(z)

     openw,3,'gamma_data'

     printf,3, 'max (2*gamma_S)=',max(z)
     printf,3, 'min (2*gamma_S)=',min(z)

     printf,3, 0
     for i_n=0,n_n-1 do begin
        printf,3, kt_rho[i_n]
     endfor
     printf,3, 1
     for i_n=0,n_n-1 do begin
        printf,3, z[n_r/2,i_n]
     endfor
     printf,3, 2
     for i=0, n_r-1 do begin
        printf,3, xkr_rho[i]
     endfor
     printf,3, 3
     for i=n_r/2, n_r-1 do begin
        printf,3, z[i,0]
     endfor

     printf,3, '-------spectrum ---------------------------------------------'
     printf,3 , '--------------------------------------------'
     printf,3 ,  TrEngGamstr
     printf,3 , '--------------------------------------------'
     for i_n=0,n_n-1 do begin
        printf,3 , ' n=',i_n
        printf,3 , '--------------------------------------------'
        for i=0,n_r-1 do begin
           printf,3 , 'p=',i,' z=', z[i,i_n] , ' kr_rho=',xkr_rho[i]
        endfor
     endfor
     print , '--------------------------------------------'
     printf,3, '--------------------------------------------------------------'


     close,3

  endif

  if (TrEngGam_flag eq 3) then begin

     openw,4,'ExBshear_data'

     printf,4, 0
     for i_n=0,n_n-1 do begin
        printf,4, kt_rho[i_n]
     endfor
     printf,4, 1
     for i_n=0,n_n-1 do begin
        printf,4, z[n_r/2,i_n]
     endfor
     printf,4, 2
     for i=0, n_r-1 do begin
        printf,4, xkr_rho[i]
     endfor
     printf,4, 3
     for i=n_r/2, n_r-1 do begin
        printf,4, z[i,0]
     endfor

     printf,4, '-------spectrum ---------------------------------------------'
     printf,4 , '--------------------------------------------'
     printf,4 ,  TrEngGamstr
     printf,4 ,  'x2y_flag=',x2y_flag
     printf,4 , '--------------------------------------------'
     for i_n=0,n_n-1 do begin
        printf,4 , ' n=',i_n
        printf,4 , '--------------------------------------------'
        for i=0,n_r-1 do begin
           printf,4 , 'p=',i,' z=', z[i,i_n] , ' kr_rho=',xkr_rho[i]
        endfor
     endfor
     print , '--------------------------------------------'
     printf,4, '--------------------------------------------------------------'


     close,4

  endif




  loadct,CTab

  contour,z,x,y, $
    nlevels=nlevels,$
    /fill, $
    xtitle='!3'+kr_rho_string, $
    xstyle=1,$
    xrange=[min(x)/zoom_TrEngGam,-min(x)/zoom_TrEngGam],$
;;          xrange=[0.0,max(y)],$
  ytitle='!3'+kt_rho_string, $
    ystyle=1,$
;;          yrange=[min(y),max(y)],$
  yrange=[0.0,max(y)/zoom_TrEngGam],$
    title='!3'+title

  plot_finish

  set_line_colors

  return
  
end
