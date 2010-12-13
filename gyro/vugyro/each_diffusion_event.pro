pro each_diffusion_event, each_diffusion

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_EACH_DIFFUSION

  widget_control, each_diffusion.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: goto, plot_it

     10: begin
        print, '************************************************************'
        print, '************************************************************'
        print, '* cycle momnent "chi-i", "chi-e", "D-e" then               *'
        print, '* Hit QL 2x then Hit QL-mult 3x                            *'
        print, '* cycling moments each time                                *'
        print, '* i_QL_flag= 0 i_QL_mult=0  usual NL diffusivities         *'
        print, '* i_QL_flag= 1 i_QL_mult=0  QL-weight                      *'
        print, '* i_QL_flag= 2 i_QL_mult=0  norm-spec                      *'
        print, '* i_QL_flag= 1 i_QL_mult=1  QL diffusivities               *'
        print, '* i_QL_flag= 0 i_QL_mult=2  NL-QL weight                   *'
        print, '* i_QL_flag= 0 i_QL_mult=3  ratio  QL to NL-QL weights     *'
        print, '* i_QL-g-phi = 0 default  phi**2 norm                      *'
        print, '*   for LINDIFF_METHOD = 4(or5) input diff_QL_n.out        *'
        print, '************************************************************'
        print, '* Hit QL-g-phi to get                                      *'   
        print, '* i_QL-g-phi = 1  g**2 norm                                *'
        print, '*   for LINDIFF_METHOD = 6      input diff_QL_n.out        *'
        print, '* i_QL-g-phi = 2  phi(g)**2 norm (same as phi**2 norm)     *'
        print, '*   for LINDIFF_METHOD = 7      input diff_QL_n.out        *'
        print, '* i_QL-g-phi = 3  g*phi(g)  norm                           *'
        print, '*   for LINDIFF_METHOD = 8      input diff_QL_n.out        *'
        print, '************************************************************'
        print, '************************************************************'

        goto, plot_it
     end


     11: begin
        i_QL_flag= 1+i_QL_flag
        if (i_QL_flag eq 1) and (exists_diff_QL_n eq 0) then begin
           ;; i_QL_flag = 2
           i_QL_flag = 1
        endif
        if (i_QL_flag eq 2) and (exists_phi_squared_QL_n eq 0) then begin
           ;; i_QL_flag = 3
           i_QL_flag = 2
        endif
        if (i_QL_flag eq 3) then begin
           i_QL_flag = 0
        endif

        print, ' i_QL_flag=',i_QL_flag
        goto, plot_it
     end
     
     
     12: begin
        i_QL_mult = 1+i_QL_mult
        if (i_QL_mult eq 4) then begin
           i_QL_mult = 0
        endif
        
        print, ' i_QL_mult=',i_QL_mult
        goto, plot_it
     end

     13: begin
        i_QL_g_phi = 1+i_QL_g_phi

        if (i_QL_g_phi eq 4) then begin
           i_QL_g_phi = 0
        endif

        print, ' i_QL_g_phi=',i_QL_g_phi
        goto, plot_it
     end


     2: begin
        i_spec = i_spec+1
        if (i_spec eq n_kinetic) then i_spec = 0
        goto, plot_it
     end   

     3: begin
        i_moment = 1-i_moment
        goto, plot_it
     end   

     4: begin
        i_f = i_f+1
        if (i_f ge n_field) then i_f = -1
        goto, plot_it
     end

     5: begin
        frac_flag = 1-frac_flag
        goto, plot_it
     end

     6: begin
        bar_plot_flag = 1-bar_plot_flag
        goto, plot_it
     end

     7: begin
        plot_mode = 2
        goto, plot_it
     end

     8: begin
        plot_export = 1
        goto, plot_it
     end

     9: widget_control, each_diffusion.top, /destroy

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

 ;; Declare quasilinear variables
  diff_QLwt_n_ave      = fltarr(n_kinetic,n_field,2,n_n)
  diff_NLQLwt_n_ave    = fltarr(n_kinetic,n_field,2,n_n)
  phi_squared_QL_n_ave = fltarr(n_n)
  g_squared_QL_n_ave   = fltarr(3,n_n)

  print, '-------------------------------------------------------'
  print, ' i_QL_flag=',i_QL_flag,' i_QL_mult=',i_QL_mult 
  print, ' i_QL_g_phi=',i_QL_g_phi
  print, ' i_spec=',i_spec,' i_moment=',i_moment

  print, '-------------------------------------------------------'

  yn = fltarr(n_n)
  
  y  = fltarr(n_time)
  
  smf_tag,title,pname,ytitle

  if (i_QL_flag eq 1)  then begin
     title = title+' QL-'
  endif
  
  if (i_QL_flag eq 0) and (i_QL_mult eq 2) then begin
     title = title+' NLQL-'
  endif
  
  if (i_QL_flag eq 0) and (i_QL_mult eq 3) then begin
     title = title+' ratio QL/NLQL-'
  endif

  title = title+' Diffusion'

  if (i_QL_flag eq 1) and (i_QL_mult eq 0) then begin
     title = title+' weight'
  endif
  
  if (i_QL_flag eq 0) and (i_QL_mult eq 2) then begin
     title = title+' weight'
  endif

  if (i_QL_flag eq 0) and (i_QL_mult eq 3) then begin
     title = title+' weights'
  endif

  if (i_QL_flag eq 2) then begin
     if (i_QL_g_phi eq 0) then begin  
        title = 'phi_squared_QL'
     endif
     if (i_QL_g_phi eq 1) then begin
        title = 'g_squared_QL'
     endif
     if (i_QL_g_phi eq 2) then begin
        title = 'phi(g)_squared_QL'
     endif
     if (i_QL_g_phi eq 3) then begin
        title = 'g-phi_QL'
     endif
  endif
  

  if i_QL_flag eq 0 then begin
     pname = 'diff_n-'+pname
  endif
  if (i_QL_flag eq 1) and (i_QL_mult eq 0) then begin
     pname = 'diff_QL_wt_n-'+pname
  endif
  if (i_QL_flag eq 1) and (i_QL_mult eq 1) then begin
     pname = 'diff_QL_n-'+pname
  endif
  if (i_QL_flag eq 0) and (i_QL_mult eq 2) then begin
     pname = 'diff_NLQL_wt-'+pname
  endif
  if (i_QL_flag eq 0) and (i_QL_mult eq 3) then begin
     pname = 'diff_ratioQL_NLQL_wt-'+pname
  endif

  if i_QL_flag eq 2 then begin
     if (i_QL_g_phi eq 0) then begin
        pname = 'phi_squared_QL_n-'+pname
     endif
     if (i_QL_g_phi eq 1) then begin
        pname = 'g1_squared_QL_n-'+pname
     endif
     if (i_QL_g_phi eq 2) then begin
        pname = 'g2_squared_QL_n-'+pname
     endif
     if (i_QL_g_phi eq 3) then begin
        pname = 'g3_squared_QL_n-'+pname
     endif
  endif

  if (i_QL_flag eq 2) and (i_QL_mult eq 0)  then begin
     for i_n=0,n_n-1 do begin
        if (i_QL_g_phi eq 0) then begin
           y[*] = phi_squared_QL_n[i_n,*]
        endif
        if (i_QL_g_phi eq 1) then begin
           y[*] = g_squared_QL_n[0,i_n,*]
        endif
        if (i_QL_g_phi eq 2) then begin
           y[*] = g_squared_QL_n[1,i_n,*]
        endif
        if (i_QL_g_phi eq 3) then begin
           y[*] = g_squared_QL_n[2,i_n,*]
        endif

        diff_stat_fast,y,it1,it2,y_ave
        yn[i_n] = y_ave
        
        ;; test
        ;; print, 'i_n=',i_n
        ;; print, 'y[n_time-1]=',y[n_time-1], 'yn[i_n]=',yn[i_n]

        if (i_QL_g_phi eq 0) then begin
           phi_squared_QL_n_ave[i_n] = y_ave
 ;;;     print , ' n=',i_n, ' phi_squared_QL_n_ave=',phi_squared_QL_n_ave[i_n]
        endif

        if (i_QL_g_phi eq 1) then begin
           g_squared_QL_n_ave[0,i_n] = y_ave
 ;;;    print , ' n=',i_n, ' g1_squared_QL_n_ave=',g_squared_QL_n_ave[0,i_n]
        endif

        if (i_QL_g_phi eq 2) then begin
           g_squared_QL_n_ave[1,i_n] = y_ave
 ;;;    print , ' n=',i_n, ' g2_squared_QL_n_ave=',g_squared_QL_n_ave[1,i_n]
        endif

        if (i_QL_g_phi eq 3) then begin
           g_squared_QL_n_ave[2,i_n] = y_ave
 ;;;    print , ' n=',i_n, ' g3_squared_QL_n_ave=',g_squared_QL_n_ave[2,i_n]
        endif

     endfor
  endif

  if i_f ge 0 then begin

     if i_QL_flag eq 0 then begin
        for i_n=0,n_n-1 do begin
           y[*] = diff_n[i_spec,i_f,i_moment,i_n,*]
           diff_stat_fast,y,it1,it2,y_ave
           yn[i_n] = y_ave
        endfor 
     endif
     if (i_QL_flag eq 1) and (i_QL_mult eq 0)  then begin 
        for i_n=0,n_n-1 do begin
           y[*] = diff_QL_n[i_spec,i_f,i_moment,i_n,*]
           diff_stat_fast,y,it1,it2,y_ave
           yn[i_n] = y_ave
           diff_QLwt_n_ave[i_spec,i_f,i_moment,i_n] = yn[i_n]
        endfor
     endif

  endif else begin

     if i_QL_flag eq 0 then begin
        for i_n=0,n_n-1 do begin
           y[*] = diff_n[i_spec,0,i_moment,i_n,*]+ $
             diff_n[i_spec,1,i_moment,i_n,*]
           diff_stat_fast,y,it1,it2,y_ave
           yn[i_n] = y_ave
        endfor
     endif
     if (i_QL_flag eq 1) and (i_QL_mult eq 0)  then begin
        for i_n=0,n_n-1 do begin
           y[*] = diff_QLwt_n[i_spec,0,i_moment,i_n,*]+ $
             diff_QLwt_n[i_spec,1,i_moment,i_n,*]
           diff_stat_fast,y,it1,it2,y_ave
           yn[i_n] = y_ave
           diff_QLwt_n_ave[i_spec,i_f,i_moment,i_n] = yn[i_n]
        endfor
     endif

  endelse

  if (i_QL_mult eq 1) and (i_QL_flag eq 1) then begin
     print, ' i_QL_mult=',i_QL_mult,' i_QL_flag=',i_QL_flag
   ;;; print,  ' diff_QLwt_n_ave=',diff_QLwt_n_ave[i_spec,i_f,i_moment,i_n]

     if (i_QL_g_phi eq 0) then begin
        for i_n=0,n_n-1 do begin
           ;;  print, ' n=',i_n, ' phi_squared_QL_n_ave=',phi_squared_QL_n_ave[i_n]
           yn[i_n] = diff_QLwt_n_ave[i_spec,i_f,i_moment,i_n]*phi_squared_QL_n_ave[i_n]
   ;;;  print, ' n=',i_n,' diff_QLwt_n_ave X phi_squared_QL_n_ave=',yn[i_n]
        endfor
     endif

     if (i_QL_g_phi eq 1) then begin
        for i_n=0,n_n-1 do begin
   ;;;  print, ' n=',i_n, ' g_squared_QL_n_ave=',g_squared_QL_n_ave[0,i_n]
           yn[i_n] = diff_QLwt_n_ave[i_spec,i_f,i_moment,i_n]*g_squared_QL_n_ave[0,i_n]
   ;;;  print, ' n=',i_n,' diff_QLwt_n_ave X g_squared_QL_n_ave=',yn[i_n]
        endfor
     endif

     if (i_QL_g_phi eq 2) then begin
        for i_n=0,n_n-1 do begin
   ;;;  print, ' n=',i_n, ' g_squared_QL_n_ave=',g_squared_QL_n_ave[1,i_n]
           yn[i_n] = diff_QLwt_n_ave[i_spec,i_f,i_moment,i_n]*g_squared_QL_n_ave[1,i_n]
   ;;;  print, ' n=',i_n,' diff_QLwt_n_ave X g_squared_QL_n_ave=',yn[i_n]
        endfor
     endif

     if (i_QL_g_phi eq 3) then begin
        for i_n=0,n_n-1 do begin
   ;;;  print, ' n=',i_n, ' g_squared_QL_n_ave=',g_squared_QL_n_ave[2,i_n]
           yn[i_n] = diff_QLwt_n_ave[i_spec,i_f,i_moment,i_n]*g_squared_QL_n_ave[2,i_n]
   ;;;  print, ' n=',i_n,' diff_QLwt_n_ave X g_squared_QL_n_ave=',yn[i_n]
        endfor
     endif

  endif

  if (i_QL_mult eq 2) and (i_QL_flag eq 0) then begin
   ;;;  print, ' i_QL_mult=',i_QL_mult,' i_QL_flag=',i_QL_flag

     if (i_QL_g_phi eq 0) then begin
        for i_n=1,n_n-1 do begin
   ;;; print, ' n=',i_n, ' phi_squared_QL_n_ave=',phi_squared_QL_n_ave[i_n]
   ;;; print,  ' diff_n_ave=',yn[i_n]
           yn[i_n] = yn[i_n]/phi_squared_QL_n_ave[i_n]
           diff_NLQLwt_n_ave[i_spec,i_f,i_moment,i_n] = yn[i_n]
   ;;; print, ' n=',i_n,' diff_n_ave / phi_squared_QL_n_ave=',yn[i_n]
        endfor
        yn[0] = 0.
        diff_NLQLwt_n_ave[i_spec,i_f,i_moment,0] = 0.
     endif

     if (i_QL_g_phi eq 1) then begin
        for i_n=1,n_n-1 do begin
   ;;; print, ' n=',i_n, ' g_squared_QL_n_ave=',g_squared_QL_n_ave[0,i_n]
   ;;; print,  ' diff_n_ave=',yn[i_n]
           yn[i_n] = yn[i_n]/g_squared_QL_n_ave[0,i_n]
           diff_NLQLwt_n_ave[i_spec,i_f,i_moment,i_n] = yn[i_n]
   ;;; print, ' n=',i_n,' diff_n_ave /g_squared_QL_n_ave=',yn[i_n]
        endfor
        yn[0] = 0.
        diff_NLQLwt_n_ave[i_spec,i_f,i_moment,0] = 0.
     endif

     if (i_QL_g_phi eq 2) then begin
        for i_n=1,n_n-1 do begin
   ;;; print, ' n=',i_n, ' g_squared_QL_n_ave=',g_squared_QL_n_ave[1,i_n]
   ;;; print,  ' diff_n_ave=',yn[i_n]
           yn[i_n] = yn[i_n]/g_squared_QL_n_ave[1,i_n]
           diff_NLQLwt_n_ave[i_spec,i_f,i_moment,i_n] = yn[i_n]
   ;;; print, ' n=',i_n,' diff_n_ave /g_squared_QL_n_ave=',yn[i_n]
        endfor
        yn[0] = 0.
        diff_NLQLwt_n_ave[i_spec,i_f,i_moment,0] = 0.
     endif

     if (i_QL_g_phi eq 3) then begin
        for i_n=1,n_n-1 do begin
   ;;; print, ' n=',i_n, ' g_squared_QL_n_ave=',g_squared_QL_n_ave[2,i_n]
   ;;; print,  ' diff_n_ave=',yn[i_n]
           yn[i_n] = yn[i_n]/g_squared_QL_n_ave[2,i_n]
           diff_NLQLwt_n_ave[i_spec,i_f,i_moment,i_n] = yn[i_n]
   ;;; print, ' n=',i_n,' diff_n_ave /g_squared_QL_n_ave=',yn[i_n]
        endfor
        yn[0] = 0.
        diff_NLQLwt_n_ave[i_spec,i_f,i_moment,0] = 0.
     endif

  endif

  if (i_QL_flag eq 0) and (i_QL_mult eq 3) then begin
     for i_n=1,n_n-1 do begin
   ;;; print, ' n=',i_n, ' diff_QLwt_n_ave=',diff_QLwt_n_ave[i_spec,i_f,i_moment,i_n]
   ;;; print, ' n=',i_n, ' diff_NLQLwt_n_ave=',diff_NLQLwt_n_ave[i_spec,i_f,i_moment,i_n]
        yn[i_n] = diff_QLwt_n_ave[i_spec,i_f,i_moment,i_n]/diff_NLQLwt_n_ave[i_spec,i_f,i_moment,i_n]
   ;;; print, ' n=',i_n, ' diff_QLwt_n_ave/diff_NLQLwt_n_ave=',yn[i_n]
     endfor
     yn[0] = 0.
  endif


  yn_tot = total(yn)

  plot_def_new,pname

  x = kt_rho[*]
  xtitle = '!3'+kt_rho_string

  if frac_flag eq 1 then begin
     ytitle='!3fractional contribution per mode'
     yn[*] = yn[*]/yn_tot
  endif else begin
     ytitle='!3absolute contribution per mode'
  endelse

  xmin = min(x)
  xmax = max(x)+0.5*(x[1]-x[0])
  ymin = min(yn)
  ymax = max(yn)*1.3
  if min(yn) gt 0.0 then ymin = 0.0
  if min(yn) lt 0.0 then ymin = min(yn)*1.3

  if (i_QL_flag eq 0) and (i_QL_mult eq 3) then begin
     ymax = 10.0
  endif

  ;; temp
  ;; ymax = 7.0
  plot,[0],[0],$
    /nodata, $
    xtitle=xtitle+' with '+t_string, $
    xminor=1, $
    xstyle=1,xrange=[xmin,xmax],$
    ytitle=ytitle, $
    ystyle=1,yrange=[ymin,ymax],$
    title=title, $
    color=line

  oplot,[xmin,xmax],[0,0],linestyle=1

  if (bar_plot_flag eq 0) then begin
     bar_oplot,x,yn,color_vec[0]
  endif else begin
     oplot,x,yn,psym=8,color=line
     oplot,x,yn,color=color_vec[0]
  endelse

  ;;----------------------------------------------------
  ;; Labels for separate ITG-ETG diffusivity

  d_itg   = 0.0
  d_etg_1 = 0.0
  d_etg_2 = 0.0
  for i_n=0,n_n-1 do begin
     if (x[i_n] gt 2.0) then d_etg_2 = d_etg_2+yn[i_n]
     if (x[i_n] gt 1.0) then d_etg_1 = d_etg_1+yn[i_n]
     if (x[i_n] le 1.0) then d_itg = d_itg+yn[i_n]
  endfor

  x0 = xmin + 0.1*(xmax-xmin)

  str = strmid(strtrim(string(d_itg),2),0,5)
  y0 = ymax - 0.06*(ymax-ymin)
  xyouts,x0,y0,xtitle+' < 1 (ITG) ['+str+']'

  str = strmid(strtrim(string(d_etg_1),2),0,5)
  y0 = ymax - 0.11*(ymax-ymin)
  xyouts,x0,y0,xtitle+' > 1 (ETG) ['+str+']'

  str = strmid(strtrim(string(d_etg_2),2),0,5)
  y0 = ymax - 0.16*(ymax-ymin)
  xyouts,x0,y0,xtitle+' > 2 (ETG) ['+str+']'
  ;;----------------------------------------------------

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for i_n=0,n_n-1 do begin
        printf,1,x[i_n],yn[i_n]
     endfor
     printf,1,'Average:',t_min,t_max
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  endif
  ;;---------------------------------------------------

  plot_finish

  return

end
