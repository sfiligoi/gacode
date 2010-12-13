pro midplane_deriv_plot

  common GLOBAL
  common MIDPLANE_DATA
  common PLOT_VARIABLES
  
  ;;=============================
  ;; radial averages
  ;;=============================
  
  yrt  = complexarr(n_r,n_time)
  f    = complexarr(n_r,n_time)
  fx   = complexarr(n_r,n_time)   
  yt   = fltarr(n_time)         
  y    = fltarr(n_n)         

  title = 'RMS '+tag[i_pwr]+': order '+strtrim(string(order),2)
  pname = 'dr'+strtrim(string(order),2)+'_ave_'+tag[i_pwr]

  y[*] = 0.0

  for i_n=0,n_n-1 do begin
     
     case (order) of

        0: begin

           yrt[*,*] = pwr[i_pwr,*,i_n,*]
           
        end

        1: begin

           f[*,*] = pwr[i_pwr,*,i_n,*]
           vector_deriv,r,f,yrt,boundary_method
           
        end

        2: begin

           f[*,*] = pwr[i_pwr,*,i_n,*]
           vector_deriv,r,f,fx,boundary_method
           vector_deriv,r,fx,yrt,boundary_method
           
        end

     endcase

     for i_time=it1,it2 do begin
        yt[i_time] = sqrt(total(abs(yrt[*,i_time])^2))/n_r     
     endfor
     diff_stat_fast,yt,it1,it2,ave_y

     y[i_n] = ave_y

  endfor

  plot_def_new,pname
  ytitle = '!3<f>!dRMS!n'

  ;; SCALE 
  x = kt_rho

  xtitle = kt_rho_string+' with'+t_string

  case (log_flag) of

     0: begin

        ;; LINEAR

        plot,[0],[0],$
          /nodata, $
          xtitle=xtitle, $
          xminor=1, $
          xstyle=1,xrange=[0,max(x)],$
          ytitle=ytitle, $
          ystyle=1,yrange=[0,max(y)],$
          title=title, $
          color=line

        oplot,x,y,psym=8,color=line
        oplot,x,y,color=color_vec[0]

     end

     1: begin

        ;; LOG

        plot_io,[0],[1],$
          /nodata, $
          xtitle=xtitle, $
          xminor=1, $
          xstyle=1,xrange=[0,max(x)],$
          ytitle=ytitle, $
          ystyle=1,yrange=[min(y),max(y)],$
          title=title, $
          color=line

        oplot,x,y,psym=8,color=line
        oplot,x,y,color=color_vec[0]

     end


  endcase

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for i=0,n_n-1 do begin
        printf,1,x[i],y[i]
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
