pro midplane_fluc_plot

  common GLOBAL
  common MIDPLANE_DATA
  common PLOT_VARIABLES
  common POLOIDAL_DATA

  if (c_table_max > 0.0) then begin

     clevels = c_table_min+$
       findgen(nlevels)*(c_table_max-c_table_min)/(nlevels-1.0)

  endif

  cn = fltarr(n_n)
  
  cn[*] = 2.0
  cn[0] = 1.0

  title = tag[i_pwr]+$
    ' fluctuations (t='+$
    strtrim(string(t[t_c]),2)+')'

  pname = 'pwr_n_'+tag[i_pwr]

  n_y = mUt*n_n

  y = fltarr(2*n_r-1,n_y)
  y[*,*] = 0.0

  alpha = findgen(n_y)/(n_y-1)

  evec = complexarr(n_n,n_y)

  for i_n=0,n_n-1 do begin
     evec[i_n,*] = exp(-2*!pi*complex(0,1)*i_n*alpha[*])*cn[i_n]
  endfor
  
  z = complexarr(n_n,2*n_r-1)
  z[*,0] = pwr[i_pwr,0,*,t_c]
  re = r[0]+findgen(2*n_r-1)*(r[n_r-1]-r[0])/(2*n_r-2)

  for ii=1,n_r-1 do begin
     z[*,2*ii] = pwr[i_pwr,ii,*,t_c]
     z[*,2*ii-1] = 0.5*(pwr[i_pwr,ii,*,t_c]+pwr[i_pwr,ii-1,*,t_c])
  endfor  

  for i=0,2*n_r-2 do begin
     for i_y=0,n_y-1 do begin
        y[i,i_y] = y[i,i_y]+$ 
          total(float(z[*,i]*evec[*,i_y]))
     endfor
  endfor

  plot_def_new,pname
  
  xtitle = '!3r/a'

  loadct,CTab
  contour,y,re,alpha, $
    nlevels=nlevels,$
    levels=clevels,$
    xtitle=xtitle, $
    xstyle=1,$
    xrange=[min(r),max(r)],$
    ytitle='!4a!3', $
    ystyle=1,$
    yrange=[0,1],$
    title='!3'+title,$
;    charsize=c_size*c_scale,$
    /fill

  set_line_colors
  plot_finish

  return

end
