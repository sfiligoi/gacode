pro gread

  id = ncdf_open("a.nc")

  info = ncdf_inquire(id)

  print,info

  ncdf_varget,id,"r",r
  ncdf_varget,id,"q",q
  ncdf_varget,id,"theta",theta
  ncdf_varget,id,"t",t
  ncdf_varget,id,"n",n
  ncdf_varget,id,"diff_t",diff_t
  ncdf_varget,id,"potential-r",u
  ncdf_varget,id,"potential-theta",u

  i_field = 0
  i_time = 300

  xsize=350
  ysize=350

  print,n

  i_n = 0
  window,1,xsize=xsize,ysize=ysize
  plot,r,u[0,*,i_n,i_time],xstyle=1,xrange=[min(r),max(r)]
  oplot,r,u[1,*,i_n,i_time],linestyle=1

  i_n = 1
  window,2,xsize=xsize,ysize=ysize
  plot,r,u[0,*,i_n,i_time],xstyle=1,xrange=[min(r),max(r)]
  oplot,r,u[1,*,i_n,i_time],linestyle=1

;  i_n = 2
;  window,3,xsize=xsize,ysize=ysize
;  plot,r,u[*,i_field,i_n,i_time,0],xstyle=1,xrange=[min(r),max(r)]
;  oplot,r,u[*,i_field,i_n,i_time,1],linestyle=1
 
  window,3,xsize=xsize,ysize=ysize
  plot,r,diff_t[*,i_time],xstyle=1,xrange=[min(r),max(r)]

  n_time=n_elements(t)
  ave = fltarr(n_time)
  for i_time=0,n_time-1 do begin
  ave[i_time] = total(diff_t[*,i_time])/n_elements(r)
  endfor

  window,4,xsize=xsize,ysize=ysize
  plot,t,ave[*],xstyle=1,xrange=[min(t),max(t)]

end
