pro bar_oplot,x,y,color

  n = n_elements(x)

  dx = (x[1]-x[0])*0.45

  for i=0,n-1 do begin
     oplot,(x[i]-dx)*[1,1],[0,y[i]],color=color
     oplot,(x[i]+dx)*[1,1],[0,y[i]],color=color
     oplot,[x[i]-dx,x[i]+dx],y[i]*[1,1],color=color
     oplot,[x[i]],[y[i]],psym=8
  endfor

end
