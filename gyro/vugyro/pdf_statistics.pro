pro pdf_statistics,y,i1,i2,y_bin,pdf

  common GLOBAL

  pdf = fltarr(nbin)
  y_bin = fltarr(nbin)

  ymin = min(y[i1:i2])
  ymax = max(y[i1:i2])

  d_y = (ymax-ymin)/(nbin-1)

  y_bin = ymin+findgen(nbin)*d_y

  for ibin=0,nbin-1 do begin
     pdf[ibin] = 0.0 
     for i=i1,i2-1 do begin

        dt = t[i+1]-t[i]
        f  = 0.5*(y[i+1]+y[i])
        
        f_min = y_bin[ibin]-0.5*d_y
        f_max = y_bin[ibin]+0.5*d_y 
        
        if (f ge f_min and f lt f_max) then begin
           pdf[ibin] = pdf[ibin]+dt
        endif
        
     endfor
  endfor

  pdf[*] = pdf[*]/total(pdf)

  return

end

   

   
