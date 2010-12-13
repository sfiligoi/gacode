pro time_phi_labels

  common GLOBAL
  common PLOT_VARIABLES

  ;;----------------------------------;
  ;; Some common plotting definitions:;
  
  dotsize = 1.5
  
  if (i_f eq 0) then begin
    ytitle  = ['!4u!3(r,!4h)!3 (Real and Imaginary)', $
               '!3<!4u!3(r)>!dRMS!n', $
               '!3<!4u(h)!3>!dRMS!n', $
               '!3<!4u!3>!dRMS!n']
  endif else begin
    ytitle  = ['!3A!d!9#!3!n(r,!4h)!3 (Real and Imaginary)', $
               '!3<A!d!9#!3!n(r)>!dRMS!n', $
               '!3<A!d!9#!3!n(!4h)!3>!dRMS!n', $
               '!3<A!d!9#!3!n>!dRMS!n']
  endelse
    
  title = ''
  
  ;;----------------------------------;

end
