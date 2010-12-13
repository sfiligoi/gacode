
;========================================================================
 function envelope,input
;========================================================================

 n = n_elements(input)
 case n-(n/2)*2 of
    0:begin
      DATA = reverse(input(0:*))
      output = sqrt(DATA^2+float(hilbert(DATA))^2)
      output = reverse(output(0:*))
      DATA = reverse(DATA(0:*))
      output = output + sqrt(DATA^2+float(hilbert(DATA))^2)
      output = output / 2
      end
    1:begin
      DATA = reverse(input(1:*))
      output = sqrt(DATA^2+float(hilbert(DATA))^2)
      output = reverse(output(0:*))
      DATA = reverse(DATA(0:*))
      output = output + sqrt(DATA^2+float(hilbert(DATA))^2)
      hold = [0.0,output/2]
      DATA = reverse(input(0:n-2))
      output = sqrt(DATA^2+float(hilbert(DATA))^2)
      output = reverse(output(0:*))
      DATA = reverse(DATA(0:*))
      output = output + sqrt(DATA^2+float(hilbert(DATA))^2)
      output = [output/2,0.0]
      output = (hold + output) / 2
      output([0,n]) = output([0,n])*2
      end
 endcase

 return,output
 end
