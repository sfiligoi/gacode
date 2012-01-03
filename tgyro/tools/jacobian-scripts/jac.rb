#!/usr/bin/ruby
h=Hash.new
filename=Array.new
r0=0.1875

[0,1,2,3].each do 
  |i|
  h.merge!("#{i}"=>i)
end

puts h

["e", "i"].each do 
  |x|
  filename.push("jac-#{x}.out")
end


jac=Array.new
jac_i=Array.new
jac_e=Array.new
elem_i = Array.new
elem_e=Array.new
ct = 0
filename.each do 
  |x|
  file= File.new("data/#{x}","r")
  if x=~/i/ 
    ct=0
    while(a=file.gets)
      ct+=1
      elem_i.push(a.split(/\s/)[-1])
      if(ct==2)
        jac_i.push(elem_i)
        ct=0
        elem_i = Array.new
      end
    end
  end

  if x=~/e/ 
    ct=0
    while(a=file.gets)
      ct+=1
      elem_e.push(a.split(/\s/)[-1])
      if(ct==2)
        jac_e.push(elem_e)
        ct=0
        elem_e= Array.new
      end
    end
  end
#    end       

    0.upto(jac_i.length-1) do |i| 
      jac.push(jac_i[i], jac_e[i])
    end

    file_out=File.open("output/jac.out","w")
    file_out.puts jac
    file_out.close
#  end

  jac=Array.new
  file.close
end  








