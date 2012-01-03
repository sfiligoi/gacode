#!/usr/bin/ruby
h=Hash.new
filename=Array.new
r0=0.1875

[0,1,2,3].each do 
  |i|
  h.merge!("#{i}"=>i)
end

#puts h

["e", "i"].each do 
  |x|
  filename.push("Q-z-#{x}.out")
end


z_i=Array.new
q_i=Array.new
q_e=Array.new
z_e=Array.new
elem_q = Array.new
elem_z=Array.new
ct = 0
filename.each do 
  |x|
  file= File.new("data/#{x}","r")
  if x=~/i/ 
    ct=0
    while(a=file.gets)
      ct+=1
      elem_q.push(a.split(/\s+/)[-1])
      elem_z.push(a.split(/\s+/)[-2])
      if(ct==2)
        q_i.push(elem_q[0])
        q_e.push(elem_q[1])
        z_i.push(elem_z[0])
        ct=0
        elem_q = Array.new
        elem_z = Array.new
      end
    end
  end



        elem_q = Array.new
        elem_z = Array.new

  if x=~/e/ 
    ct=0
    while(a=file.gets)
      ct+=1
       elem_z.push(a.split(/\s+/)[-2])
      if(ct==2)
        z_e.push(elem_z[0])
        ct=0
        elem_z = Array.new
      end
    end
  end

  file_out=File.open("output/vec.out","w")
  file_out.puts z_i
  file_out.puts q_i
  file_out.puts z_e
  file_out.puts q_e

  file_out.close
  file.close


end








