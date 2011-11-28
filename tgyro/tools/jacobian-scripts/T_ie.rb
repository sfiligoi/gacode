#!/usr/bin/ruby
file = File.new("data/profile.out", "r")
file1 = File.new("output/profiles.out","w")

while (b=file.gets)
  if b=~/\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ 

    r= $1
    ti = $4
    te = $5
      a=Array.new   

    [r, ti, te].each do 
      |x|
      if(x=~/^[a-z]/i)
        a.push("#"+x)
      else 
        a.push(x)
      end
    end
    file1.printf("%15s%15s%15s\n",a[0], a[1], a[2])

  end
end

file1.close
file.close
