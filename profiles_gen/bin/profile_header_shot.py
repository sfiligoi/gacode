import sys
import string

inputfile = sys.argv[1]

shot = '0'
time = '0'
with open(inputfile,'r') as f:
    datalines = f.readlines()
for line in datalines:
   if 'SHOT' in line:
      line = line.strip()
      if len(line) > 1:
         x = line.split()
         shot = x[4]

   if 'EFITD' in line:
       line = line.strip()
       if len(line) > 1:
          x = line.split()
          y = x[4].split('ms')
          time = y[0]
          if shot == '0':
             shot = x[3].split('#')[1]

with open('profile_shot','w') as outfile:
   outfile.write(shot+'\n')
   outfile.write(time+'\n')
