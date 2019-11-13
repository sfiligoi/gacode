import sys
import string

inputfile = sys.argv[1]
outfile = open('profile_shot','w')

index = 0
for line in open(inputfile,'r').readlines():
   if 'SHOT' in line:
      index = 1

   if index > 0:
      line = string.strip(line)
      if len(line) > 1:
         x = line.split()
         outfile.write(' '+x[4]+'\n')
      else:
         break

index = 0
for line in open(inputfile,'r').readlines():
   if 'EFITD' in line:
      index = 1

   if index > 0:
      line = string.strip(line)
      if len(line) > 1:
         x = line.split()
         x[4] = x[4].strip('ms')
         outfile.write(' '+x[4]+'\n')
      else:
         sys.exit()

