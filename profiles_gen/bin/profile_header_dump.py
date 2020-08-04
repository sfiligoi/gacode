import sys
import string

inputfile = sys.argv[1]

with open('profile_header','w') as outfile:

   index = 0
   with open(inputfile,'r') as fin:
      in_lines = fin.readlines()
   for line in in_lines:
      if index > 0:
         index = index+1
         line = line.strip()
         if len(line) > 1:
            x = line.split()
            outfile.write(x[2]+' '+x[3]+' '+x[4]+'\n')
         else:
            sys.exit()

      if 'IONS' in line:
         index = 1


