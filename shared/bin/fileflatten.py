# This is a utility file to take an ASCII data file with arbitrary 
# column lengths and create a new single-column file.
#
# The flattened file is given the aded extension '.conv'

import sys
import string

ifile = open(sys.argv[1],'r')
ofile = open(sys.argv[1]+'.conv','w')

for line in ifile.readlines():
    x = string.splitfields(line)
    for i in range(len(x)):
        ofile.write('%+.4E\n' % float(x[i]))

