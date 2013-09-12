import numpy as np
import os
import string
import sys


# ---- BEGIN input parsing ----

x = sys.argv[1]
var1 = x.split(',')[0]
var2 = x.split(',')[1] 

x = sys.argv[2]
min1 = float(x.split(',')[0])
min2 = float(x.split(',')[1])

x = sys.argv[3]
max1 = float(x.split(',')[0])
max2 = float(x.split(',')[1])

x = sys.argv[4]
n1 = int(x.split(',')[0])
n2 = int(x.split(',')[1])

run = int(sys.argv[5])

dir = sys.argv[6]

# ---- END input parsing ----

vec1 = min1+(max1-min1)*np.arange(n1)/float(n1-1)
vec2 = min2+(max2-min2)*np.arange(n2)/float(n2-1)

if run == 0:
    print var1,vec1
    print var2,vec2
    sys.exit()

# IMPORTANT: Set working directory
os.chdir(dir)

os.system('cp input.neo input.neo.template')

outfile = open('out.neo.mrun_z','w')

for i2 in vec2:
    for i1 in vec1:
        a = var1+'='+str(i1)
        b = var2+'='+str(i2)
        os.system('cp input.neo.template input.neo')
        os.system('echo '+a+' >> input.neo')
        os.system('echo '+b+' >> input.neo')
        os.system('neo -p '+dir+' -e . > out')

        data = np.loadtxt('out.neo.transport')
        print 'INFO: (neo_mrun) '+a+' ; '+b
        for i in range(len(data)):
            outfile.write(str(data[i])+' ')
        outfile.write('\n')

outfile.close()

# X:
np.savetxt('out.neo.mrun_x',vec1)
# Y:
np.savetxt('out.neo.mrun_y',vec2)

os.system('rm out')
os.system('mv input.neo.template input.neo')

