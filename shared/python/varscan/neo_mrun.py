import numpy as np
import os
import string


os.system('cp out.neo.localdump input.neo.template')

var1 = 'DLNTDR_1'
vec1 = [1.0,1.1]

var2 = 'DLNNDR_1'
vec2 = [1.0,1.1]

outfile = open('out.neo.array','w')
ovec = np.array([0,5,6,7,13,14,15])
for i1 in vec1:
    for i2 in vec2:
        a = var1+'='+str(i1)
        b = var2+'='+str(i2)
        os.system('cp input.neo.template input.neo')
        os.system('echo '+a+' >> input.neo')
        os.system('echo '+b+' >> input.neo')
        os.system('neo -e . > out')

        data = np.loadtxt('out.neo.transport')

        print a+' ; '+b
        for i in ovec:
            outfile.write(str(data[i])+' ')
        outfile.write('\n')

os.system('rm out')
