import os
import numpy as np
import sys
import string

workdir = 'bdir'
tools   = os.environ['GACODE_ROOT']+'/neo/tools/'

if len(sys.argv) < 6:
   print "python neo_boot.py <rmin> <rmaj> <q> <nuee> <tau>"
   sys.exit()

rmin = sys.argv[1]
rmaj = sys.argv[2]
q    = sys.argv[3]
nuee = sys.argv[4]
tau  = sys.argv[5]

# Prepare simulation directory
os.system('rm -rf '+workdir)
os.system('mkdir '+workdir)
os.system('cp '+tools+'input.neo.neo_boot '+workdir+'/input.neo')

list = ['DLNNDR_1',
        'DLNNDR_2',
        'DLNNDR_3',
        'DLNTDR_1',
        'DLNTDR_2',
        'DLNTDR_3']

# Open input.neo, append parameters, close
# T_norm=T_e (species 1), m_norm=m_i (species 2)
neoin = open(workdir+'/input.neo','a') 

# Set input: r_min
neoin.write('RMIN_OVER_A='+rmin+'\n')

# Set input: r_maj
neoin.write('RMAJ_OVER_A='+rmaj+'\n')

# Set input: q
neoin.write('Q='+q+'\n')

# Set input: nu_ee/(cs/a)
neoin.write('NU_1='+nuee+'\n')

# Set input: tau=Ti/Te
neoin.write('TEMP_2='+tau+'\n')
neoin.write('TEMP_3='+tau+'\n')

neoin.close()

cneo = []
csauter = []

# Run NEO
for i in range(6):
   neoin = open(workdir+'/input.neo','a') 
   # Overlay gradients
   neoin.write('# '+str(i)+'\n')
   for j in range(6):
      if j == i:
         neoin.write(list[j]+'=1\n')
      else:        
         neoin.write(list[j]+'=0\n')
   neoin.close()
   os.system('neo -e '+workdir)

   # Harvest output
   neoout = np.loadtxt(workdir+'/out.neo.transport') 
   jneo=neoout[2]

   neoout = open(workdir+'/out.neo.diagnostic_geo','r') 
   line = neoout.readlines()[0]
   ipsi = float(string.splitfields(line,'=')[1].rstrip())

   neoout = np.loadtxt(workdir+'/out.neo.equil') 
   rhostar=neoout[3]

   neoout = np.loadtxt(workdir+'/out.neo.theory') 
   jsauter=neoout[10]

   print 'jneo        '+str(jneo)
   print 'jsauter     '+str(jsauter)
   print 'I*Psi*rho_* '+str(ipsi*rhostar)

   cneo.append(jneo/(ipsi*rhostar))
   csauter.append(jsauter/(ipsi*rhostar))

print list
print cneo
print csauter
