# Usage:
#
#   python gacode_ufile_tool.py <datafile> <time>
#
#   NOTES: 
#   - Interpolate data to t=<time>.
#   - If time is absent, print list of tags.
#   - Operates on 1d and 2D datafiles
#   - DOC: http://tokamak-profiledb.ccfe.ac.uk/DOCS/PR08MAN/pdbman.html

import sys
import numpy as np

def extractcom(infile,t0):

    # Write the following strings:
    # <file rootname>
    # <tokamak>
    # <shot/pulse number>
    # <time>

    f=open('out.com','w')
    f.write(infile.split('_com')[0]+'\n')

    run   = ''
    pulse = ''

    for line in open(infile,'r').readlines():
        if line.count("Tokamak:") == 1:
            tok=line.split(":")[1].strip()
        if line.count("Run number:") == 1:
            run=line.split(":")[1].strip()
        if line.count("Pulse") == 1:
            pulse=line.split(":")[1].strip()

    f.write(tok+'\n')
    if len(run) == 0:
        f.write(pulse+'\n')
    else:
        f.write(run+'\n')

    # Also write the time 
    f.write(t0+'\n')
    f.close()


def extract1d(infile,t0):

    data_region = 0
    for line in open(infile,'r').readlines():

        if line.count("END-OF-DATA") == 1:
            data_region = 0
            if t0 >= vt[0] and t0 <= vt[-1]:
                # Compute average if time-point in range
                print 'Converting '+var+'  '+str(vt[0])+' <= '+str(t0)+ ' <= '+str(vt[-1])
                yave = np.zeros(1)
                yave[0] = np.interp(t0,vt,vy)
                # Write the averaged data for current profile (var)
                np.savetxt('out.'+var+'.ave',yave,fmt='%1.6e')
            else:
                print 'INFO: (ufile_tool) Time window: '+'t=['+str(vt[0])+','+str(vt[-1])+']'
                return

        if data_region == 0:

            if line.count("-DEP") == 1:
                # Extract current variable name
                var=line.split('  ')[0].strip()
             
            if line.count("OF PTS") == 1:
                # Get length of radial grid
                nt=int(line.split(";-# OF PTS")[0].strip())
                vt = np.zeros(nt)
                vy = np.zeros(nt)
                # Signal header end
                data_region = 1
                it = 0
                iy = 0

        else:

            # Here we are reading numbers

            # Number of columns 
            n = line.count("E")

            if it < nt:
                # Read the time grid [nt points]
                for i in range(n):
                    vt[it]=float(line[1+13*i:1+13*i+13])
                    it = it+1
            else:
                # Read the radial-time grid [(nx+1)*ny points]
                for i in range(n):
                    vy[iy]=float(line[1+13*i:1+13*i+13])
                    iy = iy+1

def extract2d(infile,t0):

    data_region = 0
    for line in open(infile,'r').readlines():
        #print line
        if line.count("END-OF-DATA") == 1:
            data_region = 0
            if t0 >= vt[0] and t0 <= vt[-1]:
                # Compute average if time-point in range
                fxt  = vy.reshape((nx,nt),order='F')
                yave = np.zeros(nx)
                print 'Converting '+var+'  '+str(vt[0])+' <= '+str(t0)+ ' <= '+str(vt[-1])
                for i in range(nx):
                    yave[i] = np.interp(t0,vt,fxt[i,:])
                    # Write the averaged data for current profile (var)
                    np.savetxt('out.'+var+'.ave',np.transpose((vx,yave)),fmt='%1.6e')
            else:
                print 'INFO: (ufile_tool) Time window: '+'t=['+str(vt[0])+','+str(vt[-1])+']'
                return

        if data_region == 0:

            if line.count("-DEP") == 1:
                # Extract current variable name
                var=line.split('  ')[0].strip()
             
            if line.count("X PTS") == 1:
                # Get length of radial grid
                nx=int(line.split(";-# OF X")[0].strip())
                vx = np.zeros(nx)

            if line.count("Y PTS") == 1:
                # Get size of time grid
                nt=int(line.split(";-# OF Y")[0].strip())
                vy = np.zeros(nx*nt)
                vt = np.zeros(nt)
                # Signal header end
                data_region = 1
                ix = 0
                it = 0
                iy = 0
                # Write dimensions
                if var=="NE":
                    f=open('out.dim','w')
                    f.write(str(nx)+'\n')
                    f.close()

        else:

            # Here we are reading numbers

            # Number of columns 
            n = line.count("E")

            if ix < nx:
                # Read the radial grid [nx points]
                for i in range(n):
                    vx[ix]=float(line[1+13*i:1+13*i+13])
                    ix = ix+1
            elif it < nt:
                # Read the time grid [nt points]
                for i in range(n):
                    vt[it]=float(line[1+13*i:1+13*i+13])
                    it = it+1
            else:
                # Read the radial-time grid [(nx+1)*ny points]
                for i in range(n):
                    vy[iy]=float(line[1+13*i:1+13*i+13])
                    iy = iy+1

#----------------------------------------------------------------------------
# Manage input parameters

# <datafile>
try:
    infile=sys.argv[1]
    infileval='0'
    if infile.split("_")[-1]=='1d.dat':
        infileval='1'
    if infile.split("_")[-1]=='2d.dat':
        infileval='2'
except:
    print 'Usage: python split.py <datafile> <time>'
    sys.exit()

if infile.split("_")[-1]=='com.dat':
    extractcom(infile,sys.argv[2])
    sys.exit()

# <time> Time-point for desired radial output
try:
    t0=float(sys.argv[2])
except:
    # GENERATE DIAGNOSTICS IF NO TIME SPECIFIED
    varlist = []
    t0 = -1.0
    # Generate list of included profiles (tags)
    for line in open(infile,'r').readlines():
        if line.count("-DEP") == 1:
            # Extract current variable name
            var=line.split('  ')[0].strip()
            varlist.append(var)
    # Print list of included profiles (tags) in blocks of 10
    print 'INFO: (ufile_tool) '+infileval+'D tags ->'
    for i in np.arange(start=0,stop=len(varlist),step=10):
        print '       '+' '.join(varlist[i:i+10])

if infile.split("_")[-1]=='1d.dat':
    extract1d(infile,t0)

if infile.split("_")[-1]=='2d.dat':
    extract2d(infile,t0)


