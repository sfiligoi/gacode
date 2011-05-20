#!/usr/bin/env python
"""
  Simple example showing how to use readInputProfiles.py and
  readPfile.py
"""
import re, sys, os, glob, filecmp, datetime, shutil
from string_to_number import *
import numpy
import readPfile
import readInputProfiles
import scipy
from matplotlib import pylab


# Read in the data.  This is where I hardcoded values.
pf=readPfile.pfileData("p132014_3000_f8099.txt")
idb=readInputProfiles.inputProfileData("input.profiles")

# Interpolate pfile data onto the input.profiles grid and translate units
# The differences in units by formats developed by guys across the 
#   hall from each other is something you have to just smile at.
#   (and we have the same problems at Tech-X)
xidb=idb.polflux.values
xmin=xidb[0];       xmax=xidb[xidb.shape[0]-1]
xidb=(xidb-xmin)/(xmax-xmin)
teinterp=scipy.interp(xidb, pf.te.grid, pf.te.values)
tiinterp=scipy.interp(xidb, pf.ti.grid, pf.ti.values)
neinterp=scipy.interp(xidb, pf.ne.grid, pf.ne.values) * 10.
niinterp=scipy.interp(xidb, pf.ni.grid, pf.ni.values) * 10.
ominterp=scipy.interp(xidb, pf.omeg.grid, pf.omeg.values)

# Assign interpolated values 
idb.ne.values=neinterp
idb.ni_1.values=niinterp
idb.Te.values=teinterp
idb.Ti_1.values=tiinterp

# Now write out the values
idb.writeInputFile("input.profiles_wPfileMod")

# Read in the new values 
idm=readInputProfiles.inputProfileData("input.profiles_wPfileMod")
xidm=idm.polflux.values
xmin=xidm[0];       xmax=xidm[xidm.shape[0]-1]
xidm=(xidm-xmin)/(xmax-xmin)
tedm=idm.Te.values
tidm=idm.Ti_1.values
nedm=idm.ne.values
nidm=idm.ni_1.values
# Compare against pfile and old values
xtepf=pf.te.grid;     tepf=pf.te.values
xtipf=pf.ti.grid;     tipf=pf.ti.values
xnepf=pf.ne.grid;     nepf=pf.ne.values * 10.
xnipf=pf.ni.grid;     nipf=pf.ni.values * 10.
idb=readInputProfiles.inputProfileData("input.profiles")
tedb=idb.Te.values
tidb=idb.Ti_1.values
nedb=idb.ne.values
nidb=idb.ni_1.values

#sys.exit()

# Read in the new values and compare via plots
pylab.subplot(211)
pylab.plot(xidb,tedb,marker="1",markevery=4,label="Te: old",ms=20.0)
pylab.plot(xidb,tidb,marker="2",markevery=4,label="Ti: old",ms=20.0)
pylab.plot(xidm,tedm,marker="o",markevery=4,label="Te: input.profiles")
pylab.plot(xtepf,tepf,marker="*",markevery=4,label="Te: pfile")
pylab.plot(xidm,tidm,marker="x",markevery=4,label="Ti: input.profiles")
pylab.plot(xtipf,tipf,marker="+",markevery=4,label="Ti: pfile")
pylab.legend()
pylab.subplot(212)
pylab.plot(xidb,nedb,marker="1",markevery=4,label="ne: old",ms=20.0)
pylab.plot(xidb,nidb,marker="2",markevery=4,label="ni: old",ms=20.0)
pylab.plot(xidm,nedm,marker="o",markevery=4,label="ne: input.profiles")
pylab.plot(xnepf,nepf,marker="*",markevery=4,label="ne: pfile")
pylab.plot(xidm,nidm,marker="x",markevery=4,label="ni: input.profiles")
pylab.plot(xnipf,nipf,marker="+",markevery=4,label="ni: pfile")
pylab.legend()
pylab.show()
