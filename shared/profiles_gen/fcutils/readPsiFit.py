#!/usr/bin/env python
#
# $Id: readPsiFit.py 3057 2010-04-27 03:11:48Z ammar $
#

import sys, os
import numpy
import optparse
# Append to path to get other methods
fcutilsexecdir = os.path.dirname(sys.argv[0])
import fcaddpaths
fcaddpaths.fcaddpaths(fcutilsexecdir)
import txppImport
import tables

# read in the options
cpsr = optparse.OptionParser()
cpsr.add_option('-m', '--magnetic-geometry', dest='mapFile', action='store',
                help='Name of fluxgrid created magnetic geometry file', default=None)
cpsr.add_option('-o', '--output', dest='output', action='store',
		help='Name of output file', default=None)
cpsr.add_option('-r', '--rename', dest='rename', action='store',
		help='Name of output array', default=None)
cpsr.add_option('-c', '--conversion-factor', dest='conversion', action='store',
		help='Multiplication factor for output data', default='1.0')
cpsr.add_option('--h5output', dest='h5out', default='',
                help='Name of HDF5 file to write data')

opts, args = cpsr.parse_args(sys.argv)

killScript = False
if opts.mapFile == None:
    print "Must specify the name of magnetic geometry file created by fluxgrid"
    killScript = True
if opts.output == None:
    print "Must specify the name of output file to write"
    killScript = True
if opts.rename == None:
    print "Must specify the name of output array"
    killScript = True

if killScript:
    exit(1)

fileName = args[1]
# open file and read data
data = numpy.loadtxt(fileName, skiprows=2)

# read metadata on top of file
fieldName = open(fileName).readline()

psi = data[:,0]
psiVals = float(opts.conversion)*data[:,1]

# now open magnetic geometry file to read psi arrays
rhoList, rhoBool = txppImport.getVar(opts.mapFile, 'rho')
rho = numpy.array(rhoList).astype(numpy.double)
psiFgList, psiFgBool = txppImport.getVar(opts.mapFile, 'NormPolFlux')
psiFg = numpy.array(psiFgList).astype(numpy.double)

psiValsFg = numpy.interp(psiFg, psi, psiVals)

# now write data
outFile = open(opts.output, 'w')
outFile.writelines('rho = %s\n' % str(rho.tolist()))
outFile.writelines('%s = %s\n' % (opts.rename, str(psiValsFg.tolist())))

if opts.h5out != '':
  fh = tables.openFile(opts.h5out, "w")
  fh.createGroup("/", "profiles")
  fh.createArray("/profiles", "rho", rho)
  fh.createArray("/profiles", opts.rename, psiValsFg)

outFile.close()
