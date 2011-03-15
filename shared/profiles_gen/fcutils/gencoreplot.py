#!/usr/bin/env python
#
# $Id: gencoreplot.py 2598 2009-07-22 19:44:28Z cary $
import tables
import pylab
import optparse
import sys

# read in the options
cpsr = optparse.OptionParser()
cpsr.add_option('-i', '--input', dest='inputFile', action='store',
                help='Simulation name', default='')
cpsr.add_option('-d', '--dump', dest='dumpNumber', action='store',
		help='Dump number', default=0)
cpsr.add_option('-c', '--component', dest='component', action='store',
		help='Name of top-level component', default='facets')
cpsr.add_option('-v', '--variable', dest='plotVar', action='store',
                help='Variable to plot. One of density, temperature_electron or'
		' temperature_H2p1',
                default='temperature_H2p1')
cpsr.add_option('-y', '--ylimits', dest='ylimits', action='store',
		help='Limits for Y axis')

opts, args = cpsr.parse_args(sys.argv)

# open file to read data from
fh = tables.openFile('%s_%d.h5' % (opts.inputFile, int(opts.dumpNumber)), 'r')
# read appropriate data-set
dat = eval('fh.root.%s.%s' % (opts.component, opts.plotVar))

outName_mp_d = opts.inputFile + "_" + opts.plotVar
outName_mp = "%s_%05d.png" % (outName_mp_d, int(opts.dumpNumber))
outName_mp_dat = "%s_%05d.dat" % (outName_mp_d, int(opts.dumpNumber))

# get mesh
meshNm = eval('fh.root.%s.%s._v_attrs.vsMesh' % (opts.component, opts.plotVar))
# get mesh object
cells = eval('fh.root.%s.%s._v_attrs.vsNumCells' % (opts.component, meshNm))
lower = eval('fh.root.%s.%s._v_attrs.vsLowerBounds' % (opts.component, meshNm))
upper = eval('fh.root.%s.%s._v_attrs.vsUpperBounds' % (opts.component, meshNm))

dx = (upper[0]-lower[0])/cells[0]
rho = pylab.linspace(lower[0]+0.5*dx, upper[0]-0.5*dx, cells[0])

# plot data
fig1 = pylab.figure()
pylab.plot(rho, dat[:,0])
if opts.ylimits:
  ylimits = eval(opts.ylimits)
  fig1.gca().set_ylim( [ylimits[0], ylimits[1]] )


# put titles
pylab.xlabel('rho')
pylab.ylabel('%s' % opts.plotVar)
pylab.title('Time %g' % fh.root.timeData._v_attrs.vsTime)

# save the figure
fig1.savefig(outName_mp)
# now write the data
dataFp = open(outName_mp_dat, 'w')
# print some header information
dataFp.writelines('# Outboard midplane %s\n' % opts.plotVar)
# write data
for i in range(rho.shape[0]):
  dataFp.writelines('%g %g\n' % (rho[i], dat[i]))
dataFp.close()

