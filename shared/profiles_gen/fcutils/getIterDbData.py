#!/usr/bin/env python
"""
  Use the readIterDB module to get specific names
  Typical usage:
    getIterDbData.py --print iterdb_118898_03400_e0120 | more
    <select field names>
    getIterDbData.py -f qbeame,qbeami -o iterdb118898.03400.pre iterdb_118898_03400
"""

import sys, os
import re
import numpy
import optparse
import tables

# Append to path to get other methods
fcutilsexecdir = os.path.dirname(sys.argv[0])
import fcaddpaths
fcaddpaths.fcaddpaths(fcutilsexecdir)

print "sys.path =", sys.path

import readIterDB
import txppImport

def writeNoInterp(myidb,fp,allFields,options):
    """
    Just write the fields as given in the iterDB file to the filehandle fp
    """
    
    # Check if we are to write an h5 output file
    h5fp = None
    if options.h5out != '':
      # open file and create a group inside it
      h5fp = tables.openFile(options.h5out, "w")
      h5fp.createGroup("/", "profiles")

    ifield=-1
    allFields=["rho"]+allFields
    for field in allFields:
      if not field=="rho": ifield=ifield+1
      if myidb.variables.has_key(field):
        if not myidb.variables[field].has_key('values'):
          print "Field "+field+" does not have values from parsing"
          continue
        values=myidb.variables[field]['values']
        if len(options.convertFields)>0 and not field=="rho":
          values=values*float(options.convertFields[ifield])
        if len(options.renameFields)>0 and not field=="rho":
          fieldName=options.renameFields[ifield]
        else:
          fieldName=field
        fp.writelines(fieldName+" = %s\n" % str(values.tolist()))
	if h5fp:
	  h5fp.createArray("/profiles", fieldName, values)
      else:
        print "Field not found in iterdb file: ", field
    return

def writeWithInterp(myidb,fp,allFields,mapFile,options):
    """
    Interpolate and then write
    """

    # Check if we are to write an h5 output file
    h5fp = None
    if options.h5out != '':
      # open file and create a group inside it
      h5fp = tables.openFile(options.h5out, "w")
      h5fp.createGroup("/", "profiles")

    # Get the grid that we will want to interpolate too
    rhoList, FAIL = txppImport.getVar(mapFile,"rho")
    if FAIL:
      print "ERROR: Unable to get rho variable from pre file: ", mapfile
      return FAIL
    rhoMap=numpy.array(rhoList).astype(numpy.double)

    # Get the rho the iterDB fields are on
    if myidb.variables.has_key("rho"):
      rho=myidb.variables['rho']['values']
    else:
      print "ERROR: Unable to get rho from iterDB file "
      return 1

    # Interpolate all of the fields specified to the new rho
    mappedFields={}
    mappedFields['rho']=rhoMap
    for field in allFields:
      if myidb.variables.has_key(field):
        if not myidb.variables[field].has_key('values'):
          print "Field "+field+" does not have values from parsing"
          continue
        values=myidb.variables[field]['values']
        mappedFields[field]=numpy.interp(rhoMap,rho,values)
      else:
        print "Field not found in iterdb file: ", field

    # Write out to the pre file
    ifield=-1
    for field in ['rho']+allFields:
      if not field=="rho": ifield=ifield+1
      values=mappedFields[field]
      if len(options.convertFields)>0 and not field=="rho":
        values=values*float(options.convertFields[ifield])
      if len(options.renameFields)>0 and not field=="rho":
        fieldName=options.renameFields[ifield]
      else:
        fieldName=field
      fp.writelines(fieldName+" = %s\n" % str(values.tolist()))
      if h5fp:
	h5fp.createArray("/profiles", fieldName, values)
    return


def main():

    parser = optparse.OptionParser(usage="%prog [options] iterDBfile")
    parser.add_option('-p', '--print', dest='doPrint',
                      help='Print all possible fields to standard out along with description',
                      action='store_true')
    parser.add_option('-i', '--input', dest='input',
                      help='Name of ITER DB input file if not specified as argument.',
                      default='')
    parser.add_option('-o', '--output', dest='outputFile',
                      help='Name of output file (.pre format).',
                      default='')
    parser.add_option('-f', '--field', dest='fields',
                      help='Name of field(s) to get (rho always included).\
                      For multiple fields, use comma-delimited list. \
                      Default="qbeame,qbeami"',
                      default='qbeame,qbeami')
    parser.add_option('-r', '--rename', dest='renameFields',
                      help='Name of field(s) to appear in prefile.\
                      For multiple fields, use comma-delimited list. \
                      The names must be in the same order as the fields"',
                      default='')
    parser.add_option('-c', '--convert', dest='convertFields',
                      help='Conversion factors between the iterDB fields \
                      and pre-file fields: \
                        preField=(Conversion Factor)*Field           \
                      For multiple fields, use comma-delimited list. \
                      The names must be in the same order as the fields"',
                      default='')
    parser.add_option('-m', '--map', dest='mapFile',
                      help='Name of pre file containing a rho variable to interpolate the fields onto.',
                      default='')
    parser.add_option('--h5output', dest='h5out', default='',
                      help='Name of HDF5 file to write data')

    options, args = parser.parse_args()

    # Too many arguments
    if len(args) > 1:
      parser.print_usage()
      return
    elif len(args) == 1:
      inputFile=args[0]
    # No arguments
    else:
      if options.input == '':
        print "Must specify an iterDB file"
        return
      else:
        inputFile=options.input

    myidb=readIterDB.iterDB(inputFile)

    if options.doPrint:
      myidb.printDbInfo()
      return

    allFields=options.fields.split(",")

    # Do sanity checks on fields related inputs
    #parser.add_option('-r', '--rename', dest='renameFields',
    #parser.add_option('-c', '--convert', dest='convertFields',
    if len(options.renameFields)>0:
      options.renameFields=options.renameFields.split(",")
      if len(options.renameFields) != len(allFields):
        print "rename option is incompatibile with fields"
        print len(options.renameFields), len(options.fields)
        return
    if len(options.convertFields)>0:
      options.convertFields=options.convertFields.split(",")
      if len(options.convertFields) != len(allFields):
        print "convert option is incompatibile with fields"
        return

    # If output not written, then just write to stdout:
    if options.outputFile == '':
      fp=sys.stdout
    else:
      fp = open(options.outputFile, "w")

    # Need a header to specify the profiles
    if len(options.renameFields)>0:
        fp.writelines("profiles = "+str(options.renameFields).replace("'","")+"\n")
    else:
        fp.writelines("profiles = "+str(allFields).replace("'","")+"\n")

    # ------------------------------------------------------------
    #  Now that the input file parsed and output file initialized,
    #  do the work
    # ------------------------------------------------------------
    # Rho is special because it is the grid for everything.
    # Part of our problem is that iterDB does not include the
    # normalizaiton of rho
    if myidb.variables['rho']['units'][0:1]=='m':
      # Dimensional rho
      myidb.variables['rhoDim']=myidb.variables['rho']
      # We want rho to be dimensionless
      rtmp=myidb.variables['rho']['values']
      myidb.variables['rho']['values']=rtmp/rtmp[len(rtmp)-1]

    if options.mapFile == '':
      writeNoInterp(myidb,fp,allFields,options)
    else:
      writeWithInterp(myidb,fp,allFields,options.mapFile,options)

if __name__ == '__main__': main()
