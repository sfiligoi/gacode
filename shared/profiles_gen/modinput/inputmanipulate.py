#!/usr/bin/env python
""" 'iterdb2gyro -d' manipulates data in INPUT_profiles

The input manipulation mode replaces a profile 'a' in the file INPUT_profiles file 'A' by the 
profile s * 'b', where the profile 'b' is extracted from file 'B' and scaled by a factor s. The 
result is written in a new file C.
File 'A' can be the same file as 'B' but neither of them can be the same as C. Also, 'a' can be the 
same profile as 'b'.

The plasma composition mode genetartes ion density profiles based on the electron density profile 
found in the the input file, given ion charges, Zeff and if three ions are considered the third to 
first ion density ratio (ni_3/ni_1).

Usage: 

    input manipulaion mode:  
         
        iterdb2gyro -d <options>

    plasma composition mode:
    
        iterdb2gyro -d -comp <options>

optional arguments for the input manipulaion mode:

    -inf1        -- INPUT_profile file to be modified (file A)
                        (Default: INPUT_profiles)
    -inf2        -- the INPUT_profiles file to extract data from (file B)
                        (Default: INPUT_profiles) 
    -outf        -- output file. (file C) 
                        (Default: INPUT_profiles_mod)
    -scal        -- scaling factor (the scaling factor s)
                        (Default: 1.0)
    -pin         -- this data column will be read from the file 'inf2' (profile b) 
                        (Default: ni_1)
    -pout        -- this data column of 'inf1' will be replaced by the new data (profile a)
                        (Default: ni_2)

optional arguments for the plasma composition mode:

    -inf1        -- the input INPUT_profile file 
                        (Default: INPUT_profiles)
    -outf        -- the output INPUT_profile file
                        (Default: INPUT_profiles_mod)
    -Z1          -- Charge of the first ion species (not necessarily integer)
                        (Default: 1)
    -Z2          -- Charge of the second ion species (not necessarily integer)
                        (Default: 6)
    -Z3          -- Charge of the third ion species (not necessarily integer)
                        (Default: 2)
    -Zeff        -- Effective ion charge
                        (Default: 2.0)
    -n3n1        -- ni_3 / ni_1 
                        (Default: 0.  ; no third ion species)


Examples:

    iterdb2gyro -d -inf1 firstdir/INPUT_profiles -inf2 seconddir/INPUT_profiles -scal 0.85 -pin Te -pout Ti_2

    The command above will generate the file 'INPUT_profiles_mod' which 
    contains the same data as 'firstdir/INPUT_profiles' except that the 
    Ti_2 profile is equal to 0.85*Te, where Te is from 
    'seconddir/INPUT_profiles'.

    iterdb2gyro -d -comp -Z2 8 -Zeff 2.2 -n3n1 0.1 -outf INPUT_o

    Based on the electron density profile ni_1, ni_2 and ni_3 are "
    generated for a plasma with Zeff=2.2, Z1=1, Z2=8, Z3=4 and"
    ni_3 / ni_1 = 0.1. The new data is written in 'INPUT_ox'.

NOTES: 

    The the input file(s) and the output file cannot be the same! 

    In plasma composition mode it is not controlled that the parameters are consistent (wrong input values
    can lead to negative densities).

    It is assumed that an INPUT_profiles file contains 8 data blocks with five data colums which are 
    16 characters wide. The number of blocks and the column width can be set in insertcol.sh and 
    extractcol.sh. Furthermore it is assumed that the lines of scalar data are not started with space or "-".

(I. Pusztai ; 10/14/2010)
"""
import os
import sys
import getopt
import shlex
import subprocess
import math

# Variables and default values
#########################################################################
path="${GYRO_DIR}/tools/modinput"

params = {
    'inputfile1': 'INPUT_profiles',         # the INPUT_profiles file to be modified
    'inputfile2': 'INPUT_profiles',         # the INPUT_profiles file to extract data from
    'outputfile': 'INPUT_profiles_mod',     # the output file
    'scaling': '1.0',                       # scaling factor multiplying the data column
    'paramin': 'ni_1',                      # input parameter in 'inputfile1'
    'paramout': 'ni_2',                     # output parameter in 'outputfile'
    'extractname': 'extract_',              # name of the extract column file
    'paramout': 'ni_2',
    'compflag': '0', 
    'Z1': '1', 
    'Z2': '6', 
    'Z3': '2', 
    'Zeff': '2.0',
    'n3n1': '0.0'}              

def usage():
    """Print the doc string."""
    print __doc__

def get_parameters():
    """Get command line options and arguments."""
    global params

       # Parse command-line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:i:j:o:s:p:r", \
                  ["modhelp","inf1=", "inf2=", "outf=", "scal=", "pin=",\
                   "pout=", "comp", "Z1=", "Z2=", "Z3=", "Zeff=", "n3n1="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--modhelp"):
            usage()
            sys.exit()
        elif opt in ("-i", "--inf1"):
            params['inputfile1'] = arg
        elif opt in ("-j", "--inf2"):
            params['inputfile2'] = arg
        elif opt in ("-o","--outf"):
            params['outputfile'] = arg
        elif opt in ("-s", "--scal"):
            params['scaling'] = arg
        elif opt in ("-p", "--pin"):
            params['paramin'] = arg
        elif opt in ("-r", "--pout"):
            params['paramout'] = arg
        if opt in ("--comp"):
            params['compflag'] = '1'
        elif opt in ("--Z1"):
            params['Z1'] = arg
        elif opt in ("--Z2"):
            params['Z2'] = arg
        elif opt in ("--Z3"):
            params['Z3'] = arg
        elif opt in ("--Zeff"):
            params['Zeff'] = arg
        elif opt in ("--n3n1"):
            params['n3n1'] = arg

def print_variable_values(pipe=sys.stdout):
    """Print the variables' values to a pipe.
    
    Keyword arguments:
    pipe  --  where to print the info, eg a file, default: sys.stdout

    """
    print >> pipe
    print >> pipe, "First input file          = " + params['inputfile1']
    print >> pipe, "Second input file         = " + params['inputfile2']
    print >> pipe, "Output file               = " + params['outputfile']
    print >> pipe, "Scaling factor            = " + params['scaling']
    print >> pipe, "Input parameter           = " + params['paramin']
    print >> pipe, "Output parameter          = " + params['paramout']
    print >> pipe

def print_variable_values_comp(pipe=sys.stdout):
    """Print the variables' values to a pipe in plasma composition mode.
    
    Keyword arguments:
    pipe  --  where to print the info, eg a file, default: sys.stdout

    """
    print >> pipe
    print >> pipe, "iterdb2gyro -d Plasma composition mode"
    print >> pipe, "Input file   = " + params['inputfile2']
    print >> pipe, "Output file  = " + params['outputfile']
    print >> pipe, "Zeff         = " + params['Zeff']
    print >> pipe, "Z1           = " + params['Z1']
    print >> pipe, "Z2           = " + params['Z2']
    if params['n3n1'] != '0.0':
        print >> pipe, "Z3           = " + params['Z3']
        print >> pipe, "ni_3/ni_1    = " + params['n3n1']

    print >> pipe

def extractcolsh():
    """ Extracts a given data column from an INPUT_profiles file using the script 'extractcol.sh'. """
    
    os.system('sh ' + path + '/extractcol.sh ' + params['paramin'] + ' ' + params['inputfile2'])

def scalecol():
    f = open(params['extractname'] + params['paramin'], 'r') 
    os.system('echo > ' + params['extractname'] + params['paramin'] + '_n')
    nf = open(params['extractname'] + params['paramin'] + '_n', 'w') 
    datacol = f.readlines()
    for line in datacol: 
        achar = line.strip()           
        a = float(achar),
        b = a[0]
        scale=float(params['scaling'])
        b = scale * b
        result = str('% 1.7E\n' % b)
        nf.write(result)

    f.close()
    nf.close()

def insertcolsh():
    """ Inserts a given data column into an INPUT_profiles file using the script 'insertcol.sh'. """

    os.system('sh ' + path + '/insertcol.sh ' + params['paramout'] + ' ' + params['extractname'] + \
                   params['paramin'] + '_n ' + params['inputfile1'] + ' ' + params['outputfile'])

def calculate_n1():
    """ Calculates the first ion density and sets the scaling factor """
    Z1 = float(params['Z1'])
    Z2 = float(params['Z2'])
    Z3 = float(params['Z3'])
    Zeff = float(params['Zeff'])
    C = float(params['n3n1'])
    scalenum = (Zeff - Z2) / (Z1*Z1 + C * Z3*Z3 - Z1 * Z2 - C * Z2 * Z3)

    params['scaling']=str('%1.8E\n' % scalenum)

def calculate_n2():
    """ Calculates the second ion density and sets the scaling factor """
    Z1 = float(params['Z1'])
    Z2 = float(params['Z2'])
    Z3 = float(params['Z3'])
    Zeff = float(params['Zeff'])
    C = float(params['n3n1'])
    scalenum =  (Zeff * (Z1 + C * Z3) - Z1*Z1 - C * Z3*Z3) / (Z2 * (Z1 * Z2 + C * Z2 * Z3 - Z1*Z1 - C * Z3*Z3))

    params['scaling']=str('%1.8E\n' % scalenum)

def calculate_n3():
    """ Calculates the third ion density and sets the scaling factor """
    Z1 = float(params['Z1'])
    Z2 = float(params['Z2'])
    Z3 = float(params['Z3'])
    Zeff = float(params['Zeff'])
    C = float(params['n3n1'])
    scalenum = C * (Zeff - Z2) / (Z1*Z1 + C * Z3*Z3 - Z1 * Z2 - C * Z2 * Z3)

    params['scaling']=str('%1.8E\n' % scalenum)

def main():
    """ main part of the code """
    get_parameters()
    if params['compflag'] == '0':      # input profile manipulation mode
        if (os.path.isfile(params['inputfile1']) and os.path.isfile(params['inputfile2'])): # do the input files exist? 
            print_variable_values()    # print info
            extractcolsh()             # extract the data 'pin' from the file 'inf2'
            scalecol()                 # do the scaling
            insertcolsh()              # replace the data 'pout' in the file 'inf1'
            os.system('rm ' + params['extractname'] + params['paramin'])  # delete temporary files
            os.system('rm ' + params['extractname'] + params['paramin'] + '_n')
        else:                          # in case of missing input
            print "ERROR:    Input file does not exist."
    else:                              # plasma composition mode
        params['inputfile2']=params['inputfile1']      # there is only one input file
        if os.path.isfile(params['inputfile2']):       # check if the input file exists
            print_variable_values_comp()               # print info             
            finaloutputfile = params['outputfile']     # saves the output file name into a variable
            params['paramin'] = 'ne'       # The ion densities are calculated from the electron density
# first ion
            calculate_n1()                 # Modifies the scaling factor for the proper value for ni_1
            extractcolsh()
            scalecol()                     # Does the scaling for ni_1
            params['paramout'] = 'ni_1'    # Parameter to be replaced
            params['outputfile'] = 'INPUT_profile_ni_1'  # temporary output with the new ni_1
            insertcolsh()
            os.system('rm ' + params['extractname'] + params['paramin'])  # remove temporary files
            os.system('rm ' + params['extractname'] + params['paramin'] + '_n')
# second ion                              
            params['inputfile1'] = 'INPUT_profile_ni_1'  # here the temporary file will be the new inputfile
            calculate_n2()                 # Modifies the scaling factor for the proper value for ni_2
            extractcolsh()
            scalecol()                     # Does the scaling for ni_2 
            params['paramout'] = 'ni_2'    # Parameter to be replaced 
            params['outputfile'] = 'INPUT_profile_ni_2'  # temporary output with new ni_1 and ni_2
            insertcolsh()
            os.system('rm ' + params['extractname'] + params['paramin'])  # remove temporary files
            os.system('rm ' + params['extractname'] + params['paramin'] + '_n')
# third ion                              
            params['inputfile1'] = 'INPUT_profile_ni_2'  # here the temporary file will be the new inputfile
            calculate_n3()               # Modifies the scaling factor for the proper value for ni_3
            extractcolsh()
            scalecol()                   # Does the scaling for ni_3 
            params['paramout'] = 'ni_3'  # Parameter to be replaced 
            params['outputfile'] = finaloutputfile       # final output file with the new densities
            insertcolsh()
            os.system('rm ' + params['extractname'] + params['paramin'])  # remove temporary files
            os.system('rm ' + params['extractname'] + params['paramin'] + '_n')
            os.system('rm INPUT_profile_ni_1')
            os.system('rm INPUT_profile_ni_2')
        else:
            print "ERROR:    Input file does not exist."

if __name__ == "__main__":
    main()

