#!/usr/bin/env python
"""transp2gyro.py: Create input.profiles from TRANSP server.

transp2gyro.py writes a trxpl script (.ind) from user input, which
it uses to connect to a TRANSP server and write a plasma 
state file (.cdf) and an efit G-EQDSK file (.geq). transp2gyro then
calls profiles_gen to produce input.profiles (and, optionally, 
input.profiles.geo), for use with gacode (GYRO, etc.).

Usage: transp2gyro.py <options>

Optional Arguments:

    -k (--tok) tokamak       -- the TRANSP id for the tokamak (eg NSTX)
    -y (--year) run_year     -- the run year for the shot (eg 08)
    -i (--shot) transp_id    -- TRANSP run id in tree (eg 129354A01)
    -t (--time) time         -- time slice of interest in ms
    -a (--avg) ave_time      -- time-averaging in ms
    -o (--trout) output_name -- name of output folder
                                   (default: "machine_transp_id_time")
    --res                    -- don't use default resolution parameters
    -g (--geo)               -- also create input.profiles.geo for
                                   numeric equilibrium
    --prgen arguments        -- pass arguments to profiles_gen formatted as:
                                   arg1;arg2;arg3;arg4
    -h (--trhelp)            -- display this message

File flow:
    TRANSP -> trxpl -> plasma state and efit files ->
    prgen -> input.profiles

Output:
    output_name/                  -- directory with shot files
        input.profiles            -- gacode-format input.profiles
        input.profiles.geo        -- numeric equilibrium file (optional)
        input.profiles.log        -- log file from profiles_gen
        output_name.log           -- log file from transp2gyro.py
        trxpl.log                 -- log file from trxpl
        transp2gyro_files/        -- directory with helper files
            output_name.cdf       -- plasma state file
            output_name.ind       -- trxpl script
            output_name.mdescr    -- machine description file
            output_name.sconfig   -- shot configuration file
            output_name.geq       -- shot geq file

Examples:
    transp2gyro.py -g -k NSTX -y 07 --shot 124901A02 -t 300 -a 2
    transp2gyro.py -r
    transp2gyro.py --tokamak ITER --year 06 --shot 20100P07

Common Errors:
    - trxpl and/or profiles_gen are not executible in user's $PATH
        - Make sure user has access to MDS+ server defined by 
          $MDS_TRANSP_SERVER
        - At GA see:
          http://www.pppl.gov/~jpeterso/documents/fusion_grid_access_ga.notes
        - At PPPL see:
          http://www.pppl.gov/~jpeterso/documents/fusion_grid_access_pppl.notes
    - Requested TRANSP id is not in database
    - Specified same output_name for different shots
        - transp2gyro.py will not overwrite a pre-existing directory, but
          will create top-level output_name.log before quitting.
    - Cannot create input.profiles.geo (option -g errors)
        - make sure gato.f is compiled correctly and 'profiles_gen -g' works
    - Python errors
        - you need to be running at least python version 2.5

--------
Author: Luc Peterson (jpeterso@pppl.gov)

Version History:
    0.1  --  July 28, 2009   : Parses input and writes trxpl script
    0.2  --  August 11, 2009 : Runs trxpl and iterdb2gyro
    0.3  --  August 12, 2009 : Cleaned up, debugged log files and input
    0.4  --  August 20, 2009 : Added INPUT_profiles.geo option
    0.5  --  August 26, 2009 : Added trxpl.log and changed for MDS+ access
    0.6  --  October 8, 2009 : Executible from iterdb2gyro
    0.7  --  April 22, 2011  : Use gacode name conventions (profiles_gen)

"""
import os
import sys
import getopt
import subprocess

# Variables and default values
#########################################################################

params = { 
    'machine': '',                   # Machine name (eg NSTX)
    'run_year': '',                  # Shot year (eg 08)
    'transp_id': '',                 # TRANSP run identification 
                                     #    (eg 129354A01)
    'time': '',                      # Time slice of interest in ms 
                                     #    (eg 240)
    'ave_time': '',                  # Smoothing time for data in ms 
                                     #    (eg 2)
    'output_name': '',               # Output file name
                                     #    default: machine_transp_id_time
    'theta_points': '128',           # Num theta points for 2D spline
    'r_points': '128',               # Num r points for overlay grid
    'z_points': '128',               # Num z points for overlay grid
    'B_tor_orientation': '0',        # Toroidal Magnetic Field Orient.
                                     #     1: counter-clockwise
                                     #    -1: clockwise
                                     #     0: read from TRANSP
    'I_pol_orientation': '0',        # Poloidal Plasma Current Orient.
    'number_of_radial_zones': '50',  # Number of radial zones in state 
                                     #    profile
    'use_default_resolutions': True, # Whether to use default 
                                     #    resolutions
    'create_numeric_equil': False,   # Whether to make input.profiles.geo
    'I_name': 'from TRANSP',         # Name of I_pol flag
    'B_name': 'from TRANSP',         # Name of B_tor flag
    'prgenargs': ''}                 # Arguments passed to profiles_gen
         
##########################################################################

# Functions
##########################################################################

def usage():
    """Print the doc string."""
    print __doc__

def print_welcome_message(pipe=sys.stdout):
    """Print version info to a pipe.
    
    Keyword arguments:
    pipe  --  where to print the info, eg a file, default: sys.stdout

    """
    print >> pipe, "-----------------------------------------------------"
    print >> pipe, "transp2gyro.py  : J. Luc Peterson (jpeterso@pppl.gov)"
    print >> pipe, "Version 0.7"
    print >> pipe, "-----------------------------------------------------"

def print_variable_values(pipe=sys.stdout):
    """Print the variables' values to a pipe.
    
    Keyword arguments:
    pipe  --  where to print the info, eg a file, default: sys.stdout

    """
    print >> pipe
    print >> pipe, "Machine                   = " + params['machine']
    print >> pipe, "Run Year                  = " + params['run_year']
    print >> pipe, "TRANSP Run ID             = " + params['transp_id']
    print >> pipe, "Output Directory          = " + params['output_name']
    print >> pipe, "Time (ms)                 = " + params['time']
    print >> pipe, "Averaging Time (ms)       = " + params['ave_time']
    print >> pipe, "Number of theta pts       = " + params['theta_points']
    print >> pipe, "Number of r pts           = " + params['r_points']
    print >> pipe, "Number of z points        = " + params['z_points']
    print >> pipe, "B_toroidal orientation    = " + params['B_name']
    print >> pipe, "I_poloidal orientation    = " + params['I_name']
    print >> pipe, "Number of radial zones    = " + \
                   params['number_of_radial_zones']
    print >> pipe, "Create input.profiles.geo = " + \
                   str(params['create_numeric_equil'])
    print >> pipe

def write_trxpl_script():
    """Write the trxpl script as params["output_name"].ind."""
    # Open file
    script_name = params["output_name"] + "/" + params['output_name'] + \
                  ".ind"
    script_file = open(script_name, 'w')

    # Times are passed in seconds
    time = str(float(params['time']) / 1000.0)
    avg = str(float(params['ave_time']) / 1000.0)

    # Write the data
    script_file.write('; ' + script_name + '\n')
    script_file.write('; Invoke with trxpl @' + script_name + '\n');
    script_file.write('; Created by transp2gyro.py\n')
    script_file.write('*S*\n')
    script_file.write('m\nd\nd\n')
    script_file.write(params['machine'] + '.')
    script_file.write(params['run_year'] + '\n')
    script_file.write(params['transp_id'] + '\n')
    script_file.write('a\n')
    script_file.write(time + '\n')
    script_file.write(avg + '\n')
    script_file.write(params['theta_points'] + '\n')
    script_file.write(params['r_points'] + '\n')
    script_file.write(params['z_points'] + '\n')
    script_file.write(params['B_tor_orientation'] + '\n')
    script_file.write(params['I_pol_orientation'] + '\n')
    script_file.write('y\nx\nN\n')
    script_file.write(params['number_of_radial_zones'] + '\n')
    script_file.write('W\n')
    script_file.write(params['output_name'] + '\n')
    script_file.write('q\nq\n')

    # Done!
    script_file.close()

def get_parameters():
    """Get command line options and arguments."""
    global params

    # Parse command-line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hk:y:i:o:t:a:g:", \
                     ["trhelp", "tokamak=", "year=", "shot=", "trout=", \
                      "time=", "aveg=", "res", "geo=", "prgen="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--trhelp"):
            usage()
            sys.exit()
        elif opt in ("-k", "--tokamak"):
            params['machine'] = arg
        elif opt in ("-y", "--year"):
            params['run_year'] = arg
        elif opt in ("-i", "--shot"):
            params['transp_id'] = arg
        elif opt in ("-o", "--trout"):
            params['output_name'] = arg
        elif opt in ("-t", "--time"):
            params['time'] = arg
        elif opt in ("-a", "--avg"):
            params['ave_time'] = arg
        elif opt in ("--res"):
            params['use_default_resolutions'] = False
        elif opt in ("-g", "--geo"):
            params['create_numeric_equil'] = True
        elif opt in ("--prgen"):
            params['prgenargs'] = params['prgenargs'] + arg

    # Prompt user for data if not given as arguments

    # Identification Information
    if not params['machine']:
        params['machine'] = raw_input("Please enter the machine ID: ")
    if not params['run_year']:
        params['run_year'] = raw_input("Please enter the shot run year: ")
    if not params['transp_id']:
        params['transp_id'] = raw_input("Please enter the TRANSP run ID: ")

    # Get time slice information
    if not params['time']:
        params['time'] = raw_input("Please enter the time of interest " + \
                                   "in ms: ")
    if not params['ave_time']:
        params['ave_time'] = raw_input("Please enter the averaging " + \
                                       "time in ms: ")

    # Create default string for file names
    if not params['output_name']:
        params['output_name'] = params['machine'] + '_' + \
                                params['transp_id'] + '_' + params['time']

    # Get resolutions if not using defaults
    success = False
    while not success:
          answer = raw_input("Are default mapping resolutions ok? (y/n) ")
          if answer == 'y' or answer == 'yes' or answer == 'Y':
              success = True
          elif answer == 'n' or answer == 'no' or answer == 'N':
              params['use_default_resolutions'] = False
              success = True
          else:
              print "Error: Please enter 'y' or 'n'"

    if not params['use_default_resolutions']:
        params['theta_points'] = raw_input("Please enter the number " + \
                                           "of theta points for 2D " + \
                                           "spline: ")
        params['r_points'] = raw_input("Please enter the number of " + \
                                       "radial points for grid overlay: ")
        params['z_points'] = raw_input("Please enter the number of " + \
                                       "z points for grid overlay: ")
        success = False
        while not success:
            params['B_tor_orientation'] = raw_input("B_toroidal " + \
                                          "orientation:\n" + \
                                          "   1: counter-clockwise \n" + \
                                          "  -1: clockwise \n" + \
                                          "   0: read from TRANSP \n" + \
                                          "Selection: ")
            if params['B_tor_orientation'] == '1':
                params['B_name'] = "counter-clockwise"
                success = True
            elif params['B_tor_orientation'] == '-1':
                params['B_name'] = "clockwise"
                success = True
            elif params['B_tor_orientation'] == '0':
                params['B_name'] = "from TRANSP"
                success = True
            else:
                print "Invalid Input Value! Options: [1, -1, 0]"
                success = False

        success = False
        while not success:
            params['I_pol_orientation'] = raw_input("I_poloidal " + \
                                  "orientation:\n" + \
                                  "   1: counter-clockwise \n" + \
                                  "  -1: clockwise \n" + \
                                  "   0: read from TRANSP \n" + \
                                  "Selection: ")
            if params['I_pol_orientation'] == '1':
                params['I_name'] = "counter-clockwise"
                success = True
            elif params['I_pol_orientation'] == '-1':
                params['I_name'] = "clockwise"
                success = True
            elif params['I_pol_orientation'] == '0':
                params['I_name'] = "from TRANSP"
                success = True
            else:
                print "Invalid Input Value! Options: [1, -1, 0]"
                success = False

    # Prompt for Numeric Equilibrium
    if not params['create_numeric_equil']:
        success = False
        while not success:
            entry = raw_input("Shall I make input.profiles.geo for " + \
                              "numeric equilibria? (y/n) ")
            if entry == 'y' or entry == 'Y' or entry == 'yes':
                params['create_numeric_equil'] = True
                success = True
            elif entry == 'n' or entry == 'N' or entry == 'no':
                params['create_numeric_equil'] = False
                success = True
            else:
                print 'Invalid entry! Options: [y, n]'


def run_profiles_gen(log):
    """Call profiles_gen to create input.profiles from plasma state.

    Keyword Arguments:
    log  --  file handle to print the output of profiles_gen

    """
    # Files we're working with
    directory = params["output_name"]
    plasma_state_file = directory + ".cdf"
    g_file = directory + ".geq"

    # Create list of arguments to pass to profiles_gen
    prgen_command = "profiles_gen;-i;" + plasma_state_file + params['prgenargs']
    if params['create_numeric_equil']:
        prgen_command = prgen_command + ";-g;" + g_file

    # Run profiles_gen and write output to log
    try:
        prgen = subprocess.Popen(prgen_command.split(";"), \
                                  stdout = subprocess.PIPE, \
                                  stderr = subprocess.STDOUT, \
                                  cwd = directory)
        prgen.wait()
        output = prgen.communicate()[0]
        log.write(output)
    except OSError, e:
        msg = "profiles_gen execution failed: " + e
        print_pipe_and_screen(msg, log)

    # Check for success by looking for input.profiles
    if not os.access(directory + "/input.profiles", os.F_OK):
        msg = "Unable to generate input.profiles.\n" + "Quitting..."
        print_pipe_and_screen(msg, log)
        sys.exit(1)

    # Check for success by looking for input.profiles.geo
    if params["create_numeric_equil"]:
        if not os.access(directory + "/input.profiles.geo", os.F_OK):
            msg = "Unable to generate input.profiles.geo.\n" + "Quitting..."
            print_pipe_and_screen(msg, log)
            sys.exit(1)

def run_trxpl_script(log):
    """Call trxpl on .ind script to create plasma state file.

    Writes output of trxpl to trxpl.log.

    Keyword Arguments:
    log    --   The log file for transp2gyro for error messages.

    """
    # Run trxpl and write output of trxpl to trxpl.log
    # Because trxpl can have a lot of output subprocess can
    # lock up (buffer overflow) if not writing to separate file
    command_name = "@" + params["output_name"] + ".ind"
    trxpl_log = open(params["output_name"] + "/trxpl.log", "w")
    try:
        trxpl = subprocess.Popen(["trxpl", command_name], \
                cwd = params["output_name"], stdout = trxpl_log, \
                stderr = trxpl_log)
        trxpl.wait()
    except OSError, e:
        print_pipe_and_screen("trxpl execution failed: " + e, log)
        print_pipe_and_screen("Quitting", log)
        sys.exit(1)

    # Test for success by looking for .cdf and .geq files
    plasma_state_files = [ params["output_name"] + "/" + \
                           params["output_name"] + ".cdf", \
                           params["output_name"] + "/" + \
                           params["output_name"] + ".geq" ]
    for plasma_state_file in plasma_state_files:
        if not os.access(plasma_state_file, os.F_OK):
            msg = "Unable to generate " + plasma_state_file + ".\n" + \
                  "Check shot parameters and TRANSP/trxpl settings.\n" + \
                  "trxpl.log may shed some light on what went wrong.\n" + \
                  "Quitting..."
            print_pipe_and_screen(msg, log)
            sys.exit(1)
    trxpl_log.close()

def does_program_exist(program):
    """Check to see if 'program' is executible in user's path.

    Keyword Arguments:
    program  --  the program to check

    """
    for path in os.environ.get('PATH', '').split(':'):
        if os.access(os.path.join(path, program), os.X_OK):
            return True
    return False

def print_pipe_and_screen(message, pipe):
    """Print message to both screen and a pipe."""
    print message
    print >> pipe, message

def main():
    """Create input.profiles from TRANSP server."""
    print_welcome_message()

    # Test to see if we can run profiles_gen and trxpl
    if not does_program_exist('profiles_gen'):
        print "Unable to execute profiles_gen! Check $PATH..."
        print "Quitting..."
        sys.exit(1)
    if not does_program_exist('trxpl'):
        print "Unable to execute trxpl! Check $PATH..."
        print "Quitting..."
        sys.exit(1)

    # Read in parameters
    get_parameters()
    print_variable_values()

    # Set up log file
    log_name = params["output_name"] + ".log"
    log = open(log_name, "w")
    print_welcome_message(log)
    print_variable_values(log)

    # Make new working directory
    print_pipe_and_screen("Making new directory " + \
                          params["output_name"] + "...", log)
    try:
        subprocess.check_call(['mkdir', params["output_name"]])
    except OSError, e:
        print_pipe_and_screen("Unable to create new directory: " + e, log)
        print_pipe_and_screen("Quitting", log)
        sys.exit(1)
    except subprocess.CalledProcessError, e:
        print_pipe_and_screen(e, log)
        print_pipe_and_screen("Quitting", log)
        sys.exit(1)

    # Write the trxpl script from the parameters
    msg = "Writing trxpl script " + params["output_name"] + ".ind..."
    print_pipe_and_screen(msg, log)
    write_trxpl_script()

    # Run trxpl script, connect to TRANSP server, etc
    msg = "Connecting to TRANSP server w/ trxpl..."
    print_pipe_and_screen(msg, log)
    run_trxpl_script(log)

    # Run profiles_gen to create input.profiles
    print_pipe_and_screen("Running profiles_gen...", log)
    run_profiles_gen(log)

    # Cleanup files
    print_pipe_and_screen("Cleaning files...", log)
    subprocess.Popen(["rmdir", "tmp"], cwd = params["output_name"])
    subprocess.Popen(["mkdir", "transp2gyro_files"], \
                     cwd = params["output_name"])
    docs = [params["output_name"] + ".cdf", \
            params["output_name"] + ".geq", \
            params["output_name"] + ".ind", \
            params["output_name"] + ".mdescr", \
            params["output_name"] + ".sconfig"]
    for doc in docs:
        subprocess.Popen(["mv", doc, "transp2gyro_files"], \
                         cwd = params["output_name"])

    # We made it! Close and move log file.
    print_pipe_and_screen("Success!", log)
    log.close()
    subprocess.Popen(["mv", log_name, params["output_name"]])

if __name__ == "__main__":
    main()
