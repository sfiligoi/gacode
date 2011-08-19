"""This script gives help for the requested GYROData data object."""

import sys
import pprint
from pyrats.gyro.data import GYROData

if len(sys.argv) == 1:
    print "Please specify a directory containing GYRO files."
else:
    sim = GYROData(sys.argv[1])
    if len(sys.argv) == 2:
        print "Please specify a GYROData method.  The options are:"
        print
        pprint.pprint(sorted(sim.__dict__.keys()))

    else:

        msg_str       = 'Array has dimensions: '
        n_kinetic_str = 'n_kinetic='    + str(sim.profile['n_kinetic'])    +', '
        n_field_str   = 'n_field='      + str(sim.profile['n_field'])      +', '
        n_x_str       = 'n_x='          + str(sim.profile['n_x'])          +', '
        n_n_str       = 'n_n='          + str(sim.profile['n_n'])          +', '
        n_theta_str   = 'n_theta_plot=' + str(sim.profile['n_theta_plot']) +', '
        n_kinetic_str = 'n_kinetic='    + str(sim.profile['n_kinetic'])    +', '
        n_energy_str  = 'n_energy='     + str(sim.profile['n_energy'])     +', '
        n_lambda_str  = 'n_lambda='     + str(sim.profile['n_lambda'])     +', '

        n_time_str    = 'n_time='    + str(sim.t['n_time'])

        if sys.argv[2] == 't':
            for key in sorted(sim.t.keys()):
                print key.rjust(16),type(sim.t[key])
        elif sys.argv[2] == 'freq':
            for key in sorted(sim.freq.keys()):
                print key.rjust(16),type(sim.freq[key])
        elif sys.argv[2] == 'profile':
            for key in sorted(sim.profile.keys()):
                print key.rjust(18),type(sim.profile[key])
        elif sys.argv[2] == 'diff':
            print msg_str+n_kinetic_str+n_field_str+'2, '+n_time_str+')'
        elif sys.argv[2] == 'diff_i':
            print msg_str+n_kinetic_str+n_field_str+'2, '+n_x_str+n_time_str+')'
        elif sys.argv[2] == 'diff_n':
            print msg_str+n_kinetic_str+n_field_str+'2, '+n_n_str+n_time_str+')'
        elif sys.argv[2] == 'gbflux':
            print msg_str+n_kinetic_str+n_field_str+'4, '+n_time_str+')'
        elif sys.argv[2] == 'gbflux_i':
            print msg_str+n_kinetic_str+n_field_str+'4, '+n_x_str+n_time_str+')'           
        elif sys.argv[2] == 'gbflux_n':
            print msg_str+n_kinetic_str+n_field_str+'4, '+n_n_str+n_time_str+')'         
        elif sys.argv[2] == 'moment_u':
            print msg_str+n_theta_str+n_x_str+n_field_str+n_n_str+n_time_str+')'
        elif sys.argv[2] == 'moment_n':
            print msg_str+n_theta_str+n_x_str+n_kinetic_str+n_n_str+n_time_str+')'
        elif sys.argv[2] == 'moment_e':
            print msg_str+n_theta_str+n_x_str+n_kinetic_str+n_n_str+n_time_str+')'
        elif sys.argv[2] == 'moment_v':
            print msg_str+n_theta_str+n_x_str+n_kinetic_str+n_n_str+n_time_str+')'
        elif sys.argv[2] == 'moment_zero':
            print msg_str+n_x_str+n_kinetic_str+'2, '+n_time_str+')'
        elif sys.argv[2] == 'flux_velocity':
            print msg_str+n_energy_str+n_lambda_str+n_kinetic_str+n_field_str+'2, '+n_n_str+n_time_str
        elif sys.argv[2] == 'k_perp_squared':
            print msg_str+n_n_str+n_time_str
        elif sys.argv[2] == 'docstrings':
            help(GYROData)
        else:
            print "Incorrect GYROData method.  The options are:"
            print "t"
            print "freq"
            print "profile"
            print "geometry"
            print "diff"
            print "diff_i"
            print "diff_n"
            print "gbflux"
            print "gbflux_i"
            print "gbflux_n"
            print "moment_u"
            print "moment_n"
            print "moment_e"
            print "moment_v"
            print "moment_zero"
            print "flux_velocity"
            print "k_perp_squared"
            print "docstrings (not a method, instead prints the docstrings for GYROData)"          
