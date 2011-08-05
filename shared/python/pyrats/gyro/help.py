"""This script gives help for the requested GYROData data object."""

import sys
from pyrats.gyro.data import GYROData

if len(sys.argv) == 1:
    print "Please specify a directory containing GYRO files."
else:
    sim1 = GYROData(sys.argv[1])
    if len(sys.argv) == 2:
        print "Please specify a GYROData method.  The options are:"
        print "t"
        print "freq"
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
    else:
        if sys.argv[2] == 't':
            for key in sorted(sim1.t.keys()):
                print key
        elif sys.argv[2] == 'freq':
            for key in sorted(sim1.freq.keys()):
                print key
        elif sys.argv[2] == 'diff':
            print "Array has dimensions of: (n_kinetic=" + str(sim1.profile['n_kinetic']) + ", n_field=" + str(sim1.profile['n_field']) + ", 2, n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'diff_i':
            print "Array has dimensions of: (n_kinetic=" + str(sim1.profile['n_kinetic']) + ", n_field=" + str(sim1.profile['n_field']) + ", 2, n_x=" + str(sim1.profile['n_x']) + ", n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'diff_n':
            print "Array has dimensions of: (n_kinetic=" + str(sim1.profile['n_kinetic']) + ", n_field=" + str(sim1.profile['n_field']) + ", 2, n_n=" + str(sim1.profile['n_n']) + ", n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'gbflux':
            print "Array has dimensions of: (n_kinetic=" + str(sim1.profile['n_kinetic']) + ", n_field=" + str(sim1.profile['n_field']) + ", 4, n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'gbflux_i':
            print "Array has dimensions of: (n_kinetic=" + str(sim1.profile['n_kinetic']) + ", n_field=" + str(sim1.profile['n_field']) + ", 4, n_x=" + str(sim1.profile['n_x']) + ", n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'gbflux_n':
            print "Array has dimensions of: (n_kinetic=" + str(sim1.profile['n_kinetic']) + ", n_field=" + str(sim1.profile['n_field']) + ", 4, n_n=" + str(sim1.profile['n_n']) + ", n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'moment_u':
            print "Array has dimensions of: (n_theta_plot=" + str(sim1.profile['n_theta_plot']) + ", n_x=" + str(sim1.profile['n_x']) + ", n_field=" + str(sim1.profile['n_field']) + ", n_n=" + str(sim1.profile['n_n']) + ", n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'moment_n':
            print "Array has dimensions of: (n_theta_plot=" + str(sim1.profile['n_theta_plot']) + ", n_x=" + str(sim1.profile['n_x']) + ", n_kinetic=" + str(sim1.profile['n_kinetic']) + ", n_n=" + str(sim1.profile['n_n']) + ", n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'moment_e':
            print "Array has dimensions of: (n_theta_plot=" + str(sim1.profile['n_theta_plot']) + ", n_x=" + str(sim1.profile['n_x']) + ", n_kinetic=" + str(sim1.profile['n_kinetic']) + ", n_n=" + str(sim1.profile['n_n']) + ", n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'moment_v':
            print "Array has dimensions of: (n_theta_plot=" + str(sim1.profile['n_theta_plot']) + ", n_x=" + str(sim1.profile['n_x']) + ", n_kinetic=" + str(sim1.profile['n_kinetic']) + ", n_n=" + str(sim1.profile['n_n']) + ", n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'moment_zero':
            print "Array has dimensions of: (n_x=" + str(sim1.profile['n_x']) + ", n_kinetic=" + str(sim1.profile['n_kinetic']) + ", n_moment=" + str(sim1.profile['n_moment']) + ", n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'flux_velocity':
            print "Array has dimensions of: (n_energy=" + str(sim1.profile['n_energy']) + ", n_lambda=" + str(sim1.profile['n_lambda']) + ", n_kinetic=" + str(sim1.profile['n_kinetic']) + ", n_field=" + str(sim1.profile['n_field']) + ", 2, n_n=" + str(sim1.profile['n_n']) + ", n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'k_perp_squared':
            print "Array has dimensions of: (n_n=" + str(sim1.profile['n_n']) + ", n_time=" + str(sim1.t['n_time']) + ")"
        elif sys.argv[2] == 'docstrings':
            help(GYROData)
        else:
            print "Incorrect GYROData method.  The options are:"
            print "t"
            print "freq"
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
