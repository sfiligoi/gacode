"""This script gives help for the requested GYROData data object."""

import sys
from pyrats.gyro.data import GYROData

sim1 = GYROData(sys.argv[1])

if sys.argv[2] == 't':
    print sorted(sim1.t.keys())
elif sys.argv[2] == 'freq':
    print sorted(sim1.freq.keys())
elif sys.argv[2] == 'diff':
    print "Array has dimensions of: (n_kinetic, n_field, 2, n_time)"
elif sys.argv[2] == 'diff_i':
    print "Array has dimensions of: (n_kinetic, n_field, 2, n_x, n_time)"
elif sys.argv[2] == 'diff_n':
    print "Array has dimensions of: (n_kinetic, n_field, 2, n_n, n_time)"
elif sys.argv[2] == 'gbflux':
    print "Array has dimensions of: (n_kinetic, n_field, 4, n_time)"
elif sys.argv[2] == 'gbflux_i':
    print "Array has dimensions of: (n_kinetic, n_field, 4, n_x, n_time)"
elif sys.argv[2] == 'gbflux_n':
    print "Array has dimensions of: (n_kinetic, n_field, 4, n_n, n_time)"
elif sys.argv[2] == 'moment_u':
    print "Array has dimensions of: (2, n_theta_plot, n_x, n_field, n_n, n_time)"
elif sys.argv[2] == 'moment_n':
    print "Array has dimensions of: (2, n_theta_plot, n_x, n_kinetic, n_n, n_time)"
elif sys.argv[2] == 'moment_e':
    print "Array has dimensions of: (2 ,n_theta_plot, n_x, n_kinetic, n_n, n_time)"
elif sys.argv[2] == 'moment_v':
    print "(2, n_theta_plot, n_x, n_kinetic, n_n, n_time)"
elif sys.argv[2] == 'moment_zero':
    print "Array has dimensions of: (n_x, n_kinetic, n_moment, n_time)"
elif sys.argv[2] == 'flux_velocity':
    print "Array has dimensions of: (n_energy, n_lambda, n_kinetic, n_field, 2, n_n, n_time)"
elif sys.argv[2] == 'k_perp_squared':
    print "Array has dimensions of: (n_n, n_time)"
elif sys.argv[2] == 'docstrings':
    help(GYROData)
