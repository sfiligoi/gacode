#! /usr/bin/env python

"""examples.py

    Example scripts for running pyrats.

    Usage:

    >>>from pyrats.examples import test_tgyro
    >>>test_tgyro()

    Contents:

    test_tgyro - sanity test that creates treg01 object, prints info
    tgyro_stab - diagnostics from tgyro stability mode, sreg02

"""

def test_tgyro(sim='$GACODE_ROOT/tgyro/tools/input/treg01'):
    """Create TGYROData object from sim and print diagnostics.

    Keywords:
    sim -- simulation directory to try, default is treg01
    """
        

    from pyrats.tgyro.data import TGYROData

    treg01 = TGYROData(sim)

    print "Testing TGYROData"
    print "directory = ", treg01.directory_name

    n_iter = treg01.n_iterations

    print "Total iterations:", n_iter

    print "Iteration --- Global Residual --- Number of Calls to Flux Driver"

    for a in range(n_iter):
        print '   ', a, '\t\t', treg01.get_global_res(a), '\t\t\t', treg01.get_flux_count(a)

    print "---- done test_tgyro ---"

def tgyro_stab(sim='$GACODE_ROOT/tgyro/tools/input/sreg02', radius=0):
    """Output linear stability analysis data.

    Prints most unstable ion and electron direction modes,
    and full spectra at specified radial point.

    Keywords:
    sim      -- simulation directory, Default: sreg02
    radius   -- radial index to print full spectrum, Default: 0
                should be an integer between 0 and n_radial-1
    """

    from numpy import array_str
    from pyrats.tgyro.data import TGYROData

    tgyro_data = TGYROData(sim)

    r = tgyro_data.get_r()

    print
    print "------ TGYRO linear stability analysis  ------"
    print
    print "file name:", tgyro_data.directory_name

    # Most unstable modes

    print "r/a     most unstable ion mode"
    print "          (k,  omega,  gamma)"
    for r_i in range(tgyro_data.n_radial):
         print r[r_i], '     ', tgyro_data.get_most_unstable_at_radius(r_i, 'ion')

    print
    print "r/a     most unstable electron mode"
    print "          (k,  omega,  gamma)"
    for r_i in range(tgyro_data.n_radial):
         print r[r_i], '     ', tgyro_data.get_most_unstable_at_radius(r_i, 'elec')

    # Full spectrum

    print
    print "growth spectra at radius:", r[radius]
    print
    ki, gamma_i = tgyro_data.get_stability_at_radius(radius, 'i', 'ion')
    omega_i = tgyro_data.get_stability_at_radius(radius, 'r', 'ion')[1]
    ke, gamma_e = tgyro_data.get_stability_at_radius(radius, 'i', 'elec')
    omega_e = tgyro_data.get_stability_at_radius(radius, 'r', 'elec')[1]

    print "ions"
    print "k    ", array_str(ki, precision=3)
    print "gamma", array_str(gamma_i, precision=3)
    print "omega", array_str(omega_i, precision=3)

    print
    print "electrons"
    print "k    ", array_str(ki, precision=3)
    print "gamma", array_str(gamma_e, precision=3)
    print "omega", array_str(omega_e, precision=3)

