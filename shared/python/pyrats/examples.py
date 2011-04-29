#! /usr/bin/env python

"""examples.py

    Example scripts for running pyrats.

    Usage:

    >>>from pyrats.examples import test_tgyro
    >>>test_tgyro

    Contents:

    test_tgyro - sanity test that creates treg01 object, prints info

"""

def test_tgyro(sim='$GACODE_ROOT/tgyro/tools/input/treg01'):
    """Create TGYROData object from sim and print diagnostics.

    Keywords:
    sim -- simulation directory to try, default is treg01
    """
        

    from pyrats.data import TGYROData

    treg01 = TGYROData(sim)

    print "Testing TGYROData"
    print "directory = ", treg01.directory_name

    n_iter = treg01.n_iterations

    print "Total iterations:", n_iter

    print "Iteration --- Global Residual --- Number of Calls to Flux Driver"

    for a in range(n_iter):
        print a, treg01.get_global_res(a), treg01.get_flux_count(a)

    print "---- done test_tgyro ---"
