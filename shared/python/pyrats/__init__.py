#! /usr/bin/env python
# pyrats.py:
"""PYRATS:  PYthon Routines for Analyzing Transport Simulations

    This package provides routines for analyzing gacode transport
    simulations. It is intended to interface simulation data
    with numpy and matplotlib, free python-based plotting and 
    analysis tools, which are required for these tools to work.

    For more information see:

    python:        http://www.python.org/
    numpy:         http://numpy.scipy.org/
    matplotlib:    http://matplotlib.sourceforge.net/

    -------------------
    Author(s):

    Luc Peterson (jdpeters@princeton.edu)

    -------------------
    Version History:

    0.1 -- 4/28/2011 -- Initial attempt

    -------------------
    Contents:

    data
        TGYROData      TGYRO Data Object
    examples
        tgyro_test.py  test for treg01

    -------------------
    Example Usage:

    Basic Test:

    python
    >>>from pyrats.examples import test_tgyro
    >>>tgyro_test()

    Plotting Test:
    python
    >>>from matplotlib import pyplot
    >>>from pyrats.data import TGYROData
    >>>sim1 = TGYROData('$GACODE_ROOT/tgyro/tools/input/treg01')
    >>>pyplot.plot(sim1.get_r(), sim1.get_Te())
    >>>pyplot.show()

    For Help:
    python
    >>>from pyrats import *
    >>>help(data)
    >>>help(examples)

"""

# Package Contents
__all__ = ["data", "examples"]

# Test for required modules
modules = ["numpy", "matplotlib"]

import imp

for module in modules:
    try:
        imp.find_module(module)
    except ImportError:
        print "Error: pyrats needs the module", module
        print "Make sure this is installed correctly."
