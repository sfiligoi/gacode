# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
import sys
from pyrats.neo.data import NEOData
import matplotlib.pyplot as plt
import math

#First are several functions that perform error checking that comes up frequently.
def opendir(directory):
    try:
        sim = NEOData(directory)
    except IOError:
        print
