#!/usr/bin/env python
from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# expromake input parameters
x.add('Z1','1.0')
x.add('Z2','1.0')
x.add('Z3','1.0')
x.add('B_REF','2.0')
x.add('A_RHO','1.0')
x.add('KAPPA','1.0')
x.add('DELTA','0.0')
x.add('TE_AXIS','5.0')
x.add('ALTE','2.0')
x.add('TI_AXIS','5.0')
x.add('ALTI','2.0')
#
x.add('SET_B_REF','0')
x.add('SET_A_RHO','0')
x.add('SET_KAPPA','0')
x.add('SET_DELTA','0')
x.add('SET_TE_AXIS','0')
x.add('SET_ALTE','0')
x.add('SET_TI_AXIS','0')
x.add('SET_ALTI','0')

# Perform the parsing
x.read_input('input.expromake')

x.printmsg()        

sys.exit(x.error)

