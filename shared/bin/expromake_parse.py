#!/usr/bin/env python
from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# expromake input parameters
x.add('Z1','1.0')
x.add('Z2','1.0')
x.add('Z3','1.0')

# Perform the parsing
x.read_input('input.expromake')

x.printmsg()        

sys.exit(x.error)

