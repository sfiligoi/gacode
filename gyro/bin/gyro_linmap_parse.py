#!/usr/bin/env python
from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# lin_linmap input parameters

x.add('DET_TOLERANCE','1e-03')
x.add('ERROR_TOLERANCE','1e-03')
x.add('WI_STABLE_LIMIT','1e-02')

# Deprecated parameters

x.dep('MYPARAM','Deprecated.')
  
# Perform the parsing
x.read_input('input.linmap')

x.printmsg()        

sys.exit(x.error)


