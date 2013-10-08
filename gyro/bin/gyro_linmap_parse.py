#!/usr/bin/env python
from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# lin_linmap input parameters

x.add('DET_TOLERANCE','1e-03')
x.add('ERROR_TOLERANCE','1e-03')
x.add('WI_STABLE_LIMIT','1e-02')

x.add('startWR','-2.0E-01')
x.add('startWI','2.0E-01')
x.add('RADIUSstart','0.5')
x.add('RADIUSmin','0.4')
x.add('RADIUSmax','0.6')
x.add('RADIUSstep','0.05')
x.add('L_Ystart','0.5')
x.add('L_Ymin','0.4')
x.add('L_Ymax','0.6')
x.add('L_Ystep','0.05')

# Deprecated parameters

x.dep('MYPARAM','Deprecated.')
  
# Perform the parsing
x.read_input('input.linmap')

x.printmsg()        

sys.exit(x.error)


