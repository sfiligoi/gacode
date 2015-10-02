from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# VGEN input parameters
x.add('ER_METHOD','1')
x.add('VEL_METHOD','2')
x.add('ERSPECIES_INDX','2')
x.add('NTHETA_MIN','16')
x.add('NTHETA_MAX','22')
x.add('TYPE','DC')

# Perform the parsing
x.read_input('input.vgen')

x.printmsg()        

sys.exit(x.error)


