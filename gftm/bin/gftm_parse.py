from gacodeinput import *
import sys
import gftm_defaults

# fill x with the GFTM defaults
x = gftm_defaults.set_defaults()

x.set_extension('.gen')

# Perform the parsing
x.read_input('input.gftm')

x.printmsg()        

sys.exit(x.error)

