from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# PROFILES_3d input parameters
x.add('PSI_NORM','0.45')
x.add('NTHETA_FOURIER','6')
x.add('A_METERS','0.6')
x.add('IMPURITY_FLAG','0')
x.add('Z_IMP','6')
x.add('RZTEST_MODEL','0')

# Perform the parsing
x.read_input('input.profiles_3d')

x.printmsg()

sys.exit(x.error)
