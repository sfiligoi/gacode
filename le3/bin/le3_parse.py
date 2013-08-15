from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# LE3 input parameters
x.add('NTHETA','11')
x.add('NPHI','11')
x.add('NTHETAS','4')
x.add('NPHIS','4')
x.add('NTHETAP','64')
x.add('NPHIP','32')
x.add('RMIN','0.1')
x.add('RMAJ','1.0')
x.add('SHIFT','0.0')
x.add('KAPPA','1.0')
x.add('S_KAPPA','0.0')
x.add('DELTA','0.0')
x.add('S_DELTA','0.0')
x.add('ZETA','0.0')
x.add('S_ZETA','0.0')
x.add('ZMAG','0.0')
x.add('DZMAG','0.0')
x.add('HMIN','0.0')
x.add('Q','1.842')
x.add('M0','2')
x.add('N0','2')
x.add('TOL','1e-10')

# Perform the parsing
x.read_input('input.le3')

x.printmsg()        

sys.exit(x.error)


