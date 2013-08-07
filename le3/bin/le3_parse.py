from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# LE3 input parameters
x.add('NTHETA','21')
x.add('NPHI','21')
x.add('NTHETAS','40')
x.add('NPHIS','40')
x.add('RMIN','0.5')
x.add('RMAJ','3.0')
x.add('HMIN','0.05')
x.add('Q','1.842')
x.add('M0','2')
x.add('N0','4')
x.add('TOL','1e-10')
x.add('RESTART_FLAG','0')
x.add('SOLVE_METHOD','1')

# Perform the parsing
x.read_input('input.le3')

x.printmsg()        

sys.exit(x.error)


