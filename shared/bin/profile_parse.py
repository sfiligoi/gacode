from gacodeinput import *
import sys

x = ProfileInput()

x.set_extension('.gen')

# input.profiles scalar parameters

x.add('N_EXP','41')
x.add('BT_EXP','-2.0')
x.add('ARHO_EXP','0.3')
 
# Perform the parsing
x.read_profile('input.profiles')

x.printmsg()        

sys.exit(x.error)

