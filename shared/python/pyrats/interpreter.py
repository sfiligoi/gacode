"""interpreter.py contains command line interpretation to feed into data.py.

    Contents:
        

"""

import sys
import os
import matplotlib.pyplot as plt
from pyrats.data import TGYROData

fignum = 1
sim1 = TGYROData(sys.argv[1])
for arg in sys.argv:
    if arg == '-ps':
        fig = plt.figure(fignum)
        fignum = fignum + 1
        ax1 = fig.add_subplot(221)
        ax1.set_ylabel('ne (1/cm^3)')
        ax1.set_xlabel('r/a')
        ax1.set_title('Electron Density vs. Radius')
        line1, = ax1.plot(sim1.get_r(), sim1.get_ne())
        ax2 = fig.add_subplot(222)
        ax2.set_ylabel('ni (1/cm^3)')
        ax2.set_xlabel('r/a')
        ax2.set_title('Ion Density vs. Radius')
        line2, = ax2.plot(sim1.get_r(), sim1.get_ni())
        ax3 = fig.add_subplot(223)
        ax3.set_ylabel('te (keV)')
        ax3.set_xlabel('r/a')
        ax3.set_title('Electron Temperature vs. Radius')
        line3, = ax3.plot(sim1.get_r(), sim1.get_Te())
        ax4 = fig.add_subplot(224)
        ax4.set_ylabel('ti (keV)')
        ax4.set_xlabel('r/a')
        ax4.set_title('Ion Temperature vs. Radius')
        line4, = ax4.plot(sim1.get_r(), sim1.get_Ti())
        plt.show()
    elif arg == '-gs':
        fig = plt.figure(fignum)
        fignum = fignum + 1
    elif arg == '-cs':
        #plot stuff
    elif arg == '-fs':
        #plot stuff
    elif arg = '-ms':
        #plot stuff
