"""interpreter.py contains command line interpretation to feed into data.py.

    Contents:
        

"""

import sys
import os
import matplotlib.pyplot as plt
from pyrats.tgyro.data import TGYROData

fignum = 1
sim1 = TGYROData(sys.argv[1])
n1=sys.argv[2]
n2=sys.argv[3]
verbose=bool(int(sys.argv[4]))
legend=bool(int(sys.argv[5]))

for arg in sys.argv:
    if arg == 'ps':
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
    elif arg == 'gs':
        fig = plt.figure(fignum)
        fignum = fignum + 1
        ax1 = fig.add_subplot(221)
        ax1.set_ylabel('a/Lni')
        ax1.set_xlabel('r/a')
        ax1.set_title('a/Lni vs. Radius')
        ax1.plot(sim1.get_r(), sim1.get_gradient_factor('a/Lni'))
        ax2 = fig.add_subplot(222)
        ax2.set_ylabel('a/Lne')
        ax2.set_xlabel('r/a')
        ax2.set_title('a/Lne vs. Radius')
        ax2.plot(sim1.get_r(), sim1.get_gradient_factor('a/Lne'))
        ax3 = fig.add_subplot(223)
        ax3.set_ylabel('a/LTi')
        ax3.set_xlabel('r/a')
        ax3.set_title('a/LTi vs. Radius')
        ax3.plot(sim1.get_r(), sim1.get_gradient_factor('a/LTi'))
        ax4 = fig.add_subplot(224)
        ax4.set_ylabel('a/LTe')
        ax4.set_xlabel('r/a')
        ax4.set_title('a/LTe vs. Radius')
        ax4.plot(sim1.get_r(), sim1.get_gradient_factor('a/LTe'))
        plt.show()
    elif arg == 'cs':
        fig = plt.figure(fignum)
        fignum = fignum + 1
        ax1 = fig.add_subplot(221)
        ax1.set_ylabel('chi_e_turb')
        ax1.set_xlabel('r/a')
        ax1.set_title('chi_e_turb vs. Radius')
        ax1.plot(sim1.get_r(), sim1.get_chi_e_turb())
        ax2 = fig.add_subplot(222)
        ax2.set_ylabel('chi_e_neo')
        ax2.set_xlabel('r/a')
        ax2.set_title('chi_e_neo vs. Radius')
        ax2.plot(sim1.get_r(), sim1.get_chi_e_neo())
        ax3 = fig.add_subplot(223)
        ax3.set_ylabel('chi_i_turb')
        ax3.set_xlabel('r/a')
        ax3.set_title('chi_i_turb vs. Radius')
        ax3.plot(sim1.get_r(), sim1.get_chi_i_turb())
        ax4 = fig.add_subplot(224)
        ax4.set_ylabel('chi_i_neo')
        ax4.set_xlabel('r/a')
        ax4.set_title('chi_i_neo vs. Radius')
        ax4.plot(sim1.get_r(), sim1.get_chi_i_neo())
        plt.show()
    elif arg == 'fs':
        fig = plt.figure(fignum)
        fignum = fignum + 1
        ax1 = fig.add_subplot(221)
        ax1.set_ylabel('flux_e_turb')
        ax1.set_xlabel('r/a')
        ax1.set_title('flux_e_turb vs. Radius')
        ax1.plot(sim1.get_r(), sim1.get_flux_e_turb())
        ax2 = fig.add_subplot(222)
        ax2.set_ylabel('flux_e_neo')
        ax2.set_xlabel('r/a')
        ax2.set_title('flux_e_neo vs. Radius')
        ax2.plot(sim1.get_r(), sim1.get_flux_e_neo())
        ax3 = fig.add_subplot(223)
        ax3.set_ylabel('flux_i_turb')
        ax3.set_xlabel('r/a')
        ax3.set_title('flux_i_turb vs. Radius')
        ax3.plot(sim1.get_r(), sim1.get_flux_i_turb())
        ax4 = fig.add_subplot(224)
        ax4.set_ylabel('flux_i_neo')
        ax4.set_xlabel('r/a')
        ax4.set_title('flux_i_neo vs. Radius')
        ax4.plot(sim1.get_r(), sim1.get_flux_i_neo())
        plt.show() 
    elif arg == 'ms':
        fig = plt.figure(fignum)
        fignum = fignum + 1
        ax1 = fig.add_subplot(221)
        ax1.set_ylabel('mflux_e_turb')
        ax1.set_xlabel('r/a')
        ax1.set_title('mflux_e_turb vs. Radius')
        ax1.plot(sim1.get_r(), sim1.get_mflux_e_turb())
        ax2 = fig.add_subplot(222)
        ax2.set_ylabel('mflux_e_neo')
        ax2.set_xlabel('r/a')
        ax2.set_title('mflux_e_neo vs. Radius')
        ax2.plot(sim1.get_r(), sim1.get_mflux_e_neo())
        ax3 = fig.add_subplot(223)
        ax3.set_ylabel('mflux_i_turb')
        ax3.set_xlabel('r/a')
        ax3.set_title('mflux_i_turb vs. Radius')
        ax3.plot(sim1.get_r(), sim1.get_mflux_i_turb())
        ax4 = fig.add_subplot(224)
        ax4.set_ylabel('mflux_i_neo')
        ax4.set_xlabel('r/a')
        ax4.set_title('mflux_i_neo vs. Radius')
        ax4.plot(sim1.get_r(), sim1.get_mflux_i_neo())
        plt.show()

