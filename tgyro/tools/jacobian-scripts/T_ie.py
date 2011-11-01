#!/usr/bin/python
#try to color some plots
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.text as text
import matplotlib.mlab as mlab
import matplotlib.cbook as cbook

def plot_T_ie():
    prop = matplotlib.font_manager.FontProperties(size=14)

    fig = plt.figure(figsize=(6,4))
    plot = fig.add_subplot(111)
    
    data = np.genfromtxt('output/profiles.out',usecols=(0,1,2))
    radius = data[:,0]
    ti = data[:,1]
    te = data[:,2]

    color = []

    plot.set_xlabel(r'r/a')
    plot.set_ylabel(r'$T_e$, $T_i$ (kev)')
    plt.title('Electron, ion temperature')
    plot.grid(True)

    plot.plot(radius,te,'-o',color='darkgreen', lw=2.5, label="Electron")
    plot.plot(radius,ti,'-o', color='red', lw=2.5, label="Ion")
    handles, labels=plot.get_legend_handles_labels()
    plot.legend(labels)

#plt.show()

    plt.savefig('plots/T_ie.pdf')
#a=np.array([1,2,3,4])


if __name__=="__main__":
    print('Plotting T_i, T_e')
    plot_T_ie()
