"""Generate a plot of a profile function"""

import sys
import string
import numpy as np
from gacodeplotdefs import *
from profiles_gen.data import profiles_genData

#---------------------------------------------------------------
def fancytag(tag):
 
    if tag == 'bunit':
        fancy = r'$B_\mathrm{unit}$'
    elif tag == 'rmin':
        fancy = r'$r$'
    elif tag == 'rmaj':
        fancy = r'$R_0$'
    elif tag == 'drmaj':
        fancy = r'$dR_0/dr$'
    elif tag == 'kappa':
        fancy = r'$\kappa$'
    elif tag == 'skappa':
        fancy = r'$r \, d\kappa/dr$'
    elif tag == 'delta':
        fancy = r'$\delta$'
    elif tag == 'sdelta':
        fancy = r'$r \, d\delta/dr$'
    elif tag == 'Te':
        fancy = r'$T_\mathrm{e}$'
    elif tag == 'ne':
        fancy = r'$n_\mathrm{e}$'
    elif tag == 'Ti_1':
        fancy = r'$T_\mathrm{i1}$'
    elif tag == 'Ti_2':
        fancy = r'$T_\mathrm{i2}$'
    elif tag == 'Ti_3':
        fancy = r'$T_\mathrm{i3}$'
    elif tag == 'Ti_4':
        fancy = r'$T_\mathrm{i4}$'
    elif tag == 'Ti_5':
        fancy = r'$T_\mathrm{i5}$'
    elif tag == 'ni_1':
        fancy = r'$n_\mathrm{i1}$'
    elif tag == 'ni_2':
        fancy = r'$n_\mathrm{i2}$'
    elif tag == 'ni_3':
        fancy = r'$n_\mathrm{i3}$'
    elif tag == 'ni_4':
        fancy = r'$n_\mathrm{i4}$'
    elif tag == 'ni_5':
        fancy = r'$n_\mathrm{i5}$'
    elif tag == 'omega0':
        fancy = r'$\omega_0$'
    elif tag == 'w0p':
        fancy = r'$d \omega_0/dr$'
    elif tag == 'jbs_neo':
        fancy = r'$\mathrm{jbs}_\mathrm{neo}$'
    elif tag == 'jbs_sauter':
        fancy = r'$\mathrm{jbs}_\mathrm{sauter}$'
    elif tag == 'jbs_nclass':
        fancy = r'$\mathrm{jbs}_\mathrm{nclass}$'
    elif tag == 'pow_e':
        fancy = r'$\mathrm{P}_\mathrm{e}\;\mathrm(MW)$'
    elif tag == 'pow_i':
        fancy = r'$\mathrm{P}_\mathrm{i}\;\mathrm(MW)$'
    elif tag == 'pow_e_aux':
        fancy = r'$\mathrm{P}_\mathrm{e,aux}\;\mathrm(MW)$'
    elif tag == 'pow_i_aux':
        fancy = r'$\mathrm{P}_\mathrm{i,aux}\;\mathrm(MW)$'
    elif tag == 'pow_ei':
        fancy = r'$\mathrm{P}_\mathrm{ei}\;\mathrm(MW)$'
    elif tag == 'pow_e_line':
        fancy = r'$\mathrm{P}_\mathrm{e,rad}\;\mathrm(MW)$'
    else:
        fancy = r'$\mathrm{'+tag+'}$'

    return fancy
#---------------------------------------------------------------

rvar    = sys.argv[1]
infiles = sys.argv[2]
plots   = sys.argv[3]
ftype   = sys.argv[4]
loc     = int(sys.argv[5])
t       = sys.argv[6]
title   = sys.argv[7]

plotvec = string.splitfields(plots,',')
filevec = string.splitfields(infiles,',')

n = len(plotvec)

for j in range(n):

    if len(t.split(',')) == n:
        tlabel = r'$\;\mathrm{'+t.split(',')[j]+'}$'
    else:
        tlabel = ''
    

    if filevec[j] == '.':
        filevec[j] = filevec[j-1]

    if plotvec[j] == '.':
        plotvec[j] = plotvec[j-1]

    prof = profiles_genData(filevec[j])
    tag  = plotvec[j]
    keys  = sorted(prof.data.keys())

    success = 0
    for i in range(len(keys)):
        if tag == keys[i].split()[0]:
            success = 1
            fulltag = keys[i]

    if success == 0:
        print "ERROR: (profiles_gen_plot) Bad profile = "+tag
        sys.exit()

    if j==0:
        fig = plt.figure(figsize=(12,7))
        ax  = fig.add_subplot(111)

        if title != 'null':
            ax.set_title(r'$'+title+'$',size=18)

        if rvar == "r":
            ax.set_xlabel(r"$r \, \mathrm{(m)}$",fontsize=GFONTSIZE)
        if rvar == "r/a":
            ax.set_xlabel(r"$r/a$",fontsize=GFONTSIZE)
        if rvar == "rho":
            ax.set_xlabel(r"$\rho$",fontsize=GFONTSIZE)

        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")

    if rvar == "r":
        x = prof.data['rmin']

    if rvar == "r/a":
        x = prof.data['rmin']
        x = x/max(x)

    if rvar == "rho":
        x = prof.data['rho']

    ftag = fancytag(tag)
    ax.plot(x,prof.data[fulltag],'o-',label=ftag+tlabel)

ax.legend(loc=loc)

if ftype == 'screen':
    plt.show()
else:
    print "Saving plot to "+ftype
    outfile = ftype
    plt.savefig(outfile)
