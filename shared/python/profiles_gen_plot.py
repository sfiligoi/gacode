import sys
import matplotlib.pyplot as plt
from profiles_genData import profiles_genData

fignum = 1
prof1 = profiles_genData(sys.argv[1])
if len(sys.argv) > 2:
    if sys.argv[2] == '-options':
        keys = []
        for k, v in prof1.data.iteritems():
            s = 0
            for x in range(len(v)):
                s = s + float(v[x])
            if s != 0:
                lflag = 0
                for letter in k:
                    if letter == '(':
                        lflag = 1
                    if letter == '#':
                        lflag = 1
                if lflag == 0:
                    keys.append(k)
        keys.sort()
        for k in keys:
            print k
    else:
        args = sys.argv[2:]
        for arg in args:
            fig = plt.figure(fignum)
            fignum = fignum + 1
            ax = fig.add_subplot(111)
            for k in prof1.data.iterkeys():
                if (arg + ' ') == k[:len(arg) + 1]:
                    ax.set_ylabel(k)
            ax.set_xlabel('rho (-)')
            ax.set_title(arg + ' vs. rho')
            line, = ax.plot(prof1.get('rho'), prof1.get(arg), 'k')
        plt.show()
if len(sys.argv) <= 2:
    print "ERROR: Please specify data to be plotted, or ask for"
    print "options."
from profiles_genData import profiles_genData
prof1 = profiles_genData()
