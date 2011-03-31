#!/bin/env python

import pylab, numpy


def plotVars(rGrid,gridVarName,data,dVars,dataLabel,dataTitle,dMarker):
   pylab.xlabel(gridVarName)
   pylab.ylabel(dataLabel)
   for iv in range(data.shape[0]):
     pylab.plot(rGrid,data[iv,:],label=dVars[iv],markevery=2,ms=10.0,marker=dMarker[iv])
   pylab.legend()
   pylab.title(dataTitle)
   #pylab.show()


echarge=numpy.float(1.602e-19)
bt_exp=numpy.float(2.011546195E+00)
b_norm=bt_exp

possibleMarkers=['+','^','*','<','>','d',',','o','H','1','2','3','4','D','_','.','h','p','s','v','x','|']
equilVars=["rho","Er_norm","q","rho_star","rmaj","omega_norm","shear_norm","n_norm","T_norm","v_norm/a"]
equilSpeciesVars=["nnorm","tnorm","a/Ln","a/LT","nu"]
for var in equilSpeciesVars: equilVars.append(var+"_i")
for var in equilSpeciesVars: equilVars.append(var+"_e")

# Load the data.  Note that pylab loads it by row, when we want it by
# column to make it easier to loop over the names
equilData=pylab.loadtxt("equil.out").transpose()
nradial=equilData.shape[1]

#Most commonly used
n_norm=equilData[7]*1.e19
t_norm=equilData[8]
v_norm=equilData[9]
#a_norm=transportExpData[0,0]/transportData[0,0] #a_norm
print "n_norm = ",n_norm
print "T_norm = ",t_norm
print "v_norm = ",v_norm

##========================================================================
#   Plotting
##========================================================================
if __name__ == "__main__":
   #tecFileWrite("jbs.dat",jbsVars,jbs,gridVars,grid)
   #tecFileWrite("efluxe.dat",efluxeVars,efluxe,gridVars,grid)
   #tecFileWrite("norm.dat",normVars,norm,gridVars,grid)
   # Choose which type of grid: ["rho","r/a","r","Rout"]
   pylab.subplot(221)
   dTitle="Densities"
   plotData= numpy.zeros((2,nradial), numpy.float)
   plotData[0]=equilData[10]*n_norm
   plotData[1]=equilData[15]*n_norm
   dMarkers=possibleMarkers[0:plotData.shape[0]]
   pVars=["n_i","n_e"]
   plotVars(equilData[0,:],"r/a",plotData,pVars,"Density",dTitle,dMarkers)

   pylab.subplot(222)
   dTitle="Temperature"
   plotData[0]=equilData[11]*t_norm
   plotData[1]=equilData[16]*t_norm
   dMarkers=possibleMarkers[0:plotData.shape[0]]
   pVars=["T_i","T_e"]
   plotVars(equilData[0,:],"r/a",plotData,pVars,"Temperatures",dTitle,dMarkers)

   pylab.subplot(223)
   dTitle="Density Scale Lengths"
   plotData[0]=equilData[12]
   plotData[1]=equilData[17]
   dMarkers=possibleMarkers[0:plotData.shape[0]]
   pVars=["a/Ln_i","a/Ln_e"]
   plotVars(equilData[0,:],"r/a",plotData,pVars,"a/Ln",dTitle,dMarkers)

   pylab.subplot(224)
   dTitle="Temperature Scale Lengths"
   plotData[0]=equilData[13]
   plotData[1]=equilData[18]
   dMarkers=possibleMarkers[0:plotData.shape[0]]
   pVars=["a/LT_i","a/LT_e"]
   plotVars(equilData[0,:],"r/a",plotData,pVars,"a/LT",dTitle,dMarkers)

   pylab.show()
