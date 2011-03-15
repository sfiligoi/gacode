#!/bin/env python

import pylab, numpy


def myPlotVars(rGrid,gridVarName,data,dVars,dataLabel,dataTitle,dMarker):
   pylab.xlabel(gridVarName)
   pylab.ylabel(dataLabel)
   for iv in range(data.shape[0]):
     pylab.plot(rGrid,data[iv,:],label=dVars[iv],markevery=2,ms=10.0,marker=dMarker[iv])
   pylab.legend()
   pylab.title(dataTitle)
   #pylab.show()

# This writes the output file into a form that Eric likes
def tecFileWrite(fname,dataVars,varData,gridVars,grid):
      datFile=open(fname,"w")
      datFile.write("VARIABLES = ")
      for ivar in range(len(gridVars)):
         datFile.write(gridVars[ivar]+",")
      for ivar in range(len(dataVars)-1):
         datFile.write(dataVars[ivar]+",")
         datFile.write(dataVars[-1]+"\n")
      for iradial in range(varData.shape[1]):
         for ivar in range(len(gridVars)):
            datFile.write('%g ' % (grid[ivar][iradial]))
         for ivar in range(len(dataVars)):
            datFile.write('%g ' % (varData[ivar][iradial]))
         datFile.write("\n")
      datFile.close()


echarge=numpy.float(1.602e-19)
bt_exp=numpy.float(2.011546195E+00)
b_norm=bt_exp

possibleMarkers=['+','^','*','<','>','d',',','o','H','1','2','3','4','D','_','.','h','p','s','v','x','|']
theoryVars=["rho","nflux_hh","efluxi_hh","efluxe_hh","Jbs_hh","ki_hh","Vpar_hh","Vpol_hh","efluxi_ch","efluxi_chmod","efluxi_tg","efluxi_tgmod","Jbs_S","ki_S","Vpar_S","Vpol_S","ePhi_HR","nfluxi_hs","efluxi_hs","nfluxe_hs","efluxe_hs","nfluxi_hsmod","efluxi_hsmod","nfluxe_hsmod","efluxe_hsmod"]
equilVars=["rho","Er_norm","q","rho_star","rmaj","omega_norm","shear_norm","n_norm","T_norm","v_norm/a","norm_density","norm_temp","norm_dens_scale","norm_temp_scale","norm_collision_freq"]
transportVars=["r/a","phi_1","J_BS","V_phi","Vpar_0","nfluxi_neo","efluxi_neo","momrfluxi_neo","vpari_neo","k1i_neo","K1i_neo","Vpoli_neo","Vtori_neo","nfluxe_neo","efluxe_neo","momrfluxe_neo","vpare_neo","k1e_neo","K1e_neo","Vpole_neo","Vtore_neo"]
transportExpVars=transportVars
transportExpVars[0]="r"   # Not exactly the same

# Load the data.  Note that pylab loads it by row, when we want it by
# column to make it easier to loop over the names
theoryData=pylab.loadtxt("theory.out").transpose()
equilData=pylab.loadtxt("equil.out").transpose()
# Transport data does it by order so it is a bit different
tempTransportExpData=pylab.loadtxt("transport_exp.out").transpose()
tempTransportData=pylab.loadtxt("transport.out").transpose()
norder=tempTransportData.shape[1]/theoryData.shape[1]
nradial=tempTransportData.shape[1]/norder
nvars=tempTransportData.shape[0]
#transportData=numpy.zeros((nvars,norder,nradial), numpy.float)
transportData=tempTransportData.reshape(nvars,norder,nradial)
transportExpData=tempTransportExpData.reshape(nvars,norder,nradial)
# Sum over the orders to allow plotting of best solution
transportDataSum=transportData.sum(1)
transportExpDataSum=transportExpData.sum(1)


##
#   Reorganize the data to make it easier to manipulate
###
## 
# 
gridVars=["rho","r/a","r","Rout"]
grid = numpy.zeros((len(gridVars),theoryData.shape[1]), numpy.float)
grid[0]=theoryData[0]
grid[1]=transportData[0,0]
grid[2]=transportExpData[0,0]
grid[3]=equilData[4]  # Renormalize down below

###
##  Normalization
# 
n_norm=equilData[7]*1.e19  #n_norm
t_norm=equilData[8]  #t_norm
v_norm=equilData[9]  # vth_norm/a
a_norm=transportExpData[0,0]/transportData[0,0] #a_norm

# Also renormalize R_maj
grid[3]=equilData[4]*a_norm + grid[2]


###
##  Electronstatic potential
# 
phiVars=["HR","neo1","neosum"]
phi_norm=1.
phiMKSfac=transportExpData[1,0,4]/transportData[1,0,4]
phi_mult=phiMKSfac
phi = numpy.zeros((len(phiVars),theoryData.shape[1]), numpy.float)
phi[0]=theoryData[16]*phi_mult
phi[1]=transportData[1,0]*phi_mult
phi[2]=transportDataSum[1]*phi_mult
#phi[2]=transportExpData[1,0]
#phi[2]=transportExpDataSum[1]

###
##  Bootstrap current  -- dimensionalize as we go along
# 
jbsVars=["hh","S","neo","neosum"]
# This puts it into A/m^2 T  which is kind of funny
jbs_norm=echarge*n_norm*v_norm
jbsMKSfac=transportExpData[2,0,4]/transportData[2,0,4]
jbs_mult=jbsMKSfac
#jbs[2]=transportExpData[2,0]
jbs = numpy.zeros((len(jbsVars),theoryData.shape[1]), numpy.float)
jbs[0]=theoryData[4]*jbs_mult 
jbs[1]=theoryData[12]*jbs_mult 
jbs[2]=transportData[2,0]*jbs_mult 
jbs[3]=transportDataSum[2]*jbs_mult 
#jbs[2]=transportExpData[2,0]
#jbs[3]=transportExpDataSum[2]

###
##  Electron energy flux -- dimensionalize as we go along
# 
efluxeVars=["hh","hs","hsmod","neo1","neosum"]
efluxe_norm=n_norm*v_norm*t_norm
efluxeMKSfac=transportExpData[14,0,4]/transportData[14,0,4]
efluxe_mult=efluxeMKSfac
efluxe = numpy.zeros((len(efluxeVars),theoryData.shape[1]), numpy.float)
efluxe[0]=theoryData[3]*efluxe_mult 
efluxe[1]=theoryData[20]*efluxe_mult 
efluxe[2]=theoryData[24]*efluxe_mult 
efluxe[3]=transportData[14,0]*efluxe_mult 
efluxe[4]=transportDataSum[14]*efluxe_mult 
#efluxe[3]=transportExpData[14,0]*efluxe_mult 
#efluxe[4]=transportExpDataSum[14]*efluxe_mult 

###
##  Ion energy flux -- dimensionalize as we go along
# 
efluxiVars=["hh","ch","chmod","tgch","tgchmod","hs","hsmod","neo","neosum"]
efluxi_norm=n_norm*v_norm*t_norm
efluxiMKSfac=transportExpData[6,0,4]/transportData[6,0,4]
efluxi_mult=efluxiMKSfac
efluxi = numpy.zeros((len(efluxiVars),theoryData.shape[1]), numpy.float)
efluxi[0]=theoryData[2]*efluxi_mult 
efluxi[1]=theoryData[8]*efluxi_mult 
efluxi[2]=theoryData[9]*efluxi_mult 
efluxi[3]=theoryData[10]*efluxi_mult 
efluxi[4]=theoryData[11]*efluxi_mult 
efluxi[5]=theoryData[18]*efluxi_mult 
efluxi[6]=theoryData[22]*efluxi_mult 
efluxi[7]=transportData[6,0]*efluxi_mult 
efluxi[8]=transportDataSum[6]*efluxi_mult 
#efluxi[7]=transportExpData[6,0]*efluxi_mult 
#efluxi[8]=transportExpDataSum[6]*efluxi_mult 



##========================================================================
#   Now write out to tecplot-like files
##========================================================================
if __name__ == "__main__":
   #tecFileWrite("jbs.dat",jbsVars,jbs,gridVars,grid)
   #tecFileWrite("efluxe.dat",efluxeVars,efluxe,gridVars,grid)
   #tecFileWrite("norm.dat",normVars,norm,gridVars,grid)


   # Choose which type of grid: ["rho","r/a","r","Rout"]
   igrid=0
   dTitle="Analytic calculations of bootstrap current"
   dMarkers=possibleMarkers[0:jbs.shape[0]]
   pylab.subplot(221)
   ylabel=numpy.str("J_BS (A/m^2)")
   myPlotVars(grid[igrid],gridVars[igrid],jbs,jbsVars,ylabel,dTitle,dMarkers)

   pylab.subplot(222)
   dTitle="<e Phi/T_norm>"
   dMarkers=possibleMarkers[0:phi.shape[0]]
   ylabel="<e Phi/T> (V^2)"
   myPlotVars(grid[igrid],gridVars[igrid],phi,phiVars,ylabel,dTitle,dMarkers)

   pylab.subplot(223)
   dTitle="Analytic calculations of electron energy flux"
   dMarkers=possibleMarkers[0:efluxe.shape[0]]
   ylabel="E_flux_e (e19 m-2 s-1)"
   myPlotVars(grid[igrid],gridVars[igrid],efluxe,efluxeVars,ylabel,dTitle,dMarkers)

   pylab.subplot(224)
   dTitle="Analytic calculations of ion energy flux"
   dMarkers=possibleMarkers[0:efluxi.shape[0]]
   ylabel="E_flux_i (e19 m-2 s-1)"
   myPlotVars(grid[igrid],gridVars[igrid],efluxi,efluxiVars,ylabel,dTitle,dMarkers)
   pylab.show()


