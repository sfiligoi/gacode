import sys
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
from scipy import interpolate

def prgen_contour(g,mult,nc,psinorm,narc,quiet):

   #-------------------------------------------------------------------------
   # Extract data from gfile
   nx    = g['NW']
   ny    = g['NH']
   psi0  = g['SIMAG']
   psi1  = g['SIBRY']
   efitp = g['PRES']
   efitq = g['QPSI']
   efitf = g['FPOL']
   lx    = g['RDIM']
   ly    = g['ZDIM']

   npsi = len(efitp)
   psi_efit = g['PSIRZ']
   #-------------------------------------------------------------------------

   #-------------------------------------------------------------------------
   # Pixel representation
   xv = np.linspace(0,1,nx)
   yv = np.linspace(0,1,ny)

   # spline interpolation (quadratic, not cubic, to avoid oscillation)
   spl = interpolate.RectBivariateSpline(yv,xv,psi_efit,kx=2,ky=2)

   if mult > 1:
      mx = mult*nx
      my = mult*ny

      xv = np.linspace(0,1,mx)
      yv = np.linspace(0,1,my)

      psi_efit = spl(yv,xv)

   else:
      mx = nx
      my = ny

   # Pixel values of magnetic axis
   kx = int((g['RMAXIS']-g['RLEFT'])/g['RDIM']*(mx-1))
   ky = int((g['ZMAXIS']-g['ZMID']+g['ZDIM']/2)/g['ZDIM']*(my-1))
   #-------------------------------------------------------------------------

   #-------------------------------------------------------------------------
   # Find separatrix using bisection
   dpsi = psi1-psi0
   z0 = psi1-0.5*dpsi
   dz = dpsi/3
   tol = 1e-14
   err = 1.0
   while dz > tol:
      z0 = z0+dz
      contours = measure.find_contours(psi_efit,z0)
      found = 0
      for contour in contours:
         if contour[0,0] == contour[-1,0] and contour[0,1] == contour[-1,1]:
            xc0 = contour[:,1] ; yc0 = contour[:,0]
            if min(xc0) < kx < max(xc0) and min(yc0) < ky < max(yc0):
               found = 1
               xc = xc0 ; yc = yc0
      if found == 0:
         z0 = z0-dz
         dz = dz/2

   psi_sep = z0
   #-------------------------------------------------------------------------

   psi_efit_min = np.min(psi_efit[yaxis-8:yaxis+8,xaxis-8:xaxis+8])
   print('INFO: (prgen_contour) psi_min (% efit) = {:.3f}'.format(100*psi_efit_min/psi0))
   print('INFO: (prgen_contour) psi_sep (% efit) = {:.3f}'.format(100*psi_sep/psi1))

   #-------------------------------------------------------------------------
   # Compute interior contours
   #
   rv = np.zeros([narc,nc])
   zv = np.zeros([narc,nc])
   psic = np.zeros(nc)

   # Move inside of separatrix by amount psinorm
   dpsi = psinorm*(psi1-psi0)/nc

   for i in range(nc):
      z0 = psi0+(i+1.0)*dpsi
      psic[i] = z0
      contours = measure.find_contours(psi_efit,z0)
      for contour in contours:
         if contour[0,0] == contour[-1,0] and contour[0,1] == contour[-1,1]:
            xc = contour[:,1] ; yc = contour[:,0]
            if min(xc) < xaxis < max(xc) and min(yc) < yaxis < max(yc):
               # Arc length
               dl = np.sqrt(np.diff(xc)**2+np.diff(yc)**2)
               larc = np.zeros([len(xc)]) ; larc[1:] = np.cumsum(dl)

               # Cubic interpolation from fine contour mesh to coarse t-mesh
               t = np.linspace(0,1,narc)*larc[-1]
               cs = interpolate.splrep(larc,xc,per=True) ; r=interpolate.splev(t,cs) 
               cs = interpolate.splrep(larc,yc,per=True) ; z=interpolate.splev(t,cs)

               # Shift elements so that first index at max(R).
               s = np.argmax(r)
               r[0:-1] = np.roll(r[0:-1],-s) ; r[-1] = r[0] 
               z[0:-1] = np.roll(z[0:-1],-s) ; z[-1] = z[0]

               # Map back from pixels to physical dimensions
               rv[:,i] = g['RLEFT']+g['RDIM']/(mx-1)*r
               zv[:,i] = g['ZMID']-g['ZDIM']/2+g['ZDIM']/(my-1)*z
   #-------------------------------------------------------------------------

   #-------------------------------------------------------------------------
   # Radial profiles of p,f,q
   #
   efitpsi = np.linspace(psi0,psi1,npsi)
   cs = interpolate.interp1d(efitpsi,efitp,kind='quadratic') ; out_p = cs(psic)
   cs = interpolate.interp1d(efitpsi,efitf,kind='quadratic') ; out_f = cs(psic)
   cs = interpolate.interp1d(efitpsi,efitq,kind='quadratic') ; out_qe = cs(psic)

   # Recalculate q based on definition (and some identities)
   loopint = np.zeros(nc)
   for i in range(narc-1):
      loopint[:] = loopint[:]+(rv[i+1,:]-rv[i,:])*(zv[i+1,:]+zv[i,:])/(rv[i+1,:]+rv[i,:])

   cs = interpolate.splrep(psic,loopint) ; out_q = interpolate.splev(psic,cs,der=1)
   out_q = out_f*out_q/(2*np.pi)
   #-------------------------------------------------------------------------

   return rv,zv,psic,out_p,out_f,out_qe,out_q
