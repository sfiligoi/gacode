import numpy as np
from scipy import interpolate, integrate, ndimage
import matplotlib.path as mp
import matplotlib._contour as _contour
import matplotlib.pyplot as plt

def contourPaths(x, y, Z, levels, remove_boundary_points=False, smooth_factor=1):
    '''
    :param x: 1D x coordinate
    :param y: 1D y coordinate
    :param Z: 2D data
    :param levels: levels to trace
    :param remove_boundary_points: remove traces at the boundary
    :param smooth_factor: smooth contours by cranking up grid resolution
    :return: list of segments
    '''

    sf = int(round(smooth_factor))
    if sf > 1:
        x = ndimage.zoom(x,sf)
        y = ndimage.zoom(y,sf)
        Z = ndimage.zoom(Z,sf)

    [X,Y]=np.meshgrid(x,y)
    contour_generator = _contour.QuadContourGenerator(X, Y, Z, None, True, 0)

    mx=min(x)
    Mx=max(x)
    my=min(y)
    My=max(y)

    allsegs = []
    for level in levels:
        segs=contour_generator.create_contour(level)
        if not remove_boundary_points:
            segs_ = segs
        else:
            segs_ = []
            for segarray in segs:
                x_ = segarray[:,0]
                y_ = segarray[:,1]
                valid = []
                for i in range(len(x_)-1):
                    if np.isclose(x_[i],x_[i+1]) and (np.isclose(x_[i],Mx) or np.isclose(x_[i],mx)):
                        continue
                    if np.isclose(y_[i],y_[i+1]) and (np.isclose(y_[i],My) or np.isclose(y_[i],my)):
                        continue
                    valid.append((x_[i],y_[i]))
                    if i==len(x_):
                        valid.append(x_[i+1],y_[i+1])
                if len(valid):
                    segs_.append(np.array(valid))

        segs=list(map(mp.Path,segs_))
        allsegs.append(segs)

    return allsegs

def paraboloid(x, y, z):
    '''
    z = ax*x^2 + bx*x + ay*y^2 + by*y + c

    NOTE: This function uses only the first 5 points of the x, y, z arrays
    to evaluate the paraboloid coefficients

    :return: ax,bx,ay,by,c
    '''
    if np.any(np.isnan(x.flatten())) or np.any(np.isnan(y.flatten())) or np.any(np.isnan(z.flatten())):
        raise(np.linalg.LinAlgError('paraboloid could not be fit with x=%s y=%s z=%s'%(x,y,z)))
    A=[]
    for k in range(5):
        A.append([x[k]**2,x[k],y[k]**2,y[k],1])
    A=np.array(A)
    ax,bx,ay,by,c=np.dot(np.linalg.inv(A),np.array(z[:5]))
    return ax,bx,ay,by,c

class RectBivariateSplineNaN:
    def __init__(self, Z, R, Q, *args, **kw):
        tmp=Q.copy()
        nan=np.isnan(tmp)
        self.thereAreNaN=False
        if nan.any():
            self.thereAreNaN=True
            tmp[nan]=0
            self.mask=interpolate.RectBivariateSpline(Z, R, nan ,kx=1, ky=1)
        self.spline=interpolate.RectBivariateSpline(Z, R, tmp, *args, **kw)

    def __call__(self,*args):
        tmp=self.spline(*args)
        if self.thereAreNaN:
            mask=self.mask(*args)
            tmp[mask>0.01]=np.nan
        return tmp

    def ev(self,*args):
        tmp=self.spline.ev(*args)
        if self.thereAreNaN:
            mask=self.mask.ev(*args)
            tmp[mask>0.01]=np.nan
        return tmp

def prgen_contour(geqdsk,nrz,levels,psinorm,narc,quiet):

    Rin      = np.linspace(0,geqdsk['RDIM'],geqdsk['NW'])+geqdsk['RLEFT']
    Zin      = np.linspace(0,geqdsk['ZDIM'],geqdsk['NH'])-geqdsk['ZDIM']/2.0+geqdsk['ZMID']
    psirz    = geqdsk['PSIRZ']
    efitpsi0 = geqdsk['SIMAG']
    efitpsi1 = geqdsk['SIBRY']
    efitp    = geqdsk['PRES']
    efitq    = geqdsk['QPSI']
    efitf    = geqdsk['FPOL']
    slfR0    = geqdsk['RMAXIS']
    slfZ0    = geqdsk['ZMAXIS']
    slfrlim  = geqdsk['RLIM']
    slfzlim  = geqdsk['ZLIM']

    if any(np.isnan(slfrlim)) or any(np.isnan(slfzlim)):
        print('ERROR: (prgen_contour) rlim/zlim arrays contain NaNs')
        return

    #-----------------------------------------------------------------
    # Change resolution
    #-----------------------------------------------------------------
    if not quiet:
        print('INFO: (prgen_contour) Levels based on psi')

    nrzz = int(nrz*(min(Zin)-max(Zin))/(min(Rin)-max(Rin)))
    
    slfR   = np.linspace(min(Rin),max(Rin),nrz)
    slfZ   = np.linspace(min(Zin),max(Zin),nrzz)
    slfPSI = RectBivariateSplineNaN(Zin,Rin,psirz)(slfZ,slfR)

    if not quiet:
       dres=np.sqrt((slfR[1]-slfR[0])**2+(slfZ[1]-slfZ[0])**2)
       print('INFO: (prgen_contour) Grid diagonal resolution [m] = {:.5f}'.format(dres))

    #-----------------------------------------------------------------
    # Crop 
    #-----------------------------------------------------------------
    if len(slfrlim) and len(slfzlim):
        #if not quiet:
        #    print('INFO: (prgen_contour) Cropping tables')
        bbox=[min(slfrlim),max(slfrlim),min(slfzlim),max(slfzlim)]
        limits=[max([np.argmin(abs(slfR-bbox[0]))-1,0]),
                min([np.argmin(abs(slfR-bbox[1]))+1,len(slfR)-1]),
                max([np.argmin(abs(slfZ-bbox[2]))-1,0]),
                min([np.argmin(abs(slfZ-bbox[3]))+1,len(slfZ)-1])]
        slfPSI=slfPSI[limits[2]:limits[3],limits[0]:limits[1]]
        slfR=slfR[limits[0]:limits[1]]
        slfZ=slfZ[limits[2]:limits[3]]

    #-----------------------------------------------------------------
    # Find axis
    #-----------------------------------------------------------------
    if not quiet:
        print('INFO: (prgen_contour) Finding magnetic axis and separatrix')

    Raxis = slfR0 ; Zaxis = slfZ0
    dmax = (slfR[-1] - slfR[0]) / 2.0 / 1.5

    RR,ZZ = np.meshgrid(slfR-Raxis,slfZ-Zaxis)
    DD = np.sqrt(RR**2+ZZ**2) < dmax
    tmp = slfPSI.copy()
    tmp[np.where(DD == 0)] = np.nan

    # figure out sign
    ax = 0
    for k in range(1,int(min([len(slfR)//2,len(slfZ)//2])))[::-1]:
        ri = (len(slfR)//2 + np.array([-k, 0, +k, 0, 0])).astype(int)
        zi = (len(slfZ)//2 + np.array([0, 0, 0, +k, -k])).astype(int)
        try:
            ax, bx, ay, by, c = paraboloid(slfR[ri], slfZ[zi], tmp[zi, ri])
            break
        except np.linalg.LinAlgError:
            pass
    if ax > 0:
        # look for the minimum
        m = np.nanargmin(tmp)
    else:
        # look for the maximum
        m = np.nanargmax(tmp)
        
    Zmi = int(m / slfPSI.shape[1])
    Rmi = int(m - Zmi * slfPSI.shape[1])

    # pick center points based on the grid
    slfR0 = slfR[Rmi]
    slfZ0 = slfZ[Zmi]

    # fit paraboloid in the vicinity of the grid-based center
    ri = (Rmi + np.array([-1, 0, +1, 0, 0])).astype(int)
    zi = (Zmi + np.array([0, 0, 0, +1, -1])).astype(int)
    ax, bx, ay, by, c = paraboloid(slfR[ri], slfZ[zi], slfPSI[zi,ri])
    slfR0_interp = -bx/(2*ax)
    slfZ0_interp = -by/(2*ay)

    # set as the center value of PSI (based on R0 and Z0)
    psi0 = RectBivariateSplineNaN(slfZ, slfR, slfPSI).ev(slfZ0_interp,slfR0_interp)
   
    #-----------------------------------------------------------------
    # Find separatrix
    #-----------------------------------------------------------------
    accuracy = 9
        
    # separatrix is found by looking for the largest closed path enclosing the magnetic axis

    flxm=np.nanmin(slfPSI)
    flxM=np.nanmax(slfPSI)
        
    kdbgmax=100
    psi1 = None
    for kdbg in range(kdbgmax):

        flx=0.5*(flxM+flxm)
    
        line=[]
        paths=contourPaths(slfR,slfZ,slfPSI,[flx])[0]
        for path in paths:
            if not np.isnan(path.vertices[:]).any() and np.allclose(path.vertices[0,0],path.vertices[-1,0]) and np.allclose(path.vertices[0,1],path.vertices[-1,1]):
                path.vertices[0,0] = path.vertices[-1,0] = (path.vertices[0,0]+path.vertices[-1,0])*0.5
                path.vertices[0,1] = path.vertices[-1,1] = (path.vertices[0,1]+path.vertices[-1,1])*0.5
                simplePath=mp.Path(path.vertices[::len(path.vertices[:,0])//10+1,:])
                if np.max(simplePath.vertices[:,0]) > slfR0 and np.min(simplePath.vertices[:,0]) < slfR0 and np.max(simplePath.vertices[:,1]) > slfZ0 and np.min(simplePath.vertices[:,1]) < slfZ0:
                    if simplePath.contains_point((slfR0,slfZ0)):
                        dR = path.vertices[1,0]-path.vertices[0,0]
                        dZ = path.vertices[1,1]-path.vertices[0,1]
                        orientation = int(np.sign((path.vertices[0,1]-slfZ0)*dR-(path.vertices[0,0]-slfR0)*dZ))
                        line=line=path.vertices[::orientation,:]
                        break
            
        if len(line):
            try:
                # stop condition
                np.testing.assert_array_almost_equal(sep/slfR0,line/slfR0,accuracy)
                break
            except Exception:
                pass
            finally:
                sep = line
                psi1 = flx
                flxm = flx
        else:
            open_sep=paths
            flxM=flx

    if kdbg==kdbgmax-1:
        print('Finding of last closed flux surface aborted after %d iterations!'%kdbgmax)

    print('INFO  (prgen_contour) dpsi = {:.9f} [EFIT] {:.9f} [new]'.format(efitpsi1-efitpsi0,psi1-psi0))

    # Separatrix hits edges of computation domain
    if ((np.abs(np.min(sep[:,0])-np.min(slfR)) < 1e-3) or
        (np.abs(np.max(sep[:,0])-np.max(slfR)) < 1e-3) or
        (np.abs(np.min(sep[:,1])-np.min(slfZ)) < 1e-3) or
        (np.abs(np.max(sep[:,1])-np.max(slfZ)) < 1e-3)):
        print("WARNING: (prgen_contour) New separatrix hits computation boundary")
            
    #-----------------------------------------------------------
    # Find surfaces
    #-----------------------------------------------------------
        
    if not quiet:
        print('INFO: (prgen_contour) Contour levels = {:d} | n_arc = {:d} | nrz = {:d}'.format(levels,narc,nrz))

    # absolute psi levels to interpolate (psinorm is normalized)
    # packsep < 1.0 increases separatrix packing
    packsep = 0.8
    
    out_psi = (np.linspace(0,psinorm,levels))**packsep*(psi1-psi0)+psi0

    CS = contourPaths(slfR,slfZ,slfPSI,out_psi)

    RI = np.zeros([narc,levels])
    ZI = np.zeros([narc,levels])
        
    for k,item1 in enumerate(CS):
        if k==0:
            # axis
            RI[:,k] = slfR0*np.ones(narc) ; ZI[:,k] = slfZ0*np.ones(narc)
        else:
            # all others
            path=item1[-1]
            # Reverse order (compared to original OMFIT order)
            r=path.vertices[::-1,0]
            z=path.vertices[::-1,1]
            if any(np.isnan(r*z)):
                print('ERROR: (prgen_contour) NaN encountered')
            r[0]=r[-1]=(r[0]+r[-1])*0.5
            z[0]=z[-1]=(z[0]+z[-1])*0.5
            n0 = len(r)
            # Arc length
            larc = np.zeros([n0])
            for i in range(n0-1):
               larc[i+1] = larc[i]+np.sqrt((r[i+1]-r[i])**2+(z[i+1]-z[i])**2)
               
            # Cubic interpolation from fine t0-mesh to coarse t-mesh
            t =np.linspace(0,1,narc)*larc[-1]
            cs=interpolate.CubicSpline(larc,r,bc_type='periodic') ; RI[:,k]=cs(t) 
            cs=interpolate.CubicSpline(larc,z,bc_type='periodic') ; ZI[:,k]=cs(t)

    efitpsi = np.linspace(out_psi[0],out_psi[-1],len(efitp))
    cs = interpolate.interp1d(efitpsi,efitp,kind='quadratic') ; out_p = cs(out_psi)
    cs = interpolate.interp1d(efitpsi,efitq,kind='quadratic') ; out_q = cs(out_psi)
    cs = interpolate.interp1d(efitpsi,efitf,kind='quadratic') ; out_f = cs(out_psi)

    # Recalculate q based on definition (and some identities)
    loopint = np.zeros([levels])
    for k in range(levels):
       for i in range(narc-1):
          loopint[k] = loopint[k]+(RI[i+1,k]-RI[i,k])*(ZI[i+1,k]+ZI[i,k])/(RI[i+1,k]+RI[i,k])
                    
    loopint = loopint/(2*np.pi)
    new_q = out_f*np.gradient(loopint,out_psi)

    return RI,ZI,out_psi,new_q,out_p,out_f

