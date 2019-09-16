import sys
import string
import numpy as np

def minmax(x,y):

    n = len(x)

    # High accuracy estimate of minimum x given initial guess as
    # minimum of array x.

    # Maximum x
    i = np.argmax(x)

    ip = np.mod(i+1,n-1)
    im = np.mod(i-1,n-1)
    x1 = x[im] 
    x2 = x[i]
    x3 = x[ip]
    y1 = y[im]
    y2 = y[i]
    y3 = y[ip]

    if np.abs(x1-x2) <= 1e-11 and np.abs(x2-x3) <= 1e-11:
        print('ERROR: (minmax) x1=x2=x3 error')
        sys.exit()

    d12 = y1-y2
    d13 = y1-y3
    d23 = y2-y3
    d3 = d13*d23
    d2 = -d12*d23
    d1 = d12*d13
    y0 = (y1+y2)*x3/d3+(y1+y3)*x2/d2+(y2+y3)*x1/d1
    y0 = y0/(2*x3/d3+2*x2/d2+2*x1/d1)
    x0 = (y0-y1)*(y0-y2)*x3/d3+(y0-y1)*(y0-y3)*x2/d2+(y0-y2)*(y0-y3)*x1/d1

    return x0,y0

def intersect(x,y,y0):

    n = len(x)

    x_c = 0.5*(np.min(x)+np.max(x))

    for i in range(n-1):

        ya = y[i]
        yb = y[i+1]
        xa = x[i]
        xb = x[i+1]

        if (ya-y0)*(yb-y0) < 0.0:

            # Linear approximation (unused)
            # x0 = (y0-ya)*xb/(yb-ya)+(y0-yb)*xa/(ya-yb)

            if np.abs(ya-y0) < np.abs(yb-y0):
                if i == 0:
                    i0 = n-2
                else:
                    i0 = i-1

                y1 = y[i0]
                y2 = ya
                y3 = yb
                x1 = x[i0]
                x2 = xa
                x3 = xb
            else:
                if i == n-2:
                    i0 = 1
                else:
                    i0 = i+2

                y1 = ya
                y2 = yb
                y3 = y[i0]
                x1 = xa
                x2 = xb
                x3 = x[i0]

            d12 = y1-y2
            d13 = y1-y3
            d23 = y2-y3
            d3 = d13*d23
            d2 = -d12*d23
            d1 = d12*d13

            # Quadratic approximation:

            x0 = (y0-y1)*(y0-y2)*x3/d3+(y0-y1)*(y0-y3)*x2/d2+(y0-y2)*(y0-y3)*x1/d1

            if x0 > x_c:
                xp = x0
            else:
                xm = x0

    return xp,xm

def prgen_fshape(rd,zd,nf):

    nd = len(rd)

    # Find the centroid
    s1 = s2 = s3 = 0.0
    for i in range(nd-1):
        ip = i+1
        dz = zd[ip]-zd[i]       ; dr = rd[ip]-rd[i]
        r0 = 0.5*(rd[ip]+rd[i]) ; z0 = 0.5*(zd[ip]+zd[i])

        # A = Int R dZ = -Int Z dR
        s1 = s1+dz*r0    # Area
        s2 = s2-dr*z0*r0 # Centroid (R)
        s3 = s3+dz*r0*z0 # Centroid (Z)

    r_c = s2/s1 ; z_c = s3/s1

    # Find major radii (rp,rm) at the height of the centroid.

    rp,rm = intersect(rd,zd,z_c)

    rmin = 0.5*(rp-rm) ; rmaj = 0.5*(rp+rm)

    # Find rightmost R (r0)

    r0,z0 = minmax(rd,zd)

    imax = np.argmax(rd)
    r2 = rd[imax]
    z2 = zd[imax]

    if z2 > z0:
        dl0 = np.sqrt((r2-r0)**2+(z2-z0)**2)
    else:
        dl0 = -np.sqrt((r2-r0)**2+(z2-z0)**2)

    l_tot = 0.0
    for i in range(nd-1):
        dl = np.sqrt((rd[i+1]-rd[i])**2+(zd[i+1]-zd[i])**2)
        l_tot = l_tot+dl

    # Construct equally-spaced poloidal angle
    theta  = np.linspace(0,1,nd)*2*np.pi
    dtheta = theta[1]-theta[0]

    ar = np.zeros(nf+1)
    br = np.zeros(nf+1)
    az = np.zeros(nf+1)
    bz = np.zeros(nf+1)

    ds = dtheta/np.pi

    # Trapezoidal integration spectrally-accurate
    for i in range(nd-1):
        for n in range(nf+1):
            y = np.cos(n*theta[i])*rd[i] ; ar[n] = ar[n]+ds*y
            y = np.sin(n*theta[i])*rd[i] ; br[n] = br[n]+ds*y
            y = np.cos(n*theta[i])*zd[i] ; az[n] = az[n]+ds*y
            y = np.sin(n*theta[i])*zd[i] ; bz[n] = bz[n]+ds*y

    return ar,br,az,bz

