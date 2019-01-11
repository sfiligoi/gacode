import gapy
import numpy as np
from scipy import integrate,optimize,interpolate

data=np.loadtxt('out.locpargen')
gapy.geo.geo_q_in = data[0]
gapy.geo.geo_rmin_in = data[1]
gapy.geo.geo_rmaj_in = data[2]
gapy.geo.geo_drmaj_in = data[3]
gapy.geo.geo_kappa_in = data[4]
db = data[5]

n = 32
     
theta=np.linspace(0,2*np.pi,n)-np.pi
gapy.geo.geo_interp(theta,True)

b_gapy = gapy.geo.geo_b
gsin_gapy = gapy.geo.geo_gsin
gsin_gapy[0] = 0.0 ; gsin_gapy[-1] = 0.0
gt_gapy = gapy.geo.geo_g_theta
gr_gapy = gapy.geo.geo_grad_r
    
b_s    = interpolate.CubicSpline(theta,b_gapy,bc_type='periodic')
gsin_s = interpolate.CubicSpline(theta,gsin_gapy,bc_type='periodic')
gt_s   = interpolate.CubicSpline(theta,gt_gapy,bc_type='periodic')
gr_s   = interpolate.CubicSpline(theta,gr_gapy,bc_type='periodic')

def func(t):

   b = b_s(t)
   f = (1+l*b/2)/np.sqrt(1.0-l*b)/b*gsin_s(t)*gt_s(t)*gr_s(t)
   return f

def bounce(t):
   return 1-l*b_s(t)

lb = 1.0/b_s(np.pi)
lmax = 1.0/b_s(0.0)
l = lb*0.999

if l < lb:
   t0 = np.pi
elif l < lmax:
   t0 = optimize.fsolve(bounce,0.01)[0]
else:
   print 'l too large'
   sys.exit()

y,err = integrate.quad(func,0,t0)

print '* orbit width/a = {:.3f}'.format(np.abs(y*db))

