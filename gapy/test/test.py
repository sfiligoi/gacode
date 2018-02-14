import numpy as np
import gapy

#print gapy.p_bound_deriv.__doc__
#print gapy.p_cub_spline_deriv.__doc__

r = np.arange(0.1,1,0.1)
f = r**2
n = len(r)
df = np.zeros(n)

print r
print f

gapy.p_bound_deriv(df,f,r)

print df/(2*r)

gapy.p_cub_spline_deriv(r,f,df)

print df/(2*r)
