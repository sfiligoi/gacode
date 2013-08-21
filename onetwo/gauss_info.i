
c     ngauss_max sets the maximum number of points used in the
c     integrations involved in the time dependent beam calcs.
c     The actuall number of points used is user determinable and
c     is given by ngauss.
         integer  
     .         ngauss,ngauss_max
      parameter (ngauss_max = 100 )
c
      real *8
     .        tau0_x(ngauss_max),tau0_w(ngauss_max)
      common /Gauss_Leg/
     .       tau0_x,tau0_w,ngauss
