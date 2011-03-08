c@ut.m 28-Aug-00 G.M. Staebler
c 1-22-01 jek added rho_xp(mxgrid+1)
cc
      integer mxgrd
      parameter (mxgrd=200)
      integer nr,nut
      real*8 rho_xp(mxgrd+1),ut(ntmax),u1d_t(ntmax,n1d)
     > ,u2d_t(mxgrd,ntmax,n2d)
      common /ut/ u2d_t,u1d_t,ut,rho_xp,nr,nut
