      real*8 ft
      real*8 wdx(nxm),k2x(nxm)
      real*8 p0x(nxm),Bx(nxm)
      real*8 hxn(ns,nxm),hxp1(ns,nxm),hxp3(ns,nxm)
      real*8 hxr11(ns,nxm),hxr13(ns,nxm),hxr33(ns,nxm)
      real*8 gxn(ns,nxm),gxp1(ns,nxm),gxp3(ns,nxm)
      real*8 gxr11(ns,nxm),gxr13(ns,nxm),gxr33(ns,nxm)
      common /x_grid/
     > ft,
     > wdx,k2x,p0x,Bx,
     > hxn,hxp1,hxp3,
     > hxr11,hxr13,hxr33,
     > gxn,gxp1,gxp3,
     > gxr11,gxr13,gxr33

