      real*8 gradB
      real*8 ave_wd(nb,nb),ave_modwd(nb,nb)
      real*8 ave_gradB(nb,nb),ave_lnB(nb,nb)
      real*8 ave_k2(nb,nb)
      real*8 ave_kpar(nb,nb),ave_kpar2(nb,nb)
      real*8 ave_modkpar(nb,nb),ave_kparinv(nb,nb)
      real*8 ave_p0(nb,nb),ave_p0inv(nb,nb)
      common /avetheta/
     > gradB,ave_wd,ave_modwd,
     > ave_gradB,ave_lnB,
     > ave_k2,
     > ave_kpar,ave_kpar2,
     > ave_modkpar,ave_kparinv,
     > ave_p0,ave_p0inv

