      real*8 gradBgp1,gradBgp3
      real*8 gradBgr11,gradBgr13,gradBgr33
      real*8 gradBgu1,gradBgu3,gradBgu33
      real*8 ave_gradBgp1(ns,nb,nb),ave_gradBgp3(ns,nb,nb)
      real*8 ave_gradBgr11(ns,nb,nb),ave_gradBgr13(ns,nb,nb)
      real*8 ave_gradBgr33(ns,nb,nb)
      real*8 ave_gradBgu1(ns,nb,nb),ave_gradBgu3(ns,nb,nb)
      real*8 ave_gradBgu33(ns,nb,nb)
      common /gradBg/
     > gradBgp1,gradBgp3,gradBgr11,gradBgr13,gradBgr33,
     > gradBgu1,gradBgu3,gradBgu33,
     > ave_gradBgp1,ave_gradBgp3,
     > ave_gradBgr11,ave_gradBgr13,ave_gradBgr33,
     > ave_gradBgu1,ave_gradBgu3,ave_gradBgu33 
  
