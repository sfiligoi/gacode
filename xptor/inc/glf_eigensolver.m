      integer neq,iar
      parameter(neq = 15*ns,iar=neq*nb)
      integer matz,ieq,iur,nroot
      real*8 fv1(iar),fv2(iar),fv3(iar)
      real*8 rr(iar), ri(iar)
      real*8 ar(iar,iar), ai(iar,iar)
      real*8 vr(iar,iar), vi(iar,iar)
      common / eigen /
     > matz,ieq,iur,nroot,
     > fv1,fv2,fv3,rr,ri,ar,ai,vr,vi
