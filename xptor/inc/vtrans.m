c@vtrans.m 28-jan-02 J. Kinsey, General Atomics
c 28-jan-02 added wt
c 18-jan-01 added vrho3 after nu_p
c 9-14-99 added conv3 and nu_p for Staebler's version of PTOR
c---------------------------------------------------------------------
      integer mxfds, mxgd
      parameter (mxfds=5, mxgd=300)
c
      real*8 dt_implicit
      real*8 work(6*mxfds-2,mxfds,mxgd)
      real*8 DIFF(mxfds,mxfds,mxgd)
      real*8 nu(mxfds,mxfds,mxgd)
      real*8 nu_p(mxfds,mxfds,mxgd)
      real*8 nu_2(mxfds,mxfds,mxgd)
      real*8 nu_pt(mxfds,mxfds,mxgd)
      real*8 nu_2t(mxfds,mxfds,mxgd)
      real*8 Told(mxfds,mxgd),Tstart(mxfds,mxgd)
      real*8 vrho(mxfds,mxgd)
      real*8 vrho3(mxfds,mxfds,mxgd)
      real*8 S(mxfds,mxgd),St(mxfds,mxgd)
      real*8 INTEGRAL_RHS(mxfds,mxgd)
      real*8 INTEGRAL_LHS(mxfds,mxgd)
      real*8 grow(mxfds,mxgd)
      real*8 CONV(mxfds,mxgd)
      real*8 CONV3(mxfds,mxfds,mxgd)
      real*8 Tnew(mxfds,mxgd)
      real*8 work2(mxgd,2)
      real*8 wp(3)
c
      integer nwork(mxgd*mxfds)
      integer nfields
c
      common /vtranscm/ work, DIFF, nu, nu_p, nu_2 
     & , nu_pt, nu_2t,CONV3
     & , vrho3, Told, Tstart, vrho, S, St
     & , INTEGRAL_RHS, INTEGRAL_LHS, grow, CONV
     & , Tnew, work2, wp, dt_implicit, nwork, nfields
