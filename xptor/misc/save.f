ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ptor_dv
c
c "PTor" Predict Torus:  Predict Tokamak Profiles and Performance
c by doing simple 1-D power balance calculations,
c using the IFS-PPPL transport model.
c
c Written by Greg Hammett, PPPL, March 17, 1995.
c 

cbasis R.E. Waltz modified for basis March 27, 1995
cbasis converted from equal space rmin to equal space rho->r
cbasis    so elong is function of rho and grad_rho_sq from exp.
cbasis modifing subroutine ptor, ptorinit (ptorexit), trcoef
cbasis and leaving others unchangeds much as possible.
cbasis ..added dinpro ion density dinpro=denpro*ni_exp/ne_exp
cbasis .. changed "rho" weights to "vrho"
cbasis .. added ifix to set vprime(1)->vprime(2) or 
cbasis    vprime(2)->vprime(1) to bring conformity between
cbasis    ptor and  modelt...testing chie=chii=1.
c
      implicit none
c
      include 'input.m'
      include 'ptor.m'
      include 'vtrans.m'
      include 'tport.m'
      include 'model.m'
      include 'glf.m'
c
c LOCAL variables:
      integer i,j,k,m,istep,ierror,isave_pt1
      integer iter,unit,ireject
      real*8 dt_imp,mindtmult,dtmult,maxramp, mindiag
      real*8 maxdelta, mindelta, maximp, deltat, start
      real*8 normS, normT, normDT, save_dt
      real*8 sumnormS,sumnormT,sumnormDT,savenorm,sumnormD
      real*8 test1,test2, rescale, normf, normdf, df
      real*8 diffmin,dx,numax
c
c**********************************************************************
c
c Read input data and initialize:
c
      unit=6
      if(ilog.eq.1)then
        unit=16
        open(unit,file='runlog',status='unknown')
      endif
      j=0
      do i=1,mxfields
        if(itport_pt(i).ne.0)then
          j=j+1
          smult(j)=s0_pt(i)
          if (i_proc .eq. 0)write(unit,*)smult(j)
        endif
        do k=1,mxgrid
          s(i,k) = 0.D0
          flux(i,k) = 0.D0
          flux_s(i,k) = 0.D0
          do m=1,mxfields
            conv3(i,m,k) = 0.D0
            diff(i,m,k)=0.D0
            nu(i,m,k)=0.D0
            nu_p(i,m,k)=0.D0
          enddo
        enddo
      enddo
      nfields = j
      ngrid=jout_m
      if (i_proc .eq. 0)write(unit,20)nfields,ngrid,nsteps_v
 20   format('  nfields = ',i3,'  ngrid = ',i3,
     >   '  nsteps_v = ',i3)
      jin_m=1      
      if(xparam_pt(1).eq.0.D0)xparam_pt(1)=0.5D0
      if(xparam_pt(2).eq.0.D0)xparam_pt(2)=0.3D0
      if(xparam_pt(9).eq.0.D0)xparam_pt(9)=2.D0
      if(iparam_pt(4).eq.0)iparam_pt(4)=10
      mindelta=1.0D-12
      mindtmult=1/2.D0**iparam_pt(4)
      dtmult=1.D0
c
      call ptorinit_dv
c
      if (i_proc .eq. 0) write(unit,*) 'PTOR_DV running ...'
      if (i_proc .eq. 0) write(unit,*) ' '
      dt_imp=dt_v
c
c      if (i_proc .eq. 0) write(6,10) istep,te_m(0),ti_m(0)
c 10   format('step=',i3,', Te0=',f8.4,', Ti0=',f8.4)
c
c get initial sources and transport coefficients
c
c      dt_imp=10000.D0
      call trcoef_dv(dt_imp)
      sumnormS=0.D0
      sumnormT=0.D0
      sumnormDT=0.D0
      do i=1,nfields
        normS=0.D0
        normT=0.D0
        normDT=0.D0
        do k=1,ngrid-1
          flux_s(i,k)=flux(i,k)
          normS = normS + s(i,k)*s(i,k)
          normT = normT + Told(i,k)*Told(i,k)
          normDT = normDT + DABS(s(i,k)*Told(i,k))
        enddo
        normS = DSQRT(normS)
        sumnormS = sumnormS + normS
        normDT = normDT/DMAX1(normT,delt_v)
        sumnormDT = sumnormDT + normDT
      enddo
        normDT = sumnormDT
      istep=0
      if (i_proc .eq. 0)write(unit,11)istep,mindtmult,sumnormS,normDT
 11   format(i4,3(4x,1pE12.4))
      if((endtime_pt.le.xp_time))then
        if (i_proc .eq. 0)
     >  write(unit,*)"Finding steady state at ",xp_time
        savenorm=sumnormS
c
c************************************************
c
c Start of steady state loop
c
c ************************************************
      do istep=1,nsteps_v
        ireject = 0
        deltat=dtmult
        call vtrans_dv(dt_imp,ierror)
        call advancefields(Tnew,deltat)
        if(deltat.lt.mindelta)then
          if(i_proc.eq.0) write(6,*)"exit",deltat,"<",mindelta
          exit
        endif
 10   continue
c recompute the sources
      isave_pt1=iparam_pt(1)
      iparam_pt(1)=-1
      call trcoef_dv(dt_imp)
      iparam_pt(1)=isave_pt1
c compute normS
      sumnormS=0.D0
      sumnormT=0.D0
      sumnormDT=0.D0
      rescale = 1.D0
      test1=0.D0
      do i=1,nfields
        normS=0.D0
        normT=0.D0
        normDT=0.D0
        normf=0.D0
        normdf=0.D0
        do k=1,ngrid-1
          normS = normS + s(i,k)*s(i,k)
          normT = normT + Told(i,k)*Told(i,k)
          normDT = normDT + DABS(s(i,k)*Told(i,k))
          normf = normf + flux_s(i,k)*flux_s(i,k)
          df = flux(i,k)-flux_s(i,k)
          normdf = normdf + df*df
        enddo
        normS = DSQRT(normS)
        sumnormS = sumnormS + normS
        normDT = normDT/DMAX1(normT,delt_v)
        sumnormDT = sumnormDT + normDT
        test1 = DMAX1(test1,DSQRT(normdf/DMAX1(normf,DFLOAT(ngrid-1))))
      enddo
        normDT = sumnormDT
        test2=DSQRT(sumnormS/savenorm)
c
        if (i_proc .eq. 0)write(unit,11)istep,deltat,sumnormS,normDT
        if(test2.lt.1.D0/(1.D0+2.D0*dtmult))then
          if(dtmult.le.0.5)dtmult=2.D0*dtmult
        endif
        if(test1.gt.xparam_pt(9).or.test2.gt.xparam_pt(9))then
          ireject = ireject+1
          rescale = 1.D0/DMAX1(test1,test2)
          if(i_proc.eq.0)write(unit,*)
     >    " ireject = ",ireject," rescale = ",rescale
          if(ireject.gt.iparam_pt(4))exit
c reject this step
          save_dt = rescale*deltat
          deltat = -deltat+save_dt
          call advancefields(Tnew,deltat)
          deltat = save_dt
          dtmult=0.7071D0*dtmult
          if(dtmult.le.mindtmult)exit
          go to 10
        endif
        if(sumnormS.gt.100*savenorm)then
c
c  adjust stepsize for reduction of normS
c
         do iter = 1,iparam_pt(4)
          dtmult=dtmult/2.D0
          if(dtmult.le.mindtmult)exit
          deltat=-DABS(deltat)/2.D0
          call advancefields(Tnew,deltat)
          isave_pt1=iparam_pt(1)
          iparam_pt(1)=-1
          call trcoef_dv(dt_imp)
          iparam_pt(1)=isave_pt1
          sumnormS=0.D0
          do i=1,nfields
            normS=0.D0
            do k=1,ngrid-1
              normS = normS + s(i,k)*s(i,k)
            enddo
            normS = DSQRT(normS)
            sumnormS = sumnormS + normS
          enddo
          if (i_proc .eq. 0)write(unit,11)istep,deltat,sumnormS
          if(sumnormS.lt.savenorm)exit
         enddo
c
         if(dtmult.le.mindtmult)then
           if(iparam_pt(5).eq.0)exit
         endif
c
        endif
c
c recompute transport coefficients
c
        call trcoef_dv(dt_imp)
c 
        savenorm=sumnormS
        deltat=DABS(deltat)
        do k=1,ngrid-1
          mask_r(k)=0
        do i=1,nfields
          Told(i,k) = Told(i,k)+deltat*Tnew(i,k)
          flux_s(i,k)=flux(i,k)
          if(DABS(deltat*Tnew(i,k)).gt.delt_v)mask_r(k)=1
        enddo
        enddo
         if(ave_dv.gt.0.0)then
c average fields
           call ave_fields
           do i=1,nfields
             k=1
             Told(i,k) = (1.D0-ave_dv)*Told(i,k)+
     >       ave_dv*(Told(i,k+1)+Told(i,k))/2.D0
            do k=2,ngrid-1
              Told(i,k) = (1.D0-ave_dv)*Told(i,k)+
     >        ave_dv*(Told(i,k+1)+Told(i,k-1))/2.D0
            enddo
           enddo
         endif
c
        if(sumnormS .lt. mindelta)exit
        if(normDT .lt. 0.001D0)exit
       enddo
      endif
c
c  save steady state profiles for retarts
c
        ntime_t=0
        time=xp_time
        if(restart_pt.ge.xp_time)time=restart_pt
        time_t(ntime_t)=time
        call timeprofiles_dv
c        if(i_proc.eq.0)then
c          write(*,*)(s(1,j),j=1,ngrid-1)
c          write(*,*)(diff(1,1,j),j=1,ngrid-1)
c        endif
c************************************************
c
c Start of time-step loop
c
c ************************************************
       if(endtime_pt.gt.xp_time)then
        dtmult=mindtmult
        xparam_pt(1)=1.D0
        xparam_pt(2)=1.D0
        maxramp=0.D0
        do i=1,mxfields
         if(itport_pt(i).ne.0)then
          maxramp=DMAX1(maxramp,sramp_pt(i))
         endif
        enddo
        start=0.D0
        if(restart_pt.ge.xp_time)then
         maxramp=maxramp+restart_pt
         start = restart_pt
        endif
        deltat=dt_v
c
c main timestep loop
c
       do istep=1,nsteps_v
        ireject=0
        if(time.lt.maxramp)then
          j=0
          do i=1,mxfields
            if(itport_pt(i).ne.0)then
              j=j+1
              smult(j)=s0_pt(i)
              if(time+dtmult*deltat-start.le.sramp_pt(i))then
                smult(j)=smult(j)+(s1_pt(i)-s0_pt(i))*
     >          (time+dtmult*deltat-start)/sramp_pt(i)
              else
                smult(j)=s1_pt(i)
              endif
              if (i_proc .eq. 0)write(unit,*)smult(j)
            endif
          enddo
         endif
c  make  full timestep
        call vtrans_dv(deltat,ierror)
        call advancefields(Tnew,dtmult)
 100    continue
c recompute the sources
      isave_pt1=iparam_pt(1)
      iparam_pt(1)=-1
      call trcoef_dv(deltat)
      iparam_pt(1)=isave_pt1
cdebug      if(i_proc.eq.0)
cdebug     >write(*,*)"flux=",flux(1,25),"flux_s=",flux_s(1,25),
c check errors
        sumnormS=0.D0
        sumnormT=0.D0
        sumnormDT=0.D0
        rescale=1.D0
        test1=0.D0
       do i=1,nfields
        normS=0.D0
        normT=0.D0
        normDT=0.D0
        normf=0.D0
        normdf=0.D0
        do k=2,ngrid-1
          normS = normS + s(i,k)*s(i,k)
          normT = normT + Told(i,k)*Told(i,k)
          normDT = normDT + DABS(s(i,k)*Told(i,k))
          df = flux(i,k)-flux_s(i,k)
          normf = normf + flux_s(i,k)*flux_s(i,k)
          normdf = normdf + df*df
        enddo
        normS = DSQRT(normS)
        sumnormS = sumnormS + normS
        normDT = normDT/DMAX1(normT,delt_v)
        sumnormDT = sumnormDT + normDT
        test1 = DMAX1(test1,DSQRT(normdf/DMAX1(normf,DFLOAT(ngrid-1))))
cdebug        if(i_proc.eq.0)
cdebug     >   write(*,*)i,"test1=",test1,"normdf=",normdf,"normf=",normf
       enddo
        normDT = sumnormDT
        if (i_proc .eq. 0)write(unit,11)istep,dtmult,sumnormS,normDT
        if(dtmult*deltat.lt.mindelta)exit
        if(test1.gt.xparam_pt(9))then
          ireject = ireject+1
          rescale = 1.D0/test1
          if(i_proc.eq.0)
     >    write(unit,*)" ireject = ",ireject," rescale = ",rescale
          if(ireject.gt.5)exit
c reject this step
          save_dt=dtmult*rescale
          dtmult=-dtmult+save_dt
          call advancefields(Tnew,dtmult)
c try smaller step
          dtmult = save_dt
          go to 100
        endif
        if((normDT*dtmult*deltat .gt. 1.D0)
     >      .and.(dtmult.gt.mindtmult))then
c reduce timestep multiplier dtmult
          dtmult=dtmult/2.D0
c  back up the fields in time
         save_dt = -dtmult
         call advancefields(Tnew,save_dt)
          go to 100
        endif
c advance fields in time 
         do i=1,nfields
          do k=1,ngrid-1
            Told(i,k) = Told(i,k)+dtmult*Tnew(i,k)
          enddo
         enddo
         if(ave_dv.gt.0.0)then
c average fields
           call ave_fields
           do i=1,nfields
             k=1
             Told(i,k) = (1.D0-ave_dv)*Told(i,k)+
     >       ave_dv*(Told(i,k+1)+Told(i,k))/2.D0
            do k=2,ngrid-1
              Told(i,k) = (1.D0-ave_dv)*Told(i,k)+
     >        ave_dv*(Told(i,k+1)+Told(i,k-1))/2.D0
            enddo
           enddo
         endif
c advance time and save profiles
         time=time+dtmult*deltat
         ntime_t=ntime_t+1
         time_t(ntime_t)=time
         if(time_series.ne.0)call transfer_exp(time)
         call timeprofiles_dv
         if (i_proc .eq. 0) write(unit,*)"time = ",time
         if(time.gt.endtime_pt)exit
c adjust timestep upward
        if(normDT*dtmult*deltat .lt. 0.25D0)then
          dtmult=DMIN1(0.5D0,2.D0*dtmult)
        endif
        if(dtmult.lt.mindtmult)dtmult=mindtmult
c recompute the transport coefficients
        call trcoef_dv(deltat)
        dx = 1.D0
        do k=1,ngrid-1
         diffmin=0.D0
         numax=0.D0
        do i=1,nfields
          flux_s(i,k) = flux(i,k)
          do j=1,nfields
           if(diff(i,j,k).gt.diffmin)diffmin=diff(i,j,k)
           if(nu(i,j,k).gt.numax)numax=nu(i,j,k)
          enddo
        enddo
         if(diffmin/numax.lt.dx)dx=diffmin/numax
        enddo
        if(i_proc.eq.0)
     >  write(*,*)"dx = ",dx,dr(1,1)
       enddo
       endif
c**************************debug*********************
       if(i_proc.eq.0)then
       open (9,file='data',status='unknown')
       do k=0,mxgrid
cdebug       write(6,*)nu(1,1,k),theta_exp(k),diff(1,1,k),zptheta_exp(k)
cdebug        write(6,13)k,s(1,k),diff(1,1,k),vphiflux(k),gamma_p_m(k)
cdebug 13     format(2x,i3,2x,1pe14.7,2x,1pe14.7,2x,1pe14.7,2x,1pe14.7)
        write(9,13) k,rho(k),te_exp(k),te_m(k),ti_exp(k),ti_m(k),
     >     egamma_exp(k),egamma_m(k)
 13     format(i4,7(2x,1pe14.7))
       enddo
      close (9)
cdebug      do k=1,ngrid
cdebug      write(unit,*)k,s(1,k),s(2,k),s(3,k),s(4,k)
cdebug      enddo
      if(lprint_gf.eq.10)then
       write(6,14)rlte_gf
       write(6,14)rlti_gf
       write(6,14)rlne_gf
       write(6,14)rlni_gf
       write(6,14)rlnimp_gf
       write(6,14)q_gf
       write(6,14)shat_gf
       write(6,14)alpha_gf
       write(6,14)betae_gf
       write(6,14)taui_gf
       write(6,14)xnu_gf
       write(6,14)rmaj_gf
       write(6,14)rmin_gf
       write(6,14)elong_gf
       write(6,14)amassgas_gf
       write(6,14)dil_gf
       write(6,14)apwt_gf
       write(6,14)aiwt_gf
       write(6,14)gamma_star_gf
       write(6,14)gamma_mode_gf
       write(6,14)xalpha
       write(6,14)alpha_e_gf
       write(6,14)alpha_p_gf
       write(6,14)gamma_e_gf
       write(6,14)gamma_p_gf
       write(6,14)chie_gf
       write(6,14)chii_gf
       write(6,14)xky0_gf
       write(6,14)rms_theta_gf
       write(6,14)xwell_gf
       write(6,14)park_gf
       write(6,14)ghat_gf
       write(6,14)gchat_gf
       write(6,14)adamp_gf
       write(6,14)xkdamp_gf
       write(6,14)amassimp_gf
       write(6,14)zimp_gf
       do k=1,24
         write(6,14)xparam_gf(k)
       enddo
 14    format(2x,1pe18.11)
       write(6,15)eigen_gf
       write(6,15)nroot_gf
       write(6,15)lprint_gf
       do k=1,5
        write(6,15)iflagin_gf(k)
       enddo
 15    format(2x,i5)
       endif
       endif
c ************************************************
c End of Debug output
c
c ************************************************
c
c End of Time Step Loop.
c
c ************************************************
c
      if(unit.eq.16)close(unit)
      return
      end
