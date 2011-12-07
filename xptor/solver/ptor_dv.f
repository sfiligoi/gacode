cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      use tglf_tg
c
      implicit none
c
      include '../inc/input.m'
      include '../inc/ptor.m'
      include '../inc/vtrans.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/glf.m'
c
c LOCAL variables:
      integer i,j,k,m,istep,ierror,isave_pt1
      integer iendstep,iretry,nloop,gridcount
      integer iter,unit,ireject,ifail,iphase,mask_p(301)
      integer tglf_restarts
      real*8 dt_imp,dt_newton,dtmult,maxramp, mindiag
      real*8 maxdelta, mindelta, maximp, deltat, start
      real*8 normS, normDT, normG, save_dt, normS_last
      real*8 STEST,Stest2,sumS
      real*8 test1,test2, rescale, normf, normdf
      real*8 save_s(5,301),df(5,301),normT(5)
      real*8 gamma_net_s(301),gamma_net,dT(5,301),deltap
      real*8 errf,ddf,savenorm,dt_min,dt_max,startnormS
      real*8 dt_courant,test
      real*8 save_nb_tg
      LOGICAL save_new_eikonal_tg
      parameter(STEST=1.0D-4)
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
        INTEGRAL_LHS(i,0)=0.0
        INTEGRAL_RHS(i,0)=0.0
        do k=1,mxgrid
          s(i,k) = 0.D0
          flux(i,k) = 0.D0
          flux_s(i,k) = 0.D0
          Tstart(i,k)=0.D0
          St(i,k) = 0.D0
          INTEGRAL_LHS(i,k)=0.0
          INTEGRAL_RHS(i,k)=0.0
          do m=1,mxfields
            conv3(i,m,k) = 0.0
            diff(i,m,k) = 0.0
            stiff(i,m,k) = 0.0
            nu(i,m,k) = 0.0
            nu_p(i,m,k) = 0.0
            nu_2(i,m,k) = 0.0
            nu_pt(i,m,k) = 0.0
            nu_2t(i,m,k) = 0.0
            vrho3(i,m,k) = 0.0
          enddo
        enddo
      enddo
      nfields = j
      ngrid=jout_m
      Stest2 = STEST/(2.D0*DBLE(nfields)*DSQRT(DBLE(ngrid)))
      nloop=iparam_pt(2)
      if (i_proc .eq. 0)write(unit,20)nfields,ngrid,nsteps_v
 20   format('  nfields = ',i3,'  ngrid = ',i3,
     >   '  nsteps_v = ',i3)
c      jin_m=1      
      if(xparam_pt(1).eq.0.D0)xparam_pt(1)=0.5D0
      if(xparam_pt(2).eq.0.D0)xparam_pt(2)=0.5D0
      if(xparam_pt(9).eq.0.D0)xparam_pt(9)=2.D0
      dtmult=1.D0
      if(xparam_gf(27).ne.0.D0)dtmult=xparam_gf(27)
      if(i_proc.eq.0) write(*,*)"dtmult = ",dtmult
c      dt_min=1.D-7
c
c  dt_imp is the full implicit timestep in seconds,
c  dtmult is the maximum fraction of dt_imp that the fields will be advanced
c  deltat is the adaptive fraction of dt_imp wich is actually used
c          to advance the fields
c
c      errf=xparam_pt(9)
      errf=1.D-8
      do k=1,ngrid-1
        mask_r(k)=1
        mask_p(k)=0
      enddo
      if (i_proc .eq. 0) write(unit,*) 'PTOR_DV running ...'
      if (i_proc .eq. 0) write(unit,*) ' '
c
      call ptorinit_dv
c
      do i=1,nfields
      do k=1,ngrid
        Tstart(i,k) = Told(i,k)
      enddo
      enddo
c
c start up TGLF 
      if(imodel.eq.82)then
        call tglf_startup
        save_new_eikonal_tg=new_eikonal_tg
        if(new_eikonal_tg.eqv. .FALSE.)new_eikonal_tg=.TRUE.
        save_nb_tg=nbasis_max_tg
c        nbasis_max_tg=1
c        nbasis_min_tg=1   ! possible side effect here 
        tglf_restarts=0
        if(i_proc.eq.0)
     >  write(unit,*)"startin TGLF with nbasis = ",nbasis_max_tg
      endif
c
c get initial sources and transport coefficients
c
      dt_imp = DABS(dt_v)
      iparam_pt(1)=1
c      if(nsteps_v.eq.0)iparam_pt(1)=-1
 100  call trcoef_dv(dt_imp)
      if(imodel.eq.82)then
       new_eikonal_tg=save_new_eikonal_tg
       if(i_proc.eq.0)write(*,*)"finished initializing eikonal"
      endif
      call norm_dv(normS,normDT,normT,ngrid)
        if(i_proc.eq.0)write(unit,*)"zero",normS,normDT
      do k=1,ngrid-1
        gamma_net_s(k) = anrate_m(k)-alpha_e_gf*DABS(egamma_m(k))
        do i=1,nfields
          flux_s(i,k)=flux(i,k)
          save_s(i,k)=s(i,k)
          df(i,k) = 1.D0
        enddo
      enddo
c
c  compute the courant timestep  
c
      dt_courant=0.D0
      do i=1,nfields
       do k=1,ngrid
        test = DABS(diff(i,i,k)/(dr(k,1)*dr(k,2)))
        if(test.gt.dt_courant)then
         dt_courant = test 
        endif
       enddo
      enddo
      dt_courant = 1.0/MAX(dt_courant,1.0)
c      write(*,*)"dt_courant = ",dt_courant
      dt_max = dt_imp
      if(dt_v.gt.0.D0)dt_max=DMIN1(0.2/normDT,dt_imp)
c      dt_min = dt_max/256.D0
      dt_min = dt_courant/10.D0
      dt_imp=dt_max
      if(i_proc.eq.0)
     >  write(*,*)"starting timestep = ",dt_imp," minimum = ",dt_min
      istep=0
      ireject=0
      if (i_proc .eq. 0)
     > write(unit,11)istep,ireject,dtmult,normS,normDT
 11   format(2(i4),3(4x,1pE12.4))
c      if((endtime_pt.le.xp_time))then
c        if (i_proc .eq. 0)
c     >  write(unit,*)"Finding steady state at ",xp_time
        savenorm=normS
c
c************************************************
c
c Start of main loop
c
c ************************************************
       ntime_t=0
       if(restart_pt.ge.0.D0)xp_time=xp_time+restart_pt
c use time relative to xp_time 
       time=0.D0
       time_t(ntime_t)=time
       call timeprofiles_dv
       if (i_proc .eq. 0)
     >   write(unit,*)"ntime = ",ntime_t,"time = ",time
       if(nsteps_v.eq.0)then
         norm_s=normS
         norm_dt=normDT
         return
       endif
c iendstep = 0 still trying
c iendstep = 1 converged
c iendstep = 2 failed
      iendstep=0
      iretry=0
      startnormS = normS
      deltat=xparam_pt(1)
 101  continue
        istep=istep+1
        ireject = 0
        ifail=0
        deltat=xparam_pt(1)
        dt_newton=dt_imp*dtmult
        call vtrans_dv(dt_newton,ierror)
c        call vtrans_dv(dt_imp,ierror)
        if(ifilter.eq.1)call filter_dv
        call advancefields(Told,Tnew,deltat,0)
c recompute the fluxes and adjust deltat
 31     continue
      isave_pt1=iparam_pt(1)
      iparam_pt(1)=-1
      call trcoef_dv(dt_imp)
      iparam_pt(1)=isave_pt1
      call norm_dv(normS,normDT,normT,ngrid)
        if (i_proc .eq. 0)
     >  write(unit,11)istep,ireject,deltat,normS,savenorm
c check for pre-convergence of TGLF with nbasis_min
       if(imodel.eq.883)then
         if(normS.lt.100*STEST)then
           if(tglf_restarts.eq.0)then
            if(nbasis_max_tg.ne.save_nb_tg)then
              nbasis_max_tg = save_nb_tg
              if(i_proc .eq. 0)
     >        write(unit,*)"switching to nbasis = ",nbasis_max_tg
              tglf_restarts=1
            endif
            if(new_eikonal_tg.eqv. .FALSE.)then
              new_eikonal_tg=.TRUE.
              if(i_proc.eq.0)write(*,*)"recomputing eikonal"
              tglf_restarts=1
            endif
            if(tglf_restarts.eq.1)then
              call trcoef_dv(dt_imp)
              call norm_dv(normS,normDT,normT,ngrid)
              new_eikonal_tg = save_new_eikonal_tg
              savenorm=normS
              if (i_proc .eq. 0)
     >        write(unit,11)istep,ireject,deltat,normS,savenorm
              go to 101
            endif
          endif
         endif
       endif    
c check for convergence
      if(normS.lt.DMAX1(normDT*STEST,STEST))then
        if(i_proc.eq.0)write(*,*)"converged ",normS
        iendstep=1
        go to 10
      endif 
c check to see if going uphill
c      if(ireject.gt.0.and.normS.gt.normS_last)then
c choose scale between last two and exit reject loop
c        deltat=1.5D0*deltat
c        call advancefields(Told,Tnew,deltat,1)
c        if(i_proc.eq.0)
c     >   write(unit,*)"uphill ",normS,norms_last
c        ireject=0
c        go to 10
c      endif
      if(ireject.lt.iparam_pt(3).and.
     >    normS.gt.xparam_pt(9)*savenorm)then
c do a line search along the newton direction
        ireject = ireject+1
c        if (i_proc .eq. 0)
c     >  write(unit,*)"reject",istep,ireject,deltat,normS,savenorm
        rescale = deltat/2.D0 - deltat
        call advancefields(Told,Tnew,rescale,1)
        deltat=deltat/2.D0        
c        iparam_pt(1)=isave_pt1
        normS_last=normS
        go to 31
      endif
      if(dt_v.lt.0.D0.and.istep.ge.nloop)iendstep=1
 10   continue
c      if (i_proc .eq. 0)
c     > write(unit,11)istep,ireject,deltat,normS,normDT
      if(dt_imp.lt.dt_min)then
        iendstep=2  
        dt_imp=dt_min
        if(i_proc.eq.0)
     >   write(unit,*)"timestep failed, resetting dt = ",dt_imp
      endif
      if(iendstep.gt.0)then
        if(i_proc.eq.0)
     >   write(unit,*)"end of timestep"
        do i=1,nfields
        do k=1,ngrid
          Tstart(i,k) = Told(i,k)
        enddo
        enddo
        if(istep.eq.1.and.dt_imp.lt.dt_max)then
         dt_imp=1.1D0*dt_imp
        if(i_proc.eq.0)
     >    write(unit,*)"increasing timestep = ",dt_imp
        endif
      endif
      if((iendstep.eq.0.and.istep.ge.nloop))then
c try again from the beginning
        istep=0
        ireject=0
        if((normS.gt.normDT/10.D0.or.iretry.gt.0))then
         normS=startnormS
         dt_imp=dt_imp/2.D0
         do i=1,nfields
         do k=1,ngrid-1
           Tnew(i,k) = Told(i,k) - Tstart(i,k)
         enddo
         enddo
         deltat=-1.D0
         call advancefields(Told,Tnew,deltat,1)
         deltat=xparam_pt(1)
         if(i_proc.eq.0)
     >    write(unit,*)"retry with timestep = ",dt_imp
         iretry=0
        else
         iretry=1
        endif
      endif
c
c recompute transport coefficients
c
        gridcount=ngrid-1
cc        do k=1,ngrid-1
cc          mask_r(k)=1
cc            sumS=0.D0
cc            do i=1,nfields
cc              sumS = sumS + DABS(S(i,k))/normT(i)
cc            enddo
cc            if(sumS.lt.DMAX1(normDT*Stest2,Stest2))then
cc              mask_r(k)=0
cc              gridcount=gridcount-1
cc            endif
cc        enddo
cc        if(i_proc.eq.0)write(*,*)"active gridpoints = ",gridcount
c
      call trcoef_dv(dt_imp)
c        write(*,*)"third call", nglf23
      call norm_dv(normS,normDT,normT,ngrid)
      savenorm=normDT
      if(dt_v.lt.0.D0)savenorm=normS
c
      if(iendstep.eq.0)go to 101
c
c advance time and save profiles
         time=time+dt_imp
         ntime_t=ntime_t+1
         time_t(ntime_t)=time
         if(time_series.ne.0)call transfer_exp(time)
         call timeprofiles_dv
         j=0
         do i=1,mxfields
           if(time.lt.sramp_pt(i))then
             if(itport_pt(i).ne.0)then
              j=j+1
              smult(j)= smult(j)+(s1_pt(i)-s0_pt(i))*time/sramp_pt(i)
              if(i_proc.eq.0)write(unit,*)i,j,"smult",smult(j),time
             endif
           endif
         enddo
         if (i_proc .eq. 0)then
           write(unit,11)istep,ireject,deltat,normS,normDT
           norm_s=normS
           norm_dt=normDT
           write(unit,*)"ntime = ",ntime_t,"time = ",time
         endif
c         if(time.gt.endtime_pt)return
      istep=0
      iendstep=0
      iretry=0
      savenorm=normDT
      if(dt_v.gt.0.D0)dt_max= DMIN1(0.2D0/normDT,dt_v)
      dt_courant=0.D0
      do i=1,nfields
       do k=1,ngrid
        test = DABS(diff(i,i,k)/(dr(k,1)*dr(k,2)))
        if(test.gt.dt_courant)then
         dt_courant = test 
        endif
       enddo
      enddo
      dt_courant = 1.0/dt_courant
      dt_min=dt_courant/10.D0
      if(dt_imp.lt.dt_min)then
       dt_imp = dt_min
       write(*,*)"increasing dt = ",dt_imp
      endif
      deltat=xparam_pt(1)
c reset nbasis to 1 for TGLF
c      if(imodel.eq.82)then
c        nbasis_max_tg=1
c        if(i_proc.eq.0)
c     >  write(unit,*)"starting TGLF with nbasis = ",nbasis_max_tg
c      endif
      if(ntime_t.lt.nsteps_v.and.normDT.gt.1.D-3)go to 101
c ***************************************************
c
c End of main Loop.
c
c**************************debug*********************
       if(i_proc.eq.0.and.lprint_gf.eq.10)then
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
      if(unit.eq.16)close(unit)
      return
      end
