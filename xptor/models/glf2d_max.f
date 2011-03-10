      subroutine glf2d_max(rms_theta_max)
      implicit none
      include '../inc/glf_dimensions.m'
      include '../inc/glf_spec.m'
      include '../inc/glf.m'
      include '../inc/glf_ave_theta.m'
      include '../inc/glf_ave_h.m'
c
      real*8 tp,ta,tb,tm,dt,ga,gb,gm,dg,fm
      real*8 test_wd0,test_wd1
      real*8 test_ws0,test_ws1
      real*8 rms_theta_save,tmax,tmin,gamma_max
      real*8 nbsave,freqs, save_park, dgamma
      real*8 park_factor, drift
      integer i,nt,imax,ntm,nbmax
c      parameter(nt=41)
c      parameter(nt=13)
      parameter(nt=20,nbmax=4)
      real*8 dg1,dg2,dgmin, rms_theta_max
      real*8 gamma_n(nt),freq_n(nt),rms_theta_n(nt)
      integer gmax(nt),igamma,branch
c
      branch = ibranch_gf
      save_park =park_gf
      nbsave = nbasis_gf
      if(iflagin_gf(1).ne.0)nbasis_gf = iflagin_gf(1)
      rms_theta_save = rms_theta_gf
      rms_theta_gf=DABS(rms_theta_gf)
      drift = shat_gf-alpha_gf-0.5D0 
      park_factor=1.D0 + xparam_gf(4)*DABS(drift)
      if(drift .gt.0.D0)then
        tmax = DLOG10(rms_theta_gf*park_factor)
c        tmax = DLOG10(DABS(rms_theta_gf - 0.3D0*drift))
      else 
        tmax = DLOG10(rms_theta_gf)
      endif
      tmin = DLOG10(0.2D0)
c
      open(10,file = 'rms_theta.dat',status='unknown')
      dt = (tmax-tmin)/FLOAT(nt-1)
      tp = tmin
      rms_theta_gf = 10.D0**tp
      park_gf = save_park*park_factor
      freqs=0.D0
      test_wd0 = 0.D0
c      test_ws0 = 0.D0
      do i=1,nt
c        write(*,*)"rms_theta_gf = ",rms_theta_gf
c       write(*,*)"nbasis_gf = ",nbasis_gf
       if(xparam_gf(4).gt.0) write(*,*)"park_gf = ",park_gf
       call glf2d
       if(ibranch_gf.eq.0)then
        branch = 1
        if(gamma_gf(2).gt.gamma_gf(1))branch = 2
       endif
       test_wd1 = rms_theta_gf*DABS(ave_wd(1,1))
       write(10,*)i,rms_theta_gf,gamma_gf(branch),freq_gf(branch)
c       write(*,*)i,rms_theta_gf,gamma_gf(branch),freq_gf(branch)
       rms_theta_n(i)=rms_theta_gf
       gamma_n(i) = gamma_gf(branch)
       freq_n(i) = freq_gf(branch)
c       if(test_wd1.lt.test_wd0)exit
c       if(test_ws1.lt.test_ws0)exit
       test_wd0 = test_wd1
c       test_ws0 = test_ws1
        tp = tp + dt
       rms_theta_gf = 10.D0**tp
       park_gf = save_park*park_factor
      enddo
      ntm = MIN0(i,nt)
      if(ntm.lt.3)ntm=3
      close(10)
c find local maxima gamma
      dg1 = gamma_n(2)-gamma_n(1)
      imax=0
c      if(dg1.lt.0.D0)then
c        imax = 1
c        gmax(imax)=1
c      endif
      do i=3,ntm
        dg2 = gamma_n(i)-gamma_n(i-1)
        if(dg1.gt.0.D0.and.dg2.lt.0.D0)then
          imax = imax+1
          gmax(imax)=i-1
        endif
        dg1 = dg2
      enddo
      if(dg1.gt.0.D0)then
         imax = imax+1
         gmax(imax) = ntm
      endif
      if(imax.eq.0)then
         gmax(1)=ntm
         imax = 1
      endif
c      write(*,*)"**** local maxima ****"
c      do i=1,imax
c        write(*,200)gmax(i),gamma_n(gmax(i)),freq_n(i),
c     >     rms_theta_n(gmax(i))
c      enddo
c
       igamma = gmax(1)
       if(imax.gt.1)then
c find top maxima
          gamma_max = -1.D0
          do i=1,imax
            if(gamma_n(gmax(i)).gt.gamma_max)then
              gamma_max = gamma_n(gmax(i))
              igamma = gmax(i)
            endif
          enddo
      endif
       rms_theta_gf = rms_theta_n(igamma)
c       write(*,*)"**** maximum gamma ****"
c       write(*,*)"gamma_gf = ",gamma_n(igamma)
c       write(*,*)"freq_gf = ",freq_n(igamma)
c       write(*,*)"rms_theta_gf = ",rms_theta_gf
       gamma_gf(branch) = gamma_n(igamma)
       freq_gf(branch) = freq_n(igamma)
       rms_theta_max=rms_theta_gf
       if(gamma_n(igamma).ne.0.D0.and.nbasis_gf.ne.nbsave)then
c refine eigenvalue with more basis functions
        nbasis_gf=nbsave
       park_gf = save_park*park_factor
       park_gf = 1.D0 + xparam_gf(5)
 100    call glf2d
c
        if(ibranch_gf.eq.0)then
         branch = 1
         if(gamma_gf(2).gt.gamma_gf(1))branch = 2
        endif
c
        dgamma = (gamma_gf(branch)-gamma_n(igamma))**2/
     >           DMAX1(gamma_n(igamma)**2,0.001D0)
        dgamma = dgamma +
     >   (freq_gf(branch)-freq_n(igamma))**2/
     >    DMAX1(freq_n(igamma)**2,0.001D0)
        dgamma = DSQRT(dgamma/2.D0)
c        write(*,*)"nbasis_gf = ",nbasis_gf
c        write(*,*)"gamma_gf = ",gamma_gf(branch),
c     >            " freq_gf = ",freq_gf(branch)
c        write(*,*)"rms change in eigenvalue = ",dgamma
        if(dgamma.gt.xparam_gf(1).and.nbasis_gf.lt.nbmax)then
         gamma_n(igamma) = gamma_gf(branch)
         freq_n(igamma) = freq_gf(branch)
         nbasis_gf = nbasis_gf+2
c         write(*,*)"retry with nbasis_gf = ",nbasis_gf
         go to 100
        endif     
c        write(*,*)" maximum gamma for nbasis = ",nbasis_gf
c        write(*,*)"gamma_gf(1) = ",gamma_gf(1)
c        write(*,*)"freq_gf(1) = ",freq_gf(1)
c        write(*,*)"gamma_gf(2) = ",gamma_gf(2)
c        write(*,*)"freq_gf(2) = ",freq_gf(2)
c        gamma_n(igamma)=gamma_gf(branch)
c        freq_n(igamma) = freq_gf(branch)
       endif
c
        park_gf = save_park
        nbasis_gf = nbsave
c
 200  format(i2,2x,0p6f10.5)
      return
      end   
