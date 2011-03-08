c@rungks.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c ivar=1 gamma vs k
c ivar=2 gamma vs a/L_T
c ivar=4 gamma vs beta
c ivar=8 gamma vs magnetic shear
c ivar=13 gamma vs alpha
c
c jvar=1 vary k
c jvar=2 vary a/L_T
c jvar=6 vary Ti/Te
c jvar=7 vary q
c jvar=10 vary L_n
c jvar=12 vary gamma_p
c
c Example: gamma vs a/L_T   L_n=1.75,1.5...0.25
c
c    ivar=2;jvar=10
c    irunmax=31;jrunmax=7
c    xvarmin=3.0;xvarmax=.01
c    zvarmin=1.75,zvarmax=0.25
c    rungks
c
c This runs GKS from a/L_T=0.01 to 3.0 (31 values) at various
c L_n from 0.25 to 1.75 (7 cases)
c
c Note: need mprint from tport.m
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine rungks
c
      use gks_var
      implicit none
      include 'input.m'
      include 'glf.m'
      include 'data_exp.m'
c
      integer irunmx, jrunmx
      parameter (irunmx=100,jrunmx=100)
      integer i, j, irun, jrun
      real xp(irunmx), yp(irunmx,jrunmx)
      real xvar(irunmx),yvar(irunmx,20)
      real zvar(jrunmx)
      real pk_h,shat_h,tprim1_h,tprim3_h,aky1_h,uprim1_h
      real beta_h,temp3_h,fprim1_h,shift_h
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... holders
c
      write(*,*) 'ivar = ',ivar
      write(*,*) 'jvar = ',jvar
      write(*,*) 'irunmax =',irunmax
      write(*,*) 'jrunmax =',jrunmax
      write(*,*) 'xvarmin = ',xvarmin
      write(*,*) 'xvarmax = ',xvarmax
      write(*,*) 'zvarmin = ',zvarmin
      write(*,*) 'zvarmax = ',zvarmax
      write(*,*) 'igraph = ',igraph
      write(*,*) 'mprint5 = ',mprint(5)
c
      write(*,*)
      write(*,*) 'pk = ',pk
      write(*,*) 'shat = ',shat
      write(*,*) 'tprim1 = ',tprim1
      write(*,*) 'tprim3 = ',tprim3
      write(*,*) 'aky1 = ',aky1
      write(*,*) 'uprim1 = ',uprim1
      write(*,*) 'beta = ',beta
      write(*,*) 'temp3 = ',temp3
      write(*,*) 'fprim1 =',fprim1
      write(*,*) 'shift = ',shift
c
      pk_h=pk
      shat_h=shat
      tprim1_h=tprim1
      tprim3_h=tprim3
      aky1_h=aky1
      uprim1_h=uprim1
      beta_h=beta
      temp3_h=temp3
      fprim1_h=fprim1
      shift_h=shift
c
      do jrun=1,jrunmax
        zvar(jrun)=zvarmin+(jrun-1)*(zvarmax-zvarmin)/
     &             (jrunmax-1+1.e-6)
        if(jvar.eq.1) aky1=sqrt(2.0)*sqrt(1.0/temp3)*zvar(jrun)
        if((jvar.eq.1).and.(ivar.eq.6)) aky1=sqrt(2.0)*zvar(jrun)
        if(jvar.eq.2) then
          tprim1=zvar(jrun)
          tprim3=zvar(jrun)
        endif
        if(jvar.eq.4) beta=1/temp3*zvar(jrun)
        if(jvar.eq.6) temp3=1/zvar(jrun)
        if(jvar.eq.7) pk=0.3333*2.0/zvar(jrun)
        if(jvar.eq.8) shat=zvar(jrun)
        if(jvar.eq.10) fprim1=zvar(jrun)
        if(jvar.eq.12) uprim1=2.0/sqrt(2.0)*sqrt(temp3)*zvar(jrun)
c
        do irun=1,irunmax
          xvar(irun)=xvarmin+(irun-1)*(xvarmax-xvarmin)/
     &               (irunmax-1+1.e-6)
          if(ivar.eq.1) aky1=sqrt(2.0)*sqrt(1.0/temp3)*xvar(irun)
          if(ivar.eq.2) then
            tprim1=xvar(irun)
            tprim3=xvar(irun)
          endif
c
          if(ivar.eq.4) beta=1/temp3*xvar(irun)
          if(ivar.eq.6) temp3=1/xvar(irun)
          if((jvar.eq.1).and.(ivar.eq.6)) aky1=aky1*sqrt(1.0/temp3)
c
          if(ivar.eq.7) pk=0.3333*2.0/xvar(irun)
          if(ivar.eq.8) shat=xvar(irun)
          if(ivar.eq.10) fprim1=xvar(irun)
          if(ivar.eq.12) uprim1=2.0/sqrt(2.0)*sqrt(temp3)*xvar(irun)
          if(ivar.eq.13) shift=xvar(irun)
          if(ivar.eq.14) beta=1/temp3*xvar(irun) 
          if(ivar.eq.14) shift=beta/(.010)*2.0*(0.333333/pk)**2/2.0
c
          call gstotal
c
          yvar(irun,1)=agammas(1)
          yvar(irun,2)=dgammas(1)
c
          xp(irun)=xvar(irun)
          if(igraph.gt.0) yp(irun,jrun) = yvar(irun,1)
c
        enddo
c
        write(*,*)
        do j=1,irunmax
          write(*,100) j, xvar(j), yvar(j,1), yvar(j,2)
        enddo
      enddo
c
c... end of do loop
c
      if(mprint(5).eq.-10) then
        write(*,*)
        if(igraph.gt.0) then
          do i=1,irunmax
            do j=1,jrunmax
              write(*,150) i, j, xp(i), yp(i,j)
            enddo
          enddo
        endif
      endif 
c
      pk=pk_h
      shat=shat_h
      tprim1=tprim1_h
      tprim3=tprim3_h
      aky1=aky1_h
      uprim1=uprim1_h
      beta=beta_h
      temp3=temp3_h
      fprim1=fprim1_h
      shift=shift_h
c
 100  format(2x,i2,2x,1p6e14.6)
 150  format(2x,i2,2x,i2,2x,1p6e14.6)
c
      return
      end
