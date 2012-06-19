ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine advancefields(T,dT,scale,iflag)
c
c     advance the transported fields in time 
c
      implicit none
c
      include '../inc/input.m'
      include '../inc/ptor.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/glf.m'
c
      real*8 T(mxflds,mxgrd),dT(mxflds,mxgrd)
      real*8 scale,test,rescale
      real*8 normT, normDT
      real*8 x13,x14
      integer i,is,j,k,iflag
c
      x13 = xparam_pt(13)
      x14 = xparam_pt(14)
c
      if(iflag .eq. 1)go to 100
c
c first check for too large of a change in fields
c
      test=0.D0
      j=0
      if(itport_pt(1).ne.0)then
        j=j+1
        do i=1,ngrid-1
          test = dmax1(test,dabs(scale*dT(j,i)/ne_m(i)))
        enddo
        if(test.gt.xparam_pt(2))scale=xparam_pt(2)*scale/test
        test=0.D0
      endif
      if(itport_pt(2).ne.0)then
        j=j+1
        do i=1,ngrid-1
          test = dmax1(test,dabs(scale*dT(j,i)/te_m(i)))
        enddo
        if(test.gt.xparam_pt(2))scale=xparam_pt(2)*scale/test
        test=0.D0
      endif
      if(itport_pt(3).ne.0)then
        j=j+1
        do i=1,ngrid-1
          test = dmax1(test,dabs(scale*dT(j,i)/ti_m(i)))
        enddo
        if(test.gt.xparam_pt(2))scale=xparam_pt(2)*scale/test
        test=0.D0
      endif
      if(itport_pt(4).ne.0)then
        j=j+1
        normT=0.D0
        normDT=0.D0
        do i=1,ngrid-1
          normT = normT + vexb_m(i)*vexb_m(i)*ne_m(i)*ne_m(i)
          normDT = normDT + dT(j,i)*dT(j,i)
c          test = DMAX1(test,scale*ABS(dT(j,i))
c     >          /DMAX1(1.0,ABS(ne_m(i)*vexb_m(i))))
        enddo
        normT = DSQRT(DMAX1(DFLOAT(ngrid-1),normT))
        normDT = DSQRT(normDT)
        test = DMAX1(test,scale*normDT/normT)
        if(test.gt.xparam_pt(2))scale=xparam_pt(2)*scale/test
        test=0.D0
      endif  
      if(itport_pt(5).ne.0)then
        j=j+1
        normT=0.D0
        normDT=0.D0
        do i=1,ngrid-1
          normT = normT + vpol_m(i)*vpol_m(i)*ne_m(i)*ne_m(i)
          normDT = normDT +  dT(j,i)*dT(j,i)
        enddo
        normDT = DSQRT(normDT)
        normT = DSQRT(DMAX1(DFLOAT(ngrid-1),normT))
        test = DMAX1(test,scale*normDT/normT)
        if(test.gt.xparam_pt(2))scale=xparam_pt(2)*scale/test
        test=0.D0
      endif
 100  continue
c
c advance the fields
c
      j=0
      if(itport_pt(1).ne.0)then
c the fast ion and impurity ion density are taken from the data
c so only the electron and main ion densities are changed 
        j=j+1
        do i=1,ngrid-1
          T(j,i) = T(j,i)+scale*dT(j,i)
          ne_m(i) = T(j,i)
          ni_m(i) = fi_m(i)*ne_m(i)
          nz_m(i) = fz_m(i)*ne_m(i)
        enddo
        ne_m(0)=ne_m(1)
        ni_m(0)=ni_m(1)
        nz_m(0)=nz_m(1)
      endif
      if(itport_pt(2).ne.0)then
        j=j+1
        do i=1,ngrid-1
          T(j,i) = T(j,i)+scale*dT(j,i)
          te_m(i) = T(j,i)
        enddo
        te_m(0)=te_m(1)
      endif
      if(itport_pt(3).ne.0)then
        j=j+1
        do i=1,ngrid-1
          T(j,i) = T(j,i)+scale*dT(j,i)
          ti_m(i) = T(j,i)
        enddo
        ti_m(0) = ti_m(1)
      endif
      if(itport_pt(4).ne.0)then
        j=j+1
        do i=1,ngrid-1
          T(j,i) = T(j,i)+x14*scale*dT(j,i)
          vexb_m(i)= T(j,i)/ne_m(i)
        enddo
        vexb_m(0)=vexb_m(1)
      endif
      if(itport_pt(5).ne.0)then
        j=j+1
        do i=1,ngrid-1
          T(j,i) = T(j,i)+x14*scale*dT(j,i)
          vpol_m(i)= T(j,i)/ne_m(i)
        enddo
        vpol_m(0)=vpol_m(1)
      endif
c
      if(alpha_dia.ne.0.0)then
c        call advance_neo_flows(xparam_pt(13))
        call neo_flows(ngrid,vneo_new,vdia_new)
        do i=1,nspecies
          do k=1,ngrid-1
            vneo_m(i,k) = x13*vneo_m(i,k)+(1.0-x13)*vneo_new(i,k)
            vdia_m(i,k) = x13*vdia_m(i,k)+(1.0-x13)*vdia_new(i,k)
          enddo
          vneo_m(i,0) = vneo_m(i,1)
          vdia_m(i,0) = vdia_m(i,1)
        enddo
      endif
c
      do k=1,ngrid-1
        if(itport_pt(4).eq.0)then
          if(irotstab.eq.1)then
c keep vphiz_m fixed and update vexb_m
           vexb_m(k) = (-c_per(k)*(vpol_m(k)+vneo_m(3,k))
     >     +vphiz_m(k))/c_tor(k) -vdia_m(3,k)
          endif
          if(irotstab.eq.2)then
c keep vpar_m fixed and update vexb_m
           vexb_m(k) =(-a_pol(k)*(vpol_m(k)+vneo_m(2,k))
     >     +vpar_m(k))/a_tor(k) - vdia_m(2,k)    
          endif
         endif
c carry the toroidal and parallel velocities along for the ride
          vpare_m(k)=a_pol(k)*(vpol_m(k)+vneo_m(1,k))
     >      +a_tor(k)*(vdia_m(1,k)+vexb_m(k))
          vpar_m(k)=a_pol(k)*(vpol_m(k)+vneo_m(2,k))
     >      +a_tor(k)*(vdia_m(2,k)+vexb_m(k))
          vparz_m(k)=a_pol(k)*(vpol_m(k)+vneo_m(3,k))
     >      +a_tor(k)*(vdia_m(3,k)+vexb_m(k))
          vphie_m(k)=c_per(k)*(vpol_m(k)+vneo_m(1,k))
     >      +c_tor(k)*(vexb_m(k)+vdia_m(1,k))
          vphi_m(k)=c_per(k)*(vpol_m(k)+vneo_m(2,k))
     >      +c_tor(k)*(vexb_m(k)+vdia_m(2,k))
          vphiz_m(k)=c_per(k)*(vpol_m(k)+vneo_m(3,k))
     >      +c_tor(k)*(vexb_m(k)+vdia_m(3,k))
      enddo
c
      vexb_m(0)=vexb_m(1)
      vphi_m(0)=vphi_m(1)
      vphiz_m(0)=vphiz_m(1)
      vphie_m(0)=vphie_m(1)
      vpar_m(0)=vpar_m(1)
      vparz_m(0)=vparz_m(1)
      vpare_m(0)=vpare_m(1)
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advance_neo_flows(x13)
c
c  advance the neoclassical flows with relaxation
c
      implicit none
c
      include '../inc/input.m'
      include '../inc/ptor.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/glf.m'
c
      real*8 x13
      integer i,k
c
      do k=1,ngrid-1
        nem = (ne_m(k+1)+ne_m(k))/2.D0
        tim = (ti_m(k+1)+ti_m(k))/2.D0
        tem = (te_m(k+1)+te_m(k))/2.D0
        fim = (fi_m(k+1)+fi_m(k))/2.D0
        fzm = (fz_m(k+1)+fz_m(k))/2.D0
        nim = fim*nem
        nzm = fzm*nem
        vexbm= 0.0
        vpolm= 0.0
        gradnem = (ne_m(k+1)-ne_m(k))/dr(k,2)
        gradtim = (ti_m(k+1)-ti_m(k))/dr(k,2)
        gradtem = (te_m(k+1)-te_m(k))/dr(k,2)
        gradvexbm = 0.0
        gradvpolm = 0.0
        gradfim = (fi_m(k+1)-fi_m(k))/dr(k,2)
        gradfzm = (fz_m(k+1)-fz_m(k))/dr(k,2)
        gradnim = fim*gradnem + nem*gradfim
        gradnzm = fzm*gradnem + nem*gradfzm
        jm=k
        call neoclassical
c
        do i=1,nspecies
          vneo_m(i,k+1) = x13*vneo_m(i,k+1)+(1.0-x13)*vneo(i)
          vdia_m(i,k+1) = x13*vdia_m(i,k+1)+(1.0-x13)*vdia(i)
        enddo
      enddo
c
      do i=1,nspecies
        vneo_m(i,1) = vneo_m(i,2)
        vdia_m(i,1) = vdia_m(i,2)
        vneo_m(i,0) = vneo_m(i,1)
        vdia_m(i,0) = vdia_m(i,1)
      enddo
c
      return
      end
