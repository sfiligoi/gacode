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
      real*8 x13,x14,reference,limit
      real*8 rminm,rhom,doppler_shearm
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
          test = dmax1(test,dabs(scale*dT(j,i)/T(j,i)))
        enddo
        if(test.gt.xparam_pt(2))scale=xparam_pt(2)*scale/test
        test=0.D0
      endif
      if(itport_pt(2).ne.0)then
        j=j+1
        do i=1,ngrid-1
          test = dmax1(test,dabs(scale*dT(j,i)/T(j,i)))
        enddo
        if(test.gt.xparam_pt(2))scale=xparam_pt(2)*scale/test
        test=0.D0
      endif
      if(itport_pt(3).ne.0)then
        j=j+1
        do i=1,ngrid-1
          test = dmax1(test,dabs(scale*dT(j,i)/T(j,i)))
        enddo
        if(test.gt.xparam_pt(2))scale=xparam_pt(2)*scale/test
        test=0.D0
      endif
      if(itport_pt(4).ne.0)then
        j=j+1
c        normT=0.D0
c        normDT=0.D0
c        do i=1,ngrid-1
c          normT = normT + T(j,i)**2
c          normDT = normDT + dT(j,i)**2
c        enddo
c        normT = DSQRT(DMAX1(DFLOAT(ngrid-1),normT))
c        normDT = DSQRT(normDT)
c        test = DMAX1(test,scale*normDT/normT)
        do i=1,ngrid-1
          test = dmax1(test,dabs(scale*dT(j,i)/DMAX1(1.0,DABS(T(j,i)))))
        enddo
        if(test.gt.xparam_pt(2))scale=xparam_pt(2)*scale/test
        test=0.D0
      endif  
      if(itport_pt(5).ne.0)then
        j=j+1
        normT=0.D0
        normDT=0.D0
        do i=1,ngrid-1
          normT = normT + T(j,i)**2
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
          ne_m(i) = T(j,i)/h_m(i)
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
          te_m(i) = T(j,i)/h_m(i)
        enddo
        te_m(0)=te_m(1)
      endif
      if(itport_pt(3).ne.0)then
        j=j+1
        do i=1,ngrid-1
          T(j,i) = T(j,i)+scale*dT(j,i)
          ti_m(i) = T(j,i)/h_m(i)
        enddo
        ti_m(0) = ti_m(1)
      endif
      if(itport_pt(4).ne.0)then
        j=j+1
        do i=1,ngrid-1
          T(j,i) = T(j,i)+scale*dT(j,i)
          vexb_m(i)= T(j,i)/h_m(i)
        enddo
        vexb_m(0)=vexb_m(1)
      endif
      if(itport_pt(5).ne.0)then
        j=j+1
        do i=1,ngrid-1
          T(j,i) = T(j,i)+scale*dT(j,i)
          vpol_m(i)= T(j,i)/h_m(i)
        enddo
        vpol_m(0)=vpol_m(1)
      endif
c
c advance ExB velocity shear and doppler shear with relaxation
      do k=1,ngrid-1
        rminm=(rmin_exp(k+1)+rmin_exp(k))/2.0
        rhom=arho_exp*(rho(k+1)+rho(k))/2.0
        doppler_shearm = -rminm/rhom*drhodr(k)*theta_exp(k)*
     >  (vexb_m(k+1)-vexb_m(k))/dr(k,2)
c        reference = 9.79D5*DSQRT(ti_exp(k)*1.0D3/amassgas_exp)
c     >   /(rmaj_exp(k)*100.0)/cv
        reference = anrate_m(k)*csda_m(k)/cv
        limit = 1.0
        test = 0.0
        if(reference.gt.1.0D-10)then
           test = ABS(doppler_shearm-doppler_shear_m(k))/reference
        endif
        if(test.gt.x14)limit = x14/test
        doppler_shear_m(k) = limit*doppler_shearm
     >   +(1.0-limit)*doppler_shear_m(k)
      enddo
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
          vneo_m(i,ngrid) = vneo_m(i,ngrid-1)
          vdia_m(i,ngrid) = vdia_m(i,ngrid-1)
        enddo
      endif
c
      do k=1,ngrid-1
        if(itport_pt(4).eq.0)then
          if(irotstab.eq.1)then
c keep cer ion toroidal velocity fixed and update vexb_m
           if(cer_ion_exp.eq.1)then
           vexb_m(k) = (-c_per(k)*(vpol_m(k)+vneo_m(2,k))
     >     +vphi_m(k))/c_tor(k) -vdia_m(2,k)
           elseif(cer_ion_exp.eq.2)then
             vexb_m(k) = (-c_per(k)*(vpol_m(k)+vneo_m(3,k))
     >       +vphiz_m(k))/c_tor(k) -vdia_m(3,k)
           endif
          endif
          if(irotstab.eq.2)then
c keep cer ion parallel velocity fixed and update vexb_m
            if(cer_ion_exp.eq.1)then
              vexb_m(k) =(-a_pol(k)*(vpol_m(k)+vneo_m(2,k))
     >        +vpar_m(k))/a_tor(k) - vdia_m(2,k)  
            elseif(cer_ion_exp.eq.2)then
              vexb_m(k) =(-a_pol(k)*(vpol_m(k)+vneo_m(3,k))
     >        +vparz_m(k))/a_tor(k) - vdia_m(3,k)  
           endif  
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
