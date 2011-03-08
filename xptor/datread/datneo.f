c@datneo.f
c jek 18-Jan-11
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... calls kapisn neoclassical using exp data in _d variables
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine datneo
c
      implicit none
      include 'mpif.h'
      include '../inc/tport.m'
      include '../inc/data.m'
      include '../inc/model.m'
      include '../inc/input.m'
      include '../inc/glf.m'
c
      integer j
      real*8 xkapi_nc(nj), xnstare_nc(nj), xnstari_nc(nj), 
     &       zfluxlim_nc(nj), zrhoi_nc(nj), ztaui_nc(nj), 
     &       delta_nc(nj), xkmneo_d(nj)
      real*8 rhob_nc(nj,3)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      if(i_proc.eq.0) write(*,30) use_xneo_m
c
c... simple artificial chi at small r
c
      if (use_xneo_m.eq.-1) then
        if(i_proc.eq.0) write(*,20)
        do j=1,nj_d
          xkineo_d(j)=dexp(-(rho(j)/xchi)**chiaddexp)
        enddo
      endif
c
c... kapisn model
c
      ng_nc=1                     ! number of hyd. species
      numzones_nc=nj_d            ! number of zones
      btf_nc=DABS(btor_d)         ! tf (tesla) fld at cntr of outer flux
      drshaf_nc=0.D0              ! shaf shift of outer zone bdr. (cm)
      rminor_nc=rminavnpsi_d(nj_d)*1.D2
      rmajor_nc=rmajavnpsi_d(nj_d)*1.D2
c
        do j=1,ng_nc
         aplasm_nc(j)=amassgas_exp
        enddo
c
        do j=1,nj_d
          rhoel_nc(j)=dabs(ene_d(j)*1.D-6) !electron density (cm**-3)
          rhob_nc(j,1)=dabs(en_d(j,1)*1.D-6) !hyd. spec. den s (cm**-3)
          rhi_nc(j)=dabs(en_d(j,2)*1.D-6)  !z.c. array of av. impurity s density
          rhoi_nc(j)=rhi_nc(j)+rhob_nc(j,1) !z.c. array of total ion density
          te_nc(j)=dabs(te_d(j)*1.D3)      !z.c. array of Te (ev)
          ti_nc(j)=dabs(ti_d(j)*1.D3)      !z.c. array of Ti (ev)
          zeff_nc(j)=zeff_d(j)             !z.c. array of plasma zeff
          q_nc(j)=q_d(j)                   !z.c. array of safety factor
          aimp_nc(j)=pimpa                 !mass of impurity
          xzimp_nc(j)=pimpz                !charge of impurity
        enddo
c
        call kapisn(
     >          nkimod_nc,         !*kapai model nr desired
     >          aimp_nc,           !*atomic mass of av. impurity
     >          xzimp_nc,          !*atomic number of av. impurity
     >          aplasm_nc,         !array of atomic masses of hyd. species
     >          ng_nc,             !number of hyd. species
     >          nj,                !radial domain size
     >          rhoel_nc,          !zone centered electron density (cm**-3)
     >          rhob_nc,           !z.c. array of hyd. spec. den s (cm**-3)
     >          rhi_nc,            !z.c. array of av. impurity s density  
     >          rhoi_nc,           !z.c. array of total ion density       
     >          te_nc,             !z.c. array of Te (ev)
     >          ti_nc,             !z.c. array of Ti (ev)
     >          zeff_nc,           !z.c. array of plasma zeff
     >          numzones_nc,       !number of zones 
     >          q_nc,              !z.c. array of safety factor
     >          btf_nc,            !tf (tesla) fld at cntr of outer flux 
     >          drshaf_nc,         !shaf shift of outer zone bdr. (cm)
     >          rminor_nc,         !plasma minor radius (cm)
     >          rmajor_nc,         !major radius (cntr of outer flux) (cm)
     >          istringer_nc,      !Stringer correctioon
     >          xkapi_nc,          !o Neo-Class Ion thermal diff. (cm**2/sec)
     >          xnstari_nc,        !o nu-star-ions, 
     >          xnstare_nc,        !o nu-star-elecs, 
     >          ztaui_nc,          !o ion collision time
     >          zrhoi_nc,          !o ion poloidal gyro-radius
     >          zfluxlim_nc)       !o flux lim flow max temp grad length
c
       do j=1,nj_d
         xkapi_nc(j)=xkapi_nc(j)*1.D-4*en_d(j,1)/elongx_d(j) ! 1/(m*s)
c        write(*,'(i3,2x,0p1f4.2,1pe13.5)') j, rho_d(j), xkapi_nc(j)
         delta_nc(j)=rminavnpsi_d(j) / rmajor_d
         zrhoi_nc(j)=zrhoi_nc(j)*1.D-2 ! convert to meters
         zptineomax(j-1)=r_d(jmaxm+1)*100.D0/zfluxlim_nc(j)
         xkineo_d(j)=xkapi_nc(j)
       enddo
c
c... calculate conductivity and diffusivity
c
       if ((xkineo_d(nj_d).eq.0).or.(imyneoclass.eq.1)) then
         if (i_proc.eq.0) then
            write(6,*) 'Neoclassical xkineo calculated by kapisn'
         endif
         do j=1,nj_d
           xkineo_d(j)=xkapi_nc(j)
c           xkmneo_d(j)=0.1D0*delta_nc(j)**2.D0*
c     >                 (zrhoi_nc(j)**2.D0/ztaui_nc(j))
           xkmneo_d(j)=delta_nc(j)**1.5D0*xkineo_d(j)
c           write(*,60) j, rho(j), delta_nc(j), zrhoi_nc(j), ztaui_nc(j),
c     >                 xkmneo_d(j)/en_d(j,1)/grho2npsi_d(j),
c     >                 xkineo_d(j)/en_d(j,1)/grho2npsi_d(j)
         enddo
       endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
 20   format(' Using artificial chi for neoclassical transport')
 30   format(' Using KAPISN for neoclassical transport, ',
     >       'use_xneo_m =',i2)
 60   format(2x,i2,2x,0p1f10.6,1p6e15.7)
c
       end

