c@tenspline.f
c jek 18-Jan-11
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Smoothes data using CURVSS tension spline routine
c      ne : ismoo_ne
c      ni : ismoo_ni
c      nf : ismoo_nf
c      zeff : ismoo_zeff
c      q : ismoo_q
c      qnb : ismoo_qnb
c      vrot : ismoo_vrot
c      delta  : ismoo_delta
c      grho1,grho2 : ismoo_grho 
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c 
      subroutine tenspline
c
      implicit none
      include 'mpif.h'
      include '../inc/tport.m'
      include '../inc/data.m'
      include '../inc/glf.m'
c
      integer j, ier
      real*8 eps, sigmas
      real*8 wt(nj), ys(nj), ysp(nj), td(nj), tsd1(nj), hd(nj),
     &       hsd1(nj), hsd2(nj), rd(nj), rsd1(nj), rsd2(nj), v(nj)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      do j=1,nj_d
        rho_d(j)=rho(j-1)
        wt(j)=1.D0
      enddo
c
      eps=1.D-6
      sigmas=0.D0
      if (ismoo_ne.gt.0) then
        if (i_proc.eq.0) write(*,85) ismoo_ne
        wt(1)=1.D-3     ! do not smooth central and edge values
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,ene_d*1.D-19,wt,0,ismoo_ne,eps,ys,ysp,
     >              sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j), 
     >                          ene_d(j), ys(j)*1.D19
          ene_d(j)=ys(j)*1.D19
        enddo
      endif
c
      if (ismoo_ni.gt.0) then
        if (i_proc.eq.0) write(*,86) ismoo_ni
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,en_d(1:nj_d,1)*1.D-19,wt,0,ismoo_ni,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j), 
     >                          en_d(j,1), ys(j)*1.D19
          en_d(j,1)=ys(j)*1.D19
        enddo
      endif
c
      if (ismoo_nf.gt.0) then
        if (i_proc.eq.0) write(*,90) ismoo_nf
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,enbeam_d(1:nj_d)*1.D-19,wt,0,ismoo_nf,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j), 
     >                          enbeam_d(j), ys(j)*1.D19
          enbeam_d(j)=ys(j)*1.D19
        enddo
      endif
c
      if (ismoo_zeff.gt.0) then
        if (i_proc.eq.0) write(*,87) ismoo_zeff
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,zeff_d,wt,0,ismoo_zeff,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j),
     >                          zeff_d(j), ys(j)
          zeff_d(j)=ys(j)
        enddo
      endif
c
      if (ismoo_q.gt.0) then
        if (i_proc.eq.0) write(*,88) ismoo_q
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,q_d,wt,0,ismoo_q,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j),
     >                          q_d(j), ys(j)
          q_d(j)=ys(j)
        enddo
      endif
c
      if (ismoo_qnb.gt.0) then
        if (i_proc.eq.0) write(*,95) ismoo_qnb
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,qbeame_d*1.D-6,wt,0,ismoo_qnb,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j),
     >                          qbeame_d(j), ys(j)*1.D6
          qbeame_d(j)=ys(j)*1.D6
        enddo
        call curvss(nj_d,rho_d,qbeami_d*1.D-6,wt,0,ismoo_qnb,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j),
     >                          qbeami_d(j), ys(j)*1.D6
          qbeami_d(j)=ys(j)*1.D6
        enddo
        call curvss(nj_d,rho_d,sbeam_d*1.D-19,wt,0,ismoo_qnb,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j),
     >                          sbeam_d(j), ys(j)*1.D19
          sbeam_d(j)=ys(j)*1.D19
        enddo
      endif
c
      eps=1.D-6
      sigmas=0.D0
c
      if (ismoo_vrot.gt.0) then
        if (i_proc.eq.0) write(*,89) ismoo_vrot
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,angrot_d*1.D-4,wt,0,ismoo_vrot,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j),
     >                          angrot_d(j), ys(j)*1.D4
          angrot_d(j)=ys(j)*1.D4
        enddo
      endif
c
      if (ismoo_delta.gt.0) then
        if (i_proc.eq.0) write(*,91) ismoo_delta
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,deltax_d,wt,0,ismoo_delta,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j),
     >                          deltax_d(j), ys(j)
          deltax_d(j)=ys(j)
        enddo
      endif
c
      if (ismoo_grho.gt.0) then
        if (i_proc.eq.0) write(*,91) ismoo_grho
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,grho1npsi_d,wt,0,ismoo_grho,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j),
     >                          grho1npsi_d(j), ys(j)
          grho1npsi_d(j)=ys(j)
        enddo
        wt(1)=1.D-3     ! do not smooth central and edge values  
        wt(nj_d)=1.D-3
        call curvss(nj_d,rho_d,grho2npsi_d,wt,0,ismoo_grho,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j),
     >                          grho2npsi_d(j), ys(j)
          grho2npsi_d(j)=ys(j)
        enddo
        do j=1,nj_d
          wt(j)=0.02D0  ! don't smooth rmin as much
        enddo
        wt(1)=1.D-4     ! do not smooth central and edge values  
        wt(nj_d)=1.D-4
        call curvss(nj_d,rho_d,rminavnpsi_d,wt,0,ismoo_grho,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j),
     >                          rminavnpsi_d(j), ys(j)
          rminavnpsi_d(j)=ys(j)
          wt(j)=1.D0     ! reset weights if smoothing later
        enddo
        do j=1,nj_d
          wt(j)=0.02D0  ! don't smooth rmaj as much
        enddo
        wt(1)=1.D-4     ! do not smooth central and edge values  
        wt(nj_d)=1.D-4
        call curvss(nj_d,rho_d,rmajavnpsi_d,wt,0,ismoo_grho,
     >       eps,ys,ysp,sigmas,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ier)
        if (ier.gt.0 .and. i_proc.eq.0) write(*,*) 'ier = ',ier
        do j=1,nj_d
          if (lprint_smoo.gt.0) write(*,50) j, rho_d(j),
     >                          rmajavnpsi_d(j), ys(j)
          rmajavnpsi_d(j)=ys(j)
          wt(j)=1.D0     ! reset weights if smoothing later
        enddo
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
 50   format(i2,2x,0p1f4.2,0p6e13.5)
 85   format(' Smoothing electron density profile ...',0p1f4.2)
 86   format(' Smoothing ion density profile ...',0p1f4.2)
 87   format(' Smoothing Zeff profile ...',0p1f4.2)
 88   format(' Smoothing q-profile ...',0p1f4.2)
 89   format(' Smoothing angrot profile ...',0p1f4.2)
 90   format(' Smoothing fast ion density profile ...',0p1f4.2)
 91   format(' Smoothing triangularity profile ...',0p1f4.2)
 95   format(' Smoothing Qnb profile ...',0p1f4.2)
c
       end
