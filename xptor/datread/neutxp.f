c@neutxp.f
c jek 16-feb-10
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c... Calls ONETWO neutrals package using experimental profiles
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  
      subroutine neutxp
c
      implicit none
c
      include '../inc/data.m'
      include '../inc/tport.m'
      include '../inc/input.m'
      include '../inc/model.m'
      include '../inc/ptor.m'
c
      integer j
      real*8 distr, drm, dvoldr_p, dvoldr_m, nim_p, nim_m, taupm
      real*8 srecom_dm(jmaxmm,2), sion_dm(jmaxmm,2), 
     >       qione_dm(jmaxmm), qioni_dm(jmaxmm), qcx_dm(jmaxmm)
c
      write(*,*) 'Calling NEUT neutrals package ...'
      pi_m=atan2(0.0D0,-1.0D0)
c
      distr=0.D0
      do j=2,nj_d
        drm=r_d(j)-r_d(j-1)
        dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
        dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)
        nim_p=1.D-19*en_d(j,1)
        nim_m=1.D-19*en_d(j-1,1)
        distr=distr+kevdsecpmw*1.D19*0.5D0*(dvoldr_p*nim_p+
     >        dvoldr_m*nim_m)*drm
        taupi_m(j)=distr/(flow_exp(j-1)+1.D-10)
c        write(*,100) j, rho(j-1), nim_m, distr, taupi_m(j)
      enddo
      taupm=taupi_m(nj_d)
      write(*,'(a7,2x,1pe13.5)') 'taup = ',taupm
c
      call neut(taupm,srecom_dm,sion_dm,qione_dm,qioni_dm,qcx_dm)
c
c      do j=1,nj_d-1
c        write(*,100) j, rho(j), powe_ion_exp(j), powi_ion_exp(j),
c     >               powi_cx_exp(j)
c      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
 100  format(i2,2x,0p1f4.2,1p8e13.5)
c
      end
