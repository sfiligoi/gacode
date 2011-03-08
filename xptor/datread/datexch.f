c@datexch.f
c jek 20-Jan-11
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c... Compute collisional exchange using _d data
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  
      subroutine datexch
c
      implicit none
      include '../inc/data.m'
      include '../inc/tport.m'
c
      integer j
      real*8 zbrac, alamda, tau_e
c
c NRL formulary, watts/meter**3

       if(qdelt_d(nj_d).eq.0.or.iexp_exch.ge.1) then
         do j=1,nj_d
           zbrac=en_d(j,1)/ene_d(j)+apgasa/
     >           apimpa*en_d(j,2)/ene_d(j)*apimpz**2
           alamda=24.D0-dlog(dsqrt(ene_d(j)*1.D-6)/(te_d(j)*1.D3))
           tau_e=3.44D5*(te_d(j)*1.D3)**1.5D0/(ene_d(j)*1.D-6)/alamda
           qdelt_d(j)=-10.D0**6*1.6022D-19*zbrac*3.D0/(1836.D0*apgasa)/
     >                 tau_e*ene_d(j)*1.D-6*(te_d(j)-ti_d(j))*1.D3
         enddo
       endif
c
       if(iexp_exch.eq.-1) then
        do j=1,nj_d
         qdelt_d(j)=0.D0
        enddo
       endif
c
      end
