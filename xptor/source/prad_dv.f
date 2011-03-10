cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prad_dv
c
c same as prad only using dv method variables
c
c Calculates the Brehmsstrahlung radiation
c from the NRL plasma formulary and the synchrotron
c radiation (Trubnikov, JETP Letters 16, 25 (1972).
c
      implicit none
c
      include '../inc/ptor.m'
      include '../inc/tport.m'
      include '../inc/data.m'
      include '../inc/glf.m'
c
      real*8 volint  ! volume integrating function
      real*8 zpi, charge, xmasse, cee
      real*8 wbx, teerg, wpsq, chi, phi, vdotsq, ene, cyclo
      integer k
c
      cee = 2.99792458D10    ! speed of light (cm/s)
      charge = 4.803207D-10  ! electron charge (Stat-Coul)
      xmasse = 9.10939e-28   ! electron mass (gm)
      zpi = atan2 ( 0.0D0, -1.0D0 )
c
c Brehmsstrahlung 
c
      do k=1,ngrid
         Pradb(k)=1.69D-32*1.D26*ne_m(k)**2*zeff_exp(k)
     &            *sqrt(1.D3*te_m(k))
c         write(*,*) k, Pradb(k)
      enddo
c
      Pradb_tot=volint(Pradb)
c
c Synchrotron radiation 
c
      if (radref.gt.0.) then
        wbx = charge*abs(btor_d*1.D4) / (xmasse*cee)
        do k=1,ngrid
           teerg = te_m(k)*1.6D-9
           ene = ne_m(k)*1.e13
           wpsq = 4.D0*zpi*ene*charge**2.D0/xmasse
           chi = (r_d(nj_d)*100.D0)/(rmajor_d*1.D2)*
     &           dsqrt(xmasse*cee**2.D0/teerg)
           phi = 60.D0*(teerg/(xmasse*cee**2.D0))**1.5D0*
     &           dsqrt(cee*wbx/(r_d(nj_d)*100.D0*wpsq)*
     &           radref*(1.D0+chi))
           vdotsq = 2.D0*wbx**2.D0*teerg/xmasse
           cyclo = ene*1.5D0*charge**2.D0/cee**3.D0*
     &             vdotsq*phi/1.6D-9
           Prads(k)=cyclo*1.6022D-10*1.D-6
c          write(*,100) k, r_d(k+1)*100., teerg, ene, cyclo
        enddo
c
      Prads_tot=volint(Prads)
      endif
c      if(i_proc.eq.0) write(*,*) 'Computed Prad = ',Pradb_tot
c
 100  format(i2,2x,0p1f8.4,1p6e13.5)
      return
      end

