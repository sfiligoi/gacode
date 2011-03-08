c@pohmic.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c     Calculates the Ohmic Heating power
c     using either Neoclassical or Spitzer resistivity and
c     current density profile
c        iohm=1 Neoclassical resistivity (ohm-cm)
c        iohm=2 Spitzer resistivity (ohm-cm)
c     For the neoclassical resistivity, we use
c        Hirshman, et al., Nucl. Fusion 17, 3 (1977)
c        Hirshman, et al., Nucl. Fusion 21, 9 (1981)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine pohmic
c
      implicit none
c
      include '../inc/ptor.m'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
c	   
      integer k
      real*8 volint       ! volume integrating function
      real*8 cur(mxgrid)  ! current density profile (amps/cm**2)
      real*8 eta(mxgrid)  ! resistivity profile (ohm-cm)
      real*8 ipsum, ip_tot, lnlam, ds
      real*8 delt, xft, xmasse, charge, charg4, vthe, xlam,
     &       ene, xnue, xnuse, xnusem, zf, zr, cr, cra,
     &       zeta, xkap11, ftad, fftrap, etasp, etapar
c
      data xmasse / 9.1093897D-28 /   ! electron mass (gm)
      data charge / 4.8032067D-10 /   ! electron charge (esu)
      data charg4 / 5.3226161D-38 /   ! charge**4 (esu**4)
c
      if (iohm.eq.1) then
c
c...Neoclassical resisitivity
c
      do k=1,ngrid
         ene=denpro(k)*1.D14
         delt=r(k,1)*100.D0/(dsqrt(kappa_d)*rmajor_d*100.D0)
         xft = 1.D0 - (1.D0 - delt)**2.D0 / dsqrt( 1.D0 - delt**2.D0) /
     &         (1.D0 + 1.46D0*dsqrt(delt))
         vthe = dsqrt(2.D0*tepro(k)*1.6D-9/xmasse)
         xlam=24.D0-dlog( dsqrt(ene)/(1.D3*tepro(k)) )
         xnue=(4.D0/3.D0)*dsqrt(2.D0*pi)*ene*zeffpro(k)*
     &        charg4*xlam / dsqrt(xmasse*(tepro(k)*1.6D-9)**3.D0)
         xnuse=dsqrt(2.D0)*rmajor_d*100.D0*dabs(q_exp(k))*xnue/
     &         (delt**(3./2.)*vthe)
c
         zf=zeffpro(k)
         xnusem=xnuse/zf
         zr=((0.222D0*zf+1.198D0)*zf+1.D0) / 
     &      ((0.753D0*zf+2.966D0)*zf+1.D0)
         cr=0.56D0*(3.D0-zf)/(3.D0+zf)/zf
         zeta=0.58D0+0.20D0*zf
         xkap11=1.D0/zr
         etasp=xmasse*xnue/(ene*charge**2.D0)/xkap11  ! Spitzer resistivity (s)
         if (zf.lt.3.D0) then
           cra=cr
         else
           cra=0.D0
         endif
         ftad=xft/(1.D0+zeta*xnusem)
         fftrap=ftad*(1.D0+cra*(1.D0-ftad))
         etapar=etasp / (1.D0-fftrap)   ! parallel resisitivity
c        write(*,100) k, r(k,1)*100.D0, etasp, etapar
         eta(k)=etapar * 9.D9 * 100.D0  ! secs -> ohm-cm
         cur(k)=1.D0/eta(k)
      enddo
c
      else
c
c Spitzer resistivity 
c NRL formulary (ohm-cm)
c Note : factor of 0.5 for parallel vs. perp resisitivity
c
      do k=1,ngrid
         lnlam=24.D0-dlog(sqrt(1.D14*denpro(k))/1.D3/Tepro(k))
         lnlam=dmax1(lnlam,1.D0)
         eta(k) = 0.5D0*1.03D-2*Zeffpro(k)*lnlam/(1000.D0*tepro(k))**1.5D0
c         write(*,100) k, rho(k), eta(k)/1.e2/9.e9  ! eta (secs)
         cur(k) = 1.D0/eta(k)
      enddo
c
      endif
c
c D3D formula to relate Ip (in MegAmps) to q_95, etc., (obtained from
c Jim Deboo in Feb. 1995., who says it usually agrees with efit to within
c 5-10%):
c
      ip_tot=(5.D0*amin**2*Btor/Rmaj/qa)*(1.D0+1.5D0*(amin/Rmaj)**2)
     &     *(1.D0+kappa**2)/2.D0
c
c calculate area-integrated current to normalize the current-density:
c
      ipsum=0.D0
      do k=1,ngrid
         ds=vprime(k,1)/(2.D0*pi*Rmaj)*dr(k,1)
         ipsum=ipsum+ds*cur(k)
      enddo
c
c 1.e6 convert Iptot from MA to A, 1.e4 to convert cur from 
c amps/m**2 to amps/cm**2:
c
      do k=1,ngrid
         cur(k)=cur(k)*ip_tot*1.D6/ipsum/1.D4
      enddo
c
      if ((iohm.eq.2).and.(epulse_pt.ne.0)) then
        if(time.lt.time1_pt) then
          do k=1,ngrid
            curden_exp(k)=cur(k)
          enddo
        else
          do k=1,ngrid
            cur(k)=curden_exp(k)
          enddo
        endif
      endif
c
c Ohmic power density in MW/m**3
c eta (ohm*cm), cur (A/cm**2)
c Note: qohm_d (W/m**3), curden_d (A/m**2) from data
c
      do k=1,ngrid
         Pohpro(k)=eta(k)*cur(k)**2
c        write(*,100) k, r(k,1)*100.D0, eta(k), qohm_d(k),
c    &                eta(k)*cur(k)**2,
c    &                eta(k)*(curden_d(k)/1.D4)**2
      enddo
c
      Pohmic_tot=volint(Pohpro)
c
 100  format(i4,2x,0pf10.5,1p6e12.4)
c
      return
      end
