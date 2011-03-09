      subroutine neuden(nneu,kion,jmaxmm,jmaxm,rm,hcap,volfac,volsn,
     &                  sfarea,dn1,dn2,dnv,wn1,wn2,wnv,eirate,fluxn,
     &                  cx12r,en,enn,ennw,tn)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c NEUDEN calculates the normalizing factors that
c yield particle conservation under present
c plasma conditions.  The neutral density profiles
c are linearly proportional to the magnitude of their
c respective sources.
c After normalizing the single source neutral densities
c with the appropriate factors, the total neutral
c density and average neutral energy are computed.
c
c If a scrape-off layer is being modelled (njs>nj),
c only sources within the main plasma (j< = nj), where
c particle conservation is maintained, are
c used in calculating the normalization.
c
c Note: flxmod is the number of particles in error,
c       which adjusts the recycling rate in SOURCE
c     include 'param.i'
c     include 'geom.i'
c     include 'mesh.i'
c     include 'neut.i'
c     include 'numbrs.i'
c     include 'solcon.i'
c     include 'soln.i'
c     include 'sourc.i'
c     include 'verbose.i'
c
      integer nneu, kion, jmaxm, jmaxmm, in
      real*8 flxnum(2)
      real*8 volfac, sfarea
      real*8 sionvn, svoln, sionwn, fvn, f1w, 
     &          fluxn(2), flxmod(2)
      real*8 rm(jmaxmm), hcap(jmaxmm), cx12r(jmaxmm)
      real*8 dn1(jmaxmm,2), dn2(jmaxmm,2), 
     &          dnv(jmaxmm,2), wn1(jmaxmm,2), wn2(jmaxmm,2),
     &          wnv(jmaxmm,2), volsn(jmaxmm,2)
      real*8 en(jmaxmm,kion)
      real*8 ennv(jmaxmm,2), ennw(jmaxmm,2), scx(jmaxmm,2)
      real*8 eirate(jmaxmm), enn(jmaxmm,2), tn(jmaxmm,2)
c
c... some settings
c
      in = 1
      flxmod(1) = 0.
      flxmod(2) = 0.
      dtt = 0.
c
      if (nneu .eq. 2)  go to 2390
c
c          One neutral species present
c          Volume sources - scale dnv to make total ionization source
c          rate = total neutral volume source rate.
c     total electron impact ionization,#/sec,due to neutral density
c     established by volume sources is sionvn:
c
      call totsrc (jmaxm,rm,hcap,volfac,eirate,dnv(1,in), sionvn)
c
c     the volume source that established dnv is just the recombination
c     rate  and charge exchange with beam neutrals .
c     The total,integrated source,#/sec,is volfac*svoln:
      call trapv (rm,volsn(1,in),hcap,jmaxm,svoln)
c
      fvn = 0.0
      if (sionvn .ne. 0.0) fvn = volfac*svoln/sionvn
      do 2382 j=1,jmaxm
 2382 ennv(j,in) = fvn*dnv(j,in)
c
c          Wall sources - scale dn1 to make total ionization source
c          rate = total inward neutral flux at surface (which includes
c          a correction to attain particle conservation).
c
c     sionwn is the #/sec of electron impact ionization events due
c     to neutral density dn1 (dn1 was established by wall sources):
c
      call totsrc(jmaxm,rm,hcap,volfac,eirate,dn1(1,in),sionwn)
      flxnum(in) = sfarea*fluxn(in)
****  if ( dt .ne. 0.0)  flxnum(in) = flxnum(in) + flxmod(in)/dt
      if (dtt .ne. 0.0)
     .      flxnum(in) = flxnum(in) + flxmod(in)/dtt     ! HSJ 2/9/96
      f1w = flxnum(in)/sionwn
      do j=1,jmaxm
        ennw(j,in) = f1w*dn1(j,in)
        enn (j,in) = ennv(j,in) + ennw(j,in)
        tn  (j,in) = (ennw(j,in)*wn1(j,in)
     .              + ennv(j,in)*wnv(j,in))/enn(j,in)
      end do
      if (neucgvb .gt. 0) then
         do j=1,jmaxm
            if (ennw(j,in) .lt. 0)  write (6,'("NEUDEN:,j,ennw(j,in) =",
     .                                         i5,1pe16.4)')j,ennw(j,in)
            if (ennv(j,in) .lt. 0)  write (6,'("NEUDEN:,j,ennv(j,in) =",
     .                                         i5,1pe16.4)')j,ennv(j,in)
            jneg=jneg+1
         end do
         write (6, '("wall source normalization factor f1w =",
     .                             1pe12.4)')f1w
         write (6, '("volume source normalization factor fvn =",
     .                             1pe12.4)')fvn
      end if
c
      return
c
c          Two neutral species present
c          Volume sources - scale the dnv profiles so that the total
c          source rate of each ion = total neutral volume source rate
c          of that ion.  This requires solving for the scaling factors
c          f1v and f2v such that:
c               f1v * (sionv1+scxv1) - f2v * (scxv2)        = svoln1
c               f1v * (-scxv1)       + f2v * (sionv2+scxv2) = svoln2
c
 2390 call totsrc(jmaxm,rm,hcap,volfac,eirate,dnv(1,1),sionv1) ! total, #/sec ioniz. rate, species 1
      call totsrc(jmaxm,rm,hcap,volfac,eirate,dnv(1,2),sionv2) ! total, #/sec ioniz. rate, species 2
      do 2392 j=1,jmaxm
      scx(j,1) = en(j,2)*dnv(j,1)*cx12r(j)
 2392 scx(j,2) = en(j,1)*dnv(j,2)*cx12r(j)
      call trapv (r,scx(1,1),hcap,jmaxm,scxv1)
      scxv1 = volfac*scxv1
      call trapv (r,scx(1,2),hcap,jmaxm,scxv2)
      scxv2 = volfac*scxv2
      call trapv (r,volsn(1,1),hcap,jmaxm,svoln1)! due to recombination and
      svoln1 = volfac*svoln1                  ! fast ion charge exchange
      call trapv (r,volsn(1,2),hcap,jmaxm,svoln2)
      svoln2 = volfac*svoln2
      call solve2(sionv1+scxv1,-scxv2,svoln1,
     .            -scxv1,sionv2+scxv2,svoln2,f1v,f2v)
      do 2395 j=1,jmaxm
      ennv(j,1) = f1v*dnv(j,1)
 2395 ennv(j,2) = f2v*dnv(j,2)
      fvn = 0.5*(f1v+f2v)
c
c          Wall sources - scale dn1 and dn2 (adjust recycle flux)
c          so that for each species the total ion source rate =
c          total neutral flux at surface (which includes a correction
c          to attain particle conservation).  This involves scaling
c          fluxn(i) by the factors fiw such that:
c             f1w*(sion11+scxw1) + f2w*(sion12-scxw2) = flxnum(1)
c             f1w*(sion21-scxw1) + f2w*(sion22+scxw2) = flxnum(2)
c
c     recall that dn1(j,1) is density of species 1 due to edge flux
c     of species 1
c     recall that dn1(j,2) is density of species 1 due to edge flux
c     of species 2
c     and the same for species2,dn2(j,1),dn2(j,2)
c     below sionmn is total electron impact ionization rate,#/sec,
c     for neutral density  species m ,which was established by edge flux
c     of species n,(m=1,2,n=1,2)
c
      call totsrc(jmaxm,rm,hcap,volfac,eirate,dn1(1,1),sion11)
      call totsrc(jmaxm,rm,hcap,volfac,eirate,dn1(1,2),sion12)
      call totsrc(jmaxm,rm,hcap,volfac,eirate,dn2(1,1),sion21)
      call totsrc(jmaxm,rm,hcap,volfac,eirate,dn2(1,2),sion22)
c
c     now get the charge exchange rates due to these individually
c     sourced neutral densities.
c     scx(j,i) is the net gain of species i due to charge exchange
c     at spatial grid point j,due to species 1 and 2 neutral
c     densities sourced by edge flux i:
c
      do 2400 j=1,jmaxm
      scx(j,1) = (en(j,2)*dn1(j,1) - en(j,1)*dn2(j,1))*cx12r(j)
 2400 scx(j,2) = (en(j,1)*dn2(j,2) - en(j,2)*dn1(j,2))*cx12r(j)
      call trapv (r,scx(1,1),hcap,jmaxm,scxw1)
      scxw1 = volfac*scxw1
      call trapv (r,scx(1,2),hcap,jmaxm,scxw2)
      scxw2 = volfac*scxw2
      do i=1,2
        flxnum(i) = sfarea*fluxn(i)
****    if ( dt .ne. 0.0) flxnum(i) = flxnum(i) + flxmod(i)/dt
        if (dtt .ne. 0.0)
     .      flxnum(i) = flxnum(i) + flxmod(i)/dtt       ! HSJ 2/9/96
      end do
      call solve2 (sion11+scxw1,sion12-scxw2,flxnum(1),
     .             sion21-scxw1,sion22+scxw2,flxnum(2),f1w,f2w)
      do j=1,jmaxm
        scx(j,1) = 0.0
        scx(j,2) = 0.0
        ennw(j,1) = f1w*dn1(j,1) + f2w*dn1(j,2)
        ennw(j,2) = f1w*dn2(j,1) + f2w*dn2(j,2)
        enn(j,1) = ennw(j,1) + ennv(j,1)
        enn(j,2) = ennw(j,2) + ennv(j,2)
        tn(j,1) = (f1w*dn1(j,1)*wn1(j,1) + f2w*dn1(j,2)*wn1(j,2)
     .              + ennv(j,1)*wnv(j,1))/enn(j,1)
        tn(j,2) = (f1w*dn2(j,1)*wn2(j,1) + f2w*dn2(j,2)*wn2(j,2)
     .              + ennv(j,2)*wnv(j,2))/enn(j,2)
      end do
      if (neucgvb .gt. 0) then
         jneg=0
         do j=1,jmaxm
            if (ennw(j,1) .lt. 0)  write (6, '("NEUDEN,j,ennw(j,1) =",
     .                                         i5,1pe16.4)')j,ennw(j,1)
            if (ennv(j,1) .lt. 0)  write (6, '("NEUDEN,j,ennv(j,1) =",
     .                                         i5,1pe16.4)')j,ennv(j,1)
            if (ennw(j,2) .lt. 0)  write (6, '("NEUDEN,j,ennw(j,2) =",
     .                                         i5,1pe16.4)')j,ennw(j,2)
            if (ennv(j,1) .lt. 0)  write (6, '("NEUDEN,j,ennv(j,2) =",
     .                                         i5,1pe16.4)')j,ennv(j,2)
         end do
         write (6, '("neutral wall source normalization factor",
     .                           "  f1w,f2w =", 2(1x,1pe12.4))')f1w,f2w
         write (6, '("total number of neutrals required:" /
     .               "  species 1 and 2 flxnum(1),flxnum(2) =",
     .                  2(1x,1pe12.4))')flxnum(1),flxnum(2)
         write (6, '("sfarea,fluxn(1),flxmod(1),dtt =",4(x,1pe12.4))')
     .                sfarea,fluxn(1),flxmod(1),dtt
         write (6, '("sfarea,fluxn(2),flxmod(2),dtt =",4(x,1pe12.4))')
     .                sfarea,fluxn(2),flxmod(2),dtt
         write (6, '("sion11,scxw1 =",2(1x,1pe12.4))')sion11,scxw1
         write (6, '("sion12,scxw2 =",2(1x,1pe12.4))')sion12,scxw2
         write (6, '("sion21,scxw1 =",2(1x,1pe12.4))')sion21,scxw1
         write (6, '("sion22,scxw2 =",2(1x,1pe12.4))')sion22,scxw2
         write (6, '(/// "neutral volume source normalization factor",
     .                   "  f1v,f2v =", 2(1x,1pe12.4))')f1v,f2v
         write (6, '("volume source of neutrals,species 1 and 2:",
     .                2(1x,1pe12.4))')svoln1,svoln2
        write (6, '("sionv1 =",1x,1pe12.4)')sionv1
        write (6, '("sionv2 =",1x,1pe12.4)')sionv2
        write (6, '("scxv1  =",1x,1pe12.4)')scxv1
        write (6, '("scxv2  =",1x,1pe12.4)')scxv2
      end if
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine totsrc (jmaxm, rm, hcap, volfac, eirate, dn, stot)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     TOTSRC computes the total rate at which the given neutral
c     density profile is producing ions via electron ionization (1/sec)
c
c     include 'param.i'
c     include 'geom.i'
c     include 'mesh.i'
c     include 'neut.i'
c     include 'numbrs.i'
c
      dimension  dn(*), siont(jmaxm), eirate(*), rm(*), hcap(*)
c
      do j=1,jmaxm
        siont(j) = dn(j)*eirate(j)
      end do
      call trapv (rm, siont, hcap, jmaxm, stot)
      stot = volfac * stot
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine solve2 (a11, a12, c1, a21, a22, c2, f1, f2)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     SOLVE2 solves a set of 2 linear equations in 2 unknowns:
c
c       a11*f1 + a12*f2 = c1
c       a21*f1 + a22*f2 = c2
c
c     Zero determinants are not tolerated.
c
c     include 'param.i'
c     include 'io.i'
c
      det = a11*a22 - a12*a21
      if (det .eq. 0.0)  go to 100
      f1  = (a22*c1 - a12*c2)/det
      f2  = (a11*c2 - a21*c1)/det
      return
  100 write  (*, 8000)
 8000 format (/ ' ERROR (nneu = 2 neutral transport): det=0.'     /
     .          ' Probably due to zero ion or neutral densities.' /
     .          ' subroutine SOLVE2: determinant is zero.' /)
c     call STOP ('subroutine SOLVE2: determinant is zero', 44)
      return
c
      end
