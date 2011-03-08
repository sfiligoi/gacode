c@etaw17a.f
c jek 15-sep-98 initialized zgm matrix elements to zero
c               added subroutine abortb at end for incorporation into MLT
c pis 12-jun-98 removed  if ( zfns .lt. zepsmach) zgns = 0
c pis  5-may-98 removed some uneccessary loops and unused variables.
c pis  5-may-98 relabelled routine etaw17a.tex and introduced tomsqz instead of
c               tomslz. QZ algorithm more stable and exact than LZ.
c pis 29-apr-98 added wexb shearing rate to argument list, currently this
c               is implemented in two different ways defined by the value of
c                   letain(15) = 0 (default) matrices redifened
c                   letain(15) = 1 shearing rate subtracted from growth rates
c                                  (this implementationwill produce err mess.)
c pis 28-apr-98 replaced nag library routines with tomsqz (wrapper for lzhes
c               and lzit)
c jek 30-nov-97 added modified FLR parameter to include elongation
c               coeff in Alfven frequency is cetain(25) (default=0.D0)
c rgb 17-may-97 added if ( zfns .lt. zepsmach) zgns = 0
c   multiplied zamr(11,9) by zimp
c jek 12-apr-96 added cross term to electromagnetic version with
c               parallel ion motion
c jek 25-mar-96 added parallel ion motion in strong ballooning limit
c               to nine equation model (H=S/2q). cetain(12) (default=0.D0)
c               added to turn on/off
c rgb 13-feb-96 input gnein=-R(d n_e / d r)/n_e rather than gnsin
c   rearranged order of argument list to put gnein before gnhin
c rgb 12-feb-96 zgp* = zgp* + zgn*  -->  zgp* = zgt* + zgn*
c   (this should have no effect on the results)
c   corrected header for eleven eqns
c rgb 08-jan-96 diagnostic printout of frequencies, fluxes, and phases
c rgb 03-jul-95 computed zerrmax and used it for error control
c   implemented control of imatrx by letain(2)
c rgb 02-jul-95 added impurity parallel ion motion 11 equations
c rgb 01-jul-95 when computing fluxes, loop over all eigenvalues ieq
c rgb 29-jun-95 protected zalp, zalf, and kps when zgnh+zgth.lt.zepsqrt
c   corrected zamr(3,9), zamr(3,10), zbmr(3,10)
c rgb 28-jun-95 cleaned up complex matrix option
c rgb 24-jun-95 converted from eps*in to normalized gradients g*in
c   use zgnh instead of zgne in definition of k1 and k2
c   removed vef = 0.10
c   zamr(7,1) = zhalf*zgte - zone  inserted in the 7 eqn model
c   replaced 1./kpc with zanorm = cetain(11)
c   added the following equation to the 10 eqn set
c   zamr(9,8) = kps * ( zone - zfnz - zfns ) / ( zone - zft )
c   Pass nout through the argument list just after neq
c 04-apr-95 added vef, ion motion, finite beta
c 04-jan-95 added vef to argument list
c 06-dec-94 include disp9 from Weiland
c
c  THIS ROUTINE IS A MODIFICATION OF THE LINEAR PART OF ETAWN6
c  WRITTEN BY GLENN BATEMAN. INCLUDING COLLISIONS TO THE
c  TRAPPED ELECTRONS YIELDS A 7 EQUATION SYSTEM. WHEN PARALLEL
c  ION MOTION IS INCLUDED THERE ARE 8 EQUATIONS. THE 9 EQUATION
c  SYSTEM INCLUDES BOTH COLLISIONS AND PARALLEL ION MOTION.
c  WITH PARALLEL ION MOTION IN THE STRONG BALLOONING LIMIT.
c  THE 9 EQUATION SYSTEM IS MODIFIED TO INCLUDE SHEAR EFFECTS.
c  THE NEW ELECTROMAGNETIC SYSTEM INCORPORATES THE PARALLEL
c  VECTOR POTENTIAL AND ITS TIME DERIVATIVE AS NEW VARIABLES.
c  THIS LEADS TO A SYSTEM OF 10 EQUATIONS WHEN COLLISIONS AND
c  PARALLEL ION MOTION ARE INCLUDED. WITH PARALLEL ION MOTION
c  OF IMPURITIES, THE SYSTEM IS EXTENDED TO 11 EQUATIONS.
c  MOST PARAMETERS ARE TRANSFERRED THROUGH COMMON BLOCKS LIKE GRAD,
c  IMP AND ETAWN6. NOTE THE INVERSE DEFINITION OF ZTAUH AND ZTAUZ !
c
c   The variables are ordered as: e phi/Te, Th, Nh, Te, Nz, Tz, F, Vp,
c   Av, K where F is due to collisions, Vp is due to parallel ion motion,
c   Av is the vector potential and K is its time derivative.
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine etaw17a (letain, cetain, lprintin, neq, nout
     & , gnein, gnhin, gnzin
     & , gtein, gthin, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, betahin, betazin, ftrapein
     & , vef, q, shear, kappa, ekyrhoin, ekparlin, wexb
     & , ndim, omega, gamma, difthi, velthi, chieff
     & , nmodes, perform, nerr )
c
c..see the table of argument list variable names and definitions
c    at the end of this document
c
c  gnein  = - R ( d n_e / d r ) / n_e
c  gnhin  = - R ( d n_H / d r ) / n_H
c  gnzin  = - R ( d n_Z / d r ) / n_Z
c  gtein  = - R ( d T_e / d r ) / T_e
c  gthin  = - R ( d T_H / d r ) / T_H
c  gtzin  = - R ( d T_Z / d r ) / T_Z
c  tauhin = T_H / T_e
c  tauzin = T_Z / T_e
c
c  This version of etaw17a is intended for use on workstations
c
c  Note that two external subroutines need to be provided:
c
c  abortb (nout,text)  Prints line of text on output unit nout and
c                      then stops the run (it calls abort or stop)
c
c  tomsqz   Code wrapper for ACM/TOMS routine 535 implementing the QZ algorithm
c           complex generailized eigenvalue problem ( eigenvalues and eigenvecto
c           The algorithm itself consists of the three routines
c           cqzhes, cqzval and cqzvec
c
c  Compile this routine  and routines that call it with a compiler option
c  such as -r8  to convert real to double precision when used on workstations.
c
      implicit none
      include '../inc/glf.m'
c
      integer idp, ndim
      parameter ( idp = 15 )
c
      logical inital
      data inital /.true./
c
      dimension letain(32), cetain(32)
     &  , omega(*), gamma(*), difthi(ndim,*), velthi(*)
     &  , chieff(*), perform(*)
c
      real*8 cetain
     & , gnein, gnhin, gnzin
     & , gtein, gthin, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, betahin, betazin, ftrapein
     & , vef, q, shear, kappa, H, ekyrhoin, ekparlin, wexb
     & , omega, gamma, difthi, velthi, chieff, perform
c
      integer letain, lprintin, lprint, neq, nmodes, nerr, nout
     & , ieq, idim, imatrx, j1, j2, j, jd
c
c ndim  = first dimension of the 2-D array difthi
c         and the maximum number of unstable modes allowed
c nmodes = number of unstable modes
c ieq   = number of equations
c
c imatrx = the number of elements computed along each row and
c column of the transport matrix
c
c  Note:  Normally the transport matrix is  ieq-1 by ieq-1
c    however, if there are 6 equations, compute only a 4 by 4 matrix
c    since the impurity temperature equation is not used by most
c    transport codes including BALDUR
c
      real*8 zamr(idp,idp), zami(idp,idp), zbmr(idp,idp), zbmi(idp,idp)
      real*8 zamrt(idp,idp),zamit(idp,idp),zbmrt(idp,idp),zbmit(idp,idp)
      real*8 zvr(idp,idp), zvi(idp,idp), zomega(idp), zgamma(idp)
      real*8 zalfr(idp), zalfi(idp), zbeta(idp)
c
      integer ifail
c
c ( zamr(i,j), zami(i,j) ) = matrix A
c ( zbmr(i,j), zbmi(i,j) ) = matrix B
c
c Note that the eigenvalues are
c   zomega(j) = zalfr(j) / zbeta(j)
c   zgamma(j) = zalfi(j) / zbeta(j) - i wexb *letain(15)
c and the eigenvectors are
c   zevec(j) = cmplx( zvr(j), zvi(j) )
c Here, zbeta(j) will be 0.0 in the case of an infinite eigenvalue
c
      complex  zevec(idp,idp)
c
      real*8  zepsmach, zepsqrt
     & , zone, ztwo, zthree, zfour, zfive, zhalf, zquarter, ztvr, zftr
     & , ztwohlf
     & , zgne, zgnh, zgnz, zgns, zgte, zgth, zgtz
     & , ztauh, ztauz, zft
     & , zimp, zfnz, zmass, zfns, zflh, zflz, zgamax
     & , zetae, zetah
c
c zepsmach = machine epsilon
c zepsqrt  = sqrt ( machine epsilon )
c zone     = 1.0
c ztwo     = 2.0
c zthree   = 3.0
c zfour    = 4.0
c zfive    = 5.0
c zhalf    = 0.5
c zquarter = 0.25
c ztwohalf = 2.5
c ztvr     = 2.0 / 3.0
c zftr     = 5.0 / 3.0
c
      real*8 zflxph, zflxnh, zflxpe, zflxnz, zflxpz
     &  ,  zphsph, zphsnh, zphspe, zphsnz, zphspz
     &  ,  ztemp1, ztemp2, zreal, zimag
     &  ,  ztempa(idp), ztempb(idp), zerreal, zerimag, zerrmax
c
c  These are the effective thermal and particle diffusivities
c  computed in different ways for comparison.
c
      real*8 zflxm(idp,idp), zchim(idp,idp), zgm(idp,idp), zdg
c
c  zflxm(j,jd) are the normalized transport fluxes of:
c  zflxm(1,jd)   n_H T_H
c  zflxm(2,jd)   n_H
c  zflxm(3,jd)   n_e T_e
c  zflxm(4,jd)   n_Z
c  zflxm(5,jd)   n_Z T_Z
c
c  zgm(j,jd)   is a matrix of gradient scale lengths
c  zgm(1,jd) = zgth
c  zgm(2,jd) = zgnh
c  zgm(3,jd) = zgte
c  zgm(4,jd) = zgnz
c
c  Here jd=1 uses the input values zgph, zgnh, zgpe, zgnz
c  jd=2 uses zgph + zdg, zgnh, zgpe, zgnz
c  jd=3 uses zgph, zgnh + zdg, zgpe, zgnz
c  jd=4 uses zgph, zgnh, zgpe + zdg, zgnz
c  jd=5 uses zgph, zgnh, zgpe, zgnz + zdg
c  That is, the gradients are imcemented one at a time.
c
c..local variables added 4-jan-95
c
      real*8 bt, bt1,  zeni, k1, k2
      real*8 zkpsh, zkpsz, zalp, zalf, zanorm, zbetae, zrav
c
c
      save idim, zepsmach, zepsqrt
     &  , zone, ztwo, zthree, zfour, zfive, zhalf, zquarter, ztvr, zftr
     &  , ztwohlf
     &  , inital
c
c..initialize variables
c
      lprint = lprintin
c
      if ( nout .lt. 1  .or.  nout .gt. 99 ) nout = 6
c
      ieq = max ( 2, neq )
c
      vef = cetain(15) * vef
      imatrx = min ( ieq - 1, 4 )
      if ( letain(2) .gt. 0 ) imatrx = min ( imatrx, letain(2)-1 )
c
      if ( inital ) then
c
        idim = idp
c
        zone   = 1.D0
        ztwo   = 2.D0
        zthree = 3.D0
        zfour  = 4.D0
        zfive  = 5.D0
c
        zhalf = zone / ztwo
        ztwohlf  = ztwo + zhalf
        zquarter = zone / zfour
        ztvr = ztwo  / zthree
        zftr = zfive / zthree
c
        zepsmach = zhalf
  2     if ( zhalf * zepsmach + zone .gt. zone ) then
          zepsmach = zhalf * zepsmach
          go to 2
        endif
c
        zepsqrt = sqrt ( zepsmach )
c
c..Print Header
c
        if (i_proc .eq. 0) then
        write (nout,*)
        write (nout,*)
     &   'Weiland-Nordman eigenvalue equations, subroutine etaw17a'
        write (nout,*)
c
        if ( letain(6) .gt. 1 ) then
          write (nout,*)
     &     ' Eigenvalues computed using IMSL routine gvcrg'
        else
          write (nout,*)
     &     ' Eigenvalues computed using ACM/TOMS routine 535'
        endif
c
        if ( ieq .eq. 11) then
          write (nout,*)
          write (nout,*)
     &     'Eleven eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z,
     &     F, Vp_H, Av, K, and Vp_Z'
        elseif ( ieq .eq. 10) then
          write (nout,*)
          write (nout,*)
     &     'Ten eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z,
     &     F, Vp, Av, and K'
        elseif ( ieq .eq. 9) then
          write (nout,*)
          write (nout,*)
     &     'Nine eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z,
     &     F, Av, K, and Vp in strong ballooning limit'
        elseif ( ieq .eq. 8) then
          write (nout,*)
          write (nout,*)
     &     'Eight eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, and Vp'
        elseif ( ieq .eq. 7) then
          write (nout,*)
          write (nout,*)
     &     'Seven eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z, and F'
        elseif ( ieq .eq. 6 ) then
          write (nout,*)
          write (nout,*)
     &     'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
	       elseif ( ieq .eq. 5 ) then
          write (nout,*)
          write (nout,*)
     &     ' Five eqns for e phi/T_e, T_H, n_H, T_e, and Vp'
        elseif ( ieq .eq. 4 ) then
          write (nout,*)
          write (nout,*) ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
        elseif ( ieq .eq. 2 ) then
           write (nout,*)
           write (nout,*) ' Two eqns for e phi/T_e and T_H'
        else
          write (nout,*)
          write (nout,*) ieq,' = ieq in sbrtn etaw17a'
          call abortb(nout
     &     ,'the value of ieq is wrong in sbrtn etaw17a')
        endif
        endif
c
        inital = .false.
c
      endif
c
c..end of initialization
c
      nerr = 0
c
c..print header
c
      if ( lprint .gt. 0 ) then
c
        write (nout,*)
        write (nout,*)
     & 'Weiland-Nordman eigenvalue equations, subroutine etaw17a'
        write (nout,*) '(all frequencies normalized by omega_{De})'
        write (nout,*) '(all diffusivities normalized by '
     &    ,'omega_{De} / k_y^2'
c
        write (nout,108) 'gnein', gnein
        write (nout,108) 'gnhin', gnhin
        write (nout,108) 'gnzin', gnzin
        write (nout,108) 'gtein', gtein
        write (nout,108) 'gthin', gthin
        write (nout,108) 'gtzin', gtzin
        write (nout,108) 'tauhin', tauhin
        write (nout,108) 'tauzin', tauzin
        write (nout,108) 'ftrapein', ftrapein
        write (nout,108) 'fnzin', fnzin
        write (nout,108) 'czin', czin
        write (nout,108) 'azin', azin
        write (nout,108) 'fnsin', fnsin
        write (nout,108) 'betaein', betaein
        write (nout,108) 'betahin', betahin
        write (nout,108) 'betazin', betazin
        write (nout,108) 'vefin', vef
        write (nout,108) 'qin', q
        write (nout,108) 'shearin', shear
        write (nout,108) 'kappain', kappa
        write (nout,108) 'ekyrhoin', ekyrhoin
        write (nout,108) 'ekparlin', ekparlin
        write (nout,*) 'nmodes = ', nmodes
 108    format (1x,a8,' = ',1pe14.6,',')
c
      endif
c
      zgne = gnein
      zgnh = gnhin
      zgnz = gnzin
      zgth = gthin
      zgte = gtein
      zgtz = gtzin
c
c..check validity of input data
c
      if ( neq .lt. 2 ) call abortb (nout
     & ,' neq .lt. 2 in sbrtn etaw17a')
c
      if ( ndim .gt. idim ) call abortb (nout
     & ,' ndim .gt. idim in sbrtn etaw17a')
c
c
c..initialize arrays
c
      do j1=1,ndim
        omega(j1)  = 0.D0
        gamma(j1)  = 0.D0
        chieff(j1) = 0.D0
        velthi(j1) = 0.D0
        perform(j1) = 0.D0
        do j2=1,ndim
          difthi(j1,j2) = 0.D0
          zchim(j1,j2)  = 0.D0
          zflxm(j1,j2)  = 0.D0
        enddo
      enddo
c
c...JEK 9/15/98 Error fixed: reset zgm to zero
c
      do j1=1,idp
      do j2=1,idp
        zgm(j1,j2) = 0.D0
      enddo
      enddo
c
c..set up initial gradients
c
      zimp   = czin
      zmass  = azin
      zfnz   = fnzin * zimp
      zfns   = fnsin
c
c..Save zgns = - R ( d n_s / d r ) / n_e
c    where n_s is the fast ion density
c
      zgns = zgne - zgnh * ( zone - zfnz - zfns ) - zgnz * zfnz
c
cccccpis      if ( zfns .lt. zepsmach) zgns = 0
c
      if ( lprint .gt. 1 ) then
	  write (nout,109) 'gnsin', zgns
 109      format (1x,a8,' = ',1pe14.6
     &       ,', computed from quasi-neutrality')
      endif
c
c
c..set up gradient matrix zgm(j1,jd)
c
c  zdg = the finite difference used to construct the zgm matrix
c
      zdg = cetain(30)
      if ( abs(zdg) .lt. zepsqrt ) zdg = 1.D-2
c
c  input values for zgm matrix
c
      do jd=1,imatrx+1
        zgm(1,jd) = zgth
        zgm(2,jd) = zgnh
        zgm(3,jd) = zgte
        if ( ieq .gt. 4 ) then
          zgm(4,jd) = zgnz
          zgm(5,jd) = zgtz
        endif
      enddo
c
c  incremental values
c  Note:  if zgm(jd,jd) is negative, increment in negative direction
c
      do jd=1,imatrx+1
        zgm(jd,jd+1) = zgm(jd,jd+1) + sign ( zdg, zgm(jd,jd+1) )
      enddo
c
      if ( lprint .gt. 8 ) then
c
        write (nout,*)
        write (nout,*) '  zgm(j1,jd) matrix,  -> jd'
        do j1=1,ieq-1
          write (nout,132) (zgm(j1,jd),jd=1,ieq)
        enddo
c
      endif
c
c
c..loop over gradient matrix
c
      do 50 jd=1,imatrx+1
c
c..set variables
c
      do j1=1,ndim
        zomega(j1)  = 0.D0
        zgamma(j1)  = 0.D0
        zalfr (j1)  = 0.D0
        zalfi (j1)  = 0.D0
        zbeta (j1)  = 0.D0
        do j2=1,ndim
          zamr (j1,j2) = 0.D0
          zami (j1,j2) = 0.D0
          zbmr (j1,j2) = 0.D0
          zbmi (j1,j2) = 0.D0
          zvr  (j1,j2) = 0.D0
          zvi  (j1,j2) = 0.D0
        enddo
      enddo
c
c..compute gradient scale lengths from the zgm(j1,jd) matrix
c
      zgth = zgm(1,jd)
      zgnh = zgm(2,jd)
      zgte = zgm(3,jd)
      zgnz = zgm(4,jd)
      zgtz = zgm(5,jd)
c
c..Since the density gradients may have changed,
c  reconstruct the normalized electron density gradient.
c
      zgne = zgnh * ( zone - zfnz - zfns ) + zgnz * zfnz + zgns
c
c..compute the rest of the dimensionless variables needed
c
      ztauh  = tauhin
c
c  check for incompatibilities:  Weiland's definitions
c
      zft    = ftrapein
      zflh   = ekyrhoin**2.D0
c
      zeni = zhalf * zgne
c
      zbetae = max ( cetain(20) * betaein, zepsqrt )
      zimp   = max ( czin, zone )
      zmass  = max ( azin, zone )
      zfnz   = fnzin * zimp
c
c  ******  NOTE THE INVERSE DEFINITION OF ZTAUZ ! ******
c
c  end of Weiland's definitions
c
      zft    = ftrapein
      zflh   = ekyrhoin**2.D0
c
      if ( neq .gt. 4 ) then
        ztauz  = tauzin / czin
        zflz   = zmass * zflh / zimp**2.D0
      else
        ztauz  = zone
        zflz   = 0.D0
      endif
c
c..Note:
c
c  ztauz = T_Z / ( Z T_e )
c  zfnz  = Z n_Z / n_e
c  zimp  = Z
c  zmass = m_Z / m_H
c  zflh  = k_y^2 \rho_{sH}^2
c  zflz  = k_y^2 \rho_{sZ}^2
c
      zetae = zgte / zgne
      zetah = zgth / zgne
c
c..diagnostic output
c
      if ( lprint .gt. 2 ) then
        write (nout,*)
        write (nout,*) '--------------------------------------'
        write (nout,*)
        write (nout,*) ' jd = ',jd
       if ( lprint .gt. 2 .or. ( lprint .gt. 0 .and. jd .eq. 1 ) ) then
        write (nout,*)
c
        write (nout,*) zgnh,' = zgnh'
        write (nout,*) zgne,' = zgne'
        write (nout,*) zgnz,' = zgnz'
        write (nout,*) zgns,' = zgns'
        write (nout,*) zgth,' = zgth'
        write (nout,*) zgte,' = zgte'
        write (nout,*) zgtz,' = zgtz'
        write (nout,*) zetae,' = zetae'
        write (nout,*) zetah,' = zetah'
        write (nout,*)
c
        write (nout,*) ztauh,' = ztauh'
        write (nout,*) ztauz,' = ztauz'
        write (nout,*)
        write (nout,*) zft,' = zft'
        write (nout,*) zimp,' = zimp'
        write (nout,*) zmass,' = zmass'
        write (nout,*) fnzin,' = fnz'
        write (nout,*) zfnz,' = zfnz'
        write (nout,*) zfns,' = zfns'
        write (nout,*) zflh,' = zflh'
        write (nout,*) zflz,' = zflz'
        write (nout,*)
        write (nout,*) zepsqrt,' = zepsqrt'
        write (nout,*) zepsmach,' = zepsmach'
       endif
      endif
 
      if ( ieq .gt. 4 ) then
c
c..constants from Nilsson and Weiland, NF 34 (1994) 803  section 2
c  and from Weiland and Hirose NF 32 (1992) 151
c
c  ** Note:  zanorm = 1./kpc is used to normalize A_\parallel and K
c
        zanorm = cetain(11)
c
        bt  = 1.5D0
        bt1 = bt - ztwohlf
c
c  Expressions from Weiland and Hirose, NF 32 (1992) 151  Eq. (6)
c  Note that k1 rather than k1**2 appears under sqrt in zalp.
c  This is the result of a numerical fit rather than analytic solution.
c  Also, only free electrons are used from zbetae in zalf
c  Hence, zbetae -> zbetae * (zone-zft).
c
        if ( zgnh + zgth .lt. zepsqrt ) then
c
          zalp = 0.D0
          zalf = 0.D0
          zkpsh = 0.D0
          zkpsz = 0.D0
          zrav = 1.D0
c
        else
c
          k1 = zquarter * q * q * zflh *
     &         sqrt( zhalf * (zgnh + zgth) * ztauh / (zone-zft) )
          k2 = q * q * zflh * zflh * zhalf * (zgnh + zgth) *
     &         ztauh / ( zone-zft)
          zalp = zhalf * ( k1 + sqrt( k1 + shear * shear * k2) )
          zalf = zalp / (ztwo * zflh * q * q * zbetae * (zone-zft))
          zkpsh = cetain(10) * zhalf * sqrt( zalp / zflh ) / q
          zrav = zone
     &      + cetain(25) * ( zquarter * ( ztwo * shear - zone
     &      + kappa * kappa * ( shear - zone ) * ( shear - zone ) )
     &        / zalp )
          if ( zmass .gt. zepsqrt ) then
            zkpsz = zkpsh / sqrt ( zmass )
          else
            zkpsz = zkpsh
          endif
c
        endif
c
        if (lprint .gt. 5 ) then
        write (nout,*)
        write (nout,*) k1,'  = k1'
        write (nout,*) k2,'  = k2'
        write (nout,*) zalp,' = zalp'
        write (nout,*) zalf,' = zalf'
        write (nout,*) zkpsh,' = zkpsh'
        write (nout,*) zkpsz,' = zkpsz'
        write (nout,*) kappa,' = kappa'
        write (nout,*) zrav,' = rav'
        write (nout,*)
        write (nout,*) bt1,' = bt1'
        write (nout,*) zanorm,' = zanorm'
        endif
c
      endif
c
c
      if (ieq .eq. 11) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..Eleven equations with impurities, trapped electrons, FLR,
c  collisional effects, parallel hydrogenic and impurity ion motion,
c  electromagnetic (finite beta) effects
c
c  equations for
c  e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Vph, Av, K, Vpz
c
      if (lprint .gt. 5 ) then
      write(nout,*)
     &  'Ten eqns for '
     &  ,'e phi / T_e, T_H, n_H, T_et, n_Z, T_Z, F, Vph, Av, K, Vpz'
      endif
 
c
c  hydrogenic density
c
      zamr(1,1) = - zone + zhalf * (zgnh - zflh * ztauh * (zgnh + zgth))
      zamr(1,2) = - ztauh
      zamr(1,3) = - ztauh
      zamr(1,8) = zkpsh
c
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
c
c  hydrogenic energy
c
      zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
      zamr(2,2) = - ztauh * zftr
c
      zbmr(2,2) = zone
      zbmr(2,3) = - ztvr
c
c  trapped electron density expressed through quasineutrality
c
      zamr(3,1) = zft * zeni - zone
      zami(3,1) = vef * (zone-zft)
      zamr(3,3) = zone - zfnz - zfns
      zami(3,3) = - vef * (zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = - vef * zfnz
      zami(3,7) = vef * zft
      zamr(3,9) =  - ( zone - zft ) * zeni * zanorm
      zami(3,9) = vef * ( zone - zft ) * zeni * zanorm
      zamr(3,10) = ( zone - zft ) * (zone + zeni) * zanorm
      zami(3,10) = - vef * ( zone - zft ) * zanorm
c
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
      zbmr(3,10) = ( zone - zft ) * zanorm
c
c  trapped electron energy
c
      zamr(4,1) = zft * zhalf * ( zgte - ztvr * zgne )
      zami(4,1) = vef * ztvr * (bt - ztwohlf * (zone-zft))
      zami(4,3) = - vef * ztvr * bt1 * (zone - zfnz - zfns)
      zamr(4,4) = zft * zftr
      zami(4,5) = - vef * ztvr * bt1 * zfnz
      zami(4,7) = - zftr * vef * zft
c
      zbmr(4,1) = (zone - zft) * ztvr
      zbmr(4,3) = - (zone - zfnz - zfns) * ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * ztvr
c
c  impurity density
c
      zamr(5,1) = - zone + zhalf*(zgnz - zflz * ztauz * (zgnz + zgtz))
      zamr(5,5) = - ztauz
      zamr(5,6) = - ztauz
      zamr(5,11) = zkpsz
c
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
c
c  impurity energy
c
      zamr(6,1) = zhalf * (zgtz - ztvr * zgnz)
      zamr(6,6) = - ztauz * zftr
c
      zbmr(6,5) = - ztvr
      zbmr(6,6) = zone
c
c  variable F
c
      zamr(7,1) = zhalf * zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = - vef
c
      zbmr(7,1) = - zone
      zbmr(7,7) = zone
c
c  hydrogenic parallel ion motion Vp = Vpi/Cs
c
      zamr(8,1) = zkpsh
      zamr(8,2) = zkpsh * ztauh
      zamr(8,3) = zkpsh * ztauh
      zamr(8,9) = - zkpsh * zhalf * ( zgth + zgnh ) * ztauh
c
      zbmr(8,8) = zone
      zbmr(8,9) = zkpsh
c
c  electromagnetic parallel vector potential Av = e A_par /Te
c
      zamr(9,1) = zeni
      zamr(9,8) = zkpsh * ( zone - zfnz - zfns ) / ( zone - zft )
      zamr(9,9) = (zhalf * (zgne + zgte) - zalf * zflh * zrav)*zanorm
      zamr(9,10) = - zeni * zanorm
      zamr(9,11) = zkpsz * zfnz / ( zone - zft )
c
      zbmr(9,1) = zone
      zbmr(9,9) = zanorm
      zbmr(9,10) = - zanorm
c
c  time derivative of Av, K = omega * Av
c
      zamr(10,10) = zone
c
      zbmr(10,9) = zone
c
c  impurity parallel ion motion
c
      zamr(11,1) = zimp * zkpsz
      zamr(11,5) = zimp * zkpsz * ztauz
      zamr(11,6) = zimp * zkpsz * ztauz
      zamr(11,9) = - zimp * zkpsz * zhalf * ( zgtz + zgnz ) * ztauz
c
      zbmr(11,9) = zimp * zkpsz
      zbmr(11,11) = zone
c
c
      elseif (ieq .eq. 10) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..Ten equations with impurities, trapped electrons, FLR,
c  collisional effects, parallel ion motion, and electromagnetic effects
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Vp, Av, K
c
      if (lprint .gt. 5 ) then
      write(nout,*)
     &'Ten eqns for e phi/T_e, Th, n_h, Tet, n_z, T_z, F, Av, K'
      endif
c
c  hydrogenic density
c
      zamr(1,1) = - zone + zhalf * (zgnh - zflh * ztauh * (zgnh + zgth))
      zamr(1,2) = - ztauh
      zamr(1,3) = - ztauh
      zamr(1,8) = zkpsh
c
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
c
c  hydrogenic energy
c
      zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
      zamr(2,2) = - ztauh * zftr
c
      zbmr(2,2) = zone
      zbmr(2,3) = - ztvr
c
c  trapped electron density expressed through quasineutrality
c
      zamr(3,1) = zft * zeni - zone
      zami(3,1) = vef * (zone-zft)
      zamr(3,3) = zone - zfnz - zfns
      zami(3,3) = - vef * (zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = - vef * zfnz
      zami(3,7) = vef * zft
      zamr(3,9) =  - ( zone - zft ) * zeni * zanorm
      zami(3,9) = vef * ( zone - zft ) * zeni * zanorm
      zamr(3,10) = ( zone - zft ) * (zone + zeni) * zanorm
      zami(3,10) = - vef * ( zone - zft ) * zanorm
c
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
      zbmr(3,10) = ( zone - zft ) * zanorm
c
c  trapped electron energy
c
      zamr(4,1) = zft * zhalf * ( zgte - ztvr * zgne )
      zami(4,1) = vef * ztvr * (bt - ztwohlf * (zone-zft))
      zami(4,3) = - vef * ztvr * bt1 * (zone - zfnz - zfns)
      zamr(4,4) = zft * zftr
      zami(4,5) = - vef * ztvr * bt1 * zfnz
      zami(4,7) = - zftr * vef * zft
c
      zbmr(4,1) = (zone - zft) * ztvr
      zbmr(4,3) = - (zone - zfnz - zfns) * ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * ztvr
c
c  impurity density
c
      zamr(5,1) = - zone + zhalf*(zgnz - zflz * ztauz * (zgnz + zgtz))
      zamr(5,5) = - ztauz
      zamr(5,6) = - ztauz
c
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
c
c  impurity energy
c
      zamr(6,1) = zhalf * (zgtz - ztvr * zgnz)
      zamr(6,6) = - ztauz * zftr
c
      zbmr(6,5) = - ztvr
      zbmr(6,6) = zone
c
c  variable F
c
      zamr(7,1) = zhalf * zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = - vef
c
      zbmr(7,1) = - zone
      zbmr(7,7) = zone
c
c  parallel ion motion Vp = Vpi/Cs
c
      zamr(8,1) = zkpsh
      zamr(8,2) = zkpsh * ztauh
      zamr(8,3) = zkpsh * ztauh
      zamr(8,9) = - zkpsh * zhalf * ( zgth + zgnh ) * ztauh
c
      zbmr(8,8) = zone
      zbmr(8,9) = zkpsh
c
c  electromagnetic parallel vector potential Av = e A_par /Te
c
      zamr(9,1) = zeni
      zamr(9,8) = zkpsh * ( zone - zfnz - zfns ) / ( zone - zft )
      zamr(9,9) = (zhalf * (zgne + zgte) - zalf * zflh * zrav)*zanorm
      zamr(9,10) = - zeni * zanorm
c
      zbmr(9,1) = zone
      zbmr(9,9) = zanorm
      zbmr(9,10) = - zanorm
c
c  time derivative of Av, K = omega * Av
c
      zamr(10,10) = zone
c
      zbmr(10,9) = zone
c
c
      elseif (ieq .eq. 9) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..Nine equations with impurities, trapped electrons, FLR,
c  collisional effects, and electromagnetic effects
c  and parallel ion motion in strong ballooning limit
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Av, K
c
      if (lprint .gt. 5 ) then
      write(nout, *)
     &'Nine eqns for e phi/T_e, Th, n_h, Tet, n_z, T_z, F, Av, K'
      endif
c
      H = cetain(12) * zhalf * abs( shear ) / max ( q, zepsqrt )
c
c  ion continuity
c
      zamr(1,1) = - zone + zhalf * (zgnh - zflh * ztauh * (zgnh + zgth))
      zami(1,1) = - H
      zamr(1,2) = - ztauh
      zami(1,2) = - ztauh*H
      zamr(1,3) = - ztauh
      zami(1,3) = - ztauh*H
c
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
c
c  ion energy
c
      zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
      zamr(2,2) = - ztauh * zftr
c
      zbmr(2,2) = zone
      zbmr(2,3) = - ztvr
c
c  total electron density expressed through quasineutrality
c
      zamr(3,1) = zft * zeni - zone
      zami(3,1) = vef * (zone-zft)
      zamr(3,3) = zone - zfnz - zfns
      zami(3,3) = - vef * (zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = - vef * zfnz
      zami(3,7) = vef * zft
      zamr(3,8) = - ( zone - zft ) * zeni * zanorm
      zami(3,8) = vef * ( zone - zft ) * zeni * zanorm
      zamr(3,9) = ( zone - zft ) * (zone + zeni) * zanorm
      zami(3,9) = - vef * ( zone - zft ) * zanorm
c
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
      zbmr(3,9) = ( zone - zft ) * zanorm
c
c  trapped electron energy
c
      zamr(4,1) = zft * zhalf * ( zgte - ztvr * zgne )
      zami(4,1) = vef * ztvr * (bt - ztwohlf * (zone-zft))
      zami(4,3) = - vef * ztvr * bt1 * (zone - zfnz - zfns)
      zamr(4,4) = zft * zftr
      zami(4,5) = - vef * ztvr * bt1 * zfnz
      zami(4,7) = - zftr * vef * zft
c
      zbmr(4,1) = (zone - zft) * ztvr
      zbmr(4,3) = - (zone - zfnz - zfns) * ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * ztvr
c
c  impurity density
c
      zamr(5,1) = -zone+zhalf*(zgnz-zflz*ztauz*(zgnz+zgtz))
      zamr(5,5) = -ztauz
      zamr(5,6) = -ztauz
c
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
c
c  impurity energy
c
      zamr(6,1) = zhalf*(zgtz-ztvr*zgnz)
      zamr(6,6) = -ztauz*zftr
c
      zbmr(6,5) = -ztvr
      zbmr(6,6) = zone
c
c  variable F
c
      zamr(7,1) = zhalf*zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = -vef
c
      zbmr(7,1) = -zone
      zbmr(7,7) = zone
c
c  electromagnetic parallel vector potential Av = e A_par /Te
c
      zamr(8,1) = zeni
      zamr(8,8) = ( zhalf*(zgne+zgte)-zalf*zflh*zrav ) * zanorm
      zamr(8,9) = - zeni * zanorm
c
      zbmr(8,1) = zone
      zbmr(8,8) = zanorm
      zbmr(8,9) = - zanorm
c
c   time derivative of Av
c
      zamr(9,9) = zone
c
      zbmr(9,8) = zone
c
c
      elseif (ieq .eq. 8) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..Eight equations with impurities, trapped electrons, FLR,
c  collisional effects, and parallel ion motion
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Vp
c
      if (lprint .gt. 5 ) then
      write(nout,*)
     &'Eight eqns for e phi/T_e, Th, n_h, Tet, n_z, T_z, F, Vp'
      endif
c
c  ion continuity
c
      zamr(1,1) = -zone + zhalf*(zgnh-zflh*ztauh*(zgnh+zgth))
      zamr(1,2) = -ztauh
      zamr(1,3) = -ztauh
      zamr(1,8) = zkpsh
c
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
c
c  ion energy
c
      zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
      zamr(2,2) = -ztauh * zftr
c
      zbmr(2,2) = zone
      zbmr(2,3) = -ztvr
c
c  total electron density expressed through quasineutrality
c
      zamr(3,1) = zft*zeni - zone
      zami(3,1) = vef*(zone-zft)
      zamr(3,3) = zone - zfnz- zfns
      zami(3,3) = -vef*(zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = -vef*zfnz
      zami(3,7) = vef*zft
c
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
c
c  trapped electron energy
c
      zamr(4,1) = zft*zhalf*(zgte-ztvr*zgne)
      zami(4,1) = vef*ztvr*(bt-ztwohlf*(zone-zft))
      zami(4,3) = -vef*ztvr*bt1*(zone - zfnz - zfns)
      zamr(4,4) = zft*zftr
      zami(4,5) = -vef*ztvr*bt1*zfnz
      zami(4,7) = -zftr*vef*zft
c
      zbmr(4,1) = (zone-zft)*ztvr
      zbmr(4,3) = - (zone - zfnz - zfns)*ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = -zfnz*ztvr
c
c  impurity density
c
      zamr(5,1) = -zone+zhalf*(zgnz-zflz*ztauz*(zgnz+zgtz))
      zamr(5,5) = -ztauz
      zamr(5,6) = -ztauz
c
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
c
c  impurity energy
c
      zamr(6,1) = zhalf*(zgtz-ztvr*zgnz)
      zamr(6,6) = -ztauz*zftr
c
      zbmr(6,5) = -ztvr
      zbmr(6,6) = zone
c
c  variable F
c
      zamr(7,1) = zhalf*zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = -vef
c
      zbmr(7,1) = -zone
      zbmr(7,7) = zone
c
c  parallel ion motion Vp = Vpi/Cs
c
      zamr(8,1) = zkpsh
      zamr(8,2) = zkpsh*ztauh
      zamr(8,3) = zkpsh*ztauh
c
      zbmr(8,8) = zone
c
c
      elseif ( ieq .eq. 7 ) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..Seven equations with impurities, trapped electrons, FLR,
c  and collisions
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z and F
c  Here, F is defined as F = GM*e phi/T_e
c  where GM=1+etae/(epsn*(omega-1+i*vef))
c
      if ( lprint .gt. 5 ) then
      write (nout,*)
      write (nout,*)
     & 'Seven eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z and F'
      endif
c
c  hydrogen density
c
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
c
c  hydrogen energy
c
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(2,2) = zone
        zbmr(2,3) = - ztvr
c
c  trapped electron density
c
        zamr(3,1) = - zone + zft * zeni
        zami(3,1) = vef*(zone-zft)
        zamr(3,3) = zone - zfnz - zfns
        zami(3,3) = -vef*(zone - zfnz - zfns)
        zamr(3,4) = zft
        zamr(3,5) = zfnz
        zami(3,5) = -vef*zfnz
        zami(3,7) = vef*zft
c
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft*zhalf*(zgte-ztvr*zgne)
        zami(4,1) = vef*ztvr*(bt-ztwohlf*(zone-zft))
        zami(4,3) = -vef*ztvr*bt1*(zone - zfnz - zfns)
        zamr(4,4) = zft * zftr
        zami(4,5) = -vef*ztvr*bt1*zfnz
        zami(4,7) = -zftr*vef*zft
c
        zbmr(4,1) = ( zone - zft ) *ztvr
        zbmr(4,3) = - ( zone - zfnz - zfns ) *ztvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * ztvr
c
c  impurity density
c
        zamr(5,1) = - zone+zhalf*(zgnz-zflz*ztauz*(zgnz+zgtz))
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
c
        zbmr(5,1) = zflz
        zbmr(5,5) = zone
c
c  impurity energy
c
        zamr(6,1) = zhalf*(zgtz-ztvr*zgnz)
        zamr(6,6) = - ztauz * zftr
c
        zbmr(6,5) = - ztvr
        zbmr(6,6) = zone
c
c  variable F
c
        zamr(7,1) = zhalf*zgte - zone
        zami(7,1) = vef
        zamr(7,7) = zone
        zami(7,7) = -vef
c
        zbmr(7,1) = -zone
        zbmr(7,7) = zone
c
c
      elseif ( ieq .eq. 6 ) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..Six equations with impurities, trapped electrons, and FLR
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, and T_Z
c
      if ( lprint .gt. 1 .and. jd .eq. 1 ) then
        write (nout,*)
        write (nout,*)
     &   'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
      endif
c
c  hydrogen density
c
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
c
c  hydrogen energy
c
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(2,2) = zone
        zbmr(2,3) = - ztvr
c
c  trapped electron density
c
        zamr(3,1) = zft*zeni - zone
        zamr(3,3) = zone - zfnz - zfns
        zamr(3,4) = zft
        zamr(3,5) = zfnz
c
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft * zhalf * (zgte - ztvr*zgne)
        zamr(4,4) = zft * zftr
c
        zbmr(4,1) = ( zone - zft ) * ztvr
        zbmr(4,3) = - ( zone - zfnz - zfns ) * ztvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * ztvr
c
c  impurity density
c
        zamr(5,1) = -zone + zhalf * (zgnz - zflz*ztauz*(zgnz+zgtz))
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
c
        zbmr(5,1) = zflz
        zbmr(5,5) = zone
c
c  impurity energy
c
        zamr(6,1) = zhalf * (zgtz - ztvr*zgnz)
        zamr(6,6) = - ztauz * zftr
c
        zbmr(6,5) = - ztvr
        zbmr(6,6) = zone
c
c
      elseif ( ieq .eq. 5 ) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..5 equations with trapped electrons, FLR effects, and parallel ion motion
c
c  equations for e phi/T_e, T_H, n_i, T_e, and Vp
c
       if ( lprint .gt. 1 ) then
         write (nout,*)
         write (nout,*)
     &    ' Five eqns for e phi/T_e, T_H, n_H, T_e, and Vp'
       endif
c
c  ion continuity
c
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
        zamr(1,5) = zkpsh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
c
c  ion energy
c
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(2,2) =   zone
        zbmr(2,3) = - ztvr
c
c  trapped electron continuity
c
c   Calculates the total electron density perturbation and replaces it
c   by the ion density perturbation.
c   The dilution factor 1-zfnz has now been added.
c
        zamr(3,1) = - zone + zft * zhalf * zgne
        zamr(3,3) = zone - zfnz - zfns
        zamr(3,4) = zft
c
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
c
c  trapped electron energy
c
        zamr(4,1) = zft * ( zgte - ztvr * zgne ) * zhalf
        zamr(4,4) = zft * zftr
c
        zbmr(4,1) = ( zone - zft ) * ztvr
        zbmr(4,3) = - zone * ztvr
        zbmr(4,4) = zft
c
c
c  Parallel ion motion Vpi/Cs
c
        zamr(5,1) = zkpsh
        zamr(5,2) = zkpsh*ztauh
        zamr(5,3) = zkpsh*ztauh
c
        zbmr(5,5) = zone
c
      elseif ( ieq .eq. 4 ) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..4 equations with trapped electrons and FLR effects
c
c  equations for e phi/T_e, T_H, n_i, and T_e
c
       if ( lprint .gt. 1 ) then
         write (nout,*)
         write (nout,*) ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
       endif
c
c  ion continuity
c
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
c
c  ion energy
c
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(2,2) =   zone
        zbmr(2,3) = - ztvr
c
c  trapped electron continuity
c
c   Calculates the total electron density perturbation and replaces it
c   The dilution factor 1-zfnz has now been added.
c
        zamr(3,1) = - zone + zft * zhalf * zgne
        zamr(3,3) = zone - zfnz - zfns
        zamr(3,4) = zft
c
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
c
c  trapped electron energy
c
        zamr(4,1) = zft * ( zgte - ztvr * zgne ) * zhalf
        zamr(4,4) = zft * zftr
c
        zbmr(4,1) = ( zone - zft ) * ztvr
        zbmr(4,3) = - zone * ztvr
        zbmr(4,4) = zft
c
      elseif ( ieq .eq. 2 ) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..two equations when trapped particles and FLR effects omitted
c
c  equations for e phi/T_e and T_H
c
       if ( lprint .gt. 1 .and. jd .eq. 1 ) then
         write (nout,*)
         write (nout,*) ' Two eqns for e phi/T_e and T_H'
       endif
c
        zamr(1,1) = zhalf * zgnh - ztauh - zone
        zamr(1,2) = - ztauh
        zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(1,1) = zone
        zbmr(1,2) = 0.D0
        zbmr(2,1) = - ztvr
        zbmr(2,2) = zone
c
c
      else
c
        write (nout,*)
        write (nout,*) ieq,' = ieq in sbrtn etaw17a'
        call abortb(nout
     & ,'the value of ieq is wrong in sbrtn etaw17a')
c
      endif
c
c..find the eigenvalues and eigenvectors using ACM/TOMS routine 535
c
      ifail = -1
c
c  use ACM/TOMS routine 535 for complex matricies
c
       if ( lprint .gt. 1 .and. jd .eq. 1 ) then
         write (nout,*)
         write (nout,*)
     &    ' Eigenvalues computed using complex ACM/TOMS routine 535'
       endif
 
c..fill complex matrices zamc and zbmc
c
      do j1=1,ieq
       do j2=1,ieq
        zamr(j1,j2)=zamr(j1,j2)+abs(wexb)*zbmi(j1,j2)
        zami(j1,j2)=zami(j1,j2)-abs(wexb)*zbmr(j1,j2)
        zamrt(j1,j2)= zamr(j1,j2)
        zamit(j1,j2)= zami(j1,j2)
        zbmrt(j1,j2)= zbmr(j1,j2)
        zbmit(j1,j2)= zbmi(j1,j2)
       enddo
      enddo
 
c
c..diagnostic output
c
      if ( lprint .gt. 6 ) then
c
        write (nout,*)
        write (nout,*) ieq
c
        write (nout,*)
        write (nout,*) ' zamr(j1,j2)  j2 ->'
        do j1=1,ieq
          write (nout,192) (zamr(j1,j2),j2=1,ieq)
        enddo
c
        write (nout,*)
        write (nout,*) ' zami(j1,j2)  j2 ->'
        do j1=1,ieq
          write (nout,192) (zami(j1,j2),j2=1,ieq)
        enddo
c
        write (nout,*)
        write (nout,*) ' zbmr(j1,j2)  j2->'
        do j1=1,ieq
          write (nout,192) (zbmr(j1,j2),j2=1,ieq)
        enddo
c
        write (nout,*)
        write (nout,*) ' zbmi(j1,j2)  j2->'
        do j1=1,ieq
          write (nout,192) (zbmi(j1,j2),j2=1,ieq)
        enddo
 192  format (1p10e12.4)
      endif
 
      call r8tomsqz(idim,ieq,zamr,zami,zbmr,zbmi,ZALFR,ZALFI,ZBETA,
     &            ZVR,ZVI,IFAIL)
 
      do j=1,ieq
        ztemp1 = zbeta(j)
        if ( abs(zbeta(j)) .lt. zepsqrt ) ztemp1 = zepsqrt
        zomega(j) = zalfr(j)  / ztemp1
        zgamma(j) = zalfi(j)  / ztemp1
        do j1=1,ieq
          zevec(j1,j) = cmplx ( zvr(j1,j), zvi(j1,j))
        end do
      enddo
 
      nerr = ifail
 
c
c..save growth rates and frequencies during first element of jd loop
c
      if ( jd .eq. 1 ) then
        zgamax = 0.D0
        do j=1,ieq
          omega(j) = zomega(j)
          gamma(j) = zgamma(j)
          zgamax = max ( zgamax, zgamma(j) )
        enddo
        if ( zgamax .lt. zepsqrt ) go to 90
      endif
c
c..check the eigenfunctions
c
      if ( lprint .gt. 12 ) then
        write (nout,*)
        write (nout,*) ' Checking eigenfunctions'
      endif
c
c  Real and imaginary parts
c
      zerrmax = 0.D0
c
      do j=1,ieq
c
        do j1=1,ieq
c
          ztempa(j1) = 0.D0
          ztempb(j1) = 0.D0
c
          do j2=1,ieq
              zerreal =
     &            zamrt(j1,j2) * real(zevec(j2,j))
     &          - zamit(j1,j2) * aimag(zevec(j2,j))
     &          - zbmrt(j1,j2) * real(zevec(j2,j))  * zomega(j)
     &          + zbmrt(j1,j2) * aimag(zevec(j2,j)) * zgamma(j)
     &          + zbmit(j1,j2) * real(zevec(j2,j))  * zgamma(j)
     &          + zbmit(j1,j2) * aimag(zevec(j2,j)) * zomega(j)
              zerimag =
     &            zamrt(j1,j2) * aimag(zevec(j2,j))
     &          + zamit(j1,j2) * real(zevec(j2,j))
     &          - zbmrt(j1,j2) * aimag(zevec(j2,j)) * zomega(j)
     &          - zbmrt(j1,j2) * real(zevec(j2,j))  * zgamma(j)
     &          + zbmit(j1,j2) * aimag(zevec(j2,j)) * zgamma(j)
     &          - zbmit(j1,j2) * real(zevec(j2,j))  * zomega(j)
c
              ztempa(j1) = ztempa(j1) + zerreal
              ztempb(j1) = ztempb(j1) + zerimag
 
            enddo
c
            zerrmax = max ( zerrmax, abs(ztempa(j1)), abs(ztempb(j1)) )
c
          enddo
c
          if ( lprint .gt. 12 ) then
            write (nout,*)
            write (nout,*) ' LHS - RHS for j =  ',j
            do j1=1,ieq
              write (nout,142) ztempa(j1), ztempb(j1)
            enddo
          endif
c
        enddo
c
        if ( lprint .gt. 0 ) then
          if (abs(zerrmax) .gt. zepsqrt) nerr = max (nerr , 1)
          write (nout,*)
          write (nout,*) zerrmax,' = zerrmax', nmodes
          write (nout,*) nerr,' = nerr, error in eigenvalue in etaw17a'
        endif
c
 142  format (1p10e12.4)
c
c
c..compute effective diffusivities directly from eigenvalues
c  assume eigenvectors are arranged in the order of
c  e\phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z
c
        nmodes = 0
c
        do j=1,ieq
c
          ztemp1 =  real(zevec(1,j)) *  real(zevec(1,j))
     &           + aimag(zevec(1,j)) * aimag(zevec(1,j))
c
          if ( zgamma(j) .gt. zepsqrt
     &      .and. ztemp1 .gt. zepsqrt
     &       ) then
c
            nmodes = nmodes + 1
c
c
            zreal =  real(zevec(2,j)) +  real(zevec(3,j))
            zimag = aimag(zevec(2,j)) + aimag(zevec(3,j))
c
            zphsph = - ( zimag * real(zevec(1,j))
     &          - zreal * aimag(zevec(1,j)) ) / ztemp1
c
            zflxph = - 2.D0 * ( abs ( zgamma(j) ) )**2.D0
     &        * ( zimag * real(zevec(1,j))
     &          - zreal * aimag(zevec(1,j)) ) / ztemp1
c
            zflxm(1,jd) = zflxm(1,jd) + zflxph
c
c            if ( lprint .gt. 7 ) then
c              write (nout,*)
c              write (nout,*) zflxph,' = zflxph for j = ',j
c            endif
c
            zphsnh = - ( aimag(zevec(3,j)) * real(zevec(1,j))
     &          - real(zevec(3,j)) * aimag(zevec(1,j)) ) / ztemp1
c
 
            zflxnh = - 2.D0 * ( abs ( zgamma(j) ) )**2.D0
     &        * ( aimag(zevec(3,j)) * real(zevec(1,j))
     &          - real(zevec(3,j)) * aimag(zevec(1,j)) ) / ztemp1
c
            zflxm(2,jd) = zflxm(2,jd) + zflxnh
c
c            if ( lprint .gt. 7 ) then
c              write (nout,*)
c              write (nout,*) zflxnh,' = zflxnh for j = ',j
c            endif
c
            if ( ieq .gt. 3 ) then
c
              zreal =  real(zevec(4,j))
     &          + ( zone - zfnz - zfns ) * real(zevec(2,j))
     &          + zfnz * real(zevec(5,j))
              zimag = aimag(zevec(4,j))
     &          + ( zone - zfnz - zfns ) * aimag(zevec(2,j))
     &          + zfnz * aimag(zevec(5,j))
c
c  Note, the electron heat flux is reduced by the fraction of
c    trapped electrons
c
              zphspe = - zft * ( zimag * real(zevec(1,j))
     &          - zreal * aimag(zevec(1,j)) ) / ztemp1
c
              zflxpe =
     &          - 2.D0 * zft * (abs( zgamma(j) ))**2.D0
     &        * ( zimag * real(zevec(1,j))
     &          - zreal * aimag(zevec(1,j)) ) / ztemp1
c
              zflxm(3,jd) = zflxm(3,jd) + zflxpe
c
c              if ( lprint .gt. 7 ) then
c                write (nout,*)
c                write (nout,*) zflxpe,' = zflxpe for j = ',j
c              endif
            endif
c
            if ( ieq .gt. 4 ) then
c
              zphsnz = - ( aimag(zevec(5,j)) * real(zevec(1,j))
     &          - real(zevec(5,j)) * aimag(zevec(1,j)) ) / ztemp1
c
              zflxnz = - 2.D0 * (abs( zgamma(j) ))**2.D0
     &        * ( aimag(zevec(5,j)) * real(zevec(1,j))
     &          - real(zevec(5,j)) * aimag(zevec(1,j)) ) / ztemp1
c
              zflxm(4,jd) = zflxm(4,jd) + zflxnz
c
c              if ( lprint .gt. 7 ) then
c                write (nout,*)
c                write (nout,*) zflxnz,' = zflxnz for j = ',j
c              endif
            endif
c
            if ( ieq .gt. 5 ) then
c
              zreal =  real(zevec(5,j)) +  real(zevec(6,j))
              zimag = aimag(zevec(5,j)) + aimag(zevec(6,j))
c
              zphspz = - ( zimag * real(zevec(1,j))
     &          - zreal * aimag(zevec(1,j)) ) / ztemp1
c
              zflxpz = - 2.D0 * (abs( zgamma(j) ))**2.D0
     &        * ( zimag * real(zevec(1,j))
     &          - zreal * aimag(zevec(1,j)) ) / ztemp1
c
              zflxm(5,jd) = zflxm(5,jd) + zflxpz
c
c              if ( lprint .gt. 7 ) then
c                write (nout,*)
c                write (nout,*) zflxpz,' = zflxpz for j = ',j
c              endif
            endif
c			
c..header for diagnostic printout of frequencies, fluxes, and phases
c
        if ( letain(29) .gt. 0 .and.  jd .eq. 1 ) then
          write (nout,134)
c
 134  format (//' Diagnostic printout of frequencies, fluxes and phases'
     &  /' Note: fluxph = flux of hydrogenic thermal energy'
     &  ,' = 2.0 * gamma^2 * phaseph'
     &  /9x,' (phaseph is related to the phases of the perturbations)'
     &  /7x,'fluxnh = flux of hydrogenic ions = 2.0 * gamma^2 * phasenh'
     &  /7x,'fluxpe = flux of electron thermal energy'
     &  ,' = 2.0 * gamma^2 * phasepe'
     &  /7x,'fluxnz = flux of impurity ions = 2.0 * gamma^2 * phasenz'
     &  /7x,'fluxpz = flux of impurity thermal energy'
     &  ,' = 2.0 * gamma^2 * phasepz'
     &  //1x,'radius',t10,'omega',t20,'gamma'
     &  ,t30,'fluxph',t40,'phaseph',t50,'fluxnh',t60,'phasenh'
     &  ,t70,'fluxpe',t80,'phasepe',t90,'fluxnz',t100,'phasenz'
     &  ,t110,'fluxpz',t120,'phasepz  #m')
c
c..diagnostic printout of frequencies and fluxes mode by mode
c
        write (nout,135) cetain(29), zomega(j), zgamma(j)
     &    , zflxph, zphsph, zflxnh, zphsnh, zflxpe, zphspe
     &    , zflxnz, zphsnz, zflxpz, zphspz
        endif
c
 135  format (0pf7.3,1p12e10.2,' #m')
c
          endif
c
        enddo
c
c..compute effective total diffusivities
c
        zchim(1,jd) = zflxm(1,jd)
     &   / sign ( max ( abs ( zgth ), zepsqrt ),  zgth )
        zchim(2,jd) = zflxm(2,jd)
     &   / sign ( max ( abs ( zgnh ), zepsqrt ),  zgnh )
        zchim(3,jd) = zflxm(3,jd)
     &   / sign ( max ( abs ( zgte ), zepsqrt ),  zgte )
        zchim(4,jd) = zflxm(4,jd)
     &   / sign ( max ( abs ( zgnz ), zepsqrt ),  zgnz )
        zchim(5,jd) = zflxm(5,jd)
     &   / sign ( max ( abs ( zgtz ), zepsqrt ),  zgtz )
c
      if ( lprint .gt. 2 ) then
c
c..print eigenvalues and eigenfunctions
c
        write (nout,121)
        do j=1,ieq
          write (nout,122) zomega(j), zgamma(j)
        enddo
 121    format (/' Solution of the eigenvalue equations'
     &   /t4,'zomega',t18,'zgamma')
 122    format (1p2e14.5,i5)
c
        write (nout,*)
        write (nout,*) ' Effective diffusivities'
     &    ,' normalized by omega_{De} / k_y^2'
c
        write (nout,130)
        write (nout,132) (zchim(j1,jd),j1=1,ieq-1)
c
      endif
c
      if ( lprint .gt. 99 ) then
c
        write (nout,*)
        write (nout,*) ' Eigenvectors zevec(j1,j2) j2->'
        do j1=1,ieq
          write (nout,124) (zevec(j1,j2),j2=1,ieq)
 124      format (2(1p12e11.3))
        enddo
c
        write (nout,*)
      endif
c
      if ( ifail .gt. 0 ) call abortb ( nout
     & ,'ifail .gt. 0 after call f02gjf in sbrtn etaw17a')
c
c
c..end of loop over diffusivity matrix
c
  50  continue
c
c
c..save effective diffusivities
c
      do j1=1,min(ieq-1,5)
        chieff(j1) = zflxm(j1,1) /
     &    sign( max( abs(zgm(j1,1)), zepsqrt), zgm(j1,1) )
      enddo
c
      if ( imatrx .lt. 1 ) go to 90
 
c
c..computation of difthi(j1,j2) from zflxm(j1,jd) and zgm(j1,jd)
c
      do j1=1,imatrx
c
        ztemp2 = 0.D0
c
        do jd=1,imatrx
c
          ztemp1 = zgm(jd,jd) - zgm(jd,1)
          if ( abs(ztemp1) .lt. zepsqrt )
     &       ztemp1 = zdg * sign(zone,zgm(jd,1))
c
          difthi(j1,jd) = ( zflxm(j1,jd+1) - zflxm(j1,1) ) / ztemp1
          ztemp2 = ztemp2 + difthi(j1,jd) * zgm(jd,1)
c
        enddo
c
c
        if ( letain(7) .eq. 0 ) then
c
c..compute convective velocities
c
          velthi(j1) = zflxm(j1,1) - ztemp2
c
        else
c
c..alternatively, set the convective velocities to zero
c    and rescale the diffusivity matrix
c
          velthi(j1) = 0.D0
          if ( abs(ztemp2) .lt. zepsqrt )
     &      ztemp2 = sign ( zepsqrt, ztemp2 )
          do jd=1,imatrx
            difthi(j1,jd) = difthi(j1,jd) * zflxm(j1,1) / ztemp2
          enddo
c
        endif
c
      enddo
c
c..diagnostic printout
c
      if ( lprint .gt. 6 ) then
c
        write (nout,*)
        write (nout,*) ' matrix zflxm(j1,jd)  -> jd'
        do j1=1,ieq-1
          write (nout,132) (zflxm(j1,jd),jd=1,ieq-1)
        enddo
c
      endif
c
  90  continue
c
      if ( lprint .gt. 0 ) then  ! changed for output
c
        write (nout,*)
        write (nout,*) '--------------------------------------'
     &    ,'  Final results from sbrtn etaw17a'
c
        write (nout,*)
        write (nout,*) ' diffusion matrix from eigenvectors'
     &    ,' normalized by omega_{De} / k_y^2'
        write (nout,*)
        write (nout,*) ' difthi(j1,j2)  -> j2'
        do j1=1,imatrx+1
          write (nout,132) (difthi(j1,j2),j2=1,imatrx+1)
        enddo
c
        write (nout,*)
        write (nout,*) ' convective velocities'
     &    ,' normalized by omega_{De} / k_y'
        do j1=1,imatrx
          write (nout,132) velthi(j1)
        enddo
c
        write (nout,*)
        write (nout,*) ' Effective diffusivities from eigenvectors'
     &    ,' normalized by omega_{De} / k_y^2'
c
        write (nout,130)
        write (nout,132) (chieff(j1),j1=1,imatrx)
c
      endif
c
c
      return
c
 130    format (/t3,'chi_i',t16,'dif_h',t29,'chi_e'
     &    ,t42,'dif_Z',t55,'chi_Z')
 132    format (1p11e12.4)
c
      end
c
      subroutine abortb (nout,text)
c
c   Subrtn abortb prints out text and stops run
c
      character *(*) text
c
      write (nout,*) text
c
      stop
      end
