cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine vtrans_dv(dt,ierror)
c
      implicit none
c
      include '../inc/ptor.m'
      include '../inc/vtrans.m'
c
c INPUT variables:
c
      real*8 dt                 ! time step
                                ! or dV/dr at zone outer-boundary (j=2)
c WORK arrays:
      real*8 deltaT(nfields,ngrid)
      real*8 ab(6*nfields-2,nfields,ngrid)
c     integer nwork(ngrid*nfields)    ! Work array for lapack in vtrans.
c
c LOCAL variables:
      integer nhi               ! HIGHEST GRID POINT TO SOLVE FOR (=ngrid-1)
      integer ldab              ! leading dimension of ab
      real*8 rhop,rhom,up,um,ap,am,ca
      integer i,j,k,klower,kupper,ioff,ioffp,ioffm,nsize,ierror,isp
c
c      ca = 2.D0/3.D0
      ca = 1.D0
      ap = (1.D0 + ca)/2.D0
      am = (1.D0 - ca)/4.D0
      do k=1,mxgrid
        do i=1,mxfields
          Tnew(i,k) = 0.D0
        enddo
      enddo
c
c  Set up coefficients of an implicit equation of the form Ax=b,
c  where A is stored in LAPACK/NAG_s banded matrix format.  I.e.,
c       A_{i,j} is stored in AB(kupper+klower+1+i-j,j)
c
c The indexing scheme I_m using is kind of complicated.  Not only is
c the i index in A_{i,j} compressed to handle just the banded part,
c but the j index is split into two parts, one for the species index
c and the other for the radial zone index.  One has to think through
c all of this very slowly and carefully...
c
c Regarding the banded matrix storage, one must be careful because
c LAPACK uses slightly different storage schemes depending on what
c will be done with the banded matrix.  For our purposes, which involves
c an LU decomposition, the proper layout is described in the source
c code for sgbtrf, or in the documentation for NAG (LAPACK is also
c distributed as part of NAG, as well as being available as freeware
c through http:www.netlib.org):

c  AB      (input/output) REAL array, dimension (LDAB,N)
c          On entry, the matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c  The band storage scheme is illustrated by the following example, when
c  M = N = 6, KL = 2, KU = 1:
c
c  On entry:                       On exit:
c
c      *    *    *    +    +    +       *    *    *   u14  u25  u36
c      *    *    +    +    +    +       *    *   u13  u24  u35  u46
c      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
c     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
c     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
c     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
c
c  Array elements marked * are not used by the routine; elements marked
c  + need not be set on entry, but are required by the routine to store
c  elements of U because of fill-in resulting from the row interchanges.
c
      klower=2*nfields-1        ! # of upper diagonals
      kupper=2*nfields-1        ! # of lower diagonals
      ldab=6*nfields-2          !includes extra storage for LAPACK
      nhi=ngrid-1
c
c  First the right-hand-side of Ax=b (where Tnew is used to store b):
c
      do k=1,nhi
         do i=1,nfields
            deltaT(i,k)=S(i,k)*dt
         enddo
      enddo
c
c boundary condition on deltaT
c
      k=ngrid
         do i=1,nfields
            deltaT(i,k)=0.D0
         enddo
c
c Initialize A to 0:
c
      do k=1,ngrid
         do i=1,nfields
            do ioff=1,ldab
               ab(ioff,i,k)=0.D0
            enddo
         enddo
      enddo
c
c d/dt term gives:
c
      rhop = 1.D0
      rhom = 1.D0
      do k=2,nhi
        do i=1,nfields
        do j=1,nfields
          ioff=kupper+klower+1+i-j
          ioffp=ioff-nfields
          ioffm=ioff+nfields
          ab(ioffp,j,k+1) = ab(ioffp,j,k+1) + am*vrho3(i,j,k)
          ab(ioff,j,k)    = ab(ioff,j,k) + ap*vrho3(i,j,k)
          ab(ioffm,j,k-1) = ab(ioffm,j,k-1) + am*vrho3(i,j,k)
        enddo
        enddo
      enddo
c boundary condition rhop*Tnew(0)=rhop*Tnew(1)
      k=1
      do i=1,nfields
      do j=1,nfields
         ioff=kupper+klower+1+i-j
         ioffp=ioff-nfields
         ab(ioffp,j,k+1) = ab(ioffp,j,k+1)+am*vrho3(i,j,k)
         ab(ioff,j,k)    = ab(ioff,j,k)+(ap+am)*vrho3(i,j,k)    
      enddo
      enddo
c boundary condition is Tnew(ngrid)=0 but need non-zero matrix element here
      k=ngrid
      do i=1,nfields
      do j=1,nfields
         ioff=kupper+klower+1+i-j
         ioffm=ioff+nfields
         ab(ioff,j,k)    = ab(ioff,j,k)+ap*vrho3(i,j,k)
         ab(ioffm,j,k-1) = ab(ioffm,j,k-1)+am*vrho3(i,j,k)
      enddo
      enddo
c
c nu_ij equilibration term gives:
c
      do k=2,nhi
         do i=1,nfields
            do j=1,nfields
               ioff=kupper+klower+1+i-j
               ioffp=ioff-nfields
               ioffm=ioff+nfields
               ab(ioffp,j,k+1)=ab(ioffp,j,k+1)+
     >          am*rhop*nu(i,j,k)*dt
               ab(ioff,j,k)=ab(ioff,j,k)+
     >          ap*rhom*nu(i,j,k)*dt
               ab(ioffm,j,k-1)=ab(ioffm,j,k-1)+
     >          am*rhom*nu(i,j,k)*dt
            enddo
         enddo
      enddo
c
c boundary condition is  nu(0)*Tnew(0)=nu(1)*Tnew(1)
c
      k=1
         do i=1,nfields
            do j=1,nfields
               ioff=kupper+klower+1+i-j
               ioffp=ioff-nfields
               ab(ioffp,j,k+1)=ab(ioffp,j,k+1)+
     >          am*rhop*nu(i,j,k)*dt
               ab(ioff,j,k)=ab(ioff,j,k)+
     >         (ap+am)*rhom*nu(i,j,k)*dt
            enddo
         enddo
c
c nu_p reaction term with gradient
c gradient assumed to vanish at k=1
c
      do k=2,nhi
         do i=2,nfields
            do j=1,nfields
               ioff=kupper+klower+1+i-j
               ioffp=ioff-nfields
               ioffm=ioff+nfields
               ab(ioffp,j,k+1)=ab(ioffp,j,k+1)
     >           + wp(3)*rhop*nu_p(i,j,k+1)*dt/dr(k,2)
     >           + wp(1)*rhop*nu_p(i,j,k)*dt/(dr(k,1)+dr(k,2))
               ab(ioff,j,k)=ab(ioff,j,k)
     >           - wp(3)*rhop*nu_p(i,j,k+1)*dt/dr(k,2)
     >           + wp(2)*rhom*nu_p(i,j,k)*dt/dr(k,1)
               ab(ioffm,j,k-1)=ab(ioffm,j,k-1)
     >           - wp(2)*rhom*nu_p(i,j,k)*dt/dr(k,1)
     >           - wp(1)*rhop*nu_p(i,j,k)*dt/(dr(k,1)+dr(k,2))
            enddo
         enddo
      enddo
c
c boundary condition is rhom*T(0)=rhop*T(1)
c
       k=1
         do i=1,nfields
            do j=1,nfields
               ioff=kupper+klower+1+i-j
               ioffp=ioff-nfields
               ab(ioffp,j,k+1)=ab(ioffp,j,k+1)
     >           + wp(3)*rhop*nu_p(i,j,k+1)*dt/dr(k,2)
     >           + wp(1)*rhop*nu_p(i,j,k)*dt/(dr(k,1)+dr(k,2))
               ab(ioff,j,k)=ab(ioff,j,k)
     >           - wp(3)*rhop*nu_p(i,j,k+1)*dt/dr(k,2)
     >           - wp(1)*rhop*nu_p(i,j,k)*dt/(dr(k,1)+dr(k,2))
            enddo
         enddo
c
c nu_2 reaction term 
c
      do k=2,nhi
         do i=2,nfields
            do j=1,nfields
               ioff=kupper+klower+1+i-j
               ioffp=ioff-nfields
               ioffm=ioff+nfields
               ab(ioffp,j,k+1)=ab(ioffp,j,k+1)
     >           + wp(3)*rhop*nu_2(i,j,k+1)*dt/2.D0
     >           + wp(1)*rhop*nu_2(i,j,k)*dt/2.D0
               ab(ioff,j,k)=ab(ioff,j,k)
     >           + wp(3)*rhom*nu_2(i,j,k+1)*dt/2.D0
     >           + wp(2)*rhop*nu_2(i,j,k)*dt/2.D0
               ab(ioffm,j,k-1)=ab(ioffm,j,k-1)
     >           + wp(2)*rhom*nu_2(i,j,k)*dt/2.D0
     >           + wp(1)*rhom*nu_2(i,j,k)*dt/2.D0
            enddo
         enddo
      enddo
c
c boundary condition is rhom*T(0)=rhop*T(1)
c
       k=1
         do i=1,nfields
            do j=1,nfields
               ioff=kupper+klower+1+i-j
               ioffp=ioff-nfields
               ab(ioffp,j,k+1)=ab(ioffp,j,k+1)
     >           + wp(3)*rhop*nu_2(i,j,k+1)*dt/2.D0
     >           +wp(1)* rhop*nu_2(i,j,k)*dt/2.D0
               ab(ioff,j,k)=ab(ioff,j,k)
     >           + wp(3)*rhop*nu_2(i,j,k+1)*dt/2.D0
     >           + wp(1)*rhop*nu_2(i,j,k)*dt/2.D0
     >           + wp(2)*rhop*nu_2(i,j,k)*dt
            enddo
         enddo
c
c Diffusion terms give:
c
      do k=2,nhi
         do i=1,nfields
            do j=1,nfields
               ioff=kupper+klower+1+i-j
               ioffp=ioff-nfields
               ioffm=ioff+nfields
               ab(ioffp,j,k+1)=ab(ioffp,j,k+1)-dt/vprime(k,1)/dr(k,1)*(
     &              vprime(k,2)*rhop*diff(i,j,k)/dr(k,2))
               ab(ioff,j,k)=ab(ioff,j,k)+dt/vprime(k,1)/dr(k,1)*(
     &              vprime(k,2)*rhom*diff(i,j,k)/dr(k,2)
     &              +vprime(k-1,2)*rhop*diff(i,j,k-1)/dr(k-1,2))
               ab(ioffm,j,k-1)=ab(ioffm,j,k-1)-dt/vprime(k,1)/dr(k,1)*(
     &              +vprime(k-1,2)*rhom*diff(i,j,k-1)/dr(k-1,2))
            enddo
        enddo
      enddo
c
c boundary condition at origin is that there is no flux through r=0,
c which can be thought of as vprime(0,2)=0:
c
      k=1
      do i=1,nfields
         do j=1,nfields
            ioff=kupper+klower+1+i-j
            ioffp=ioff-nfields
            ab(ioffp,j,k+1)=ab(ioffp,j,k+1)-dt/vprime(k,1)/dr(k,1)*(
     &           vprime(k,2)*rhop*diff(i,j,k)/dr(k,2))
            ab(ioff,j,k)=ab(ioff,j,k)+dt/vprime(k,1)/dr(k,1)*(
     &           vprime(k,2)*rhom*diff(i,j,k)/dr(k,2))
         enddo
      enddo
c
c Convection terms (no upwind differencing):
c
      rhop = 1.D0
      rhom = 1.D0
      do k=2,nhi
         do i=1,nfields
           do j=1,nfields
cgms added ioff,ioffp,ioffm these were missing in original code.
            ioff=kupper+klower+1+i-j
            ioffp=ioff-nfields
            ioffm=ioff+nfields
            up=0.5D0*conv3(i,j,k)
            um=0.5D0*conv3(i,j,k-1)
            ab(ioffp,j,k+1)=ab(ioffp,j,k+1)+dt/vprime(k,1)/dr(k,1)*(
     &          vprime(k,2)*up*rhop)
            ab(ioff,j,k)=ab(ioff,j,k)+dt/vprime(k,1)/dr(k,1)*(
     &           vprime(k,2)*up*rhom
     &          -vprime(k-1,2)*um*rhop)
            ab(ioffm,j,k-1)=ab(ioffm,j,k-1)-dt/vprime(k,1)/dr(k,1)*(
     &           vprime(k-1,2)*um*rhom)
           enddo
         enddo
      enddo
c
c Convection terms affecting the innermost zone (which has zero flux
c through the axis at its inner boundary):
c
      k=1
      do i=1,nfields
        do j=1,nfields
            ioff=kupper+klower+1+i-j
            ioffp=ioff-nfields
            up=0.5*conv3(i,j,k)
            ab(ioffp,j,k+1)=ab(ioffp,j,k+1)+dt/vprime(k,1)/dr(k,1)*(
     &          vprime(k,2)*up*rhop)
            ab(ioff,j,k)=ab(ioff,j,k)+dt/vprime(k,1)/dr(k,1)*(
     &          vprime(k,2)*up*rhom)
        enddo
      enddo
c
c Finally ready to call LAPACK routines to solve Ax=b:
c
      nsize=nfields*(ngrid)
c     call sgbtrf(nsize,nsize,klower,kupper,ab,ldab,nwork,ierror)
      call dgbtrf(nsize,nsize,klower,kupper,ab,ldab,nwork,ierror)
c (nag only has the double precision version, d.....)
      if(ierror .ne. 0) then
         write(6,*)'Fatal Error in vtrans calling sgbtrf, ierror=',
     &        ierror
         return
      endif
c     call sgbtrs('N',nsize,klower,kupper,1,ab,ldab,nwork,deltaT,nsize,
c    &     ierror)
      call dgbtrs('N',nsize,klower,kupper,1,ab,ldab,nwork,deltaT,nsize,
     &     ierror)
      do i=1,nfields
        do k=1,ngrid-1
          Tnew(i,k) = deltaT(i,k)
        enddo
      enddo
c
      return
      end                           
