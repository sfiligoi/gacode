!  gfortran -march=native -O3 -fno-stack-arrays -fimplicit-none -fdefault-real-8 landau.F90  -c
!intel:
!ifort -stand f15 -warn all -march=native -O3 -heap-arrays 10 -implicitnone -real-size 64 landau.F90 -c

#ifdef LANDAU_PREC16
!at least in the past couldn't use quad precision with frumpy pgi compiler
!for comparison purposes
#define HIREAL real(16)
#define HRS _16
#else
#define HIREAL real
#define HRS
#endif

module landau
  HIREAL, parameter,private :: pi1=atan(1.HRS)*4
  ! real, parameter :: intnorm=pi1**(-1.5)*4*pi1
  ! int gphys(v) d^3v=1=int g(x) x^2 dx *intnor
  ! Maybe one should normalise this once for a given xmax and density 1?
  real, parameter :: vunit=real(sqrt(2.HRS))
  ! v*vunit is v in units of sqrt(T/m)
  ! v itself is v in units of sqrt(2T/m) ("collision operator units")
  real, parameter :: intnorm=real(4/sqrt(pi1))
  real, parameter :: normcol=real(sqrt(8/pi1))
  ! Normalisation of ctest/emat cfield/emat for the self-collisions of the
  ! reference species without caa or cab.
  ! Note: must take into account muref^3=sqrt(1/8), since mu^2=m/(2T)
  real, parameter :: normcolcgyro=real(4/sqrt(pi1))

  integer :: verbose=0
contains
  real function tau_ab_inv(mua,mub,nb,za,zb,Ta,Tb)
    real,intent(in) :: mua,mub,nb,Za,Zb,Ta,Tb
    optional :: mua,mub,nb,Za,Zb,Ta,Tb
    real mua1,mub1,nb1,Za1,Zb1,Ta1,Tb1
    ! mua,mub are the factors such that exp(-(mua*v)**2) is the gaussian factor
    !     (relative to the reference mu, which may be the electrons or the
    !      ions, i.e., mua=mua(real)/muref)
    ! Ta,Tb are the Temperatures relative to the reference temperature.
    ! in fact  (tau=Ta/Tb)
    !  sqrt(8/pi)/3 * cab/(sqrt(ma Ta^3)) * nb * (1+ma/mb) * (1+ma/mb/tau)^1.5
    ! we return the ratio
    !        ***** tau_ab_inv/c_aaref*sqrt(m_ref)*Tref^1.5  *****
    !  tau_ab_inv taken from Eq. (47) in
    !      Hinton, "Collisional Transport in Plasma" (1983)
    ! c_aaref:= 2pi e_ref^4 ln(Lambda)
    ! c_ab/c_aaref=za^2*zb^2
    ! ma/maref=mua**2*Ta     ;   x:=v/(sqrt(2)*vref)
    ! ==> sqrt(m_ref)*Tref^1.5/[ sqrt(ma)*Ta^1.5 ] = 1/(mua*Ta**2)
    mua1=1.
    mub1=1.
    nb1=1.
    Za1=1.
    Zb1=1.
    Ta1=1.
    Tb1=1.
    if (present(Ta)) Ta1=Ta
    if (present(Tb)) Tb1=Tb
    if (present(Za)) Za1=Za
    if (present(Zb)) Zb1=Zb
    if (present(nb)) nb1=nb
    if (present(mua)) mua1=mua
    if (present(mub)) mub1=mub
    tau_ab_inv=real(sqrt(8./pi1))/3*za1**2*zb1**2*nb1/(mua1*Ta1**2)*&
         (1+mua1**2*Ta1/(mub1**2*Tb1))/(1+mua1**2/(mub1**2))**1.5
  end function tau_ab_inv
    
  subroutine gentestkernel(n,a1,b1,c1,xmax,beta,gp,gw,ngauss,lor_int,dif_int,addcutoff,t1t2_int)
    implicit none
    integer, intent(in) :: n,ngauss
    logical, intent(in) :: addcutoff
    real, intent(in) :: a1(n),b1(n),c1(n),& !poly rec-coeffis
         gp(ngauss),gw(ngauss),& ! gauss points and weights ngauss> n
         xmax,beta !integration range and relative factor in exponential
    real, dimension(n,n), intent(out) :: lor_int,dif_int,t1t2_int
    optional :: t1t2_int ! only calculate that term if it is present.
    optional :: addcutoff
    !allocate(lor_int(n,n),dif_int(n,n))
    real p0,p1,pi,q0,q1,qi,r0,r1,ri,s0,s1,si,x,xmax1
    real f1,f2
    integer i,j,k
    HIREAL xmaxx,x1
    logical t1t2,addcutoff1
    real t1,t2 ! for timing

    xmaxx=xmax
    t1t2=present(t1t2_int)
    addcutoff1=.true.
    if (present(addcutoff)) addcutoff1=addcutoff

2   format (A,3G25.16)
    
    if (verbose>2) then
       print 2,'1,lo,hi',fct1lo(xmaxx),fct1hi(xmaxx),fct1(xmaxx)
       print 2,'2,lo,hi',fct2lo(xmaxx),fct2hi(xmaxx),fct2(xmaxx)
    end if

    ! lor_int and dif_int are symmetric in i,j
    ! t1t2_int is not symmetric (but not antisymmetric.)
    ! The first index is in that case the output index

    ! regarding all the signs of the collision operator.
    ! we view it as being negative signed, because on the right hand side it
    ! reduces the dist.fct
    ! the Matrix element is thus
    ! -int g1 g2 (1-uhat uhat)/u (grad phi-grad psi)^2 dv1^3 dv2^3
    ! The fieldop=intkernel+deltaintegral should therefore be
    ! +2* int g1g2 grad phi grad psi (1-uhat uhat)/u dv1 dv2
    ! (1-uhat uhat)/u is positive and corresponds therefore
    !  - nabla nabla /laplace^2.
    ! I.e. the fieldop is
    ! -2* int g1 grad phi nabla nabla/lap^2 g2 grad psi  dv1
    ! I.e. the fieldop is then
    ! +2* int phi nabla g1 nabla nabla/lap^2 g2 nabla psi dv1
    ! This integral contains a positive deltaintegral, if the nablas are
    ! commuted towards the laplacian.
    ! i.e., intkernel calculates the *negativ* other part.

    ! But what should be the true normalisation?
    ! When the scattering partner has a more compact distribution function
    ! as indicated by an increase in beta, that means it is becoming more
    ! massive. The deflection rate of the lighter particle is hardly affected
    ! by this. The diffusion rate is however.
    ! fct1(w) is for large w proportional to 1/w. I.e. the overall
    ! normalisation should be proportional beta.
    ! fct2(w) is for large w proportional to 1/w as well. But diffusion should
    ! be proportional to 1/beta., so one beta less.
    ! Now if we look at the deflection rate for ions off electrons it should
    ! decrease for decreasing beta. And it does. But see, we have altered time.
    ! if we hold the ion speed fixed, and change the electron mass, the
    ! electrons become relatively faster by beta. They carry out more
    ! collisions and so on. So there is the additional factor beta.
    ! for intkernel it's the same. First you've got to add the beta**3 instead
    ! of the beta and then add the time shifting. Maybe one should look at the
    ! electrons as basis, or we can look at the middle ground, the geometrical
    ! mean of electron and ion speed.

    ! Now the test Operator for phi is
    ! + int g1 (grad phi)^2 nabla nabla/lap^2 g2  dv1
    ! This is very roughly
    !   int g1 (grad phi)^2 1/lap g2  dv1, which is dominantly negative.
    ! howver lor_int and dif_int in the following are dominantly positive.
    ! i.e. we take their negatives.
    ! for a scattering of a single species, we have to take also two times
    ! this test op. If scatter two species, take the two contributions.

    ! calculate test operator matrix elements
    ! concretely we calculate the two integrals
    ! lor_int=int_0...xmax exp(-x**2) l(l+1)p^2*fct1(v*beta)*beta^-2 dv
    ! and
    ! dif_int=int_0...xmax exp(-x**2) (p')^2 *fct2(v*beta)*beta^-3 dv
    ! where fct1(w)=-1/w ddw 1/lap_w^2 exp(-w^2)
    !       = 1/(8w^3)(exp(-w^2) w-(1-2w^2) I(w))
    ! and   fct2(w)=-w^2 d2dw2 1/lap_w^2 exp(-w^2)=
    !       =1/(4w)(-w exp(-w^2)+I(w))
    ! with I(w)=int_0^w exp(-x^2)dx=sqrt(pi)/2 erf(x)
    ! fct1=nulor*v^2 ; fct2=nudif
    !fct1(w)=.125*w**(-3)*(exp(-w**2)*w+(2*w**2-1)*.5*sqrt(pi1)*erf(w))
    !fct2(w)=.25*w**(-1)*(.5*sqrt(pi1)*erf(w)-w*exp(-w**2))
    ! ** problem: fct1 and fct2 suffer terrible cancellation near w=0.
    ! this reduces accuracy dramatically.
    !  can improve with series expansion.
    ! we calculate for the moment separately the two integrals.
    ! w== x1
    ! v==x
    ! pi=phi
    ! qi=phi'

    call cpu_time(t1)
    lor_int=0
    dif_int=0
    if (t1t2) t1t2_int=0
    if (addcutoff1 .and. beta>1) then
       xmax1=xmax/beta
    else
       xmax1=xmax
    end if
    do k=1,ngauss
       x=gp(k)*xmax1
       x1=x*beta
       if (addcutoff1) then
          f1=real(fct1lo(x1))
          f2=real(fct2lo(x1))
       else
          f1=real(fct1(x1))
          f2=real(fct2(x1))
       endif
       do j=1,n
          if (j==1) then
             q1=0
             qi=0
             p1=0
             pi=c1(1)*gw(k)*exp(-x**2)*xmax1
          else
             q0=qi
             qi=c1(j)*pi+(c1(j)*x-a1(j))*qi-b1(j)*q1
             q1=q0
             p0=pi
             pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
             p1=p0
          endif
          i=1
          r1=0
          ri=c1(1)
          s1=0
          si=0
          lor_int(j,i)=lor_int(j,i)+pi*ri*f1
          if (t1t2) t1t2_int(j,i)=t1t2_int(j,i)+2*qi*ri*real(x1)*f2
          ! ^^ vgl 1.8.19 (2)
          do i=2,n
             s0=si
             si=c1(i)*ri+(c1(i)*x-a1(i))*si-b1(i)*s1
             s1=s0
             r0=ri
             ri=(c1(i)*x-a1(i))*ri-b1(i)*r1
             r1=r0
             lor_int(j,i)=lor_int(j,i)+pi*ri*f1
             dif_int(j,i)=dif_int(j,i)+qi*si*f2
             if (t1t2) t1t2_int(j,i)=t1t2_int(j,i)+2*qi*ri*real(x1)*f2
          enddo
       enddo
    enddo
    if (addcutoff1 .and. beta>1) then
       do k=1,ngauss
          x1=xmax*(1+(beta-1)*gp(k))
          x=real(x1)/beta
          f1=real(fct1hi(x1))
          f2=real(fct2hi(x1))

          do j=1,n
             if (j==1) then
                q1=0
                qi=0
                p1=0
                pi=c1(1)*gw(k)*exp(-x**2)*xmax*(1-1/beta)
             else
                q0=qi
                qi=c1(j)*pi+(c1(j)*x-a1(j))*qi-b1(j)*q1
                q1=q0
                p0=pi
                pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
                p1=p0
             endif
             i=1
             r1=0
             ri=c1(1)
             s1=0
             si=0
             lor_int(j,i)=lor_int(j,i)+pi*ri*f1
             if (t1t2) t1t2_int(j,i)=t1t2_int(j,i)+2*qi*ri*real(x1)*f2
             ! ^^ vgl 1.8.19 (2)
             do i=2,n
                s0=si
                si=c1(i)*ri+(c1(i)*x-a1(i))*si-b1(i)*s1
                s1=s0
                r0=ri
                ri=(c1(i)*x-a1(i))*ri-b1(i)*r1
                r1=r0
                lor_int(j,i)=lor_int(j,i)+pi*ri*f1
                dif_int(j,i)=dif_int(j,i)+qi*si*f2
                if (t1t2) t1t2_int(j,i)=t1t2_int(j,i)+2*qi*ri*real(x1)*f2
             enddo
          enddo
       enddo
    end if

    !1/v ddv 1/lap_v^2 exp(-beta**2 x**2) *beta = beta**-2 fct1(beta*x) *beta
    ! = beta**-1 1/lap_w exp(-w^2)
    !w=beta*v, lap_w=beta**-2 *lap_v
    !1/v ddv 1/lap_v^2 propto v^2 =beta**-2 w^2 propto 1/w ddw 1/lap_w^2
    ! or other way:
    ! large v,w : 1/v ddv 1/lap_v^2 exp(-beta**2 x**2) propto 1/v 1/beta**3
    !   1/lap_w exp(-w^2) propto 1/w   ==> 1/v 1/beta^3 / (1/w) = 1/beta^2.
    !
    ! v^2 ddv^2 1/lap_v^2 exp(-beta**2 x**2) *beta= beta**2/w^2 beta**-2 /lap_w
    ! exp(-w) *beta =
    ! v^2 ddv^2 1/lap_v^2 propto v^4 = beta**-4 w^4 propto ..
    !                       =beta/w^2 /lap_w exp(-w^2)
    !may need to split the integral into a part up to xmax/beta
    ! and one up to xmax for large beta, due to singular behavior.
    lor_int=lor_int*beta**(-1) ! ***<--- must be multiplied with lphys(lphys+1)
    ! this is the only l dependence.
    dif_int=dif_int*beta**(-3)
    !*** Note: we have rescaled the above with beta, to account for the
    ! changing norm of g2=exp(-beta**2 v2**2)
    ! i.e. we multiply with beta!!
    if (t1t2) t1t2_int=t1t2_int*beta**(-4) ! ***<--- must be multiplied with (Ta/Tb-1)
    call cpu_time(t2)

    if (verbose>0) print '(A,G13.3)','gentestkernel took ',t2-t1
  contains
    elemental HIREAL function fct1(w)
      implicit none
      HIREAL,intent(in) :: w
      fct1=.125HRS*w**(-3)*(exp(-w*w)*w+(2*w*w-1)*.5HRS*sqrt(pi1)*erf(w))
    end function fct1
    elemental HIREAL function fct2(w)
      implicit none
      HIREAL,intent(in) :: w
      fct2=.25HRS*w**(-1)*(.5HRS*sqrt(pi1)*erf(w)-w*exp(-w**2))
    end function fct2
    elemental HIREAL function fct1lo(w)
      implicit none
      HIREAL,intent(in) :: w
      fct1lo=fct1(w)-exp(-xmaxx**2)/6  ! for w<xmax
    end function fct1lo
    elemental HIREAL function fct1hi(w)
      implicit none
      HIREAL,intent(in) :: w
      !  fct1hi(w)=w**(-3)*((1./24)*exp(-xmaxx**2)*xmaxx*(3-6*w**2+2*xmaxx**2)+(sqrt(pi1)/16)*(2*w**2-1)*erf(xmaxx))
      fct1hi=-w**(-3)/48*((-2*exp(-xmaxx**2)*xmaxx*(3+2*xmaxx**2)+3*sqrt(pi1)*erf(xmaxx))+&
           w**2*(12*exp(-xmaxx**2)*xmaxx-6*sqrt(pi1)*erf(xmaxx)))
      !Inaccuracy for
      !./hhv 5 .5 0 50 1 1
      !was due to inaccuracies in this routine, likely due to cancellations in the fcts.
    end function fct1hi
    elemental HIREAL function fct2lo(w)
      implicit none
      HIREAL,intent(in) :: w
      fct2lo=fct2(w)-w*w*exp(-xmaxx**2)/6 !for w<xmax
    end function fct2lo
    elemental HIREAL function fct2hi(w)
      implicit none
      HIREAL,intent(in) :: w
      fct2hi=-w**(-1)*((1.HRS/12)*xmaxx*exp(-xmaxx**2)*(3+2*xmaxx**2)-(sqrt(pi1)/8)*erf(xmaxx))
    end function fct2hi
  end subroutine gentestkernel

  subroutine genintkernel(n,lmax,a1,b1,c1,xmax,beta,t1t2ratio,gp,gw,ngauss,gp2,gw2,ng2,addcutoff,intkernel,intlokernel)
    implicit none
    optional :: intlokernel
    logical,intent(in) :: addcutoff
    integer, intent(in) :: n,lmax,ngauss,ng2  ! lmax=max l-**index**, lphys=l-1
    real, intent(in) :: a1(n),b1(n),c1(n),& !poly rec-coeffis
         gp(ngauss),gw(ngauss),& ! gauss points and weights ngauss> n
         gp2(ng2),gw2(ng2),& ! gauss points and weights for small kernel ng2~8
         xmax,beta,t1t2ratio !integration range and relative factor in exponential
    real, dimension(n,n,lmax), intent(out) :: intkernel(n,n,lmax),intlokernel(n,n,lmax) !intlokernel only used if t1t2ratio/=1
    real, parameter :: epsmax=1e-17 !epsmax is for cutoff of b1integral
    real, dimension(:,:), allocatable :: b1integral
    real, dimension(:,:,:), allocatable :: v2integral
    real, dimension(:,:,:,:), allocatable :: delta35integral,delta35cutoff
    real, dimension(:), allocatable ::  pimax ! only needed for addcutoff
    real mymin,mymax,v1max,v2max,x,p0,p1,pi
    real val,val1,val2
    real s12,s22
    real t1,t2 !for timing
    integer i,j,k,k1,l,lphys,m1,m2
    if (n==0 .or. lmax==0) return
    call cpu_time(t1)
    ! This routine calculates the main part of the field operator for the Background distribution functions
    ! g1=exp(-x^2) and g2=exp(-x^2*beta^2)*beta
    ! To get the integral for g1=exp(-x^2*alpha^2)*alpha and g2=exp(-x^2*beta^2)*beta with polynomials Pi(x*alpha) ,Pj(x*beta)
    ! one can use the fact that multiplication of all the velocity arguments with a factor f causes a change in the
    ! integral by a factor f^-3*f^2=1/f.
    ! Therefore fcode(alpha,beta)=fcode(1,beta/alpha)/alpha=f(alpha/beta,1)/beta.
    ! So fij(1,lambda)=fji(1,1/lambda)/lambda.
    ! Otherwise the actual field operator is
    ! Fij=-cab*32/pi *Fcode_complete (1,beta/alpha)/(beta^1.5*mu^2.5) ; alpha=sqrt(ma) beta=sqrt(mb)
    ! The same normalisation applies to the Testoperator.
    ! Fij_code_complete(1,lambda)=Fij_c(1,lambda)+Fji_c(1,1/lambda)

    ! Ok, now for m1, m2 and one p_i and one p_j and one beta calculate
    ! #=int v1=0..xmax v2=0..xmax p_i(v1) p_j(v2) exp(-v1**2 -v2**2) v1^m1 v2^m2 theta(beta*v1-v2) dv1 dv2.
    ! for this to be valid, m2>0, and m1+m2>-2. But we have assurance that m1+m2>=3

    ! In the following this is rescaled with a prefactor. to determine it, let
    ! us use w2=v2/beta v2=w2*beta as integration variable. We obtain for #:
    ! #=int v1=0^xmax w2=0^xmax/beta p_i(v1) p_j(w2 beta ) *
    !     *exp(-v1**2 -beta**2 w2**2) * v1^m1 w2^m2 beta^m2 *
    !     *theta(v1-w2) dv1  dw2 *beta
    ! ***Important: We rescale this with the factor beta^(-m2) . Then we get
    ! something similar to the quantities in NEO
    ! #rescaled =
    !    int v1=0..xmax v2=0..xmax p_i(v1) p_j(v2) exp(-v1**2 -v2**2) v1^m1 v2^m2 beta^-m2 theta(beta*v1-v2) dv1 dv2  

    ! compare the following with the expression in deltaintegrate.
    ! intkernel_rescaled(i,j,beta)=
    !   =int v1=0^xmax w2=0^xmax/beta p_i(v1) p_j(w2 beta ) *
    !     *exp(-v1**2 -beta**2 w2**2) * v1^m1 w2^m2 * theta(v1-w2) dv1  dw2 beta
    ! This is for i=1,j=1 equal to jeffs integral F(beta^2,m1,m2)*beta

    ! In the end delta35integral(j,i,idx,l)=
    ! int v1=0^xmax v2=0^xmax p_i(v1) p_j(v2) exp(-v1**2 -v2**2) v1^m1 (v2/beta)^m2 theta(beta*v1-v2) dv1 dv2
    ! with m1=1-l m2=l+2 for idx=1 and m2=l+4 for idx=2 (l: physical l = idx_l-1)
    !
    ! To get the other half by a call to this routine with 1/beta:
    ! 2-half rescaled=int v1=0^xmax w2=0^xmax/beta p_i(v1) p_j(w2 beta ) *
    !     *exp(-v1**2 -beta**2 w2**2) * v1^m2 w2^m1 * theta(w2-v1) dv1  dw2 *beta
    ! again with m2>0 and m1 maybe <0
    ! we reintroduce v2=w2*beta and replace v1=w1/beta <-> w1=v1*beta
    ! 2-half resc=int w1=0^xmax*beta w2=0^xmax p_i(w1/beta) p_j(v2) *
    !      * exp(-v1**2 beta**-2 -v2**2) * w1^m2/beta^m2 v2^m1/beta^m1 *
    !      *theta(v2-w1) dw1/beta dv2/beta *beta
    !   =int w1=0^xmax*beta w2=0^xmax p_i(w1/beta) p_j(v2) *
    !      * exp(-v1**2 beta**-2 -v2**2) * w1^m2 v2^m1 *
    !      *theta(v2-w1) dw1 dv2/ beta^(m1+m2+1)
    ! =intkernel_rescaled(j,i,1/beta)/beta^(m1+m2)
    ! Note m1+m2 is either 3 or 5
    ! Note however, that within the intkernel there is a factor 1/(s12*s22) propto
    !  beta**2 for the terms with m1+m2=5 and 1/s12 =const or 1/s22 prop beta**2
    !  for m1+m2=3.
    !  the factor 1/(s12*s22) flips from beta**2 to beta**-2
    ! ==> taking this into account we get
    ! 2-half resc=intkernel_rescaled(j,i,1/beta)/beta^(5)*beta**4
    ! 2-half resc=intkernel_rescaled(j,i,1/beta)/beta
    ! in addition



    ! beta>1 ==> split integral over v1 into
    !   A= int dv1_0^xmax/beta  dv2_0^xmax p_i(v1) p_j(v2) exp(-v1**2 -v2**2) v1^m1 v2^m2 theta(beta*v1-v2) dv1 dv2
    !  and  B= int_xmax/beta..xmax p_i(v1) v1^m1*beta^-m2 exp(-v1**2) dx * int_0^xmax p_j(v2) v2^m2 exp(-v2**2) dx
    !  in general the v1 integral in B is with m1<0 and needs to be cut off at (xmax/beta)*(1e16)^[1/(1-|m1|)]
    !  then it contains only about 64 polynomials or so.
    !  The 2nd part of the B-integral is trivial more or less, since we have
    !      the steen vertices? Yes.
    !  For the A integral we carry out the v2 integral for all v1 vertices, that we need for a nice Gauss rule.
    ! beta<1 ==> can replace the v2-boundary xmax by beta*xmax, and need only
    !       compute values of v2-integral for beta* the Gauss points of the v1-integral.
    !   ***** in case beta<1 we rescale the v2 integral components with beta^(-m2)
    ! for stepping between.
    ! for each l=0...lmax, we need (m1,m2)=(1-l,l+4) (delta=5)
    !   (m1,m2)=(1-l,l+2) (delta=3)
    !   (m1,m2)=(3-l,l+2) (delta=5)  (only for l>=2).
    ! i.e. we need delta=3 and delta=5.
    ! Leading coeffi for delta=3 m2>=2 for delta=5 m2>=4 ,i.e.
    ! m1=1,0,-1,...-(lmax-1), i.e. lmax+1 values. Lmax=nxi-1 ==>nxi values.
    ! m2=2..lmax+4 ==> lmax+3 values.
    if (beta>1) then
       !calculate first half of B-integrals for all n orthonormal polynomials
       !*** we originally rescaled b1integral by beta**m1
       !    now we rescale by mymin**(-m1)=(xmax/beta)**(-m1)=
       !       =(beta/xmax)**m1   to avoid overflow
       allocate(b1integral(n,lmax+1))
       b1integral=0
       mymax=xmax
       mymin=xmax/beta
       do l=1,lmax+1
          m1=2-l
          if (m1<-1) mymax=min(xmax,mymin*epsmax**(1./(m1+1)))
          if (verbose>2) print '(I3,A,3I3)',l,'ngauss',ngauss,mymin,mymax
          do k=1,ngauss
             x=gp(k)*(mymax-mymin)+mymin
             p1=0
             pi=c1(1)*(gw(k)*exp(-x**2)*(x/mymin)**m1)*(mymax-mymin) ! ** <--rescaling of v1 here
             ! here all constant factors.
             j=1
             b1integral(j,l)=b1integral(j,l)+pi
             do j=2,n
                p0=pi
                pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
                p1=p0
                b1integral(j,l)=b1integral(j,l)+pi
             enddo
          enddo
       enddo
!!$2      format (*(G0,"  "))
!!$         if (verbose>2) then
!!$            print 2,'Check b1integral:'
!!$            print 2,'n=1 (Pn=const)'
!!$            do l=1,lmax+1
!!$               m1=2-l
!!$               val=c1(1)*.5*(gsl_sf_gamma_inc(.5*(1+m1),mymin**2)&
!!$                    &-gsl_sf_gamma_inc(.5*(1+m1),xmax**2))
!!$               if (beta>1) val=val*(beta/xmax)**m1 ! ** rescaling of v1
!!$               print 2,'l=',l,'m1=',m1,b1integral(1,l),val,val/b1integral(1,l)-1
!!$            enddo
!!$            print 2,'n=2 (Pn=const)'
!!$            do l=1,lmax+1
!!$               m1=2-l
!!$               val1=c1(1)*.5*(gsl_sf_gamma_inc(.5*(1+m1),mymin**2)&
!!$                    &-gsl_sf_gamma_inc(.5*(1+m1),xmax**2))
!!$               val2=c1(2)*c1(1)*.5*(gsl_sf_gamma_inc(.5*(1+m1+1),mymin**2)&
!!$                    &-gsl_sf_gamma_inc(.5*(1+m1+1),xmax**2))
!!$               val=val2-a1(2)*val1
!!$               if (beta>1) val=val*(beta/xmax)**m1 ! ** rescaling of v1
!!$               print 2,'l=',l,'m1=',m1,b1integral(2,l),val,val/b1integral(2,l)-1
!!$            enddo
!!$         endif
    endif
    !chacka!
    ! now calculate all necessary v2 integrals for the gausspoints of v1
    ! we do this with a secondary gauss rule with lower point number in between
    ! the v1 gausspoints. ng2 is relatively small, maybe 8 or so. hope this'
    ! nuff.
    ! The last igauss index is for the full v2max integral.
    ! we rescale each v2integral(n,l,ig) by 1/(v2max*gp(ig))**m2
    !  = 1/(xmax*min(1,beta)*gp(ig))**m2
    ! instead of rescaling all at once by beta**-m2 if beta<0.

    allocate(v2integral(n,lmax+3,ngauss+1))
    v2integral=0
    v2max=xmax*min(1.,beta)
    mymin=0.
    do k=1,ngauss+1
       do l=1,lmax+3
          m2=l+1
          if (k<=ngauss) then
             mymax=gp(k)*v2max
          else
             mymax=v2max
          endif
          ! Correctly scale down from the earlier gausspoint:
          if (k>1) v2integral(:,l,k)=v2integral(:,l,k-1)*(mymin/mymax)**m2

          do k1=1,ng2
             x=(mymax-mymin)*gp2(k1)+mymin
             p1=0
             ! rescaling by mymax**(-m2)=1/(xmax*min(1,beta)*gp(ig))**m2 :
             pi=c1(1)*(gw2(k1)*exp(-x**2)*(x/mymax)**m2)*(mymax-mymin)

             j=1
             v2integral(j,l,k)=v2integral(j,l,k)+pi
             do j=2,n
                p0=pi
                pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
                p1=p0
                v2integral(j,l,k)=v2integral(j,l,k)+pi
             enddo
          enddo
       enddo
       mymin=mymax
    enddo
    if (.true. .and. verbose>2) then
       print '(A,G23.16)','Check complete v2integral up to v2=',v2max
       do i=1,min(n,lmax+3)
          l=i
          m2=l+1
          val=0
          do k=1,ngauss
             x=v2max*gp(k)
             p1=0
             pi=c1(1)*(gw(k)*exp(-x**2)*(x/v2max)**m2)*v2max
             j=1
             do j=2,i
                p0=pi
                pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
                p1=p0
             enddo
             val=val+pi
          enddo

4         format (A,2I3,5(G25.16))
          print 4,'v2int',i,l,v2integral(i,l,ngauss+1),val,v2integral(i,l,ngauss+1)/val-1,&
               maxval(abs(v2integral(i,l,:)))/v2integral(i,l,ngauss+1)
          ! agreement here perfect for beta=0.02 e.g.
          ! agreement here less than perfect for n=50 beta=1.
          ! even though the integral is *not* zero??
          ! well the integral somehow must be practically zero exactly for beta=1.
          ! but not otherwise.
          ! however this is probably ok because the results should be practically zero
          ! because we have practical orthogonality of the higher polynomials
       enddo
    endif
    ! Compare end result for beta>=1
    if (verbose>2 .and. beta >=1 .and. .false.) then
       do l=1,lmax+3,2
          print '(G23.16)',v2max
          m2=l+1
          do i=2,n,2
             val2=0
             do k=1,ngauss
                x=v2max*gp(k)
                p1=0
                pi=c1(1)
                j=1
                do j=2,i
                   p0=pi
                   pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
                   p1=p0
                enddo
                val2=val2+gw(k)*pi*exp(-x**2)*(x/v2max)**m2*v2max
             enddo

             print 4,'v2compare',i,l,v2integral(i,l,ngauss+1),val2,v2integral(i,l&
                  &,ngauss+1)/val2-1
             !so for i>m2+1=l+3 this value should be zero.
             !looks good.
          enddo
       enddo
    endif

    !now calculate all the v1-integrals and later add the bintegrals.
    allocate(delta35integral(n,n,2,lmax+1))
    delta35integral=0
    ! delta35integral(v2-index,v1-index,1 or 2,lindex)
    !index 1 is with delta=3 and index 2 with delta=5
    v1max=xmax*min(1.,1./beta)  ! in case beta>0 we must add later the missing
    ! piece.
    do k=1,ngauss
       x=v1max*gp(k)
       do l=1,lmax+1
          m1=2-l
          m2=l+1
          p1=0
          ! we have rescaled v2integral(i,l,ig) by a factor
          !   1/(v2max*gp(ig))**m2, wereby m2=l+1
          !  = 1/(xmax*min(1,beta)*gp(ig))**m2
          ! note, that v2max= beta*v1max
          ! now we must multiply each (among other things by
          ! x**m1 = (v1max*gp(ig))**m1
          ! i.e. to cancel the scaling and to multiply at the same time
          ! we multiply with (v1max*gp(ig))**m1*(v2max*gp(ig))**m2
          ! =(v1max*gp(ig))**(m1+m2) *beta**m2
          !*** now at this point we apply the scaling beta **-m2
          ! so we multiply only with
          ! (v1max*gp(ig))**(m1+m2)

          pi=c1(1)*(gw(k)*exp(-x**2)*x**(m1+m2))*v1max !m1+m2=3

          ! hmm we may get an overflow here for very small x.
          ! for too many gausspoints.
          ! of course this is then cancelled by the v2integral
          ! but that in turn may have had an underflow.
          ! so you want to rescale the v2integral by it's boundary
          ! x. i.e. divide it by x**m2.
          i=1
          delta35integral(:,i,1,l)=delta35integral(:,i,1,l)+pi*v2integral(:,l,k)
          delta35integral(:,i,2,l)=delta35integral(:,i,2,l)+pi*v2integral(:,l+2,k)*x**2 ! rescaling since m2 is here +2.
          do i=2,n
             p0=pi
             pi=(c1(i)*x-a1(i))*pi-b1(i)*p1
             p1=p0
             delta35integral(:,i,1,l)=delta35integral(:,i,1,l)+pi*v2integral(:,l,k)
             delta35integral(:,i,2,l)=delta35integral(:,i,2,l)+pi*v2integral(:&
                  &,l+2,k)*x**2 ! rescaling since m2 is here +2.
          enddo
       enddo
    enddo

    ! it seems that for high l there is massive cancellation between
    ! hmm this cannot be the only problem for beta!=1.
    ! for beta=1 everything is fine for all l.
    ! b1integral and delta35integral.
    !print 2,'!!!warning leaving out bintegral'
    if (verbose>2 .and. beta>1) then
       k=ngauss+1
       do l=1,n
5         format (A,I3,5(G25.16))
         print 5,'d35l',l,delta35integral(l,l,1,l),b1integral(l,l)&
               &*v2integral(l,l,k)*v1max**3,delta35integral(l,l,1,l)+b1integral(l,l)&
               &*v2integral(l,l,k)*v1max**3
          ! for testing, e.g. do a
          !./hhv 5 10 0|grep d35l
          !and
          !./hhv 5 20 0|grep d35l
          ! whereby the last column does not change,
          ! but the first two do, because of different integral splits.
       enddo
       do l=1,n
          print 5,'d35l',l,delta35integral(l,l,2,l),b1integral(l,l)&
               &*v2integral(l,l+2,k)*v1max**5,delta35integral(l,l,2,l)+b1integral(l,l)&
               &*v2integral(l,l+2,k)*v1max**5
       enddo
    endif
    if (beta>1) then
       ! add back missing integral:
       do l=1,lmax+1
          m1=2-l
          m2=l+1
          do i=1,n
             ! b1integral(i,l) has been rescaled by (beta/xmax)**m1 m1=2-l
             ! we have rescaled v2integral(i,l,ig) by a factor
             !   1/(v2max*gp(ig))**m2
             !  = 1/(xmax*min(1,beta)*gp(ig))**m2, wereby m2=l+1
             ! note, that v2max= beta*v1max always
             ! v2integral(i,l,ngauss+1) has been rescaled by 1/v2max**m2
             ! because that is the full integral.
             ! since in this section beta>1 v2max=xmax, v1max=xmax/beta
             ! in toto we have to rescale by (v1max)**m1*(v2max)**m2
             ! = v1xmax**(m1+m2) beta**(m2)
             ! now at this point we apply the rescaling beta**-m2
             ! so we multiply only with v1max**(m1+m2)
             delta35integral(:,i,1,l)=delta35integral(:,i,1,l)+b1integral(i,l)&
                  &*v2integral(:,l,ngauss+1)*v1max**(m1+m2)  !m1+m2=3
             delta35integral(:,i,2,l)=delta35integral(:,i,2,l)+b1integral(i,l)&
                  &*v2integral(:,l+2,ngauss+1)*v1max**(m1+m2+2) !m1+m2+2=5
          enddo
       enddo
    endif
    ! if (beta > 1) then !** need to correct scaling by beta**(-m1-m2)
    !! because we rescaled the v1 component integrals above for this case.
    !! for delta=3: -m1-m2= -3 for delta=5: -m1-m2=-5
    !   delta35integral(:,:,1,:)=delta35integral(:,:,1,:)*beta**(-3)
    !  delta35integral(:,:,2,:)=delta35integral(:,:,2,:)*beta**(-5)
    !endif
    if (verbose>2) then
       do l=1,lmax+1
          print 5,'lend',l,delta35integral(1,1,1,l),delta35integral(1,1,1,l)/c1(1)**2,&
               delta35integral(1,1,2,l),delta35integral(1,1,2,l)/c1(1)**2
       enddo
    endif

    if (addcutoff) then
       ! adding Ia,Ib,Iab to delta35integral(:,:,2,:)
       ! j=2...lmax+4 B(text)(max,j)=v2integral(:,j-1,ngauss+1)
       ! j=-1,...lmax-1 AN(text)(j)=b1integral(:,j+2)
       ! j=-1..lmax-1 I>(text)(j,j+3/5)=delta35integral(:,:,1/2,j+2)
       allocate(pimax(n),delta35cutoff(n,n,3,lmax+1))
       pi=c1(1)*exp(-xmax**2)
       pimax(1)=pi
       p1=0
       do i=2,n
          p0=pi
          pi=(c1(i)*xmax-a1(i))*pi-b1(i)*p1
          p1=p0
          pimax(i)=pi
       end do
       do l=1,lmax+1
          m1=l-2
          m2=m1+5 !l+3
          do i=1,n
             do j=1,n
                !delta35integral(j,i,2,l)=I>(l-2,l+3) for Pi-1,Pj-1
                !v2integral(j,l+2,ngauss+1)=B(v1max,l+3) for Pj-1
                !b1integral(i,l)=AN(l-2)
                !==>
                !I>a(l-2,l+3)=B(v1max,l+3)xmax^(l+3-(l-2)-1) psia
                !           = .... xmax^(5-1) psia
                !I>b(l-2,l+3)=An(l-2)v1max^(5-1)phib
                !I>ab(l-2,l+3)=psia phib*(1/beta)^(l+3-1)*xmax^(5-2)
                !mua=1/2
                delta35cutoff(j,i,1,l)=&
                     v2integral(j,l+2,ngauss+1)*pimax(i)*xmax**(5-1)*(v1max/xmax)**(l+3)
                !pimax(j)=psia(xmax)
                !I>a=B(max,m2)vamax^(m2-m1-1) psia(xmax)
                if (beta>1) then
                   !phia(xmax/beta)=pimax(i)
                   !I>b=AN(m1)phia(xmax/beta)vbmax^(m2-m1-1)
                   !mub=1/(2beta)
                   delta35cutoff(j,i,2,l)=&
                        b1integral(i,l)*pimax(j)*v1max**(5-1)*beta
                   delta35cutoff(j,i,3,l)=pimax(i)*pimax(j)*beta**(-l-2)*xmax**(5-2)*beta
                else
                   delta35cutoff(j,i,2,l)=0
                   if (beta==1) then
                      delta35cutoff(j,i,3,l)=.5*pimax(i)*pimax(j)*beta**(-l-2)*xmax**(5-2)*beta
                   else
                      delta35cutoff(j,i,3,l)=0
                   endif
                endif
             end do
          end do
       end do
    end if

    ! now adding the three kernels together
    ! note: in "lowerkernel in gl2check1.nb and in the writeups
    ! we have used v1<v2. but here v2<v1 ==> s12 and s22
    ! are interchanged.
    ! I.e. we calculate here the "hikernel" of gl2check1.nb
    ! and: sigma1^2=.5, sigma2^2=.5/beta^2

    ! this does of course not matter if beta=1. and t1=t2
    s12=.5
    s22=.5/beta**2
    !allocate(intkernel(n,n,lmax))
    do l=1,lmax
       lphys=l-1 ! difference between index and actual l
       intkernel(:,:,l)=-(lphys+1.)*(lphys+2)/(2*(1+2*lphys)*(3+&
            &2*lphys)*s12*s22)*delta35integral(:,:,2,l)
       if (lphys>=2) then
          intkernel(:,:,l)=intkernel(:,:,l)+ 1./(1+2*lphys)* lphys*(lphys-1.)&
               /(2*(2*lphys-1.)*s12*s22)* delta35integral(:,:,2,l-2)
       endif
       if (t1t2ratio/=1) then
          intlokernel(:,:,l)=intkernel(:,:,l)*t1t2ratio+1./(1+2*lphys)*&
               ((1+lphys)/(s12/t1t2ratio)-lphys/(s22))*delta35integral(:,:,1,l)
          if (addcutoff) then
             intlokernel(:,:,l)=intlokernel(:,:,l)-(lphys+1.)*(lphys+2)/(2*(1+2*lphys)*(3+&
                  &2*lphys))*(delta35cutoff(:,:,1,l)/s22&
                  +delta35cutoff(:,:,2,l)/s12*t1t2ratio&
                  +delta35cutoff(:,:,3,l))
             if (lphys>=2) then
                intlokernel(:,:,l)=intlokernel(:,:,l)+ 1./(1+2*lphys)* lphys*(lphys-1.)&
                     /(2*(2*lphys-1.))*(delta35cutoff(:,:,1,l-2)/s22&
                     +delta35cutoff(:,:,2,l-2)/s12*t1t2ratio&
                     +delta35cutoff(:,:,3,l-2))
             end if

             do i=1,n
                intlokernel(:,i,l)=intlokernel(:,i,l)+1./(1+2*lphys)*&
                     ((1+lphys)*v2integral(:,l,ngauss+1)*pimax(i)*xmax**(3-1)*(v1max/xmax)**(l+1))
                if (beta>1) &
                     intlokernel(:,i,l)=intlokernel(:,i,l)+1./(1+2*lphys)*(&
                     -lphys*b1integral(i,l)*pimax(:)*(v1max)**(3-1))  *beta
             enddo
          endif
       endif
       intkernel(:,:,l)=intkernel(:,:,l)/t1t2ratio+1./(1+2*lphys)*&
            ((1+lphys)/(s12)-lphys/(s22*t1t2ratio))*delta35integral(:,:,1,l)
       if (addcutoff) then
          intkernel(:,:,l)=intkernel(:,:,l)-(lphys+1.)*(lphys+2)/(2*(1+2*lphys)*(3+&
               &2*lphys))*(delta35cutoff(:,:,1,l)/(s22*t1t2ratio)&
               +delta35cutoff(:,:,2,l)/s12&
               +delta35cutoff(:,:,3,l))
          if (lphys>=2) then
             intkernel(:,:,l)=intkernel(:,:,l)+ 1./(1+2*lphys)* lphys*(lphys-1.)&
                  /(2*(2*lphys-1.))*(delta35cutoff(:,:,1,l-2)/(s22*t1t2ratio)&
                  +delta35cutoff(:,:,2,l-2)/s12&
                  +delta35cutoff(:,:,3,l-2))
          end if
          do i=1,n
             intkernel(:,i,l)=intkernel(:,i,l)+1./(1+2*lphys)*&
                  ((1+lphys)*v2integral(:,l,ngauss+1)*pimax(i)*xmax**(3-1)*(v1max/xmax)**(l+1))
             if (beta>1) &
                  intkernel(:,i,l)=intkernel(:,i,l)+1./(1+2*lphys)*(&
                  -lphys*b1integral(i,l)*pimax(:)*(v1max)**(3-1))  *beta
          enddo
       endif
    enddo
    !search for near cancellations:
    if (verbose>2) then
       do l=1,lmax
          lphys=l-1
          do i=1,n
             do j=1,n
                val=intkernel(i,j,l)
                val1=-(lphys+1.)*(lphys+2)/(2*(1+2*lphys)*(3+&
                     &2*lphys)*s12*s22)*delta35integral(i,j,2,l)
                val2=1./(1+2*lphys)*&
                     &((1+lphys)/s22-lphys/s12)*delta35integral(i,j,1,l)
                if (abs(val)<abs(val1)*1e-5 .or. abs(val)<abs(val2)*1e-5) then
6                  format (A,3I3,5(G25.16))
                   print 6,'nc',i,j,l,val,val1,val2
                endif
             end do
          end do
       end do

       do l=1,min(3,lmax)
          lphys=l-1 ! difference between index and actual l
          do i=1,5
             print 4,'ik',l,i,intkernel(i,1:5,l)
             print 4,'ik1',l,i,-(lphys+1.)*(lphys+2.)/(2.*(1.+2*lphys)*(3.+2*lphys)*s12&
                  &*s22)&
                  &*delta35integral(i,1:5,2,l)
             print 4,'ik2',l,i,1./(1+2*lphys)*&
                  &((1+lphys)/s22-lphys/s12)*delta35integral(i,1:5,1,l)
             if (lphys>=2) then
                print 4,'ik3',l,i,1./(1+2*lphys)*&
                     &lphys*(lphys-1.)/(2*(2*lphys-1.)*s12*s22)*&
                     &delta35integral(i,1:5,2,l-2)
             endif
          enddo
       enddo
    endif
    !ts unfortunately, we have calculated the transposed version
    ! not easy to turn around.
    ! hence:

    do l=1,lmax
#ifdef DIMATCOPY
#ifdef __INTEL_COMPILER
       !use MKL
       call mkl_dimatcopy('c','t',n,n,1.,intkernel(:,:,l),n,n)
       if (t1t2ratio/=1) call mkl_dimatcopy('c','t',n,n,1.,intlokernel(:,:,l),n,n)
#else
       !use openBLAS
       call dimatcopy('c','t',n,n,1.,intkernel(:,:,l),n,n)
       if (t1t2ratio/=1) call dimatcopy('c','t',n,n,1.,intlokernel(:,:,l),n,n)
#endif
#else
         !alternative to dimatcopy
         do i=1,n
            do j=1,i-1
               val=intkernel(i,j,l)
               intkernel(i,j,l)=intkernel(j,i,l)
               intkernel(j,i,l)=val
               if (t1t2ratio/=1) then
                  val=intkernel(i,j,l)
                  intlokernel(i,j,l)=intlokernel(j,i,l)
                  intlokernel(j,i,l)=val
               endif
            enddo
         enddo
#endif
    enddo
    !** rescale everything with beta**2, since we did not completely balance
    !the shrinking of the integral of the exp(-beta**2 x**2) in **3d-space**
    !(only 1d)
    !  intkernel=intkernel*beta**2
    call cpu_time(t2)
    if (verbose>0) print '(A,G13.3)','genintkernel took',t2-t1
  end subroutine genintkernel

  subroutine deltaintegrate(n,a1,b1,c1,xmax,beta,gp,gw,ngauss,deltaintegral)
    ! ??? This should be recoded with Steen points not Gauleg points.
    implicit none
    integer, intent(in) :: n,ngauss
    real, intent(in) :: a1(n),b1(n),c1(n),& !poly rec-coeffis
         gp(ngauss),gw(ngauss),& ! gauss points and weights
         xmax,beta !integration range and relative factor in exponential
    real, intent(out) :: deltaintegral(n,n)
    !allocate(deltaintegral(n,n))

    real v2max,x,x1,p0,p1,pi,q0,q1,qi
    integer i,j,k
    real t1,t2 ! timing
    ! Now for the last integral:  (v2=w2*beta, w2=v2/beta)
    ! deltaintegral(k,j,beta):=
    ! int pk(v1) pj(w2*beta) exp(-v1**2-beta**2*w2**2) delta(v1-w2) v1**2 *beta dv1 dw2
    ! = int pk(v1) pj(v2) exp(-v1**2-v2**2) delta(v1-v2/beta) v1**2 dv1 dv2
    ! =int pk(v2/beta) pj(v2) exp(-(1+beta^-2) v2**2) v2**2 dv2 *beta**(-2)
    !   v2<xmax,xmax*beta
    ! now replace again w2=v2/beta
    ! =int pk(w2) pj(v2) exp(-(1+beta^-2) v2**2) w2**2 dv2
    ! =int pk(w2) pj(w2*beta) exp(-(1+beta^2) w2**2) w2**2 dw2*beta
    ! =deltaintegral(j,k,1/beta)/beta
    ! ** in the rescaling we did not take into account, that the Gaussian
    ! is a *3D* function and should be rescaled by beta**3 and not just beta
    ! I.e. we rescale di(i,k,beta)=beta**2 diold(i,k,beta)
    ! i.e. di(i,k,beta)=beta**2 diold(k,i,1/beta)/beta
    ! and diold(k,i,1/beta)*1/beta**2=di(k,i,1/beta)
    !     diold(k,i,1/beta)=beta**2 di(k,i,1/beta)
    ! ==> di(i,k,beta)=beta**3 di(k,i,1/beta)
    v2max=xmax*min(1.,beta) ! same as already defined before

    call cpu_time(t1)

    deltaintegral=0
    do i=1,ngauss
       x=gp(i)*v2max  ! vertex for v2
       x1=x/beta      ! vertex for v1==w2
       do j=1,n
          if (j==1) then
             p1=0
             pi=c1(1)*gw(i)*exp(-(1+beta**(-2))*x**2)*x1**2*v2max
          else
             p0=pi
             pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
             p1=p0
          endif
          q1=0
          qi=c1(1)*pi
          k=1
          deltaintegral(k,j)=deltaintegral(k,j)+qi
          do k=2,n
             q0=qi
             qi=(c1(k)*x1-a1(k))*qi-b1(k)*q1
             q1=q0
             deltaintegral(k,j)=deltaintegral(k,j)+qi
          enddo
       enddo
    enddo
    call cpu_time(t2)
    if (verbose>0) print '(A,G13.3)','deltaintegrate took',t2-t1

    if (verbose>2) then
5      format (A,I3,5(G25.16))
       do i=1,min(n,6)
          print 5,'di',i,deltaintegral(i,:min(6,n))
       enddo
       do i=max(1,n-4),n
          print 5,'di',i,deltaintegral(i,max(1,n-4):)
       enddo
       do i=1,min(n,6)
          print 5,'di*beta',i,deltaintegral(i,:min(6,n))*beta
       enddo
       do i=max(1,n-4),n
          print 5,'di*beta',i,deltaintegral(i,max(1,n-4):)*beta
       enddo
    endif
    ! new rescaling:
    !  deltaintegral=deltaintegral*beta**2
  end subroutine deltaintegrate

  subroutine genenergymatrix(n,a1,b1,c1,xmax,gp,gw,ngauss,deltaintegral)
    ! calculates int v^2 exp(-v^2) P_i(v)P_j(v) dv.
    ! adapted from deltaintegrate.
    implicit none
    integer, intent(in) :: n,ngauss
    real, intent(in) :: a1(n),b1(n),c1(n),& !poly rec-coeffis
         gp(ngauss),gw(ngauss),& ! gauss points and weights
         xmax !integration range and relative factor in exponential
    real, intent(out) :: deltaintegral(n,n)
    !allocate(deltaintegral(n,n))
    !Note on normalisation:
    !With cgyro tauab_inv:
    ! with lambda=va/vb =sqrt(Ta/Tb * mb/ma)
    ! Changerate_a=4*pi^(-1.5) tauab_inv emat_a^-1 *lambda^2 ( Test_ab(1,lambda)+ ma/mb*field_ab(1,lambda))
    ! We have  na*<a|changerate_a_due to field> = nb*<b|changerate_b_due to field> hopefully.
    ! and tauab_inv=tauba_inv*nb/na*mb^2/ma^2*vb^3/va^3=tauba_inv*nb/na va/vb Tb^2/Ta^2
    ! changerate_b_dtf=tauab_inv*na/nb vb/va Ta^2/Tb^2 mb/ma*1/lambda^2 field_ab(1,1/lambda)
    !                 = ....          ".....................*1/lambda (field_ab(1,lambda)^T
    !                 =tauab_inv*na/nb va^2/vb^2 ma/mb * (field_av(1,lambda))^T
    !                 =tauab_inv*na/nb lambda^2 ma/mb * (field_av(1,lambda))^T ok.

    real x,p0,p1,pi,q0,q1,qi
    integer i,j,k

    deltaintegral=0
    do i=1,ngauss
       x=gp(i)*xmax
       do j=1,n
          if (j==1) then
             p1=0
             pi=c1(1)*gw(i)*exp(-x**2)*x**2*xmax
          else
             p0=pi
             pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
             p1=p0
          endif
          q1=0
          qi=c1(1)*pi
          k=1
          deltaintegral(k,j)=deltaintegral(k,j)+qi
          do k=2,n
             q0=qi
             qi=(c1(k)*x-a1(k))*qi-b1(k)*q1
             q1=q0
             deltaintegral(k,j)=deltaintegral(k,j)+qi
          enddo
       enddo
    enddo
  end subroutine genenergymatrix
  subroutine genfieldkernel(n,lmax,a1,b1,c1,xmax,beta,t1t2ratio,gp,gw,ngauss&
       ,gp2,gw2,ng2,intkernel,intkernel2)
    ! driver routine.
    ! generates field particle kernel in intkernel dependent on beta and
    ! t1t2ratio.
    implicit none
    optional :: intkernel2
    integer, intent(in) :: n,lmax,ngauss,ng2  ! lmax=max l-**index**, lphys=l-1
    real, intent(in) :: a1(n),b1(n),c1(n),& !poly rec-coeffis
         gp(ngauss),gw(ngauss),& ! gauss points and weights ngauss> n
         gp2(ng2),gw2(ng2),& ! gauss points and weights for small kernel ng2~8
         xmax,beta,t1t2ratio !integration range and relative factor in exponenti
    real, dimension(n,n,lmax), intent(out) :: intkernel(n,n,lmax)&
         ,intkernel2(n,n,lmax)
    real, allocatable,dimension(:,:,:) :: intkernel1,intlokernel,intlokernel1
    real, allocatable,dimension(:,:) :: deltaintegral
    integer i,j,l
    real  t1,t2

    call cpu_time(t1)

    if (t1t2ratio/=1 .and. .not. present(intkernel2)) then
       print '(A,G23.16,A)','genfieldkernel: t1t2ratio=',t1t2ratio,'and optional arg intker&
            &nel2 missing.'
    end if

    if (t1t2ratio==1) then
       call genintkernel(n,lmax,a1,b1,c1,xmax,beta,1.,gp,gw,ngauss,gp2&
            ,gw2,ng2,.true.,intkernel)
       if (beta==1) then
          do i=1,n
             do j=1,i-1
                intkernel(i,j,:)=intkernel(i,j,:)+intkernel(j,i,:)
                intkernel(j,i,:)=intkernel(i,j,:)
             enddo
             intkernel(i,i,:)=intkernel(i,i,:)*2
          enddo
       else
          allocate(intkernel1(n,n,lmax))
          call genintkernel(n,lmax,a1,b1,c1,xmax,1./beta,1.,gp,gw&
               ,ngauss,gp2 ,gw2,ng2,.true.,intkernel1)
          do i=1,n
             do j=1,n
                intkernel(i,j,:)=intkernel(i,j,:)+intkernel1(j,i,:)/beta
             enddo
          enddo
          deallocate(intkernel1)
       end if
    else
       allocate(intlokernel(n,n,lmax),intlokernel1(n,n,lmax))

       call genintkernel(n,lmax,a1,b1,c1,xmax,beta,t1t2ratio,gp,gw,ngauss,gp2&
            ,gw2,ng2,.true.,intkernel,intlokernel)

       call genintkernel(n,lmax,a1,b1,c1,xmax,1./beta,1/t1t2ratio,gp,gw&
            ,ngauss,gp2 ,gw2,ng2,.true.,intkernel2,intlokernel1)
       do i=1,n
          do j=1,n
             intkernel(i,j,:)=intkernel(i,j,:)+intlokernel1(j,i,:)/beta
             intkernel2(i,j,:)=intkernel2(i,j,:)+intlokernel(j,i,:)*beta
          enddo
       enddo
       deallocate(intlokernel,intlokernel1)
    end if
    allocate(deltaintegral(n,n))
    call deltaintegrate(n,a1,b1,c1,xmax,beta,gp,gw,ngauss,deltaintegral)
    do l=1,lmax
       intkernel(:,:,l)=deltaintegral-intkernel(:,:,l)
    end do
    if (t1t2ratio/=1) then
       do i=1,n
          do j=1,n
             intkernel2(i,j,:)=deltaintegral(j,i)*beta-intkernel2(i,j,:)
          end do
       end do
    else
       if (present(intkernel2)) then
          do i=1,n
             do j=1,n
                intkernel2(i,j,:)=intkernel(j,i,:)*beta
             end do
          end do
       endif
    endif
!!$2   format (*(G0,"  "))
!!$      block
!!$        real di1(n,n),val
!!$        call deltaintegrate(n,a1,b1,c1,xmax,1/beta,gp,gw,ngauss,di1)
!!$        print 2,'Check deltaintegral self adjointness:'
!!$        print 2,deltaintegral(2,1),di1(1,2)/beta
!!$        val=0
!!$        do i=1,n
!!$           do j=1,n
!!$              val=max(val,abs(deltaintegral(i,j)-di1(j,i)/beta))
!!$           end do
!!$        end do
!!$        print 2,'max nonself',val
!!$      end block
    deallocate(deltaintegral)
    call cpu_time(t2)
    if (verbose>1) print '(A,G13.3)','genfieldkernel took',t2-t1
  end subroutine genfieldkernel
end module landau
