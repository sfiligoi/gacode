!gfortran -march=native -g -fcheck=all -fno-stack-arrays -fimplicit-none -fdefault-real-8 gyrotransformation.f90 -c
!gfortran -march=native -O3 -fno-stack-arrays -fimplicit-none -fdefault-real-8 gyrotransformation.f90 -c
!gfortran -march=native -O3 -fno-stack-arrays -fimplicit-none -fdefault-real-8 gyrotransformation.f90 landau.f90 half_hermite.f90 gauss_legendre.f90 -c
!intel:
!ifort -stand f15 -warn all -march=native -O3 -heap-arrays 10 -implicitnone -real-size 64 gyrotransformation.f90 -c
!ifort -stand f15 -warn all -xCORE-AVX512 -O3 -heap-arrays 10 -implicitnone -real-size 64 gyrotransformation.f90 -c

module gyrotransformation
  real, parameter,private :: epspullback=1e-13,epsgyrocolmat=1e-15
  real, parameter,private :: pi1=atan(1.)*4
  integer :: verbose=0
contains
  integer function est_mpullback(kbounda,kboundb,eps)
    implicit none
    !real :: kbounda=16.
    real,intent(in) :: kbounda,kboundb, &   !usually krho*xmax with maximum
         eps                                !krho planned.
    optional :: kboundb,eps
    real val,val1,kba,kbb,epsp
    integer mmax
    ! The maximum mpullback (with two species) relies on m being a conserved
    ! quantum number in the collisions. Therefore we can use the product of
    ! the smallnesses of the Bessel functions as indicator of the global
    ! smallness. (For large m, the Bessel functions are asymptotically equal
    ! to Jn[x]~(2/m)^(1/3) Ai[(2/m)^(1/3)(n-x)] increase slowly up to their
    ! 1st maximum, which is around n+1.01879*(n/2)^(1/3) and on the order of
    ! (2/m)^(1/3)*0.535656656015.)  So mpullback depends on kbounda and
    ! kboundb
    
    if (verbose>0) print '(A)','Warning, need abs(kbound) est_mpullback!'
    if (present(kboundb)) then
       kbb=kboundb
    else
       kbb=kbounda
    end if
    !kbb=kbounda
    if (present(eps)) then
       epsp=eps
    else
       epsp=epspullback
    endif
    if (verbose>0) print '(A)','WARNING: using max of kb for both species in est_mpullback please improve me.'
    kba=kbounda
    kba=max(kbounda,kbb)
    kbb=kba
    val=0
    mmax=max(200,ceiling(2*max(kbounda,kbb)))
    do est_mpullback=mmax,1,-1
       val1=val
       val=val+abs(bessel_jn(est_mpullback,kba)*bessel_jn(est_mpullback,kbb))
       if (val>epsp) exit
    end do
!    est_mpullback=est_mpullback*2
!!$    if (verbose>0) print "(3(A,G0),4(A,ES11.4),2ES11.4)",'mpullback based on Jn is ',est_mpullback,' initial guess was ',ceiling(2 *max(kbounda&
!!$         ,kbb)),'k1,k2',kbounda,kbb,' Bessel_j(' ,est_mpullback,',',kbounda,')=',bessel_jn(est_mpullback ,kbounda),' squaresum=',val,' sqrt=' &
!!$         ,sqrt(val) ,val1,sqrt(val1)
    if (verbose>0) print '(2(A,I0),A,2G23.16,A,I0,4(A,G23.16),5(G24.16))','mpullback based on Jn is ',&
         est_mpullback,' initial guess was ',ceiling(2*max(kbounda,kbb)),' k1,k2',kbounda,kbb,&
         ' Bessel_j(' ,est_mpullback,',',kbounda,')=',&
         bessel_jn(est_mpullback ,kbounda),' squaresum=',val,' sqrt=',sqrt(val) ,val1,sqrt(val1)
    if (est_mpullback==max(200,ceiling(2*max(kbounda,kbb)))) then
       print '(A)','mpullback wrong!'
       stop
    end if
  end function est_mpullback

  integer function est_extradegree(kbound,eps,mmode)
    ! kbound is necessary only for the species for which the extradegree is to
    ! be calculated. The extradegree may be different for two species.
    ! OBSERVATION: est_extradegree near identical est_k_sampling
    implicit none
    !real :: kbound=16.
    real,intent(in) :: kbound, &   !usually krho*xmax with maximum krho planned.
         eps
    integer,intent(in),optional :: mmode ! azimuthal mode number m used for
                                         ! estimation, normally should be 0, the worst
                                         ! case.
    optional :: eps
    real epsp
    real val,x
    integer i,j,k,mm
    real, allocatable :: sample(:)

    if (verbose>0) print '(A)','Warning, need abs(kbound) est_extradegree!'
    mm=0
    if (present(mmode)) mm=mmode
    if (present(eps)) then
       epsp=eps
    else
       epsp=epspullback
    endif
    k=ceiling(kbound)+20
    allocate(sample(2*k))
    do i=1,2*k
       x=(2*i-1)*(.25*pi1/k)
       !Now if we want to deal with polynomials up to degree k,
       !we need at least k vertices. With 2k vertices we can deal with
       !polynomials up to degree 2k.
       !x should from pi/4k up to pi-pi/4k

       sample(i)=bessel_jn(mm,(cos(x)+1)*(.5*kbound))-bessel_jn(mm,0.)
       ! subtracting bessel(0) here cures numerical precision issues for krho~0
       !! sample(i)=gsl_sf_bessel_jl(mm,(sin(x)+1)*.5*kbound)-gsl_sf_bessel_jl(mm,0.)
    end do
    do j=k,1,-1
       val=0
       do i=1,2*k
          x=(2*i-1)*(.25*pi1/k)
          val=val+sample(i)*cos(j*x)
       end do
       val=val/k
       if (verbose>4) print '(A,I5,G24.16)','deg',j,val
       if (abs(val)>epsp) exit
    end do
!!$    if (j==0) then
!!$       print 2,'Did not find proper degree to sample Bessel fct!'
!!$       stop
!!$    end if
    if (verbose>0) print '(A,I0,A,G23.16)',&
         'Extra degree needed for sampling of Bessel fct=',j,' at k=',kbound
    est_extradegree=j
  end function est_extradegree

  integer function est_k_sampling(kb1,eps)
    ! Estimate number of Chebyshev points needed to represent exp(i k)
    ! in the range k=0...kbound with at most eps error.
    ! This is necessary because we do not want to calculate the gyro-
    ! transformation for *every* individual kperp occuring in the code,
    ! which would cause much double work due to many near-identities.
    
    implicit none
    !real :: kb1=16.
    real,intent(in) :: kb1, &   !usually max(kperp)*xmax*2 for worst species
         eps
    optional :: eps
    real epsp
    real x,kbound
    complex val
    integer i,j,k

    kbound=abs(kb1)
    if (present(eps)) then
       epsp=eps
    else
       epsp=epspullback
    endif
    k=ceiling(kbound)+20
!!$    chebysample:  block
!!$      complex*16 sample(2*k),val
!!$      do i=1,2*k
!!$         x=(i-.5)*(.5_16*pi1/k)
!!$         ! with n=2k
!!$         ! sampling exp(i kbound*y) at  y=.5*(1+cos(pi/n*(2*i-1)))
!!$         sample(i)=exp(sin(.5*x)**2*((0,1)*kbound))-1
!!$      end do
    do j=k,1,-1
!!$         val=0
!!$         do i=1,2*k
!!$            x=(i-.5)*(.5*pi1/k)
!!$            val=val+sample(i)*cos(j*x)
!!$         end do
!!$         val=val/k
!!$         if (verbose>4) print 2,'deg',j,val,'should',&
!!$              exp((0,.5)*kbound)*(0,1.)**j*2*bessel_jn(j,kbound*.5)
       val=exp((0,.5)*kbound)*(0,1.)**j*2*bessel_jn(j,kbound*.5)
       if (verbose>4) print '(A,I5,": (",G23.16,",",G23.16,")")','deg',j,val
       if (abs(val)>epsp) exit
    end do
!!$    end block chebysample
    j=j+1 ! j is the min. *frequency* - we need one more interpolation point
    if (verbose>0) print '(A,I0,A,G23.16,A)','Found k-sample-number needed for gyro phases=',&
         j,' at k=',kbound,' WARNING should depend on both species'
    est_k_sampling=j
  end function est_k_sampling
  
  subroutine calc_projleg(projleg,lmax2,lmax0,gpl,gwl)
    ! calculate half of the projection matrix for the Legendrepolynomials.
    ! only half is needed, because the other half can be infered from parity.
    ! the complete set of vertices is therefore lmax2*2, which is even.
    ! actually--> this is part of calc_projassleg. Should use that one.
    implicit none
    integer,intent(in) :: lmax2,lmax0
    real,intent(out) :: projleg(lmax2,lmax0)
    real,intent(in) :: gwl(lmax2),gpl(lmax2) ! gauss weights and points
    ! ^^may be defined also for more than lmax2, no problemo
    real p1,p0,pi,x
    integer i,k,j,l
    
    ! projleg might be improvable by using an eigenvector routine or real*16 variables
    k=0
    do i=1,lmax2
       x=gpl(i)
       p1=0
       pi=sqrt(gwl(i))
       j=1
       projleg(i,j)=pi*sqrt(j-1+.5)
       do j=2,lmax0
          p0=pi
          pi=((2*(j-2)+1)*x*pi-(j-2)*p1)/(j-1)
          p1=p0
          projleg(i,j)=pi*sqrt(j-1+.5)
       enddo
       if (verbose>4) then
          do l=1,lmax0
             if (projleg(i,l)==0) then
                print *,'projleg zero??',i,l,projleg(i,l)
                k=k+1
             else if(abs(projleg(i,l))<epsgyrocolmat) then
                print *,'projleg near zero',i,l,projleg(i,l)
                k=k+1
             end if
          end do
       end if
    enddo
    if (verbose>4) then
       print '(A,I5,A)','projleg had',k,' strange coefficients.'
       do l=1,lmax0
          print 2,'projleg l norm=1',l,sum(projleg(:,l)**2)
       enddo
       do l=1,lmax0-1
          print 2,'projleg l norm=0',l,sum(projleg(:,l)*projleg(:,l+1))
       enddo
       do l=1,lmax2
          print 2,'projleg l norm=1',l,sum(projleg(l,:)**2)
       enddo
       do l=1,lmax2-1
          print 2,'projleg l norm=0',l,sum(projleg(l,:)*projleg(l+1,:))
       enddo
    end if
2   format (A,I5,G25.16)
    
  end subroutine calc_projleg


  subroutine gyrotrafo(gyrocolmat,nmax0,lmax0,colmat,ncolmat,lcolmat,projsteen,sp,nsteen,xmax,krhoa,krhob,eps)
    implicit none
    integer,intent(in) :: lmax0,nmax0,ncolmat,lcolmat,nsteen
    real,intent(in) :: xmax,krhoa,krhob,eps !for a and b species
    optional :: krhob,eps !default is equal to krhoa

    ! In the first iteration we use equal nmax0 and nmaxpoly for all species.
    ! For unequal n we also need to update landau.f90
    ! In that case we need two projsteen and sp, and beskl (see below)
    ! We may however compute the actually needed nmaxpoly and possibly nmax0.

    ! We take the angular functions square normalised to unity (not as the Ylm).
    ! As does energymatrix() and the other landau routines.

    real,intent(out) :: gyrocolmat(nmax0,lmax0,nmax0,lmax0)
    ! output calculated gyrokinetic collision matrix.
    ! a-species: second pair of indices, b: first pair

    real,intent(in) :: colmat(ncolmat,ncolmat,lcolmat) !input Landau collision matrix
    real,intent(in) :: projsteen(nsteen,nsteen),sp(nsteen)
    !projection matrix (i.e. sqrt(weight(i))*P(sp(i)) and vertices.
    real krb
    real,dimension(:),allocatable :: gpl,gwl ! Gauss points
    real,dimension(:,:),allocatable :: projleg,v_theta_space,v_l_space,steen0_theta_space,tmp_steen_space
    real,dimension(:,:,:),allocatable :: beskla,besklb,steen_l_space
    real,dimension(:,:,:,:),allocatable :: projassleg 
    real val
    real epsp
    integer mpullback,lmax,extradegree,nmaxpoly
    integer i,j,k,l,m,m2max,l2
    integer oe,moe !odd even angular momentum =1: even =2:odd, moe=0/1: even/odd m numbers
    integer mphys,lphys  ! real mode numbers, otherwise indices are used, which are those +1
    integer m1,m2,o,q !used for verbose output in projassleg
    integer, external ::idamax

    ! originally in block "timeblock":
    real tbes,t(11),tto,t1,t2,cost(11)
    real,allocatable,dimension(:,:,:) :: max_n_comp,max_l_comp

    ! originally in block "lminblock":
    integer lminb
    real lbval,lminbval

    ! originally in certain diagnostic blocks:
    integer idx,a(4)
    real v1
    integer a1(4),a2(4)
    real div,siv,sym,sym1

    real,allocatable,dimension(:) :: maxn,maxl

    if (present(krhob)) then
       krb=krhob
    else
       krb=krhoa
    end if

    if (present(eps)) then
       epsp=eps
    else
       epsp=epspullback
    endif

    ! First a few sanity checks and estimates


    if (verbose>0 .and. mod(lcolmat,2)==1) then
2      format (A,9(I0,A))
       print 2,'Warning in gyrotrafo: lmax is assumed to be even, but input collision matrix has l=',lcolmat,'.'
       print 2,'Cannot use this extra information. Because lmax is needed to be even.'
    end if

    mpullback=est_mpullback(krhoa*xmax,kboundb=krb*xmax,eps=eps)
    if (verbose>0) print 2,'mpullback must be >=',mpullback
    lmax=lmax0+mpullback
    ! make lmax an even number
    if (mod(lmax,2)==1) then
       mpullback=mpullback+1
       lmax=lmax0+mpullback
       if (verbose>0) print 2,'Corrected mpullback to make lmax even mp.=',mpullback
    end if
    if (verbose>0) print 2,'Using lmax=',lmax,'. (lmax0=',lmax0,')'


    if (lcolmat<lmax) then
       print 2,'Error in gyrotrafo: lcolmat=',lcolmat,'<lmax. Cannot continue.'
       stop
    end if

    extradegree=est_extradegree(max(abs(krhoa),abs(krb))*xmax,eps=eps)

    nmaxpoly=nmax0+extradegree
    if (verbose>0) print 2,'Need extradegree',extradegree,'-> max. polynomial degree=',nmaxpoly-1
    if (nmaxpoly>ncolmat) then
       print *,'Error in gyrotrafo: Max. polynomial degree in colmat='&
            ,ncolmat-1,', but need',nmaxpoly-1
       stop
    end if
    if (nmaxpoly>nsteen) then
       print *,'Error in gyrotrafo: Max. polynomial degree in projsteen/sp='&
            ,nsteen-1,', but need',nmaxpoly-1
       stop
    end if

    ! So now for the pullback:
    ! we need first the pullback projection matrices.
    ! we define a maximum L, that is considered. Everything above Lmax is
    ! projected out in real particle space.

    ! First we set up separate Steen and Legendre Projectionmatrices: (We
    ! normalise the Legendre polynomials by the factor sqrt(n+.5) so that the
    ! squareintegral is 1).

    !projleg might be improvable by using an
    !eigenvector routine or real*16 variables

    if (verbose>4) then
       print *,'Steen function projsteen(nmax0,nmaxpoly)=',projsteen(nmax0&
            ,nmaxpoly),'at sp=',sp(nmaxpoly),'nmax0=',nmax0
       print *,'Steen function projsteen(nmaxpoly,nmaxpoly)=' &
            ,projsteen(nmaxpoly,nmaxpoly),'at gp=',sp(nmaxpoly),'nmaxpoly='&
            ,nmaxpoly
!!$  do i=1,nmaxpoly
!!$     print *,'projsteen norm=1',i,sum(projsteen(i,:)**2)
!!$  enddo
!!$  do i=1,nmaxpoly-1
!!$     print *,'projsteen norm=0',i,sum(projsteen(i,:)*projsteen(i+1,:))
!!$  enddo

    end if

    allocate(gwl(lmax),gpl(lmax))
    call gauss_legendre(-1.,1.,gpl,gwl,lmax)
    ! If lmax is even, we need only half of the vertices.
    allocate(projleg(lmax/2,lmax0)) ! Convert from lmax0-L-space to lmax-xi-space
    call calc_projleg(projleg,lmax/2,lmax0,gpl,gwl)

    !Now calculate projector for associated Legendre functions at Gauss points gpl.
    !These functions are only polynomials for even m.
    !Parity: P^m_l is an even function, if m+l even.
    !Therefore often only half sums have to be calculated.
    ! And we assume lmax to be even.

    allocate(projassleg(lmax/2,lmax/2,2,mpullback+1))
    call calc_projassleg(projassleg,lmax/2,lmax/2,mpullback+1,gpl,gwl)
    !Note: projleg(i,j)=projassleg(i,(j+1)/2,mod(j+1,2),1)
    !here: m=mphys+1, j=jphys+1: (later m=mphys; mpullback=phys.)
    !oe=mod(j+m,2)+1
    !projassleg-old(i,j,m)=projassleg(i,(j-1)/2+1,mod(j+m,2)+1,m)
    !                     =projassleg(i,(j+1)/2,mod(j+m,2)+1,m)
    !projassleg(i,j2,oe,m)=projassleg-old(i,j2*2-mod(j,2),m)
    !                     =projassleg-old(i,j2*2-mod(oe+m-1,2),m)
    !                     =projassleg-old(i,j2*2+mod(oe+m,2)-1,m)
    !projassleg-old is even, if j+m is even, or if oe is 1 (and not 2).

    !old:
    !projassleg(i,j,m)=Pmj(gpl(i))* sqrt(   (2lphys+1) (lphys-mphys)! /(2(lphys+mphys)!) )
    !  
    ! nor(j,m)=sqrt(   (2j-1) (j-m)! /(2(j+m-2)!) )
    ! nor(m,m)=sqrt( (2m-1)/(2 (2m-2)!)
    ! nor(m,m)/nor(m-1,m-1)= sqrt((2m-1)/(2m-3)/(2m-2)/(2m-3))=sqrt((2m-1)/(2m-2)) / (2m-3)
    !            =sqrt((2m-1)(2m-4)/(2m-3)/(2m-2))
    ! nor(j,m)/nor(j-1,m)=sqrt( (2j-1)/(2j-3)*(j-m)/(j+m-2) )

    ! projassleg might be improvable by using an eigenvector routine or real*16 variables
!!$    m1=0;m2=0;q=0;o=20 !only for verbose
!!$    projassleg=0
!!$    j=0
!!$    do i=1,lmax/2
!!$       x=gpl(i)
!!$       do m=1,mpullback+1
!!$          ! we start at pmm and move to higher L
!!$          p1=0
!!$          mphys=m-1
!!$          if (m==1) then
!!$             !P00=1*sqrt(1/2)
!!$             pi=sqrt(.5*gwl(i))
!!$             projassleg(i,1,1,1)=pi
!!$          else
!!$             pi=-sqrt(1-x**2)*projassleg(i,m/2,1,m-1)*sqrt((2*m-1.)/(2*m-2.))
!!$             projassleg(i,(m+1)/2,1,m)=pi
!!$          endif
!!$          do l=m+1,lmax
!!$             lphys=l-2
!!$             p0=pi
!!$             pi=( (2*lphys+1)*x*pi-(lphys+mphys)*p1 )/(lphys-mphys+1)
!!$             pi=pi*sqrt( (2*l-1)/(2*l-3.)* ((l-m)/(l+m-2.)) )
!!$             p1=p0*sqrt( (2*l-1)/(2*l-3.)* ((l-m)/(l+m-2.)) )
!!$             projassleg(i,(l+1)/2,mod(l+m,2)+1,m)=pi
!!$          enddo
!!$          if (verbose>0) then
!!$             do l2=1,lmax/2
!!$                do oe=1,2
!!$                   if (l2*2+mod(oe+m,2)-1<m) cycle
!!$                   j=j+1
!!$                   if (projassleg(i,l2,oe,m)==0) then
!!$                      print 2,'projassleg zero??',i,x,l2,oe,m,projassleg(i,l2,oe,m)
!!$                      m1=m1+1
!!$                   else if(abs(projassleg(i,l2,oe,m))<epsgyrocolmat) then
!!$                      !print 2,'projassleg near zero',i,x,l2,oe,m,projassleg(i,l2,oe,m)
!!$                      m2=m2+1
!!$                   end if
!!$                end do
!!$             end do
!!$          end if
!!$       enddo
!!$    enddo
!!$    print 2,'projassleg had',m1,'0 and',m2,'near 0 coefficients of',j,'total',(mpullback+1)*lmax/2*lmax/2*2
    ! Now should better check the normalisation
!!$  m=(lmax+1)/2
!!$  m=mpullback+1
!!$  do l=m,lmax
!!$     print 2,'m,l norm=1',m,l,sum(projassleg(:,l,m)**2)
!!$  enddo
!!$  do l=m,lmax-1
!!$     print 2,'m,l norm=0',m,l,sum(projassleg(:,l,m)*projassleg(:,l+1,m))
!!$  enddo
    ! expansion of exp(i k*v_perp) = exp(i k_perp*v) in spherical harmonics:
    !    sum l=0 ... (2l+1)i^l Jspher_l(kvperp) P_l(kvperp*sin(theta))
    ! On the other hand Y0l()=sqrt((2l+1)/4pi) P_l(cos theta)
    !    ==> sum l=0.. sqrt((2l+1)*4pi) Y0l
    ! However we are not normalising our P_l and Pml correctly but so that
    ! <P_l^2>=2pi
    !  so  exp(i k*v_perp) = sum l=0 ... (2l+1)i^l Jspher_l(kvperp) sqrt(2/(2l+1))myP_l(kvperp*sin(theta))
    !  so  exp(i k*v_perp) = sum l=0 ... i^l Jspher_l(kvperp) sqrt(2*(2l+1))myP_l(kvperp*sin(theta))
    ! Now loop over steen index and L index

    allocate(beskla(nsteen,lmax/2,0:mpullback))
    if (krb/=krhoa) allocate(besklb(nsteen,lmax/2,0:mpullback))
    allocate(steen0_theta_space(nmax0,lmax/2))
    allocate(v_theta_space(nsteen,lmax/2))
    allocate(v_l_space(nsteen,lmax/2))
    allocate(steen_l_space(nmaxpoly,0:mpullback/2,lmax/2),tmp_steen_space(nmaxpoly,0:mpullback/2))
    ! Wir zerlegen alle Funktionen in die m,L und setzen sie wieder zusammen.
    ! Dann schauen wir ob was am RMS fehlt.
    !    print 2,'sumprojsteen',sum(projsteen),sum(projassleg),sum(projleg)
    gyrocolmat=0
    val=0
    if (verbose>4) then
       allocate(max_n_comp(nmaxpoly,nmax0,lmax0),max_l_comp(lmax,nmax0,lmax0))
       max_n_comp=0
       max_l_comp=0
       q=0;o=20 ! for error messages
    end if
    tbes=0
    t=0
    cost=0
    call cpu_time(tto)
    call cpu_time(t1)
    do m=mpullback,0,-1
       do l2=1,lmax/2
          beskla(:,l2,m)=bessel_jn(m,krhoa*sp*sqrt(1-gpl(l2)**2))
          if (krb/=krhoa) then
             besklb(:,l2,m)=bessel_jn(m,krb*sp*sqrt(1-gpl(l2)**2))
          end if
          ! don't like sqrt(1-gpl(l)**2) for sin theta. Should use something more accurate. ???
       enddo
       !         print 2,'m=',m,'sumbeskl',sum(beskla(:,1:lmax/2,m))
    end do
    call cpu_time(t2)
    tbes=tbes+(t2-t1)


    jloop: do j=1,lmax0
       oe=mod(j+1,2)+1 !oe=1: even in z oe=2: odd in z
       lminb=0
       lminbval=0
       !m=mphys l=lphys+1 here
       !l>=m+1, l<=lmax and oe+m+l is even.
       !==> l>=m+oe. since lmax is even l<=lmax-mod(oe+m,2)
       !the correct arguments are projassleg(,(l+1)/2,oe,m+1)
       !==> (l+1)/2 >=(m+oe+1)/2  (l+1)/2<=(lmax-mod(oe+m,2)+1)/2=lmax/2 +(1-mod(oe+m,2))/2=lmax/2
!!$          do l2=(m+oe+1)/2,lmax/2 !n_elements=lmax/2-(m+oe+1)/2+1=(lmax+2-m-oe)/2
!!$             lbval=sum(beskla(nmaxpoly,:)*projleg(:lmax/2,j)&
!!$                  *projassleg(:,l2,oe,m+1))
!!$             !if (m==1 .and. j==30) print 2,'m,j,l,lbval',m,j,l,lbval
!!$             if (abs(lbval)>epsp) then
!!$                lminb=l2*2+mod(oe+m,2)-1  +1
!!$                lminbval=lbval
!!$             end if
!!$          end do
!!$          print 2,'m',m,'j',j,'lminb',lminb,'lminbval',lminbval
       !^^ can use this boundary value as estimate of max l needed.
       ! most conservative value for m=0.
       ! can also reduce lmax by that. But will not be greatly different from lmax.

       iloop: do i=1,nmax0
          moeloop: do moe=0,1 ! do odd and even m separately
             ! we do now all odd or even m in one bunch to be able to reuse the read in collision matrix.
             ! costs a little more memory.
             mloop:  do m=mpullback-mod(moe+mpullback,2),moe,-2 ! Leaving out negative m's here. That must be accounted for when doing sums.
                !number of loop traversals=(mpullback-mod(moe+mpullback,2)-moe+2)/2
                !                         =(mpullback-moe+2-mod(mpullback-moe+2,2))/2
                !                         =(mpullback-moe+2)/2  (carefull here: (-1)/2=0)
                ! max. traversals occur for moe=0 and are (mpullback+2)/2=mpullback/2+1 since mpullback>=0.
                call cpu_time(t1)
                do l2=1,lmax/2
                   do k=1,nsteen
                      v_theta_space(k,l2)=beskla(k,l2,m)*projsteen(i,k)*projleg(l2,j)
                      ! should be same as *projassleg(l2,j/2,oe,1)
                      ! Note: sum_i=-infty,+infty  bessel_jn(i,x)^2 = 1 for any x.
                   enddo
                enddo
                call cpu_time(t2)
                t(1)=t(1)+(t2-t1)
                cost(1)=cost(1)+nsteen*lmax/2
                !cost=lmax*n maybe lmax/2*n
                ! Now transform into l space. For this we have to multiply with the Pml for l>=m.
                ! for even/odd jphys: l+m must be even/odd.
                ! Smallest possible lphys:
                !  for even jphys (oe=1): lphys=mphys
                !  for odd  jphys (oe=2): lphys=mphys+1
                !     (sometimes cannot be mphys, because that has wrong parity.)
                ! ==> lphysmin=mphys+oe-1 always.
                ! l2min=(lphysmin)/2+1=(mphys+oe-1)/2+1=(mphys+oe+1)/2
                ! Since lmax is even (!):
                ! number of l2=lmax/2+1-(m+oe+1)/2=(lmax+2)/2-(m+oe+1)/2=(lmax+2-m+oe)/2
                ! This is independent for every k.
                ! Matrix multiply with upper triangular matrix projassleg:
                !             v_l_space(:,1:m)=0
                if (verbose>4) v_l_space=1e300
                call cpu_time(t1)
                call dgemm('n','n',nsteen,(lmax-m+2-oe)/2,lmax/2,1.,v_theta_space,nsteen,&
                     projassleg(:,(m+oe+1)/2:lmax/2,oe,m+1),lmax/2,0.,&
                     v_l_space(:,(m+oe+1)/2:lmax/2),nsteen)
                call cpu_time(t2)
                t(2)=t(2)+(t2-t1)
                cost(2)=cost(2)+nsteen*((lmax-m+2-oe)/2)*lmax/2
                ! Only l remain here, for which m+j+l=even.
                ! or j+oe=even.
                ! Of course also on the way back parity is conserved.
                ! We have therefore separate even l and odd l gyrocolmat.
                ! This makes total sense.
                ! Pml and projassleg is an even function for 
                ! Now we work on v_l_space with the steen polynomials and the collision operator.
                ! we multiply over the left index of v_l_space
                if (verbose>4) steen_l_space(:,m/2,:)=1e300
                call cpu_time(t1)
                call dgemm('n','n',nmaxpoly,(lmax-m+2-oe)/2,nsteen,1.,projsteen,nsteen&
                     ,v_l_space(:,(m+oe+1)/2:lmax/2),nsteen,0.,&
                     steen_l_space(:,m/2,(m+oe+1)/2),nmaxpoly*(1+mpullback/2))
                call cpu_time(t2)
                t(3)=t(3)+(t2-t1)
                cost(3)=cost(3)+nmaxpoly*((lmax-m+2-oe)/2)*nsteen
                !could conceivably cut off this matrix product
                if (verbose>4) then
                   do k=1,nmaxpoly
                      max_n_comp(k,i,j)=max(max_n_comp(k,i,j),maxval(steen_l_space(k,m/2,(m+oe+1)/2:lmax/2)))
                   end do
                   do k=(m+oe+1)/2,lmax/2
                      l=k*2-mod(oe+m,2)
                      max_l_comp(l,i,j)=max(max_l_comp(l,i,j),maxval(steen_l_space(:,m/2,k)))
                   end do
                end if
                ! epspullback-small whenever l>j+mpullback or l<j-mpullback (bec of spherical Bessels)
                ! on the otherhand: l>m
                ! for backtransform small, if
                ! j<l-mpullback or j> l+mpullback
                ! could be checked how far out needed, also |v| cut-off
                ! if there is no redistribution of l,
                ! it is sufficient to have sqrt(epspullback)??
                ! only lmax cutoff, not mpullback
                ! 
                ! Now we do whatever we want with steen_l_space.

                !  here
                ! can use dgemm_batch of mkl library.
                if (verbose>4) then
                   steen_l_space(:,m/2,:(m+oe+1)/2-1)=1e300
                end if
             end do mloop
             if (verbose>4 .and. (mpullback-mod(moe+mpullback,2))/2<mpullback/2) then
                steen_l_space(:,mpullback/2,:)=1e300
             end if
             call cpu_time(t1)
             ! This is really a costly procedure. dgemm cannot reuse colmat matrix elements for different l,l2.
             ! The dominant cost is reading through the colmat matrices.
             ! This is on an ivy-bridge 10x slower than the other operations per flop.
             ! Maybe one should do several i,j at once.
             ! The smallest l2 for given m is lphysmin/2+1=(m+oe-1)/2+1=(m+oe+1)/2.
             ! The smallest l2 for *all* m in this loop is (moe+oe+1)/2, since the smallest m is moe.
             ! For that l2 we have at least moe as a possible m.
             do l2=(moe+oe+1)/2,lmax/2
                l=l2*2-mod(oe+moe,2) !--> comments in calc_projassleg
                ! l=lphys+1
                ! this should be equal (l2-(moe+oe+1)/2)*2+lmin,
                !            lmin(moe,oe)=moe+1+mod(oe+1,2)=moe+oe
                ! i.e. l=(l2-(moe+oe+1)/2)*2+moe+oe
                !        =l2*2-(moe+oe+1)/2*2+moe+oe
                !        =l2*2+moe+oe+1-((moe+oe+1)/2*2-1
                !        =l2*2+mod(moe+oe+1,2)-1
                !       =l2*2-mod(moe+oe,2)   ok!

                ! for small l2 we do not have to do all m. The condition on m for given l2 is
                ! l2>=(m+oe+1)/2 <=> 2*l2>=(m+oe+1)/2*2 <=> 2*l2+1>=m+oe+1 <=> m<=2*l2-oe>=moe
                ! check: oe=1 ==> (moe+oe+1)/2*2-oe=moe/2*2+1=1
                !        oe=2 ==> (moe+oe+2)/2*2-oe=(moe+1)/2*2=moe ok.
                ! In addition mod(m+moe,2)==0 ==> m<=2*l2-oe-mod(2*l2-oe+moe,2)=
                !          =2*l2-oe-mod(oe+moe,2)
                ! # ==> m/2<=(2*l2-oe-mod(oe+moe,2))/2=l2-oe+(oe-mod(oe+moe,2))/2
                !    =l2-oe+(moe+oe-moe-mod(oe-moe,2))/2
                !    =l2-oe+(oe-moe)/2
                ! check: case distinction for #:
                ! oe=1,moe=0 ==> (2*l2-1-mod(1,2))/2=(2*l2-2)/2=l2-1
                ! oe=1,moe=1 ==> (2*l2-1)/2=l2-1 (l2>=1)
                ! oe=2,moe=0 ==> (2*l2-2)/2=l2-1
                ! oe=2,moe=1 ==> (2*l2-2-1)/2=(2*l2-3)/2 and l2>=2 in that case
                !            ==> (2*l2-4)/2=l2-2
                ! in total l2-1 excebt if oe+moe=3 then l2-2. ==> l2-(oe+moe+1)/2
                ! ok clear. For every m+=2 we gain another l2 and vice versa.
                ! in addition m<=mpullback/2.
                m2max=min(l2-(oe+moe+1)/2,mpullback/2)
                if (verbose>4 .and. m2max<mpullback/2) then
                   if (maxval(abs(steen_l_space(:,m2max+1:mpullback/2,l2)-1e300))/=0) then
                      print *,'Error in heart of gyrotransformation with',&
                           lmax,oe,moe,(oe+moe+1)/2,'l2',l2,'l',l,'m2max',m2max,'mmax',m2max*2+moe,&
                           'l2min(m2max)=',((m2max)*2+moe+oe+1)/2,&
                           'l2min(m2max+1)=',((m2max+1)*2+moe+oe+1)/2,'mpullback=',mpullback,&
                           steen_l_space(1,m2max+1,l2),&
                           maxval(abs(steen_l_space(:,m2max+1:mpullback/2,l2)-1e300))
                      stop
                   end if
                end if
                call dgemm('n','n',nmaxpoly,1+m2max,nmaxpoly,1.,colmat(:,:,l),ncolmat,&
                     steen_l_space(:,0:m2max,l2),nmaxpoly,0.,&
                     tmp_steen_space(:,0:m2max),nmaxpoly)
                call dcopy(nmaxpoly*(1+m2max),tmp_steen_space,1,steen_l_space(:,0:m2max,l2),1)
                ! steen_l_space(:,0:m2max,l2)=tmp_steen_space(:,0:m2max)
             end do
             call cpu_time(t2)
             t(11)=t(11)+(t2-t1)
             m2max=min(lmax/2-(oe+moe+1)/2,mpullback/2)
             cost(11)=cost(11)+nmaxpoly**2*(&
                  ((m2max+1)*m2max)/2+&
                  (lmax/2-(oe+moe+1)/2-m2max)*(mpullback/2+1))

!!$                 !Achtung:
!!$                 !vergleiche:
!!$                 ksh@gadget:~/w> ./gt 5 5 0 1 1 .1
!!$ Dabei ist hier der Fehler sogar viel zu klein!!!
!!$ evtl. könnte man bei nmax noch eine maximal rechts lokalisierte Funktion ausrechnen.
!!$                 und was nicht geht: zu großer error in den diagonalen.
!!$                 ksh@gadget:~/w> ./gtt 5 5 0 1 1 .1

             mloop2:  do m=mpullback-mod(moe+mpullback,2),moe,-2 ! Leaving out negative m's here. That must be accounted for when doing sums.

                ! ..... and then:
                if (verbose>4) v_l_space=1e300
                call cpu_time(t1)
                call dgemm('t','n',nsteen,(lmax-m+2-oe)/2,nmaxpoly,1.,projsteen,nsteen&
                     ,steen_l_space(:,m/2,(m+oe+1)/2),nmaxpoly*(1+mpullback/2),0.,&
                     v_l_space(:,(m+oe+1)/2:lmax/2),nsteen)
                call cpu_time(t2)
                t(4)=t(4)+(t2-t1)
                cost(4)=cost(4)+nsteen*((lmax-m+2-oe)/2)*nmaxpoly

                if (verbose>4) v_theta_space=1e300
                call cpu_time(t1)
                call dgemm('n','t',nsteen,lmax/2,(lmax-m+2-oe)/2,1.,&
                     v_l_space(:,(m+oe+1)/2:lmax/2),nsteen,&
                     projassleg(:,(m+oe+1)/2:lmax/2,oe,m+1),lmax/2,0.,&
                     v_theta_space,nsteen)
                call cpu_time(t2)
                t(5)=t(5)+(t2-t1)
                cost(5)=cost(5)+lmax/2*((lmax-m+2-oe)/2)*nsteen

                call cpu_time(t1)
                if (krb/=krhoa) then
                   v_theta_space=v_theta_space*besklb(:,:,m)
                else
                   v_theta_space=v_theta_space*beskla(:,:,m)
                end if
                call cpu_time(t2)
                t(6)=t(6)+(t2-t1)
                cost(6)=cost(6)+nsteen*lmax/2

                if (verbose>4) steen0_theta_space=1e300
                call cpu_time(t1)
                call dgemm('n','n',nmax0,lmax/2,nsteen,1.,projsteen,nsteen,v_theta_space,nsteen,0.,steen0_theta_space,nmax0)
                call cpu_time(t2)
                t(7)=t(7)+(t2-t1)
                cost(7)=cost(7)+nsteen*lmax/2*nmax0
                !             ! SUM UP gyrocolmat (!) ==>
                call cpu_time(t1)
                if (m==0) then
                   !projleg delivers 1<=l<=lmax0. thereby oe=1: 1,3,... oe=2: 2,4,...  max(l)=lmax0-mod(oe+lmax0,2)
                   ! n(l)=(max(l)-oe)/2+1=(lmax0-mod(oe+lmax0,2)-oe+2)/2=(lmax0-oe-mod(lmax0-oe,2)+2)/2=(lmax0-oe+2)/2 correct.
                   call dgemm('n','n',nmax0,(lmax0-oe)/2+1,lmax/2,4.,steen0_theta_space&
                        ,nmax0,projleg(:,oe:),lmax,1.,gyrocolmat(:,oe:,i,j),nmax0*2)
                   !factor 8 because of lmax/2 and in total 4 L projectors in a row.
                   !note dimensions: projleg(lmax/2,lmax0) gyrocolmat(nmax0,lmax0,nmax0,lmax0)
                   !note further: parity conservation requires gyrocolmat(_,i,_,j) to be zero for odd i+j.
                   !lmax and nmax0*2 because of parity conservation
                else
                   call dgemm('n','n',nmax0,(lmax0-oe)/2+1,lmax/2,8.,steen0_theta_space,nmax0,&
                        projleg(:,oe:),lmax,1.,gyrocolmat(:,oe:,i,j),nmax0*2)
                endif
                call cpu_time(t2)
                t(8)=t(8)+(t2-t1)
                cost(8)=cost(8)+lmax/2*nmax0*((lmax0-oe)/2+1)


                ! Then we transform everything back into normal space.
                ! steen_l_space is the nice pullbacked matrix.
                ! In principle the full specification is
                ! steen_l_space(1..n,m/2,l=m..lmax [,i=1..nmax0,j=1..lmax0,m=1..mpullback+1])   l>=m
                !    [..] suppressed (loop variables)
                ! it may be useful to actually caclulate this one.
                ! ok, we probably don't have to store the m's simultaneously.
                ! the outer loop should be over m.
                ! totally clear, that the m sum is outside. and we store the bessel_jn(m, ... k ... l)
                ! independent from the rest.
                ! with that the storage is reduced to n*lmax for bessel, n*nmax0 for projsteen
                ! lmax*lmax0 for projleg and (m*) lmax*mpullback.
             enddo mloop2
          end do moeloop
       end do iloop
    enddo jloop
    call cpu_time(t1)
    if (verbose>4) then
       !find index of max absolute element of gyrocolmat:
       idx=idamax((nmax0*lmax0)**2,gyrocolmat,1)
       a=Mod((idx-1)/(/nmax0*lmax0**2,nmax0*lmax0,lmax0,1/),(/nmax0,lmax0,nmax0,lmax0/))+(/1,1,1,1/)
       v1=gyrocolmat(a(1),a(2),a(3),a(4))
       if (abs(v1)>val) then
          val=abs(v1)
          print "(A,I3,A,4I4,A,ES25.17)",'At m=',m,' component maximum at i,j,k,l=',a,' value=' ,v1
       end if
       m1=0
       do i=1,nmax0
          do j=1,lmax0
             do k=1,nmax0
                do l=1,lmax0
                   if (mod(j+l,2)==1) then
                      if (gyrocolmat(i,j,k,l)/=0) then
                         print *,'gyrocolmat parity??',i,j,k,l,gyrocolmat(i,j,k,l)
                         m1=m1+1
                      end if
                   else
                      if (gyrocolmat(i,j,k,l)==0) then
                         if (m==0 .or. krhoa /=0 .and. krb/=0) then
                            !in case for some reason kperp=0 we do not print the m>0.
                            print *,'gyrocolmat zero??',i,j,k,l,gyrocolmat(i,j,k,l)
                            m1=m1+1
                         end if
                      else if (abs(gyrocolmat(i,j,k,l))<1e-16*epsgyrocolmat) then
                         print '(A,4I5,G25.16)','extraneous (near) zero at',i,j,k,l,gyrocolmat(i,j,k,l)
                      end if
                   end if
                   if ((i<k .or. i<=k .and. j<l) .and. abs(gyrocolmat(i,j,k,l)-gyrocolmat(k,l,i,j))>epsgyrocolmat*(1 &
                        +abs(gyrocolmat(i,j,k,l)))) then
!!$                           if(q<o) print 2,'gyrocolmat symmetry??',m,i,j,k,l,gyrocolmat(i,j,k,l),gyrocolmat(k,l,i,j) &
!!$                                ,gyrocolmat(i,j,k ,l)-gyrocolmat(k ,l,i,j)
!!$                           q=q+1
                      m1=m1+1
                   endif
                end do
             end do
          end do
       end do
       if (m1/=0) print 2,'gyrocolmat at m=',m,' had',m1,' issues.'
    end if
    call cpu_time(t2)
    t(9)=t(9)+(t2-t1)

    call cpu_time(t2)
    tto=t2-tto
    if (verbose>0) then
       print '(A,G24.16)','tbes:',tbes
       cost(9)=1
       do i=1,11
          print "(I2,9(' ',A,ES9.2))",i,'t(i):',t(i),'cost(i):',cost(i),'t(i)/cost(i):',t(i)/cost(i),'flops/cycle(i):',2*cost(i)&
               /(t(i)*3.9d9)
       end do
       print '(A,2(G24.16))','tot',tto,tbes+sum(t)
    end if
    if (verbose>4) then
       allocate(maxn(nmax0:nmaxpoly),maxl(lmax0:lmax))
       do i=nmax0,nmaxpoly
          maxn(i)=maxval(max_n_comp(i,:,:))
       end do
       do i=lmax0,lmax
          maxl(i)=maxval(max_l_comp(i,:,:))
       end do
       print '(A,10(ES9.2))','maxl',maxl
       print '(A,10(ES9.2))','maxn',maxn
       ! e.g. ./gt 5 .1 1 1 1 100
       deallocate(maxn,maxl,max_n_comp,max_l_comp)
    end if
    ! Now go through gyrocolmat and check for larger than normal discrepancies
    if (verbose>4) then
       o=20
       q=0
       div=0
       siv=0
       sym=0
       do i=1,nmax0
          do j=1,lmax0
             do k=1,nmax0
                do l=1,lmax0
                   if (mod(j+l,2)==1) then
                      if (gyrocolmat(i,j,k,l)/=0) then
                         print *,'gyrocolmat parity??',i,j,k,l,gyrocolmat(i,j,k,l)
                      end if
                   else
                      if (gyrocolmat(i,j,k,l)==0) then
                         !if (abs(j-l)<=2 .and. abs(i-k)<=5) &
                         print *,'gyrocolmat zero??',i,j,k,l,gyrocolmat(i,j,k,l)
                         ! this vanishes if xmax<=3:
                         ! ksh@napoleon:~/w> ./gt 15 5 2 1 1 6
                         ! ksh@napoleon:~/w> ./gt 5 30 2 1 1 1
                         ! produces zero??
                         ! ksh@napoleon:~/w> ./gt 15 3 2 1 1  10
                         ! ksh@napoleon:~/w> ./gt 5 1 2 1 1 30
                         ! produces nothing.

!!$                   else if (abs(gyrocolmat(i,j,k,l))<epsgyrocolmat) then
!!$                      print 2,'extraneous (near) zero at',i,j,k,l,gyrocolmat(i,j,k,l)
                      end if
                   end if
                   if ((i<k .or. i<=k .and. j<l)) then
                      sym1= abs(gyrocolmat(i,j,k,l)-gyrocolmat(k,l,i,j))
                      if (sym1>sym) then
                         sym=sym1
                         if(q<o .and. sym>epsgyrocolmat) &
                              print *,'gyrocolmat symmetry??',i,j,k,l,gyrocolmat(i,j ,k,l)&
                              ,gyrocolmat(k,l,i,j),gyrocolmat(i,j,k ,l)-gyrocolmat(k,l,i,j)
                         q=q+1
                      end if
                   endif
!!$                     if (i==k .and. j==l) then
!!$                        if (abs(gyrocolmat(i,j,k,l)-1)>div) then
!!$                           div=abs(gyrocolmat(i,j,k,l)-1)
!!$                           a1=(/i,j,k,l/)
!!$                           print '("gyrocolmat diagonal error",i5,i3,ES24.16)',i,j,gyrocolmat(i,j,k,l)-1
!!$                        endif
!!$                     else
!!$                        if (abs(gyrocolmat(i,j,k,l))>siv) then
!!$                           siv=abs(gyrocolmat(i,j,k,l))
!!$                           a2=(/i,j,k,l/)
!!$                           print '("gyrocolmat side error",4i3,ES24.16)',i,j,k,l,gyrocolmat(i,j,k,l)
!!$                        end if
!!$                     endif
                end do
             end do
          end do
       end do
!!$         if (div>epsp*4 .or. siv>epsp*4) then
!!$            print 2,'something went wrong'
!!$         else
!!$            print 2,'everything ok'
!!$         endif
       print "(A,ES12.2)",'max symmetry error',sym
!!$         print "(A,ES12.2,2I3,ES25.16)",'max error on diagonal',div,a1(1:2),gyrocolmat(a1(1),a1(2),a1(3),a1(4))
!!$         print "(A,ES9.2,4I3,2ES25.16)",'max error on sides',siv,a2&
!!$              ,gyrocolmat(a2(1),a2(2),a2(3),a2(4)),gyrocolmat(a2(3),a2(4)&
!!$              ,a2(1),a2(2))
!!$       do i=2,nmax0-1
!!$          do j=1,lmax0
!!$             print 2,'testm',i,j,gyrocolmat(i-1:i+1,j,i,j)
!!$          end do
!!$       enddo
!!$       do i=2,nmax0-1
!!$          do j=3,lmax0
!!$             print 2,'testmj-2',i,j,gyrocolmat(i-1:i+1,j-2,i,j)
!!$          end do
!!$       enddo
!!$       do i=2,nmax0-1
!!$          do j=1,lmax0-2
!!$             print 2,'testmj+2',i,j,gyrocolmat(i-1:i+1,j+2,i,j)
!!$          end do
!!$       enddo
    end if
  end subroutine gyrotrafo

  subroutine calc_projassleg(projassleg,ng,lmax2,mmax,gpl,gwl)
    !Calculate associated Legendre functions.
    !They are normalised, such that their square integrals over [-1,1] is 1.
    !Due to known parity, normally only half the values are needed.
    !Note also: P-ml=Pml*(-)^m for this normalisation.
    implicit none
    integer,intent(in) :: lmax2,mmax,ng
    real,intent(out) :: projassleg(ng,lmax2,2,mmax)
    real,intent(in) :: gpl(ng),gwl(ng)
    ! ^^may be defined also for more or less than lmax2, no problemo

    real p1,pi,p0,x
    integer i,j,k,l,l2,m,oe
    integer lphys,mphys
    integer m1,m2,q,o !only for verbose

    !allocate(projassleg(lmax2,lmax2,2,mmax))
    !use projassleg(i_vertex,l2=jphys/2+1,oe,mphys+1)
    !here: m=mphys+1, j=jphys+1: (later m=mphys; mpullback=phys.)
    !oe=mod(j+m,2)+1=mod(jphys+mphys,2)+1
    !  => jphys=jphys/2*2+mod(jphys,2)=2*(l2-1)+mod(oe-1+mphys,2)=2*l2-2+mod(oe+1+mphys,2)=
    !          =2*l2-2+1-mod(oe+mphys,2)=2*l2-1-mod(oe+mphys,2)
    ! ==> j=2*l2-mod(oe+mphys,2)
    
    !projassleg-old(i,j,m)=projassleg(i,(j-1)/2+1,mod(j+m,2)+1,m)
    !                     =projassleg(i,(j+1)/2,mod(j+m,2)+1,m)
    !projassleg(i,j2,oe,m)=projassleg-old(i,j2*2-mod(j,2),m)
    !                     =projassleg-old(i,j2*2-mod(oe+m-1,2),m)
    !                     =projassleg-old(i,j2*2+mod(oe+m,2)-1,m)
    !projassleg-old is even, if j+m is even, or if oe is 1 (and not 2).

    !old:
    !projassleg(i,j,m)=Pmj(gpl(i))* sqrt(gwl(i))*sqrt(   (2lphys+1) (lphys-mphys)! /(2(lphys+mphys)!) )
    ! so that int(x=-1..1) projassleg**2/gwl dx=1 or sum_i projassleg**2 =1.
    ! possibly not all Gauss points need be calculated due to parity, then ng is reduced to about half.
    !  
    ! nor(j,m)=sqrt(   (2j-1) (j-m)! /(2(j+m-2)!) )
    ! nor(m,m)=sqrt( (2m-1)/(2 (2m-2)!)
    ! nor(m,m)/nor(m-1,m-1)= sqrt((2m-1)/(2m-3)/(2m-2)/(2m-3))=sqrt((2m-1)/(2m-2)) / (2m-3)
    !            =sqrt((2m-1)(2m-4)/(2m-3)/(2m-2))
    ! nor(j,m)/nor(j-1,m)=sqrt( (2j-1)/(2j-3)*(j-m)/(j+m-2) )
    ! projassleg might be improvable by using an eigenvector routine or real*16 variables
    
    m1=0;m2=0;q=0;o=20 !only for verbose>4
    projassleg=0
    j=0
    do i=1,ng
       x=gpl(i)
       do m=1,mmax
          ! we start at pmm and move to higher L
          p1=0
          mphys=m-1 ! the mphys we are calculating now.
          if (m==1) then
             !P00=1*sqrt(1/2)
             pi=sqrt(.5*gwl(i))
             projassleg(i,1,1,1)=pi
          else
             pi=-sqrt(1-x**2)*projassleg(i,m/2,1,m-1)*sqrt((2*m-1.)/(2*m-2))
             ! P(l+1,l+1)=-(2l+1)*sqrt(1-x^2)*P(l,l) (l=lphys)
             ! Normalised:
             ! P(l+1,l+1)=-(2l+1)*sqrt(1-x^2)*P(l,l)*sqrt((2l+3)/[(2l+1)*(2l+2)(2l+1)])
             !    = -sqrt(1-x^2)*P(l,l)*sqrt((2l+3)/(2l+2))
             !here: lphys=mphys=m-1 => l2=(m+1)/2=(lphys+2)/2=lphys/2+1
             projassleg(i,(m+1)/2,1,m)=pi
          endif
          do l=m+1,2*lmax2
             lphys=l-2
             p0=pi
             pi=( (2*lphys+1)*x*pi-(lphys+mphys)*p1 )/(lphys-mphys+1)
             ! for (m,l)=phys
             ! (l-m+1)P(m,l+1)=(2l+1)xP(m,l)-(l+m)P(m,l-1)
             ! Normalised:
             ! Pnor(m,l)=P(m,l)*sqrt((2l+1)(l-m)!/[2*(l+m)!])
             ! (l-m+1)Pn(m,l+1)/sqrt((2l+3)(l+1-m)/(l+m+1))=
             !   =(2l+1)xPn(m,l)/sqrt(2l+1)-(l+m)Pn(m,l-1)/sqrt((2l-1)(l+m)/(l-m))
             ! Pn(m,l+1)*sqrt((l-m+1)*(l+m+1)/(2l+3))=
             !   =xPn(m,l)*sqrt(2l+1)-Pn(m,l-1)*sqrt((l+m)(l-m)/(2l-1))
             ! 
             pi=pi*sqrt( (2*l-1)/(2*l-3.)* ((l-m)/(l+m-2.)) )
             p1=p0*sqrt( (2*l-1)/(2*l-3.)* ((l-m)/(l+m-2.)) )
             projassleg(i,(l+1)/2,mod(l+m,2)+1,m)=pi
          enddo
          if (verbose>4) then
             do l2=1,lmax2
                do oe=1,2
                   if (l2*2+mod(oe+m,2)-1<m) cycle
                   j=j+1
                   if (projassleg(i,l2,oe,m)==0) then
                      print *,'projassleg zero??',i,x,l2,oe,m,projassleg(i,l2,oe,m)
                      m1=m1+1
                   else if(abs(projassleg(i,l2,oe,m))<epsgyrocolmat) then
                      !print 2,'projassleg near zero',i,x,l2,oe,m,projassleg(i,l2,oe,m)
                      m2=m2+1
                   end if
                end do
             end do
          end if
       enddo
    enddo
    if (verbose>4) print '(8(A,I3))','projassleg had',m1,'zero and',m2,'near zero coefficients of',j,'total',(mmax)*ng*lmax2*2
    ! Now should better check the normalisation
!!$  m=(lmax+1)/2
!!$  m=mpullback+1
!!$  do l=m,lmax
!!$     print 2,'m,l norm=1',m,l,sum(projassleg(:,l,m)**2)
!!$  enddo
!!$  do l=m,lmax-1
!!$     print 2,'m,l norm=0',m,l,sum(projassleg(:,l,m)*projassleg(:,l+1,m))
!!$  enddo
!!$ 2   format (5(G0,"  "))
  end subroutine calc_projassleg
  
  subroutine calc_projasslegm(projasslegm,ng,lmax,mmode,gpl,gwl)
    !Calculate associated Legendre functions.
    !Calculate only one specific m-mode with l=mmode .. lmax
    !They are normalised, such that their square integrals over [-1,1] is 1.
    !That means their RMS is 1/sqrt(2) (over the sphere or over [0,1] or over [-1,1]).
    !Note also: P-ml=Pml*(-)^m for this normalisation.
    !Due to known parity, normally only half the values are needed.
    implicit none
    ! do not distinguish any longer between odd and even momenta. To difficult here.
    integer,intent(in) :: lmax,mmode,ng
    real,intent(out) :: projasslegm(ng,lmax-mmode+1)
    real,intent(in) :: gpl(ng),gwl(ng)
    ! ^^may be defined also for more than lmax2, no problemo

    real p1,pi,p0,x
    integer i,j,k,l,m
    integer lphys,mphys
    integer m1,m2,q,o !only for verbose>4

    if (lmax > 2*ng) then
       print '(A,I3,A,I3,A)','calc_projasslegm error 2ng=',2*ng,'lmax=',lmax,'-> scalar products incorrect.'
       stop
    end if
    
    m1=0;m2=0;q=0;o=20;j=0 !only for verbose
    do i=1,ng
       x=gpl(i)
       !       do m=1,mmax
       m=mmode
       ! we start at pmm and move to higher L
       p1=0
       ! for phys l:
       ! Pll=(2l-1)!! (-sqrt(1-x^2))**(l)
       ! Necessary normalisation factor: sqrt[ (2l+1) (l-m)!/[2(l+m)! ]
       !     = sqrt[ (2l+1)/[2(2l)! ]
       pi=sqrt(.5*gwl(i))
       do m=2,mmode
          pi=-sqrt(1-x**2)*pi*sqrt((2*m-1.)/(2*m-2))
       end do
       projasslegm(i,1)=pi
       m=mmode
       mphys=m-1
       p1=0
       do l=mmode+1,lmax !l: l idx of poly to be calculated
          lphys=l-2 !lphys of poly. stored in pi currently
          p0=pi
          pi=( (2*lphys+1)*x*pi-(lphys+mphys)*p1 )/(lphys-mphys+1)
          ! for (m,l)=phys
          ! (l-m+1)P(m,l+1)=(2l+1)xP(m,l)-(l+m)P(m,l-1)
          ! Normalised:
          ! Pnor(m,l)=P(m,l)*sqrt((2l+1)(l-m)!/[2*(l+m)!])
          ! (l-m+1)Pn(m,l+1)/sqrt((2l+3)(l+1-m)/(l+m+1))=
          !   =(2l+1)xPn(m,l)/sqrt(2l+1)-(l+m)Pn(m,l-1)/sqrt((2l-1)(l+m)/(l-m))
          ! Pn(m,l+1)*sqrt((l-m+1)*(l+m+1)/(2l+3))=
          !   =xPn(m,l)*sqrt(2l+1)-Pn(m,l-1)*sqrt((l+m)(l-m)/(2l-1))
          ! 
          pi=pi*sqrt( (2*l-1)/(2*l-3.)* ((l-m)/(l+m-2.)) )
          p1=p0*sqrt( (2*l-1)/(2*l-3.)* ((l-m)/(l+m-2.)) )
          projasslegm(i,l-mmode+1)=pi
       enddo
       if (verbose>4) then
          do l=mmode,lmax
             j=j+1
             if (projasslegm(i,l-mmode+1)==0) then
                print '(A,I3,G25.16,I3,I3,G25.16)','projasslegm zero??',i,x,l,mmode,projasslegm(i,l-mmode+1)
                m1=m1+1
             else if(abs(projasslegm(i,l-mmode+1))<epsgyrocolmat) then
                !print 2,'projasslegm near zero',i,x,l,mmode,projasslegm(i,l-mmode+1)
                m2=m2+1
             end if
          end do
       end if
    enddo
    if (verbose>4) then
       print '(8(A,I3))','projasslegm had',m1,'zero and',m2,'near zero coefficients of',j,'total',ng*(lmax-mmode+1)
       ! Now should better check the normalisation
       ! ACTUALLY: Norm is not so accurate and looses 4 digits. Maybe there is a
       ! numerically better way to calculate these polynomials
       do l=max(mmode,lmax-5),lmax
          print '(A,2I3,G25.16)','l norm=1',mmode,l,sum(projasslegm(:,l-mmode+1)**2)
       enddo
       do l=max(mmode,lmax-7),lmax-2
          print '(A,2I3,G25.16)','l norm=0',mmode,l,sum(projasslegm(:,l-mmode+1)*projasslegm(:,l+2-mmode+1))
       enddo
    end if
  end subroutine calc_projasslegm

  subroutine gyroproj(fout,nmax0,lmax0,lmode,mmode,projsteen,sp,nsteen,xmax,krho,eps)
    ! Gyro-project a given distribution function with given lmode and mmode into the gyrokinetic
    ! space. This is useful for creating particle, heat, momentum sources etc..  Different from
    ! gyrotrafo, this time the normalisation of the orthogonal angular functions matters,
    ! because they are returned.  We want the angular coefficients to be equal to the angular
    ! RMS they cause, when transformed into real space. That way, the exact geometrical area of
    ! validity (sphere, [-1,1], [0,2pi]) does not matter.
    ! At krho=0 we expect fout to be exactly the unit matrix in the radial polynomial index for
    ! lmode=l and mmode=1 (physical mode numbers by 1 smaller).
    ! mphys=mmode-1
    ! For k<0 we get output=(-1)^mphys*output for -k due to the Bessel function.
    ! Negative charge corresponds to negative krho.
    ! Not sure whether the code works for k<0.
    ! For mphys<0 we have Yl,-m=(-1)^m Yl,m^* = (-1)^m [Yl,m*Exp(-imphi)] exp(imphi)
    ! J_-m=(-1)^m J_m
    ! Pl,-m=(-1)^m Pl,m ==> Identical output.


    implicit none
    integer,intent(in) :: nmax0,lmax0,lmode,mmode,nsteen
    real,intent(in) :: xmax,krho,eps
    !for a and b species
    optional :: eps

    ! In the first iteration we use equal nmax0 and nmaxpoly for all species.
    ! For unequal n we also need to update landau.f90
    ! In that case we need two projsteen and sp, and beskl (see below)
    ! We may however compute the actually needed nmaxpoly and possibly nmax0.

    real,intent(out) :: fout(nmax0,lmax0,nmax0) ! result is real,
    ! since we multiply with the exp(-i*m*theta)-term of the Bessel expansion of a plane wave.

    real,intent(in) :: projsteen(nsteen,nsteen),sp(nsteen)
    !projection matrix (i.e. sqrt(weight(i))*P(sp(i)) and vertices.
    real,dimension(:),allocatable :: gpl,gwl !,tmp_steen_space ! Gauss points
    real,dimension(:,:),allocatable :: projleg,pltest,projasslegm,v_theta_space,beskl,steen0_theta_space
    real,dimension(:,:,:,:),allocatable :: paltest
    real val
    real epsp
    integer mpullback,lmax,extradegree,nmaxpoly
    integer i,j,k,l,l2
    integer oe !odd even angular momentum =1: even =2:odd.
    integer mphys  ! real mode numbers, otherwise indices are used, which are those +1
    integer o,q !used for verbose output in projassleg
    integer ng2 ! half the Gauss points for Legendre.
    logical timing
    integer, external ::idamax
    integer ngs

    ! originally in block "timeblock":
    real tbes,t(11),tto,t1,t2,cost(11)

    ! originally in block "lminblock"
    integer lminb
    real lbval,lminbval

    timing=verbose>0


    if (present(eps)) then
       epsp=eps
    else
       epsp=epspullback
    endif

    ! First a few sanity checks and estimates

1   format (9(A))
    if (lmode>lmax0) then
       print 1,'Warning in gyroproj: lmode is larger than lmax0. ',&
            'Projection will work, but result does not at all represent the original.'
    end if

    if (mmode>lmode .or. mmode<1 .or. lmode<1) then
2      format (A,30(I3))
       print 2,'Illegal combination of lmode,mmode:',lmode,mmode
       stop
    end if

    mpullback=est_mpullback(krho*xmax,eps=eps**2)
    if (verbose>0) print 2,'mpullback must be >=',mpullback
    !mpullback is what the plane wave expansion requires in terms of m and l.
    !for gyroprojection we only need m=mmode, but l=mmode ... mpullback+1
    !from the plane wave expansion.

    if (mpullback+1<mmode) then
       ! can never produce any non-zero m=0 component.
       fout=0
       return
    end if

    ! Following short discussion in physical l,m numbers:
    ! For a given real space mode (l,m) the product with a mode (lp,mp) in the plane wave expansion
    ! can only produce something with m=0 if m=-mp.
    ! Parity conservation requires for the output mode mod(lo,2)=mod(l+m,2),
    ! since the parity of the plane wave is even.
    ! Since the parity of all plane wave expansion modes is even, mod(lp+mp,2)=0.
    ! The maximum considered lp and mp is mpullback.
    ! lo<=l+lp-mod(lp,2) (parity conservation and lp constraint)
    ! lo>=|l-lp|+mod(lp,2)
    ! lo>=mod(m+l,2) due to parity
    ! m<=lp<=mpullback
    ! ==> lo<=l+mpullback-mod(mpullback,2) and that is sharp if m<=mpullback.
    ! >>> incomplete discussion <<< needs to be corrected/refined later.


    ! The parity of the output mode (lo,0) is mod(lo,2) and must be the same as the one
    ! of the input mode, mod(l+m,2) since the plane wave is even.

    ! The product of lmode,mmode and lplanewave,mmode, lplanewave<=mpullback+1, and
    ! mod(lplanewave+mmode,2)=0 since the plane wave has even parity, has parity
    ! mod(lmode+mmode+lplanewave+mmode,2)=mod(lmode+lplanewave,2)=mod(lmode+mmode,2).
    ! The minimum l index is potentially abs(lmode-lplanewave)+1, but only if it has the
    ! same parity as (lmode,mmode), i.e., if mod(abs(lmode-lplanewave)+lmode+mmode,2)=0.
    ! <=> mod(lplanewave+mmode,2)=0.
    ! That means only l indices (lphys+1) from
    !     abs(lmode-lplanewave+mod(lplanewave+mmode,2))+1
    !    to (lmode+lplanewave-mod(lplanewave+mmode,2)-1)
    ! in steps of 2 can result in fout.

    ! So the maximum l index at all, that could result is:
    lmax=lmode+mpullback
    ! The minimum l index at all is max(0,lmode-lplanewave+1)
    if (verbose>0) print 2,'Minimum physical L after projection is',max(0,lmode-mpullback)
    if (lmax0>lmax) then
3      format (A,20(I3,A))
       print 3,'warning gyroproj: lmax0=',lmax0,'>lmax=',lmax,'=lmode(',lmode,&
            ')+mpullback(',mpullback,')'
       print 3,'=> the coefficients from lmax+1 .. lmax0 will be < epspullback. Reduce lmax0.'
    end if
    !the maximum m for projassleg is actually mmode+mpullback and not lmax.
    !The parity of the mode is mod(lmode+mmode,2) (0=even 1=odd)
    !The parity of a mode in the plane wave expansion is
    ! make lmax an even number
    if (mod(lmax,2)==1) then
       mpullback=mpullback+1
       lmax=lmode+mpullback
       print 2,'Corrected mpullback to make lmax even mp.=',mpullback
    end if
    print 3,'Using lmax=',lmax,'. (lmode=',lmode,'lmax0=',lmax0,')'

    extradegree=est_extradegree(abs(krho)*xmax,eps=eps)

    nmaxpoly=nmax0+extradegree
    print 3,'Need extradegree',extradegree,'-> max. polynomial degree=',nmaxpoly-1

    if (nmaxpoly>nsteen) then
       print 3,'Error in gyroproj: Max. polynomial degree in projsteen/sp='&
            ,nsteen-1,', but need',nmaxpoly-1
       stop
    end if


    ng2=max(lmax/2,(lmax0+1)/2)
    allocate(gwl(ng2*2),gpl(ng2*2))
    call gauss_legendre(-1.,1.,gpl,gwl,ng2*2)
    ! If lmax is even, we need only half of the vertices.
    !Therefore often only half sums have to be calculated.
    ! And we assume lmax to be even.
    allocate(projleg(ng2,lmax0)) ! Convert from lmax0-L-space to lmax-xi-space
    call calc_projasslegm(projleg,ng2,lmax0,1,gpl,gwl)
    if (verbose>4) then
       allocate(pltest(ng2,lmax0))
       call calc_projleg(pltest,ng2,lmax0,gpl,gwl)
       print 2,'Max deviation of projleg,pltest:',maxval(abs(projleg-pltest))
    end if

    allocate(projasslegm(ng2,lmax-mmode+1))
    call calc_projasslegm(projasslegm,ng2,lmax,mmode,gpl,gwl)

    if (verbose>4) then
       allocate(paltest(ng2,lmax/2,2,mmode)) ! need actually only one Pml for mmode and lmode.
       !    call calc_projassleg(paltest,lmax/2,(lmode+1)/2,mmode,gpl,gwl)

       call calc_projassleg(paltest,ng2,lmax/2,mmode,gpl,gwl)
       print 5,'Comparison of new projasslegm with old projassleg:'
       ngs=max(1,ng2-4)
       do l=mmode,lmax
5         format (A,I3,A,G25.16)
          print 5,'l=',l,'projasslegm(ngs:ng2,l)=',projasslegm(ngs:ng2,l-mmode+1)
          print 5,'l=',l,'    paltest(ngs:ng2,l)=',paltest(ngs:ng2,(l+1)/2,mod(l+mmode,2)+1,mmode)
          print 5,'l=',l,'      delta(ngs:ng2,l)=',paltest(ngs:ng2,(l+1)/2,mod(l+mmode,2)+1,mmode)&
               -projasslegm(ngs:ng2,l-mmode+1)
       end do
    end if

    allocate(beskl(nsteen,ng2))
    allocate(steen0_theta_space(nmax0,ng2))
    allocate(v_theta_space(nsteen,ng2))


    val=0
    if (verbose>4) then
       q=0;o=20 ! for error messages
    end if
    if (timing) then
       tbes=0; t=0; cost=0
       call cpu_time(tto)
    end if
    mphys=mmode-1 !m is mphys here.

    if (timing) call cpu_time(t1)
    if (verbose>4) beskl=1e300
    do l2=1,ng2
       beskl(:,l2)=bessel_jn(mphys,krho*sp*sqrt(1-gpl(l2)**2))
       ! The plane wave exp(iky)=sum_(m in Z) bessel_jn(m,k|r|)exp(im phi)
       ! r=(x,y)
       ! kx*rho_x=-kx*mv_y/qB for B=Bz positive.
       ! i.e.
       ! exp(i kx rho_x)=exp(-i kx*m/(qB) vy) =sum_m bessel_jn(m,-krho |v|) exp(imphi)
       ! call kx*m/(qB)=krho.
       ! multiplication of real space distribution function with exp(ikx*rhox)
       ! gives the gyrocenter distribution function. (Both Fourier transformed)
       ! Thereafter we determine the l,m=0 components == gyro averaging.
       ! We take the -mphys component of this expansion to get a non-zero contribution
       ! from mphys=(mmode-1)
       ! Bessel_jn(-m,-krho |v|) exp(-imphi)=Bessel_jn(m,krho|v|) exp(-imphi)
       ! Note bessel_jn(-m)=bessel_jn(m)*(-1)^m
       !      Bessel_jn(m,-x)=(-1)^m Bessel_Jn(m,x)=Bessel_jn(-m,x)
       ! must use -m, but this gives just a (-1)^m sign.
       ! don't like sqrt(1-gpl(l)**2) for sin theta. Should use something more accurate. ???
    enddo
    !         print 5,'m=',m,'sumbeskl',sum(beskla(:,1:ng2))
    if (timing) then
       call cpu_time(t2)
       tbes=tbes+(t2-t1)
    end if
    j=lmode
    oe=mod(lmode+mmode,2)+1 !oe=1: even in z oe=2: odd in z
    lminb=0
    lminbval=0
    !m=mphys l=lphys+1 here
    !l>=m+1, l<=lmax and oe+m+l is even.
    !==> l>=m+oe. since lmax is even l<=lmax-mod(oe+m,2)
    !the correct arguments are projassleg(,(l+1)/2,oe,m+1)
    !==> (l+1)/2 >=(m+oe+1)/2  (l+1)/2<=(lmax-mod(oe+m,2)+1)/2=ng2 +(1-mod(oe+m,2))/2=ng2
!!$          do l2=(m+oe+1)/2,ng2 !n_elements=ng2-(m+oe+1)/2+1=(lmax+2-m-oe)/2
!!$             lbval=sum(beskla(nmaxpoly,:)*projleg(:ng2,j)&
!!$                  *projassleg(:,l2,oe,m+1))
!!$             !if (m==1 .and. j==30) print 2,'m,j,l,lbval',m,j,l,lbval
!!$             if (abs(lbval)>epsp) then
!!$                lminb=l2*2+mod(oe+m,2)-1  +1
!!$                lminbval=lbval
!!$             end if
!!$          end do
!!$          print 2,'m',m,'j',j,'lminb',lminb,'lminbval',lminbval
    !^^ can use this boundary value as estimate of max l needed.
    ! most conservative value for m=0.
    ! can also reduce lmax by that. But will not be greatly different from lmax.

    iloop: do i=1,nmax0
       if (timing) call cpu_time(t1)
       do l2=1,ng2
          do k=1,nsteen
             v_theta_space(k,l2)=beskl(k,l2)*projsteen(i,k)*projasslegm(l2,lmode-mmode+1)
             ! Now we are at RMS=1/sqrt(2) due to the normalisation of projasslegm.
             !write down mode ampl*plane wave with -m.
             ! Note: sum_i=-infty,+infty  bessel_jn(i,x)^2 = 1 for any x.
          enddo
       enddo
       if (timing) then
          call cpu_time(t2)
          t(1)=t(1)+(t2-t1)
          cost(1)=cost(1)+nsteen*ng2
          !cost=lmax*n maybe ng2*n
       end if
       if (verbose>4) steen0_theta_space=1e300
       if (timing) call cpu_time(t1)
       call dgemm('n','n',nmax0,ng2,nsteen,1.,projsteen,nsteen,v_theta_space,nsteen,&
            0.,steen0_theta_space,nmax0)
       ! angular RMS unchanged, just transform into orthonormal radial polynomials.
       if (timing) then
          call cpu_time(t2)
          t(7)=t(7)+(t2-t1)
          cost(7)=cost(7)+nsteen*ng2*nmax0
       end if
       if (verbose>4) then
          fout(:,3-oe::2,i)=1e100 ! these are actually zero due to parity and should not be needed later.
          fout(:,oe::2,i)=1e300
       end if
       if (timing) call cpu_time(t1)
       !projleg delivers 1<=l<=lmax0. thereby oe=1: 1,3,... oe=2: 2,4,...  max(l)=lmax0-mod(oe+lmax0,2)
       ! n(l)=(max(l)-oe)/2+1=(lmax0-mod(oe+lmax0,2)-oe+2)/2=(lmax0-oe-mod(lmax0-oe,2)+2)/2=(lmax0-oe+2)/2 correct.
       if (lmax0>=oe) call dgemm('n','n',nmax0,(lmax0-oe)/2+1,ng2,2.,steen0_theta_space,nmax0,projleg(1&
            ,oe),ng2*2,0.,fout(1,oe,i),nmax0*2)
       ! Now we have integrated over the domain [0,1] with a function of RMS=1/sqrt(2).
       ! i.e. ==> with the RMS=1/sqrt(2.) before we are now in principle at 1/2 and must
       ! multiply 2 back.
       ! In other words: multiply a factor of 2 because we have exactly one legendre/ass.-legendre projection.
       if (timing) then
          call cpu_time(t2)
          t(8)=t(8)+(t2-t1)
          cost(8)=cost(8)+ng2*nmax0*((lmax0-oe)/2+1)
       end if
    enddo iloop
    if (timing) then
       call cpu_time(t2)
       tto=t2-tto
       print '(A,ES9.2)','tbes:',tbes
       cost(9)=1
       do i=1,11
          if (t(i)/=0) &
               print "(I2,9(' ',A,ES9.2))",i,'t(i):',t(i),'cost(i):',cost(i),'t(i)/cost(i):',&
               t(i)/cost(i),'flops/cycle(i):',2*cost(i)/(t(i)*3.9d9) ! for I7-3770
       end do
       print '(A,2ES9.2)','tot',tto,tbes+sum(t)
    end if

  end subroutine gyroproj

  subroutine gen_polys(polys,names,unitsenergy,npolys,emat,projsteen,sp,a1,b1,c1,xmax,krho,&
       nsteen,nmax0,lmax0,eps)
    ! generate unit source terms for standard thermondyn quantities and their energies
    ! useful for test programs that calculate the standard fluid transport coefficients.
    use landau, only: vunit,intnorm
    implicit none
    real, dimension(:,:,:), allocatable :: polys
    character(12),dimension(:), allocatable :: names
    real, dimension(:), allocatable :: unitsenergy
    integer npolys,nsteen,nmax0,lmax0
    real :: projsteen(nsteen,nsteen),emat(:,:)
    real,dimension(nsteen) :: sp,a1,b1,c1
    real xmax,krho,eps
    intent (out) :: polys,names,unitsenergy,npolys
    intent (in) :: emat,projsteen,sp,a1,b1,c1,xmax,krho,nsteen,nmax0,lmax0,eps
    optional :: eps

    real, dimension(:,:,:,:,:), allocatable :: fout
    real, dimension(:), allocatable :: den,dens,en,temp,temps,mom,moms,enfl,hfl,hfls,&
         pre,pxx,pxxs,pzz,pzzs,pxx_yy,pxx_yys,pyz,pyzs
    real densener,tempsener,momsener,hflsener,pzzsener,pxxsener,pxx_yysener,pyzsener
    real v,w
    integer l,m

    if (nmax0<4) then
       print '(A)','gen_pols: nmax0 must be >=4, otherwise the heat flux cannot be represented.'
       stop
    end if
    allocate(fout(nmax0,lmax0,nmax0,3,3))
    do l=1,3
       do m=1,l
          if (l==5 .and. m==3) then
             call gyroproj(fout(:,:,:,l,m),nmax0,lmax0,l,m,projsteen,sp,nsteen,xmax,krho,eps)!,verbose)
          else
             call gyroproj(fout(:,:,:,l,m),nmax0,lmax0,l,m,projsteen,sp,nsteen,xmax,krho,eps)
          end if
       end do
    end do
    allocate(den(nmax0),dens(nmax0),temp(nmax0),temps(nmax0),en(nmax0),mom(nmax0),moms(nmax0))
    allocate(enfl(nmax0),hfl(nmax0),hfls(nmax0))
    allocate(pre(nmax0),pzz(nmax0),pxx_yy(nmax0),pyz(nmax0),pzzs(nmax0),pxx_yys(nmax0),pyzs(nmax0))
    ! P(n)=(c1(n)*x-a1(n))*P(n-1)-b1(n)*P(n-2)
    ! c1(n)x*P(n-1)=P(n)+a1(n)*P(n-1)+b1(n)*P(n-2)
    ! c1(n+1)x*P(n)=P(n+1)+a1(n+1)*P(n)+b1(n+1)*P(n-1)
    ! x*P(n)=(P(n+1)+a1(n+1)*P(n)+b1(n+1)*P(n-1))/c1(n+1)
    ! the 3 first v powers:
    den=0
    den(1)=1/c1(1)
    mom=0
    mom(1:2)=[a1(2),1.]*den(1)/c1(2)
    mom=mom*vunit
    en=0
    en(2:3)=mom(1:2)/c1(2:3)
    en(1:2)=en(1:2)+mom(1:2)*a1(2:3)/c1(2:3)
    en(1)=en(1)+mom(2)*b1(3)/c1(3)
    en=.5*en*vunit
    enfl=0
    enfl(2:4)=en(1:3)/c1(2:4)
    enfl(1:3)=enfl(1:3)+en(1:3)*a1(2:4)/c1(2:4)
    enfl(1:2)=enfl(1:2)+en(2:3)*b1(3:4)/c1(3:4)
    enfl=enfl*vunit

    densener=1/(den(1)*emat(1,1)*den(1)*intnorm)
    dens=den*densener  ! -> now <dens|emat|den>=1. <dens|emat|dens=densener
    !dens: source of unit density, den: measurement "bra" vector.
5   format (A,(G23.16,4G25.16))
    print 5,'density energy=',densener,dens(1)**2*emat(1,1)*intnorm
    momsener=1/(dot_product(mom(1:2),matmul(emat(1:2,1:2),mom(1:2)))*2./6.*intnorm)
    ! Factor 2./6: 2 Y functions contribute, each of them 1/sqrt(6)^2
    ! or: x^2 averaged over the sphere is 1/3.
    moms=mom/momsener
    print 5,'momenergy=',momsener,dot_product(moms(1:2),matmul(emat(1:2,1:2),moms(1:2)))*2./6.*intnorm,&
         dot_product(mom(1:2),matmul(emat(1:2,1:2),moms(1:2)))*intnorm
    ! energy of momentum "1" (*2)
    v=dot_product(en(1:3),emat(1:3,1))*dens(1)/densener*intnorm
    temps=en-dens*v
    w=1.5/(dot_product(temps(1:3),matmul(emat(1:3,1:3),en(1:3)))*intnorm)
    temps=temps*w
    print 6,'temps=en*',w,'-dens*',v*w
    tempsener=dot_product(temps(1:3),matmul(emat(1:3,1:3),temps(1:3)))*intnorm ! should be 1.5
    temp=temps/tempsener
    print 6,'tempener',tempsener,'entemps',dot_product(temps(1:3),matmul(emat(1:3,1:3),en(1:3)))*intnorm,&
         'dentemps',den(1)*sum(emat(1,1:3)*temps(1:3))*intnorm
    ! energy of temperature "1" (*2)
    v=dot_product(enfl(1:4),matmul(emat(1:4,1:2),moms(1:2)))/momsener*2./6*intnorm
    hfls=enfl-moms*v
    w=1/(dot_product(hfls(1:4),matmul(emat(1:4,1:4),enfl(1:4)))*2./6.*intnorm)
    hfls=hfls*w
6   format (4(A,G23.16))
    print 6,'hfls=enfl*',w,'-moms*',v*w
    hflsener=dot_product(hfls(1:4),matmul(emat(1:4,1:4),hfls(1:4)))*2./6*intnorm
    hfl=hfls/hflsener
    print 6,'hflener',hflsener,'enfl*hfls',dot_product(hfls(1:4),matmul(emat(1:4,1:4),enfl(1:4)))*2./6*intnorm,&
         'momhfls',dot_product(hfls(1:4),matmul(emat(1:4,1:2),mom(1:2)))*2./6*intnorm
    ! energy of heat flux "1" (*2)

    ! (2zz-xx-yy)/3 stress:
    pre=2*en
    w=1/(dot_product(pre(1:3),matmul(emat(1:3,1:3),pre(1:3)))*4./5/9*intnorm)
    pzzs=pre*w
    pzzsener=dot_product(pzzs(1:3),matmul(emat(1:3,1:3),pzzs(1:3)))*4./5/9*intnorm
    pzz=pzzs/pzzsener
    w=1/(dot_product(pre(1:3),matmul(emat(1:3,1:3),pre(1:3)))*(2./30)*intnorm)
    pyzs=pre*w
    pyzsener=dot_product(pyzs(1:3),matmul(emat(1:3,1:3),pyzs(1:3)))*(2./30)*intnorm
    pyz=pyzs/pyzsener
    pxx_yys=pyzs
    pxx_yysener=pyzsener
    pxx_yy=pyz
    pxxs=pzzs
    pxxsener=pzzsener
    pxx=pzz
    print 6,'pzzener',pzzsener,'pyzener=',pyzsener,'pxx_yyener',pxx_yysener

    npolys=11
    allocate(polys(nmax0,lmax0,npolys),names(npolys),unitsenergy(npolys))

    ! actually it might be less confusing to write (*i) instead of (/i) here.
    names=[character(len=12) :: 'n','T','vy (=1/i)','vz','qy (=1/i)','qz',&
         'pizz','piyz (=1/i)','pixx-yy/2','pixx','piyy']
    unitsenergy=[densener,tempsener,momsener,momsener,hflsener,hflsener,pzzsener,pyzsener,pxx_yysener,pxxsener,pxxsener]
    do l=1,lmax0
       if (mod(l+1,2)==0) then
          ! even parity quantities
          polys(:,l,1)=matmul(fout(:,l,:,1,1),dens)
          polys(:,l,2)=matmul(fout(:,l,:,1,1),temps)
          polys(:,l,3)=matmul(fout(:,l,:,2,2),moms)*2*sqrt(1./6) !*i mphys=+-1
          polys(:,l,4)=0
          polys(:,l,5)=matmul(fout(:,l,:,2,2),hfls)*2*sqrt(1./6) !*i
          polys(:,l,6)=0
          polys(:,l,7)=matmul(fout(:,l,:,3,1),pzzs)*sqrt(4./5)/3
          polys(:,l,8)=0
          polys(:,l,9)=matmul(fout(:,l,:,3,3),pxx_yys)*(2/sqrt(30.)) ! (xx-yy)/2 stress. gyrostress
          polys(:,l,10)=matmul(fout(:,l,:,3,1),pxxs)*(-1/sqrt(5.)/3)+&
               matmul(fout(:,l,:,3,3),pxxs)*(2*sqrt(3./10.)/3) ! xx stress.  gyrostress
          polys(:,l,11)=matmul(fout(:,l,:,3,1),pxxs)*(-1/sqrt(5.)/3)+&
               matmul(fout(:,l,:,3,3),pxxs)*(-2*sqrt(3./10.)/3) ! yy stress.  gyrostress
       else
          ! odd parity quantities
          polys(:,l,1)=0
          polys(:,l,2)=0
          polys(:,l,3)=0
          polys(:,l,4)=matmul(fout(:,l,:,2,1),moms)*(-sqrt(1./3))
          polys(:,l,5)=0
          polys(:,l,6)=matmul(fout(:,l,:,2,1),hfls)*(-sqrt(1./3))
          polys(:,l,7)=0
          polys(:,l,8)=matmul(fout(:,l,:,3,2),pyzs)*2/sqrt(30.)  !*i yz stress odd parity can never interact.
          polys(:,l,9)=0
          polys(:,l,10)=0
          polys(:,l,11)=0
       end if
       !zero is the: xy stress,
    end do

  end subroutine gen_polys
end module gyrotransformation
