  MODULE curden_terms
  USE nrtype, ONLY : DP,I4B
  USE dep_var, ONLY : rbp
  IMPLICIT NONE
  SAVE
  PUBLIC
  REAL(DP),DIMENSION(:),ALLOCATABLE :: curden,curohm,curboot,       &
              curbeam,currf,curdri,bp,q,dqdr,shearp,curpar,             &
              xden,xdet,sl31,sl32,sl34,eta,xjbte,xjbne,xjbnf,       &
              dfion,dfast,curbe,curbet,curbi
  REAL(DP),DIMENSION(:,:),ALLOCATABLE :: xdin,xdit,xjbti,xjbni
  REAL(DP) tot_cur,totohm_cur,totbeam_cur,totrf_cur,totboot_cur,    &
           q0_max,q0_mult,q0_radius,rbsaxis
  INTEGER(I4B) jboot,ibcur,irfc,curtype,beam_cur_comp

  !q on axis control with these parameters is not active in gcnm code:
  DATA  q0_max,q0_mult,q0_radius,curtype /-1.0_DP,1.0_DP,0.0_DP,0_I4B /

  DATA beam_cur_comp /0/  ! individual components of beam driven current not
                          ! avaialable.
  CONTAINS

     SUBROUTINE allocate_cur_parms
 
#ifdef GCNMP
       USE ions_gcnmp,                  ONLY : nion
#elif defined NFREYA
       USE ions_gcnmp,                  ONLY : nion
#else
       USE numbrs,                      ONLY : nion
#endif

       USE grid_class,                  ONLY : nj

       USE common_constants,            ONLY : zeroc

       IMPLICIT NONE
       IF(.NOT. ALLOCATED(xjbte))then
          ALLOCATE(xjbte(nj),xjbne(nj))
          ALLOCATE(xjbni(nj,nion),xjbti(nj,nion))
       ENDIF
       IF(.NOT. ALLOCATED(curpar))ALLOCATE(curpar(nj))
       IF(.NOT. ALLOCATED(xdin))ALLOCATE(xdin(nion,nj-1))
       IF(.NOT. ALLOCATED(xdit))ALLOCATE(xdit(nion,nj-1)) 
       IF(.NOT. ALLOCATED(xden))ALLOCATE(xden(nj-1))      
       IF(.NOT. ALLOCATED(xdet))ALLOCATE(xdet(nj-1)) 
       IF(.NOT. ALLOCATED(sl31))ALLOCATE(sl31(nj-1))
       IF(.NOT. ALLOCATED(sl32))ALLOCATE(sl32(nj-1))
       IF(.NOT. ALLOCATED(sl34))ALLOCATE(sl34(nj-1))
       IF(.NOT. ALLOCATED(bp))ALLOCATE(bp(nj))
       IF(.NOT. ALLOCATED(shearp))ALLOCATE(shearp(nj))  
       IF(.NOT. ALLOCATED(xjbne))ALLOCATE(xjbne(nj))      
       IF(.NOT. ALLOCATED(dfion))ALLOCATE(dfion(nj))          
       IF(.NOT. ALLOCATED(dfast))ALLOCATE(dfast(nj))  
       IF(.NOT. ALLOCATED(xjbnf))ALLOCATE(xjbnf(nj))   
       IF(.NOT. ALLOCATED(curbe))ALLOCATE(curbe(nj))
       IF(.NOT. ALLOCATED(curbet))ALLOCATE(curbet(nj))
       IF(.NOT. ALLOCATED(curbi))ALLOCATE(curbi(nj)) 
       IF(.NOT. ALLOCATED(currf))ALLOCATE(currf(nj))  
       IF(.NOT. ALLOCATED(curdri))ALLOCATE(curdri(nj))  
       IF(.NOT. ALLOCATED(eta))ALLOCATE(eta(nj-1)) 
       IF(.NOT. ALLOCATED(q))ALLOCATE(q(nj)) 
       IF(.NOT. ALLOCATED(dqdr))ALLOCATE(dqdr(nj-1)) 

       IF(.NOT. ALLOCATED(curden))ALLOCATE(curden(nj))
       IF(.NOT. ALLOCATED(curohm))ALLOCATE(curohm(nj))
       IF(.NOT. ALLOCATED(curbeam))ALLOCATE(curbeam(nj))
       IF(.NOT. ALLOCATED(curboot))ALLOCATE(curboot(nj))
       IF(.NOT. ALLOCATED(rbp))ALLOCATE(rbp(nj))

       xdin(:,:)    = zeroc    ; xdit(:,:)  = zeroc
       xjbte(:)     = zeroc    ; xjbne(:)   = zeroc
       xjbni(:,:)   = zeroc    ; xjbti(:,:) = zeroc
       xden(:)      = zeroc    ; xdet(:)    = zeroc
       xjbnf(:)     = zeroc    ; eta(:)     = zeroc
       q(:)         = zeroc    ; shearp(:)  = zeroc
       sl31(:)      = zeroc    ; sl32(:)    = zeroc  
       sl34(:)      = zeroc    ; bp(:)      = zeroc
       curbet(:)    = zeroc    ; curbi(:)   = zeroc
       curohm(:)    = zeroc    ;
       curdri(:)    = zeroc    ; curbe(:)   = zeroc
       curden(:)    = zeroc    ; curpar(:)  = zeroc
       dfion(:)     = zeroc    ; dfast(:)   = zeroc
       curboot(:)   = zeroc   
   
       !curbeam(:) = zeroc  ;  currf(:)   = zeroc these are read in from iterdb
     
       
       RETURN
     END      SUBROUTINE allocate_cur_parms
 

      SUBROUTINE curcalc 
!
!
! ----------------------------------------------------------------------
! --- subroutine calculates the current density ,<Jphi*R0/R>,using Ampere's
! --- law in the form:
!           <Jphi*R0/R> = (1/u0)*(1.0/(h*r)*d(g*h*r*bp)/dr
! --- The parallel current density
!     < Jtotal  dot B/Bt0> = (1/u0)*(1./(F**2H*r)*d(FGHrbp)/dr
! --- it is assumed that d<Jphi*R0/R>/dr = 0 at the magnetic axis.
! --- We enforce this condition below by extrapolating the cubic polynomial
! --- fit of curden near the magnetic axis (see subroutine CUBICEXTRP)
!
! --- input
!  rbp(j)                      j = 1,..nj,rbp(j) = f(j)*g(j)*h(j)*r(j)*bp(j)
!                              in Tesla m
!
!  fcap(j)                     j = 1,2..nj,fcap(j) = f(psilim)/f(psi)
!
!  hcap(j)                     j = 1,2..nj,hcap(j) = fcap(j)/<R0**2/R**2>
!  gcap(j)                     j = 1,2,..nj,gcap(j) = <(grad rho R0/R)**2>
!
!  r(j)                        j = 1,2..nj,r(j) = rho value (in cm)
!
!  nj                          the grid size
!
!
! --- output
!  curden(j)                j = 1,2..nj, the current density,amps/m**2,
!                           on the r(j) grid,with zero gradient at r = 0
!  curpar(j)           parallel current, see below
!
! ------------------------------------------------------------------ HSJ
!

      USE nrtype,                 ONLY : DP,I4B

      USE common_constants,       ONLY :  Permeability

      USE dep_var,                ONLY : rbp
      
      USE grid_class,             ONLY : fcap,gcap,hcap,nj,r



      IMPLICIT  NONE
      INTEGER(I4B) j
      REAL(DP)   ghrbp (nj), drbpdr(nj), dghrbpdr(nj),const

      const = 1.0 /Permeability    ! u0 = Permeability 
!
! --- get g*h*rho*bp
!
      ghrbp(:) = rbp(:)/fcap(:)

!
! --- differentiate ghrbp
!
      CALL difydx(r,ghrbp,dghrbpdr,nj)

!
! --- get curden in amps/m**2
!

      DO j=2,nj
        curden(j) = (dghrbpdr(j)/(hcap(j)*r(j)))*const
      END DO

!
! --- get current density at magnetic axis,consistent with
! --- d<Jphi*R0/R>/dr  = 0 at r=0
!
      CALL cubicextrp (curden(2), curden(3), curden(4), &
                       r(2), r(3), r(4), curden(1), 2)
!     make sure curden is normalized to the total current
!     (uses rbp(nj):
      CALL curnorm
!
!     get the parallel current , curpar,given  by
!     curpar =
!     <J*B/Bt0> = < J*R0/R>/fcap +(1/u0)*(gcap/fcap**2)*d/dr(fcap)
!

      CALL difydx (r, rbp, drbpdr, nj)
      DO j=2,nj
        curpar(j) = const/(fcap(j)*fcap(j)*hcap(j)*r(j)) &
                                                       *drbpdr(j)
      END DO
      curpar(1)=curpar(2)
 
      RETURN
!
      END SUBROUTINE curcalc
 




      SUBROUTINE curnorm
! ----------------------------------------------------------------------------
!     normalize curden  to total toroidal current
!     NOTe that I(rho(J)) = Integral_0^rho(j) { 2 pi H rho drho < J Ro/R> } 
!                      = (2pi/u0) *G H rho Bp0 = (2pi/u0)* rbp(j)
!     recall  that curden is defined as <J R0/R> 
!-----------------------------------------------------------------------------
 
      USE nrtype,                             ONLY : DP,I4B

      USE grid_class,                         ONLY : nj,r,hcap

      USE common_constants,                   ONLY : twopi,Permeability

      USE dep_var,                            ONLY : rbp

      IMPLICIT  NONE

      REAL(DP) xdum(nj),ydum(nj)
      REAL(DP) xnorm

         xdum(:)  = hcap(:) * r(:) * curden(:)
         CALL trap2(r,xdum,ydum,nj)                ! here ydum(nj) = I_o/2pi
                                                   ! where I_o is the total 
                                                   ! current  obtained by
                                                   ! integrating up curden

         xnorm  = rbp(nj)/(Permeability*ydum(nj))  ! now modify curden so that it
                                                   ! yields the  desired current I_b
                                                   ! Note that the boundary 
                                                   ! condition input value of I_b 
                                                   ! is given by 
                                                   ! I_b = (2 pi/u0) * rbp(nj) 
                                                 
         curden = curden * xnorm
         tot_cur = twopi *rbp(nj) / Permeability
      RETURN
      END       SUBROUTINE curnorm



  END MODULE curden_terms
