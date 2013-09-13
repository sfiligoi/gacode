  MODULE paleocl
! --------------------------------------------------------------------------
! Compute  the paleoclassical transport coeffiecient for 
! electron heating. See Callen,"Paleoclassical transport in low collisionality
! toroidal plasmas" to be published. (UW-CPTC 04-3(revised))
! -----------------------------------------------------HSJ-04/09/07----------



  USE param,     ONLY : kj
  USE aid_newton,ONLY : freeze_xte
  USE constnts,  ONLY : pi,rot2pi,charg4,xmasse,charge,u0,root2,cee
  USE extra,     ONLY : q
  USE geom ,     ONLY : hcap
  USE io,        ONLY : ddebug,nout
  USE ions,      ONLY : zeff
  USE mesh ,     ONLY : ra,r,dr
  USE machin ,   ONLY : rmajor 
  USE numbrs ,   ONLY : nj
  USE soln,      ONLY : te,ene
  USE sourc,     ONLY : eta
  USE tcoef,     ONLY : chie_paleo,d,xke_paleo
  USE mesh,      ONLY : grho2_mesh


  IMPLICIT NONE
  INTEGER j,include_paleo  
  REAL *8 rmajloc,vthe,wpe,deltae,dqdrho,tea,enea,          &
          xlam,teaerg,qa,etapc,eta0,nmax,zeffa,             &
          paleo_mult,etaa,paleo_pinch_mult,M,tem,           &
          tp1,tm1,tm0,hcapa,chipc,fac0,                     &
          r_paleo_pinch_cut,paleo_pinch_mult2,grho2a
  REAL *8, ALLOCATABLE,DIMENSION(:) :: lame,lmax,fac,xnue
  REAL *8, ALLOCATABLE,DIMENSION(:) :: xm


  CONTAINS


    SUBROUTINE paleocl_tr(wrt,nion)
! -------------------------------------------------------------------
!   Determined chieo_paleo the paleoclassical diffusion 
!   coefficient in units of cm**2/sec. Add it to d.
! --------------------------------------------------------------------

  INTEGER wrt,nion
  REAL *8 nmaxp,qdp,qdp2(nj),temp
  LOGICAL qdpset

  IF( .NOT. ALLOCATED(lame))THEN 
     ALLOCATE(lame(nj-1))
     ALLOCATE(lmax(nj-1))
     ALLOCATE(fac(nj-1))
     ALLOCATE(xnue(nj-1))
     ALLOCATE(xm(nj-1))
  ENDIF


  rmajloc = rmajor
  IF (rmajloc .LT. 10.) rmajloc = rmajloc*100.0  !make sure rmajloc is in cm
  qdpset = .FALSE.
  DO j = 1,nj-1                                     !zone center quantities
    IF(freeze_xte .EQ. 0)THEN
      dqdrho  = ABS(q(j+1)-q(j))/(r(j+1)-r(j))
      zeffa   = 0.5*(zeff(j+1)+zeff(j))
      tea     = 0.5*(te(j+1)+te(j))                   !kev
      enea    = 0.5*(ene(j+1)+ene(j))                 !#/cm**3
      qa      = 0.5*ABS(q(j+1)+q(j))
      vthe    = 1.325e9*SQRT(tea)                     !cm/sec
      wpe     = 5.64e4*SQRT(enea)                     !rad/sec
      grho2a  = 0.5*(grho2_mesh(j+1)+grho2_mesh(j))
      deltae  = cee/wpe                               !cm
      teaerg  = 1.6e-9 * tea                          !erg
      xlam    = 24.0 - LOG (SQRT (enea)/(1.0e3*tea))  !for e-e collisions
      xnue(j) = 1.33333*rot2pi*enea*zeffa*charg4*xlam  &
                       / SQRT (xmasse*teaerg**3)      !1/sec
      lame(j) = vthe/xnue(j)                          !cm
      nmax   = 1./SQRT(pi*deltae*dqdrho)
                !nmaxp correction added 03/22/05  as per Callen mail :
      IF(ra(j)/r(nj) .LT. r_paleo_pinch_cut) THEN
         !second derivative of q, non uniform grid:
         IF( .NOT. qdpset)CALL dif2ydx(dr, q, qdp2, nj)
         qdpset = .TRUE.
         !qdp = 0.5*(qdp2(j) + qdp2(j+1))
         !use the mag axis value for all j (Callen email) :
         qdp = qdp2(1)
         nmaxp   = 1. / ((pi*deltae)**2* ABS(qdp) )**0.3333
         nmax    = MIN(nmax,nmaxp)
      ENDIF
      fac0    = pi*rmajloc                            !cm
      fac(j)  = fac0*qa                               !cm,assumes Rbar = rmajor
      lmax(j) = fac(j)*nmax                           !cm
      eta0    = xmasse*xnue(j)*8.98755179e11/(enea*charge**2)   !ohm cm      
      etaa    = eta(j)*8.98755179e11                            !ohm cm 
      etapc   = 100.*etaa/u0                                    ! cm**2/sec
      temp    = fac(j)*(lmax(j)+lame(j))
      IF(temp .ne. 0.0)THEN
!         xm(j) = (lmax(j)+lame(j))/temp
          xm(j) = (lmax(j)*lame(j))/temp               !changed 04/09/07 HSJ 
          ! 08/29/2011 SPS
          ! Use grho2 (<|rho|^2>) normalization to match ONETWO definitions 
         chie_paleo(j) = 1.5*(xm(j)+1.)*etapc*grho2a   !cm^2/sec
      ELSE
          xm(j) =  0.0
          chie_paleo(j) = 0.0
      ENDIF
      !update xke_paleo only if freeze_xte = 0.
      !if freeze_xte =1 then use old value for xke_paleo
      xke_paleo(j) =  enea*chie_paleo(j)
    ENDIF
        !set d for electron energy equation to include paleoclassical transport
        d(nion+1,nion+1,j) = d(nion+1,nion+1,j)+ paleo_mult*xke_paleo(j)
  ENDDO


  RETURN
  END SUBROUTINE paleocl_tr
  

  SUBROUTINE paleocl_output(nout,jprt)
    INTEGER,INTENT(in)  :: nout,jprt
    INTEGER j,j1prt
     WRITE(nout,'(" --------------Paleoclassical results -----")')
     WRITE(nout,'("    j","     ra      ", "     te    ",               &
        &   "      ene   ","    q    ","   dqdr   ",                     &
        &   "   chiepc  ","   lame    ",                                 &
        &   "   lmax   ","    pi_R0_q ","     nue     ")')
     DO j=1,nj-1
        qa = ABS((q(j+1)+q(j)))*0.5
        tea = 0.5*(te(j+1)+te(j))
        enea = 0.5*(ene(j+1)+ene(j))
        dqdrho  = ABS(q(j+1)-q(j))/(r(j+1)-r(j))
        j1prt = ((j-1)/jprt)*jprt
        IF (j1prt .NE. j-1 .AND. j .NE. nj-1) go to 371
        WRITE(nout,1) j,ra(j),tea,enea,qa,dqdrho,chie_paleo(j),            &
                    lame(j),lmax(j),fac(j),xnue(j)
371     CONTINUE
     ENDDO
1    FORMAT(i4,'.5',10(x,1pe10.2))
  END SUBROUTINE paleocl_output

  END MODULE paleocl
