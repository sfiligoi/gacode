


      double precision function enorm(n,x)
      integer n
      double precision x(n)
c     **********
c
c     function enorm
c
c     given an n-vector x, this function calculates the
c     euclidean norm of x.
c
c     the euclidean norm is computed by accumulating the sum of
c     squares in three different sums. the sums of squares for the
c     small and large components are scaled so that no overflows
c     occur. non-destructive underflows are permitted. underflows
c     and overflows do not occur in the computation of the unscaled
c     sum of squares for the intermediate components.
c     the definitions of small, intermediate and large components
c     depend on two constants, rdwarf and rgiant. the main
c     restrictions on these constants are that rdwarf**2 not
c     underflow and rgiant**2 not overflow. the constants
c     given here are suitable for every known computer.
c
c     the function statement is
c
c       double precision function enorm(n,x)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an input array of length n.
c
c     subprograms called
c
c       fortran-supplied ... dabs,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i
      double precision agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs,
     *                 x1max,x3max,zero
      data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant/floatn
      do 90 i = 1, n
         xabs = dabs(x(i))
         if (xabs .gt. rdwarf .and. xabs .lt. agiant) go to 70
            if (xabs .le. rdwarf) go to 30
c
c              sum for large components.
c
               if (xabs .le. x1max) go to 10
                  s1 = one + s1*(x1max/xabs)**2
                  x1max = xabs
                  go to 20
   10          continue
                  s1 = s1 + (xabs/x1max)**2
   20          continue
               go to 60
   30       continue
c
c              sum for small components.
c
               if (xabs .le. x3max) go to 40
                  s3 = one + s3*(x3max/xabs)**2
                  x3max = xabs
                  go to 50
   40          continue
                  if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2
   50          continue
   60       continue
            go to 80
   70    continue
c
c           sum for intermediate components.
c
            s2 = s2 + xabs**2
   80    continue
   90    continue
c
c     calculation of norm.
c
      if (s1 .eq. zero) go to 100
         enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
         go to 130
  100 continue
         if (s2 .eq. zero) go to 110
            if (s2 .ge. x3max)
     *         enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
            if (s2 .lt. x3max)
     *         enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
            go to 120
  110    continue
            enorm = x3max*dsqrt(s3)
  120    continue
  130 continue
      return
c
c     last card of function enorm.
c
      end



      subroutine eval_non_lin_eq(neq,x,residd)
c-----------------------------------------------------HSJ-----
c----This routine supplies all the transport related information required
c    by the non-linear equation solver routines
c INPUT:
c     neq # equations and unknowns
c     x(1...neq) the value of the unknowns at which the transport 
c     finite difference equations are to be evaluated

c
c OUTPUT:
c     residd(neq)   the residual of each equation
c     iflag is normally not changed or used in this routine but
c     iflag can be set to a negative # in this routine  to indicate that 
c     the solution process should be terminated
c--------------------------------------------------------------------

      USE param
      USE aid_newton
      USE     solcon
      USE     soln
      USE mhdpar
      USE numbrs
      USE mesh
      USE flags
      USE tordlrot
      USE soln2d
      USE tcoef
      USE nonlin
      USE ifs
      USE staebler
      implicit  integer (i-n), real*8 (a-h, o-z)
      real *8 x(*),residd(*)
      integer j,k,nu

c      include 'param.i'    !kk,
c      include 'flags.i'    !itran(kk)
c      include 'mesh.i'
c      include 'numbrs.i'   !nk,nkt,nj
c      include 'soln.i'     !u(kk,kj),rbp,te,ti,en
c      include 'soln2d.i'   !xi_include
c      include 'solcon.i'   !cparam, note n is in common in solcon !!!
c      include 'tordlrot.i' !iangrot,angrot
c      include 'nonlin.i'   !xscale,resid,nrmax,iteration_method
c      include 'tcoef.i'    !temporary for debug only 
c      include 'mhdpar.i'    !temporary for debug only 
c      include 'ifs.i'      !temporary for debug only 
c      include 'staebler.i' !temporary for debug only 
      data icalln /0/

      

c     load u with current values
      nu = 0
      do j =1,nrmax
         do k=1,nk
            if(itran(k) .gt. 0)then
c              skip rbp at the origin (its value is fixed at 0.0)
               if(j .eq. 1 .and. k .eq. nk-iangrot)go to 10
               if(j .ge. te_index .and. k .eq. nk-iangrot-2)go to 10
               if(j .ge. ti_index .and. k .eq. nk-iangrot-1)go to 10
               if(j .ge. rot_index .and. k .eq. nk .and. 
     .                                          iangrot .eq. 1)go to 10
               if(j .ge. ni_index .and. k .eq. 1)go to 10 !jmp.den
               nu = nu+1
               if(k .le. nk-iangrot)then
chsj                u(k,j) = (xscale(nu)*x(nu))**2 !forces non-negative value
                    u(k,j) = x(nu)**2 !forces non-negative value
               else
chsj                u(k,j) = xscale(nu)*x(nu)
                    u(k,j) = x(nu)
               endif
 10         endif
         enddo
      enddo
c update sets en,te,ti,rbp,angrot to the values in vector u:
      call update (u, en, te, ti, rbp, nk, nj, kj, kk, iangrot, angrot)

c evaluate the terms required to form the equations:
c if iteration_method .eq. 1 then xchi is not changed for this
c solution cycle (note that default iteration_method = 0):
      if( iteration_method .eq. 0) call diffus (xi_include)
      if(freeze_source .eq. 0)  call source 
      call fluxx
c generate the equations by calling solve(solve in turn calls residual):
      call solve
c the required residual vector, resid(1...neq) is now known  in nonlin.i:
      call copya(resid,residd,nu) !resid ,residd required due to code structure
      icalln = icalln+1

      if(steady_state .lt. 1.e-5) 
     .    call savesol (u, usave, nk, nj, kk) !HSJ 8/14/04 copy u to usave,ene to enesav

 
      return
      end





      subroutine get_vars(maxscale,xvset,typx,
     .           nameu_gp,grid_point,x,xval,xval_old,nu)
c------------------------------------------------------------HSJ------------
c     Get info on variables at grid points that are to be calculated
c     For variables that have physical meaning only when positive
c     (or at least non-negative),
c     we return xval(i) = sqrt(v) where v is one of the dependent
c     variables at grid point i. Thus if the solution method
c     produces a negative xval(i) the actual value used in the
c     equations is xval(i)**2 .
c     This transformation is not done  for toroidal rotation since
c     it can be negative. (Negative poloidal B field is not allowed
c     so it is done for Faraday's law)
c    
c  INPUT (from argument list - there is additional  input from include files)
c   nrmax            #(in nonlin.i) of grid points for which a solution is sought.
c                    nrmax is no necessarily equal tp nj.
c   xval_old(1...nu) previous version of xval
c   xvset            logical, set to true if xval is to be averaged with
c                    xval_old
c   maxscale         logical, selects scaling option, use maxscale =flase
c  OUTPUT
c    nu              total number of variables for which a solution is sought
c    x(1..nk)        the scale factors for each dependent variable
c                    (includes ones that are not in simulation mode)
c    nameu_gp(1..nu) name of variables
c    grid_point(1..nu) grid point # associated with nameu_gp)
c    typx(1..nu)    typical size of variable
c    xval(1...nu)   
c---------------------------------------------------------------------------

      USE param
      USE     solcon  !get n (number of time step) from here
      USE     soln
      USE     ions
      USE numbrs
      USE mesh
      USE flags
      USE tordlrot
      USE verbose,only : newtonvb
      USE nonlin
      implicit  integer (i-n), real*8 (a-h, o-z)


      real *8 x(*),typx(*),xval(*),xval_old(*)
      integer *4 grid_point(*)
      integer *4 nuu,nu
      logical maxscale,xvset
      character *(*) nameu_gp(*)
c     find the scale factor for the physical variables.
      call maxu (u, itran, nk, nj, kk, x, iangrot)
      
      nu = 0
      do j =1,nrmax
         nuu =0
         do k=1,nk

            if(itran(k) .gt. 0)then
               nuu = nuu +1
c              skip rbp at the origin (its value is fixed at 0.0):
               if(j .eq. 1 .and. k .eq. nk-iangrot .and.
     .                              diffeq_methd .ne. 1)go to 10
c              skip boundary values:
               if(j .ge. te_index  .and. k .eq. nk-iangrot-2)go to 10
               if(j .ge. ti_index .and. k .eq. nk-iangrot-1)go to 10
               if(j .ge. rot_index .and. k .eq. nk .and. 
     .                                          iangrot .eq. 1)go to 10
               if(j .ge. ni_index .and. k .eq. 1)go to 10 !jmp.den
               
               nu = nu+1
               nameu_gp(nu) = nameu(k)
               grid_point(nu) = j
               if(k .le. nk-iangrot)then !eni,te,ti,rbp, must be non negative
                  xscale(nu) = SQRT(x(nuu))
                 if(k .eq. nk-iangrot .and. j .eq. 1)
     .            xscale(nu) = 1.0d0
chsj             xval(nu) = sqrt(u(k,j))/xscale(nu)
                 xval(nu) = sqrt(u(k,j))
                 if(maxscale)then
                      typx(nu) = xscale(nu)
                 else
                      typx(nu) = xval(nu) 
                 if(k .eq. nk-iangrot .and. j .eq. 1)
     .                typx(nu) = 1.0d0
                 endif
               else                      ! toroidal rotation can be negative
                 xscale(nu) = x(nuu)   
chsj             xval(nu) = u(k,j)/xscale(nu)
                 xval(nu) = u(k,j)
                 if(maxscale )then 
                    typx(nu) = xscale(nu)
                 else
                     typx(nu) = ABS(xval(nu))    !this assumes that initial input is
                 endif                       !about right for the scale of the solution
               endif
               if(xvset)then !we cut the previous time step in half
                   xval(nu)= 0.5*(xval(nu)+xval_old(nu))
               endif
 10         endif
         enddo
      enddo


      if(n .lt. 2 .and. newtonvb .ne.0)then ! print variable map info (variable n is in solcon.i ,
                       ! which is very baaad), n = # time steps taken so far
          do j=1,nu
                   print *,'j,gridpt,name',j,grid_point(j),nameu_gp(j)
                   print *,'j,typx,xscale,xval ',j,typx(j),xscale(j),
     .                                                         xval(j)
           enddo
           print *,'x =',x(1:nk)
      endif


      if(nu .le. 0)
     .   call STOP ('sub get_vars: no. equations to solve = 0', 0)

      return
      end



      subroutine nl_driver_12(n,x0,FVEC,global,analjac,cheapF,
     .                 bw_low,bw_up,sparse_jac,Jc,
     .                 factsec,typF,Fdigits,fvectol,steptol,
     .                 maxstep,itnlimit,delta,printcode,get_typf,
     .                 solmethod,noscale,xf,termcode,gradtol,
     .                 jmaxresid,jmaxgrad,
     .                 frozen,freeze_type,ssqrmin,
     .                 zero_step,first_step)
c--------------------------------------------------------------------------
c INPUT
c n          # unknowns
c x0(1...n)  initial guess at unknowns
c FVEC       External subroutine name for evaluating F(1...n)
c            and sum of squares of residulas
c global  =   1 ! line search 
c global  =   2 ! Hookstep trust region
c global  =   3 Dogleg trust region 
c analjac    logical true  if anlaytic Jacobian is to be used
c            analjac = true results in a call to ???
c cheapf     used only if analjac = .flase.
c              set cheapf = .true.  if f is cheap to evaluate
c             (uses finite differences if cheapf = .true.,
c              uses secant updates if cheapf = .false.)
c            if /  uses secants otherwise
c factsec    used only if analjac = .false. and cheapf = .false.
c            if factsec = .true. use factored secant approximation
c            to update a QR factorization of the Jacobian.
c            factsec = .false. is not implemented !
c typx(1...n)  gives typical magnitude of unknown i
c noscale 
c typx (from ai_newton)
c typF
c Fdigits
c fvectol
c gradtol    if gradient of 0.5*ssq resid is les than this number then
c            the problem is considered converged
c steptol,
c maxstep
c itnlimit    max iterations
c printcode   Fortranunit io unit  for error diagnostics
c get_typf    if 0 means do one call to residual to get the scale of
c             the equations, then exit.
c             if 1 then typf is assumed set and we go for the solution.
c nameu_gp(j)  from aid_newton,character *8 array gives name of variable j
c grid_point(j) gives rho grid point associated with index j (aid_newton)
c
c ssqrmin       minimum sum of squares of residuals to be achieved
c               (if possible)
c
c jac_known    logical value set to .true. if jacobian at current point is 
c              avaialble in JC   . If jacobian at this point is not
c              known then set jac_known = .false. and it will be calculated.
c              (from  aid_newton)
c jac_skip     in aid_newton the jacobian will be calculated fresh each time
c              Mod(iters,jac_skip) = 0. 
c              for time dependent problems if the time step is small
c              enough it may be possible to save execution time by 
c              setting jac_skip to a number greater than 1.
c              most likely jac_skip =1,2,3, will be most useful but
c              you can try larger numbers. note that jac_skip =1 means 
c              calculate the jacobian on every iteration.
c    
c
c OUTPUT:
c
c xf(i)      i=1,..n solution vector
c termcode   solution indicator
c ssqr       sum of squares of residuals
c itncount    # of iterations done 

c------------------------------------------------------------HSJ-----------
      USE aid_newton
      USE allocate_err
      USE COM  
      USE copy
      USE terminate
      USE verbose
      implicit none
c      include 'verbose.i'
!INPUT:
      INTEGER 
     .       n,global,Fdigits,itnlimit,printcode,frozen,freeze_type,
     .       bw_low,bw_up,get_typf,jmaxresid,jmaxgrad,
     .       printcode1,first_step
      REAL *8
     .       x0(n),typf(n),Jc(n,n),
     .       fvectol,steptol,maxstep,delta,
     .       macheps,eta,gradtol,ssqrmin
      LOGICAL
     .       analjac,cheapf,factsec,noscale,sparse_jac,
     .       zero_step
      CHARACTER *(*) solmethod
      EXTERNAL FVEC
c     EXTERNAL JAC           !name of analytic jacobian subroutine
c LOCAL:
       REAL *8 fc,fp, FVc(n),gc(n),M2(n),M(n,n),
     .       sumd,delta_prev,mu,phi,phip,bandw,smallestv,
     .       largestv
       REAL *8 ,dimension(:)  ,allocatable :: xc,Sn,xp,Hc
       INTEGER i,j,retcode,istat,consecmax,nn,kbw,ibw,
     .         ibwmax,jbwmax
       LOGICAL 
     .       restart,maxtaken,testing
!OUTPUT:
      INTEGER termcode
      REAL*8
     .     xf(n)
 


!may be allocated with different size from previous call 
      IF( allocated(Sf) )                      
     .  deallocate (Sf, STAT = istat)
        allocate (Sf(n),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("Sf, nedriver",0,istat)

      IF( allocated(Sx) )                      
     .  deallocate (Sx, STAT = istat)
        allocate (Sx(n),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("Sx, nedriver",0,istat)

      IF(allocated(FVp))
     .  deallocate (FVp,STAT = istat)
        allocate (FVp(n),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("FVp, nedriver",0,istat)




      IF(allocated(xc))
     .  deallocate (xc,STAT = istat)
        allocate (xc(n),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("xc, nedriver",0,istat)

      IF( allocated(xp))
     .  deallocate (xp,STAT = istat)
        allocate (xp(n),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("xp, nedriver",0,istat)

      IF( allocated(Sn))
     .  deallocate (Sn,STAT = istat)
        allocate (Sn(n),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("Sn, nedriver",0,istat)

      nn = (n*(n+1))/2
      IF( allocated(Hc))
     .  deallocate (Hc,STAT = istat)
        allocate (Hc(nn),STAT = istat)
        IF(istat .ne. 0)
     .          call allocate_error("Hc, nedriver",0,istat)
c         macheps = 1.d-14   !test HSJ 11/18/02
c        if(ABS(macheps) .lt. 1.d-25)
c     . CALL MACHINE_EPSD(macheps)
       CALL MACHINE_EPSD(macheps)   !hsj 7/20/02
c     set Sx = 1./typx, sf =1./typf,eta,dleta,termcode = 0 or -1:
      CALL NEINCK(n,macheps,x0,typx,typf,Fdigits,fvectol,
     .                   steptol,maxstep,itnlimit,printcode,
     .                   delta, global,analjac,cheapf,factsec,
     .                   Sx,Sf,eta,termcode)

      IF(termcode .lt. 0)then
         call copya(x0,xf,n)      !copy x0 to xf
         write(printcode,'("error in input to nedriver",/, 
     .   "termcode = ",i5     )')termcode
         CALL STOP('subroutine NEDRIVER: input error', 6)
      ENDIF 
      itncount = 0
c     get 0.5* sum of squares of residuals of nonlinear equations:
      CALL nefn(n,x0,fc,FVEC,NOSCALE)  !fc = 0.5*Norm(Sf*Fvp)**2,FVp loaded in common

c      CALL nestop0(n,FVp,Sf,fvectol,termcode,consecmax)
       termcode = 0 !dont use nestop0 to set termcode HSJ 4/30/02
      !if we are just trying to discover the scale of the equations
      !then get_typf = 0:
      if(get_typf .eq. 0)termcode = 1 !force exit below
      IF(termcode .gt. 0)then
         CALL copya(x0,xf,n)          !nestop0 returned termcode = 1
      ELSE  !get inital jacobian      !copy x0 to xf and exit below
         if(.not. jac_known)then
             If(analjac) then
c                CALL JAC(n,x0,Jc)
             ELSE IF(sparse_jac)then
                 CALL fdjac_sparse(n,x0,FVp,FVEC,Jc,bw_low,bw_up,eta,Sx)
             ELSE
                 CALL fdjac(n,x0,FVp,FVEC,Sx,eta,Jc)
             ENDIF
             DO i =1,n
                sumd = 0.0d0
                DO j=1,n
                   sumd =sumd+Jc(j,i)*FVp(j)*Sf(j)**2
                ENDDO
                gc(i) = sumd
             ENDDO
             !bandwidth not otherwise used -----------------------------
             testing = .false.
             if(testing)then
             kbw =0
             DO  j=1,n
                 bandw = 0.0
                 smallestv  =0.0
                 largestv =0.0
                 ibw = 0
                 DO i =1,n
                     bandw = bandw +ABS(jc(j,i)) ! *Sf(j)) 
                     smallestv = min(smallestv,ABS(Jc(j,i))) !*Sf(j)))
                     largestv = max(largestv,ABS(Jc(j,i))) !*Sf(j)))
                     if(ABS(Jc(j,i)) .gt. 1.e-10)ibw =ibw+1
                     if(ABS(Jc(j,i)) .gt. 1.e-10 .and. j .eq. 21)
     .               print *,'jc(21,i),i =',jc(21,i),i
                 ENDDO
                 print *,'# non zero ,j =',ibw,j
                 bandw = bandw/n   
                 ibwmax = MAX(ibw,ibwmax)
                 if(ibw .eq. ibwmax)jbwmax=j  !equation j has max elements ne 0
                 ibw = 0
                 do i=1,n                 
                      if(ABS(Jc(j,i)) .gt. bandw)ibw =ibw+1
                 enddo !bandwidth not otherwise used
                 print *,'ibw,j =',ibw,j 
                 print *,'smalst,largest abs value=',smallestv,largestv
                 kbw = MAX(kbw,ibw)
             ENDDO 
             print *,'kbw =',kbw  
             print *,'equation # ,jac(j,i) :',jbwmax
             do i=1,n
                print *,jc(jbwmax,i),nameu_gp(i)
             enddo
             !end bandwidth calcs ----------------
             ENDIF
          else
             !jacobian is known from previous solution point.
             print *,'nl_driver jac_known =',jac_known
          endif
             CALL copya(FVp,FVc,n)
             CALL copya(x0,xc,n)
      ENDIF
      restart =.true.  !used only for factored broyden method

c     iteration section:
c     repeatedly solve the linearized problem using 
c     (globally convergent) inexact Newton methods
c     until max iterations are used up or problem is converged
      first_step = 0          !flags.i used if dv method selected,and ddebug(50)=1
      DO WHILE (termcode .eq. 0)
         itncount =itncount+1
            !First calculate the Newton step,Sn, using nemodel or nemodelfac
            !Also factor
            !the Jacobian into M and Hc(the new solution point would 
            !be xc+Sn. Which is formed in linesearch,hookstep, or dogleg.
         IF(analjac .or. cheapf .or. .NOT. factsec)THEN

            CALL nemodel(n,FVc,Jc,gc,Sf,Sx,macheps,global,
     .                   M,Hc,Sn)

         ELSE
            CALL nemodelfac(n,FVc,gc,Sf,Sx,macheps,global,
     .                      restart,M,M2,Jc,Hc,Sn) 

         ENDIF
         IF(global .eq. 1)THEN
c            print *,'entering linesearch'
c             print *,'maxstep =',maxstep
            CALL linesearch(n,xc,fc,FVEC,gc,Sn,Sx,noscale,Sf,maxstep,
     .                      steptol,retcode,xp,fp,FVp,maxtaken)
c           linesearch returns fp(= SSQR) as determined in FVEC
c           with no additional scaling.
c            call STOP('subroutine line seacrh 616',1)

         ELSEIF(global .eq. 2)THEN
            CALL hookdriver(n,xc,fc,FVEC,gc,M,Hc,Sn,Sx,noscale,Sf,
     .                    FVp,maxstep,steptol,itncount,macheps,delta,
     .                    mu,delta_prev,phi,phip,retcode,xp,fp,
     .                    maxtaken,printcode)
c           hookdriver (through trust_region) returns fp(= SSQR) 
c           as determined in FVEC
c           with no additional scaling.

         ELSE   !global = 3
            CALL drive_dogleg(n,xc,fc,FVEC,gc,M,Sn,Sx,noscale,Sf,
     .                        FVp,maxstep,steptol,delta,
     .                        retcode, xp,fp,maxtaken,
     .                        printcode)
c           drive_dogleg(through trust_region) returns fp(= SSQR) 
c           as determined in FVEC
c           with no additional scaling.

         ENDIF

c         print *,'retcode ,rstart,analjac,cheapftermcode=',
c     .   retcode,restart,analjac,cheapf,termcode


         !retcode = 1 means satisfactory new step not found
         !otherwise retcode = 0
 10      IF(retcode .eq. 1 .and. (.not. restart) .and. (.not. analjac)
     .         .and. (.not. cheapf))THEN 
             !this section is for factsec = true cases 
             !( .not. analjac and .not. cheapf ==> factsec = true)
             IF(sparse_jac)then
                CALL fdjac_sparse(n,xc,FVc,FVEC,Jc,bw_low,
     .                            bw_up,eta,sx)
             ELSE
                CALL fdjac(n,xc,FVc,FVEC,Sx,eta,Jc)
             ENDIF
            DO i =1,n
               sumd = 0.0d0
               DO j=1,n
                  sumd =sumd+Jc(j,i)*FVc(j)*Sf(j)**2
               ENDDO
               gc(i) = sumd
            ENDDO
            IF(global .eq. 2 .or. global .eq. 3)delta = -1
            restart = .true.
            termcode = 0    !HSJ added
c            jac_known = .true.       !set for next time subroutine nl_driver_12 is called
         ELSE   !complete the iteration by evaluating the Jacobian at 
                !the new point xp

            !use old Jacobian if called for
            if(MOD(itncount,jac_skip).eq. 0)then
                IF(analjac)then
c                   CALL JAC(n,xp,Jc)
                ELSEIF(cheapf)then
                   IF(sparse_jac)then
                       CALL fdjac_sparse(n,xp,FVp,FVEC,Jc,bw_low,
     .                            bw_up,eta,sx)
                   ELSE
                       CALL fdjac(n,xp,FVp,FVEC,Sx,eta,Jc)
                   ENDIF
                ELSEIF(factsec)then
                   print *,'calling broyfac'
                   CALL broyfac(n,xc,xp,FVc,Fvp,eta,Sx,Sf,Jc,
     .                          M,M2)
                ELSE !factsec =.false.
                   CALL broyunfac(n,xc,xp,FVc,Fvp,eta,Sx,Jc)
                ENDIF
             endif
  

            IF(factsec)THEN      !get gc ( gradient) from QR factorization
               !gc = Jc*Sf*FVp     !note that Jc contains Q transpose here
               Do i =1,n  
                  sumd =0.0d0
                  DO j=1,n
                     sumd =sumd+Jc(i,j)*FVp(j)*Sf(j)
                  ENDDO
                  gc(i) = sumd
               ENDDO
               Do i =n,1,-1
                  sumd =0.0d0
                  DO j=1,i
                     sumd =sumd+M(j,i)*gc(j)
                  ENDDO
                  gc(i) = sumd
               ENDDO
            ELSE                 !get gc ( gradient) from definition
                                 !in scaled form
               DO i =1,n
                  sumd = 0.0d0
               DO j=1,n
                  sumd =sumd+Jc(j,i)*FVp(j)*Sf(j)**2
               ENDDO
                  gc(i) = sumd
               ENDDO

            ENDIF
            zero_step = .false.
            printcode1 = printcode
            if(newtonvb .eq. 0)printcode1 = 0
            CALL check_convergence(n,xc,xp,FVp,gc,fp,Sx,Sf,retcode,
     .            fvectol,gradtol,steptol,noscale,itncount,itnlimit,
     .            maxtaken,consecmax,termcode,printcode1,gradmax,ssqr,
     .            solmethod,jmaxresid,jmaxgrad,nameu_gp,grid_point,
     .            frozen,freeze_type,ssqrmin)
            print *,'termcode after  check_conver',termcode
c           if step failed set zero_step
            if(termcode .eq. 2 .or. termcode .eq. 3 )
     .                                       zero_step = .true.
            IF(termcode .eq. 2 .and. (.not. restart) .and. 
     .                (.not. analjac) .and. (.not. cheapf))THEN
               !restart
               retcode = 1

               go to 10
            ELSEIF (termcode .gt. 0)THEN
               call copya(xp,xf,n)
            ELSE
               restart = .FALSE.
            ENDIF
            CALL copya(xp,xc,n)
            fc = fp
            FVc(:) = FVp(:)
         ENDIF
c         jac_known = .true.  !set for next time subroutine nl_driver_12 is called
      ENDDO  ! termcode = 0 DO WHILE 




c     deallocation required, see comments  in driver routine
c      deallocate(Sf,STAT = istat)
c      deallocate(Sx,STAT = istat)
      deallocate(Fvp,STAT = istat)
      deallocate(xc,STAT = istat)
      deallocate(xp,STAT = istat)
      deallocate(Sn,STAT = istat)
      deallocate(Hc,STAT = istat)



      RETURN
      END



      subroutine output_test(nu,xval,restart,termcode,
     .                       ipturb,resid,frozen,ind_test,
     .                        nameu_gp,typf,grid_point,ssqrmin,
     .                        NOSCALE,iendit)
c-----------------------------------------------------------
          USE aid_newton, only :  typx
          USE param
          USE mhdpar
          USE com
          USE numbrs
          USE mesh
          USE machin
          USE tcoef
          USE ifs
          implicit  integer (i-n), real*8 (a-h, o-z)
c          include 'param.i'
c          include 'mhdpar.i' !kpsi for use in ifs.i
c          include 'ifs.i'    !xke_ifs,xki_ifs,xkang_ifs (used for glf23 storage)
c          include 'numbrs.i' ! nion
c          include 'machin.i' !rmajor
c          include 'mesh.i'   !te,ti,rot_index
c          include 'tcoef.i'  !d(i,j,k),chii_ms,chie_ms,etaphi_ms
          integer restart,termcode,frozen,
     .    grid_point(nu),ientry
          real *8 xval(nu),resid(nu),typf(nu),ssqrmin,
     .    dumy,sum,dumy2,dumy3
          LOGICAL  exists, opened, NOSCALE
          character *3 status
          character *(*),nameu_gp(nu)
          data ientry / 0 /
          iounit = 92
          inquire (file = 'nwt.test', iostat = iostat,
     .              exist = exists, opened = opened)
          if (opened) then   ! file is open
             !do nothing
          else
             if (exists )          ! old file exists, trash it
     .               call DESTROY ('nwt.test')
             !open the file
             open (unit = iounit, file = 'nwt.test',
     .            status = 'NEW', iostat = iostat)
          endif
           ientry = ientry+1


c          the values in vector resid are stale, update them here:
           call eval_non_lin_eq(nu,xval,resid)

           sum =0.0d0
          do j =1,nu
              sum = sum +(resid(j)/typf(j))**2
          enddo
          write(iounit,1)ientry,grid_point(nu),restart,termcode,ipturb,
     .                   frozen,ind_test,sum,te_index,ti_index,
     .                   rot_index,nu
 1         format(7(2x,i8),2x,1pe14.8,4(2x,i4))

           do j =1,nu
              dumy = xval(j)
              k = grid_point(j)
              if(nameu_gp(j)(1:2) .ne. 'wr')dumy =dumy*dumy
              if(nameu_gp(j)(1:2) .eq. 'te') then
                   dumy2 = d(nion+1,nion+1,k)
c                   dumy3 = chie_ms(k)*1.e4
                   dumy3 = xke_ifs(k)
              endif
              if(nameu_gp(j)(1:2) .eq. 'ti')then 
                   dumy2 = d(nion+2,nion+2,k)
c                   dumy3 = chii_ms(k)*1.e4
                   dumy3 = xki_ifs(k)
              endif
              if(nameu_gp(j)(1:2) .eq. 'bp')then 
                   dumy2 = d(nion+3,nion+3,k)
                   dumy3 = 0.0
              endif
              if(nameu_gp(j)(1:2) .eq. 'wr')then 
                   dumy2 = d(nion+4,nion+4,k)
                   dumy3 = xkang_ifs(k)*dmassden(k)*rmajor*rmajor
              endif
              if(nameu_gp(j)(1:2) .eq. 'ni')then  !'ni' is not a vaild name 
                   dumy2 = d(1,1,k)
                   dumy3 = 0.0
              endif
              write(iounit,2)dumy,resid(j)/typf(j),dumy2,dumy3,
     .        j,k,nameu_gp(j)(1:8)
 2            format(4(2x,1pe16.8),2x,i4,2x,i3,2x,a)
           enddo
           if(iendit .eq. 1)then
              write(iounit,FMT ='("the end")')
              write(iounit,3)(typx(j),j=1,nu)
              write(iounit,3)(typf(j),j=1,nu)
              write(iounit,3)(xval(j),j=1,nu)
 3            format(5(2x,1pe14.6))
           endif


c           print *,'iendit in out =',iendit
           return
           end





      subroutine residual(jrgrid,a,b,c,g)
c---------------------------------------------------HSJ---------------
c      This subroutine calculates the residual for the finite differenced
c      form of the transport equations. This residual can be used in
c      various non linear solver schemes.
c      Each radial grid point, designated by jrgrid,introduces a set of
c      nkt unknowns where nkt is the number of itran(k) values that
c      are non zero (eg the number of dependent variables run in
c      simulation mode).
c      OUTPUT
c      nequations    current equation number (the total number of equations
c                    and the total number of unknowns must be equal of course).
c      resid         residual of equation number neq,in nonlin.i
c      typf          typical magnitude of  of equation Nequation
c                    returned if get_typf =0 (in nonlin.i)
c---------------------------------------------------------------------
c
c
      USE param
      USE     soln
      USE numbrs
      USE mesh
      USE flags
      USE tordlrot
      USE nonlin
      implicit  integer (i-n), real*8 (a-h, o-z)
      integer l,k1,k2,jrgrid 
c      include 'param.i'    !kk,
c      include 'flags.i'    !itran(kk)
c      include 'numbrs.i'   !nk,nkt,nj
c      include 'mesh.i'     !te,ti,rot_index
c      include 'nonlin.i'   !resid,typf
c      include 'soln.i'     !u(kk,kj)
c      include 'tordlrot.i' !iangrot
      real *8 sum, a(kk,kk),b(kk,kk),c(kk,kk),g(kk)
      !temp HSJ
      logical testing
      REAL*8 aa,bb,cc,xx,ss,gg,terma,termb,termc
      common /temporary / aa,bb,cc,xx(51),ss,gg,terma,termb,termc
      !end temp HSJ

      if(jrgrid .eq. nj)return
      
      k1=0
      do k = 1,nk
         if( itran(k) .gt. 0)then
            k1 =k1 +1
            scale_eqn = 1.e0
            sum = 0.0d0                            ! nkt is the number of simulation variables
            if(jrgrid .eq. 1 .and. k .eq. nk-iangrot)
     .                           go to 30
            if(jrgrid .ge. te_index .and. k .eq. 
     .                         nk-2-iangrot)go to 30
            if(jrgrid .ge. ti_index  .and. k .eq. 
     .                         nk-1-iangrot)go to 30
            if(jrgrid .ge. rot_index .and. k .eq. nk .and. 
     .                                   iangrot .eq. 1)go to 30
            if(jrgrid .ge. ni_index .and. k .eq. 1)go to 30 !jmp.den
chsj            if(k .eq. nk-iangrot)scale_eqn = 1.e-2
chsj            if(k .eq. nk-2-iangrot .or. k .eq. nk-1-iangrot)
chsj     .                                        scale_eqn = 1.e-15
            nequations = nequations+1     
            k2 =0
            do l =1,nk
               if(itran(l) .gt. 0)then
                  k2 = k2+1 
                  if(jrgrid .eq. 1)go to 10
                  sum =sum +a(k1,k2)*u(l,jrgrid-1)
                  terma = a(k1,k2)*u(l,jrgrid-1)
c                  write(*,'("terma = ",1pe12.4)')terma
 10               sum =sum + b(k1,k2)*u(l,jrgrid)
                  termb = b(k1,k2)*u(l,jrgrid)
c                  write(*,'("termb = ",1pe12.4)')termb
                  if(jrgrid .eq. nj) go to 20
                  sum = sum + c(k1,k2)*u(l,jrgrid+1)
                  termc = c(k1,k2)*u(l,jrgrid+1)
c                  write(*,'("termc = ",1pe12.4)')termc
c                  write(*,'("jrgrid,k1,k2 =",3(1x,i5))')jrgrid,k1,k2
c                  print *,"k1,k2,l,jrgrid =",k1,k2,l,jrgrid
c                  if( jrgrid .gt. 1 .and. jrgrid .lt. nj)then
c                     print *,"u(l,jrgrid-1,jrgrid,jrgrid+1) =",
c     .               u(l,jrgrid-1),u(l,jrgrid),u(l,jrgrid+1)
c                  else if(jrgrid .eq. 1) then
c                     print *,"u(l,jrgrid,jrgrid+1) =",
c     .               u(l,jrgrid),u(l,jrgrid+1)
c                  else
c                     print *,"u(l,jrgrid-1,jrgrid) =",
c     .               u(l,jrgrid-1),u(l,jrgrid)
c                  endif



 20            endif
            enddo
chsj        resid(nequations)  = scale_eqn*(sum -g(k1))
            
            resid(nequations)  = (sum -g(k1))
c            write(*,'("sum,g(k1)  ",2(1pe25.14,x))')sum,g(k1)
 
!            if(get_typf .eq. 0)typf(nequations) = ABS(resid(nequations)) !hsj 10/3/05
!            if(get_typf .eq. 0)typf(nequations) = ABS(sum)  !previously used this
             if(get_typf .eq. 0)typf(nequations) = ABS(.01*sum)
c            print *,'get_typf =',get_typf
c            print *,'nequations,resid(nequations) =',nequations,
c     .                                  resid(nequations)
c            print *,"sum,g(k1)=",sum,g(k1)
 30      endif
         testing = .false.  !set in sub fdjac  as well
         IF(testing)then
             do l=1,nk
              if(itran(l) .gt. 0)then
              do j=1,51   !note 51 
              xx(j)  = u(l,j)
              enddo
              endif
             enddo       !temp for debug
          ENDIF
       enddo
       return
       end

      subroutine set_freeze (freeze_type,itte,itti,itangrot,itran,
     .                freeze_xte,freeze_xti,freeze_xwr,
     .                freeze_xrbp,freeze_xni,itran_save,
     .                itti_save,itte_save,itangrot_save,
     .                init,frozen,nk,iangrot,jac_known)
c-------------------------------------------------------------------


      use flags, only : iten !jmp.den
      implicit none

      integer freeze_type,itte,itti,nk,itangrot,itran(nk),
     .                freeze_xte,freeze_xti,freeze_xwr,
     .                freeze_xrbp,freeze_xni,itran_save(nk),
     .                itti_save,itte_save,itangrot_save,init,
     .                frozen,iangrot
       logical jac_known

      if(init .eq. 1)then   !make sure problem is un-frozen start of 
         frozen = 0         !each time step
         freeze_xte = 0
         freeze_xti =0
         freeze_xrbp =0
         freeze_xwr =0
         freeze_xni = 0
         itran_save(:) = itran(:)
         itte = itran(nk-iangrot-2)
         itti = itran(nk-iangrot-1)
c         itxj = itran(nk-iangrot)
         if(iangrot .gt. 0)itangrot = itran(nk)
         itti_save = itti
         itte_save = itte
         itangrot_save = itangrot
         return
       endif

      if(init .eq. -1)then
         itran(:) = itran_save(:)
         itte = itte_save
         itti = itti_save
c         itxj = itxj_save
c         itangrot_save = itangrot
         itangrot = itangrot_save !hsj 10/19/05
         return
       endif

c      Freeze something according to freeze_type
c      -5 to -1 freezes the profiles so they drop out of the
c      problem
c      1 to 5 freezes the associated anomalous diffusion coefficient
c      the corresponding profiles are still part of the unknown
c      set we are solving for. 

           if(freeze_type .eq. -5 .and. itte*itti*itangrot .eq. 1)then
               itran(nk-iangrot-2) =  0
               itte = 0
               itran(nk-iangrot-1) = 0
               itti = 0
               itran(nk) = 0
               itangrot = 0
               frozen = 1
           else if(freeze_type .eq. -4 .and. itte*itti .eq. 1)then 
               itran(nk-iangrot-2) =  0
               itte = 0
               itran(nk-iangrot-1) = 0
               itti = 0
               frozen = 1
           else if(freeze_type .eq. -3 .and. itte .eq. 1)then 
               itran(nk-iangrot-2) =  0
               itte = 0
               frozen = 1
           else if(freeze_type .eq. -2 .and. itti .eq. 1)then 
               itran(nk-iangrot-1) =  0
               itti = 0
               frozen = 1
           else if(freeze_type .eq. -1 .and. itangrot .eq. 1)then 
               itran(nk) =  0
               itangrot = 0
               frozen = 1
c          else if(freeze_type .eq. 0  .and. itangrot .eq. 1)then 
c                print *,'no freezing  method selected'
           else if(freeze_type .eq. 1 .and. itangrot .eq. 1)then 
               freeze_xwr = 1
               frozen = 1
           else if(freeze_type .eq. 2 .and. itti .eq. 1)then 
               freeze_xti = 1
               frozen = 1
           else if(freeze_type .eq. 3 .and. itte .eq. 1)then 
               freeze_xte = 1
               frozen = 1
           else if(freeze_type .eq. 4 .and. itte*itti .eq. 1)then 
               freeze_xte = 1
               freeze_xti = 1
               frozen = 1
           else if(freeze_type .eq. 5 .and. itte*itti*itangrot
     .                                                 .eq. 1)then 
               freeze_xte = 1
               freeze_xti = 1
               freeze_xwr = 1
               frozen = 1
           else if(freeze_type .eq. 6) then !jmp.den
               freeze_xte = 1               !jmp.den
               freeze_xti = 1               !jmp.den
               freeze_xwr = 1               !jmp.den
               freeze_xni = 1               !jmp.den
               frozen = 1                   !jmp.den
           endif
           
           if ((iten.eq.1).and.(frozen.eq.1)) then !jmp.den
               freeze_xni = 1               !jmp.den
           endif
           
           jac_known = .false.
           print *,'freeze_type,frozen,itte ,itti,itangrot,iangrot ='
           print *,freeze_type,frozen,itte,itti,itangrot,iangrot
           print *,'init , itte_save,itti_save,itran ='
           print *,init,itte_save,itti_save,itran
           print *,'itran_save =',itran_save
       return
       end



      subroutine solve_newton(info_conv)
c-------------------------------------------------------------HSJ-------
c     initialize the non-linear solver routines
c     This is the driver for the newton based methods
c     use allocatable arrays so that we dont grab core if this solution
c     method is not being used.
c-----------------------------------------------------------------------

      USE     param
      USE     aid_newton
      USE     ions
      USE     io
      USE     solcon
      USE     soln
      USE     COM
      USE numbrs
      USE mesh
      USE verbose
      USE flags
      USE tordlrot
      USE soln2d
      USE tcoef
      USE nonlin
      implicit  integer (i-n), real*8 (a-h, o-z)

      integer termcode,global, printcode,term3,xval_global_set,
     .        eff_bandwidth,jmaxresid,jmaxgrad,ind_test,
     .        ipturb,ipturbmax,frozen,restart,iostat,
     .        init_nl,imonit,imonit_old,ictr,itry1000,bandwidth_save  
      real *8  maxstep,jacmax,dumy,x(kk),xmonit(2)


       integer *4 ,dimension(:),allocatable,save ::itran_global
      real *8 , dimension(:),allocatable,save::  change,xval_global
      real *8 , dimension(:),allocatable,save::xscale_global,typf_global
      real *8 , dimension(:),allocatable,save::typx_global
      real *8 ,dimension(:,:),allocatable,save :: Jc
      real *8 ,dimension(:),allocatable,save ::xval,xval_in
      real *8 ,dimension(:),allocatable,save :: xval_old
      integer *4 ,dimension(:),allocatable,save ::itran_save
      real *4 ,dimension(:),allocatable ::xcRandom
      character *16,solmethod
      LOGICAL  noscale,analjac,cheapf,factsec,sparse_jac,
     .         rescale,maxscale,zero_step,uniform_weight

      external eval_non_lin_eq       


      if( .not. allocated(nameu_gp))then
        allocate (nameu_gp(nkt*nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("nameu_gp, solve_non_lin",0,istat)
        nameu_gp(:) = ' '
      endif

      if( .not. allocated(itran_save))then
        allocate (itran_save(nk),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("itran_save, solve_non_lin",0,istat)
      endif
      if( .not. allocated(itran_global))then
        allocate (itran_global(kk),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("itran_global, solve_non_lin",0,istat)
      endif
      if( .not. allocated(grid_point))then
        allocate (grid_point(nkt*nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("grid_point, solve_non_lin",0,istat)
        grid_point(:) = 0
      endif
      if( .not. allocated(xval))then
        allocate (xval(nkt*nj),STAT = istat)
        if(istat .ne. 0)
     .          call allocate_error("xval, solve_non_lin",0,istat)
      endif
      if( .not. allocated(xval_in))then
        allocate (xval_in(nkt*nj),STAT = istat)
        if(istat .ne. 0)
     .          call allocate_error("xval_in, solve_non_lin",0,istat)
      endif
      if( .not. allocated(xval_old))then
        allocate (xval_old(nkt*nj),STAT = istat)
        if(istat .ne. 0)
     .          call allocate_error("xval_old, solve_non_lin",0,istat)
      endif
      if( .not. allocated(xval_global))then
        allocate (xval_global(nkt*nj),STAT = istat)
        if(istat .ne. 0)
     .       call allocate_error("xval_global, solve_non_lin",0,istat)
      endif
      if( .not. allocated(xscale_global))then
        allocate (xscale_global(nkt*nj),STAT = istat)
        if(istat .ne. 0)
     .       call allocate_error("xscale_global, solve_non_lin",0,istat)
      endif
      if( .not. allocated(typf_global))then
        allocate (typf_global(nkt*nj),STAT = istat)
        if(istat .ne. 0)
     .       call allocate_error("typf_global, solve_non_lin",0,istat)
      endif
      if( .not. allocated(typx_global))then
        allocate (typx_global(nkt*nj),STAT = istat)
        if(istat .ne. 0)
     .       call allocate_error("typx_global, solve_non_lin",0,istat)
      endif
      if( .not. allocated(typx))then
        allocate (typx(nkt*nj),STAT = istat)
        if(istat .ne. 0)
     .          call allocate_error("typx, solve_non_lin",0,istat)
      endif

      if( .not. allocated(change))then
        allocate (change(nkt*nj),STAT = istat)
        if(istat .ne. 0)
     .          call allocate_error("change, solve_non_lin",0,istat)
      endif

    
      tipts =0
      uniform_weight = .false.
      xval_global_set = 0
      SSQR_global = 1.d100            !arbitrary  large number to start with
      init_nl = 1 !make sure problem  is un-frozen start of cycle
      call set_freeze(freeze_type,itte,itti,itangrot,itran,
     .                freeze_xte,freeze_xti,freeze_xwr,
     .                freeze_xrbp,freeze_xni,itran_save,
     .                 itti_save,itte_save,itangrot_save,
     .                 init_nl,frozen,nk,iangrot,jac_known)
      init_nl = 0
c       commented ot 4/8/03 HSJ==============================
c      if(freeze_nl .eq. 1) !set in tport if dt is too small 
c         ! we have reached dt lt 2*dtmin . This means we are not
c         ! converging . Freeze something according to freeze_type
c         !-5 to -1 freezes the profiles so they drop out of the
c         !problem
c     .   call set_freeze(freeze_type,itte,itti,itangrot,itran,
c     .                freeze_xte,freeze_xti,freeze_xwr,
c     .                freeze_xrbp,freeze_xni,itran_save,
c     .                itti_save,itte_save,itangrot_save,
c     .                init_nl,frozen,nk,iangrot,jac_known)
c=========================================      

        

      
c
c      freeze model after this number of iterations
c      if this is a time dependent run(for time independent runs
c       see below):
      if(iters_freeze .eq. 0) iters_freeze = 20 !normally set in inone,default 20
      restart = 0                               !set to > tot_iters_max to turn off
      navg_chi = 0              !determines # terms used to avg chi in diffuse
      tot_iters = 0
      ihalve_nl = 0             !returned as 1 if time step needs to be cut down
      info_conv = 0             !will be set to 1 if satisfactory solution is found
      if(tot_iters_max .eq. 0) tot_iters_max = 100
      iendit = 0 
      rescale = .true. 
      maxscale = .false. 
      freeze_xnue = 0  !xnue in qdimpl is allowed to change
      freeze_source = 0 ! call sub source
c      freeze_source = 1 !temp test only HSJ JAN 18
      nrmax = nj-1 ! values on boundary are known, hence do not form part of the equation set
       bandwidth_save = bandwidth
       itry1000 = 0
 1000  continue    ! loop back  to here if we drop some equations due to
                   ! lack of convergence
       imonit_old = 0
       ictr = 0
       itry1000 = itry1000+1
       call get_vars(maxscale,xvset,typx,
     .      nameu_gp,grid_point,x,xval,xval_old,nu)
      if(nu .le. 0)
     .   call STOP ('sub solve_newton: no. equations to solve = 0', 0)

c----set some arrays to dynamic size:
      if( .not. allocated(Jc))then
        allocate (Jc(nu,nu),STAT = istat)
        if(istat .ne. 0)
     .          call allocate_error("Jc, solve_non_lin",0,istat)
      endif



c     jacobian is banded, we have nkt variables in simulation mode
c     and each term in the residual equations depends on grid
c     points j-1,j,j+1 . This give us:
      eff_bandwidth = 2*nkt -1
      if(itry1000  .ge. 2 .and. freeze_type .lt. 0)then
             bandwidth = MAX(bandwidth_save - 3,0)
      endif 
      if(bandwidth .ne. 0 )eff_bandwidth = eff_bandwidth + bandwidth
      njup = eff_bandwidth
      njlow =eff_bandwidth
      printcode = 6       !need to protect write statments 
                          !before setting this to 0


      nback =0   ! back average,used if ddebug(3) .lt. 0,in io.i
       global = non_lin_method      ! set global from namelist input
       if(global .eq. 1)then
          solmethod = 'line search'
       else if(global .eq. 2)then
          solmethod = 'hookstep'
       else if(global .eq. 3)then
          solmethod = 'dogleg'
       else
          call STOP('solve_newton, no method selected',0)
       endif
       if(global .lt. 1 .or. global .gt. 3)
     .         call STOP("SUBROUTINE solve_newton,global incorrect",0)

      analjac = .false.
      noscale = .false.         !scaling done automatically
      if(get_typf .eq. 0)typf(:) = 1.d0
      cheapf =  .true.        !cheapf is used only if analjac = .false.
      sparse_jac = .false.
      factsec = .false.        !used if analjac = .false. and cheapf = .false.
      if(abs(jacobian_type) .eq. 1 .or. jacobian_type .eq. 3) 
     .                                         cheapf = .false.
      if(jacobian_type .ge. 2 .or. jacobian_type .eq. 3) 
     .                                      sparse_jac = .true.
      if(jacobian_type .eq. 3 .or. jacobian_type .eq. 1) 
     .                                        factsec = .true.


      if(fdigits .eq. 0) fdigits = 16            !could try -1 here
      if(fvectol .eq. 0.0d0)fvectol = 1.d-6        !maximum residual in any single equation
      if(gradtol .eq. 0.0d0)gradtol = 1.e-8          !min  gradient
      if(steptol .eq. 0.0d0) steptol = 1.d-10        !min step size
      maxstep = 0.05        !max step size
      if(maxfev .eq. 0 .or. maxfev .gt. tot_iters_max) 
     .                            maxfev = tot_iters_max      
      if(ssqrmin .le. 0.0d0) ssqrmin   = 1.d-8 !max sum squares residuals
      term3 = 0
      ipass =0
      ipassmax = switch_method  !number of times we are allowed switch
                                !solution methods !before quiting.
       ipassmax = MAX(ipassmax,1)
       ipturb = 0
       ipturbmax = random_pert   !random_pert will not generally be
                                 !an issue for time dependent cases.
                                 !(instead the time step is cut in half if
                                 !we get "stuck"). For steady state
                                 !cases if we get stuck at a local minimum
                                 !we will randomly perturb the solution
                                 ! a bit to get unstuck. Up to ipturbmax
                                 ! such perturbations may be carried out.
                                 !xval_global saves the best result in this
                                 ! case and is used as the returned result.
                                 !if freeze_type .ne. 0 then the solution
                                 !returned (if convergence to the desired
                                 !degree is not achieved)depends on
                                 ! tot_iters_max. if 
                                 !the number of iterations exceeds
                                 ! tot_iters_max before ipturbmax perturbations
                                 !are done then freeze_type is ignored.
                                 !otherwise if ipturbmax perturbations have
                                 !been done but  the number of iterations 
                                 ! is less than  tot_iters_max then freeze_type
                                 !is used to freeze some aspect of the problem
                                 !and more iterations are taken. if this effect
                                 !is not wanted then simply set
                                 !freeze_type =0 in inone.





 100  ipass = ipass +1         !loop back to here if solution found
                               !is not converged
      print *,'ipass,ipassmax,term3  =',ipass,ipassmax,term3
      print *,'ipturb,ipturbmax,time =',ipturb,ipturbmax,time


      ipassif: if(ipass .le. ipassmax)then
         delta = -1.0        !means trust region radius is to be initialized
c        generate the solution: 
c         do j=1,nu
c           xval_in(j) = xval(j)
c         enddo
         maxstep = 0.05
         if(steady_state .eq. 0.0)    jac_known = .false.
         call nl_driver_12(nu,xval,eval_non_lin_eq,global,analjac,
     .              cheapF,njlow,njup,sparse_jac,Jc,
     .              factsec,typF,fdigits,fvectol,steptol,
     .              maxstep,maxfev,delta,printcode,get_typf,
     .              solmethod,
     .              noscale,xval_in,termcode,gradtol,
     .              jmaxresid,jmaxgrad,frozen,
     .              freeze_type,ssqrmin,zero_step,first_step)
         xval(:)= xval_in(:)
         tot_iters =tot_iters + itncount         !itncount is iters taken for current method
                                                 !tot_iters is total for this call to solve_newton
         print *,'tot_iters,freeze_source =',tot_iters,freeze_source
         print *,'tot_iters_max,termcode  =',tot_iters_max,termcode
         print *,'freeze_xsn ,conv_skip =',freeze_xsn,conv_skip
c        if(tot_iters .ge. tot_iters_max .and. termcode .lt. 1) go to 500
         if(tot_iters .gt. tot_iters_max .and. conv_skip .eq. 0
     .      .and. ssqr .gt. ssqrmin )then
               go to 500                                !reduce time step
         else if(tot_iters .gt. tot_iters_max .and. conv_skip .eq.1)then
               termcode = 9                       !skip to next time step
         endif
         if(get_typf .eq. 0) termcode = -1

         if(get_typf .eq. 0 .and. uniform_weight) then
            typfavg  = 0.0
            do j=1,nu
                typfavg  = typfavg + typf(j)
            enddo
            typfavg = typfavg/nu
            do j =1,nu
              typf(j) = typfavg
            enddo
         endif



         if(tot_iters .gt.iters_freeze .and. freeze_xsn .eq. 1 .and.
     .       steady_state .gt. 0.0  .and. ssqr .gt. ssqrmin )then
              freeze_xnue = 1      !xnue in qdimpl will not be changed anymore
              freeze_source = 1
              if(ABS(steady_state -1.0) .lt. 1.e-5)then        !time dep run
                     print *,'calling freeze tot_iters =',tot_iters
                     print *,'iters before freeze =',iters_freeze
                     print *,'line 1437 ,itte_save =',itte_save
                     call set_freeze(freeze_type,itte,itti,itangrot,
     .               itran,freeze_xte,freeze_xti,freeze_xwr,
     .               freeze_xrbp,freeze_xni,itran_save,
     .               itti_save,itte_save,itangrot_save,
     .               init_nl,frozen,nk,iangrot,jac_known)
                     if(itry1000 .eq. 1)go to 1000
              else if(tot_iters .gt. 0.85*tot_iters_max )then
                     print *,'tot_iters_max =',tot_iters_max
                     print *,'line 1445 ,itte_save =',itte_save
                     call set_freeze(freeze_type,itte,itti,itangrot,
     .                     itran,freeze_xte,freeze_xti,freeze_xwr,
     .                     freeze_xrbp,freeze_xni,itran_save,
     .                     itti_save,itte_save,itangrot_save,
     .                     init_nl,frozen,nk,iangrot,jac_known)
                     if(itry1000 .eq. 1)go to 1000
              endif
              xval_old(:) = xval(:)
         endif



         conv_rate = 1.e6
         imonit = tot_iters/10
         if(imonit-imonit_old .gt. 0 .and. tot_iters .ge. 30)then
              imonit_old = imonit
              ictr = ictr+1
              SELECT CASE(ictr)
                  CASE(1)
                      xmonit(ictr) = SSQR
                  CASE(2)
                      xmonit(ictr) = SSQR
                      conv_rate = 100.*ABS(xmonit(2)
     .                                       -xmonit(1))/xmonit(2)
                  CASE(3)
                      xmonit(1) =xmonit(2)
                      xmonit(2) = SSQR
                      conv_rate = 100.*ABS(xmonit(2)
     .                          -xmonit(1))/xmonit(2)
                      ictr  = 2
              END SELECT
         endif ! monit
         if(conv_rate .lt. .2 .and. SSQR .gt. SSQRMIN .and.
     .               termcode .ne. 6 .and. termcode .ne. 1)then
             termcode = 3
             term3 = 3
         endif 
         if(termcode .gt. 0 .and. SSQR .le. ssqrmin)termcode = 8

         SELECT CASE(termcode)
                CASE(-1)
                           write(printcode,'("determined typf values")')
                           get_typf = 1
                           ipass =ipass - 1   !get typf passes dont count
c                          do j=1,nu
c                               print *,'typf(j) =',typf(j)
c                               print *,'typx(j) =',typx(j)
c                          enddo
                           go to 100
                CASE(0)
                            write(printcode,
     .                      '("error nl_driver_12 returned with",/,       
     .                       " termcode = 0. Indicates bug in code")')
                            Call STOP(
     .                      'Subroutine solve_newton: error return',0)
                CASE(1)
                            print *,'case 1'
                            write(printcode,
     .                      '("Norm of scaled gradient is less than",
     .                      1pe12.5,/,
     .                      "this indicates a good solution unless",
     .                      " gradtol is set",/,
     .                      "too large for this problem")')gradtol
c                            if(SSQR .le. 10.*ssqrmin)info_conv = 1
                             info_conv = 1
                            write(printcode,
     .                     '("sum of squares of residuals =",
     .                                                1pe12.5)')SSQR
                           write(printcode,
     .                        '("info_conv =",i5)')info_conv
                CASE(2)
                           print *,'case 2'
                           write(printcode,
     .                     '("scaled distance between last two ",/,
     .                          "steps is less than ",1pe12.5,/,
     .                          "may indicate converged result if ",/,
     .                          " sum of squares of residuals is "
     .                          "sufficiently small")')steptol
                            write(printcode,
     .                      '("sum of squares of residuals =",
     .                                                 1pe12.5)')SSQR
                            if(ssqr .gt. ssqrmin)then
                               write(printcode,
     .                             '("run will be restarted from last ",
     .                        "estimate to try and reduce ssq to below "
     .                        " ssqrmin = ",1pe12.4)') ssqrmin
                              if(global .eq. 3)then
                                 global = 1
                              else
                                 global = global + 1
                              endif
                              if(global .eq. 1)solmethod = 'linesearch'
                              if(global .eq. 2)solmethod = 'hookstep'
                              if(global .eq. 3)solmethod = 'dogleg'
                              write(printcode,
     .                           '("selected new solution method ",a)')
     .                        solmethod
                              if(rescale)then
c                                get_typf = 0   ! redetermine the scaling in the equations
                                 do j=1,nu      !and in the independent variables
                                    typx(j) = ABS(xval(j))
c                                   typf(j) = 1.0
                                 enddo
                               endif
                               go to 100
                            else   ! ssqr .le. ssqrmin solution OK
                               info_conv = 1 
                            endif
                CASE(3)
                    print *,'case 3'
                    write(printcode,
     .                      '("last step failed to decrease",
     .                      "the objective function",/,  "may "
     .                      "indicate that steptol is set too large",/,
     .                      "current value of steptol =",1pe12.4)')
     .                                                          steptol
                    write(printcode,
     .                      '("sum of squares of residuals =",
     .                                                   1pe12.5)')SSQR
                    !if we quit on zero step size 
                    !three times in a row
                    !then all three methods have failed to 
                    !improve the current estimate:
                    if(itncount .le. 2)term3 = term3 +1 
                    if(ssqr .gt. ssqrmin .and. term3 .le. 2)then
                           write(printcode,
     .                             '("run will be restarted from last ",
     .                            "estimate to try and reduce ssq"
     .                            " to below ssqrmin = ",
     .                            1pe12.4)') ssqrmin
                           write(printcode,
     .                            '("we will continue the run with ",
     .                             "a different solution method")')
                           if(global .eq. 3)then
                                       global = 1
                                       solmethod = 'linesearch'
                           else
                                       global = global + 1
                           endif
                           if(global .eq. 2)solmethod ='hookstep'
                           if(global .eq. 3)solmethod = 'dogleg'
                           write(printcode,
     .                            '("selected new solution method ",a)')
     .                             solmethod
                           if(rescale)then
c                                    get_typf = 0   ! redetermine the scaling in the equations
                                     do j=1,nu      !and in the independent variables
                                       typx(j) = ABS(xval(j))
c                                      typf(j) = 1.0
                                     enddo
                            endif
                            go to 100
                    else  if( ssqr .gt. ssqrmin .and. term3 .ge. 3)then
                        term3 = 0
                        ind_test = -1
                        ipturb = ipturb+1



                        if(ipturb .le. ipturbmax)then
                             SSQR_global = MIN(SSQR_global,SSQR)
                             if(SSQR_global .eq. SSQR)then
                                xval_global(:) = 
     .                              xval(1:SIZE(xval_global))
                                xscale_global(:) = 
     .                               xscale(1:SIZE(xscale_global))
                                typf_global(:) = 
     .                                   typf(1:SIZE(typf_global))
                                typx_global(:) = 
     .                              typx(1:SIZE(typx_global))
                                itran_global(:) = itran(:)
                                 xval_global_set = 1
                             endif
                             IF(wrt_nwt .NE.  0)
     .                       call output_test(nu,xval,restart,
     .                               termcode,ipturb,resid,frozen,
     .                               ind_test,nameu_gp,typf,grid_point,
     .                                           ssqr,NOSCALE,iendit)
                             !we cant get a better solution. Freeze 
                             !something according to freeze_type:
c                            if(ABS(steady_state -1.d0) .lt. 1.e-5)then
c                            if this is a time dependent run 
c                            call set freeze
c                            if it is a time independent run 
c                            (steady_state =0)
c                            then call set_freeze if the number 
c                            of random perturbations
c                            has reached ipturbmax:
                             if((ABS(steady_state -1.d0) .lt. 1.e-5)
     .                           .or. (steady_state .lt. 1.e-5 
     .                           .and. ipturb .eq. ipturbmax))then
                             print *,'line 1636 ,itte_save =',itte_save
                                      call set_freeze(freeze_type,
     .                                        itte,itti,itangrot,itran,
     .                                freeze_xte,freeze_xti,freeze_xwr,
     .                                freeze_xrbp,freeze_xni,itran_save,
     .                                itti_save,itte_save,itangrot_save,
     .                                init_nl,frozen,nk,iangrot,
     .                                                        jac_known)


                             else
                                 !maximum gradient component will be perturbed
                                 !grad (0.5*Ftranspose*F) = Jtranspose *F
                                 !where F (maps Rn into Rn ) is 
                                 !the  set of equations to be solved. 

                                 jacmax =0.0
                                 do j= 1,nu
                                   dumy = ABS(Jc(j,jmaxgrad))/typf(j)
                                   jacmax = MAX(dumy,jacmax)
                                   if(jacmax .eq. dumy)jcmax = j 
                                 enddo
                                 if(Jc(jcmax,jmaxgrad) .gt. 0.0d0)then
                                      xval(jcmax) = 0.999*xval(jcmax)
                                 else
                                      xval(jcmax) = 1.001*xval(jcmax)
                                 endif
                                 IF( .not. allocated(xcRandom) )
     .                           allocate ( xcRandom(nu),STAT = istat)
                                 CALL Random_Number(xcRandom) 
                                 print *,'///////////////////////////',
     .                                   '///////////////////////////'
                                 print *,'         PERTURBING SOLUTION:'
                                 write(*,1011)(xcrandom(j),j=1,nu)
 1011                            format(2(2x,1pe14.4))
                                 do j=1,nu
                                   if(ssqr .gt. 1.0)then
                                     xval(j) = xval(j)*
     .                               DBLE(1.+(xcRandom(j) - 0.5)*.001)
                                   else
                                     xval(j) = xval(j)*
     .                               DBLE(1.+(xcRandom(j) - 0.5)*.001)
                                   endif
                                 enddo
                                 go to 100
                             endif              !ABS(steady_state -1.d0) case
                             if(frozen .eq. 1 .and. freeze_type .lt. 0)
     .                                                             then
                               go to 1000
                             else if(frozen .eq. 1 .and. freeze_type 
     .                                                  .gt. 0)then
                               !freeze_type .gt. 0 does not involve
                               !changing number of unknowns so we
                               !continue at label 100 
                               go to 100
                             else
                               print *,'dt,dtmin,time =',dt,dtmin,time
                               info_conv = 0
c                              Call STOP('Subroutine solve_newton:no', freeze 
c     .                                   possible',0)
c                              instead of stopping continue with 
c                              this solution: HSJ 4/29/02
c                              info_conv = 1                         
                             endif                       
                        endif      !ipturb .le. ipturbmax case
                        print *,'case 3 ipturb .gt. ipturbmax'
c                    else if(ssqr .gt. ssqrmin)then
c                        write(printcode,'("used all 3 methods without ",
c     .                  "further progress",/,"time step will be",/,
c     .                  " aborted")')
c                        info_conv =0
c                        ihalve_nl = 1
c                    else if(ssqr .le. ssqrmin)then
c                        info_conv = 1
                    else
                        Call STOP('Subroutine solve_newton:case3',0) 
                    endif



                CASE(4)
                    print *,'case 4,termcode,ipass,ipassmax',
     .              termcode,ipass,ipassmax
                    if(ipass .le. ipassmax)then
                        write(printcode,'("iteration limit per method",
     .                  " was exceeded",/,
     .                  "run is not yet converged, we will reset the",
     .                  "  iteration counter  to zero. ",/,
     .                  " max iterations per method =",i6,/,
     .                  "  total iterations  taken = ",i6
     .                                )')maxfev,tot_iters
                        write(printcode,'("run will be restarted from"
     .                        " last estimate to try and reduce ssq "
     .                        " to below ssqrmin = ",1pe12.4)') ssqrmin
                        write(printcode,'("we will continue the run "
     .                        "with a different solution method")')
                        if(global .eq. 3)then
                           global = 1
                        else
                           global = global + 1
                        endif
                        if(global .eq. 1)solmethod = 'linesearch'
                        if(global .eq. 2)solmethod = 'hookstep'
                        if(global .eq. 3)solmethod = 'dogleg'
                        write(printcode,'("selected new solution"
     .                        " method ",a)') solmethod
                        term3 =0
                        ind_test = -2
                        if(wrt_nwt .NE.  0)
     .                   call output_test(nu,xval,restart,termcode,
     .                               ipturb,resid,frozen,ind_test,
     .                   nameu_gp,typf,grid_point,ssqr,NOSCALE,iendit)
                        if(rescale)then
                            print *,'rescaling equations ,termcode =',
     .                                                     termcode
c                           get_typf = 0   ! redetermine the scaling 
                                           !in the equations
                            do j=1,nu      !and in the independent variables
                               typx(j) = ABS(xval(j))
c                              typf(j) = 1.0
                            enddo
                        endif

                        if(frozen .eq. 1 .and. freeze_type .lt. 0)
     .                                                    go to 1000
                        !freexe_type .gt. 0 does not involve
                        !changing number of unknowns so we
                        !continue at label 100 
                        if(ipass .eq. ipassmax .and. frozen .eq.1 ) 
     .                                            ipass = ipass-1
                        go to 100
                 else       !ipass .gt. ipassmax
                    write(printcode,'("run will be terminated due to ",
     .                 " maximum number of passes being restricted to",
     .                  i5,/"termcode =",2x,i5)')ipassmax,termcode 
                    info_conv =0
                 endif
c                 print *,'reached end of case 4'


              CASE(5)
                 print *,'case 5'
                 write(printcode,'(" five consecutive steps of size"
     .           "maxstep =",1pe12.4,/,
     .           "were taken. we will contine with larger value of",/,
     .           " maxstep ")')maxstep
                 maxstep = maxstep *1.5
                 if(rescale)then
c                   get_typf = 0   ! redetermine the scaling in the equations
                   do j=1,nu      !and in the independent variables
                      typx(j) = ABS(xval(j))
c                      typf(j) = 1.0
                   enddo
                 endif
                 term3 = 0
                 ipassmax = ipassmax+1 !otherwise we will exit at label 100
                 go to 100
                CASE(6)
                    print *,'case 6'
                    write(printcode,'("maximum residual of any single ",
     .                 "equation was less than ",1pe12.4,/,
     .                 "Satisfactory solution found")')fvectol 
                    info_conv = 1 
                CASE(7)
                    print *,'case 7'
                    write(printcode,'("maximum iters exceeded ",
     .                 "try to cut time step in half")') 
                    info_conv = 0 
                    ihalve_nl = 1   !cut time step in half
                CASE(8)
                    print *,'case 8'
                    write(printcode,'("converged: SSQR le ssqrmin ")') 
                    write(printcode,'("SSQr,ssqrmin =",2(2x,1pe12.4))')
     .                        SSQR,ssqrmin
                    info_conv = 1
                CASE(9)
                    print *,'case 9'
                    !conv_skip =1, so dont test convergence after
                    !tot_iters_max iterations are done.
                    write(printcode,'("convergence test skipped")') 
                    write(printcode,'("SSQr,ssqrmin =",2(2x,1pe12.4))')
     .                        SSQR,ssqrmin
                 info_conv = 1
         END SELECT


      else if(frozen .eq. 0)then ipassif !ipass branch
           write(printcode,'("did not converge in maximum",
     .              i5," passes.",/,
     .     " Will continue with xchis frozen until this time",
     .     " step is converged")')ipassmax
           xval_old(:) = xval(:)
           ipass = MAX(0,ipassmax -10)
           !we cant get a better solution. Freeze 
           !something according to freeze_type:
           print *,'line 1830 ,itte_save =',itte_save
           call set_freeze(freeze_type,itte,itti,itangrot,itran,
     .                freeze_xte,freeze_xti,freeze_xwr,
     .                freeze_xrbp,freeze_xni,itran_save,
     .                 itti_save,itte_save,itangrot_save,
     .                 init_nl,frozen,nk,iangrot,jac_known)
           if(frozen .eq. 1 .and. freeze_type .lt. 0)then
                 print *,'location d,frozen =',frozen
                 if(ipass .eq. ipassmax) ipass =ipass-1
                 if(itry1000 .eq. 1)go to 1000
           else if(frozen .eq. 1 .and. freeze_type .gt. 0)then
                           !freexe_type .gt. 0 does not involve
                           !changing number of unknowns so we
                           !continue at label 100 
                 if(ipass .eq. ipassmax) ipass =ipass-1
                 go to 100
           else if( steady_state .gt. 1.e-5 .and. conv_skip .eq.0 )then
                 Call STOP('solve_newton:no freeze possible 2',0)
           endif

      endif   ipassif               !ipass .lt. ipassmax





c     generate the final solution vector returned by nl_driver_12:

 500  if(xval_global_set .eq. 1 .and. steady_state .lt. 1.e-5 
     .            .and. SSQR .gt. SSQR_global)then
         xval(1:SIZE(xval_global)) = xval_global(:)     !for the steady state case use the
         xscale(1:SIZE(xscale_global)) = xscale_global(:) !best solution found during the
         typf(1:SIZE(typf_global)) =   typf_global(:)   !random perturbations (if there were any).
         typx(1:SIZE(typx_global)) =   typx_global(:)
         sf(:)   =   1./typf(:)
         sx(:)   =   1./typx(:)
         SSQR = SSQR_global       
                                  
         itran(:) = itran_global(:)
      endif
      nu = 0
      do j =1,nrmax
         do k=1,nk
            if(itran(k) .gt. 0)then
c              skip rbp at the origin (its value is fixed at 0.0)
               if(j .eq. 1 .and. k .eq. nk-iangrot)go to 200
               if(j .ge. te_index .and. k .eq. nk-iangrot-2)go to 200
               if(j .ge. ti_index .and. k .eq. nk-iangrot-1)go to 200
               if(j .ge. rot_index .and. k .eq. nk .and. 
     .                                          iangrot .eq. 1)go to 200
               if(j .ge. ni_index .and. k .eq. 1)go to 200 !jmp.den     
               nu = nu+1
               if(k .le. nk-iangrot)then
chsj                u(k,j) = (xscale(nu)*xval(nu))**2 !forces non-negative value
                    u(k,j) = (xval(nu))**2 !forces non-negative value
               else
chsj                u(k,j) = xscale(nu)*xval(nu)
                    u(k,j) = xval(nu)
               endif
               change(nu) = u(k,j)-usave(k,j)
 200        endif
         enddo
      enddo


c ---------------------------------------------------------------------
      !unfreeze if necessary:
      init_nl = 1 !make sure problem  is un-frozen start of cycle
cjmp.c  print *,'line 1896 ,itte_save =',itte_save
c      call set_freeze(freeze_type,itte,itti,itangrot,itran,
c     .                freeze_xte,freeze_xti,freeze_xwr,
c     .                freeze_xrbp,freeze_xni,itran_save,
c     .                 itti_save,itte_save,itangrot_save,
c     .                 init_nl,frozen,nk,iangrot,jac_known)
c ------------------------------------------------------------------------
c      get residual of final result:
      call eval_non_lin_eq(nu,xval,resid)
      resid_max =0.D0
      if(newtonvb .eq. 1)then
           nu = 0
           SSQR = 0.0
           do j=1,nrmax
              do k=1,nk
                 if(itran(k) .gt. 0)then
                    if(j .eq. 1 .and. k .eq. nk-iangrot)go to 20
                    if(j .ge. te_index .and. k .eq. nk-iangrot-2)
     .                                                   go to 20
                    if(j .ge. ti_index .and. k .eq. nk-iangrot-1)
     .                                                   go to 20
                    if(j .ge. rot_index .and. k .eq. nk .and. 
     .                                    iangrot .eq. 1)go to 20
                    if(j .ge. ni_index .and. k .eq. 1)go to 20 !jmp.den
                    nu = nu + 1 
                    SSQR = SSQR + (resid(nu)/typf(nu))**2
                    print *,'j,resid,change scalef =',j,nameu_gp(nu),
     .                         resid(nu)/typf(nu),change(nu),typf(nu)
                    resid_max = MAX(resid_max,ABS(resid(nu)))
                 endif
 20           enddo
              print *,' ' 
           enddo
           print *,' SSQR = ',SSQR 
      else
         print *,' max resid,tot iters = ',resid_max,tot_iters
         print *,'time,info convg, dt =',info_conv,dt
      endif
c update sets en,te,ti,rbp,angrot to the values in vector u:
      if(iteration_method .ne. 1)  !default iteration_method = 0
     .call update (u, en, te, ti, rbp, nk, nj, kj, kk, iangrot, angrot)
      print *,'termcode,frozen', termcode,frozen
 



c      print *,'line 1963 ,itte_save =',itte_save
      if((tot_iters .gt. tot_iters_max  .and. conv_skip .eq. 0) .and.
     .   SSQR .ge. ssqrmin )then
           !set flag to indicate time step needs to be reduced
           !if this is a time stepping case. If this is a run
           !to steady state then we will terminate the run gracefully
           if(ABS(steady_state -1.0d0) .lt. 1.e-5 .and. frozen 
     .                      .eq. 0 )then

                  call set_freeze(freeze_type,itte,itti,itangrot,itran,
     .                                freeze_xte,freeze_xti,freeze_xwr,
     .                               freeze_xrbp,freeze_xni,itran_save,
     .                               itti_save,itte_save,itangrot_save,
     .                             init_nl,frozen,nk,iangrot,jac_known)
               tot_iters = tot_iters/2  !allow half as many
               tot_iters = MAX(tot_iters,1)
               print *,'frozen,freeze_type', frozen,freeze_type
               if(frozen .eq. 1 .and. freeze_type .lt. 0)then
                print *,'location e,frozen =',frozen
                   if(itry1000 .eq. 1)go to 1000
               else if(frozen .eq. 1 .and. freeze_type .gt. 0)then
                           !freexe_type .gt. 0 does not involve
                           !changing number of unknowns so we
                           !continue at label 100 
                   go to 100
               else

              ihalve_nl = 1
              write(printcode,'("max number of iterations exceeded",/,
     .        " will cut time step in half and try again")')
              !replaced stop with above 8/26/02 HSJ
c                    Call STOP('solve_newton:no freeze possible',0)
                    !no freeze possible
               endif
           else if(ABS(steady_state -1.0d0) .lt. 1.e-5 .and. frozen 
     .                      .eq. 1 )then
              ihalve_nl = 1
              write(printcode,'("max number of iterations exceeded",/,
     .        " will cut time step in half and try again")')
           else
              write(printcode,'("max number of iterations exceeded",/,
     .              "steady state run will be terminated")')
              go to 5000
           endif
      endif
      init_nl = -1 !restore problem setup
c      print *,'line 2015 itte_save=',itte_save
      call set_freeze(freeze_type,itte,itti,itangrot,itran,
     .                freeze_xte,freeze_xti,freeze_xwr,
     .                freeze_xrbp,freeze_xni,itran_save,
     .                 itti_save,itte_save,itangrot_save,
     .                 init_nl,frozen,nk,iangrot,jac_known)
      print *,'itte =',itte
      print *,'itran_save =',itran_save
      print *,'itran =',itran
c      call stop('line 2012',1)
      init_nl = 0
      bandwidth =  bandwidth_save
      if(ihalve_nl .eq. 1)then  !save the current best guess. After the
         xval_old(:) =  xval(:) !time step is halved in tport and we 
         xvset = .true.         !renter this routine we will use xval_old
      else                      !to produce a better initial gues
         xvset = .false.
      endif
      print *,'xvset =',xvset

c     if this run is not checked for convergence then we must
c     ensure that we have not stalled because we ended up at a
c     point where no further steps could be taken due to zero steplength
c     In this case we swtich solution methods so that a different
c     step (of non zero lenght) can  be taken:
c      if(conv_skip .eq. 1 .and. zero_step)then
  
      if(conv_skip .eq. 1 )then
         if(non_lin_method .eq. 3)then
                non_lin_method = 1
         else
                non_lin_method = non_lin_method + 1
         endif

      endif
 5000  continue
       if(steady_state  .lt. 1.d-5)then
          iendit = 1
          if(wrt_nwt .NE. 0)
     .     call output_test(nu,xval,restart,
     .        termcode,ipturb,resid,frozen,
     .        ind_test,nameu_gp,typf,grid_point,
     .                     ssqr,NOSCALE,iendit)
        endif





       if(newtonvb .eq. 1)print *,'leaving solve_newton'




      return
      end
