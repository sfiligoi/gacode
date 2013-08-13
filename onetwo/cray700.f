  
      
      SUBROUTINE  broyfac(n,xc,xp,fc,fp,eta,Sx,Sf,Z,M,M2)
c------------------------------------------------------------------------
c--- Factored form of Broydens update for Jacobian
c--- INPUT:
c--- INPUT/OUTPUT:
c--- z(n,n)  contains Q transpose of A on input
c            Q+, transpose of A+ on output
c--- LM( )   LM contains R (factor of matrix a)
c            on input and R+, factor of matrix A+ on output
c            
c -- M2(n)   contains diagonals of M
c------------------------------------------------------------------------



      IMPLICIT none
c  INPUT:
      INTEGER n
      REAL *8 xc(n),xp(n),fc(n),fp(n),Sx(n),Sf(n),eta
C  LOCAL:
      INTEGER
     .        i,j
      LOGICAL skipupdate
      REAL*8
     .       s(n),t(n),w(n),sumd,denom
c INPUT/OUTPUT:
      REAL*8 M(N,N),Z(n,n),M2(n)


      DO i =1,n
         M(i,i) = M2(i)
          s(i) = xp(i)-xc(i)
      ENDDO
      skipupdate = .true.

      DO i =1,n
         sumd =0.0D0
         DO j = i,n
            sumd =sumd +M(i,j)*s(j)
         ENDDO
         t(i) = sumd   !t = R*s  R in upper triangle of M
      ENDDO
      DO i =1,n
         sumd =0.0D0
         DO j = 1,n
            sumd =sumd +Z(j,i)*t(j)
         ENDDO
         w(i) = Sf(i)*(Fp(i) -Fc(i)) - sumd
         IF(w(i) .gt. eta * Sf(i)*(ABS(Fp(i)+ABS(Fc(i)))))THEN
            skipupdate = .false.
         ELSE
            w(i) = 0.0D0
         ENDIF
      ENDDO
      

      IF(.not. skipupdate )THEN
         DO i =1,n
            sumd =0.0D0
            DO j = 1,n
               sumd =sumd +Z(i,j)*w(j)
            ENDDO
            t(i) = sumd
         ENDDO
         sumd =0.0D0
         DO i=1,n
            sumd = sumd + (Sx(i)*s(i))**2
         ENDDO
         denom = sumd
         DO i =1,n
            s(i) = Sx(i)**2*s(i)/denom
         ENDDO

         CALL QRUPDATE(n,t,s,1,Z,M)

         DO i =1,n
           M2(i) = M(i,i)
         ENDDO
      ENDIF


      RETURN
      END

         
      SUBROUTINE broyunfac(n,xc,xp,FVc,Fvp,eta,Sx,Jc)
c---------------------------------------------------------------
c --- unfactored Broyden method,not verys useful ??
c -- not implemented at this time 
c-----------------------------------------------------------HSJ
      USE terminate
      IMPLICIT NONE
!INPUT:
      INTEGER 
     .       n
      REAL *8
     .       xc(n),xp(n),Jc(n,n),FVc(n),FVp(n),
     .       Sx(n),eta
      CALL STOP("SUB BROYUNFAC NOT IMPLEMENTED",6)
      RETURN
      END
C
C




      SUBROUTINE CHECK_CONVERGENCE(N,XCUR,XNEW,FVp,GRADIENT,FNORM,
     &             SCALEX,SCALEF,RETCODE,FVECTOL,GRADTOL,STEPTOL,
     &             NOSCALE,ITNCOUNT,ITNLIMIT,MAX_TAKEN,CONSECMAX,
     &             RETURN_CODE,IOUNIT,GRADMAX,SSQR,SOLMETHOD,
     &             JMAXRESID,JMAXGRAD,NAMEU_GP,GRID_POINT,FROZEN,
     &             FREEZE_TYPE,SSQRMIN)
C-----------------------------------------------------------------------
C---SUBROUTINE DECIDES IF THE PROBLEM IS CONVERGED.
C---INPUT
C  N                    NUMBER OF ELEMENTS IN XCUR,XNEW,GRADIENT,
C                       SCALE
C  XCUR(I)              I=1,2..N,CURRENT SOLUTION POINT
C  XNEW(I)              I=1,2..N,NEW SOLUTION POINT
C  FNORM                SUM SQUARES OF RESIDUALS At XNEW
C  FVP(J)               RESIDUAL of EQUATION J
C  GRADIENT(I)          I=1,2,..N,VALUE OF GRADIENT AT XNEW
C  SCALEX(I)             I=1,2,..N,SCALE FACTORS
C  NOSCALE              BOOLEAN,NOSCALE=.TRUE. MEANS NO SCALING IS
C                       IN EFFECT
C  SCALEF(I)            TYPICAL  MAGNITUDE OF F(I)
C  RETCODE              RETCODE =0 IF VALID XNEW WAS FOUND
C                       RETCODE =1 MEANS A VALID XNEW WAS NOT FOUND
C                       (BY SUBROUTINE LINESEARCH OR DOGDRIVER). IN
C                       THIS CASE WE ASSUME THAT XCUR IS THE SOLUTION
C                       AND RETURN WITH RETURN_CODE=3 (SEE BELOW)
C  FVECTOL              MAX RESIDUAL IN ANY EQUATION
C  GRADTOL              TOLERANCE LEVEL FOR THE (SCALED) GRDIENT
C  STEPTOL              TOLERANCE LEVEL FOR THE (SCALED) STEP SIZE
C  ITNCOUNT             ITERATION  COUNTER
C  ITNLIMIT             MAX ITERATIONS ALLOWED
C  MAX_TAKEN            NUMBER OF CONSECUTIVE STEPS OF MAXSTEP TAKEN
C  SOLMETHOD            character variables contains solution method in use
c  NAMEU_GP(j)          name of umknown J
c  GRID_POINT(j)        rho grid point associated with index j
C  FROZEN               = 0 nothng frozen, =1 freeze_type gives frozen model
C  FREEZE_TYPE          -5 to 5 indicates type of freeze applied
C  SSQRMIN              RETURN  CONVERGED IF SSQ < SSQRMIN
C---OUTPUT
C  CONSECMAX            COUNTS NUMBER OF STEPS OF MAX LENGTH
C                       TAKEN IN A ROW.
C  RETURN_CODE          RETURN_CODE=0  MEANS NOT CONVERGED
C                               =1  NORM OF (SCALED) GRADIENT IS
C                                   LESS THAN GRADTOL,THIS SHOULD BE
C                                   A LOCAL MINIMIZER UNLESS GRADTOL
C                                   IS SET TOO LARGE.
C                               =2  SCALED DISTANCE BETWEEN LAST TWO STEPS
C                                   IS LESS THAN STEPTOL. COULD MEAN
C                                   WERE CONVERGED,OR STEPTOL IS TOO LARGE
C                                   OR ALGORITHM IS MAKING VERY SMALL
C                                   STEPS TOWARD A MINIMIZER
C                               =3  LAST STEP FAILED TO LOCATE A LOWER
C                                   POINT THAN XCUR,XCUR IS THE SOLUTION
C                                   OR STEPTOL IS TOO LARGE
C                               =4  ITERATION LIMIT EXCEEDED
C                               =5  FIVE CONSECUTIVE STEPS OF SIZE
C                                   MAXSTEP WERE TAKEN. MAXSTEP
C                                   IS SET TOO SMALL.
C                               =6  RESIDMAX < FVECTOL
c                               =8  SSQR < SSQRMIN
C
C-------------------------------------------------------------------HSJ
C
C
      USE aid_newton,   only : wrt_nwt,tipts,tipts2,ssqr_hist
      USE tcoef,        only : d
      USE numbrs,       only : nion,nj
      USE nonlin,       only : typf
      USE ifs,          ONLY : xke_ifs,xki_ifs,xkang_ifs,dmassden
      USE machin,       ONLY : rmajor
      USE solcon,       ONLY : time
      IMPLICIT NONE
      INTEGER N,RETURN_CODE,CONSECMAX,ITNLIMIT,ITNCOUNT,RETCODE,J,
     &        JMAXRESID,JMAXGRAD,GRID_POINT(*),JMAXSTEP,
     &        FROZEN,FREEZE_TYPE,IOSTAT,K,IOUNIT
      INTEGER IO_NWT
      REAL*8 SCALEX(*),SCALEF(*),
     &     GRADTOL,STEPTOL,GRADMAX,FMAX,DUMY,ABS,
     &     DUMY1,STEPMAX,SSQR,FVp(n),DENAXIS,TEAXIS,TIAXIS,
     .     RBPEDGEM1,ANGRAXIS,FNORM,FVECTOL,RESIDMAX,FRMAX,
     .     GRMAX,SRMAX,XCR,XCD,SSQRMIN,DUMY2,DUMY3
      REAL *8 XCUR(*),XNEW(*),GRADIENT(*)
      LOGICAL MAX_TAKEN,NOSCALE ,exists, opened
      CHARACTER *(*) SOLMETHOD,NAMEU_GP(*)
C
C
      RETURN_CODE=0
      GRADMAX=0.0d0
      RESIDMAX =0.0d0
      FMAX=MAX(FNORM,DBLE(N/2))
      DO J=1,N
          DUMY=1.0d0 
          IF(.NOT. NOSCALE)DUMY=1./SCALEX(J)
          DUMY=MAX(ABS(XNEW(J)),DUMY)
          GRMAX = ABS(GRADIENT(J))*DUMY/FMAX
          GRADMAX=MAX(GRMAX,GRADMAX)
          IF(GRADMAX .EQ. GRMAX) JMAXGRAD = J
          FRMAX = ABS(FVP(j))*SCALEF(j)
          RESIDMAX = MAX(FRMAX,RESIDMAX)
          IF(RESIDMAX .EQ. FRMAX)JMAXRESID = J
      ENDDO
      RESIDMAX = RESIDMAX*RESIDMAX
      STEPMAX=0.0
      JMAXSTEP = 1
      DO J=1,N
          DUMY=1.0
          IF(.NOT. NOSCALE)DUMY=1./SCALEX(J)
          DUMY=MAX(ABS(XNEW(J)),DUMY)
          DUMY1=ABS(XNEW(J)-XCUR(J))
          SRMAX = DUMY1/DUMY
          STEPMAX=MAX(SRMAX,STEPMAX)
          IF(STEPMAX .EQ. SRMAX)JMAXSTEP =j
      ENDDO
      SSQR=2.*FNORM  !fnorm is 0.5 * sum of squares of residuals
      SSQR_HIST = SSQR
      IF(RETCODE .EQ. 1)THEN
          RETURN_CODE=3
      ELSE IF(GRADMAX .LE. GRADTOL)THEN
          RETURN_CODE=1
      ELSE IF(STEPMAX .LE. STEPTOL)THEN
          RETURN_CODE=2
      ELSE IF(ITNCOUNT .GE. ITNLIMIT)THEN
          RETURN_CODE=4
      ELSE IF(MAX_TAKEN )THEN
          CONSECMAX=CONSECMAX+1
          IF(CONSECMAX .EQ. 5)RETURN_CODE=5
      ELSE IF( FVECTOL .GE. SQRT(RESIDMAX))THEN
          RETURN_CODE = 6
      ELSE IF( SSQRMIN  .GE. SSQR)THEN
          RETURN_CODE = 8
      ELSE
          CONSECMAX=0
      ENDIF

      CALL GET_AXIS_VALUES(DENAXIS,TEAXIS,TIAXIS,RBPEDGEM1,
     .                                          ANGRAXIS,XNEW)
      IF(IOUNIT .NE. 0)THEN
          WRITE(IOUNIT,1)ITNCOUNT,GRADMAX,GRADTOL,STEPMAX,STEPTOL,
     &                   RESIDMAX,SSQR,FROZEN,FREEZE_TYPE
    1     FORMAT(2X,'ITER,GRADMAX,GRADTOL,STEPMAX,STEPTOL,RESIDMAXSQ,',
     &    'SSQRESID :',/
     &             2X,I3,6(1PE10.3,1X),'-',I2,X,I2)
          WRITE(IOUNIT,2)DENAXIS,ANGRAXIS,TEAXIS,TIAXIS,RBPEDGEM1,
     .                   SOLMETHOD
 2        FORMAT(2X,'DEN(1) = ',1pe14.6,2x,'ANGROT(1) = ',1pe14.6,/,
     .            2x,'TE(1) = ',1pe14.6,2x,'TI(1) = ',1pe14.6,/,
     .           2X,'RBP(nj-1) = ',1pe14.6,2x,a)

          WRITE(IOUNIT,3)nameu_gp(jmaxgrad)(1:8),grid_point(jmaxgrad)
          if(nameu_gp(jmaxgrad)(1:2) .ne. 'wr')then
            xcr = XCUR(jmaxgrad)*XCUR(jmaxgrad)
            xcd = XNEW(jmaxgrad)*XNEW(jmaxgrad)
          else
            xcr = XCUR(jmaxgrad)
            xcd = XNEW(jmaxgrad)
          endif
          WRITE(IOUNIT,4)xcr,xcd
          WRITE(IOUNIT,5)nameu_gp(jmaxresid)(1:8),grid_point(jmaxresid)
          if(nameu_gp(jmaxresid)(1:2) .ne. 'wr')then
            xcr = XCUR(jmaxresid)*XCUR(jmaxresid)
            xcd = XNEW(jmaxresid)*XNEW(jmaxresid)
          else
            xcr = XCUR(jmaxresid)
            xcd = XNEW(jmaxresid)
          endif
          WRITE(IOUNIT,4)xcr,xcd
          WRITE(IOUNIT,6)nameu_gp(jmaxstep),grid_point(jmaxstep)
          if(nameu_gp(jmaxstep)(1:2) .ne. 'wr')then
            xcr = XCUR(jmaxstep)*XCUR(jmaxstep)
            xcd = XNEW(jmaxstep)*XNEW(jmaxstep)
          else
            xcr = XCUR(jmaxstep)
            xcd = XNEW(jmaxstep)
          endif
          WRITE(IOUNIT,4)xcr,xcd
 3        FORMAT(2x,'max gradient due to ',a,' grid point ',i5)
 4        FORMAT(2x,'old value ',1pe14.7,' new value ',1pe14.7)
 5        FORMAT(2x,'max residual due to ',a,' equation,',
     .                                       ' at grid point',i5)
 6        FORMAT(2x,'Max change in variable ',a,
     .           ' at grid point',i5)

      ENDIF
c
c
c     some diagnostics. Note that wrt_nwt =1 will write out the results at 
c     every iteration. This will be excessive unless the number of time steps
c     is small (eq nmax = 1)
c     if wrt_nwt = x , where x > 1 , then we save the last 
c     wrt_nwt timesteps and
c     rewrite the file evry time that number is exceeded.
c     note that each of the time steps may consist of more than 
c     one iteration and all of the ierations taken during each time step
c     are written out.
      print *,'wrt_nwt =',wrt_nwt 

      IF  (wrt_nwt >0)  THEN
           inquire (file = 'nwt.test', iostat = iostat,
     .                          exist = exists, opened = opened)
           if(wrt_nwt > 1 .and. tipts2 .gt. wrt_nwt)then
              close(unit = io_nwt)
              tipts2 = 0              !count is incremented in tport 
              opened = .false.
           endif
             IF ( .not. opened) THEN
                 io_nwt = 92
                 call getioun(io_nwt,io_nwt)
c                call giveupus(io_nwt)
                 open (unit = io_nwt, file = 'nwt.test',
     .                          status = 'UNKNOWN', iostat = iostat)
                 write(io_nwt,FMT = '(2x,"type 1 nwt file")')
             ENDIF
                 tipts = tipts +1        !uses this  when wrt_nwt = 1
                 tipts2 = tipts2 + 1     !for  wrt_nwt = > 1
                 write(io_nwt,
     .           FMT='("  nj, tot_iters,time,ssq  = ",i5,2x,i5,
     .              2(2x, 1pe16.8),i3)') nj,tipts,time, ssqr,wrt_nwt
                 WRITE(IO_NWT,16)ITNCOUNT,GRADMAX,GRADTOL,
     &                          STEPMAX,STEPTOL,RESIDMAX,
     &                          SSQR,FROZEN,FREEZE_TYPE
 16              FORMAT(
     &             2X,I3,6(1PE10.3,1X),'-',I2,X,I2)
                 do j =1,n
                       dumy = xnew(j)
                       k = grid_point(j)
                       if(nameu_gp(j)(1:2) .ne. 'wr')dumy =dumy*dumy
                       if(nameu_gp(j)(1:2) .eq. 'te') then
                           dumy2 = d(nion+1,nion+1,k)
                           dumy3 = xke_ifs(k)
                       endif
                      if(nameu_gp(j)(1:2) .eq. 'ti')then 
                           dumy2 = d(nion+2,nion+2,k)
                           dumy3 = xki_ifs(k)
                      endif
                      if(nameu_gp(j)(1:2) .eq. 'bp')then 
                           dumy2 = d(nion+3,nion+3,k)
                           dumy3 = 0.0
                      endif
                      if(nameu_gp(j)(1:2) .eq. 'wr')then 
                           dumy2 = d(nion+4,nion+4,k)
                           dumy3 = xkang_ifs(k)*dmassden(k)*
     .                                            rmajor*rmajor
                      endif
                      if(nameu_gp(j)(1:2) .eq. 'ni')then  !'ni' is currently not a vaild name 
                           dumy2 = d(1,1,k)
                           dumy3 = 0.0
                      endif
                      write(io_nwt,10)dumy,FVP(j)/typf(j),dumy2,dumy3,
     .                                            j,k,nameu_gp(j)(1:8)
 10                   format(4(2x,1pe16.8),2x,i4,2x,i3,2x,a)
                  enddo
                  write(io_nwt,FMT='("end_pass")')
      ENDIF
C
C
      RETURN
      END

C
C
C
      SUBROUTINE CHOLDECOMP(N,H,POSDEF,EPSILOND,L,ADDMAX)
C----------------------------------------------------------------------
C  CALCULATE THE PERTURBED CHOLESKY DECOMPOSITION OF THE SYMMETRIC
C  MATRIX H. H MUST BE IN UPPER TRIANGULAR SYMMETRIC STORAGE MODE.
C  (IE H(I,J),FOR J  .GE. I,IS ACCESSED AS H(K),WHERE
C       K=((I-1)(2*N-I)+2*J)/2
C   MATRIX H SHOULD BE POSITIVE DEFINITE. IF IT IS NOT THEN THE
C   DIAGONAL ELEMENS OF H ARE (IMPLICITELY) MODIFIED A MINIMAL AMOUNT UNTIL
C   H IS POSITIVE DEFINITE. NOTE THAT H IS ACTUALLY NOT CHANGED. INSTEAD
C   IT IS L THAT IS CHANGED TO ACCOUNT FOR THE NON POSITIVE DEFINITENESS
C   OF H. ON OUTPUT H WILL BE THE SAME AS ON INPUT.
C  THE LOWER TRIANGULAR MATRIX L,OF THE CHOLESKY DECOMPOSITION OF
C  H,(IE H=L*LT) IS RETURNED .
C  THIS SUBROUTINE FORMS L BY FIRST GETTING THE DIAGONAL ELEMENT
C  L(I,I) AND THEN CALCULATING THE ELEMENTS OF COLUMN I BELOW THIS
C  DIAGONAL ELEMENT (IE L(J,I),J=I+1,...N). 
C  INPUT
C  N            SQUARE SIZE OF MATRIX H
C  H            MATRIX WHICH MAY OR MAY NOT BE POSTIVE DEFINITE.
C               IF IT IS NOT THEN A DIAGONAL MATRIX D WILL BE ADDED
C               TO H TO MAKE IT POSITIVE DEFINITE. THE DIAGONAL
C               ELEMENTS OF D WILL NOT ALL BE THE SAME. THE LARGEST
C               ELEMENT OF D IS RETURNED AS ADDMAX.
C               NOTE THAT H MUST BE IN UPPER TRIANGULAR SYMMETRIC
C               STORAGE MODE AS EXPLAINED ABOVE. THE ADDITION OF
C               THE DIAGONAL ELEMENTS IS DONE IMPLICITELY,SO THAT H
C               IS ACTUALLY NEVER CHANGED.
C  POSDEF        A TOLERANCE PARAMETER.SET POSDEF=0.0 IF H IS KNOWN
C               TO BE POSITIVE DEFINITE A PRIORI.
C---OUTPUT
C  L(I,j)         LOWER TRIANGULAR CHOLESKY FACTOR OF H IN
C               Full STORAGE MODE
C  ADDMAX       THE MAXIMUM DIAGONAL ELEMENT OF MATRIX D
C               THAT WAS ADDED TO H.
C
C
C------------------------------------------------------------------HSJ
C
C
      IMPLICIT NONE
      INTEGER N,N2,I,J,IHJJ,IJI,III,K
      REAL*8 MINMU,MINMU2,POSDEF,ADDMAX,
     &     HMAX,SQRTEPS,MINLJJ,SQRT,SQRT4
      REAL*8 SUMD,H(*),L(N,N),EPSILOND
      DATA SQRTEPS,SQRT4 /0.0D0,0.0D0/
C
C
      N2=2*N
      IF(SQRTEPS .EQ. 0.0D0)SQRTEPS=SQRT(EPSILOND)
      IF(SQRT4 .EQ. 0.0D0)SQRT4=SQRT(SQRTEPS)
      MINMU2=0.0
      MINMU=SQRT4*POSDEF
      IF(POSDEF .EQ. 0.0)THEN  !H is known to be postive definite
          HMAX=-1.
          DO I=1,N
              III=(I*(N2-I+3))/2-N
              HMAX=MAX(HMAX,ABS(H(III)))
          ENDDO
          POSDEF=SQRT(HMAX)
          MINMU2=SQRTEPS*POSDEF
      ENDIF
      ADDMAX=0.0D0

c        !$OMP PARALLEL DO PRIVATE(J,IHJJ,SUMD,MINLJJ,IJI,ADDMAX)
      DO J=1,N
C         FIRST GET THE DIAGONAL ELEMENT L(ILJJ). THE NORMALIZATION
C         OF THIS ELEMENT IS DONE AFTER WE DETERMINE IF A PERTURBATION
C         IS REQUIRED.
          SUMD=0.0D0
          DO I=1,J-1
              SUMD=SUMD+L(j,i)**2
          ENDDO
          IHJJ=(J*(N2-J+3))/2-N
          L(J,J)=H(IHJJ)-SUMD
C         NEXT GET THE REMAINDER OF COLUMN J BELOW THE DIAGONAL.
C         MINLJJ KEEPS TRACK OF THE LARGEST (IN ABSOLUTE VALUE)
C         OFFDIAGONAL ELEMENT ENCOUNTERED IN THIS COLUMN:
C         THE NORMALIZATION OF L(I,J) IS DONE AFTER L(J,J) IS
C         NORMALIZED.
          MINLJJ=0.0
          DO I=J+1,N
              SUMD=0.0D0
              DO K=1,J-1
                  SUMD=SUMD+L(i,k)*L(j,k)
              ENDDO
              IJI=((J-1)*(N2-J)+2*I)/2
              L(i,j)=H(IJI)-SUMD
              MINLJJ=MAX(MINLJJ,ABS(L(I,J)))
          ENDDO
          MINLJJ=MAX(MINLJJ/POSDEF,MINMU)
C         NOW FINISH THE CALCULATION OF THE DIAGONAL TERM:
          IF(L(j,j) .GT. MINLJJ**2)THEN   !STANDARD CASE
              L(j,j)=SQRT(L(j,j))
          ELSE                            !PERTURBED CASE
              IF(MINLJJ .LT. MINMU2)MINLJJ=MINMU2
              ADDMAX=MAX(ADDMAX,MINLJJ**2-L(j,j))
              L(j,j)=MINLJJ
          ENDIF
C         FINISH THE CALCULATION OF THE REMAINDER OF COLUMN J:
          DO I=J+1,N
              L(I,J)=L(I,J)/L(J,J)
          ENDDO
      ENDDO
C
C
      RETURN
      END
C
C
C
      SUBROUTINE CHOLSOLVE(N,GRADIENT,L,STEP)
C----------------------------------------------------------------------
C   SOLVE L*LT*STEP=-GRADIENT
C  WHERE LT (IE L TRANSPOSE) IS THE UPPER TRIANGLE AND L IS THE
C  LOWER TRIANGLE OF THE CHOLESKY DECOMPOSITION OF SOME POSITIVE
C  DEFINITE MATRIX.
C  THE SOLUTION IS DONE IN TWO STEPS. DEFINE Y=-LT*STEP.
C  THEN A) SOLVE L*Y=GRADIENT FOR Y
C  AND  B) SOLVE LT*STEP=-Y FOR STEP
C  (DOING IT THIS WAY AVOIDS MODIFYING THE GRADIENT AND REQUIRES NO
C   ADDITIONAL STORAGE).
C---INPUT
C   GRADIENT(I)    I=1,2..N,THE RHS
C   L(I,J)         CHOLESKY decomposed matrix
C   N              SIZE OF PROBLEM
C---OUTPUT
C  STEP(I)         I=1,2..N THE SOLUTION
C                  NOTE THAT WE SOLVE SET OF EQUATIONS FOR THE
C                  NEGATIVE OF THE RHS,BUT DO NOT CHANGE THE SIGN
C                  OF THE RHS VECTOR.
C-----------------------------------------------------------------HSJ
C
C
      IMPLICIT NONE
      INTEGER J,I,N
      REAL *8 SUMD,L(N,N),GRADIENT(*),STEP(*)

C
C
C---SOLVE L*Y=GRADIENT,USE STEP TO STORE Y:
      Call Lsolve(n,gradient,L,step)

C
C---SOLVE LT*STEP=-Y;
      STEP(N)=-STEP(N)/L(n,n)            !NOTE THE - SIGN HERE


      DO I=N-1,1,-1
          SUMD=0.0
          DO J=I+1,N
              SUMD=SUMD+L(j,i)*STEP(J)
          ENDDO
          STEP(I)=-(SUMD+STEP(I))/L(i,i)  !NOTE REVERSAL OF SIGN
       ENDDO
C
C
      RETURN
      END
      SUBROUTINE condest(n,M,M2,est)
c--------------------------------------------------------------------------
cEstimate L1 condition number of upper triangular matrix R
cstored in upper triangle of M and diagonal M2
c------------------------------------------------------------------HSJ
      IMPLICIT none
      INTEGER 
     .       n
      REAL *8
     .       M(n,n)
c LOCAL:
       INTEGER j,i
       REAL*8
     .       est,temp,tempm,xnorm,x(n),p(n),pm(n),
     .       xp,xm

!OUTPUT:
      REAL*8
     .     M2(n)


      est =ABS(M2(1))

      DO j =2,n
         temp =0.0d0
         DO i =1,j-1
            temp =temp +ABS(M(i,j))
         ENDDO
         temp =temp +ABS(M2(j))
         est = MAX(temp,est)
      ENDDO
      x(1)=1.D0/M2(1)
      DO i =2,n
         p(i) = M(1,i)*x(1)
      ENDDO
      DO j=2,n
         xp =(1.d0-p(j))/M2(j)
         xm = -(1.d0+p(j))/M2(j)
         temp = ABS(xp)
         tempm = ABS(xm)
         DO i = j+1,n
            pm(i) = p(i) +M(j,i)*xm
            tempm = tempm + ABS(pm(i)/M2(i))
            p(i) =p(i) +M(j,i)*xp
            temp =temp +ABS(p(i)/M2(i))
         ENDDO
         IF(temp .ge. tempm)THEN
            x(j) = xp  ! ej = 1
         ELSE
            x(j) = xm
            DO i = j+1,n
               p(i) = pm(i)
            ENDDO
         ENDIF
      ENDDO
      xnorm =0.0D0
      DO j =1,n
         xnorm =xnorm +ABS(x(j))
      ENDDO
      est =est /xnorm
      CALL rsolve(n,M,M2,x)
      xnorm =0.0D0
      DO j =1,n
         xnorm =xnorm +ABS(x(j))
      ENDDO

      est=est*xnorm

      RETURN
      END
C
          SUBROUTINE DOGSTEP(N,GRADIENT,L,STEP,SCALE,NOSCALE,STEPLENGTH,
     &                       STEPLENGTH_MAX,DELTA,FIRSTDOG,CAUCHYLENGTH,
     &                       ETA,SSD,V,S,NEWTON_STEP,IOUNIT)
C-----------------------------------------------------------------------
C---HERE WE FIND AN APPROXIMATE SOLUTION TO
C    MINMIZE ( GT*S+0.5*ST*L*LT*S) SUBJECT TO THE CONSTRAINT THAT
C    THE (SCALED) STEP LENGTH IS LESS THAN DELTA.
C    THE OUTPUT STEP S(I),I=1,2..N WILL HAVE LENGTH DELTA IF THE
C    FULL NEWTON STEP COULD NOT BE TAKEN (BECAUSE IT IS OUTSIDE
C    THE TRUST REGION). IF THE FULL NEWTON STEP WAS TAKEN,INDICATED
C    BY NEWTON_STEP =TRUE,THEN THE OUTPUT STEP S(I) IS IDENTICAL
C    TO THE INPUT STEP(I).
C---INPUT
C  N                  NUMBER OF ELEMENTS IN STEP,SCALE,GRADIENT
C  GRADIENT(I)        GRADIENT OF OBJECTIVE FUNCTION AT XCUR
C  L(N,N)               CHOLESKY FACTOR OF HESSIAN
C  STEP(I)            I=1,2..N  THE CURRENT STEP VECTOR
C  SCALE(I)           DIAGONAL SCALING MATRIX
C  NOSCALE            IF .TRUE. DONT USE SCALE,
C                     IF .FALSE. SCALE CALCULATIONS USING SCALE
C  STEPLENGTH         MAGNITUDE OF STEP
C  STEPLENGTH_MAX     MAX ALLOWABLE STEP SIZE
C
C
C---INPUT AND OUTPUT:
C  FIRSTDOG           IS SET TO FALSE TO INDICATE THAT THE NEXT
C                     CALL TO THIS SUBROUTINE DOES NOT HAVE TO
C                     CALCULATE THE DOUBLE DOGLEG CURVE
C  DELTA              THE CURRENT TRUST REGION RADIUS
C  CAUCHYLENGTH       LENGTH OF THE CAUCHY STEP
C  ETA
C  SSD(I)             I=1,2..N,THE STEEPEST DESCENT STEP
C  V(I)
C
C
C---OUTPUT:
C  NEWTON_STEP        SET TO TRUE IF FULL NEWTON STEP WAS TAKEN
C                     OTHERWISE SET TO FALSE.
C  S(I)               I=1,2..N  THE NEW STEP VECTOR
C---------------------------------------------------------------------HSJ
C}}}
      USE copy
      IMPLICIT NONE
      INTEGER N,I,J,IOUNIT
      REAL*8
     &     STEPLENGTH,STEPLENGTH_MAX,DELTA,CAUCHYLENGTH,ETA,
     &     DUMY,DUMY1,DUMY2,LAMBDA,ALPHA,BETA,SCALE(*),ENORM, 
     &      L(N,N),GRADIENT(*),STEP(*),SSD(*),V(*),S(*),SDOT12,
     &     DSDOT
      LOGICAL 
     &     NOSCALE,NEWTON_STEP,FIRSTDOG
C
c      print *,'Newton length,delta =',STEPLENGTH,DELTA

      IF(STEPLENGTH .LE. DELTA)THEN  ! on first entry delta = -1 
          NEWTON_STEP=.TRUE.         !so netwon step is assumed too long
          CALL DCOPY(N,STEP,1,S,1)   !and we calculate delta below
          DELTA=STEPLENGTH
          IF(IOUNIT .NE. 0)WRITE(IOUNIT,*)'NEWTON STEP TAKEN'
      ELSE
C---NEWTON STEP TOO LONG,GET S ON DOUBLE DOGLEG CURVE:
          NEWTON_STEP=.FALSE.
          IF(FIRSTDOG)THEN
C---CALCULATE THE DOUBLE DOGLEG CURVE ON FIRST TRY
              FIRSTDOG=.FALSE.  !SET FOR NEXT CALL TO THIS ROUTINE
              IF(.NOT. NOSCALE)THEN
                  DO J=1,N
                      S(J)=GRADIENT(J)/SCALE(J)
                  ENDDO
              ELSE
                  CALL DCOPY(N,GRADIENT,1,S,1)
              ENDIF
chsj              CALL EUCLIDNORM(S,SCALE,N,.TRUE.,ALPHA)
              ALPHA = ENORM(n,s)
              ALPHA=ALPHA*ALPHA
              BETA=0.0
              DO I=1,N
                  DUMY2=0.0
                  DO J=I,N
                      DUMY=1.0
                      IF(.NOT. NOSCALE)DUMY=SCALE(J)
                      DUMY2=DUMY2+L(j,i)*GRADIENT(J)/(DUMY*DUMY)
                  ENDDO
                  BETA=BETA+DUMY2*DUMY2
              ENDDO
              DUMY=-ALPHA/BETA
              DO J=1,N
                  SSD(J)=DUMY*S(J)   ! SSD is Cauchy step 
              ENDDO
              CAUCHYLENGTH=ALPHA*SQRT(ALPHA)/BETA
              DUMY=SDOT12(GRADIENT,STEP,N)
c              DUMY=DSDOT(GRADIENT,STEP,N)
              ETA=0.2+0.8*ALPHA*ALPHA/(BETA*ABS(DUMY))
c          !$OMP PARALLEL DO PRIVATE(DUMY)
              DO J=1,N
                  DUMY=ETA
                  IF( .NOT. NOSCALE)DUMY=SCALE(J)*DUMY
                  V(J)=DUMY*STEP(J)-SSD(J)
              ENDDO
              IF(DELTA .EQ. -1.d0)
     .                   DELTA=MIN(CAUCHYLENGTH,STEPLENGTH_MAX)
          ENDIF
          IF(ETA*STEPLENGTH .LE. DELTA)THEN
C---PARTIAL LENGTH STEP IN NEWTON DIRECTION:
              DUMY=DELTA/STEPLENGTH
              DO J=1,N
                  S(J)=DUMY*STEP(J)
              ENDDO
               IF(IOUNIT .NE. 0)WRITE(IOUNIT,*)'PARTL NEWTON STEP TAKEN'
          ELSE IF(CAUCHYLENGTH .GE. DELTA)THEN
C---TAKE STEP IN STEEPEST DESCENT DIRECTION:
              DUMY=DELTA/CAUCHYLENGTH
              DO J=1,N
                  DUMY1=1.0
                  IF(.NOT. NOSCALE)DUMY1=1./SCALE(J)
                  S(J)=DUMY*SSD(J)*DUMY1
              ENDDO
               IF(IOUNIT .NE. 0)WRITE(IOUNIT,*)'STEPST DESCT STEP TAKEN'
          ELSE
C---GET THE CONVEX COMBINATION THAT HAS SCALED LENGTH DELTA:
              DUMY=SDOT12(V,SSD,N)
c              DUMY=DSDOT(V,SSD,N)
              DUMY1=SDOT12(V,V,N)
c              DUMY1=DSDOT(V,V,N)
              LAMBDA=(-DUMY+SQRT(DUMY**2-DUMY1*
     &                             (CAUCHYLENGTH**2-DELTA**2)))/DUMY1
              DO J=1,N
                  DUMY=1.
                  IF(.NOT. NOSCALE)DUMY=1./SCALE(J)
                  S(J)=DUMY*(SSD(J)+LAMBDA*V(J))
              ENDDO
                  IF(IOUNIT .NE. 0)WRITE(IOUNIT,*)'DOGLEG STEP TAKEN'
          ENDIF
      ENDIF
C
C
      RETURN
      END





      real*8 function dpmpar(i) !left here for reference , not used otherwise HSJ
      integer i
c     **********
c
c     Function dpmpar
c
c     This function provides real*8 machine parameters
c     when the appropriate set of data statements is activated (by
c     removing the c from column 1) and all other data statements are
c     rendered inactive. Most of the parameter values were obtained
c     from the corresponding Bell Laboratories Port Library function.
c
c     The function statement is
c
c       real*8 function dpmpar(i)
c
c     where
c
c       i is an integer input variable set to 1, 2, or 3 which
c         selects the desired machine parameter. If the machine has
c         t base b digits and its smallest and largest exponents are
c         emin and emax, respectively, then these parameters are
c
c         dpmpar(1) = b**(1 - t), the machine precision,
c
c         dpmpar(2) = b**(emin - 1), the smallest magnitude,
c
c         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c     Argonne National Laboratory. MINPACK Project. November 1996.
c     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'
c
c     **********
      integer mcheps(4)
      integer minmag(4)
      integer maxmag(4)
      real*8 dmach(3)
      equivalence (dmach(1),mcheps(1))
      equivalence (dmach(2),minmag(1))
      equivalence (dmach(3),maxmag(1))
c
c     Machine constants for the IBM 360/370 series,
c     the Amdahl 470/V6, the ICL 2900, the Itel AS/6,
c     the Xerox Sigma 5/7/9 and the Sel systems 85/86.
c
c     data mcheps(1),mcheps(2) / z34100000, z00000000 /
c     data minmag(1),minmag(2) / z00100000, z00000000 /
c     data maxmag(1),maxmag(2) / z7fffffff, zffffffff /
c
c     Machine constants for the Honeywell 600/6000 series.
c
c     data mcheps(1),mcheps(2) / o606400000000, o000000000000 /
c     data minmag(1),minmag(2) / o402400000000, o000000000000 /
c     data maxmag(1),maxmag(2) / o376777777777, o777777777777 /
c
c     Machine constants for the CDC 6000/7000 series.
c
c     data mcheps(1) / 15614000000000000000b /
c     data mcheps(2) / 15010000000000000000b /
c
c     data minmag(1) / 00604000000000000000b /
c     data minmag(2) / 00000000000000000000b /
c
c     data maxmag(1) / 37767777777777777777b /
c     data maxmag(2) / 37167777777777777777b /
c
c     Machine constants for the PDP-10 (KA processor).
c
c     data mcheps(1),mcheps(2) / "114400000000, "000000000000 /
c     data minmag(1),minmag(2) / "033400000000, "000000000000 /
c     data maxmag(1),maxmag(2) / "377777777777, "344777777777 /
c
c     Machine constants for the PDP-10 (KI processor).
c
c     data mcheps(1),mcheps(2) / "104400000000, "000000000000 /
c     data minmag(1),minmag(2) / "000400000000, "000000000000 /
c     data maxmag(1),maxmag(2) / "377777777777, "377777777777 /
c
c     Machine constants for the PDP-11. 
c
c     data mcheps(1),mcheps(2) /   9472,      0 /
c     data mcheps(3),mcheps(4) /      0,      0 /
c
c     data minmag(1),minmag(2) /    128,      0 /
c     data minmag(3),minmag(4) /      0,      0 /
c
c     data maxmag(1),maxmag(2) /  32767,     -1 /
c     data maxmag(3),maxmag(4) /     -1,     -1 /
c
c     Machine constants for the Burroughs 6700/7700 systems.
c
c     data mcheps(1) / o1451000000000000 /
c     data mcheps(2) / o0000000000000000 /
c
c     data minmag(1) / o1771000000000000 /
c     data minmag(2) / o7770000000000000 /
c
c     data maxmag(1) / o0777777777777777 /
c     data maxmag(2) / o7777777777777777 /
c
c     Machine constants for the Burroughs 5700 system.
c
c     data mcheps(1) / o1451000000000000 /
c     data mcheps(2) / o0000000000000000 /
c
c     data minmag(1) / o1771000000000000 /
c     data minmag(2) / o0000000000000000 /
c
c     data maxmag(1) / o0777777777777777 /
c     data maxmag(2) / o0007777777777777 /
c
c     Machine constants for the Burroughs 1700 system.
c
c     data mcheps(1) / zcc6800000 /
c     data mcheps(2) / z000000000 /
c
c     data minmag(1) / zc00800000 /
c     data minmag(2) / z000000000 /
c
c     data maxmag(1) / zdffffffff /
c     data maxmag(2) / zfffffffff /
c
c     Machine constants for the Univac 1100 series.
c
c     data mcheps(1),mcheps(2) / o170640000000, o000000000000 /
c     data minmag(1),minmag(2) / o000040000000, o000000000000 /
c     data maxmag(1),maxmag(2) / o377777777777, o777777777777 /
c
c     Machine constants for the Data General Eclipse S/200.
c
c     Note - it may be appropriate to include the following card -
c     static dmach(3)
c
c     data minmag/20k,3*0/,maxmag/77777k,3*177777k/
c     data mcheps/32020k,3*0/
c
c     Machine constants for the Harris 220.
c
c     data mcheps(1),mcheps(2) / '20000000, '00000334 /
c     data minmag(1),minmag(2) / '20000000, '00000201 /
c     data maxmag(1),maxmag(2) / '37777777, '37777577 /
c
c     Machine constants for the Cray-1.
c
c     data mcheps(1) / 0376424000000000000000b /
c     data mcheps(2) / 0000000000000000000000b /
c
c     data minmag(1) / 0200034000000000000000b /
c     data minmag(2) / 0000000000000000000000b /
c
c     data maxmag(1) / 0577777777777777777777b /
c     data maxmag(2) / 0000007777777777777776b /
c
c     Machine constants for the Prime 400.
c
c     data mcheps(1),mcheps(2) / :10000000000, :00000000123 /
c     data minmag(1),minmag(2) / :10000000000, :00000100000 /
c     data maxmag(1),maxmag(2) / :17777777777, :37777677776 /
c
c     Machine constants for the VAX-11.
c
c     data mcheps(1),mcheps(2) /   9472,  0 /
c     data minmag(1),minmag(2) /    128,  0 /
c     data maxmag(1),maxmag(2) / -32769, -1 /
c
c     Machine constants for IEEE machines.
c
      data dmach(1) /2.22044604926d-16/
      data dmach(2) /2.22507385852d-308/
      data dmach(3) /1.79769313485d+308/
c
      dpmpar = dmach(i)
      return
c
c     Last card of function dpmpar.
c
      end





C
      SUBROUTINE DRIVE_DOGLEG(N,XCUR,FCUR,FVEC,GRADIENT,L,
     &                       STEP,SCALE,NOSCALE,Sf,FVp,
     &                       STEPLENGTH_MAX,STEP_TOL,
     &                       DELTA,RETCODE,XNEW,FNEW,
     &                       MAX_TAKEN,IOUNIT)
C----------------------------------------------------------------------
C---FIND XNEW  ALONG THE DOUBLE DOGLEG CURVE SUCH THAT
C---FNEW<FCUR+ALPHA*GTRANSPOSE*P,WHERE FNEW AND FCUR ARE THE
C---OBJECTIVE FUNCTION EVALUATED AT XNEW AND XCUR RESPECTIVELY,
C---AND GTRANSPOSE IS THE GRADIENT OF OBJECTIVE FUNCTION TRANSPOSE AT XCUR
C---THE SOLUTION VECTOR P IS CONSTRAINED TO HAVE (SCALED) STEP
C---LENGTH DELTA. DELTA WILL BE INCREASED OR DECREASED AS APPROPRIATE.
C---INPUT
C  N                  NUMBER OF ELEMENTS IN XCUR,P,SCALE,XNEW
C  XCUR(I)            I=1,2...N CURRENT SOLUTION POINT
C  FCUR               VALUE OF OBJECTIVE FUNCTION AT XCUR
C  FVEC                 SUBROUTINE NAME WHICH CALCULATES FNEW
C  GRADIENT(I)        GRADIENT OF OBJECTIVE FUNCTION AT XCUR
C  STEP(I)         DESCENT VECTOR FROM XCUR
C  L(N,N)               CHOLESKY FACTOR OF HESSIAN
C  SCALE(I)           DIAGONAL SCALING MATRIX
C  NOSCALE            IF .TRUE. DONT USE SCALE,
C                     IF .FALSE. SCALE CALCULATIONS USING SCALE
C  STEPLENGTH_MAX     MAX ALLOWABLE STEP SIZE
C  STEP_TOL           THE STEP TOLERANCE
C               THE FOLLOWING FOUR STORAGE VECTORS ARE REQUIRED
C               TO HOLD INTERMEDIATE RESULTS:
C  S(I)               I=1,2..N
C  SSD(I)             I=1,2..N
C  XPREV(I)           I=1,2..N
C  V(I)               I=1,2..N
c  delta         trust region radius

C
C
C---INPUT AND OUTPUT:
C  DELTA                ON INPUT THE CURRENT TRUST REGION RADIUS
C                       ON OUTPUT THE NEW TRUST REGION RADIUS
C
C 
C---OUTPUT:
C  RETCODE              INTEGER VARIABLE WITH MEANING AS FOLLOWS
C                       RETCODE=0   VALID XNEW FOUND
C                               1   FAILED TO FIND A SATISFACTORY XNEW,
C                                   SUFFICIENTLY DISTINC FROM XCUR
C                                   XCUR IS COPIED TO XNEW. THIS
C                                   ZERO STEPLENGTH TO BE REPORTED IN
C                                   CHECK_CONVERGENCE.
C  XNEW(I)              I=1,2..N THE NEW SOLUTION POINT
C                       XNEW(I)=XCUR(I)+P(I) ON OUTPUT
C  FNEW                 VALUE OF OBJECTIVE FUNCTION AT XNEW
C  MAX_TAKEN            LOGICAL VARIABLE =.TRUE. IF MAX STEP
C                       WAS TAKEN (THIS ALLOWS US TO KEEP TRACK
C                       OF HOW MANY SUCCESSIVE STEPS OF MAXIMUM LENGTH
C                       ARE TAKEN. IF THERE ARE TOO MANY THEN THE
C                       MAX STEP LENGTH,STEPLENGTH_MAX,SHOULD BE
C                       INCREASED FOR BETTER EFFICIENCY.
C
C
C------------------------------------------------------------------HSJ
      IMPLICIT NONE
      INTEGER N,RETCODE,STEP_TYPE,IOUNIT
      REAL*8 SCALE(*),
     &     STEPLENGTH_MAX,STEP_TOL,STEPLENGTH,
     &     DELTA,CAUCHYLENGTH,ETA
      REAL *8 L(n,n),XCUR(*),GRADIENT(*),STEP(*),
     &         XNEW(*),HESSIAN(1),Sf(n),FVp(n),
     &         FCUR,FNEW,FPREV,delta_prev
      REAL*8    !these are automatic arrays (created and destroyed on each 
                !entry and exit:
     &         S(N),SSD(N),XPREV(N),V(N)
      LOGICAL MAX_TAKEN,NOSCALE,FIRSTDOG,NEWTON_STEP
      EXTERNAL FVEC
C
C
      FIRSTDOG=.TRUE.
      STEP_TYPE=2  !SELECT DOGSTEP
      RETCODE=4
C---GET THE (SCALED) LENGTH OF THE DESCENT VECTOR:
      CALL EUCLIDNORM(STEP,SCALE,N,NOSCALE,STEPLENGTH)

      DO WHILE (RETCODE .GE. 2)
C---CALCULATE THE DOUBLE DOGLEG CURVE (ON FIRST CALL ONLY) AND
C---GET THE NEW STEP VECTOR S,SUCH THAT S HAS LENGTH DELTA.
          CALL DOGSTEP(N,GRADIENT,L,STEP,SCALE,NOSCALE,STEPLENGTH,
     &                       STEPLENGTH_MAX,DELTA,FIRSTDOG,CAUCHYLENGTH,
     &                       ETA,SSD,V,S,NEWTON_STEP,IOUNIT)
C
C
          !note step_type = 2 here. Hence the HESSIAN
          !is in fact never used in the trust_region routine
          !(trust_region uses HESSIAn only if step_type = 1)
          !thus hessian does not have to be set when passed on
          !from here
          CALL TRUST_REGION(N,XCUR,FCUR,FVEC,GRADIENT,L,S,SCALE,
     &         NOSCALE,Sf,FVp,
     &         NEWTON_STEP,STEPLENGTH_MAX,STEP_TOL,STEP_TYPE,
     &         HESSIAN,DELTA,RETCODE,XPREV,FPREV,XNEW,FNEW,MAX_TAKEN)

      ENDDO
C
C
      RETURN
      END


      SUBROUTINE EUCLIDNORM(VEC,SCALE,N,NOSCALE,LENGTH)
C---------------------------------------------------------------------
C---CALCULATE THE EUCLIDEAN NORM OF VECTOR VEC.
C---IF NOSCALE= .FALSE. THEN SCALE VEC BY SCALE BEFORE
C---DOING THE CALCULATION. THIS LENGTH DOES NOT HAVE TO
C---BE VERY PRECISE SO WE DONT USE SDOT12 HERE.
C------------------------------------------------------------------HSJ
C
C
      IMPLICIT NONE
      LOGICAL NOSCALE
      INTEGER N,J
      REAL *8 VEC(*),SUM
      REAL*8 LENGTH,SCALE(*)
C
C
 
      SUM=0.0D0
      IF(.NOT. NOSCALE)THEN
c       !$OMP PARALLEL DO
c       !$OMP+ reduction(+:sum)
          DO J=1,N
              SUM=SUM+(SCALE(J)*VEC(J))**2
          ENDDO
      ELSE
c       !$OMP PARALLEL DO
c       !$OMP+ reduction(+:sum)
          DO J=1,N
              SUM=SUM+VEC(J)*VEC(J)
          ENDDO
      ENDIF
      LENGTH=SQRT(SUM)
C
C
      RETURN
      END
      SUBROUTINE fdjac(n,xc,fc,FVEC,Sx,eta,jacf)
c------------------------------------------------------------------------
c--- evaluate finite difference Jacobian


c------------------------------------------------------------------------
      USE aid_newton


      IMPLICIT none
c  INPUT:
      INTEGER n
      REAL *8 xc(n),fc(n),Sx(n),eta
C  LOCAL:
      INTEGER
     .        j,i

      REAL*8
     .       sqrteta,stepsizej,tempj,Fj(n)
      LOGICAL TESTING

c  OUTPUT:
      REAL*8 jacf(n,n)
      REAL*8 aa,bb,cc,xx,ss,gg,terma,termb,termc
      EXTERNAL FVEC
      common /temporary / aa,bb,cc,xx(51),ss,gg,terma,termb,termc

      print *,'freeze_xte,ti =',freeze_xte,freeze_xti
      sqrteta =SQRT(eta)
      DO j = 1,n !perturb each of the n variables in turn to generate a 
                 !single column of the Jacobian at a time:
            go to 10
            if(freeze_xte .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'te')
     .            go to 60
            if(freeze_xti .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'ti')
     .            go to 60
            if(freeze_xrbp .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'rbp')
     .            go to 60
            if(freeze_xwr .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'wr')
     .            go to 60
            if(freeze_xni .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'ni')
     .            go to 60
         !calculate column j of JACF
 10      stepsizej = sqrteta*MAX(ABS(xc(j)),1./Sx(j))*SIGN(1.D0,xc(j))
         tempj = xc(j)
         xc(j) = xc(j) + stepsizej
         stepsizej = xc(j)-tempj   !reduces finite precision error slightly 
         !FVEC = eval_non_lin_eq which evaluates the non
         ! linear equations and obtains the residuals in Fj .
         !eval_non_lin_eq calls solve then subrotuine residual
chsj         CALL FVEC(n,xc,Sf,Fj,ssqrc)
          call FVEC(n,xc,Fj) 
         DO i = 1, n
            jacf(i,j) = (Fj(i) -Fc(i))/stepsizej
         ENDDO
         xc(j) = tempj     !set xc back to original value
 60       CONTINUE
      ENDDO
      testing = .false.
      IF(testing)then
c-----------------------testing only
         print *,'**************** j = 19 ******'
         fj(20) = 0.0
         j=19
         stepsizej = 1.e0*
     .      sqrteta*MAX(ABS(xc(j)),1./Sx(j))*SIGN(1.D0,xc(j))
         tempj = xc(j)
         xc(j) = xc(j) + stepsizej
         call FVEC(n,xc,Fj)
         print*, 'aa,bb,cc ',aa,bb,cc
         print *,'gg,ss,termabc =',gg,ss,terma,termb,termc
         print *,'xx(20-22) =',xx(20),xx(21),xx(22)
         print *,'sqrt(xx(20)),xc(19) =',SQRT(xx(20)),xc(19)
         print *,'fj(20),fc(20)=',fj(20),fc(20)
         xc(j) = tempj
         print *,'**************** j = 20 ******'
         j=20
         stepsizej = sqrteta*MAX(ABS(xc(j)),1./Sx(j))*SIGN(1.D0,xc(j))
         tempj = xc(j)
         xc(j) = xc(j) + stepsizej
         call FVEC(n,xc,Fj)
         print *, 'aa,bb,cc ',aa,bb,cc
         print *,'gg,ss,termabc =',gg,ss,terma,termb,termc
         print *,'xx(20-22) =',xx(20),xx(21),xx(22)
         print *,'qrt(xx(21)),xc(20) =',SQRT(xx(21)),xc(20)
         print *,'fj(20),fc(20)=',fj(20),fc(20)
         xc(j) = tempj
         print *,'**************** j = 21 ******'
         j=21
         stepsizej =1.e0*
     .       sqrteta*MAX(ABS(xc(j)),1./Sx(j))*SIGN(1.D0,xc(j))
         tempj = xc(j)
         xc(j) = xc(j) + stepsizej
         call FVEC(n,xc,Fj)
         print*, 'aa,bb,cc ',aa,bb,cc
         print *,'gg,ss,termabc =',gg,ss,terma,termb,termc
         print *,'xx(20-22) =',xx(20),xx(21),xx(22)
         print *,'sqrt(xx(22)),xc(21) =',SQRT(xx(22)),xc(21)
         print *,'fj(20),fc(20)=',fj(20),fc(20)
         xc(j) = tempj
         print *,'non perturbed:'
         call FVEC(n,xc,Fj)
         print*, 'aa,bb,cc ',aa,bb,cc
         print *,'gg,ss,termabc =',gg,ss,terma,termb,termc
         print *,'xx(20-22) =',xx(20),xx(21),xx(22)
         print *,'fj(20),fc(20)=',fj(20),fc(20)
         xc(j) = tempj
         do j=1,n
           print *,xc(j)**2,xx(j)
         enddo
         call stop('testing',0)
      endif
c------------------------------------
      RETURN
      END


      subroutine fdjac_sparse(n,x,fvx,FVEC,fjac,ml,mu,
     *                        epsfcn,sx)
c-------------------------------------------------------------------
c     NOTE: this routine is 
c           slightly modified by HSJ from the originaL:
c           a) call to FVEC is changed
c           b) arguments and name of this routine are modified
c              meaning of epsfcn is changed. HSJ
c
c     original subroutine fdjac1
c
c     this subroutine computes a forward-difference approximation
c     to the n by n jacobian matrix associated with a specified
c     problem of n functions in n variables. if the jacobian has
c     a banded form, then function evaluations are saved by only
c     approximating the nonzero terms.
c
c     the original subroutine statement is
c
c       subroutine fdjac1(FVEC,n,x,fvx,fjac,ldfjac,iflag,ml,mu,epsfcn,
c                         wa1,wa2)
c
c     where
c
c       FVEC is the name of the user-supplied subroutine which
c         calculates the functions. FVEC must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine fcn(n,x,fvx,iflag) (original)
c         subroutine FVEC(N,X,fvx) (as used here)  HSJ
c         integer n,iflag
c         real*8 x(n),fvx(n),fp
c         ----------
c         calculate the functions at x and
c         return this vector in fvx.
c         return residual sum of squares in fp
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by FVEC unless
c         the user wants to terminate execution of fdjac1.
c         in this case set iflag to a negative integer.
c         new version ignores iflag
c
c       n is a positive integer input variable set to the number
c         of functions and variables.
c
c       x is an input array of length n.
c
c       fvx is an input array of length n which must contain the
c         functions evaluated at x.
c
c       fjac is an output n by n array which contains the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive integer input variable not less than n
c         which specifies the leading dimension of the array fjac.
c         NOTE: modified version uses f90, assumes fjac is exactly
c         size n by n with n changed dynamically for each new problem HSJ
c
c       iflag is an integer variable which can be used to terminate
c         the execution of fdjac1. see description of FVEC.
c
c       ml is a nonnegative integer input variable which specifies
c         the number of subdiagonals within the band of the
c         jacobian matrix. if the jacobian is not banded, set
c         ml to at least n - 1.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c       epsfcn new definition HSJ: 
c
c       mu is a nonnegative integer input variable which specifies
c         the number of superdiagonals within the band of the
c         jacobian matrix. if the jacobian is not banded, set
c         mu to at least n - 1.
c
c       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at
c         least n, then the jacobian is considered dense, and wa2 is
c         not referenced.
c
c     subprograms called
c
c       minpack-supplied ... dpmpar
c
c       fortran-supplied ... dabs,dmax1,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
c---------------------------------------------------------HSJ----------
      USE aid_newton
      integer n,ldfjac,iflag,ml,mu
      real*8 epsfcn
      real*8 x(n),sx(n),fvx(n),fjac(n,n),wa1(n),wa2(n)

      integer i,j,k,msum
      real*8 eps,epsmch,h,temp,zero
      real*8 dpmpar,sqrteta
      data zero /0.0d0/
c     **********

c
c     epsmch is the machine precision.
c
c      epsmch = dpmpar(1) !function call
      iflag = 0 
c
c      eps = dsqrt(dmax1(epsfcn,epsmch))
c       eps = SQRT(MAX(epsfcn,epsmch))
      sqrteta =SQRT(epsfcn)
      msum = ml + mu + 1
      if (msum .lt. n) go to 40
c
c        computation of dense approximate jacobian.
c
         do 20 j = 1, n  !j is variable #
            go to 11
            if(freeze_xte .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'te')
     .            go to 20
            if(freeze_xti .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'ti')
     .            go to 20
            if(freeze_xrbp .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'rbp')
     .            go to 20
            if(freeze_xwr .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'wr')
     .            go to 20
            if(freeze_xni .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'ni')
     .            go to 20
 11         temp = x(j)
c            h = eps*dabs(temp)   !original
c            if (h .eq. zero) h = eps
             h = sqrteta*MAX(ABS(temp),1./Sx(j)) ! *SIGN(1.D0,temp) !HSJ 
            x(j) = temp + h
c            call FVEC(n,x,wa1,iflag) original
            call FVEC(n,x,wa1)
            if (iflag .lt. 0) go to 30
            x(j) = temp
            do 10 i = 1, n   !i is equation #
               fjac(i,j) = (wa1(i) - fvx(i))/h
   10          continue
   20       continue
   30    continue
         go to 110




   40 continue
c
c        computation of banded approximate jacobian.
c        !$OMP PARALLEL DO
         do 90 k = 1, msum
            do 60 j = k, n, msum 
            go to 12
            if(freeze_xte .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'te')
     .            go to 60
            if(freeze_xti .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'ti')
     .            go to 60
            if(freeze_xrbp .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'rbp')
     .            go to 60
            if(freeze_xwr .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'wr')
     .            go to 60
            if(freeze_xni .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'ni')
     .            go to 60
 12         wa2(j) = x(j)
c              h = eps*dabs(wa2(j))
c              if (h .eq. zero) h = eps
               h = sqrteta*MAX(ABS(wa2(j)),1./Sx(j))  ! *SIGN(1.D0,wa2(j)) !HSJ 
               x(j) = wa2(j) + h

   60          continue

c            call FVEC(n,x,wa1,iflag)  original
c             !$OMP CRITICAL
             call FVEC(n,x,wa1)
c              !$OMP END CRITICAL
            if (iflag .lt. 0) go to 100


            do 80 j = k, n, msum
               GO TO 13
               if(freeze_xte .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'te')
     .            go to 80
               if(freeze_xti .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'ti')
     .            go to 80
               if(freeze_xrbp .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'rbp')
     .            go to 80
               if(freeze_xwr .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'wr')
     .            go to 80
               if(freeze_xni .eq. 1 .and. nameu_gp(j)(1:2) .eq. 'ni')
     .            go to 80
 13            x(j) = wa2(j)

c              h = eps*dabs(wa2(j))
c              if (h .eq. zero) h = eps
               h = sqrteta*MAX(ABS(wa2(j)),1./Sx(j)) ! *SIGN(1.D0,wa2(j)) !HSJ 
               do 70 i = 1, n
                  fjac(i,j) = zero
                  if (i .ge. j - mu .and. i .le. j + ml)
     .               fjac(i,j) = (wa1(i) - fvx(i))/h
c                    if (j .eq. 1)print *,'J,I,FJAC',J,I,FJAC(I,J)
   70          continue
   80      continue
   90   continue
  100 continue



  110 continue
      return
c
c     last card of subroutine fdjac1.
c
      end
      Subroutine get_axis_values(denaxis,teaxis,tiaxis,rbpedgem1,
     .                           angraxis,x)
c---------------------------------------------------------------------
c     get current values at magnetic axis
      USE param
      USE numbrs
      USE mesh
      USE flags
      USE tordlrot
      USE nonlin
      implicit  integer (i-n), real*8 (a-h, o-z)


c      include 'param.i'
c      include 'flags.i' !itran
c      include 'mesh.i'
c      include 'numbrs.i' !nk,
c      include 'nonlin.i' !nrmax
c      include 'tordlrot.i' !iangrot
      real *8 x(*)
      nu = 0 
      ik = 0
      denaxis = 0.0
      teaxis =  0.0
      tiaxis =  0.0
      rbpedgem1 = 0.0 
      angraxis = 0.0
chsj  j = 1
      do j=1,nrmax   !central values only
         do k=1,nk
            if(itran(k) .gt. 0)then
c              skip rbp at the origin (its value is fixed at 0.0)
               if(j .eq. 1 .and. k .eq. nk-iangrot)go to 10
               if(j .ge. te_index  .and. k .eq. nk-iangrot-2)go to 10
               if(j .ge. ti_index .and. k .eq. nk-iangrot-1)go to 10
               if(j .ge. rot_index .and. k .eq. nk .and. 
     .                                          iangrot .eq. 1)go to 10
               nu = nu+1
               if(k .le. nion)then
                  ik = ik + 1
chsj              if(ik .eq. 1)denaxis = (xscale(nu)*x(nu))**2
                  if(ik .eq. 1)denaxis = x(nu)**2
               else if(k .eq. nion+1 .and. j .eq. 1 )then
chsj              teaxis = (xscale(nu)*x(nu))**2
                  teaxis = x(nu)**2 
               else if(k .eq. nion+2 .and. j .eq. 1)then
chsj              tiaxis = (xscale(nu)*x(nu))**2
                  tiaxis = x(nu)**2 
               else if(k .eq. nion+3 .and. j .eq. nrmax)then
c                 pick up penultimate value of rbp:
chsj              rbpedgem1 = (xscale(nu)*x(nu))**2 
                  rbpedgem1 = x(nu)**2 
               else if(j .eq. 1)then
chsj              angraxis = xscale(nu)*x(nu)
                  angraxis = x(nu)
               endif
 10         endif
         enddo
         enddo


      RETURN
      END
                                             
                                                                      
                                                                        
      subroutine hookdriver(n,xcur,fcur,FVEC,gradient,L,hessian,       
     &               step,scale,noscale,Sf,FVp,steplength_max,step_tol,
     &               iterat,epsilond,delta,mu,deltaprev,        
     &               phi,phip,retcode,xnew,fnew,max_taken,iounit)      
!---------------------------------------------------------------------  
C INPUT 
c n
c xcur
c fcur
C FVEC
c gradient
c L(n,n)
c step(n)
c scale(n)                



!---the following vectors are used for local storage
!  s(i)          i=1,2..n
!  xprev(i)      i=1,2..n  (f90 automatic arrays here)
!--------------------------------------------------------------------hsj
      USE copy

      implicit none
                              ! in: number of data points
      integer n,iterat,retcode,j,                                       &
     &      i,step_type,iounit
      real *8 xcur(*),L(n,n),hessian(*),gradient(*),step(*),xnew(*),      &
     &      s(n),xprev(n),temp,beta,epsilond,fnew,fcur,fprev
      real*8 scale(*),steplength_max,step_tol,delta,                    &
     &    deltaprev,phi,phip,alpha,sqrt,steplength,                     &
     &    dumy,phipinit,Sf(n),FVp(n),enorm
      real*8 mu
      logical max_taken,noscale,firsthook,newton_step
      ! need a lower limit on steplength; this is arbitrary but a lot
      ! bigger than zero.                                               
      real*8, parameter :: steplength_min=1.0e-5
                                                                        
      external FVEC
                                   ! exception handling routine

                                                     
      retcode=4
                  
      step_type=1  !selects hookstep in trust region
      firsthook=.true.
      call euclidnorm(step,scale,n,noscale,steplength)
c       print *,'steplength in hookdriver',steplength
      !steplength is scaled length
                                                                        
      ! purpose of this section is to calculate a new delta
      ! upon initialization                                             
      if(iterat .eq. 1 .or. delta .eq. -1.)then
          mu=0.0
          if(delta .eq. -1)then ! get trust region
              if(.not. noscale)then
                  do j=1,n
                      s(j)=gradient(j)/scale(j)
                  enddo
              else
                  call dcopy(n,gradient,1,s,1)
              endif
chsj              call euclidnorm(s,scale,n,.true.,alpha)
              alpha = enorm(n,s)
chsj              ! eliminate arithmetic problems from very large
chsj              ! (and certainly pathological) alpha                      
chsj              alpha = min(alpha, 1.0d8)
              alpha=alpha*alpha
              beta=0.0d0
              do i=1,n
                  temp=0.0d0
                  do j=i,n
                      dumy=1.0d0
                      if(.not. noscale)dumy=scale(j)
                      temp=temp+L(j,i)*gradient(j)/(dumy*dumy)
                  enddo
                  beta=beta+temp*temp
              enddo
              delta=alpha*sqrt(alpha)/beta
              if(delta .gt. steplength_max)delta=steplength_max
              ! if we came in with a gradient effectively at the
              ! solution, we can have an extremely small delta;         
              ! zero is even a possibility.  this will cause troubles
              ! in hookstep.  defend against that.
chsj              if (delta .lt. steplength_min) delta=steplength_min
          endif
      endif
!                                                                       
!                                                                       
      do while (retcode .ge. 2)
          call hookstep(n,gradient,l,hessian,step,scale,noscale,        
     &              epsilond,delta,mu,deltaprev,phi,phip,firsthook,     
     &              phipinit,s,newton_step,steplength,                  
     &              iounit)
          deltaprev=delta
          call trust_region(n,xcur,fcur,FVEC,gradient,l,s,scale,        
     &        noscale,Sf,FVp,newton_step,steplength_max,step_tol,
     &        step_type,
     &        hessian,delta,retcode,xprev,fprev,xnew,fnew,max_taken)
                                    
c         if(phip .eq. 0.0)retcode=1    !hsj 5/10/93 accept solution
c         print *,'phip,retcode =',phip,retcode
      enddo
!                                                                       
                                                  
      return
      end 

          SUBROUTINE HOOKSTEP(N,GRADIENT,L,HESSIAN,STEP,SCALE,NOSCALE,
     &                EPSILOND,DELTA,MU,DELTAPREV,PHI,PHIP,FIRSTHOOK,
     &                PHIPINIT,S,NEWTON_STEP,STEPLENGTH,IOUNIT)
C---------------------------------------------------------------------
C
C
      USE copy
      IMPLICIT none
      INTEGER N,I,J,K,N2,IOUNIT,TRY
      REAL *8 TEMPVEC(n),S(*),L(N,N),HESSIAN(*),STEP(*),GRADIENT(*),
     &        EPSILOND,ENORM
      REAL *8 HI,LOW,STEPLENGTH,MU,DELTAPREV,PHI,PHIP,SCALE(*),
     &      MUUP,DELTA,PHIPINIT,DUMY,MULOW,ADDMAX,
     &      STEPLS,POSDEF
      LOGICAL NEWTON_STEP,FIRSTHOOK,NOSCALE,DONE
C
C
      N2=2*N
      HI=1.5
      LOW=0.75
      TRY = 0 
c      print *,'newton step length, delta =',STEPLENGTH,DELTA
C     the scaled step,scale*s,must have length between low*delta and hi*delta
      IF(STEPLENGTH .LE. HI*DELTA)THEN !SET S TO NEWTON STEP
          NEWTON_STEP=.TRUE.
          CALL DCOPY(N,STEP,1,S,1)
          MU=0.0
          DELTA=MIN(DELTA,STEPLENGTH) 
         IF(IOUNIT .NE. 0)THEN
             WRITE(IOUNIT,*)'NEWTON STEP ACCEPTED in HOOKSTEP'
             WRITE(IOUNIT,2)STEPLENGTH,DELTA
    2         FORMAT(2X,'STEPLENGTH,DELTA =',2(1PE14.5))
         ENDIF
      ELSE
          NEWTON_STEP=.FALSE.
          IF(MU .GT. 0.0)
     &        MU=MU-(PHI+DELTAPREV)*((DELTAPREV-DELTA)+PHI)/(DELTA*PHIP)
          PHI=STEPLENGTH-DELTA
          IF(FIRSTHOOK)THEN     !calc phipinit
              FIRSTHOOK=.FALSE.
              DO J=1,N
                  DUMY=1.0
                  IF(.NOT. NOSCALE)DUMY=SCALE(J)
                  TEMPVEC(J)=STEP(J)*DUMY*DUMY
              ENDDO

              CALL LSOLVE(N,TEMPVEC,L,TEMPVEC)
chsj              CALL EUCLIDNORM(TEMPVEC,SCALE,N,.TRUE.,PHIPINIT)
              PHIPINIT = enorm(n,tempvec)
              PHIPINIT=-PHIPINIT*PHIPINIT/STEPLENGTH
          ENDIF
          MULOW=-PHI/PHIPINIT
          DO J=1,N
              DUMY=1.0
              IF(.NOT. NOSCALE)DUMY=1./SCALE(J)
              TEMPVEC(J)=GRADIENT(J)*DUMY
          ENDDO
chsj          CALL EUCLIDNORM(TEMPVEC,SCALE,N,.TRUE.,MUUP)
          MUUP = ENORM(N,TEMPVEC)
          MUUP=MUUP/DELTA
          DONE=.FALSE.
          DO WHILE (.NOT. DONE)
              IF((MU .LT. MULOW) .OR. (MU .GT. MUUP))
     &               MU=MAX(SQRT(MULOW*MUUP),1.E-3*MUUP)
              DO I=1,N
                  K=(I*(N2-I+3))/2-N
                  DUMY=1.0d0
                  IF(.NOT. NOSCALE)DUMY=SCALE(I)
                  HESSIAN(K)=HESSIAN(K)+MU*DUMY*DUMY
              ENDDO
              POSDEF=0.0
              CALL CHOLDECOMP(N,HESSIAN,POSDEF,EPSILOND,L,ADDMAX)
c            solve (Ltranspose *L)s = -g
              CALL CHOLSOLVE(N,GRADIENT,L,S)
C             RESTORE HESSIAN :
              DO I=1,N
                  K=(I*(N2-I+3))/2-N
                  DUMY=1.0
                  IF(.NOT. NOSCALE)DUMY=SCALE(I)
                  HESSIAN(K)=HESSIAN(K)-MU*DUMY*DUMY
              ENDDO
              CALL EUCLIDNORM(S,SCALE,N,NOSCALE,STEPLS)
              PHI=STEPLS-DELTA
              DO J=1,N
                  DUMY=1.0
                  IF(.NOT. NOSCALE)DUMY=SCALE(J)
                  TEMPVEC(J)=S(J)*DUMY*DUMY
              ENDDO
              CALL LSOLVE(N,TEMPVEC,L,TEMPVEC)
chsj              CALL EUCLIDNORM(TEMPVEC,SCALE,N,.TRUE.,PHIP)
              PHIP = ENORM(N,TEMPVEC)
              PHIP=-PHIP*PHIP/STEPLS
              IF(PHIP .EQ. 0.0)DONE=.TRUE. !HSJ 5/2/94 CANT DO ANY BETTER
              IF( .NOT. DONE)THEN          !SO ACCEPT THIS AND RETURN
                  TRY = TRY + 1
                  IF((STEPLS .GE. LOW*DELTA) .AND. (STEPLS .LE.
     &               HI*DELTA) .OR. (MUUP-MULOW .LE. 0.0D0) .or. TRY  
     &                    .gt. 1000 )THEN
                      DONE=.TRUE.    !ACCEPT S AS THE STEP
                  ELSE        !S NOT ACCEPTABLE,CALCULATE NEW MU,MULOW,MUUP:
                      MULOW=MAX(MULOW,MU-PHI/PHIP)
C                     IF(PHI .GT. 0.0)MULOW=MU
                      IF(PHI .LT. 0.0)MUUP=MU
                      MU=MU-(STEPLS*PHI)/(DELTA*PHIP)
                  ENDIF
              ENDIF
          ENDDO
C         IF(IOUNIT .NE. 0)THEN
C             WRITE(IOUNIT,1)MU
C   1         FORMAT(2X,'HOOKSTEP TAKEN,MU=',2X,1PE12.4)
C             WRITE(IOUNIT,2)STEPLENGTH,DELTA
C         ENDIF
      ENDIF
C
C
      RETURN
      END

      SUBROUTINE  jackrotate(n,i,alpha,beta,method,Z,LM)
c------------------------------------------------------------------------
c--- premultiply M,Z by Jacobi rotation matrix J(i,i+1,alpha,beta)
c--- this uses upper triangular LM plus the subdiagonal(of M).
c--- for diagoanl elment (i,i) the subdiagoanl is (i+1,i). this subdiagonal
c--  element is found in postion (n*(n+1))/2 + (i*(i+1))/2
c--  in vector LM (it waas put there in routine loadRT)

c----------------------------------------------------------HSJ-----------



      IMPLICIT none
c  INPUT:
      INTEGER n,method,i
      REAL *8 alpha,beta
C  LOCAL:
      INTEGER
     .        j
      REAL*8 
     .      den,c,s ,y,w
c INPUT/OUTPUT:
      REAL*8 LM(N,N),Z(n,n)


      IF(alpha .eq. 0.0d0)THEN
         c = 0.0D0
         s = sign(1.d0,beta)
      ELSE
         den =SQRT(alpha**2+beta**2)
         c = alpha/den
         s = beta/den
      ENDIF
      DO j =i,n
         y = LM(i,j)
         w = LM(i+1,j)
         LM(i,j) = c*y -s*w
         LM(i+1,j) = s*y + c*w
      ENDDO
      
      IF(method .eq. 1)THEN
         DO j=1,n
            y=Z(i,j)
            w = Z(i+1,j)
            Z(i,j) = c*y -s*w
            Z(i+1,j) = s*y +c*w
         ENDDO
      ENDIF
      
      RETURN
      END



      SUBROUTINE LINESEARCH(N,XCUR,FCUR,FVEC,GRADIENT,P,SCALE,NOSCALE,
     &                 Sf,MAXSTEP,STEPTOL,RETCODE,XNEW,FNEW,FVp,
     &                 MAX_TAKEN)
C----------------------------------------------------------------------
C---FIND THE APPROXIMATE MINIMUM ALONG A DESCENT LINE FROM XCUR.
C---THAT IS,GIVEN A VECTOR P SUCH THAT GRADIENTTRANSPOSE*P<0
C---FIND A LAMBDA SUCH THAT XNEW=XCUR+LAMBDA*P YIELDS
C---FNEW<FCUR+ALPHA*LAMBDA*GTRANSPOSE*P,WHERE FNEW AND FCUR ARE THE
C---OBJECTIVE FUNCTION EVALUATED AT XNEW AND XCUR RESPECTIVELY.
C---INPUT
C  N                  NUMBER OF ELEMENTS IN XCUR,P,SCALE,XNEW
C  XCUR(I)            I=1,2...N CURRENT SOLUTION POINT
C  FCUR               VALUE OF OBJECTIVE FUNCTION AT XCUR
C  FVEC                 SUBROUTINE NAME WHICH CALCULATES FNEW
C  GRADIENT(I)        GRADIENT OF OBJECTIVE FUNCTION AT XCUR
C  P(I)               (ANY) DESCENT DIRECTION FROM XCUR
C  SCALE(I)           DIAGONAL SCALING MATRIX
C  MAXSTEP            MAX ALLOWABLE STEP SIZE
C  STEPTOL            STEP TOLERANCE FOR DISTINGUISHING XNEW,XCUR
C  NOSCALE            IF .TRUE. DONT USE SCALE,
C                     IF .FALSE. SCALE CALCULATIONS USING SCALE
C  Sf(n)              Typ size of F(n)
C
C---OUTPUT:
C  RETCODE              INTEGER VARIABLE WITH MEANING AS FOLLOWS
C                       RETCODE=0   VALID XNEW FOUND
C                               1   FAILED TO FIND A SATISFACTORY XNEW,
C                                   SUFFICIENTLY DISTINC FROM XCUR
C                                   XCUR IS COPIED TO XNEW IN THIS CASE.
C                                   THIS CAUSES CHECK_CONVERGENCE TO
C                                   EXIT WITH ZERO MAX STEP
C  XNEW(I)              I=1,2..N THE NEW SOLUTION POINT
C  FNEW                 VALUE OF OBJECTIVE FUNCTION AT XNEW
C  MAX_TAKEN            LOGICAL VARIABLE =.TRUE. IF MAX STEP
C                       WAS TAKEN
C
C
C------------------------------------------------------------------HSJ
      
      USE copy
      IMPLICIT NONE
      INTEGER N,RETCODE,J
      REAL*8  SCALE(n),one,
     &     MAXSTEP,STEPTOL,ALPHA,STEPLENGTH,
     &     RELLENGTH,INITSLOPE,MINLAMBDA,LAMBDA,ABS,
     &     DUMY1,DUMY2,DUMY3,LAMBDATEMP,LAMBDAPREV,
     &     DISC,A,B,LAMBDASQINV,LAMBDAPREVSQINV,
     &     XCUR(n),Sf(n),FVp(n),
     &         GRADIENT(n),P(n),XNEW(n),DDDD,
     &         FCUR,FNEW,FNEWPREV,SDOT12,DSDOT
      LOGICAL MAX_TAKEN,NOSCALE
      EXTERNAL FVEC
C
C

      MAX_TAKEN=.FALSE.
      RETCODE=2
      ALPHA=1.E-4
      alpha = 1.e-8 !hsj 5/02/02
      one = 1.d0
C---GET LENGTH OF DESCENT VECTOR:
      CALL EUCLIDNORM(P,SCALE,N,NOSCALE,STEPLENGTH)
      IF(STEPLENGTH .GT. MAXSTEP)THEN !DECREASE LENGTH IF TOO LARGE
          DDDD= MAXSTEP/STEPLENGTH
          DO J=1,N
              P(J)=P(J)*DDDD
          ENDDO
          STEPLENGTH=MAXSTEP
      ENDIF
      INITSLOPE=SDOT12(GRADIENT,P,N)
c      INITSLOPE=DSDOT(GRADIENT,P,N)
      RELLENGTH=0.0
      DO J=1,N
          DUMY1=1.0
          IF( .NOT. NOSCALE)DUMY1=1./SCALE(J)
          DUMY1=MAX(ABS(XCUR(J)),DUMY1)
          RELLENGTH=MAX(RELLENGTH,ABS(P(J))/DUMY1)
      ENDDO
      MINLAMBDA=STEPTOL/RELLENGTH
      minlambda = .01*minlambda  !HSJ 5/02/02
      LAMBDA=one
C----------------------------------------------------------------------
C---CHECK IF XNEW=XCUR+LAMBDA*P IS SATISFACTORY. IF NOT GENERATE NEW LAMBDA
C---UNTIL IT IS O.K.:
C----------------------------------------------------------------------
      DO WHILE (RETCODE .GE. 2)
          DO J=1,N
              XNEW(J)=XCUR(J)+LAMBDA*P(J)
          ENDDO
chsj           CALL FVEC(N,XNEW,SF,FVp,FNEW)
          CALL FVEC(N,XNEW,FVP)
          CALL EUCLIDNORM(FVp,Sf,N,NOSCALE,FNEW)
          FNEW=0.5*FNEW*FNEW
c           print *,'FNEW,fcur ,lambda,minlambda,initslope',
c     .            FNEW,fcur,lambda,minlambda,initslope
c           print *,'steptol,rellength=',steptol,rellength

          IF(FNEW .LT. FCUR +ALPHA*LAMBDA*INITSLOPE)THEN
              RETCODE=0   !NEW STEP IS OK
              IF(LAMBDA .EQ. one .AND. (STEPLENGTH .GT. 0.99*MAXSTEP))
     &                                                  MAX_TAKEN=.TRUE.
          ELSE IF(LAMBDA .LT. MINLAMBDA)THEN !XNEW,XCUR NOT SUFFICIENTLY
                  RETCODE=1                  !DISTINCT
                  CALL DCOPY(N,XCUR,1,XNEW,1)
          ELSE  !REDUCE LAMBDA
                  IF(LAMBDA .EQ. one)THEN !QUADRATIC FIT ON FIRST TRY
                      LAMBDATEMP=-INITSLOPE/(2.d0*(FNEW-FCUR-INITSLOPE))
                  ELSE  !CUBIC FIT THEREAFTER
                      DUMY1=FNEW-FCUR-LAMBDA*INITSLOPE
                      DUMY2=FNEWPREV-FCUR-LAMBDAPREV*INITSLOPE
                      DUMY3=1./(LAMBDA-LAMBDAPREV)
                      LAMBDASQINV=one/(LAMBDA*LAMBDA)
                      LAMBDAPREVSQINV=one/(LAMBDAPREV*LAMBDAPREV)
                      A=DUMY3*(LAMBDASQINV*DUMY1-LAMBDAPREVSQINV*DUMY2)
                      B=DUMY3*(-LAMBDAPREV*LAMBDASQINV*DUMY1
     &                                +LAMBDA*LAMBDAPREVSQINV*DUMY2)
                      DISC=B**2-3.d0*A*INITSLOPE
                      DISC = MAX(DISC,0.0d0)       ! HSJ 04/15/04
                      IF(A .EQ. 0.0)THEN
                          LAMBDATEMP=-INITSLOPE/(2.d0*B)
                      ELSE
                          LAMBDATEMP=(-B+SQRT(DISC))/(3.d0*A)
                      ENDIF
                      IF(LAMBDATEMP .GT. 0.5d0*LAMBDA)LAMBDATEMP
     &                                                  =0.5d0*LAMBDA
                  ENDIF
                  LAMBDAPREV=LAMBDA
                  FNEWPREV=FNEW
                  IF(LAMBDATEMP .LE. 0.1*LAMBDA)THEN
                      LAMBDA=0.1*LAMBDA
                  ELSE
                      LAMBDA=LAMBDATEMP
                  ENDIF
          ENDIF
      ENDDO
C
C
      RETURN
      END
C
C
      SUBROUTINE LSOLVE(N,B,L,Y)
C---------------------------------------------------------------------
C---SOLVE L*Y=B
C---WHERE L IS LOWER TRIANGULAR
C   
C--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,I,J
      REAL *8 L(N,N),Y(*),B(*),SUMD
C
C
      Y(1)=B(1)/L(1,1)

      DO I=2,N
        SUMD=0.0D0

        DO J=1,I-1
          SUMD=SUMD+L(i,j)*Y(J)
        ENDDO
        Y(I)=(B(I)-SUMD)/L(i,i)
      ENDDO
C
C
      RETURN
      END
C
C
C
C
C
C
C
C
C
      SUBROUTINE MACHINE_EPSD(EPSILOND)
C-----------------------------------------------------------------------
C---GET THE (DOUBLE PRECISION) MACHINE EPSILON
C--------------------------------------------------------------------HSJ
C
      REAL *8 EPSILOND
      EPSILOND=1.0D0
      DO WHILE (1.D0+EPSILOND .GT. 1.0D0)
         EPSILOND=EPSILOND*0.5D0
      ENDDO
C
C
      RETURN
      END


               
      SUBROUTINE nefn(n,xp,fp,FVEC,NOSCALE)
c------------------------------------------------------------------------
c     get 0.5* sum of squares of residuals of nonlinear equations
c
c INPUT:
c   n
c   xp(n)
c  NOSCALE
c
c OUTPUT:
c
c fp         = 0.5*Norm(Sf*fVp)**2


c-------------------------------------------------------------HSJ-----------
      USE COM     !Sf(n),FVp(n) brought in here
      IMPLICIT NONE

      INTEGER
     .        n,i
      REAL*8 xp(n),fp 
      Logical NOSCALE
      EXTERNAL FVEC



c     evaluate the set of equations at xp:
chsj      CALL FVEC(n,xp,Sf,FVp,fp)
      CALL FVEC(n,xp,FVp)
      CALL EUCLIDNORM(FVp,Sf,N,NOSCALE,fp)
      fp=0.5*fp*fp
c      print *,'Fvp(1..N) in nefn'
c      do i =1,n
c         print *,'i,xp(i),Sf(i),FVp(i)',i,xp(i),Sf(i),FVp(i)
c      enddo
c      print *,'fp in nefn =',fp
      RETURN
      END



      SUBROUTINE neinck(n,macheps,x0,typx,typf,Fdigits,fvectol,
     .                   steptol,maxstep,itnlimit,printcode,
     .                   delta, global,analjac,cheapf,factsec,
     .                  Sx,Sf,eta,termcode)
c--------------------------------------------------------------------

c----- check input
c
c-------------------------------------------------------------HSJ---





      IMPLICIT NONE
c INPUT:
      INTEGER
     .        n,Fdigits,itnlimit,printcode,global
      REAL *8
     .        x0(n),typx(n),typf(n),delta,
     .         macheps, fvectol,steptol,maxstep
      LOGICAL
     .        analjac,cheapf,factsec
c LOCAL:
      INTEGER i
      
      REAL*8
     .       a

c OUTPUT
      INTEGER
     .      termcode
      REAL *8
     .      Sx(n),Sf(n),eta

c-------------------------------------------------------------------

      termcode = 0
      IF(n .lt. 1)THEN
         termcode = -1
         RETURN
      ENDIF
      DO i=1,n
         Sx(i) = 1./typx(i)
         Sf(i) = 1./typf(i)
      ENDDO
      IF(Fdigits .eq. -1)THEN
         eta = macheps
      ELSE
         a = 10.d0**(-Fdigits)
         eta = MAX(macheps,a)
      ENDIF
      if(global .eq. 2 .or. global .eq. 3 .and. delta .eq. 0.0D0)
     .         delta = -1.D0
      
      RETURN
      END


       SUBROUTINE nemodel(n,fc,jacf,g,Sf,Sx,macheps,global,
     .                  M,H,Sn)                               
c------------------------------------------------------------------------
c--- solve J*dx = -F(xc)
c    Factor model Jacobian,calculate Newton step
c    If Jacobian is singular modify it 
c---
c INPUT
c n              # unknowns
c fc(1..n)       Value of n dimensional function Fc(xc)  at xc
c jacf(n,n)      Jacobian of Fc
c g(1..n)        scaled gradient
c Sf(1..n)       scale factor for Fc
c Sx(1..n)       scale factors for xc
c macheps
c global         solution method 

c global         = 1    line search 
c global         = 2    Hookstep trust region
c global         = 3    Dogleg trust region 

c OUPUT:
c M(n,n)        M is put into LOWER TRIANGULAR SYMMETRIC STORAGE
c               MODE on output.
c
c
c H             H IS IN UPPER TRIANGULAR SYMMETRIC STORAGE MODE.
c               (IE H(I,J),FOR J  .GE. I,IS ACCESSED AS H(K),WHERE
c               K=((I-1)(2*N-I)+2*J)/2
c Sn(n)         The Newton step
c-----------------------------------------------------------HSJ-------------

       USE terminate

      IMPLICIT none
c  INPUT:
      INTEGER n,global
      REAL *8 fc(n),jacf(n,n),g(n),Sf(n),Sx(n),macheps
C  LOCAL:
      INTEGER
     .        i,j,k

      REAL*8
     .       M1(n),M2(n),est,sumd,Hnorm,temp,maxadd,zero
      LOGICAL sing
c  OUTPUT:
      REAL*8 M(n,n),H(*),Sn(n) !H has  n(n+1)/2 elements


c TEMPORARY
      REAL *8 work(n)
      INTEGER info 
      zero = 0.0d0


      DO i = 1,n
         DO j = 1,n 
            M(i,j) = Sf(i)*jacf(i,j) !i is equation, j is variable number
         ENDDO                       !equation i is scaled by Sf(i)
      ENDDO


c------------------------------------------------------------
c      qrdecomp,coded from the Dennis- Schnabel algorithm
c      appears to be faulty in some way. So replace it with
c      Linpack routine dgeqrf and compensate for different
c      storage scheme:
c      CALL qrdecomp(n,M,M1,M2,sing)

      sing = .false.
      Call Dgeqrf(n,n,M,n,M1,work,n,info) ! lapack routine
c     returns qr factorization of scaled jacobian in M, presumabley, if 
c     a diagonal element of r is zero then we dont crash
c     

      do i = 1,n
         M2(i) = M(i,i)  !M1 is ony used in qrsolve
                         !which is also replaced so M1 becomes unused
         if(M2(i) .eq. 0.0d0)sing = .true.
      enddo 

      IF( .not. sing)THEN
         DO j =1 ,n
            DO i = 1,j-1
               M(i,j) = M(i,j)/Sx(j)
            ENDDO
            M2(j) = M2(j)/Sx(j)
         ENDDO
         CALL condest(n,M,M2,est)
      ELSE
         est = zero
      ENDIF  ! sing 
c       print *,'sing,est  =',sing,est 
      IF(sing .or. est .gt. 1.0d0/SQRT(macheps))THEN
         print *,'JACOBIAN SINGULAR************************************'
         print *,'sing,etst,1/sqrt(meps) ',sing,est,1.0d0/SQRT(macheps)
c        H = Jtranspose* Sf^2*J
         DO i =1,n   !perturb jacobian so it isnt singular
            DO j =i,n
               sumd = 0.0d0
               DO k =1,n
                  sumd = sumd + jacf(k,i)*jacf(k,j)*Sf(k)**2
               ENDDO
               k=((i-1)*(2*n-i)+2*j)/2  !upper triangular storage for H
               H(k) = sumd 
            ENDDO
         ENDDO
         sumd =0.0d0
         DO j =1,n
            sumd =sumd + ABS(H(j))/Sx(j)
         ENDDO
         Hnorm = sumd/Sx(1)
         DO i = 2,n
            sumd =0.0D0
            DO  j= 1,i
                K=((j-1)*(2*N-j)+2*i)/2
                sumd = sumd + ABS(H(k))/Sx(j)
            ENDDO
            temp = sumd
            sumd =0.0D0
            DO j = i+1,n
                k=((i-1)*(2*n-i)+2*j)/2
                sumd = sumd + ABS(H(k))/Sx(j)
            ENDDO
            temp = temp +sumd
            temp = temp/Sx(i)
            Hnorm = MAX(Hnorm,temp)
         ENDDO
         DO i = 1,n
            k=((i-1)*(2*n-i)+2*i)/2
            H(k) = H(k)+ SQRT(n*macheps)*Hnorm*Sx(i)**2
         ENDDO
       
         
         CALL choldecomp(n,H,zero,macheps,M,maxadd)
         CALL cholsolve(n,g,M,Sn)


      ELSE  ! jacobian not singular, get normal Newton step:
         DO j = 1,n
            DO i =1 , j-1
               M(i,j) = M(i,j)*Sx(j)
            ENDDO
            M2(j) = M2(j)*Sx(j)
         ENDDO
         DO j= 1,n
            Sn(j) = -Sf(j)*Fc(j)      
         ENDDO


c-----------------------------------------------------------------------
c        now solve J*dx = -F(xc), where dx is the Newton step increment
c        (dx = xc new - xc prev). J is in factored form from above:
c     dont use qrsolve because qrdecomp doesn't produce the
c     right factors! (see above):
c         Call qrsolve(n,M,M1,M2,Sn)
c     qrsolve doesnt change M,M1,M2 . output is Sn (eg the Newton step)
c
c      REPLACEMENT for qrsolve:
c     Recall that in call to Dgeqrf we ended up with the matrix R 
c     in the upper triangle of M.
c      To solve the system
c           J*dx = Sn  (see Sn above,defined with minus sign included)
c     we have J = Q_1 * R   (Q = (Q_1,Q-2) , Q_1,Q_2 orthogonal, but
c     Q_2 drops out because it multiplies the lower half of the factored
c     matrix which is zero)
c      hence Q_1*R*dx  = Sn is just R*dx = Q_1(transpose)*Sn
c     this upper triangular system is easily solved.
C     The reason why all this is done is to control what happens
c     when J is near singular. 
c     First get Q_1 Transpose *Sn:


         CALL DORMQR('L','T',n,1,n,M,n,M1,Sn,n,work,n,info)
         if(info .lt. 0)
     .   CALL STOP('subroutine NEMODEL: DORMQR error', 0)


c     Now solve R * dx = Sn:

         CALL  DTRTRS('U','N','N',n,1,M,n,Sn,n,info)
c        if info = i>0  then element r(i,i) is zero
c        we checked for this above after call to Dgeqrf
c        so this situation shouldnt arise here:
         if(info .gt. 0)
     .   CALL STOP('subroutine NEMODEL: R singular', 0)

c---------------------------------------------------------------------


         IF(global .eq. 2 .or. global .eq. 3) THEN
            !put R transpose into lower triangle of M:
            DO i =1,n
               M(i,i) = M2(i)
               Do j= 1,i-1
                  M(i,j)=M(j,i)
               ENDDO
            ENDDO
         ENDIF     
         IF(global .eq. 2) THEN
            !store  R transpose* R in upper triangular symmetric storage
            !mode in H: 
               DO i =1,n                       
                  sumd =0.0d0
                  DO k =1,i
                     sumd =sumd + M(i,k)**2
                  ENDDO
                  K=((i-1)*(2*n-i)+2*i)/2         !picks out diagonal elements
                  H(k) = sumd
                  DO j = i+1,n
                     sumd =0.0d0
                     DO  k = 1,i
                        sumd =sumd + M(i,k)*M(j,k)
                     ENDDO
                     K=((I-1)*(2*N-I)+2*J)/2
                     H(k) = sumd
                  ENDDO
               ENDDO 
          ENDIF
         
         
      ENDIF


      RETURN
      END


      SUBROUTINE nemodelfac(n,fc,g,Sf,Sx,macheps,global,
     .                  restart,M,M2,jacf,H,Sn)
c------------------------------------------------------------------------
c--- Factor model Jacobian,calculate Newton step,If Jacobian is singular
c--- modify it
c INPUT
c n
c fc(n)
c g(n)
c Sf(n)
c Sx(n)
c macheps
c global
c restart

c INPUT/OUTPUT
c M(n,n) M2(n),jacf(n,n)
c              If this is a restart then M contains no useful
c              information on input
c              otherwise M contains R in its upper triangle, including the
c              diagonal elements.
c              On output M contains R in its upper triangle and
c              diagonal elements are in M2.

c OUTPUT:
c H(n)              upper triangular storage mode
c Sn(n)

c------------------------------------------------------------------------



      IMPLICIT none
c  INPUT:
      INTEGER n,global
      REAL *8 fc(n),g(n),Sf(n),Sx(n),macheps
      LOGICAL restart
C  LOCAL:
      INTEGER
     .        i,j,k,info

      REAL*8
     .       M1(n),work(n),est,sumd,Hnorm,temp,maxadd,zero,
     .       jactr(n,n)
      LOGICAL sing
c  OUTPUT:
      REAL*8 M(n,n),M2(n),jacf(n,n),H(*),Sn(n)



      zero = 0.0D0


      IF(restart) THEN                   !get new QR factorization
         !we start with a two dimensional M:

         DO i =1,n
            DO j=1,n
               M(i,j)= Sf(i)*jacf(i,j)
            ENDDO
         ENDDO

c        CALL qrdecomp(n,M,M1,M2,sing)
c        M is in factored form (qform does not modify M)
c        R is in upper triangle of M but diagonal elements
c        of R are in M2. Lower triangle of M contains
c        information to generate Q.
c        next factored form of M is used to get Q transpose
c        returned in jacf:
c        Call qform(n,M,M1,jacf)
c                   HOWEVER
c        qrdecomp,coded from the Dennis- Schnabel algorithm
c        appears to be faulty in some way. So replace it with
c        Linpack routine dgeqrf and compensate for different
c        storage scheme:

        sing = .false.  !we have to assume this here Dgeqr doesnt supply info
        jacf(:,:) = M(:,:)
        Call Dgeqrf(n,n,jacf,n,M1,work,n,info)      ! lapack routines
c       jacf now contains R in its upper triangle,including the diagonal
c       Put info into M,m2 as would have been done by qrdecomp:

        do i = 1,n
           M2(i) = jacf(i,i) 
        enddo   
        M(:,:) = jacf(:,:)              !puts M into factored form
                                        !as would have been done by qrdecomp
        Call Dorgqr(n,n,n,jacf,n,M1,work,n,info)    !jacf now contains Q transpose


         jactr(:,:) = jacf(:,:)

        DO i=1,n

          DO j=1,n
               jacf(i,j) = jactr(j,i)
          ENDDO
        ENDDO



      ELSE  !restart = .false.
         sing = .false.
         DO i =1,n
            if(M(i,i) .eq. 0.0d0 )sing  = .TRUE.
         ENDDO
      ENDIF
      

      IF(.not. sing) THEN
         DO j = 1,n
            DO i = 1,j-1
                M(i,j) = M(i,j)/Sx(j)
            ENDDO
            M2(j) = M2(j)/Sx(j)
         ENDDO
         CALL condest(n,M,M2,est)
         DO j = 1,n
            DO i = 1,j-1
                M(i,j) = M(i,j)*Sx(j)
            ENDDO
            M2(j) = M2(j)*Sx(j)
         ENDDO
      ELSE
         est =0.0d0
      ENDIF


      If(sing .or. (est .gt. 1.D0/SQRT(macheps)))THEN
        !perturb Jacobian H = RT*R
         DO i =1,n
            DO j =i,n
               sumd =0.0D0
               DO k =1,i
                  sumd =sumd +M(k,i)*M(k,j)
               ENDDO
               k = ((i-1)*(2*n-i)+2*j)/2  !upper triangular storage for H
               H(k) = sumd 
            ENDDO
         ENDDO
         sumd =0.0D0
         DO j =1,n
c            sumd =sumd +ABS(H(1,j))/Sx(j)
            sumd =sumd +ABS(H(j))/Sx(j)
         ENDDO
         Hnorm = sumd/Sx(1)
         DO i =2,n
            sumd =0.0D0
            DO j =1,i
               k = ((j-1)*(2*n-j)+2*i)/2
               sumd =sumd +ABS(H(k))/Sx(j)
            ENDDO
            temp =sumd
            sumd =0.0d0
            DO j =i+1,n
               k = ((i-1)*(2*n-i)+2*j)/2
               sumd =sumd+ABS(H(k))/Sx(j)
            ENDDO
            temp = temp +sumd
            temp=temp/Sx(i)
            Hnorm = MAX(temp,Hnorm)
         ENDDO
         DO i =1,n
            k = ((i-1)*(2*n-i)+2*i)/2
            H(k) = H(k)+SQRT(n*macheps)*Hnorm*Sx(i)**2
         ENDDO


         !choldecomp takes H and outputs M (as a lower
         !triangular matrix)
         CALL choldecomp(n,H,zero,macheps,M,maxadd)
         !cholsolve takes lower triangular M and 
         !solves for Sn, the "newton' step. g is not changed:
         CALL cholsolve(n,g,M,Sn)
      ELSE
         !do normal Newton step
         !note that jacf  = Q transpose here
         !note that M contains R in upper triangle except  diagonal
         !which  is in M2
         DO i =1,n
            sumd =0.0D0
            DO j =1,n
               sumd =sumd+jacf(i,j)*Sf(j)*Fc(j)
            ENDDO
            Sn(i) = -sumd
         ENDDO
         !rsolve does not change M:
         CALL rsolve(n,M,M2,Sn)

          if(restart) print *,"****************************restart****"


         IF(global .eq. 2 .or. global .eq. 3)THEN
c            copy R transpose into lower triangle of M:
            DO i =2,n
               Do j = 1,i-1
                  M(i,j) = M(j,i)
               ENDDO
            ENDDO
         ENDIF      

          IF(global .eq. 2)THEN
c             put L*Ltranspose into H in symmetric storage mode
             DO i =1,n
                DO j =i,n
                   sumd =0.0D0
                   DO k = 1,i
                      sumd =sumd +M(i,k)*M(j,k)
                   ENDDO
                   k = ((i-1)*(2*n-i)+2*j)/2
                   H(k) =sumd
                ENDDO
             ENDDO
          ENDIF
      ENDIF

      RETURN
      END
      




      SUBROUTINE nestop0(n,F,Sf,fvectol,termcode,consecmax)
c-----------------------------------------------------------------------
c---- determine if we should stop at iteration 0
c---- because the initial guess is good enough.
c---- return termcode =0 if the guess x0 is not sufficient
c---- return termcode =1 if x0 is a sufficient approximation to the
c---- root of F(x0) =0.0
c---- condition is that maximum component is .01*fvectol or less

c--------------------------------------------------------------HSJ------

      IMPLICIT none
c  INPUT:
      INTEGER
     .        n
      REAL*8
     .        f(n),Sf(n),fvectol
c  LOCAL
      INTEGER
     .        i
      REAL*8
     .         mmax,factor
c  OUTPUT:
      INTEGER
     .        termcode,consecmax

      
      factor = 1.D-2
      consecmax = 0
      mmax =0.0d0
      DO i = 1,n
         mmax = MAX(mmax,Sf(i)*ABS(F(i)))
      ENDDO
      IF(mmax .lt. factor*fvectol)THEN
         termcode = 1
      ELSE
         termcode =0
      ENDIF


      RETURN
      END

      
      SUBROUTINE  qrupdate(n,u,v,method,Z,LM)
c------------------------------------------------------------------------
c--- This routine is part of the scheme for updateing  A to A+.
c--- Given QR factorization of matrix A, get the factorization of
c--- Q+*R+ of A+ = Q(R+u*vt) 
C--- LM contains R in upper triangular
c--- LM contains R+ on output

c------------------------------------------------------------------------



      IMPLICIT none
c  INPUT:
      INTEGER n,method
      REAL *8 u(n),v(n)
C  LOCAL:
      INTEGER
     .        i,j,k
      REAL*8 
     .       uip1,mip
c INPUT/OUTPUT:
      REAL*8 LM(N,N),Z(n,n)


      k = n
      DO WHILE (u(k) .eq. 0.0d0 .and. k .gt. 1)
         k = k -1
      ENDDO

      DO i = k-1,1,-1
         uip1 = - u(i+1)
         CALL jackrotate(n,i,u(i),uip1,method,Z,LM)
         IF(u(i) .eq. 0.0D0)THEN
            u(i) = ABS(u(i+1))
         ELSE
            u(i) = SQRT(u(i)**2 + u(i+1)**2)
         ENDIF
      ENDDO
      DO j =1,n
       LM(1,j) = LM(1,j)+U(1)*v(j)
      ENDDO
      
      DO i =1,k-1
        mip = -LM(i+1,i)
         CALL jackrotate(n,i,LM(i,i),mip,method,Z,LM)
      ENDDO
      

      RETURN
      END


      SUBROUTINE rsolve(n,M,M2,B)
c--------------------------------------------------------------------------
c  Solve Rx = b 
c  R is stored in upper triangle of M and diagonal of R is in M2
c on output B contains solution x
c------------------------------------------------------------------HSJ
      IMPLICIT none
      INTEGER 
     .       n
      REAL *8
     .       M(n,n)
c LOCAL:
       INTEGER j,i
       REAL*8
     .       sumd

!OUTPUT:
      REAL*8
     .     M2(n),B(n)


      B(n) =B(n)/M2(n)


      DO i =n-1,1,-1
         sumd =0.0D0
         DO j = i+1,n
            sumd =sumd+M(i,j)*B(j)
         ENDDO
         B(i) = (B(i) - sumd)/M2(i)
      ENDDO

      RETURN
      END

C
C
C
C
      REAL*8  FUNCTION SDOT12(A,B,N)
C----------------------------------------------------------------------
C---CALCULATE THE DOT PRODUCT OF VECTORS A AND B IN DOUBLE PRECISION
C---HERE WE UNROLL THE DO LOOP FOR BETTER PERFORMANCE
C--------------------------------------------------------------------HSJ
C
C
      IMPLICIT NONE
      REAL*8  SUMD,A(*),B(*)
      INTEGER M,N,MP1,I
C
C
      SUMD=0.0D0
      M=MOD(N,5)
      IF(M .EQ. 0)GO TO 40
      DO I=1,M
         SUMD=SUMD+A(I)*B(I)
      ENDDO
      IF(N .LT. 5)GO TO 50
   40 MP1=M+1

      DO I=MP1,N,5
         SUMD=SUMD+A(I)*B(I)+A(I+1)*B(I+1)
     &          +A(I+2)*B(I+2)+A(I+3)*B(I+3)
     &          +A(I+4)*B(I+4)
      ENDDO 
   50 SDOT12=SUMD
C
C
      RETURN
      END
C
C
C
      SUBROUTINE TRUST_REGION(N,XCUR,FCUR,FVEC,GRADIENT,L,S,SCALE,
     &         NOSCALE,Sf,FVp,
     &         NEWTON_STEP,STEPLENGTH_MAX,STEP_TOL,STEP_TYPE,HESSIAN,
     &         DELTA,RETCODE,XPREV,FPREV,XNEW,FNEW,MAX_TAKEN)
C-----------------------------------------------------------------------
C   GIVEN THE STEP S(I),I=1,2..N,PRODUCED BY THE HOOKSTEP OR DOGLEG
C   METHODS,DECIDE IF XNEW=XCUR+S IS ACCEPTABLE AS THE NEXT ITERATE.
C   IF SO THEN ADJUST THE INITIAL TRUST REGION FOR THE NEXT ITERATION
C   APPROPRIATELY. IF NOT THEN DECREASE OR INCREASE THE TRUST REGION
C   RADIUS FOR THE CURRENT ITERATION.
C---INPUT
C  N                  NUMBER OF ELEMENTS IN XCUR,GRADIENT,S,SCALE,XNEW
C  XCUR(I)            I=1,2...N CURRENT SOLUTION POINT
C  FCUR               VALUE OF OBJECTIVE FUNCTION AT XCUR
C  FVEC                 SUBROUTINE NAME WHICH CALCULATES FNEW
C  GRADIENT(I)        GRADIENT OF OBJECTIVE FUNCTION AT XCUR
C  S(I)               I=1,2..N,THE STEP TO BE TESTED.
C  L(N,N))            CHOLESKY DECOMPOSITION OF HESSIAN
C  SCALE(I)           DIAGONAL SCALING MATRIX
C  NOSCALE            IF .TRUE. DONT USE SCALE,
C                     IF .FALSE. SCALE CALCULATIONS USING SCALE
C  NEWTON_STEP        BOOLEAN,=.TRUE. ONLY IF FULL NEWTON STEP WAS TAKEN
C  STEPLENGTH_MAX     MAX ALLOWABLE STEP SIZE
C  STEP_TOL           STEP TOLERANCE LEVEL
C  STEP_TYPE          =1  FOR HOOKSTEP,=2 FOR DOUBLE DOGLEG
C  HESSIAN           in symmetric storage mode

C---INPUT AND OUTPUT:
C  DELTA                THE TRUST REGION RADIUS
C  RETCODE              INTEGER VARIABLE WITH MEANING AS FOLLOWS
C                       RETCODE=0   VALID XNEW FOUND,DELTA SET TO NEW VALUE
C                               1   FAILED TO FIND A SATISFACTORY XNEW,
C                                   SUFFICIENTLY DISTINC FROM XCUR
C                                   BECAUSE LENGTH OF XNEW-XCUR IS LESS
C                                   THAN STEP TOL. THIS MEANS THAT NO
C                                   FURTHER PROGRESS TOWARD THE SOLUTION
C                                   IS POSSIBLE. IN ALMOST ALL CASES THIS
C                                   MEANS THAT XCUR IS AN ACCEPTABLE SOLUTION.
C                               2   FNEW IS TOO LARGE. THIS MEANS THAT
C                                   THE TRUST REGION RADIUS WAS SO LARGE
C                                   THAT THE DESCENT STEP TOTALLY OVERSHOT
C                                   THE SOLUTION POINT. WE HAVE TO DECREASE
C                                   THE TRUST REGION RADIUS DELTA AND
C                                   TRY AGAIN
C                               3   FNEW IS SUFFICIENTLY SMALL BUT
C                                   THE PROBABILITY OF TAKING A LONGER
C                                   SUCCESSFUL STEP IS LARGE. HENCE
C                                   CONTINUE THE CURRENT ITERATION WITH
C                                   A LARGER TRUST REGION RADIUS DELTA.
C  XPREV
C  FPREV
C---OUTPUT:
C  XNEW(I)              I=1,2..N THE NEW SOLUTION POINT
C  FNEW                 VALUE OF OBJECTIVE FUNCTION AT XNEW
C  MAX_TAKEN            LOGICAL VARIABLE =.TRUE. IF MAX STEP
C                       WAS TAKEN
C
C
C------------------------------------------------------------------HSJ
C
C
      USE copy
      IMPLICIT NONE
      INTEGER N,RETCODE,I,J,STEP_TYPE,
     &        N2,IHIJ,IHII
      REAL*8 SCALE(*),
     &     STEPLENGTH_MAX,STEP_TOL,ALPHA,STEPLENGTH,
     &     RELLENGTH,INITSLOPE,ABS,
     &     DELTAF,DELTAF_PREDICTED,
     &     DUMY,DELTA,SUM
      REAL *8 L(N,N),XCUR(*),GRADIENT(*),SDOT12,DSDOT,
     &      S(*),HESSIAN(*),Sf(N),FVp(n),
     &      XPREV(*),XNEW(*),FNEW,FCUR,FPREV
      LOGICAL MAX_TAKEN,NOSCALE,NEWTON_STEP
      EXTERNAL FVEC
C
C
      N2=2*N
      MAX_TAKEN=.FALSE.
      ALPHA=1.E-04
      CALL EUCLIDNORM(S,SCALE,N,NOSCALE,STEPLENGTH)
c       print *,'steplength in trust region',steplength

      DO J=1,N
        XNEW(J)=XCUR(J)+S(J)
      ENDDO
C---GET THE FUNCTION VALUE AT XNEW,DERIVATIVES ARE NOT USED:
CHSJ      CALL FVEC(N,XNEW,Sf,FVp,FNEW)
      CALL FVEC(N,XNEW,FVp)
      CALL EUCLIDNORM(FVp,Sf,N,NOSCALE,FNEW)
      FNEW=0.5*FNEW*FNEW
      DELTAF=FNEW-FCUR


      INITSLOPE=SDOT12(GRADIENT,S,N)
c       INITSLOPE=DSDOT(GRADIENT,S,N)
c      print *,'fcur,fnew =',fcur,fnew
c      print *,'initslope,retcode ',initslope,retcode
      IF(RETCODE .NE. 3)FPREV=0.0
      IF((RETCODE .EQ. 3) .AND. ((FNEW .GE. FPREV) .OR. (DELTAF .GT.
     &                                         ALPHA*INITSLOPE)))THEN
          RETCODE=0
          CALL DCOPY(N,XPREV,1,XNEW,1)
          FNEW=FPREV
          DELTA=DELTA*0.5
      ELSE IF(DELTAF .GE. ALPHA*INITSLOPE)THEN
          RELLENGTH=0.0    
          DO J=1,N
              DUMY=1.0
              IF( .NOT. NOSCALE)DUMY=1./SCALE(J)
              DUMY=MAX(ABS(XNEW(J)),DUMY)
              RELLENGTH=MAX(RELLENGTH,ABS(S(J))/DUMY)
          ENDDO
          IF(RELLENGTH .LT. STEP_TOL)THEN !XNEW-XCUR TOO SMALL
              RETCODE=1
              CALL DCOPY(N,XCUR,1,XNEW,1)
          ELSE
C---REDUCE DELTA,CONTINUE GLOBAL STEP
              RETCODE=2
              DUMY=-(INITSLOPE*STEPLENGTH)/(2.0*(DELTAF-INITSLOPE))
              IF(DUMY .LT. 0.1*DELTA)THEN
                  DELTA=0.1*DELTA
              ELSE IF(DUMY .GT. 0.5*DELTA)THEN
                  DELTA=0.5*DELTA
              ELSE
                  DELTA=DUMY
              ENDIF
          ENDIF
      ELSE !FNEW IS SUFFICIENTLY SMALL
          DELTAF_PREDICTED=INITSLOPE
          IF(STEP_TYPE .EQ. 1)THEN
C---HOOKSTEP CALCULATE ST*HESSIAN*S:
              DO I=1,N
                  SUM=0.0
                  DO J=I+1,N
                      IHIJ=((I-1)*(N2-I)+2*J)/2
                      SUM=SUM+HESSIAN(IHIJ)*S(J)
                  ENDDO
                  IHII=(I*(N2-I+3))/2-N
                  DUMY=S(I)*(0.5*HESSIAN(IHII)*S(I)+SUM)
                  DELTAF_PREDICTED=DELTAF_PREDICTED+DUMY
              ENDDO
          ELSE
C---DOGLEG CALCULATE ST*L*LT*S:
              DO I=1,N
                  SUM=0.0
                  DO J=I,N
                      SUM=SUM+L(j,i)*S(J)
                  ENDDO
                  DELTAF_PREDICTED=DELTAF_PREDICTED+0.5*SUM*SUM
              ENDDO
          ENDIF
          DUMY=ABS(DELTAF_PREDICTED-DELTAF)
          IF(RETCODE .NE. 2 .AND. ((DUMY .LE. 0.1*ABS(DELTAF))
     &        .OR. (DELTAF .LE. INITSLOPE)) .AND.( .not. NEWTON_STEP)
     &              .AND. (DELTA .LE. 0.99*STEPLENGTH_MAX))THEN
C---DOUBLE DELTA AND CONTINUE THE ITERATION
              RETCODE=3 
              CALL DCOPY(N,XNEW,1,XPREV,1)
              FPREV=FNEW
              DELTA=MIN(2.*DELTA,STEPLENGTH_MAX)
          ELSE
C---ACCEPT XNEW AS NEW ITERATE,CHOOSE NEW TRUST REGION RADIUS DELTA:
              RETCODE=0
              IF(STEPLENGTH .GT. .99*STEPLENGTH_MAX)MAX_TAKEN=.TRUE.
              IF(DELTAF .GE. 0.1*DELTAF_PREDICTED)THEN
                  DELTA=0.5*DELTA
              ELSE IF(DELTAF .LE. 0.75*DELTAF_PREDICTED)THEN
                  DELTA=MIN(2.*DELTA,STEPLENGTH_MAX)
              ENDIF
          ENDIF
      ENDIF
C
C 
      RETURN
      END
