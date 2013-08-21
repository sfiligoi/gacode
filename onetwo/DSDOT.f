      REAL*8  FUNCTION DSDOT(A,B,N)
C----------------------------------------------------------------------
C---CALCULATE THE DOT PRODUCT OF VECTORS A AND B IN QUAD PRECISION
C---HERE WE UNROLL THE DO LOOP FOR BETTER PERFORMANCE
C---QUAD precision not available with pgi compiler
C--------------------------------------------------------------------HSJ
C
C
      IMPLICIT NONE
      REAL*8  A(*),B(*)

      REAL*16
     . QEXT,DSUM,ai,bi, aip1,aip2,aip3,aip4,bip1,bip2,bip3,bip4
       !problem with pgi compiler, *16 data type not accepted
      INTEGER M,N,MP1,I
C
C
      DSUM=0.0D0
      M=MOD(N,5)
      IF(M .EQ. 0)GO TO 40
      DO I=1,M
         DSUM=DSUM+QEXT(A(I))*QEXT(B(I))
      ENDDO
      IF(N .LT. 5)GO TO 50
   40 MP1=M+1
      DO I=MP1,N,5
         ai=QEXT(A(I))
         bi=Qext(B(I))
         aip1 = QEXT(A(I+1))
         bip1 = QEXT(B(I+1))
         aip2 = QEXT(A(I+2))
         bip2 = QEXT(B(I+2))
         aip3 = QEXT(A(I+3))
         bip3 = QEXT(B(I+3))
         aip4 = QEXT(A(I+4))
         bip4 = QEXT(B(I+4))
         DSUM = DSUM + ai*bi+aip1*bip1+aip2*bip2+
     &                 aip3*bip3+aip4*bip4

      ENDDO
   50 DSDOT=DSUM
C
C
      RETURN
      END
