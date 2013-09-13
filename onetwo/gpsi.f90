        
   MODULE gpsi
 
     USE nrtype,                    ONLY : I4B,DP
     USE param,                     ONLY : nwmhd,nhmhd
     IMPLICIT NONE

      REAL(DP),dimension(:,:)    :: p(nwmhd,nhmhd)
      REAL(DP)	pmin,pmax,psimx(2),xax(2),yax(2),elongax
      REAL(DP),ALLOCATABLE,DIMENSION(:) :: pppsi_eqdsk,presspsi_eqdsk, &
               fpsi_eqdsk,ffppsi_eqdsk,qpsi_eqdsk,psival_eqdsk
      INTEGER(I4B) nxeqd_eqdsk
      REAL(DP) psimag_eqdsk,psilim_eqdsk
   ENDMODULE gpsi
