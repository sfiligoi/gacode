

MODULE pelcom
   USE nrtype,                              ONLY : I4B,DP
   IMPLICIT NONE
   SAVE
!
   CHARACTER       nampel*8
!
   REAL(DP)                                            &
        pelrad, vpel, timpel(10)
   
   INTEGER(I4B)                                        &
        ipel, nbgpel, ipelet, npel, pelmod

 END MODULE pelcom 
