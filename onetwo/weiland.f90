
      MODULE weiland

	USE param,                           ONLY : kj
        USE nrtype,                          ONLY : I4B,DP
	IMPLICIT NONE
        REAL(DP)                                                        &
	                wweiland, xchie_weiland(kj), xchii_weiland(kj), &
                        d_weiland(kj), qe_weiland(kj), qi_weiland(kj),  &
                        xchie_weilandsv(kj), xchii_weilandsv(kj),       &
                        xki_weiland(kj), xke_weiland(kj), xkangwl(kj),  &
                        xkangwlsv(kj), xeffwl(kj), time_weiland   
     
        INTEGER(I4B)  include_weiland

      END MODULE weiland
