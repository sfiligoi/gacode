    MODULE neo2d
! --- flux surface quantitites, interpolated onto the transport grid,
! --- from previous (...0) and current calcs
!
    USE param, only : kj

    IMPLICIT NONE
      REAL *8        eps0(kj), xhm20(kj), xi110(kj), xi330(kj),       &
                    xips0(kj), eps(kj), xhm2(kj), xi11(kj), xi33(kj), &
                    xips(kj), ftfc(kj), depsdt(kj), dxhm2dt(kj),      &
                    dxi33dt(kj), dxi11dt(kj), dxipsdt(kj),            &
                    ftncl(kj),elong_r(kj),h88l31(kj),h88l32(kj) 



    END MODULE neo2d
