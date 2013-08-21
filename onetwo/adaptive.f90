   MODULE adaptive
      USE param, only : kj,kk
      implicit none
!
      integer           include_adaptive, include_curvature
      integer,private :: j
!
      real*8            eps_adaptive(kk), dzeta_adaptive, r_adaptive(kj), &
                       drhodt_adaptive(kj), speed_adaptive,               &
                       freeze_adaptive, curve_eps_adaptive(kk)
!
!
      data  (drhodt_adaptive(j), j=1,kj)          /kj   * 0.0/
!
!
!
   END MODULE ADAPTIVE
