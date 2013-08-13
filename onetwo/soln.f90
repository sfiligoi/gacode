      MODULE SOLN
      USE param,only : kj,kk,kion
! --- old INCLUDE file soln.i
      implicit none
!
!
      integer *4,public :: njcur 
      integer *4,public :: curtype 
      integer *4,public :: diffeq_methd
!
      real *8,public :: u(kk,kj)
      real *8,public ::usave(kk,kj)
      real *8,public :: en(kj,kion)
      real *8,public ::te(kj)
      real *8,public ::ti(kj)
      real *8,public ::rbp(kj)
      real *8,public ::ene(kj)
      real *8,public ::enesav(kj)
      real *8,public ::curden(kj)
      real *8,public ::etor(kj)
      real *8,public ::uav(kk,kj)
      real *8,public ::uav0(kk,kj)
      real *8,public ::eneav0(kj)
      real *8,public ::eneav1(kj)
      real *8,public ::curpar_soln(kj)
      real *8,public ::cur_tor_ps_soln(kj)
      real *8,public ::rbp_save(kj)
      real *8,public ::u_continue(kk,kj)
      real *8,public ::cur_seed(kj)
!
!toroidal rotation speed profile,angrot, is kept in tordlrot.i
!       curpar_soln is parallel current
!       cur_tor_ps_soln is toroidal current due to perpendicular and
!       Pfirsch Schluter current.
!       cur_seed hold current density used to control q on axis if
!       q0_max > 0.0
      end module SOLN
