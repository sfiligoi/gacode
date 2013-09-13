      MODULE SOLCON
!
      USe param, only: krf
      implicit none
!
      save
      character*8,public ::    istep
      real *8,public  :: q0_max,q0_radius,q0_mult
      real *8,public  :: cur_seed_amp,cur_seed_ctr,cur_seed_width
      real *8,public  :: time0, time, timnew, timmax, dt, dtt
      real *8,public  :: dtmin, dtmax, dtevmin, relmin, relmax, relit
      real *8,public  :: theta, dvres, timzen, erneumax,steady_state
      real *8,public  :: eqtim0, relaxsol, timav, dtsumi,cximult
      real *8,public  :: timecrit, dtmaxcrit,cparam,dcparam_in  
      real *8,public ::  cpu_time
      real *8,public  :: time_tol
!
      integer *4,public :: n, nmax, itmax, irfsw1, isecrem0
      integer *4,public :: istop, icallneu, icallnub, icallrf
      integer *4,public :: ilimdt, idterr, ifreya_old
      integer *4,public :: ineucg, ifreya, irfcnt(krf)
      integer *4,public :: glf_debug,external_beam_cur
      logical,public ::  te_range_check,write_profiles
      integer,parameter,private :: nbp0_ic = 2
      character(len=8) ,public ::  bp0_icv(nbp0_ic),bp0_ic
      data time_tol /1.e-9/
      data bp0_icv /'eqdsk','analytic'/
      end module SOLCON
