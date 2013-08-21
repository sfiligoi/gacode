  MODULE param
! --- INCLUDE file param.i
!
! --- kk is now kion + 4 instead of kion + 3,
! --- to allow for angular momentum equation
! --- kbctim changed from 50 to 100 7/25/96 HSJ
! --- kbctim changed to 150 HSJ
! --- ksplin changed from 21 to  51 12/8/97 HSJ
!
       integer,parameter :: kj =201,ktoray=kj,ktorm1 = ktoray-1
       integer,parameter :: nwmhd =129,nhmhd =129 
!
!                         !!!!!!!!!! In the toray code must change kj 
!                         !!!!!!!!!! In the gafit code must change njmax
!
!
!
!
      integer,parameter :: kprim = 3, kimp = 2, kion = kprim + kimp
      integer,parameter :: kk = kion + 4, kjm1 = kj-1, kbctim = 350
      integer,parameter :: kb = 32, ke = 3,le_mcgo = ke+1, kbe = kb*ke
! changed kb from 32 to 14 2/28/13 to be consistent with new beamlets 
!     integer,parameter :: kb = 14, ke = 3,le_mcgo = ke+1, kbe = kb*ke
      integer,parameter :: kcm = 3, kcmp1 = kcm + 1
      integer,parameter :: ksge = 20, kf = 81, kz = kf-1, krf = 30 
      integer,parameter :: kzrf = 20, krt = 10 ,kcrm = 6
      integer,parameter :: kevents = 10 + 2*krf + 7 + 2*kb
! add 2*kb to kevents t account for 32 beamlet on/off swithcing times for 
! P_Nfreya
      integer,parameter :: nap = 4 ! old  nap =10  can be reduced to 4 I think
      integer,parameter :: kjeb = kj*ke*kb,kbs = 2,kedge = kj 
!
! --- see include file fusion.i for use of ksymbp (below)
!
      integer,parameter :: maxp = 3000, kjp = 301, ksplin = MAX(51,kj)
      integer,parameter :: kar = 4,kdi = 16 
      integer,parameter :: ksymbp = (3 * kb * (3 * kb + 1)) / 2 + 1 
!
!       parameter (maxp = 300, kjp = 301, ksplin = 51,kar = 4,
!     .           kdi = 16, ksymbp = (3 * kb * (3 * kb + 1)) / 2 + 1 )
!
! --- parameter  ktab was pulled out of subroutine MIX and inserted here
!
      integer,parameter :: ktab = 25    ! formerly 51
!
   END MODULE param
!
!
!
