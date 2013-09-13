  MODULE nf_param

      INTEGER,PARAMETER :: kj =51

      INTEGER,PARAMETER :: kprim = 3, kimp = 2, kion = kprim + kimp
      INTEGER,PARAMETER :: kk = kion + 4, kjm1 = kj-1, kbctim = 350
      INTEGER,PARAMETER :: kb = 32,ke = 3,le_mcgo = ke+1, kbe = kb*ke
      INTEGER,PARAMETER :: kcm = 3, kcmp1 = kcm + 1


      INTEGER,PARAMETER :: ksge = 20, kf = 81,kz =kf-1
      INTEGER,PARAMETER :: kzrf = 20, krt = 10 ,kcrm = 6
!      INTEGER,PARAMETER :: nap = 10 !  nap can be reduced to 4 I think
      INTEGER,PARAMETER :: nap = 4 ! nap =10 in onetwo
      INTEGER,PARAMETER :: kjeb = kj*ke*kb,kbs = 2,kedge = kj 
!
! --- see include file fusion.i for use of ksymbp (below)
!
      INTEGER,PARAMETER :: maxp = 3000, kjp = 301, ksplin = MAX(51,kj)
      INTEGER,PARAMETER :: kar = 4,kdi = 16 
      INTEGER,PARAMETER :: ksymbp = (3 * kb * (3 * kb + 1)) / 2 + 1 

      INTEGER,PARAMETER :: ktab = 25    ! formerly 51
!
   END MODULE nf_param
!
!
!
