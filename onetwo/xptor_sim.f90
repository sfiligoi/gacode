MODULE xptor_sim
! stores info for for simulation of xptor runs:

       real *8,dimension(:),allocatable :: rxpt,texpt,tixpt
       real *8,dimension(:),save,allocatable :: qbeame_fixed 
       real *8,dimension(:),save,allocatable :: qbeami_fixed
       real *8,dimension(:),save,allocatable :: qrfe_fixed 
       real *8,dimension(:),save,allocatable :: qrfi_fixed
       integer    read_qbeami,read_qbeame,test_xptor
       integer    read_qrfe,read_qrfi

END MODULE
