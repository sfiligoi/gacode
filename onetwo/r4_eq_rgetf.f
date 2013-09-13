!----------------------------------------------------------
! r8real -- generated (f90 fixed form) REAL interface to:  EQ_RGETF
      SUBROUTINE R4_EQ_RGETF(
     > IVEC,R4_ZRHO,IFCN,IWANT,R4_ZVAL,IERR)

      external EQ_RGETF

! argument declarations
      INTEGER IVEC
      INTEGER IFCN
      INTEGER IWANT
      INTEGER IERR
 ! floating type, input only:
      REAL R4_ZRHO(IVEC)
  
 ! floating type, output only:
      REAL R4_ZVAL(IVEC)
  

! local (automatic array) declarations

      integer alloc_stat
      REAL*8, allocatable :: ZRHO_act(:)
      REAL*8, allocatable :: ZVAL_act(:)

! allocation of working arrays...

      ALLOCATE(ZRHO_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RGETF -- ZRHO_act ALLOCATE error!')
      ALLOCATE(ZVAL_act(IVEC)
     >   ,stat=alloc_stat)
      if(alloc_stat.ne.0) call errmsg_exit(
     >  ' ?R4_EQ_RGETF -- ZVAL_act ALLOCATE error!')

! executable code:  copy for input
      ZRHO_act=R4_ZRHO
      ZVAL_act=0
! call to original routine:  EQ_RGETF

      CALL EQ_RGETF(
     > IVEC,ZRHO_act(1),IFCN,IWANT,ZVAL_act(1),IERR)

! copy back outputs if modified.
      DEALLOCATE(ZRHO_act)
      R4_ZVAL = ZVAL_act
      DEALLOCATE(ZVAL_act)

! exit
      return
      end
