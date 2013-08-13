      subroutine eq_rgetf(ivec,zrho,ifcn,iwant,zval,ierr)
c
      use eq_xyz_cmn
c
c  f(rho) interpolation routine -- with error checking
c
      IMPLICIT NONE
c
c  input:
c
      integer ivec                      ! vector dimension
c
      REAL*8 zrho(ivec)                 ! argument
c
c  input/output:
      integer ifcn                      ! function id number
c
c  input:
      integer iwant                     ! 0:  value, 1: df/drho, 2: d2f/drho2
c
c  output:
      REAL*8 zval(ivec)                 ! result of interpolation
c
      integer ierr                      ! =0: OK, =1: zrho out of range
c
c$r8real_input:  zrho
c$r8real_output:  zval
c---------------------------
c
c  the caller provides the name of the function desired, but it also
c  provides a number which is writable.  If the caller saves this number
c  and repeats the call, it saves having to lookup the function name,
c  a performance enhancement.
c
c---------------------------
c
      integer ict(3,0:2)
      integer i
c
      data ict/
     >   1,0,0,
     >   0,1,0,
     >   0,0,1/
c
c---------------------------
c
      ierr=0
c
      if((ifcn.le.0).or.(ifcn.gt.nx_items)) then
         call eq_errmsg(' ?eq_rgetf:  invalid function id.')
         write(lunerr,*) '  ifcn=',ifcn
         ierr=1
         return
      endif
c
      call eqi_frho(ivec,zrho,ifcn,ict(1,iwant),ivec,zval,ierr)
c
      if(ierr.ne.0) then
         write(lunerr,*) ' %eq_rgetf: error flag set, ',ierr,
     >      ' points out of range?'
      endif
c
      return
      end
