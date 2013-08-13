      MODULE allocate_err
!
!
      CONTAINS
!
!
      subroutine allocate_error(var,myid,istat)
      character *(*) var
      integer istat,myid
      write(*,1)var,istat,myid
 1    format(2x,"Memory Allocation error encountered",/, &
            2x,"Unable to deallocate ",a,/, &
            2x,"status =",i5) 
      istat =0 !reset for next case
      return
      end      subroutine allocate_error
!
!
!
      subroutine deallocate_error(var,myid,istat)
      character *(*) var
      integer istat, myid
      write(*,1)var,istat,myid
 1    format(2x,"Memory DE-Allocation error encountered",/,&
            2x,"Unable to deallocate ",a,/,&
            2x,"status =",i5, " process rank =",i5)
      istat =0 !reset for next case
      return
      end      subroutine deallocate_error
!
!
!
      END MODULE allocate_err
