!ubroutine index_convert_fi(is,isb)
!
  implicit NONE
  include 'TRCOM'
!
  !  mapping from NUBEAM fast ion index to TRANSP fast ion index
!
  integer,intent(in) :: is
  integer,intent(out) :: isb
!
  integer ia,iz,itype,isi
  integer index_nbfi         ! NUBEAM fast ion index function
!
  !------------------------------------------------------
  !  internal state data
!
  integer, save :: index_table(mibs)
  real, save :: ztime_ref = -1.01e30
!
  !------------------------------------------------------
!
  if(ztime_ref.ne.ta) then
!
     !  update table if time has changed since last call
!
     index_table=0
     ztime_ref=ta
!
     do isi=1,nsbeam
!
        ia = abeams(isi)+0.5
        iz = xzbeams(isi)+0.5
!
        if(nlfprod(isi)) then
           !  ... fusion product
           if(nlfbmfpp(isi)) then
              itype=4  !  FPP
           else
              itype=3  !  Monte Carlo
           endif
        else
           !  ... beam ion
           if(nlfbmfpp(isi)) then
              itype=2  !  FPP
           else
              itype=1  !  Monte Carlo
           endif
        endif
!
        index_table(isi) = index_nbfi(iz,ia,itype)
     enddo
!
  endif
!
  isb=index_table(is)
!
!nd subroutine index_convert_fi
