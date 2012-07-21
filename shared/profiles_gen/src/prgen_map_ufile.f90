!------------------------------------------------------------
! prgen_map_iterdb.f90
!
! PURPOSE:
!  Map native iterdb data onto input.profiles standard.  
!  
! NOTES:
!  - This routine is common to both text and NetCDF formats.
!  - See map_plasmastate.f90 for analogous routine for 
!    plasmastate data.
!------------------------------------------------------------

subroutine prgen_map_ufile

  use prgen_read_globals

  implicit none

  integer :: i
  integer :: nion=1
  integer :: nb=1

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  allocate(vec(n_indx,ufile_nj))
  vec(:,:) = 0.0
  !
  vec(1,:)  = rho(:)
  vec(2,:)  = rmin(:)
  vec(3,:)  = 0.0
  vec(4,:)  = q(:)
  vec(5,:)  = kappa(:)
  vec(6,:)  = delta(:)
  vec(7,:)  = ufile_te(:)
  vec(8,:)  = 0.0
  vec(9,:)  = 0.0
  vec(11,:) = 0.0
  vec(12,:) = 0.0
  vec(13,:) = 0.0
  vec(14,:) = 0.0
  vec(15,:) = 0.0
  vec(16,:) = 0.0
  vec(17,:) = 0.0
  vec(18,:) = 0.0
  vec(19,:) = 0.0
  vec(20,:) = 0.0 

  !-----------------------------------------------------------------
  ! Construct ion densities and temperatures with reordering
  ! in general case.  Use vphi and vpol as temporary arrays.
  !
  
  do i=1,nion
     ! ni
     vec(31+i-1,:) = 0.0
     ! Ti
     vec(36+i-1,:) = 0.0
  enddo
  ! Beam ions
  do i=1,nb
     ! ni
     vec(31+i+nion-1,:) = 0.0
     ! Ti: T[keV] = (p/n)[J]/1.6022e-16[J/eV]
     vec(36+i+nion-1,:) = 0.0
  enddo

  ! reorder
  do i=1,5 
     vec(20+i,:) = vec(30+reorder_vec(i),:)
     vec(25+i,:) = vec(35+reorder_vec(i),:)
  enddo

  ! vphi
  vec(31:35,:) = 0.0

  ! Insert carbon toroidal velocity
  do i=1,5
     if (reorder_vec(i) == nion+1) then
        vec(30+i,:) = 0.0
     endif
  enddo

  ! vpol
  vec(36,:) = 0.0
  vec(37,:) = 0.0
  vec(38,:) = 0.0
  vec(39,:) = 0.0
  vec(40,:) = 0.0

end subroutine prgen_map_ufile
