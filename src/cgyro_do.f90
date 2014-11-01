!-----------------------------------------------------------------
! cgyro_do.f90
!
! PURPOSE:
!  Subroutinized main cgyro program.  
!-----------------------------------------------------------------

subroutine cgyro_do

  use timer_lib
  use mpi

  use cgyro_globals

  implicit none

  integer :: lc,lv,iv,ic,nc,nv
  integer :: ie,ix,is,ir,it
  integer, dimension(:), allocatable :: ie_v
  integer, dimension(:), allocatable :: ix_v
  integer, dimension(:), allocatable :: is_v
  integer, dimension(:), allocatable :: ir_c
  integer, dimension(:), allocatable :: it_c
  real, dimension(:,:), allocatable :: h

  n_energy  = 3
  n_xi      = 3
  n_species = 2
  n_radial  = 1
  n_theta   = 2

  nv = n_energy*n_xi*n_species
  nc = n_radial*n_theta

  allocate(ie_v(nv))
  allocate(ix_v(nv))
  allocate(is_v(nv))

  allocate(ir_c(nc))
  allocate(it_c(nc))

  lv = 0
  do ie=1,n_energy
     do ix=1,n_xi
        do is=1,n_species
           lv = lv+1
           ie_v(lv) = ie
           ix_v(lv) = ix
           is_v(lv) = is
        enddo
     enddo
  enddo

  lc = 0
  do ir=1,n_radial
     do it=1,n_theta
        lc = lc+1
        ir_c(lc) = ir
        it_c(lc) = it
     enddo
  enddo

  nv_loc = 

  allocate(h(nv,nc))

  do iv=1+i_proc,nvmn_proc
     do ic=1,nc

        h(iv,ic) = ie_v(iv)*10*0 + &
             ie_v(iv)*10**0 + &
             ix_v(iv)*10**1 + &
             is_v(iv)*10**2 + &
             ir_c(ic)*10**3 + &
             it_c(ic)*10**4

     enddo
  enddo



  print *,h

end subroutine cgyro_do
