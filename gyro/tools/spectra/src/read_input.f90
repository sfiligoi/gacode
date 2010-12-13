subroutine read_input(dir)

  use input_data

  implicit none

  character (len=*) :: dir

  integer :: i
  integer :: idum

  real :: dum
  real :: pi
  real, dimension(:), allocatable :: vec

  pi = 4.0*datan(1.0)

  open(unit=1,file=trim(dir)//'profile_vugyro.out')
  read(1,*) n_r

  allocate(r(n_r))
  allocate(vec(n_r))

  read(1,*) idum
  read(1,*) n_pass
  read(1,*) n_trap
  read(1,*) n_energy
  read(1,*) idum
  read(1,*) idum
  read(1,*) n_n

  read(1,*) idum
  read(1,*) idum
  read(1,*) idum
  read(1,*) idum
  read(1,*) n_field
  read(1,*) n_ion
  read(1,*) n_kinetic
  read(1,*) n_spec

  allocate(kx(0:n_r/2-1))
  allocate(ky(0:n_n-1))

  do i=1,4
     read(1,*) idum
  enddo
  read(1,*) r(:)
  do i=1,27+5*(n_spec-1)
     read(1,*) vec(:)
  enddo
  read(1,*) dum
  do i=1,11
     read(1,*) vec(:)
  enddo
  do i=1,n_pass+n_trap+n_energy+1
     read(1,*) dum
  enddo
  read(1,*) ky
  read(1,*) rho
  close(1)

  deallocate(vec)

  ! Define kx vector:
  do i=0,n_r/2-1
     kx(i) = i*rho*2*pi/(r(n_r)-r(1))
  enddo

end subroutine read_input
