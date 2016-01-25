!-----------------------------------------------------------------
! cgyro_mpi_grid.f90
!
! PURPOSE:
!  Subroutinized main cgyro program.  
!-----------------------------------------------------------------

subroutine cgyro_mpi_grid

  use timer_lib
  use parallel_lib
  use mpi

  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: ie,ix,is,ir,it
  integer :: d
  integer :: splitkey
  integer, external :: parallel_dim

  integer, external :: omp_get_max_threads, omp_get_thread_num

  ! Velocity-space (v) and configuration-space (c) dimensions
  nv = n_energy*n_xi*n_species
  nc = n_radial*n_theta

  call gcd(nv,nc,d)

  !-------------------------------------------------------------------------
  ! MPI diagnostics need to come early
  !
  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io,file=trim(path)//runfile_mpi,status='replace')
     write(io,*) 'Parallelization and distribution diagnostics'
     write(io,*)
     write(io,'(a,i5)') '         nv: ',nv
     write(io,'(a,i5)') '         nc: ',nc
     write(io,'(a,i5)') ' GCD(nv,nc): ',d
     write(io,*)
     write(io,*) '        [coll]    [str]     [NL]'
     write(io,*) 'n_MPI   nc_loc   nv_loc  n_split'
     write(io,*) '-----   ------   ------  -------'
     do it=1,d*n_toroidal
        if (mod(d*n_toroidal,it) == 0 .and. mod(it,n_toroidal) == 0) then
           n_proc_1 = it/n_toroidal
           nc_loc = nc/n_proc_1           
           nv_loc = nv/n_proc_1           
           write(io,'(t2,4(i5,4x))') it,nc_loc,nv_loc,1+(nv_loc*n_theta-1)/n_toroidal
        endif
     enddo
     close(io)
  endif
  !-------------------------------------------------------------------------

  allocate(ie_v(nv))
  allocate(ix_v(nv))
  allocate(is_v(nv))

  allocate(ir_c(nc))
  allocate(it_c(nc))

  allocate(ic_c(n_radial,n_theta))
  allocate(iv_v(n_energy,n_xi,n_species))

  if (test_flag == 1) return

  !-------------------------------------------------------------
  ! Check that n_proc is a multiple of n_toroidal
  !
  if (modulo(n_proc,n_toroidal) /= 0) then
     call cgyro_error('ERROR: (CGYRO) Number of processors must be a multiple of N_TOROIDAL.')
     return
  endif

  ! Assign subgroup dimensions:

  n_proc_2 = n_toroidal
  n_proc_1 = n_proc/n_toroidal

  ! Check that nv and nc are multiples of the local processor count

  if (modulo(nv,n_proc_1) /= 0 .or. modulo(nc,n_proc_1) /= 0) then
     call cgyro_error('ERROR: (CGYRO) nv or nc not a multiple of the local processor count.')
     return
  endif

  ! Local group indices:

  i_group_1 = i_proc/n_proc_1
  i_group_2 = modulo(i_proc,n_proc_1)
  !------------------------------------------------

  !-----------------------------------------------------------
  ! Split up GYRO_COMM_WORLD into groups and adjoint:
  !
  !             NEW_COMM_1  and  NEW_COMM_2
  !
  splitkey = i_proc
  call MPI_COMM_SPLIT(CGYRO_COMM_WORLD,&
       i_group_1,& 
       splitkey,&
       NEW_COMM_1, &
       i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) NEW_COMM_1 not created')
     return
  endif

  ! Local adjoint Group number

  call MPI_COMM_SPLIT(CGYRO_COMM_WORLD,&
       i_group_2,& 
       splitkey,&
       NEW_COMM_2, &
       i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) NEW_COMM_2 not created')
     return
  endif
  !
  call MPI_COMM_RANK(NEW_COMM_1,i_proc_1,i_err)
  call MPI_COMM_RANK(NEW_COMM_2,i_proc_2,i_err)
  !
  !-----------------------------------------------------------

  ! Velocity pointers
  iv = 0
  do ie=1,n_energy
     do ix=1,n_xi
        do is=1,n_species
           iv = iv+1
           ie_v(iv) = ie
           ix_v(iv) = ix
           is_v(iv) = is
           iv_v(ie,ix,is) = iv
        enddo
     enddo
  enddo

  ! Configuration pointers
  ic = 0
  do ir=1,n_radial
     do it=1,n_theta
        ic = ic+1
        ir_c(ic) = ir
        it_c(ic) = it
        ic_c(ir,it) = ic
     enddo
  enddo

  ! Linear parallelization dimensions

  ! ni -> nc
  ! nj -> nv  
  call parallel_lib_init(nc,nv,nc_loc,nv_loc,NEW_COMM_1)

  nv1 = 1+i_proc_1*nv_loc
  nv2 = (1+i_proc_1)*nv_loc

  nc1 = 1+i_proc_1*nc_loc
  nc2 = (1+i_proc_1)*nc_loc

  ! Nonlinear parallelization dimensions (returns nsplit)

  if (nonlinear_method == 1) then
     call parallel_slib_init(n_toroidal,nv_loc,nc,nsplit,NEW_COMM_2)
  else
     call parallel_slib_init(n_toroidal,nv_loc*n_theta,n_radial,nsplit,NEW_COMM_2)
  endif

  ! OMP code
  n_omp = omp_get_max_threads()

  allocate(ic_locv(nc1:nc2))
  allocate(iv_locv(nv1:nv2))
  
  ic_loc = 0
  do ic=nc1,nc2
     ic_loc = ic_loc+1
     ic_locv(ic) = ic_loc
  enddo
  iv_loc = 0
  do iv=nv1,nv2
     iv_loc = iv_loc+1
     iv_locv(iv) = iv_loc
  enddo
  
end subroutine cgyro_mpi_grid

subroutine gcd(m,n,d)

  implicit none

  integer, intent(in) :: m,n
  integer, intent(inout) :: d
  integer :: a,b,c

  a = m
  b = n

  if (a < b) then
     c = a
     a = b
     b = c
  endif

  do          
     c = mod(a, b)    
     if (c == 0) exit
     a = b         
     b = c 
  enddo

  d = b

end subroutine gcd
