!-----------------------------------------------------------
! tgyro_restart.f90
!
! PURPOSE:
!  Manage restart.
!----------------------------------------------------------

subroutine tgyro_restart

  use mpi
  use tgyro_globals

  implicit none

  integer :: i
  integer :: ip
  integer :: j
  integer :: ioerr
  real, dimension(9) :: x_read
  character(len=1) :: dummy


  if (i_proc_global == 0) then

     open(unit=1,&
          file='control.out',&
          status='old',&
          iostat=ioerr)

     if (ioerr /= 0) then

        close(1)
        loc_restart_flag = 0
        stop

     else

        read(1,*) n_r
        read(1,*) n_evolve
        read(1,*) i_tran
        close(1)

        ! Need to read both gradients and relaxation values:

        open(unit=1,file='gradient.out',status='old')
        do j=0,i_tran
           read(1,'(a)') dummy
           read(1,'(a)') dummy
           do i=1,n_r
              read(1,*) x_read(1:6)
              dlnnidr(1,i) = x_read(2)/r_min
              dlnnedr(i)   = x_read(3)/r_min
              dlntidr(1,i) = x_read(4)/r_min
              dlntedr(i)   = x_read(5)/r_min
           enddo
        enddo
        close(1)

        if (loc_ne_feedback_flag == 0) then

           open(unit=1,file='flux_target.out',status='old')
           do j=0,i_tran
              read(1,'(a)') dummy
              read(1,'(a)') dummy
              do i=1,n_r
                 read(1,*) x_read(1:5)
                 eflux_i_tot(i) = x_read(2)
                 eflux_e_tot(i) = x_read(4)
              enddo
           enddo
           close(1)

        else

           open(unit=1,file='flux_target.out',status='old')
           do j=0,i_tran
              read(1,'(a)') dummy
              read(1,'(a)') dummy
              do i=1,n_r
                 read(1,*) x_read(1:7)
                 eflux_i_tot(i) = x_read(2)
                 eflux_e_tot(i) = x_read(4)
                 pflux_e_tot(i) = x_read(6)
              enddo
           enddo
           close(1)

        endif

        open(unit=1,file='residual.out',status='old')
        do j=0,i_tran 
           read(1,'(a)') dummy
           do i=2,n_r
              read(1,*) &
                   x_read(1),(res(pmap(i,ip)),relax(pmap(i,ip)),ip=1,n_evolve)
           enddo
        enddo
        res = res*size(res)
        close(1)

     endif

  endif

  ! Reset n_evolve in case it was changed from value stored in control.out.

  n_evolve = &
       loc_ti_feedback_flag+&
       loc_te_feedback_flag+&
       loc_ne_feedback_flag+&
       loc_er_feedback_flag

  call MPI_BCAST(loc_restart_flag,&
       1,&
       MPI_INTEGER,&
       0,&
       MPI_COMM_WORLD,&
       ierr)

  call MPI_BCAST(i_tran,&
       1,&
       MPI_INTEGER,&
       0,&
       MPI_COMM_WORLD,&
       ierr)

  if (loc_restart_flag == 1) then

     call MPI_BCAST(dlnnidr,&
          size(dlnnidr),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(dlnnedr,&
          size(dlnnedr),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(dlntidr,&
          size(dlntidr),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(dlntedr,&
          size(dlntedr),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(eflux_i_tot,&
          size(eflux_i_tot),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(eflux_e_tot,&
          size(eflux_e_tot),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(pflux_e_tot,&
          size(pflux_e_tot),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(relax,&
          size(relax),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

  endif

end subroutine tgyro_restart
