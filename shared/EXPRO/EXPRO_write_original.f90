!--------------------------------------------------------
! EXPRO_write_original.f90
!
! PURPOSE:
!  Rewrite original quantities to input.profiles
!
!  NOTE:
!   q is kept positive definite.  All information
!   about sign of q is retained in IPCCW and BTCCW.  
!--------------------------------------------------------

subroutine EXPRO_write_original(tag)

  use EXPRO_globals
  use EXPRO_interface

  implicit none

  integer :: i
  integer :: ierr

  character(len=*), intent(in) :: tag
  character (len=80) :: line

  open(unit=1,file=trim(path)//'input.profiles',status='old')
  open(unit=2,file=trim(path)//'input.profiles.new',status='replace')

  write(2,'(a)') '# '//tag
  write(2,'(a)') '# '
  do 

     read(1,'(a)',iostat=ierr) line
     if (ierr < 0) exit
     write(2,'(a)') line

     if (line(2:4) == 'rho') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') &
                EXPRO_rho(i),&
                EXPRO_rmin(i),&
                EXPRO_poloidalfluxover2pi(i),&
                EXPRO_q(i),&
                EXPRO_w0(i)
        enddo
     endif

     if (line(2:5) == 'rmaj') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') &
                EXPRO_rmaj(i),&
                EXPRO_zmag(i),&
                EXPRO_kappa(i),&
                EXPRO_delta(i),&
                EXPRO_zeta(i)
        enddo
     endif

     if (line(2:3) == 'ne') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') &
                EXPRO_ne(i),&
                EXPRO_te(i),&
                EXPRO_ptot(i),&
                EXPRO_z_eff(i),&
                0.0
        enddo
     endif

     if (line(2:5) == 'ni_1') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') EXPRO_ni(1:5,i)
        enddo
     endif
     if (line(2:5) == 'ni_6') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') EXPRO_ni(6:10,i)
        enddo
     endif

     if (line(2:5) == 'Ti_1') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') EXPRO_ti(1:5,i)
        enddo
     endif
     if (line(2:5) == 'Ti_6') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') EXPRO_ti(6:10,i)
        enddo
     endif

     if (line(2:7) == 'vtor_1') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') EXPRO_vtor(1:5,i)
        enddo
     endif
     if (line(2:7) == 'vtor_6') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') EXPRO_vtor(6:10,i)
        enddo
     endif

     if (line(2:7) == 'vpol_1') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') EXPRO_vpol(1:5,i)
        enddo
     endif
     if (line(2:7) == 'vpol_6') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') EXPRO_vpol(6:10,i)
        enddo
     endif

     if (line(2:10) == 'flow_beam') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') &
                EXPRO_flow_beam(i),&
                EXPRO_flow_wall(i),&
                EXPRO_flow_mom(i),&
                0.0,&
                0.0
        enddo
     endif

     if (line(2:7) == 'pow_e(') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') &
                EXPRO_pow_e(i),&
                EXPRO_pow_i(i),&
                EXPRO_pow_ei(i),&
                EXPRO_pow_e_aux(i),&
                EXPRO_pow_i_aux(i)
        enddo
     endif

     if (line(2:10) == 'pow_e_fus') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') &
                EXPRO_pow_e_fus(i),&
                EXPRO_pow_i_fus(i),&
                EXPRO_pow_e_sync(i),&
                EXPRO_pow_e_brem(i),&
                EXPRO_pow_e_line(i)
        enddo
     endif

  enddo

  close(1)
  close(2)

end subroutine EXPRO_write_original
