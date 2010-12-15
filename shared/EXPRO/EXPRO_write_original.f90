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

subroutine EXPRO_write_original(path,tag)

  use EXPRO_interface

  implicit none

  integer :: i
  integer :: ierr

  character(len=*), intent(in) :: tag
  character(len=*) :: path
  character (len=80) :: line

  open(unit=1,file=trim(path)//'input.profiles',status='old')
  open(unit=2,file=trim(path)//'input.profiles.new',status='replace')

  write(2,'(a)') '# '//tag
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
                EXPRO_rmaj(i),&
                abs(EXPRO_q(i)),&
                EXPRO_kappa(i)
        enddo
     endif

     if (line(2:6) == 'delta') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') &
                EXPRO_delta(i),&
                EXPRO_te(i),&
                EXPRO_ne(i),&
                EXPRO_z_eff(i),&
                EXPRO_w0(i)
        enddo
     endif

     if (line(2:9) == 'flow_mom') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') &
                EXPRO_flow_mom(i),&
                EXPRO_pow_e(i),&
                EXPRO_pow_i(i),&
                EXPRO_pow_ei(i),&
                EXPRO_zeta(i)
        enddo
     endif

     if (line(2:10) == 'flow_beam') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') &
                EXPRO_flow_beam(i),&
                EXPRO_flow_wall(i),&
                EXPRO_zmag(i),&
                0.0,&
                0.0
        enddo
     endif

     if (line(2:3) == 'ni') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') EXPRO_ni(:,i)
        enddo
     endif

     if (line(2:3) == 'Ti') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') EXPRO_ti(:,i)
        enddo
     endif

    if (line(2:5) == 'vtor') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') EXPRO_vtor(:,i)
        enddo
     endif

    if (line(2:5) == 'vpol') then
        do i=1,EXPRO_n_exp
           read(1,'(a)',iostat=ierr) line
           write(2,'(5(1pe14.7,2x))') EXPRO_vpol(:,i)
        enddo
     endif
  enddo

  close(1)
  close(2)

end subroutine EXPRO_write_original
