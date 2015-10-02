subroutine cgyro_init_h

  use mpi
  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: ir,it,is,ie,ix
  real :: arg, ang

  !---------------------------------------------------------------------
  ! Check to see if we have restart data available
  !
  if (i_proc == 0) then

     open(unit=io,&
          file=trim(path)//runfile_restart_tag,&
          status='old',iostat=i_err)
     close(io)

     if (i_err == 0) then
        restart_flag = 1
     else
        restart_flag = 0
     endif

  endif

  call MPI_BCAST(restart_flag,1,MPI_DOUBLE_PRECISION,0,CGYRO_COMM_WORLD,i_err)
  !---------------------------------------------------------------------

  if (restart_flag == 1) then

     call cgyro_info('Restart data found.')
     call cgyro_read_restart

  else

     i_current = 0
     t_current = 0.0

     !-------------------------------------------------------------------------
     ! Generate analytic initial conditions
     !
     h_x(:,:) = (0.0,0.0)
     !
     iv_loc = 0
     do iv=nv1,nv2

        iv_loc = iv_loc+1

        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)

        do ic=1,nc

           ir = ir_c(ic) 
           it = it_c(ic)

           if (n == 0) then

              ! Zonal-flow initial condition

              if (zf_test_flag == 1) then
                 if (is == 1 .and. abs(px(ir)) == 1) then
                    h_x(ic,iv_loc) = 1e-6
                 endif
              else
                 ! CAUTION: Need f(p) = conjg[ f(-p) ] for n=0
                 arg = abs(px(ir))/real(n_radial)
                 h_x(ic,iv_loc) = amp*rho*exp(-arg)
                 if (ir == 1) h_x(ic,iv_loc) = (0.0,0.0)
                 ! Safer to zero out completely
                 h_x(ic,iv_loc) = 0.0
              endif

           else 

              ! Exponential in ballooning angle.

              if (n_toroidal == 1) then
                 if (is == 1) then
                    ang = theta(it)+2*pi*px(ir)
                    h_x(ic,iv_loc) = rho*exp(-(ang/2)**2) 
                 endif
              else
                 h_x(ic,iv_loc) = amp*rho
              endif

           endif

        enddo
     enddo

  endif

  call cgyro_field_c

  ! Initialize time-history of fields (-3,-2,-1) to initial field.
  field_old  = field
  field_old2 = field
  field_old3 = field

end subroutine cgyro_init_h
