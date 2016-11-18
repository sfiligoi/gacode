subroutine tgyro_write_jacobian(init)

  use tgyro_globals
  use tgyro_iteration_variables

  implicit none

  integer, intent(in) :: init
  character(len=4), dimension(n_evolve) :: vary
  character(len=4), dimension(n_evolve) :: varx
  integer :: ipt

  if (i_proc_global > 0) return

  if (init == 1) then

     open(unit=101,file='out.tgyro.jacobian',status='replace')

  else

     open(unit=101,file='out.tgyro.jacobian',status='old',position='append')

     ipt = 0
     if (loc_ti_feedback_flag == 1) then
        ipt = ipt+1
        vary(ipt) = 'dQi/'
        varx(ipt) = 'dzTi'
     endif
     if (loc_te_feedback_flag == 1) then
        ipt = ipt+1
        vary(ipt) = 'dQe/'
        varx(ipt) = 'dzTe'
     endif
     if (loc_er_feedback_flag == 1) then
        ipt = ipt+1
        vary(ipt) = 'dPi/'
        varx(ipt) = 'dw0'
     endif

     do p=1,p_max,n_evolve
        write(101,'(a,1pe13.6)') 'r/a = ',r((p-1)/n_evolve+2)/r_min
        do pp=0,n_evolve-1
           write(101,10) vary(ip+1),varx(pp+1),jf(p+pp,p+ip)
        enddo
     enddo

  endif

  close(101)

10 format(a,a,1x,1pe13.6)

end subroutine tgyro_write_jacobian

