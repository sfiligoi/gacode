subroutine tgyro_write_jacobian

  use tgyro_globals
  use tgyro_iteration_variables

  implicit none

  character(len=3), dimension(n_evolve) :: vary
  character(len=3), dimension(n_evolve) :: varx
  integer :: ipt

  open(unit=101,file='out.tgyro.jacobian',status='replace')

  ipt = 0
  if (loc_ti_feedback_flag == 1) then
     ipt = ipt+1
     vary(ipt) = 'Ti'
     varx(ipt) = 'zTi'
  endif
  if (loc_te_feedback_flag == 1) then
     ipt = ipt+1
     vary(ipt) = 'Te'
     varx(ipt) = 'zTe'
  endif
  if (loc_ne_feedback_flag == 1) then
     ipt = ipt+1
     vary(ipt) = 'ne'
     varx(ipt) = 'zne'
  endif
  if (loc_er_feedback_flag == 1) then
     ipt = ipt+1
     vary(ipt) = 'Er'
     varx(ipt) = 'w0'
  endif

  do p=1,p_max,n_evolve
     do pp=0,n_evolve-1
        write(101,10) imap(p,ip),vary(ip),varx(pp),jf(p+pp,p+ip)
     enddo
  enddo

  close(101)

10 format(i2,a,a,1pe12.5)

end subroutine tgyro_write_jacobian

