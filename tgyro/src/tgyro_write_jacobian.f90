subroutine tgyro_write_jacobian(particle)
  use tgyro_globals
  use tgyro_iteration_variables
  implicit none

  character(len=1), intent(in):: particle
     open(unit=101, file='jac-'//trim(particle)//'.out')
     open(unit=102, file='Q-z-'//trim(particle)//'.out')
     do p=1,p_max,n_evolve
        do pp=0,n_evolve-1
           write(101,*) p+pp, p+ip, jf(p+pp,p+ip)
           write(102,*) p+pp, p+ip, x_vec(p+ip), f_vec(p+pp) 
        enddo
     enddo
     close(101)
     close(102)

end subroutine tgyro_write_jacobian     

