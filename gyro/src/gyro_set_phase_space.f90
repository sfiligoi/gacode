!-----------------------------------------------------
! gyro_set_phase_space.f90
!
! PURPOSE:
!  Combine the results of make_energy_grid and 
!  make_lambda_grid -- to make phase-space weight
!  w_p(ie,k).
!
!  Also, print various weights if verbose = 1.
!-------------------------------------------------------

subroutine gyro_set_phase_space(datafile,io)

  use gyro_globals
  use math_constants

  implicit none

  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile

  character (len=8) :: species_tag

  do ie=1,n_energy
     do i=1,n_x
        do k=1,n_lambda
           ck = class(k)
           w_p(ie,i,k) = w_energy(ie)*w_lambda(i,k)*d_tau(ck)
        enddo
     enddo
  enddo

  if (i_proc == 0 .and. output_flag == 1) then

     open(unit=io,file=datafile,status='replace')

     ! Energy grid output

     write(io,'(a)') ' i    energy(i)    w_energy(i) '
     write(io,'(a)') '--- ------------  ------------ '
     do ie=1,n_energy
        write(io,'(i2,3(2x,f12.10))') ie,energy(ie),w_energy(ie)
     enddo

     ! Lambda-grid output

     do i=1,n_x,n_x-1

        write(io,*)
        write(io,'(a)') ' k   tau(k)/2pi     s_lambda(k)   lambda(k)   w_lambda(k)'
        write(io,'(a)') '--- ------------  ------------ ------------  -----------'
        do k=1,n_lambda
           write(io,'(i2,4(2x,f12.10))') & 
                k,1.0/(pi_2*omega(i,k)),s_lambda(k),lambda(i,k),w_lambda(i,k)
           if (k == n_pass) write(io,'(a,t33,f12.10)') 'lambda-TP:',lambda_tp(i)
           if (k == n_lambda) write(io,'(a,t33,f12.10)') 'lambda-MAX:',lambda_max(i) 
        enddo

     enddo

     close(io)

  endif

end subroutine gyro_set_phase_space
