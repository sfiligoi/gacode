!------------------------------------------------------
! gyro_write_efficiency.f90 [caller BigScience]
!
! PURPOSE:
!  This subroutine computes and prints the parallel
!  distribution efficiency.
!------------------------------------------------------

subroutine gyro_write_efficiency(datafile,io)

  use gyro_globals

  !--------------------------------------------------
  implicit none
  !
  integer :: msplit
  integer :: nv2
  integer :: nlit
  integer, external :: parallel_dim
  real :: x_s
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !--------------------------------------------------

  if (output_flag /= 1) return

  if (i_proc == 0) then

     open(unit=io,file=datafile,status='replace')

     write(io,*) '---------- DISTRIBUTION EFFICIENCY -------------'
     write(io,*) 'N_proc  NL iter  Eff(%)'
     write(io,*) '------  -------  ------'
     do i=n_n,n_n*n_energy*n_lambda,n_n

        nv2    = parallel_dim(n_energy*n_lambda,i/n_n)
        msplit = 1+(n_stack*nv2-1)/n_n
        nlit   = msplit*n_kinetic  
        x_s = (n_stack*n_energy*n_lambda*n_n/real(i))/(msplit*n_n)
        if (x_s >= 0.99999) then
           write(io,10) i,nlit,x_s*100,'***'
        else if (x_s >= 0.9) then
           write(io,10) i,nlit,x_s*100,'**'
        else if (x_s >= 0.8) then
           write(io,10) i,nlit,x_s*100,'*'
        else 
           write(io,10) i,nlit,x_s*100
        endif
     enddo ! i

     close(io)

  endif

10 format(t2,i6,t12,i4,t19,f6.2,1x,a)

end subroutine gyro_write_efficiency
