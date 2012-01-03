!------------------------------------------------------
! gyro_write_radial_op.f90
!
! PURPOSE:
!  Manage output for selected radial operators.
!------------------------------------------------------

subroutine gyro_write_radial_op(datafile,io)

  use gyro_globals

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !---------------------------------------------------

  select case (output_flag)

  case (1)

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='replace')

        write(io,*) ' i              w_d1               s_d1               w_d2' 
        write(io,*)
        do i_diff=-m_dx,m_dx-i_dx
           write(io,'(i3,8(2x,f11.5))') &
                i_diff,w_d1(i_diff)*d_x,s_d1(i_diff)*d_x,w_d2(i_diff)*d_x*d_x
        enddo

        if (n_1(in_1) == 0) then
           write(io,*)
           write(io,*) 'Trace of double-gyroaverage matrix:'
           write(io,*)
           do is=1,n_gk
              do i_diff=-m_gyro,m_gyro-i_gyro
                 write(io,'(i3,2(2x,f8.5))') i_diff,gyro_trace(is,i_diff)
              enddo
              write(io,*)
              write(io,'(a,f8.5)') 'Sum: ',sum(gyro_trace(is,:))
           enddo ! is
        else
           write(io,*)
           write(io,*) 'No trace available since n>0.'
           write(io,*)
        endif

        close(io)

     endif

  end select

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_write_radial_op called]'
  endif

end subroutine gyro_write_radial_op
