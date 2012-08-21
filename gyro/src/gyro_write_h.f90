!------------------------------------------------------
! gyro_write_h.f90 [caller gyro_write_master]
!
! PURPOSE:
!  Distribution function output for all kinetic species.
!-----------------------------------------------------

  subroutine gyro_write_h(datafile1,datafile2,io1,io2)

  use gyro_globals
  use gyro_pointers

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: io1
  character (len=*), intent(in) :: datafile1
  integer, intent(in) :: io2
  character (len=*), intent(in) :: datafile2



  !---------------------------------------------------

  select case (io_control)

  case(0)

     return

  case(1)

     open(unit=io1,file=trim(path)//'theta_m.out',status='replace')
     print *,'Writing to theta_m.out'
     do k=1,n_lambda
        ck = class(k) 
        do j=1,n_theta(ck)
           write(io1,*) theta_t(1,k,j)
        enddo
     enddo
     do k=1,n_lambda
        ck = class(k) 
        do j=1,n_theta(ck)
           write(io1,*) sqrt(abs(1.0-lambda(ir_norm,k)*b0_t(ir_norm,k,j)))
        enddo
     enddo
     close(io1)

     open(unit=io1,file=datafile1,status='replace')
     print *,'Opened ',datafile1
     close(io1)

     open(unit=io2,file=datafile2,status='replace')
     print *,'Opened ',datafile2
     close(io2)

  case(2)

     open(unit=io1,file=datafile1,status='old',position='append')
     open(unit=io2,file=datafile2,status='old',position='append')

     do is=1,n_kinetic

        p_nek_loc = 0

        do p_nek=1+i_proc_1,n_nek_1,n_proc_1

           p_nek_loc = p_nek_loc+1

           ie = nek_e(p_nek)  
           k  = nek_k(p_nek)   

           ck = class(k)

           if (ie == n_energy/2) then

              if (k <= n_pass) then

                 do j=1,n_theta(ck)
                    write(io1,fmtstr) real(h(m_map(ck,j,:),ir_norm,p_nek_loc,is))
                 enddo

              else

                 do j=1,n_theta(ck)
                    write(io2,fmtstr) real(h(m_map(ck,j,:),ir_norm,p_nek_loc,is))
                 enddo

              endif

           endif

        enddo

     enddo ! is

     close(io1)
     close(io2)

  case(3)

     call catch_error('ERROR: Restart not supported in write_h')

  end select

  if (debug_flag == 1) print *,'[write_h called]'

end subroutine gyro_write_h
