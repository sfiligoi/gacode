!---------------------------------------------------------
! glf23.f90
!
! PURPOSE:
!  Manage standalone GLF23 run by calling
!
!    call glf23_read_input()
!    call glf23_run() 
!
!  and then adding extra I/O. 
!---------------------------------------------------------

program glf23

  use glf23_interface

  implicit none

  integer :: i
  integer :: n
  character (len=4) :: tag(2)=(/'ion1','ion2'/)
  real :: prec

  call glf23_read_input()
  call glf23_dump_local()
  call glf23_run()
  call glf23_dump_module() 

  if (glf23_use_transport_model_in) then

     ! Output to screen

     print 20,'Gam/Gam_GB','    Q/Q_GB','  Pi/Pi_GB', '    S/S_GB'
     print 10,'elec',&
          glf23_elec_pflux_out,&
          glf23_elec_eflux_out,&
          glf23_elec_mflux_out,&
          glf23_elec_expwd_out

     prec = abs(glf23_elec_pflux_out)+&
          abs(glf23_elec_eflux_out)+&
          abs(glf23_elec_mflux_out)+&
          abs(glf23_elec_expwd_out)

     do i=1,glf23_ns_in-1
        print 10,tag(i),&
             glf23_ion_pflux_out(i),&
             glf23_ion_eflux_out(i),&
             glf23_ion_mflux_out(i),&
             glf23_ion_expwd_out(i)

        prec = prec+&
             abs(glf23_ion_pflux_out(i))+&
             abs(glf23_ion_eflux_out(i))+&
             abs(glf23_ion_mflux_out(i))
     enddo

     ! Output to file

     n = glf23_ns_in-1
     open(unit=1,file='out.glf23.gbflux',status='replace')
     write(1,'(32(1pe11.4,1x))') glf23_elec_pflux_out,glf23_ion_pflux_out(1:n),&
          glf23_elec_eflux_out,glf23_ion_eflux_out(1:n),&
          glf23_elec_mflux_out,glf23_ion_mflux_out(1:n),&
          glf23_elec_expwd_out,glf23_ion_expwd_out(1:n)
     close(1)

     open(unit=1,file='out.glf23.grid',status='replace')
     write(1,'(i2)') glf23_ns_in
     close(1)

  else

     print 10,'     ky:',glf23_ky_in
     ! Collect linear eigenvalues
     do i=1,2
        print 10,'(wr,wi):',glf23_eigenvalue_out(i)
     enddo

     prec = sum(abs(glf23_eigenvalue_out))

  endif

  open(unit=1,file=trim(glf23_path_in)//'out.glf23.prec')
  write(1,*) prec
  close(1)

10 format(a,10(1x,1pe11.4))
20 format(t7,a,t19,a,t31,a,t43,a,t55,a)

end program glf23
