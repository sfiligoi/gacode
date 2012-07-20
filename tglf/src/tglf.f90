!---------------------------------------------------------
! tglf_run.f90
!
! PURPOSE:
!  Manage standalone TGLF run by calling
!
!    call tglf_read_input()
!    call tglf_run() 
!
!  and then adding extra I/O. 
!---------------------------------------------------------

program tglf

  use tglf_interface

  implicit none

  integer :: i
  character (len=4) :: tag(5)=(/'ion1','ion2','ion3','ion4','ion5'/)
  real :: prec

  call tglf_read_input()
  call tglf_run() 

  if (tglf_use_transport_model_in)then

     print 20,'Gam/Gam_GB','    Q/Q_GB','Q_low/Q_GB','  Pi/Pi_GB', '    S/S_GB'
     print 10,'elec',&
          tglf_elec_pflux_out,&
          tglf_elec_eflux_out,&
          tglf_elec_eflux_low_out,&
          tglf_elec_mflux_out,&
          tglf_elec_expwd_out

     prec = abs(tglf_elec_pflux_out)+&
          abs(tglf_elec_eflux_out)+&
          abs(tglf_elec_eflux_low_out)+&
          abs(tglf_elec_mflux_out)

     do i=1,tglf_ns_in-1
        print 10,tag(i),&
             tglf_ion_pflux_out(i),&
             tglf_ion_eflux_out(i),&
             tglf_ion_eflux_low_out(i),&
             tglf_ion_mflux_out(i),&
             tglf_ion_expwd_out(i)

        prec = prec+&
             abs(tglf_ion_pflux_out(i))+&
             abs(tglf_ion_eflux_out(i))+&
             abs(tglf_ion_eflux_low_out(i))+&
             abs(tglf_ion_mflux_out(i))
     enddo

  else

     print 10,'     ky:',tglf_ky_in
     ! Collect linear eigenvalues
     do i=1,tglf_nmodes_in
        print 10,'(wr,wi):',tglf_eigenvalue_out(i)
     enddo

     prec = sum(abs(tglf_eigenvalue_out))

  endif

  open(unit=1,file=trim(tglf_path_in)//'out.tglf.prec')
  write(1,*) prec
  close(1)

10 format(a,6(1x,1pe11.4))
20 format(t7,a,t19,a,t31,a,t43,a,t55,a)

end program tglf
