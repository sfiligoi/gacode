program tglf

  use tglf_mpi
  use tglf_interface

  implicit none

  integer :: i
  integer :: n
  character (len=4) :: tag(5)=(/'ion1','ion2','ion3','ion4','ion5'/)
  real(8) :: prec
  integer :: ierr

  ! initialize MPI
  call MPI_INIT(ierr)

  iCommTglf = MPI_COMM_WORLD
  call MPI_COMM_RANK(iCommTglf,iProcTglf,ierr)
  call MPI_COMM_SIZE(iCommTglf,nProcTglf,ierr)

  call tglf_read_input()
  call tglf_run_mpi() 

  if (iProcTglf == iProc0Tglf) then

     print 20,'Gam/Gam_GB',' Q/Q_GB','Q_low/Q_GB',' Pi/Pi_GB', ' S/S_GB'
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

     ! Output to file

     n = tglf_ns_in-1
     open(unit=1,file='out.tglf.gbflux',status='replace')
     write(1,'(32(1pe11.4,1x))') tglf_elec_pflux_out,tglf_ion_pflux_out(1:n),&
          tglf_elec_eflux_out,tglf_ion_eflux_out(1:n),&
          tglf_elec_mflux_out,tglf_ion_mflux_out(1:n),&
          tglf_elec_expwd_out,tglf_ion_expwd_out(1:n)
     close(1)

     open(unit=1,file='out.tglf.grid',status='replace')
     write(1,'(i2)') tglf_ns_in,tglf_nxgrid_in
     close(1)

     open(unit=1,file=trim(tglf_path_in)//'out.tglf.prec')
     write(1,*) prec
     close(1)

  endif

10 format(a,10(1x,1pe11.4))
20 format(t7,a,t19,a,t31,a,t43,a,t55,a)


end program
