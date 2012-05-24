program tglf
  use tglf_mpi
  use tglf_interface

  implicit none

  integer :: i,itime0,itime1,irate
  character (len=4) :: tag(5)=(/'ion1','ion2','ion3','ion4','ion5'/)
  real(8) :: prec
  integer :: ierr

  !\---------------------------
  ! executive code starts here
  !/

  ! initialize MPI
  call MPI_INIT(ierr)

  iCommTglf = MPI_COMM_WORLD
  call MPI_COMM_RANK(iCommTglf,iProcTglf,ierr)
  call MPI_COMM_SIZE(iCommTglf,nProcTglf,ierr)

  ! read inputs    
  call tglf_read_input()

  call system_clock(itime0)
  call tglf_run_mpi() 
  call system_clock(itime1,irate)

  if(iProcTglf == iProc0Tglf) then
    print 20,'Gam/Gam_GB','    Q/Q_GB','Q_low/Q_GB','  Pi/Pi_GB'
    print 10,'elec',&
          tglf_elec_pflux_out,&
          tglf_elec_eflux_out,&
          tglf_elec_eflux_low_out,&
          tglf_elec_mflux_out

    prec = abs(tglf_elec_pflux_out)+&
           abs(tglf_elec_eflux_out)+&
           abs(tglf_elec_eflux_low_out)+&
           abs(tglf_elec_mflux_out)

    do i=1,tglf_ns_in-1
       print 10,tag(i),&
            tglf_ion_pflux_out(i),&
            tglf_ion_eflux_out(i),&
            tglf_ion_eflux_low_out(i),&
            tglf_ion_mflux_out(i)

       prec = prec+&
            abs(tglf_ion_pflux_out(i))+&
            abs(tglf_ion_eflux_out(i))+&
            abs(tglf_ion_eflux_low_out(i))+&
            abs(tglf_ion_mflux_out(i))
    enddo

    write(*,*) 'system time used:',real(itime1-itime0)/real(irate),' seconds'

    open(unit=1,file=trim(tglf_path_in)//'out.tglf.precision')
    write(1,*) prec
    close(1)
   end if

20 format(t15,a,t25,a,t43,a,t47,a)
10 format(a,5(1x,1pe13.6))

end program tglf
