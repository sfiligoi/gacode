program tglf

  use tglf_interface

  implicit none

  integer :: ierr

  include 'mpif.h'

  call MPI_INIT(ierr)

  tglf_dump_flag_in=.true.

  call tglf_read_input()
  call tglf_run() 

  print 10,'  Elec. particle flux :',tglf_elec_pflux_out
  print 10,'    Elec. energy flux :',tglf_elec_eflux_out
  print 10,'Elec. energy flux (L) :',tglf_elec_eflux_low_out
  print *
  print 10,'    Ion particle flux :',tglf_ion_pflux_out(1:tglf_ns_in-1)
  print 10,'      Ion energy flux :',tglf_ion_eflux_out(1:tglf_ns_in-1)
  print 10,'  Ion energy flux (L) :',tglf_ion_eflux_low_out(1:tglf_ns_in-1)

  call MPI_FINALIZE(ierr)

10 format(a,5(1x,1pe13.6))

end program tglf
