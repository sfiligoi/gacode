program peek

  use gyro_globals
  use peek_globals

  implicit none
  
  character (len=100) :: simdir='./'

  open(unit=1,file='/tmp/peekcfg',status='old')
  read(1,*) window
  close(1)

  ! Read GYRO data
  call read_profile_vugyro(trim(simdir))

  print '(t2,a,a)','GYRO Quick Summary: ',simdir

  call read_time(trim(simdir))

end program peek
