!ubroutine nbfinish
!
  use nbi_com
!
  !  finish up the fast ion calculations
  !  this includes bbneut; outptb2 also handles copying of 2d outputs
  !  into xplasma.
  print *,'in nbfinish,pbe =',pbe
  call outptb2
  print *,' 2 nbfinish,pbe =',pbe
  !  compute sawtooth mixing if needed
!
  call pre_sawnbi
  print *,' 3 nbfinish,pbe =',pbe
  !  load 1d profiles into xplasma (generated code).
!
  call nbo_load
!
  !  close lost orbit output files (if any)...
!
  if(NLBOUT.and.NLMCANY) then
     call nbxfin
  endif
  print *,'end nbfinish,pbe =',pbe
!nd subroutine nbfinish
