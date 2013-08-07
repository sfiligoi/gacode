program le3

  use le3_globals
  use le3_driver

  implicit none
  integer :: it,ip

  open(unit=1,file='input.le3.gen',status='old')
  read(1,*) nt
  read(1,*) np
  read(1,*) nts
  read(1,*) nps
  read(1,*) rmin
  read(1,*) rmaj
  read(1,*) hmin
  read(1,*) q
  read(1,*) m
  read(1,*) n
  read(1,*) tol
  read(1,*) restart_flag
  read(1,*) solve_method
  close(1)

  call le3_alloc(1)

  if (restart_flag == 1) then
     open(unit=1,file='out.le3.tb',status='old')
     do it=1,nt
        read(1,*) tb(it,:)
     enddo
     close(1)
  else   
     do ip=1,np
        tb(:,ip) = t(:)
     enddo

     as = 0.0
     bs = 0.0

  endif

  call le3_solver

  call le3_alloc(0)

end program le3
