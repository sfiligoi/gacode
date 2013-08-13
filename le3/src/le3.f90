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
  read(1,*) shift
  read(1,*) kappa
  read(1,*) s_kappa
  read(1,*) delta
  read(1,*) s_delta
  read(1,*) zeta
  read(1,*) s_zeta
  read(1,*) zmag
  read(1,*) dzmag
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

     if (solve_method == 1) then 
        do ip=1,np
           tb(:,ip) = t(:)
        enddo
     endif

  endif

  call le3_solver

  call le3_alloc(0)

end program le3
