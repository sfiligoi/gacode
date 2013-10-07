program le3

  use le3_globals
  use le3_write

  implicit none
  external :: le3_func

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
  close(1)

  iota = 1.0/q

  call le3_alloc(1)

  ! Set initial conditions
  as(:,:) = 0.0
  bs(:,:) = 0.0
  cs(:,:) = 0.0
  ds(:,:) = 0.0
  as(1,0) = rmin/rmaj

  ! Map (a,b,c,d) -> x
  call le3_map(xfunc,as,bs,cs,ds,nps,nts,'setx')

  ! Solve nonlinear system via MINPACK
  call hybrd1(le3_func,msize,xfunc,yfunc,tol,info,work,nwork)

  ! Map x -> (a,b,c,d)
  call le3_map(xfunc,as,bs,cs,ds,nps,nts,'setc')

  call le3_write_do
  call le3_alloc(0)

end program le3
