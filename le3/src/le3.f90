program le3

  use le3_globals

  implicit none
  integer :: it, ip, im
  external :: le3_func

  open(unit=1,file='input.le3.gen',status='old')
  read(1,*) nts
  read(1,*) nps
  read(1,*) equilibrium_model
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
  read(1,*) beta_star
  read(1,*) hmin
  read(1,*) dhmindr
  read(1,*) q
  read(1,*) s
  read(1,*) m
  read(1,*) n
  read(1,*) tol
  close(1)

  if(equilibrium_model == 1) then
     open(unit=1,file='out.profiles_3d.geo',status='old')
     read(1,*) nts_geo
     read(1,*) nps_geo
     allocate(r_geo(0:nts_geo,0:nps_geo,4))
     allocate(z_geo(0:nts_geo,0:nps_geo,4))
     allocate(rd_geo(0:nts_geo,0:nps_geo,4))
     allocate(zd_geo(0:nts_geo,0:nps_geo,4))
     do it=0,nts_geo
        do ip=0,nps_geo
           do im=1,4
              read(1,*) r_geo(it,ip,im), z_geo(it,ip,im), &
                   rd_geo(it,ip,im), zd_geo(it,ip,im)
           enddo
        enddo
     enddo
     close(1)
  endif

  nt = 2*nts+2
  if (nps > 0) then
     np = 2*nps+2
  else
     np = 1
  endif

  ! i = 1/q
  iota   = 1.0/q

  ! i' = B_unit*di/d(chi)
  iota_p = -s/(q*rmin**2)

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

  call le3_geometry
  call le3_geometry_rho

  call le3_alloc(0)
  if(equilibrium_model == 1) then
     deallocate(r_geo)
     deallocate(z_geo)
     deallocate(rd_geo)
     deallocate(zd_geo)
  endif

end program le3
