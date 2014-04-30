subroutine le3_geometry

  use le3_globals

  implicit none
  integer :: i,j,k,its,ips,ip,kt,kp
  real, dimension(:), allocatable :: vec_vdriftx
  real, dimension(:), allocatable :: vec_flux
  real, dimension(:), allocatable :: vec_upar
  real, dimension(:), allocatable :: vec_uparB 
  real, dimension(:), allocatable :: vec_fsa
  real, dimension(:), allocatable :: vec_bmag
  real, dimension(:), allocatable :: vec_thetabar
  real, dimension(:), allocatable :: vec_ntv

  print '(a,1pe12.5)','INFO: (le3) Root accuracy ->',sum(abs(yfunc))/size(yfunc)

  if (nps > 0) then
     print 30,'a_{mn}','b_{mn}','c_{mn}','d_{mn}'
     do ips=0,nps
        do its=0,nts
           print 20,its,ips,as(its,ips),bs(its,ips),cs(its,ips),ds(its,ips)
        enddo
     enddo
  else
     print 30,'a_{m}','c_{m}'
     ips=0
     do its=0,nts
        print 10,its,as(its,ips),cs(its,ips)
     enddo
  endif

  allocate(g(nt,np))
  allocate(gpp(nt,np))
  allocate(gtt(nt,np))
  allocate(gpt(nt,np))
  allocate(cosu(nt,np))
  allocate(chi1(nt,np))
  allocate(chi1p(nt,np))
  allocate(chi1t(nt,np))
  allocate(drdpt(nt,np))
  allocate(dzdpt(nt,np))
  allocate(drdt(nt,np))
  allocate(drdp(nt,np))
  allocate(dzdt(nt,np))
  allocate(dzdp(nt,np))
  allocate(btor(nt,np))
  allocate(bpol(nt,np))
  allocate(bmag(nt,np))
  allocate(dbdt(nt,np))
  allocate(dbdp(nt,np))
  allocate(bdotgrad(nt,np))
  allocate(bdotgradB_overB(nt,np))
  allocate(vdrift_x(nt,np))
  allocate(vexb_dt(nt,np))
  allocate(vexb_dp(nt,np))
  allocate(dgdp(nt,np))

  do i=1,nt
     tb(i,:) = t(i)
  enddo
  dtbdt(:,:) = 1.0
  dtbdp(:,:) = 0.0

  do j=1,np
     do i=1,nt

        do ips=0,nps
           do its=0,nts

              tb(i,j) = tb(i,j) + &
                   sinm(i,its)*(bs(its,ips)*sinn(j,ips)+as(its,ips)*cosn(j,ips)) &
                   +cosm(i,its)*(ds(its,ips)*sinn(j,ips)+cs(its,ips)*cosn(j,ips))

              dtbdt(i,j) = dtbdt(i,j) + &
                   its*cosm(i,its)*(bs(its,ips)*sinn(j,ips)+as(its,ips)*cosn(j,ips)) &
                   -its*sinm(i,its)*(ds(its,ips)*sinn(j,ips)+cs(its,ips)*cosn(j,ips))

              dtbdp(i,j) = dtbdp(i,j) + &
                   ips*sinm(i,its)*(bs(its,ips)*cosn(j,ips)-as(its,ips)*sinn(j,ips)) &
                   +ips*cosm(i,its)*(ds(its,ips)*cosn(j,ips)-cs(its,ips)*sinn(j,ips))

           enddo
        enddo

        call le3_rz(tb(i,j),&
             p(j),&
             r(i,j),&
             z(i,j),&
             drdtb(i,j),&
             drdpb(i,j),&
             dzdtb(i,j),&
             dzdpb(i,j),&
             jac(i,j),&
             rc(i,j))

     enddo
  enddo

  ! sqrt(g)
  g(:,:) = dtbdt*jac/rmin

  ! R_t
  drdt(:,:) = drdtb*dtbdt

  ! Z_t
  dzdt(:,:) = dzdtb*dtbdt

  ! R_p
  drdp(:,:) = drdpb+drdtb*dtbdp

  ! Z_p
  dzdp(:,:) = dzdpb+dzdtb*dtbdp

  ! g_pp = R^2 + R_p^2 + Z_p^2
  gpp(:,:) = r**2+drdp**2+dzdp**2

  ! g_tt = R_t^2 + Z_t^2
  gtt(:,:) = drdt**2+dzdt**2

  ! g_pt = R_t R_p + Z_t Z_p
  gpt(:,:) = drdt*drdp+dzdt*dzdp

  ! cos(u) = Z_t/sqrt(g_tt)
  cosu(:,:) = dzdt/sqrt(gtt)

  btor(:,:) = 1.0/(r*g)

  bpol(:,:) = 1.0/g*(gpt/sqrt(gtt)+iota*sqrt(gtt))

  bmag(:,:) = 1.0/g*sqrt(gpp+2.0*iota*gpt+iota**2*gtt)

  ! chi_1 = x R /sqrt(g)
  chi1(:,:) = sqrt(gtt)*r/g

  ! db/dtheta
  allocate(deriv_t(0:nt-1))
  deriv_t(0) = 0.0
  do i=1,nt-1
     deriv_t(i) = -0.5*(-1)**i/tan(0.5*t(i+1))
  enddo
  dbdt(:,:) = 0.0
  chi1t(:,:) = 0.0
  do j=1,np
     do i=1,nt
        do ip=1,nt
           k = ip-i
           if (k < 0) then
              k = k + nt
           endif
           dbdt(i,j)  = dbdt(i,j)+deriv_t(k)*bmag(ip,j)
           chi1t(i,j) = chi1t(i,j)+deriv_t(k)*chi1(ip,j)
        enddo
     enddo
  enddo

  ! db/dphi and dg/dphi
  allocate(deriv_p(0:np-1))
  deriv_p(0) = 0.0
  do i=1,np-1
     deriv_p(i) = -0.5*(-1)**i/tan(0.5*p(i+1))
  enddo
  dbdp(:,:) = 0.0
  dgdp(:,:) = 0.0
  chi1p(:,:) = 0.0
  drdpt(:,:) = 0.0
  dzdpt(:,:) = 0.0
  do j=1,nt
     do i=1,np
        do ip=1,np
           k = ip-i
           if (k < 0) then
              k = k + np
           endif
           dbdp(j,i)  = dbdp(j,i) +deriv_p(k)*bmag(j,ip)
           dgdp(j,i)  = dgdp(j,i) +deriv_p(k)*g(j,ip)
           chi1p(j,i) = chi1p(j,i)+deriv_p(k)*chi1(j,ip)
           drdpt(j,i) = drdpt(j,i)+deriv_p(k)*drdt(j,ip)
           dzdpt(j,i) = dzdpt(j,i)+deriv_p(k)*dzdt(j,ip)
        enddo
     enddo
  enddo

  ! bhat dot grad = bdotgrad * (iota d/dt + d/dp)  
  bdotgrad(:,:) = 1.0/(bmag * g)

  ! (bhat dot grad B)/B
  bdotgradB_overB(:,:) = bdotgrad * (iota * dbdt + dbdp) / bmag

  ! bhat cross grad B dot grad r / B^2
  vdrift_x(:,:) = 1/(rmin*bmag * g**2) &
       * (-dbdt * (gpp + iota * gpt) + dbdp * (gpt + iota*gtt)) / bmag**2

  ! -bhat cross grad f dot grad r / B 
  vexb_dt(:,:) = -1/(rmin*bmag * g**2) &
       * (-(gpp + iota * gpt)) 
  vexb_dp(:,:) = -1/(rmin*bmag * g**2) & 
       * (gpt + iota*gtt) / bmag

  ! construct the geo collocation matices

  matsize = 4*nts*nps+2*(nts+nps)+1

  allocate(m_indx(matsize))
  allocate(n_indx(matsize))
  allocate(itype(matsize))

  allocate(bk(nt,np))
  allocate(bkp(nt,np))
  allocate(bk_t(nt,np))
  allocate(bkp_t(nt,np))
  allocate(bk_p(nt,np))
  allocate(bkp_p(nt,np))

  allocate(mat_stream_dt(matsize,matsize))
  allocate(mat_stream_dp(matsize,matsize))
  allocate(mat_trap(matsize,matsize))
  allocate(mat_coll(matsize,matsize))
  allocate(mat_vexb_dt(matsize,matsize))
  allocate(mat_vexb_dp(matsize,matsize))
  allocate(vec_vdriftx(matsize))
  allocate(vec_flux(matsize))
  allocate(vec_upar(matsize))
  allocate(vec_uparB(matsize))
  allocate(vec_fsa(matsize))
  allocate(vec_bmag(matsize))
  allocate(vec_thetabar(matsize))
  allocate(vec_ntv(matsize))

  i=0
  do ips=0,nps
     do its=0,nts
        if(its > 0) then 
           ! amn
           i=i+1
           itype(i) = 1
           m_indx(i)     = its
           n_indx(i)     = ips
        endif
        if(ips > 0 .and. its > 0) then 
           !bmn
           i=i+1
           itype(i) = 2
           m_indx(i)     = its
           n_indx(i)     = ips
        endif
        ! cmn
        i=i+1
        itype(i) = 3
        m_indx(i)     = its
        n_indx(i)     = ips
        if(its == 0 .and. ips == 0) then
           indx_c00 = i
        endif
        if(ips  > 0) then
           ! dmn
           i=i+1
           itype(i) = 4
           m_indx(i)     = its
           n_indx(i)     = ips
        endif
     enddo
  enddo

  mat_stream_dt(:,:) = 0.0
  mat_stream_dp(:,:) = 0.0
  mat_trap(:,:)      = 0.0
  mat_coll(:,:)      = 0.0
  mat_vexb_dt(:,:)   = 0.0
  mat_vexb_dp(:,:)   = 0.0
  vec_vdriftx(:)  = 0.0
  vec_flux(:)     = 0.0
  vec_upar(:)     = 0.0
  vec_uparB(:)    = 0.0
  vec_fsa(:)      = 0.0
  vec_bmag(:)     = 0.0
  vec_thetabar(:) = 0.0
  vec_ntv(:)      = 0.0

  do i=1,matsize
     call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'d0')
     do j=1,matsize
        call le3_basis(itype(j),m_indx(j),n_indx(j),bkp(:,:),'d0')
        call le3_basis(itype(j),m_indx(j),n_indx(j),bkp_t(:,:),'dt')
        call le3_basis(itype(j),m_indx(j),n_indx(j),bkp_p(:,:),'dp')
        do kt=1,nt
           do kp=1,np
              mat_trap(i,j) = mat_trap(i,j) &
                   + bk(kt,kp) * bkp(kt,kp) * bdotgradB_overB(kt,kp) 
              mat_stream_dt(i,j) = mat_stream_dt(i,j) &
                   + bk(kt,kp) * bkp_t(kt,kp) * bdotgrad(kt,kp) * iota
              mat_stream_dp(i,j) = mat_stream_dp(i,j) &
                   + bk(kt,kp) * bkp_p(kt,kp) * bdotgrad(kt,kp) 
              mat_coll(i,j) = mat_coll(i,j) &
                   + bk(kt,kp) * bkp(kt,kp)
              mat_vexb_dt(i,j) = mat_vexb_dt(i,j) &
                   + bk(kt,kp) * bkp_t(kt,kp) * vexb_dt(kt,kp)
              mat_vexb_dp(i,j) = mat_vexb_dp(i,j) &
                   + bk(kt,kp) * bkp_p(kt,kp) * vexb_dp(kt,kp)
           enddo
        enddo
     enddo
     do kt=1,nt
        do kp=1,np
           vec_thetabar(i) = vec_thetabar(i) &
                + bk(kt,kp) * (tb(kt,kp)-t(kt))
           vec_vdriftx(i) = vec_vdriftx(i) &
                + bk(kt,kp) * vdrift_x(kt,kp)
           vec_flux(i)    = vec_flux(i) &
                + bk(kt,kp) * vdrift_x(kt,kp) * g(kt,kp)
           vec_upar(i)    = vec_upar(i) &
                + bk(kt,kp) 
           vec_uparB(i)   = vec_uparB(i) &
                + bk(kt,kp) * bmag(kt,kp) * g(kt,kp)
           vec_fsa(i)   = vec_fsa(i) &
                + bk(kt,kp) * g(kt,kp)
           vec_bmag(i)   = vec_bmag(i) &
                + bk(kt,kp) * bmag(kt,kp) 
           vec_ntv(i)   = vec_ntv(i) &
                + bk(kt,kp) * dgdp(kt,kp)
        enddo
     enddo
  enddo

  mat_stream_dt(:,:) = mat_stream_dt(:,:) / (nt*np)
  mat_stream_dp(:,:) = mat_stream_dp(:,:) / (nt*np)
  mat_trap(:,:)      = mat_trap(:,:)      / (nt*np)
  mat_coll(:,:)      = mat_coll(:,:)      / (nt*np)
  mat_vexb_dt(:,:)   = mat_vexb_dt(:,:)   / (nt*np)
  mat_vexb_dp(:,:)   = mat_vexb_dp(:,:)   / (nt*np)

  open(unit=1,file='out.le3.geomatrix',status='replace')
  do i=1,matsize
     do j=1, matsize
        write (1,'(e16.8)',advance='no') mat_trap(i,j)
        write (1,'(e16.8)',advance='no') mat_stream_dt(i,j)
        write (1,'(e16.8)',advance='no') mat_stream_dp(i,j)
        write (1,'(e16.8)',advance='no') mat_coll(i,j)
        write (1,'(e16.8)',advance='no') mat_vexb_dt(i,j)
        write (1,'(e16.8)',advance='no') mat_vexb_dp(i,j)
        write (1,*)
     enddo
  enddo
  close(1)

  ! Construct the geo vectors

  ! flux-surface d volume / dr
  vprime = 0.0
  do i=1,nt
     do j=1,np
        vprime = vprime + g(i,j)
     enddo
  enddo
  vprime = vprime / (nt*np)

  vec_thetabar(:) = vec_thetabar(:) / (nt*np)
  vec_vdriftx(:)  = vec_vdriftx(:) / (nt*np)
  vec_flux(:)     = vec_flux(:)  / (nt*np) / vprime
  vec_uparB(:)    = vec_uparB(:) / (nt*np) / vprime
  vec_upar(:)     = vec_upar(:)  / (nt*np) * (2*pi)**2
  vec_fsa(:)      = vec_fsa (:)  / (nt*np) / vprime
  vec_bmag(:)     = vec_bmag (:)  / (nt*np)
  vec_ntv(:)      = vec_ntv(:) / (nt*np) / vprime

  open(unit=1,file='out.le3.geovector',status='replace')
  do i=1,matsize
     write (1,'(e16.8)',advance='no') vec_thetabar(i)
     write (1,'(e16.8)',advance='no') vec_vdriftx(i)
     write (1,'(e16.8)',advance='no') vec_flux(i)
     write (1,'(e16.8)',advance='no') vec_uparB(i)
     write (1,'(e16.8)',advance='no') vec_upar(i)
     write (1,'(e16.8)',advance='no') vec_fsa(i)
     write (1,'(e16.8)',advance='no') vec_bmag(i)
     write (1,'(e16.8)',advance='no') vec_ntv(i)
     write (1,*)
  enddo
  close(1)

  !ips=1
  !do its=1,nt
  !   print *, t(its),vdrift_x(its,ips)*0.001
  !enddo

  open(unit=1,file='out.le3.geoscalar',status='replace')
  write (1,'(i3)') nts
  write (1,'(i3)') nps
  write (1,'(i3)') matsize
  write (1,'(i3)') indx_c00
  close(1)

  call le3_compute_theory

  deallocate(bdotgrad)
  deallocate(bdotgradB_overB)
  deallocate(vdrift_x)
  deallocate(vexb_dt)
  deallocate(vexb_dp)
  deallocate(mat_stream_dt)
  deallocate(mat_stream_dp)
  deallocate(mat_trap)
  deallocate(mat_coll)
  deallocate(mat_vexb_dt)
  deallocate(mat_vexb_dp)
  deallocate(vec_vdriftx)
  deallocate(vec_flux)
  deallocate(vec_upar)
  deallocate(vec_uparB)
  deallocate(vec_fsa)
  deallocate(vec_bmag)
  deallocate(vec_thetabar)
  deallocate(vec_ntv)

10 format('(',i2,'):',2x,2(1pe14.7,1x))
20 format('(',i2,',',i2,'):',2x,4(1pe14.7,1x))
30 format(t15,4(a,9x))

end subroutine le3_geometry

