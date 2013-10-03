module le3_write
  
  implicit none

  public :: le3_write_do
  real, dimension(:,:), allocatable :: g
  real, dimension(:,:), allocatable :: bmag, bdotgrad, bdotgradB_overB, &
       vdrift_x
  real, dimension(:,:), allocatable :: mat_stream_dt, mat_stream_dp, &
       mat_trap
  real, dimension(:), allocatable :: vec_vdriftx, vec_flux, vec_upar, &
       vec_uparB
  real :: vprime
  integer, dimension(:,:), allocatable :: mat_mindx,mat_nindx,mat_mpindx,&
       mat_npindx, mat_typel, mat_typer
  integer :: matsize

contains

  subroutine le3_write_do 
    use le3_globals
    implicit none
    integer :: i,j,k,its,ips,ip
    real :: jacs
    real :: drdtbs,dzdtbs
    real :: drdpbs,dzdpbs
    real, dimension(:), allocatable :: derivvec
    real, dimension(:,:), allocatable :: gpp, gtt, gpt
    real, dimension(:,:), allocatable :: rs,zs
    real, dimension(:,:), allocatable :: bpol, btor
    real, dimension(:,:), allocatable :: dbdt, dbdp
    integer :: int_type
 
    print '(a,1pe12.5)','INFO: (le3) Root accuracy ->',sum(abs(yfunc))/size(yfunc)
    
    print 30,'a_{mn}','b_{mn}','c_{mn}','d_{mn}'
    do ips=0,nps
       do its=0,nts
          print 20,its,ips,as(its,ips),bs(its,ips),cs(its,ips),ds(its,ips)
       enddo
    enddo
    
    deallocate(t)
    deallocate(tb)
    deallocate(dtbdt)
    deallocate(dtbdp)
    deallocate(p)
    deallocate(sinm)
    deallocate(sinn)
    deallocate(cosm)
    deallocate(cosn)
    
    allocate(t(ntp))  
    allocate(tb(ntp,npp))  
    allocate(dtbdt(ntp,npp))  
    allocate(dtbdp(ntp,npp))
    allocate(p(npp))  
    allocate(sinm(ntp,0:nts))  
    allocate(cosm(ntp,0:nts))  
    allocate(sinn(npp,0:nts))  
    allocate(cosn(npp,0:nts)) 

    do i=1,ntp
       t(i) = 2*(i-1)*pi/(ntp)
    enddo
    dt = t(2)-t(1)
    do j=1,npp
       p(j) = 2*(j-1)*pi/(npp)
    enddo
    dp = p(2)-p(1)
    
    do i=1,ntp
       tb(i,:) = t(i)
    enddo
    dtbdt(:,:) = 1.0
    dtbdp(:,:) = 0.0
    
    do its=0,nts
       do i=1,ntp
          sinm(i,its) = sin(its*t(i))
          cosm(i,its) = cos(its*t(i))
       enddo
    enddo
    do ips=0,nps
       do j=1,npp
          sinn(j,ips) = sin(ips*p(j))
          cosn(j,ips) = cos(ips*p(j))
       enddo
    enddo
    
    allocate(rs(ntp,npp))
    allocate(zs(ntp,npp))
    allocate(g(ntp,npp))
    allocate(gpp(ntp,npp))
    allocate(gtt(ntp,npp))
    allocate(gpt(ntp,npp))
    
    allocate(btor(ntp,npp))
    allocate(bpol(ntp,npp))
    allocate(bmag(ntp,npp))
    allocate(dbdt(ntp,npp))
    allocate(dbdp(ntp,npp))
    allocate(bdotgrad(ntp,npp))
    allocate(bdotgradB_overB(ntp,npp))
    allocate(vdrift_x(ntp,npp))

    do j=1,npp
       do i=1,ntp
          
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
               rs(i,j),&
               zs(i,j),&
               drdtbs,&
               drdpbs,&
               dzdtbs,&
               dzdpbs,&
               jacs)
          
          g(i,j)   = 1.0/rmin * dtbdt(i,j) * jacs
          
          gpp(i,j) = rs(i,j)**2 + (drdpbs + drdtbs * dtbdp(i,j))**2 &
               + (dzdpbs + dzdtbs * dtbdp(i,j))**2
          
          gtt(i,j) = (drdtbs**2 + dzdtbs**2) * dtbdt(i,j)**2
          
          gpt(i,j) = drdtbs * dtbdt(i,j) * (drdpbs + drdtbs * dtbdp(i,j)) &
               + dzdtbs * dtbdt(i,j) * (dzdpbs + dzdtbs * dtbdp(i,j)) 

       enddo
    enddo
    
    btor(:,:) = 1.0/(rs * g)
    
    bpol(:,:) = 1.0/g * (gpt/sqrt(gtt) + iota * sqrt(gtt))
    
    bmag(:,:) = 1.0/g * sqrt(gpp + 2.0*iota*gpt + iota**2 * gtt)
    
    ! db/dtheta
    allocate(derivvec(0:ntp-1))
    do i=1,ntp-1
       derivvec(i) = -0.5*(-1)**i/tan(0.5*t(i+1))
    enddo
    dbdt(:,:) = 0.0
    do j=1,npp
       do i=1,ntp
          do ip=1,ntp
             k = ip-i
             if(k < 0) then
                k = k + ntp
             endif
             dbdt(i,j) = dbdt(i,j) + derivvec(k) * bmag(ip,j)
          enddo
       enddo
    enddo
    deallocate(derivvec)
    
    ! db/dphi
    allocate(derivvec(0:ntp-1))
    do i=1,npp-1
       derivvec(i) = -0.5*(-1)**i/tan(0.5*p(i+1))
    enddo
    dbdp(:,:) = 0.0
    do j=1,ntp
       do i=1,npp
          do ip=1,npp
             k = ip-i
             if(k < 0) then
                k = k + npp
             endif
             dbdp(j,i) = dbdp(j,i) + derivvec(k) * bmag(j,ip)
          enddo
       enddo
    enddo
    deallocate(derivvec)
    
    ! bhat dot grad = bdotgrad * (iota d/dt + d/dp)  
    bdotgrad(:,:) = 1.0/(bmag * g)
    
    ! (bhat dot grad B)/B
    bdotgradB_overB(:,:) = bdotgrad * (iota * dbdt + dbdp) / bmag
    
    ! bhat cross grad B dot grad psi / B^2
    vdrift_x(:,:) = iota/(bmag * g**2) &
         * (-dbdt * (gpp + iota * gpt) + dbdp * (gpt + iota*gtt)) / bmag**2
    
    ! construct the geo collocation matices

    matsize = 4*nts*nps+2*(nts+nps)+1
    allocate(mat_stream_dt(matsize,matsize))
    allocate(mat_stream_dp(matsize,matsize))
    allocate(mat_trap(matsize,matsize))
    allocate(mat_mindx(matsize,matsize))
    allocate(mat_nindx(matsize,matsize))
    allocate(mat_mpindx(matsize,matsize))
    allocate(mat_npindx(matsize,matsize))
    allocate(mat_typel(matsize,matsize))
    allocate(mat_typer(matsize,matsize))
    
    mat_stream_dt(:,:) = 0.0
    mat_stream_dp(:,:) = 0.0
    mat_trap(:,:)      = 0.0

    i=0
    ! amn left collocation: sin(mt) cos(np)
    int_type = 1
    do ips=0,nps
       do its=0,nts
          if(its > 0) then 
             i=i+1
             call compute_mat(int_type,its,ips,i)
          endif
       enddo
    enddo

    ! bmn left collocation: sin(mt) sin(np)
    int_type = 2
    do ips=0,nps
       do its=0,nts
          if(ips > 0 .and. its > 0) then 
             i=i+1
             call compute_mat(int_type,its,ips,i)
          endif
       enddo
    enddo
    
    ! cmn left collocation: cos(mt) cos(np)
    int_type = 3
    do ips=0,nps
       do its=0,nts
          i=i+1
          call compute_mat(int_type,its,ips,i)
       enddo
    enddo
    
    ! dmn left collocation: cos(mt) sin(np)
    int_type = 4
    do ips=0,nps
       do its=0,nts
          if(ips  > 0) then
             i=i+1
             call compute_mat(int_type,its,ips,i)
          endif
       enddo
    enddo

    mat_stream_dt(:,:) = mat_stream_dt(:,:) / (ntp*npp)
    mat_stream_dp(:,:) = mat_stream_dp(:,:) / (ntp*npp)
    mat_trap(:,:)      = mat_trap(:,:)      / (ntp*npp)

    print *, 'matsize', matsize, i,nps,nts, mat_trap(2,2)
    
    open(unit=1,file='out.le3.geomatrix',status='replace')
    do i=1,matsize
       do j=1, matsize
           write (1,'(e16.8)',advance='no') mat_trap(i,j)
           write (1,'(e16.8)',advance='no') mat_stream_dt(i,j)
           write (1,'(e16.8)',advance='no') mat_stream_dp(i,j)
           write (1,*)
       enddo
    enddo
    close(1)

    ! Construct the geo vectors

    ! flux-surface d volume / dr
    vprime = 0.0
    do i=1,ntp
       do j=1,npp
          vprime = vprime + g(i,j)
       enddo
    enddo
    vprime = vprime / (ntp*npp)
    
    allocate(vec_vdriftx(matsize))
    allocate(vec_flux(matsize))
    allocate(vec_upar(matsize))
    allocate(vec_uparB(matsize))
    vec_vdriftx(:) = 0.0
    vec_flux(:)    = 0.0
    vec_upar(:)    = 0.0
    vec_uparB(:)   = 0.0
    int_type = 0
    its=0; ips=0; i=0    ! dummy variables
    call compute_mat(int_type,its,ips,i)
    vec_vdriftx(:) = vec_vdriftx(:) / (ntp*npp)
    vec_flux(:)    = vec_flux(:)  / (ntp*npp) / vprime
    vec_uparB(:)   = vec_uparB(:) / (ntp*npp) / vprime
    vec_upar(:)    = vec_upar(:)  / (ntp*npp) 

    open(unit=1,file='out.le3.geovector',status='replace')
    do i=1,matsize
       write (1,'(e16.8)',advance='no') vec_vdriftx(i)
       write (1,'(e16.8)',advance='no') vec_flux(i)
       write (1,'(e16.8)',advance='no') vec_uparB(i)
       write (1,'(e16.8)',advance='no') vec_upar(i)
       write (1,*)
    enddo
    close(1)

    open(unit=1,file='out.le3.geoscalar',status='replace')
    write (1,'(i3)') ntp
    write (1,'(i3)') nps
    write (1,'(i3)') matsize
    write (1,'(e16.8)') vprime
    close(1)

    !do i=1,matsize
    !   do j=1, matsize
    !      if(abs(mat_trap(i,j)) > 1e-5) then
    !         print *, mat_trap(i,j), mat_mindx(i,j), mat_nindx(i,j), &
    !              mat_mpindx(i,j), mat_npindx(i,j), &
    !              mat_typel(i,j), mat_typer(i,j)
    !      endif
    !   enddo
    !enddo
    
    open(unit=1,file='out.le3.t',status='replace')
    write(1,40) t(:)
    close(1)
    
    open(unit=1,file='out.le3.p',status='replace')
    write(1,40) p(:)
    close(1)
    
    open(unit=1,file='out.le3.tb',status='replace')
    do j=1,npp
       write(1,40) tb(:,j)-t(:)
    enddo
    close(1)
    open(unit=1,file='out.le3.r',status='replace')
    do i=1,ntp
       write(1,10) rs(i,:)
    enddo
    close(1)
    open(unit=1,file='out.le3.z',status='replace')
    do i=1,ntp
       write(1,10) zs(i,:)
    enddo
    close(1)
    open(unit=1,file='out.le3.b',status='replace')
    write(1,40) bmag(:,:)
    close(1)
    
    deallocate(rs)
    deallocate(zs)
    deallocate(g)
    deallocate(gpp)
    deallocate(gtt)
    deallocate(gpt)
    deallocate(btor)
    deallocate(bpol)
    deallocate(bmag)
    deallocate(dbdt)
    deallocate(dbdp)
    deallocate(bdotgrad)
    deallocate(bdotgradB_overB)
    deallocate(vdrift_x)
    deallocate(mat_stream_dt)
    deallocate(mat_stream_dp)
    deallocate(mat_trap)
    deallocate(mat_mindx)
    deallocate(mat_nindx)
    deallocate(mat_mpindx)
    deallocate(mat_npindx)
    deallocate(mat_typel)
    deallocate(mat_typer)
    deallocate(vec_vdriftx)
    deallocate(vec_flux)
    deallocate(vec_upar)
    deallocate(vec_uparB)

10  format(200(1pe12.5,1x))
20  format('(',i2,',',i2,'):',2x,4(1pe14.7,1x))
30  format(t15,4(a,9x))
40  format(1pe12.5)

  end subroutine le3_write_do
  
  subroutine compute_mat(int_type,its,ips,i)
    use le3_globals
    implicit none
    integer, intent(in) :: int_type,its,ips,i
    real, dimension(:,:), allocatable :: mat_func
    integer :: kt,kp,jps,jts,j
    
    allocate(mat_func(ntp,npp))
    do kt=1,ntp
       do kp=1,npp
          if(int_type == 1) then
             mat_func(kt,kp) = sinm(kt,its)*cosn(kp,ips) 
          else if(int_type == 2) then
             mat_func(kt,kp) = sinm(kt,its)*sinn(kp,ips)
          else if(int_type == 3) then
             mat_func(kt,kp) = cosm(kt,its)*cosn(kp,ips)
          else if(int_type == 4) then
             mat_func(kt,kp) = cosm(kt,its)*sinn(kp,ips)
          else
             mat_func(kt,kp) = 1.0
          endif
       enddo
    enddo
    
    j=0
    do jps=0,nps
       do jts=0,nts
          
          if(jts > 0) then
             ! amn right collocation
             j=j+1
             do kt=1,ntp
                do kp=1,npp
                   if(int_type == 0) then
                      vec_vdriftx(j) = vec_vdriftx(j) &
                           + mat_func(kt,kp) * sinm(kt,jts)*cosn(kp,jps) &
                           * vdrift_x(kt,kp)
                      vec_flux(j)    = vec_flux(j) &
                           + mat_func(kt,kp) * sinm(kt,jts)*cosn(kp,jps) &
                           * vdrift_x(kt,kp) * g(kt,kp)
                      vec_upar(j)    = vec_upar(j) &
                           + mat_func(kt,kp) * sinm(kt,jts)*cosn(kp,jps) 
                      vec_uparB(j)   = vec_uparB(j) &
                           + mat_func(kt,kp) * sinm(kt,jts)*cosn(kp,jps) &
                           * bmag(kt,kp) * g(kt,kp)
                   else
                      mat_mindx(i,j)=its
                      mat_nindx(i,j)=ips
                      mat_mpindx(i,j)=jts
                      mat_npindx(i,j)=jps
                      mat_typel=int_type
                      mat_typer=1
                      mat_trap(i,j) = mat_trap(i,j) &
                           + mat_func(kt,kp) * sinm(kt,jts)*cosn(kp,jps) &
                           * bdotgradB_overB(kt,kp) 
                      mat_stream_dt(i,j) = mat_stream_dt(i,j) &
                           + mat_func(kt,kp) * jts*cosm(kt,jts)*cosn(kp,jps) &
                           * bdotgrad(kt,kp) * iota
                      mat_stream_dp(i,j) = mat_stream_dp(i,j) &
                           + mat_func(kt,kp)*(-jps)*sinm(kt,jts)*sinn(kp,jps) &
                           * bdotgrad(kt,kp) 
                   endif
                enddo
             enddo
          endif
          
          if(jps > 0 .and. jts > 0) then
             ! bmn right collocation
             j=j+1
             do kt=1,ntp
                do kp=1,npp
                   if(int_type == 0) then
                      vec_vdriftx(j) = vec_vdriftx(j) &
                           + mat_func(kt,kp) * sinm(kt,jts)*sinn(kp,jps) &
                           * vdrift_x(kt,kp)
                      vec_flux(j)    = vec_flux(j) &
                           + mat_func(kt,kp) * sinm(kt,jts)*sinn(kp,jps) &
                           * vdrift_x(kt,kp) * g(kt,kp)
                      vec_upar(j)    = vec_upar(j) &
                           + mat_func(kt,kp) * sinm(kt,jts)*sinn(kp,jps) 
                      vec_uparB(j)   = vec_uparB(j) &
                           + mat_func(kt,kp) * sinm(kt,jts)*sinn(kp,jps) &
                           * bmag(kt,kp) * g(kt,kp)
                   else
                      mat_mindx(i,j)=its
                      mat_nindx(i,j)=ips
                      mat_mpindx(i,j)=jts
                      mat_npindx(i,j)=jps
                      mat_typel=int_type
                      mat_typer=2
                      mat_trap(i,j) = mat_trap(i,j) &
                           + mat_func(kt,kp) * sinm(kt,jts)*sinn(kp,jps) &
                           * bdotgradB_overB(kt,kp)  
                      mat_stream_dt(i,j) = mat_stream_dt(i,j) &
                           + mat_func(kt,kp) * jts*cosm(kt,jts)*sinn(kp,jps) &
                           * bdotgrad(kt,kp) * iota
                      mat_stream_dp(i,j) = mat_stream_dp(i,j) &
                           + mat_func(kt,kp) * jps*sinm(kt,jts)*cosn(kp,jps) &
                           * bdotgrad(kt,kp) 
                   endif
                enddo
             enddo
          endif
          
          ! cmn right collocation
          j=j+1
          do kt=1,ntp
             do kp=1,npp
                if(int_type == 0) then
                   vec_vdriftx(j) = vec_vdriftx(j) &
                        + mat_func(kt,kp) * cosm(kt,jts)*cosn(kp,jps) &
                        * vdrift_x(kt,kp)
                   vec_flux(j)    = vec_flux(j) &
                        + mat_func(kt,kp) * cosm(kt,jts)*cosn(kp,jps) &
                        * vdrift_x(kt,kp) * g(kt,kp)
                   vec_upar(j)    = vec_upar(j) &
                        + mat_func(kt,kp) * cosm(kt,jts)*cosn(kp,jps) 
                   vec_uparB(j)   = vec_uparB(j) &
                        + mat_func(kt,kp) * cosm(kt,jts)*cosn(kp,jps) &
                        * bmag(kt,kp) * g(kt,kp)
                else
                   mat_mindx(i,j)=its
                   mat_nindx(i,j)=ips
                   mat_mpindx(i,j)=jts
                   mat_npindx(i,j)=jps
                   mat_typel=int_type
                   mat_typer=3
                   mat_trap(i,j) = mat_trap(i,j) &
                        + mat_func(kt,kp) * cosm(kt,jts)*cosn(kp,jps) &
                        * bdotgradB_overB(kt,kp) 
                   mat_stream_dt(i,j) = mat_stream_dt(i,j) &
                        + mat_func(kt,kp) * (-jts)*sinm(kt,jts)*cosn(kp,jps) &
                        * bdotgrad(kt,kp) * iota
                   mat_stream_dp(i,j) = mat_stream_dp(i,j) &
                        + mat_func(kt,kp) * (-jps)*cosm(kt,jts)*sinn(kp,jps) &
                        * bdotgrad(kt,kp) 
                endif
             enddo
          enddo
          
          if(jps > 0) then
             ! dmn right collocation
             j=j+1
             do kt=1,ntp
                do kp=1,npp
                   if(int_type == 0) then
                      vec_vdriftx(j) = vec_vdriftx(j) &
                           + mat_func(kt,kp) * cosm(kt,jts)*sinn(kp,jps) &
                           * vdrift_x(kt,kp)
                      vec_flux(j)    = vec_flux(j) &
                           + mat_func(kt,kp) * cosm(kt,jts)*sinn(kp,jps) &
                           * vdrift_x(kt,kp) * g(kt,kp)
                      vec_upar(j)    = vec_upar(j) &
                           + mat_func(kt,kp) * cosm(kt,jts)*sinn(kp,jps) 
                      vec_uparB(j)   = vec_uparB(j) &
                           + mat_func(kt,kp) * cosm(kt,jts)*sinn(kp,jps) &
                           * bmag(kt,kp) * g(kt,kp)
                   else
                      mat_mindx(i,j)=its
                      mat_nindx(i,j)=ips
                      mat_mpindx(i,j)=jts
                      mat_npindx(i,j)=jps
                      mat_typel=int_type
                      mat_typer=4
                      mat_trap(i,j) = mat_trap(i,j) &
                           + mat_func(kt,kp) * cosm(kt,jts)*sinn(kp,jps) & 
                           * bdotgradB_overB(kt,kp)  
                      mat_stream_dt(i,j) = mat_stream_dt(i,j) &
                           + mat_func(kt,kp)*(-jts)*sinm(kt,jts)*sinn(kp,jps) &
                           * bdotgrad(kt,kp) * iota
                      mat_stream_dp(i,j) = mat_stream_dp(i,j) &
                           + mat_func(kt,kp) * jps*cosm(kt,jts)*cosn(kp,jps) &
                           * bdotgrad(kt,kp)
                   endif
                enddo
             enddo
          endif
          
       enddo
    enddo

    deallocate(mat_func)
    
  end subroutine compute_mat

end module le3_write
