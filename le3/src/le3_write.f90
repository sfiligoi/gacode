module le3_write
  
  implicit none

  public :: le3_write_do
  real, dimension(:,:), allocatable :: g
  real, dimension(:,:), allocatable :: bmag, bdotgrad, bdotgradB_overB, &
       vdrift_x
  real, dimension(:,:), allocatable :: mat_stream_dt, mat_stream_dp, &
       mat_trap
  real, dimension(:), allocatable :: vec_vdriftx, vec_flux, vec_upar, &
       vec_uparB, vec_fsa, vec_bmag, vec_thetabar
  real :: vprime
  integer :: matsize
  integer :: indx_c00

contains

  subroutine le3_write_do 
    use le3_globals
    implicit none
    integer :: i,j,k,its,ips,ip, kt, kp
    real :: jacs
    real :: drdtbs,dzdtbs
    real :: drdpbs,dzdpbs
    real, dimension(:), allocatable :: derivvec
    real, dimension(:,:), allocatable :: gpp, gtt, gpt
    real, dimension(:,:), allocatable :: rs,zs
    real, dimension(:,:), allocatable :: bpol, btor
    real, dimension(:,:), allocatable :: dbdt, dbdp
    real, dimension(:,:), allocatable :: basis, basis_prime, &
         basis_dt_prime, basis_dp_prime
    integer, dimension(:), allocatable :: m_indx, n_indx, itype
 
    print '(a,1pe12.5)','INFO: (le3) Root accuracy ->',sum(abs(yfunc))/size(yfunc)
    
    print 30,'a_{mn}','b_{mn}','c_{mn}','d_{mn}'
    do ips=0,nps
       do its=0,nts
          print 20,its,ips,as(its,ips),bs(its,ips),cs(its,ips),ds(its,ips)
       enddo
    enddo
    
    allocate(rs(nt,np))
    allocate(zs(nt,np))
    allocate(g(nt,np))
    allocate(gpp(nt,np))
    allocate(gtt(nt,np))
    allocate(gpt(nt,np))
    
    allocate(btor(nt,np))
    allocate(bpol(nt,np))
    allocate(bmag(nt,np))
    allocate(dbdt(nt,np))
    allocate(dbdp(nt,np))
    allocate(bdotgrad(nt,np))
    allocate(bdotgradB_overB(nt,np))
    allocate(vdrift_x(nt,np))

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
    allocate(derivvec(0:nt-1))
    do i=1,nt-1
       derivvec(i) = -0.5*(-1)**i/tan(0.5*t(i+1))
    enddo
    dbdt(:,:) = 0.0
    do j=1,np
       do i=1,nt
          do ip=1,nt
             k = ip-i
             if(k < 0) then
                k = k + nt
             endif
             dbdt(i,j) = dbdt(i,j) + derivvec(k) * bmag(ip,j)
          enddo
       enddo
    enddo
    deallocate(derivvec)
    
    ! db/dphi
    allocate(derivvec(0:nt-1))
    do i=1,np-1
       derivvec(i) = -0.5*(-1)**i/tan(0.5*p(i+1))
    enddo
    dbdp(:,:) = 0.0
    do j=1,nt
       do i=1,np
          do ip=1,np
             k = ip-i
             if(k < 0) then
                k = k + np
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

    allocate(m_indx(matsize))
    allocate(n_indx(matsize))
    allocate(itype(matsize))

    allocate(basis(nt,np))
    allocate(basis_prime(nt,np))
    allocate(basis_dt_prime(nt,np))
    allocate(basis_dp_prime(nt,np))

    allocate(mat_stream_dt(matsize,matsize))
    allocate(mat_stream_dp(matsize,matsize))
    allocate(mat_trap(matsize,matsize))
    allocate(vec_vdriftx(matsize))
    allocate(vec_flux(matsize))
    allocate(vec_upar(matsize))
    allocate(vec_uparB(matsize))
    allocate(vec_fsa(matsize))
    allocate(vec_bmag(matsize))
    allocate(vec_thetabar(matsize))

    i=0
    do ips=0,nps
       do its=0,nts
          if(its > 0) then 
             ! amn
             i=i+1
             itype(i) = 1
             m_indx(i)     = its
             n_indx(i)     = ips
             vec_thetabar(i) = as(its,ips)
          endif
          if(ips > 0 .and. its > 0) then 
             !bmn
             i=i+1
             itype(i) = 2
             m_indx(i)     = its
             n_indx(i)     = ips
             vec_thetabar(i) = bs(its,ips)
          endif
          ! cmn
          i=i+1
          itype(i) = 3
          m_indx(i)     = its
          n_indx(i)     = ips
          vec_thetabar(i) = cs(its,ips)
          if(its == 0 .and. ips == 0) then
             indx_c00 = i
          endif
          if(ips  > 0) then
             ! dmn
             i=i+1
             itype(i) = 4
             m_indx(i)     = its
             n_indx(i)     = ips
             vec_thetabar(i) = ds(its,ips)
          endif
       enddo
    enddo

    mat_stream_dt(:,:) = 0.0
    mat_stream_dp(:,:) = 0.0
    mat_trap(:,:)      = 0.0
    vec_vdriftx(:)  = 0.0
    vec_flux(:)     = 0.0
    vec_upar(:)     = 0.0
    vec_uparB(:)    = 0.0
    vec_fsa(:)      = 0.0
    vec_bmag(:)     = 0.0
    vec_thetabar(:) = 0.0

 
    do i=1,matsize
       call get_basis(itype(i),m_indx(i),n_indx(i),basis)
       do j=1,matsize
          call get_basis(itype(j),m_indx(j),n_indx(j),basis_prime(:,:))
          call get_basis_dt(itype(j),m_indx(j),n_indx(j),basis_dt_prime(:,:))
          call get_basis_dp(itype(j),m_indx(j),n_indx(j),basis_dp_prime(:,:))
          do kt=1,nt
             do kp=1,np
                mat_trap(i,j) = mat_trap(i,j) &
                     + basis(kt,kp) * basis_prime(kt,kp) &
                     * bdotgradB_overB(kt,kp) 
                mat_stream_dt(i,j) = mat_stream_dt(i,j) &
                     + basis(kt,kp) * basis_dt_prime(kt,kp) &
                     * bdotgrad(kt,kp) * iota
                mat_stream_dp(i,j) = mat_stream_dp(i,j) &
                     + basis(kt,kp) * basis_dp_prime(kt,kp) &
                     * bdotgrad(kt,kp) 
             enddo
          enddo
       enddo
       do kt=1,nt
          do kp=1,np
             vec_vdriftx(i) = vec_vdriftx(i) &
                  + basis(kt,kp) * vdrift_x(kt,kp)
             vec_flux(i)    = vec_flux(i) &
                  + basis(kt,kp) * vdrift_x(kt,kp) * g(kt,kp)
             vec_upar(i)    = vec_upar(i) &
                  + basis(kt,kp) 
             vec_uparB(i)   = vec_uparB(i) &
                  + basis(kt,kp) * bmag(kt,kp) * g(kt,kp)
             vec_fsa(i)   = vec_fsa(i) &
                  + basis(kt,kp) * g(kt,kp)
             vec_bmag(i)   = vec_bmag(i) &
                  + basis(kt,kp) * bmag(kt,kp) 
          enddo
       enddo
    enddo

    mat_stream_dt(:,:) = mat_stream_dt(:,:) / (nt*np)
    mat_stream_dp(:,:) = mat_stream_dp(:,:) / (nt*np)
    mat_trap(:,:)      = mat_trap(:,:)      / (nt*np)
    
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
    do i=1,nt
       do j=1,np
          vprime = vprime + g(i,j)
       enddo
    enddo
    vprime = vprime / (nt*np)
    

    vec_vdriftx(:) = vec_vdriftx(:) / (nt*np)
    vec_flux(:)    = vec_flux(:)  / (nt*np) / vprime
    vec_uparB(:)   = vec_uparB(:) / (nt*np) / vprime
    vec_upar(:)    = vec_upar(:)  / (nt*np) 
    vec_fsa(:)     = vec_fsa (:)  / (nt*np) / vprime
    vec_bmag(:)     = vec_bmag (:)  / (nt*np) 

    open(unit=1,file='out.le3.geovector',status='replace')
    do i=1,matsize
       write (1,'(e16.8)',advance='no') vec_thetabar(i)
       write (1,'(e16.8)',advance='no') vec_vdriftx(i)
       write (1,'(e16.8)',advance='no') vec_flux(i)
       write (1,'(e16.8)',advance='no') vec_uparB(i)
       write (1,'(e16.8)',advance='no') vec_upar(i)
       write (1,'(e16.8)',advance='no') vec_fsa(i)
       write (1,'(e16.8)',advance='no') vec_bmag(i)
       write (1,*)
    enddo
    close(1)

    open(unit=1,file='out.le3.geoscalar',status='replace')
    write (1,'(i3)') nts
    write (1,'(i3)') nps
    write (1,'(i3)') matsize
    write (1,'(i3)') indx_c00
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
    deallocate(vec_vdriftx)
    deallocate(vec_flux)
    deallocate(vec_upar)
    deallocate(vec_uparB)
    deallocate(vec_fsa)
    deallocate(vec_bmag)
    deallocate(vec_thetabar)
    deallocate(basis)
    deallocate(basis_prime)
    deallocate(basis_dt_prime)
    deallocate(basis_dp_prime)
    deallocate(itype)
    deallocate(m_indx)
    deallocate(n_indx)

20  format('(',i2,',',i2,'):',2x,4(1pe14.7,1x))
30  format(t15,4(a,9x))

  end subroutine le3_write_do
  
  subroutine get_basis(itype,m_indx,n_indx,basis)
    use le3_globals
    implicit none
    integer, intent(in)  :: itype, m_indx, n_indx
    real, dimension(nt,np), intent(inout) :: basis
    integer :: kt,kp
    
    do kt=1,nt
       do kp=1,np
          if(itype == 1) then
             basis(kt,kp) = sinm(kt,m_indx)*cosn(kp,n_indx) 
          else if(itype == 2) then
             basis(kt,kp) = sinm(kt,m_indx)*sinn(kp,n_indx)
          else if(itype == 3) then
             basis(kt,kp) = cosm(kt,m_indx)*cosn(kp,n_indx)
          else if(itype == 4) then
             basis(kt,kp) = cosm(kt,m_indx)*sinn(kp,n_indx)
          endif
       enddo
    enddo

  end subroutine get_basis

  subroutine get_basis_dt(itype,m_indx,n_indx,basis)
    use le3_globals
    implicit none
    integer, intent(in)  :: itype, m_indx, n_indx
    real, dimension(nt,np), intent(inout) :: basis
    integer :: kt,kp
    
    do kt=1,nt
       do kp=1,np
          if(itype == 1) then
             basis(kt,kp) = m_indx * cosm(kt,m_indx)*cosn(kp,n_indx) 
          else if(itype == 2) then
             basis(kt,kp) = m_indx * cosm(kt,m_indx)*sinn(kp,n_indx)
          else if(itype == 3) then
             basis(kt,kp) = -m_indx * sinm(kt,m_indx)*cosn(kp,n_indx)
          else if(itype == 4) then
             basis(kt,kp) = -m_indx * sinm(kt,m_indx)*sinn(kp,n_indx)
          endif
       enddo
    enddo

  end subroutine get_basis_dt

  subroutine get_basis_dp(itype,m_indx,n_indx,basis)
    use le3_globals
    implicit none
    integer, intent(in)  :: itype, m_indx, n_indx
    real, dimension(nt,np), intent(inout) :: basis
    integer :: kt,kp
    
    do kt=1,nt
       do kp=1,np
          if(itype == 1) then
             basis(kt,kp) = -n_indx * sinm(kt,m_indx)*sinn(kp,n_indx) 
          else if(itype == 2) then
             basis(kt,kp) = n_indx * sinm(kt,m_indx)*cosn(kp,n_indx)
          else if(itype == 3) then
             basis(kt,kp) = -n_indx * cosm(kt,m_indx)*sinn(kp,n_indx)
          else if(itype == 4) then
             basis(kt,kp) = n_indx * cosm(kt,m_indx)*cosn(kp,n_indx)
          endif
       enddo
    enddo

  end subroutine get_basis_dp

end module le3_write
