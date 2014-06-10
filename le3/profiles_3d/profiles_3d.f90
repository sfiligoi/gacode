program read_m3dc1
  use netcdf
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  character(len=80)  :: path= './'
  character (len=80) :: data_file_in='input.profiles_3d.rz.nc'
  character (len=80) :: data_file
  real    :: my_psi_norm
  integer :: nt_fourier
  real    :: a_meters
  integer :: impurity_flag
  integer :: zimp
  integer :: rztest_model
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: np_fourier=16 
  integer :: np_mask(0:np_fourier)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real, parameter :: pi=3.1415926535897932
  real, parameter :: temp_norm_fac   = 1602.2
  real, parameter :: charge_norm_fac = 1.6022
  real, parameter :: mass_deuterium = 3.3452   ! (x 10-27 kg)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: nr, nt, np, npsi
  integer :: verbose_flag=0
  integer :: ir, ip, it, jp, im, n
  real :: d_phi
  real :: m3d_psi0, m3d_psi1
  real, dimension(:,:,:), allocatable :: m3d_bigr, m3d_bigz
  real, dimension(:,:), allocatable :: bigr_2d, bigz_2d
  real, dimension(:), allocatable :: m3d_phi, m3d_psi, m3d_q
  integer :: nr_prof
  real, dimension(:), allocatable :: m3d_psi_prof, m3d_te, m3d_ti, m3d_tc, &
       m3d_ne, m3d_ni, m3d_nc
  real, dimension(:), allocatable :: te, ti, tc, ne, ni, nc, &
       dlntedr, dlntidr, dlntcdr, dlnnedr, dlnnidr, dlnncdr
  real, dimension(:), allocatable :: psi_norm, psi_norm_prof
  real, dimension(:,:), allocatable :: rmin, rmaj
  real, dimension(:), allocatable :: rmin_avg, rmaj_avg, dpsidr, bunit
  real, dimension(:,:,:,:), allocatable :: gvec_3d, gvec_miller
  real, dimension(:,:,:), allocatable :: gvec_2d
  real, dimension(:,:,:,:), allocatable :: amn_r, amn_z, amn_rd, amn_zd
  character (len=70) :: cdum
  integer :: ncid, varid, err
  ! local parameters
  integer :: nr_loc=1
  real, dimension(:,:,:,:), allocatable :: amn_r_loc, amn_z_loc, amn_rd_loc, &
       amn_zd_loc
  real, dimension(1) :: psi_norm_loc, rmin_avg_loc, rmaj_avg_loc,&
       q_loc, bunit_loc
  real, dimension(1) :: te_loc, ti_loc, tc_loc, ne_loc, ni_loc, nc_loc, &
       dlntedr_loc, dlntidr_loc, dlntcdr_loc, &
       dlnnedr_loc, dlnnidr_loc, dlnncdr_loc
  real, dimension(1) :: rho_loc, nu_loc
  real :: cc, loglam, vth_norm, sum, x, y, dtheta
  real, dimension(:,:), allocatable :: theta

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read the input.profiles_3d 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit=1,file=trim(path)//'input.profiles_3d.gen',status='old')
  read(1,*) my_psi_norm
  read(1,*) nt_fourier
  read(1,*) a_meters
  read(1,*) impurity_flag
  read(1,*) zimp
  read(1,*) rztest_model
  do ip=0, np_fourier
     read(1,*) np_mask(ip)
  enddo
  close(1)

  data_file=trim(path)//data_file_in

  !!!!!!!!!!!!!!!!!!!!!!!!
  ! Geometry coefficients
  !!!!!!!!!!!!!!!!!!!!!!!!

  psi_norm_loc(1) = my_psi_norm

  ! Read the data from the netcdf file
 
  err = nf90_open(data_file,NF90_NOWRITE,ncid)

  ! nr=n_radial
  err = nf90_inq_dimid(ncid,trim('nr'),varid)
  err = nf90_Inquire_Dimension(ncid,varid,len=nr)
  ! 'nt'=num toroidal=np
  err = nf90_inq_dimid(ncid,trim('nt'),varid)
  err = nf90_Inquire_Dimension(ncid,varid,len=np)
  ! 'np'=num poloidal=nt
  err = nf90_inq_dimid(ncid,trim('np'),varid)
  err = nf90_Inquire_Dimension(ncid,varid,len=nt)

  allocate(m3d_bigr(nr,np,nt+1))
  allocate(m3d_bigz(nr,np,nt+1))
  allocate(m3d_phi(np))
  allocate(m3d_psi(nr))
  allocate(m3d_q(nr))

  err = nf90_inq_varid(ncid,trim('R'),varid)
  err = nf90_get_var(ncid,varid,m3d_bigr(:,:,1:nt))
  err = nf90_inq_varid(ncid,trim('Z'),varid)
  err = nf90_get_var(ncid,varid,m3d_bigz(:,:,1:nt))
  err = nf90_inq_varid(ncid,trim('Phi'),varid)
  err = nf90_get_var(ncid,varid,m3d_phi(:))
  err = nf90_get_att(ncid,nf90_global,trim('psi_0'),m3d_psi0)
  err = nf90_get_att(ncid,nf90_global,trim('psi_1'),m3d_psi1)
  err = nf90_inq_varid(ncid,trim('psi'),varid)
  err = nf90_get_var(ncid,varid,m3d_psi(:))
  err = nf90_inq_varid(ncid,trim('q'),varid)
  err = nf90_get_var(ncid,varid,m3d_q(:))

  err = nf90_close(ncid)

  ! Convert phi to radians
  do ip=1,np
     m3d_phi(ip) = m3d_phi(ip)*pi/180.0
  enddo
  d_phi=2*pi/np

  allocate(psi_norm(nr))
  psi_norm(:) = (m3d_psi(:)-m3d_psi0)/(m3d_psi1-m3d_psi0)
  if(psi_norm_loc(1) > psi_norm(nr) .or. psi_norm_loc(1) < psi_norm(1)) then
     print *, "Error: must chose psi_norm in domain: ", psi_norm(1),",", &
          psi_norm(nr)
     stop
  endif

  print *, nr, nt, np

  ! for fluxfit, the first and last theta points must be the same
  do ir=1,nr
     do ip=1,np
        m3d_bigr(ir,ip,nt+1) = m3d_bigr(ir,ip,1)
        m3d_bigz(ir,ip,nt+1) = m3d_bigz(ir,ip,1)
     enddo
  enddo
  
  if (rztest_model == 1) then
     do ir=1,nr
        do it=1,nt+1
           x = (it-1)*2*pi/nt
           do ip=1,np
              x = (it-1)*2*pi/nt
              x = x + 0.1*sin(x)
              m3d_bigr(ir,ip,it) = 3.0 + 0.1*cos(x) &
                   + 0.001*cos(2*x-3*m3d_phi(ip))
              m3d_bigz(ir,ip,it) = 1.5*0.1*sin(x) &
                   + 0.001*sin(2*x-3*m3d_phi(ip))
           end do
        end do
     end do
  end if

  allocate(rmin(nr,np))
  allocate(rmaj(nr,np))
  allocate(gvec_miller(6,1,nr,np))
  ! first with Miller to get rmin
  do ir=1,nr
     do ip=1,np
        call fluxfit_driver(1,1,1,nt+1,&
             m3d_bigr(ir,ip,:),m3d_bigz(ir,ip,:),verbose_flag)
        open(unit=1,file=trim(path)//'fluxfit.profile',status='old')
        read(1,*) cdum
        read(1,*) gvec_miller(:,1,ir,ip)
        close(1)
        rmin(ir,ip) = gvec_miller(1,1,ir,ip)
        rmaj(ir,ip) = gvec_miller(3,1,ir,ip)
     enddo
  enddo

  ! define the minor radius as the toroidal average of rmin
  allocate(rmin_avg(nr))
  rmin_avg(:) = 0.0
  do ir=1,nr
     do ip=1,np
        rmin_avg(ir) = rmin_avg(ir) + rmin(ir,ip)
     enddo
     rmin_avg(ir) = rmin_avg(ir) / (1.0*np)
  enddo

  allocate(rmaj_avg(nr))
  rmaj_avg(:) = 0.0
  do ir=1,nr
     do ip=1,np
        rmaj_avg(ir) = rmaj_avg(ir) + rmaj(ir,ip)
     enddo
     rmaj_avg(ir) = rmaj_avg(ir) / (1.0*np)
  enddo

  ! dpsi/dr and bunit
  allocate(dpsidr(nr))
  allocate(bunit(nr))
  call bound_deriv(dpsidr,m3d_psi,rmin_avg,nr)
  bunit(:) = m3d_q(:)/rmin_avg(:)*dpsidr(:)

  ! axisymmetric component -- avg R and Z over phi
  allocate(bigr_2d(nr,nt+1))
  allocate(bigz_2d(nr,nt+1))
  do ir=1,nr
     do it=1,nt+1
        ! R
        sum = 0.0
        do ip=1,np
           sum = sum + m3d_bigr(ir,ip,it)
        enddo
        sum = sum/(1.0*np)
        bigr_2d(ir,it) = sum
        ! Z
        sum = 0.0
        do ip=1,np
           sum = sum + m3d_bigz(ir,ip,it)
        enddo
        sum = sum/(1.0*np)
        bigz_2d(ir,it) = sum
     enddo
  enddo

  ! general geometry for axisymmetric to get theta and m fourier coefficients
  allocate(gvec_2d(4,0:nt_fourier,nr))
  allocate(theta(nr,nt+1))
  do ir=1,nr
     call fluxfit_driver(2,nt_fourier,1,nt+1,&
          bigr_2d(ir,:),bigz_2d(ir,:),verbose_flag)
     open(unit=1,file=trim(path)//'fluxfit.geo',status='old')
     read(1,*) it
     read(1,*) gvec_2d(:,:,ir)
     close(1)
     open(unit=1,file=trim(path)//'fluxfit.theta',status='old')
     read(1,*) theta(ir,:)
     close(1)
  enddo
  
  open(unit=1,file=trim(path)//'out.profiles_3d.bigr',status='replace')
  do ip=1,np
     do it=1,nt+1
        write(1,*) m3d_bigr(1,ip,it)
     enddo
  enddo
  close(1)

  ! compute decomposition for non-axisymmetic components

  do ir=1,nr
     do ip=1,np
        do it=1,nt+1
           m3d_bigr(ir,ip,it)=m3d_bigr(ir,ip,it)-bigr_2d(ir,it)
           m3d_bigz(ir,ip,it)=m3d_bigz(ir,ip,it)-bigz_2d(ir,it)
        enddo
     enddo
  enddo

  allocate(gvec_3d(4,0:nt_fourier,nr,np))
  gvec_3d(:,:,:,:) = 0.0
  do ir=1,nr
     do ip=1,np
        do it=1,nt
           dtheta = theta(ir,it+1)-theta(ir,it)
           !if(abs(dtheta) > pi ) dtheta = dtheta+2*pi
           if (dtheta > pi) dtheta = dtheta-2*pi
           if (dtheta < -pi) dtheta = dtheta+2*pi
           !if(ir==1 .and. ip==1) then
           !   print *, dtheta
           !endif
           do n=0,nt_fourier
              y = 0.5*(cos(n*theta(ir,it+1))*m3d_bigr(ir,ip,it+1)&
                   +cos(n*theta(ir,it))*m3d_bigr(ir,ip,it))
              gvec_3d(1,n,ir,ip) = gvec_3d(1,n,ir,ip)+dtheta*y/pi
                 
              y = 0.5*(sin(n*theta(ir,it+1))*m3d_bigr(ir,ip,it+1)&
                   +sin(n*theta(ir,it))*m3d_bigr(ir,ip,it))
              gvec_3d(2,n,ir,ip) = gvec_3d(2,n,ir,ip)+dtheta*y/pi
                 
              y = 0.5*(cos(n*theta(ir,it+1))*m3d_bigz(ir,ip,it+1)&
                   +cos(n*theta(ir,it))*m3d_bigz(ir,ip,it))
              gvec_3d(3,n,ir,ip) = gvec_3d(3,n,ir,ip)+dtheta*y/pi
                 
              y = 0.5*(sin(n*theta(ir,it+1))*m3d_bigz(ir,ip,it+1)&
                   +sin(n*theta(ir,it))*m3d_bigz(ir,ip,it))
              gvec_3d(4,n,ir,ip) = gvec_3d(4,n,ir,ip)+dtheta*y/pi
                 
           enddo
        enddo
     enddo
  enddo

  do ip=1,np
     gvec_3d(:,:,:,ip) = gvec_3d(:,:,:,ip) + gvec_2d(:,:,:)
  enddo

  allocate(amn_r(1:nr,0:nt_fourier,0:np_fourier,4))
  allocate(amn_z(1:nr,0:nt_fourier,0:np_fourier,4))
  allocate(amn_rd(1:nr,0:nt_fourier,0:np_fourier,4))
  allocate(amn_zd(1:nr,0:nt_fourier,0:np_fourier,4))

  ! compute the 3D fourier coefficients
  amn_r(:,:,:,:) = 0.0
  amn_z(:,:,:,:) = 0.0
  ! case(1): amn sin(m theta) cos(n phi)
  ! 0 if m=0
  do ir=1,nr
     do it=1,nt_fourier
        do ip=0,np_fourier
           do jp=1,np
              amn_r(ir,it,ip,1) = amn_r(ir,it,ip,1) + gvec_3d(2,it,ir,jp) &
                   * cos(ip*m3d_phi(jp))
              amn_z(ir,it,ip,1) = amn_z(ir,it,ip,1) + gvec_3d(4,it,ir,jp) &
                   * cos(ip*m3d_phi(jp))
           enddo
           amn_r(ir,it,ip,1) = amn_r(ir,it,ip,1) * d_phi/pi
           amn_z(ir,it,ip,1) = amn_z(ir,it,ip,1) * d_phi/pi
           if(ip==0) then
              amn_r(ir,it,ip,1) = amn_r(ir,it,ip,1) * 0.5
              amn_z(ir,it,ip,1) = amn_z(ir,it,ip,1) * 0.5
           endif
        enddo
     enddo
  enddo
  ! case(2): bmn sin(m theta) sin(n phi)
  ! 0 if m=0 or n=0
  do ir=1,nr
     do it=1,nt_fourier
        do ip=1,np_fourier
           do jp=1,np
              amn_r(ir,it,ip,2) = amn_r(ir,it,ip,2) + gvec_3d(2,it,ir,jp) &
                   *sin(ip*m3d_phi(jp))
              amn_z(ir,it,ip,2) = amn_z(ir,it,ip,2) + gvec_3d(4,it,ir,jp) &
                   *sin(ip*m3d_phi(jp))
           enddo
           amn_r(ir,it,ip,2) = amn_r(ir,it,ip,2) * d_phi/pi
           amn_z(ir,it,ip,2) = amn_z(ir,it,ip,2) * d_phi/pi
        enddo
     enddo
  enddo
  ! case(3): cmn cos(m theta) cos(n phi)
  do ir=1,nr
     do it=0,nt_fourier
        do ip=0,np_fourier
           do jp=1,np
              amn_r(ir,it,ip,3) = amn_r(ir,it,ip,3) + gvec_3d(1,it,ir,jp) &
                   * cos(ip*m3d_phi(jp))
              amn_z(ir,it,ip,3) = amn_z(ir,it,ip,3) + gvec_3d(3,it,ir,jp) &
                   * cos(ip*m3d_phi(jp))
           enddo
           amn_r(ir,it,ip,3) = amn_r(ir,it,ip,3) * d_phi/pi
           amn_z(ir,it,ip,3) = amn_z(ir,it,ip,3) * d_phi/pi
           if(ip==0) then
              amn_r(ir,it,ip,3) = amn_r(ir,it,ip,3) * 0.5
              amn_z(ir,it,ip,3) = amn_z(ir,it,ip,3) * 0.5
           endif
           if(it==0) then
              ! gvec(it=0) is really a0/2
              amn_r(ir,it,ip,3) = amn_r(ir,it,ip,3) * 0.5
              amn_z(ir,it,ip,3) = amn_z(ir,it,ip,3) * 0.5
           endif
        enddo
     enddo
  enddo
  ! case(4): dmn cos(m theta) sin(n phi)
  ! 0 if n=0
  do ir=1,nr
     do it=0,nt_fourier
        do ip=1,np_fourier
           do jp=1,np
              amn_r(ir,it,ip,4) = amn_r(ir,it,ip,4) + gvec_3d(1,it,ir,jp) &
                   * sin(ip*m3d_phi(jp))
              amn_z(ir,it,ip,4) = amn_z(ir,it,ip,4) + gvec_3d(3,it,ir,jp) &
                   * sin(ip*m3d_phi(jp))
           enddo
           amn_r(ir,it,ip,4) = amn_r(ir,it,ip,4) * d_phi/pi
           amn_z(ir,it,ip,4) = amn_z(ir,it,ip,4) * d_phi/pi
           if(it==0) then
              ! gvec(it=0) is really a0/2
              amn_r(ir,it,ip,4) = amn_r(ir,it,ip,4) * 0.5
              amn_z(ir,it,ip,4) = amn_z(ir,it,ip,4) * 0.5
           endif
        enddo
     enddo
  enddo

  ! radial derivatives
  if(rztest_model == 1) then
     amn_rd = 0.0
     amn_zd = 0.0
  else
     do it=0,nt_fourier
        do ip=0,np_fourier
           do im=1,4
              call bound_deriv(amn_rd(:,it,ip,im),amn_r(:,it,ip,im),&
                   rmin_avg,nr)
              call bound_deriv(amn_zd(:,it,ip,im),amn_z(:,it,ip,im),&
                   rmin_avg,nr)
           enddo
        enddo
     enddo
  endif

  ! interpolate to get the local value
  allocate(amn_r_loc(1,0:nt_fourier,0:np_fourier,4))
  allocate(amn_z_loc(1,0:nt_fourier,0:np_fourier,4))
  allocate(amn_rd_loc(1,0:nt_fourier,0:np_fourier,4))
  allocate(amn_zd_loc(1,0:nt_fourier,0:np_fourier,4))

  do it=0,nt_fourier
     do ip=0,np_fourier
        do im=1,4
           call cub_spline(psi_norm,amn_r(:,it,ip,im),nr,&
                psi_norm_loc,amn_r_loc(:,it,ip,im),nr_loc)
           call cub_spline(psi_norm,amn_z(:,it,ip,im),nr,&
                psi_norm_loc,amn_z_loc(:,it,ip,im),nr_loc)
           call cub_spline(psi_norm,amn_rd(:,it,ip,im),nr,&
                psi_norm_loc,amn_rd_loc(:,it,ip,im),nr_loc)
           call cub_spline(psi_norm,amn_zd(:,it,ip,im),nr,&
                psi_norm_loc,amn_zd_loc(:,it,ip,im),nr_loc)
        enddo
     enddo
  enddo

  call cub_spline(psi_norm,rmin_avg,nr,psi_norm_loc,rmin_avg_loc,nr_loc)
  call cub_spline(psi_norm,rmaj_avg,nr,psi_norm_loc,rmaj_avg_loc,nr_loc)
  call cub_spline(psi_norm,m3d_q,nr,psi_norm_loc,q_loc,nr_loc)
  call cub_spline(psi_norm,bunit,nr,psi_norm_loc,bunit_loc,nr_loc)
  

  !!!!!!!!!!!!!!!!!!!!!!!!
  ! Profile data
  !!!!!!!!!!!!!!!!!!!!!!!!

  ! Read the data from the netcdf file
  err = nf90_open(data_file,NF90_NOWRITE,ncid)
  err = nf90_inq_dimid(ncid,trim('npsi'),varid)
  err = nf90_Inquire_Dimension(ncid,varid,len=nr_prof)

  allocate(m3d_te(nr_prof))
  allocate(m3d_ti(nr_prof))
  allocate(m3d_tc(nr_prof))
  allocate(m3d_ne(nr_prof))
  allocate(m3d_ni(nr_prof))
  allocate(m3d_nc(nr_prof))

  ! ele temp -- eV
  err = nf90_inq_varid(ncid,trim('Te0'),varid)
  err = nf90_get_var(ncid,varid,m3d_te(:))
  ! ion temp -- eV
  err = nf90_inq_varid(ncid,trim('Ti0'),varid)
  err = nf90_get_var(ncid,varid,m3d_ti(:))
  ! ele dens -- 1/m^3
  err = nf90_inq_varid(ncid,trim('ne0'),varid)
  err = nf90_get_var(ncid,varid,m3d_ne(:))
  ! ion dens -- 1/m^3
  err = nf90_inq_varid(ncid,trim('ni0'),varid)
  err = nf90_get_var(ncid,varid,m3d_ni(:))

  ! put in units of input.profiles (keV, 10^19/m^3)
  m3d_te(:) = m3d_te(:) / 1000.0
  m3d_ti(:) = m3d_ti(:) / 1000.0
  m3d_ne(:) = m3d_ne(:) / 1.0e19
  m3d_ni(:) = m3d_ni(:) / 1.0e19

  ! assume missing density is specified impurity
  m3d_tc(:) = m3d_ti(:)
  m3d_nc(:) = (m3d_ne(:)-m3d_ni(:))/(1.0*zimp)
  if(impurity_flag /= 1) then
     m3d_ni(:) = m3d_ne(:)
  endif

  allocate(m3d_psi_prof(nr_prof))
  err = nf90_inq_varid(ncid,trim('psi0'),varid)
  err = nf90_get_var(ncid,varid,m3d_psi_prof(:))

  allocate(psi_norm_prof(nr_prof))
  psi_norm_prof(:) = (m3d_psi_prof(:)-m3d_psi0)/(m3d_psi1-m3d_psi0)

  err = nf90_close(ncid)

  ! make a pfile
  open(unit=1,file=trim(path)//'out.profiles_3d.peq',status='replace')
  write(1,*)  nr_prof, 'psinorm  ne(10^20/m^3) dne/dpsiN'
  write(1,*)  psi_norm_prof(1)*0.0, m3d_ne(1)/10, 0.0
  do ir=1,nr_prof
     write(1,*) psi_norm_prof(ir), m3d_ne(ir)/10, 0.0
  enddo
  write(1,*)  nr_prof, 'psinorm  te(KeV) dte/dpsiN'
  write(1,*)  psi_norm_prof(1)*0.0, m3d_te(1), 0.0
  do ir=1,nr_prof
     write(1,*) psi_norm_prof(ir), m3d_te(ir), 0.0
  enddo
  write(1,*)  nr_prof, 'psinorm  ni(10^20/m^3) dni/dpsiN'
  write(1,*)  psi_norm_prof(1)*0.0, m3d_ni(1)/10, 0.0
  do ir=1,nr_prof
     write(1,*) psi_norm_prof(ir), m3d_ni(ir)/10, 0.0
  enddo
  write(1,*)  nr_prof, 'psinorm  ti(KeV) dti/dpsiN'
  write(1,*)  psi_norm_prof(1)*0.0, m3d_ti(1), 0.0
  do ir=1,nr_prof
     write(1,*) psi_norm_prof(ir), m3d_ti(ir), 0.0
  enddo
  close(1)

  ! spline onto geo grid
  allocate(te(nr))
  allocate(ti(nr))
  allocate(tc(nr))
  allocate(ne(nr))
  allocate(ni(nr))
  allocate(nc(nr))
  call cub_spline(psi_norm_prof,m3d_te,nr_prof,psi_norm,te,nr)
  call cub_spline(psi_norm_prof,m3d_ti,nr_prof,psi_norm,ti,nr)
  call cub_spline(psi_norm_prof,m3d_tc,nr_prof,psi_norm,tc,nr)
  call cub_spline(psi_norm_prof,m3d_ne,nr_prof,psi_norm,ne,nr)
  call cub_spline(psi_norm_prof,m3d_ni,nr_prof,psi_norm,ni,nr)
  call cub_spline(psi_norm_prof,m3d_nc,nr_prof,psi_norm,nc,nr)

  ! compute the gradient length scales
  allocate(dlntedr(nr))
  allocate(dlntidr(nr))
  allocate(dlntcdr(nr))
  allocate(dlnnedr(nr))
  allocate(dlnnidr(nr))
  allocate(dlnncdr(nr))
  call bound_deriv(dlntedr,te,rmin_avg,nr)
  dlntedr(:) = -1.0*dlntedr(:)/te(:)
  call bound_deriv(dlntidr,ti,rmin_avg,nr)
  dlntidr(:) = -1.0*dlntidr(:)/ti(:)
  call bound_deriv(dlntcdr,tc,rmin_avg,nr)
  dlntcdr(:) = -1.0*dlntcdr(:)/tc(:)
  call bound_deriv(dlnnedr,ne,rmin_avg,nr)
  dlnnedr(:) = -1.0*dlnnedr(:)/ne(:)
  call bound_deriv(dlnnidr,ni,rmin_avg,nr)
  dlnnidr(:) = -1.0*dlnnidr(:)/ni(:)
  call bound_deriv(dlnncdr,nc,rmin_avg,nr)
  do ir=1,nr
     if(abs(nc(ir)) > 0.0) then
        dlnncdr(:) = -1.0*dlnncdr(:)/nc(:)
     else
        dlnncdr(:) = 0.0
     endif
  enddo
        
  ! obtain local values of profile quantities

  call cub_spline(psi_norm,te,nr,psi_norm_loc,te_loc,nr_loc)
  call cub_spline(psi_norm,ti,nr,psi_norm_loc,ti_loc,nr_loc)
  call cub_spline(psi_norm,tc,nr,psi_norm_loc,tc_loc,nr_loc)
  call cub_spline(psi_norm,ne,nr,psi_norm_loc,ne_loc,nr_loc)
  call cub_spline(psi_norm,ni,nr,psi_norm_loc,ni_loc,nr_loc)
  call cub_spline(psi_norm,nc,nr,psi_norm_loc,nc_loc,nr_loc)
  call cub_spline(psi_norm,dlntedr,nr,psi_norm_loc,dlntedr_loc,nr_loc)
  call cub_spline(psi_norm,dlntidr,nr,psi_norm_loc,dlntidr_loc,nr_loc)
  call cub_spline(psi_norm,dlntcdr,nr,psi_norm_loc,dlntcdr_loc,nr_loc)
  call cub_spline(psi_norm,dlnnedr,nr,psi_norm_loc,dlnnedr_loc,nr_loc)
  call cub_spline(psi_norm,dlnnidr,nr,psi_norm_loc,dlnnidr_loc,nr_loc)
  call cub_spline(psi_norm,dlnncdr,nr,psi_norm_loc,dlnncdr_loc,nr_loc)

  ! derived quantities
  !a_meters = rmin_avg(nr)
  ! vth/a (1/s)
  vth_norm = sqrt(ti_loc(1) * temp_norm_fac &
       / (mass_deuterium)) &
       * 1.0e4 / a_meters

  ! rho
  rho_loc(1) = sqrt(ti_loc(1) * temp_norm_fac &
          * mass_deuterium) &
          / (charge_norm_fac * bunit_loc(1)) &
          * 1.0e-4 / a_meters

  ! collision frequency (main deuterium ion)
  cc = sqrt(2.0) * pi * charge_norm_fac**4 &
       * 1.0 / (4.0 * pi * 8.8542)**2 &
       * 1.0 / (sqrt(mass_deuterium) * temp_norm_fac**1.5) &
       * 1e9
  loglam = 24.0 - log(sqrt(ne_loc(1)*1e13)/(te_loc(1)*1000))
  nu_loc = cc * loglam * ni_loc / (ti_loc**1.5)
  nu_loc = nu_loc/vth_norm

  open(unit=1,file=trim(path)//'out.profiles_3d.geoall',status='replace')
  write(1,*) nt_fourier
  write(1,*) np_fourier
  do it=0,nt_fourier
     do ip=0,np_fourier
        do im=1,4
           write(1,'(1pe12.5)') amn_r_loc(1,it,ip,im)/a_meters
           write(1,'(1pe12.5)') amn_z_loc(1,it,ip,im)/a_meters
           write(1,'(1pe12.5)') amn_rd_loc(1,it,ip,im)
           write(1,'(1pe12.5)') amn_zd_loc(1,it,ip,im)
        enddo
     enddo
  enddo
  close(1)

  open(unit=1,file=trim(path)//'out.profiles_3d.geo',status='replace')
  write(1,*) nt_fourier
  write(1,*) np_fourier
  do it=0,nt_fourier
     do ip=0,np_fourier
        do im=1,4
           write(1,'(1pe12.5)') np_mask(ip)*amn_r_loc(1,it,ip,im)/a_meters
           write(1,'(1pe12.5)') np_mask(ip)*amn_z_loc(1,it,ip,im)/a_meters
           write(1,'(1pe12.5)') np_mask(ip)*amn_rd_loc(1,it,ip,im)
           write(1,'(1pe12.5)') np_mask(ip)*amn_zd_loc(1,it,ip,im)
        enddo
     enddo
  enddo
  close(1)

  open(unit=1,file=trim(path)//'out.profiles_3d.prof',status='replace')
  write (1,*)  'Data'
  write (1,20) 'geo psi_norm min=', psi_norm(1)
  write (1,20) 'geo psi_norm max=', psi_norm(nr)
  write (1,20) 'local psi_norm=', psi_norm_loc(1)
  write (1,20) 'rmin (m)=',   rmin_avg_loc(1)
  write (1,20) 'q=',      q_loc(1)
  write (1,20) 'bunit (T)=',  bunit_loc
  write (1,20) 'te (keV)=',     te_loc(1)
  write (1,20) 'ti (keV)=',     ti_loc(1)
  write (1,20) 'tc (keV)=',     tc_loc(1)
  write (1,20) 'ne (10^19/m^3)=',     ne_loc(1)
  write (1,20) 'ni (10^19/m^3)=',     ni_loc(1)
  write (1,20) 'nc (10^19/m^3)=',     nc_loc(1)
  write (1,20) 'dlntedr (1/m)=',     dlntedr_loc(1)
  write (1,20) 'dlntidr (1/m)=',     dlntidr_loc(1)
  write (1,20) 'dlntcdr (1/m)=',     dlntcdr_loc(1)
  write (1,20) 'dlnnedr (1/m)=',     dlnnedr_loc(1)
  write (1,20) 'dlnnidr (1/m)=',     dlnnidr_loc(1)
  write (1,20) 'dlnncdr (1/m)=',     dlnncdr_loc(1)
  write (1,*)
  write (1,*)  'Parameters for LE3'
  write (1,20) 'RMIN=', rmin_avg_loc(1)/a_meters
  write (1,20) 'RMAJ=', rmaj_avg_loc(1)/a_meters
  write (1,20) 'Q=', abs(q_loc(1))
  write (1,*)
  write (1,*)  'Parameters for NEO'
  write (1,20) 'RMIN_OVER_A=', rmin_avg_loc(1)/a_meters
  write (1,20) 'RHO_STAR=', abs(rho_loc(1))
  write (1,*)
  write (1,30) 'N_SPECIES=', 2+impurity_flag
  write (1,30) 'Z_1=',1
  write (1,20) 'MASS_1=',1.0
  write (1,20) 'DENS_1=',ni_loc(1)/ni_loc(1)
  write (1,20) 'TEMP_1=',ti_loc(1)/ti_loc(1)
  write (1,20) 'DLNNDR_1=',dlnnidr_loc(1)*a_meters
  write (1,20) 'DLNTDR_1=',dlntidr_loc(1)*a_meters
  write (1,20) 'NU_1=',nu_loc(1)
  write (1,*)
  write (1,30) 'Z_2=',-1
  write (1,20) 'MASS_2=',0.0002724486
  write (1,20) 'DENS_2=',ne_loc(1)/ni_loc(1)
  write (1,20) 'TEMP_2=',te_loc(1)/ti_loc(1)
  write (1,20) 'DLNNDR_2=',dlnnedr_loc(1)*a_meters
  write (1,20) 'DLNTDR_2=',dlntedr_loc(1)*a_meters
  write (1,*)
  write (1,30) 'Z_3=',zimp
  write (1,20) 'MASS_3=',zimp*1.0
  write (1,20) 'DENS_3=',nc_loc(1)/ni_loc(1)
  write (1,20) 'TEMP_3=',tc_loc(1)/ti_loc(1)
  write (1,20) 'DLNNDR_3=',dlnncdr_loc(1)*a_meters
  write (1,20) 'DLNTDR_3=',dlntcdr_loc(1)*a_meters
  write (1,*)  
  write (1,*)  'Normalizations'
  write (1,20) 'a_meters=', a_meters
  write (1,20) 'n_norm=', ne_loc(1)
  write (1,20) 't_norm=', ti_loc(1)
  write (1,20) 'vth_norm=', vth_norm*a_meters
  write (1,20) 'b_unit=', bunit(1)
  close(1)

20  format(t2,a,1pe12.5)
30  format(t2,a,i3)

end program read_m3dc1
