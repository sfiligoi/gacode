subroutine le3_geometry_rho

  use le3_globals
  use GEO_interface

  implicit none

  real, dimension(:,:), allocatable :: sys_m
  real, dimension(:), allocatable :: sys_b
  real, dimension(:,:), allocatable :: s1,s2
  real, dimension(:,:), allocatable :: c1,c2
  integer, dimension(:), allocatable :: i_piv
  integer :: i,j,k,kp
  integer :: nrow
  real :: g0,h0,g1,h1,d0
  real :: x,eta,chi_1
  real :: theta
  integer, parameter :: n_theta=128

  ! Since this is a 2x2 system, define new matsize
  nrow = 2*matsize

  ! Step 1: Generate required functions on (theta,phi) mesh
  allocate(s1(nt,np))
  allocate(s2(nt,np))
  allocate(c1(nt,np))
  allocate(c2(nt,np))

  do i=1,nt
     do j=1,np

        x = sqrt(gtt(i,j))

        eta = 1.0/rc(i,j)+cosu(i,j)/rs(i,j)

        ! chi_i -> chi_1/B_unit
        chi_1 = x*rs(i,j)/g(i,j)

        ! M0       -> g0
        ! N0       -> h0
        ! delta_M0 -> g1
        ! delta_N0 -> h1 

        g0 = rs(i,j)**2
        h0 = iota*x**2
        g1 = iota*x**2*2.0/rc(i,j)+chi_1*iota_p*x**2-g0*eta 
        h1 = 2*rs(i,j)*cosu(i,j)-h0*eta

        d0 = g(i,j)

        ! 4 pi p' D0 + (delta_N0+i delta_M0)/(chi_1 D0)
        s1(i,j) = -0.5*d0*beta_star/rmin+(h1+iota*g1)/(chi_1*d0)

        ! -delta_N0/(chi_1 D0)  
        s2(i,j) = -h1/(chi_1*d0)
        ! 
        c1(i,j) = g0/d0
        c2(i,j) = h0/d0

     enddo
  enddo

  allocate(sys_m(nrow,nrow))
  allocate(sys_b(nrow))
  allocate(i_piv(nrow))

  ! Initialize arrays
  sys_m(:,:) = 0.0
  sys_b(:)   = 0.0

  do k=1,matsize

     ! Projection basis
     call le3_basis(itype(k),m_indx(k),n_indx(k),basis(:,:),'d0')
     call le3_basis(itype(k),m_indx(k),n_indx(k),basis_dt(:,:),'dt')

     ! Construct source
     do i=1,nt
        do j=1,np

           ! (1)
           sys_b(k) = sys_b(k)+basis(i,j)*s1(i,j)

           ! (2)
           ! -basis_dt means integrated d/dt by parts
           sys_b(k+matsize) = sys_b(k+matsize)-basis_dt(i,j)*s2(i,j)

        enddo
     enddo

     do kp=1,matsize

        ! Expansion basis
        call le3_basis(itype(kp),m_indx(kp),n_indx(kp),basis_prime(:,:),'d0')
        call le3_basis(itype(kp),m_indx(kp),n_indx(kp),basis_dt_prime(:,:),'dt')

        ! Construct matrix
        do i=1,nt
           do j=1,np

              ! (1,1)
              sys_m(k,kp) = sys_m(k,kp)+ &
                   basis(i,j)*basis_dt_prime(i,j)*(-1)*c2(i,j)

              ! (1,2)
              sys_m(k,kp+matsize) = sys_m(k,kp+matsize)+ &
                   basis(i,j)*basis_prime(i,j)*(-1)*(c2(i,j)+iota*c1(i,j))

              ! (2,1)
              sys_m(k+matsize,kp) = sys_m(k+matsize,kp)- &
                   basis_dt(i,j)*basis_dt_prime(i,j)*c2(i,j)

              ! (2,2)
              sys_m(k+matsize,kp+matsize) = sys_m(k+matsize,kp+matsize)- &
                   basis_dt(i,j)*basis_prime(i,j)*c2(i,j)
              !endif

           enddo
        enddo

     enddo
  enddo


  ! Set average theta to zero
  !sys_m(1,:) = 0.0
  !sys_m(1,1) = 1.0
  !sys_b(1) = 0.0
  ! Fix I
  sys_m(matsize+1,:) = 0.0
  sys_m(matsize+1,1) = 1.0
  sys_b(matsize+1) = 0.0

  !do k=1,nrow
  !   print '(20(1pe11.4,1x))',sys_m(k,:)
  !enddo

  call DGESV(nrow,1,sys_m,nrow,i_piv,sys_b,nrow,info)

  print *,info

  open(unit=1,file='out.le3.rho',status='replace')
  do k=1,matsize
     write(1,'(2(1pe12.5,2x),3(i2,1x))') sys_b(k),sys_b(k+matsize),itype(k),m_indx(k)
  enddo
  close(1)

  GEO_model_in = 0
  GEO_rmin_in = rmin
  GEO_rmaj_in = rmaj
  GEO_drmaj_in = 0.0
  GEO_q_in = q
  GEO_s_in = s
  GEO_kappa_in = kappa
  GEO_beta_star_in = beta_star
  call GEO_alloc(1)
  call GEO_do()
  open(unit=1,file='out.miller',status='replace')
  do i=1,n_theta
     theta = (i-1)*2*pi/(n_theta-1)-pi
     call GEO_interp(theta)
     write(1,*) GEO_theta_s,GEO_psi2
  enddo
  close(1)
  call GEO_alloc(0)


  deallocate(sys_m)
  deallocate(sys_b)
  deallocate(i_piv)
  deallocate(basis)
  deallocate(basis_prime)
  deallocate(basis_dt_prime)
  deallocate(basis_dp_prime)
  deallocate(itype)
  deallocate(m_indx)
  deallocate(n_indx)

  deallocate(rs)
  deallocate(zs)
  deallocate(g)
  deallocate(gpp)
  deallocate(gtt)
  deallocate(gpt)
  deallocate(cosu)
  deallocate(btor)
  deallocate(bpol)
  deallocate(bmag)
  deallocate(dbdt)
  deallocate(dbdp)
  deallocate(dgdp)
  deallocate(deriv_t)
  deallocate(deriv_p)

end subroutine le3_geometry_rho
