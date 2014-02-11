subroutine le3_geometry_rho

  use le3_globals

  implicit none

  real, dimension(:,:), allocatable :: geo_m
  real, dimension(:), allocatable :: geo_b
  real, dimension(:,:), allocatable :: s1
  integer, dimension(:), allocatable :: i_piv
  integer :: i,j,kt,kp
  real :: g0,h0,g1,h1
  real :: x,eta,chi_1


  ! Step 1: Generate required functions on (theta,phi) mesh
  allocate(s1(nt,np))

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

     enddo
  enddo

  allocate(geo_m(matsize,matsize))
  allocate(geo_b(matsize))
  allocate(i_piv(matsize))

  ! Initialize arrays
  geo_m(:,:) = 0.0
  geo_b(:) = 0.0

  do i=1,matsize

     call le3_basis(itype(i),m_indx(i),n_indx(i),basis(:,:),'d0')

     ! Construct source
     do kt=1,nt
        do kp=1,np
           geo_b(i) = geo_b(i)+basis(kt,kp)*s1(kt,kp)
        enddo
     enddo

     do j=1,matsize

        call le3_basis(itype(j),m_indx(j),n_indx(j),basis_prime(:,:),'d0')

        ! Construct matrix
        do kt=1,nt
           do kp=1,np
              geo_m(i,j) = geo_m(i,j) &
                   + basis(kt,kp)* basis_prime(kt,kp) 
           enddo
        enddo

     enddo
  enddo

  call DGESV(matsize,1,geo_m,matsize,i_piv,geo_b,matsize,info)

  deallocate(geo_m)
  deallocate(geo_b)
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

end subroutine le3_geometry_rho
