subroutine le3_compute_theory

  use le3_globals

  implicit none

  real, dimension(:,:), allocatable :: bgradmat
  real, dimension(:), allocatable :: vec
  real, dimension(:,:), allocatable :: a, ahat, gfac, u, uhat, h, hhat
  integer :: i, j, kt, kp
  integer, dimension(:), allocatable :: i_piv
  real :: bsq_avg, bgradbsq_avg, gfac_avg, bsqahat_avg
  real :: g2fac, g1fac, g3fac, g1fac_ntv
  real :: gs1, gs2, gs3, gs4, gs5, gs6
  real, dimension(:,:), allocatable :: ubfac

  bsq_avg      = 0.0
  bgradbsq_avg = 0.0
  do i=1,nt
     do j=1,np
        bsq_avg      = bsq_avg + g(i,j) * bmag(i,j)**2
        bgradbsq_avg = bgradbsq_avg + g(i,j) &
             * (bdotgradB_overB(i,j) * bmag(i,j))**2
     enddo
  enddo
  bsq_avg      = bsq_avg / (nt*np) / vprime
  bgradbsq_avg = bgradbsq_avg / (nt*np) / vprime

  ! b dot grad matrix
  allocate(bgradmat(matsize,matsize))
  allocate(i_piv(matsize))
  allocate(vec(matsize))
  bgradmat(:,:) = mat_stream_dt(:,:) + mat_stream_dp(:,:)
  bgradmat(indx_c00,:) = 0.0 
  bgradmat(indx_c00,indx_c00) = 1.0
  call DGETRF(matsize,matsize,bgradmat(:,:),matsize,i_piv,info)

  ! Parameters for our PS theory

  ! compute ahat
  allocate(a(nt,np))
  do kt=1,nt
     do kp=1,np
        a(kt,kp) = vdrift_x(kt,kp)/ bmag(kt,kp)
     enddo
  enddo
  vec(:) = 0.0
  do i=1,matsize
     call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'d0')
     do kt=1,nt
        do kp=1,np
           vec(i) = vec(i) + bk(kt,kp) * a(kt,kp)
        enddo
     enddo
  enddo
  vec = vec / (nt*np)
  vec(indx_c00) = 0.0
  call DGETRS('N',matsize,1,bgradmat(:,:),matsize,i_piv,vec(:),matsize,info)
  allocate(ahat(nt,np))
  ahat(:,:) = 0.0
  do i=1,matsize
     call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'d0')
     do kt=1,nt
        do kp=1,np
           ahat(kt,kp) = ahat(kt,kp) + bk(kt,kp) * vec(i)
        enddo
     enddo
  enddo

  allocate(gfac(nt,np))
  gfac = ahat * (bdotgradB_overB * bmag)**2 &
       + 0.5 * a * bmag * (bdotgradB_overB * bmag)

  gfac_avg = 0.0
  do i=1,nt
     do j=1,np
        gfac_avg = gfac_avg + g(i,j) * gfac(i,j)
     enddo
  enddo
  gfac_avg = gfac_avg / (nt*np) / vprime

  bsqahat_avg = 0.0
  do i=1,nt
     do j=1,np
        bsqahat_avg = bsqahat_avg  + g(i,j) * bmag(i,j)**2 * ahat(i,j)
     enddo
  enddo
  bsqahat_avg = bsqahat_avg / (nt*np) / vprime
  g2fac = bsqahat_avg - gfac_avg * bsq_avg/bgradbsq_avg

  ! compute uhat
  allocate(u(nt,np))
  do kt=1,nt
     do kp=1,np
        u(kt,kp) = 1/bmag(kt,kp) * (gfac(kt,kp) &
             - (bdotgradB_overB(kt,kp) * bmag(kt,kp))**2 &
             / bgradbsq_avg * gfac_avg)
     enddo
  enddo
  vec(:) = 0.0
  do i=1,matsize
     call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'d0')
     do kt=1,nt
        do kp=1,np
           vec(i) = vec(i) + bk(kt,kp) * u(kt,kp)
        enddo
     enddo
  enddo
  vec = vec / (nt*np)
  vec(indx_c00) = 0.0
  call DGETRS('N',matsize,1,bgradmat(:,:),matsize,i_piv,vec(:),matsize,info)
  allocate(uhat(nt,np))
  uhat(:,:) = 0.0
  do i=1,matsize
     call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'d0')
     do kt=1,nt
        do kp=1,np
           uhat(kt,kp) = uhat(kt,kp) + bk(kt,kp) * vec(i)
        enddo
     enddo
  enddo

  g1fac=0.0
  do i=1,nt
     do j=1,np
        g1fac = g1fac + g(i,j) * bmag(i,j) * a(i,j) &
             * (u(i,j)*bmag(i,j)/(bdotgradB_overB(i,j) * bmag(i,j)) &
             - 2.0*uhat(i,j))
     enddo
  enddo
  g1fac = g1fac / (nt*np) / vprime

  g1fac_ntv=0.0
  do i=1,nt
     do j=1,np
        g1fac_ntv = g1fac_ntv + dgdp(i,j) &
             * (u(i,j)*bmag(i,j)/(bdotgradB_overB(i,j) * bmag(i,j)) &
             - 2.0*uhat(i,j))
     enddo
  enddo
  g1fac_ntv = g1fac_ntv / (nt*np) / vprime

  ! compute hhat
  allocate(h(nt,np))
  do kt=1,nt
     do kp=1,np
        h(kt,kp) = bmag(kt,kp) * ahat(kt,kp) &
             - bsqahat_avg/bsq_avg * bmag(kt,kp)
     enddo
  enddo
  vec(:) = 0.0
  do i=1,matsize
     call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'d0')
     do kt=1,nt
        do kp=1,np
           vec(i) = vec(i) + bk(kt,kp) * h(kt,kp)
        enddo
     enddo
  enddo
  vec = vec / (nt*np)
  vec(indx_c00) = 0.0
  call DGETRS('N',matsize,1,bgradmat(:,:),matsize,i_piv,vec(:),matsize,info)
  allocate(hhat(nt,np))
  hhat(:,:) = 0.0
  do i=1,matsize
     call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'d0')
     do kt=1,nt
        do kp=1,np
           hhat(kt,kp) = hhat(kt,kp) + bk(kt,kp) * vec(i)
        enddo
     enddo
  enddo

  g3fac=0.0
  do i=1,nt
     do j=1,np
        g3fac = g3fac + g(i,j) * bmag(i,j) * a(i,j) * hhat(i,j)
     enddo
  enddo
  g3fac = g3fac / (nt*np) / vprime

  ! Parameters for Simikov PoP 2009 PS theory
  ! u = 2*ahat

  ! <u B^2>
  gs1 = 2.0*bsqahat_avg

  ! <u^2 B^2>
  gs2 = 0.0
  do i=1,nt
     do j=1,np
        gs2 = gs2  + g(i,j) * bmag(i,j)**2 * ahat(i,j)**2 * 4.0
     enddo
  enddo
  gs2 = gs2 / (nt*np) / vprime

  ! <u ((bdotgradB)/B)^2>
  gs3 = 0.0
  do i=1,nt
     do j=1,np
        gs3 = gs3  + g(i,j) * bdotgradB_overB(i,j)**2 * ahat(i,j) * 2.0
     enddo
  enddo
  gs3 = gs3 / (nt*np) / vprime

  ! bdotgrad(uB^2)
  allocate(ubfac(nt,np))
  ubfac = 2.0 * a * bmag**2 + 4.0 * ahat * bmag * (bdotgradB_overB*bmag)

  ! < ((bdotgradB)/B) bdotgrad(uB^2)>
  gs4 = 0.0
  do i=1,nt
     do j=1,np
        gs4 = gs4  + g(i,j) * ubfac(i,j) * bdotgradB_overB(i,j)
     enddo
  enddo
  gs4 = gs4 / (nt*np) / vprime

  ! < (bdotgrad(uB^2) / B)^2 >
  gs5 = 0.0
  do i=1,nt
     do j=1,np
        gs5 = gs5  + g(i,j) * ubfac(i,j)**2 / bmag(i,j)**2
     enddo
  enddo
  gs5 = gs5 / (nt*np) / vprime

  ! < u ((bdotgradB)/B) bdotgrad(uB^2)>
  gs6 = 0.0
  do i=1,nt
     do j=1,np
        gs6 = gs6  + g(i,j) * ubfac(i,j) * bdotgradB_overB(i,j) &
             * ahat(i,j) * 2.0
     enddo
  enddo
  gs6 = gs6 / (nt*np) / vprime

  open(unit=1,file='out.le3.theory',status='replace')
  write (1,'(e16.8)') g1fac
  write (1,'(e16.8)') g2fac
  write (1,'(e16.8)') g3fac
  write (1,'(e16.8)') g1fac_ntv
  write (1,'(e16.8)') bsq_avg
  write (1,'(e16.8)') bgradbsq_avg
  write (1,'(e16.8)') gs1
  write (1,'(e16.8)') gs2
  write (1,'(e16.8)') gs3
  write (1,'(e16.8)') gs4
  write (1,'(e16.8)') gs5
  write (1,'(e16.8)') gs6
  close(1)

  deallocate(bgradmat)
  deallocate(vec)
  deallocate(a)
  deallocate(ahat)
  deallocate(u)
  deallocate(uhat)
  deallocate(h)
  deallocate(hhat)
  deallocate(gfac)
  deallocate(i_piv)
  deallocate(ubfac)

end subroutine le3_compute_theory
