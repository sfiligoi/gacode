!---------------------------------------------------------
! cgyro_advect_wavenumber.f90
!
! PURPOSE:
!  Manage shearing by wavenumber advection.
!---------------------------------------------------------

subroutine cgyro_advect_wavenumber(ij)

  use cgyro_globals
  use timer_lib

  implicit none

  integer, intent(in) :: ij
  integer :: ir,l,ll,j,icc,in
  integer :: irp,irm
  complex, dimension(:,:),allocatable :: h0
  complex, dimension(:,:,:),allocatable :: he
  complex, dimension(nv_loc) :: dh
  real :: scale


  call timer_lib_in('shear')
  allocate(h0(nv_loc,1-2*n_wave:n_radial+2*n_wave))
  allocate(he(n_theta,nv_loc,1-2*n_wave:n_radial+2*n_wave))

  ! Wavenumber advection ExB shear

  if (shear_method == 2) then
    !ic_c(ir,j) = j + (ir-1)*n_theta

!$omp parallel 

!$omp do private(j,ir,in,icc)
    do ir=1,n_radial
      !icc =ic_c(ir,1)-1
      icc = (ir-1)*n_theta
      do in=1,nv_loc
        do j=1,n_theta
           he(j,in,ir) = omega_eb*h_x(icc+j,in)
        enddo
      enddo
    enddo
!$omp end do nowait

   ! Zero wavenumbers outside n_radial
!$omp do private(l,ll)
    do l=n_wave,1,-1
      ll = 1-2*l
      he(:,:,ll:ll+1) = 0.0
    enddo
!$omp end do nowait

!$omp do private(l,ll)
    do l=1,n_wave
      ll = n_radial+2*l
      he(:,:,ll-1:ll) = 0.0
    enddo
!$omp end do
  ! here is an implicit barrier

!$omp do private(j,ir,in,icc,ll,l)
     do ir=1,n_radial
        !icc =ic_c(ir,1)-1
        icc = (ir-1)*n_theta
        do in=1,nv_loc
           do l=1,n_wave
              ll = 2*l-1
              do j=1,n_theta
                 rhs(icc+j,in,ij) = rhs(icc+j,in,ij)+c_wave(l)*(he(j,in,ir+ll)-he(j,in,ir-ll))
              enddo
           enddo
        enddo
     enddo
!$omp end do

!$omp end parallel 
  endif

  ! Wavenumber advection profile shear

  if (profile_shear_flag == 1) then
    !ic_c(ir,j) = j + (ir-1)*n_theta

!$omp parallel 

!$omp do private(j,ir,in,icc)
    do ir=1,n_radial
      !icc =ic_c(ir,1)-1
      icc = (ir-1)*n_theta
      do in=1,nv_loc
        do j=1,n_theta
           he(j,in,ir) = sum(omega_ss(:,icc+j,in)*field(:,icc+j))
        enddo
      enddo
    enddo
!$omp end do nowait

   ! Zero wavenumbers outside n_radial
!$omp do private(l,ll)
    do l=n_wave,1,-1
      ll = 1-2*l
      he(:,:,ll:ll+1) = 0.0
    enddo
!$omp end do nowait

!$omp do private(l,ll)
    do l=1,n_wave
      ll = n_radial+2*l
      he(:,:,ll-1:ll) = 0.0
    enddo
!$omp end do
  ! here is an implicit barrier

!$omp do private(j,ir,in,icc,ll,l)
     do ir=1,n_radial
        !icc =ic_c(ir,1)-1
        icc = (ir-1)*n_theta
        do in=1,nv_loc
           do l=1,n_wave
              ll = 2*l-1
              do j=1,n_theta
                 rhs(icc+j,in,ij) = rhs(icc+j,in,ij)+c_wave(l)*(he(j,in,ir+ll)-he(j,in,ir-ll))
              enddo
           enddo
        enddo
     enddo
!$omp end do

!$omp end parallel 
  endif

  ! Zonal damping (to reduce box-size correlation)

  if (profile_shear_flag == 1 .or. shear_method == 2) then
     if (n == 0) then
        irp = +1+1+n_radial/2
        irm = -1+1+n_radial/2
        scale = nu_global*abs(gamma_e)/2.0
!$omp parallel do private(j,h0,dh,ic)
        do j=1,n_theta

           ! p = +1
           h0(:,irp) = h_x(ic_c(irp,j),:)

           ! p = -1
           h0(:,irm) = h_x(ic_c(irm,j),:)

           ! p = +1
           dh(:) = -scale*(h0(:,irp)-h0(:,irm))
           ic = ic_c(irp,j)
           rhs(ic,:,ij) = rhs(ic,:,ij)+dh(:)

           ! p = -1
           dh(:) = -scale*(h0(:,irm)-h0(:,irp))
           ic = ic_c(irm,j)
           rhs(ic,:,ij) = rhs(ic,:,ij)+dh(:)

        enddo
     endif
  endif

  deallocate(he)
  deallocate(h0)
  call timer_lib_out('shear')

end subroutine cgyro_advect_wavenumber
