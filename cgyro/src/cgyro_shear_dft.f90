!---------------------------------------------------------
! cgyro_shear_dft.f90
!
! PURPOSE:
!  DFT shear algorithm.  DFT coefficients are shifted
!  continuously to the left or right depending upon 
!  sign of omega_eb.
!
! NOTE:
!                       k_theta*length*gamma_e
!           omega_eb  = ---------------------- 
!                               2 pi
!---------------------------------------------------------

subroutine cgyro_shear_dft

  use cgyro_globals

  implicit none

  integer :: ix,p,nxs
  integer :: ir,it
  real :: a

  include 'fftw3.f03'

  a = omega_eb*delta_t

  ! Partial-shift option (dshift > 0)
  if (dshift > 0.0) then
     gtime = gtime+a
     if (gtime > dshift) then
        a     = dshift
        gtime = gtime-a
     else if (gtime < -dshift) then
        a     = -dshift 
        gtime = gtime-a
     else
        return
     endif
  endif

  nxs = n_radial+2*shear_pad

  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     do it=1,n_theta

        fp(:) = 0.0
        do ir=1,n_radial
           ix = ir-1-n_radial/2
           if (ix < 0) ix = ix+nxs
           fp(ix) = h_x(ic_c(ir,it),iv_loc)
        enddo

        !do p=-nxs/2,nxs/2-1
        !   ix = p
        !   if (ix < 0) ix=ix+nxs
        !   print '(4(1pe12.5,1x))',fp(ix)
        !enddo

        call fftw_execute_dft(plan_p2j,fp,fj)
        do ix=0,nxs-1
           p = ix
           if (p > nxs/2-1) p=p-nxs
           fj(ix) = exp(2*pi*i_c*p*a/nxs)*fj(ix)
        enddo
        call fftw_execute_dft(plan_j2p,fj,fp)

        do ir=1,n_radial
           ix = ir-1-n_radial/2
           if (ix < 0) ix = ix+nxs
           h_x(ic_c(ir,it),iv_loc) = fp(ix)/nxs
        enddo

        !print *
        !do p=-nxs/2,nxs/2-1
        !   ix = p
        !   if (ix < 0) ix=ix+nxs
        !   print '(4(1pe12.5,1x))',fp(ix)/nxs
        !enddo
        !stop

     enddo
  enddo

  call cgyro_field_c

end subroutine cgyro_shear_dft
