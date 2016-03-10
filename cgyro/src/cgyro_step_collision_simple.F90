!------------------------------------------------------------------------------
! cgyro_step_collision_simple.F90
!
! PURPOSE:
!  Take an implicit collision step using the pre-computed collision 
!  matrix.  Effectively, we compute the new collisional cap_H: 
!
!                       H = h + ze/T G phi
!------------------------------------------------------------------------------

subroutine cgyro_step_collision_simple

  use parallel_lib
  use timer_lib

  use cgyro_globals

  implicit none

  integer :: is,ie,ix,jx,it,ir
  integer :: ivp
  complex, dimension(n_species,n_energy,n_xi) :: cvec
  complex, dimension(n_species,n_energy,n_xi) :: bvec
  real :: cvec_re,cvec_im

  !----------------------------------------------------------------
  ! Perform data tranpose from _c to _v data layouts:
  call timer_lib_in('coll_comm')
  call parallel_lib_rtrans(cap_h_c,cap_h_v)
  call timer_lib_out('coll_comm')
  !----------------------------------------------------------------

  call timer_lib_in('coll')

  do ic=nc1,nc2
     ic_loc = ic-nc1+1
     ir = ir_c(ic)
     it = it_c(ic)

     ! Set-up the RHS: H = f + ze/T G phi

     do iv=1,nv
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        cvec(is,ie,ix) = cap_h_v(ic_loc,iv)
     enddo

     ! Avoid singularity of n=0,p=0:
     if (px(ir) == 0 .and. n == 0) then
        bvec = cvec
     else

        do is=1,n_species
           do ie=1,n_energy
              bvec(is,ie,:) = (0.0,0.0)
              
              do jx=1,n_xi
                 
                 cvec_re = real(cvec(is,ie,jx),kind=kind(cmat_simple))
                 cvec_im = aimag(cvec(is,ie,jx))
                 
                 do ix=1,n_xi
                    bvec(is,ie,ix) = bvec(is,ie,ix)+ &
                         cmplx(cmat_simple(ix,jx,is,ie,it)*cvec_re, &
                         cmat_simple(ix,jx,is,ie,it)*cvec_im)
                 enddo
              enddo
           enddo
        enddo
     endif

     do iv=1,nv
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
       cap_h_v(ic_loc,iv) = bvec(is,ie,ix)
     enddo

  enddo

  call timer_lib_out('coll')

  call timer_lib_in('coll_comm')
  call parallel_lib_f(cap_h_v,cap_h_ct)
  cap_h_c = transpose(cap_h_ct)
  call timer_lib_out('coll_comm')

  call timer_lib_in('coll')

  ! Compute H given h and [phi(h), apar(h)]

  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)
     do ic=1,nc
        psi(ic,iv_loc) = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        h_x(ic,iv_loc) = cap_h_c(ic,iv_loc)-psi(ic,iv_loc)*(z(is)/temp(is))
     enddo
  enddo

  call timer_lib_out('coll')

end subroutine cgyro_step_collision_simple
