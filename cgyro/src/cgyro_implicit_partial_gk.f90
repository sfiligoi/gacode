subroutine cgyro_init_implicit_partial_gk

  use mpi
  use parallel_lib
  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: is, ir, it, ie, ix, k
  integer :: id, jt, jr, jc
  real :: rval

  if(implicit_flag /= 2) return

  if(zf_test_flag == 1 .and. ae_flag == 1) then
     call cgyro_error('ERROR: (CGYRO) ZF test with adiabatic electrons not implemented for implicit')
     return
  endif

  ! Kinetic eqn solve matrix(ic,ic,nv_loc) 
  ! gksp_mat = (1 + delta_t/2 * D_stream)

  ! initialization for umfpack
  allocate(gksp_cntl(10,nv_loc))
  allocate(gksp_icntl(20,nv_loc))
  allocate(gksp_keep(20,nv_loc))
  gksp_icntl(1,:) = 6
  gksp_icntl(2,:) = 6
  gksp_icntl(3,:) = 1
  gksp_icntl(8,:) = 0
  iv_loc = 0
  do iv=nv1,nv2
     iv_loc = iv_loc+1
     call UMZ21I(gksp_keep(:,iv_loc),gksp_cntl(:,iv_loc),&
          gksp_icntl(:,iv_loc))
  enddo
  
  gksp_nelem=(2*nup+1)*nc
  gksp_nmax=gksp_nelem*10
  allocate(gksp_mat(gksp_nmax,nv_loc))
  allocate(gksp_indx(2*gksp_nmax,nv_loc))
  gksp_mat(:,:) = (0.0,0.0)
  gksp_indx(:,:) = 0
  
  iv_loc = 0
  do iv=nv1,nv2
     
     iv_loc = iv_loc+1
     
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)
     
     k=0
     do ic=1,nc
        
        ir = ir_c(ic) 
        it = it_c(ic)
        
        rval = omega_stream(it,is)*sqrt(energy(ie))*xi(ix) 
        
        do id=-nup,nup
           jt = thcyc(it+id)
           jr = rcyc(ir,it,id)
           jc = ic_c(jr,jt)
           
           k=k+1
           gksp_mat(k,iv_loc) = abs(rval)*dtheta_up(ir,it,id) * 0.5 * delta_t
           gksp_indx(k,iv_loc) = ic
           gksp_indx(k+gksp_nelem,iv_loc) = jc 
           if(ic == jc) then
              gksp_mat(k,iv_loc) = gksp_mat(k,iv_loc) + 1.0
           endif
        enddo
        
     enddo
     
     ! Umfpack factorization 
     call UMZ2FA(nc,gksp_nelem,0,.false.,gksp_nmax,2*gksp_nmax,&
          gksp_mat(:,iv_loc),gksp_indx(:,iv_loc),&
          gksp_keep(:,iv_loc),gksp_cntl(:,iv_loc),gksp_icntl(:,iv_loc),&
          gksp_uinfo,gksp_rinfo)
     if(gksp_uinfo(1) < 0) then
        call cgyro_error('ERROR: (CGYRO) umfpack error in cgyro_init_implicit_gk')
        return
     endif
     
  enddo

  ! Extra allocations needed for solve
  allocate(gkvec(nc,nv_loc))
  allocate(gksvec(nc))
  allocate(gkwvec(2*nc))

end subroutine cgyro_init_implicit_partial_gk

subroutine cgyro_clean_implicit_partial_gk
  use cgyro_globals

  implicit none

  if(implicit_flag /= 2) return

  if(allocated(gkvec))        deallocate(gkvec)
  if(allocated(gksvec))       deallocate(gksvec)
  if(allocated(gkwvec))       deallocate(gkwvec)
  if(allocated(gksp_mat))     deallocate(gksp_mat)
  if(allocated(gksp_indx))    deallocate(gksp_indx)
  if(allocated(gksp_cntl))    deallocate(gksp_cntl)
  if(allocated(gksp_icntl))   deallocate(gksp_icntl)
  if(allocated(gksp_keep))    deallocate(gksp_keep)

end subroutine cgyro_clean_implicit_partial_gk

subroutine cgyro_step_implicit_partial_gk

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none

  integer :: is, ir, it, ie, ix
  integer :: id, jt, jr, jc
  real    :: rval

  if (implicit_flag /= 2) return
  
  call timer_lib_in('rhs_impgk')

  ! Solve the gk eqn for the part of RHS depending on old H,fields
  ! RHS = (1 - delta_t/2 * D_stream)*h_old 

  ! form the rhs

  ! Special case for n=0, p=0
  if(n == 0) then
     do ic=1,nc
        ir = ir_c(ic)
        it = it_c(ic)
        if(px(ir) == 0) then
           h_x(ic,:) = 0.0
        endif
     enddo
  endif

  gkvec(:,:) = h_x(:,:)
  iv_loc = 0
  do iv=nv1,nv2
     iv_loc = iv_loc+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc
        ir = ir_c(ic)
        it = it_c(ic)

        rval = omega_stream(it,is)*sqrt(energy(ie))*xi(ix) 

        do id=-nup,nup
           jt = thcyc(it+id)
           jr = rcyc(ir,it,id)
           jc = ic_c(jr,jt)

           gkvec(ic,iv_loc) = gkvec(ic,iv_loc) &
                + abs(rval)*dtheta_up(ir,it,id) &
                * (-0.5 * delta_t) * h_x(jc,iv_loc)

        enddo

     enddo

     ! matrix solve
     call UMZ2SO(nc,0,.false.,gksp_nmax,2*gksp_nmax,&
          gksp_mat(:,iv_loc),gksp_indx(:,iv_loc),gksp_keep(:,iv_loc),&
          gkvec(:,iv_loc),gksvec,gkwvec,&
          gksp_cntl(:,iv_loc),gksp_icntl(:,iv_loc),&
          gksp_uinfo,gksp_rinfo)
     gkvec(:,iv_loc) = gksvec(:)

  enddo

  h_x(:,:) = gkvec(:,:) 

  call cgyro_field_c

  call timer_lib_out('rhs_impgk')

end subroutine cgyro_step_implicit_partial_gk
