subroutine cgyro_init_implicit_gk

  use mpi
  use parallel_lib
  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: is, ir, it, ie, ix, k
  integer :: id, jd, idp, jdp, jt, jr, jc, ifield, jfield
  real :: rval
  real :: rfac
  complex, dimension(:,:), allocatable :: akmat, fieldmat_loc
  complex, dimension(:), allocatable   :: work_field

  if (implicit_flag /= 1) return

  if (zf_test_flag == 1 .and. ae_flag == 1) then
     call cgyro_error('ERROR: (CGYRO) Implicit ZF test with adiab. electrons not implemented')
     return
  endif

  ! Kinetic eqn solve matrix(ic,ic,nv_loc) 
  ! gksp_mat = (1 + delta_t/2 * stream)

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

  gksp_nelem=(2*nup_theta+1)*nc
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

        rval = omega_stream(it,is)*vel(ie)*xi(ix) 

        do id=-nup_theta,nup_theta
           jc = icd_c(ic,id)
           k=k+1
           gksp_mat(k,iv_loc) = (rval*dtheta(ic,id) &
                + abs(rval)*dtheta_up(ic,id)) * 0.5 * delta_t
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

  ! Field eqn matrix (ic*ifield,ic*ifield) (dense)
  allocate(idfield(nc,n_field))
  allocate(fieldmat_loc(nc*n_field,nc*n_field))
  allocate(akmat(nc,nc))
  allocate(i_piv_field(nc))
  allocate(work_field(nc))
  fieldmat_loc(:,:) = (0.0,0.0)

  id=0
  do ifield=1,n_field
     do ic=1,nc
        id = id + 1
        idfield(ic,ifield) = id
     enddo
  enddo

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     rfac = w_xi(ix)*w_e(ie)*dens(is) * z(is)**2/temp(is)

     akmat(:,:) = (0.0,0.0)
     do ic=1,nc

        ir = ir_c(ic) 
        it = it_c(ic)

        rval = omega_stream(it,is)*vel(ie)*xi(ix) 

        do id=-nup_theta,nup_theta
           jc = icd_c(ic,id)

           akmat(ic,jc) = akmat(ic,jc) &
                + (rval*dtheta(ic,id) &
                + abs(rval)*dtheta_up(ic,id)) * 0.5 * delta_t

        enddo

        akmat(ic,ic) = akmat(ic,ic) + 1.0

     enddo

     ! Lapack factorization and inverse
     ! akmat = (1 + delta_t/2 * stream)^-1
     call ZGETRF(nc,nc,akmat(:,:),nc,i_piv_field,info)
     call ZGETRI(nc,akmat(:,:),nc,i_piv_field,work_field,nc,info)

     ! int d^3 v (z G) akmat (z/T f0 G)
     ifield=1
     jfield=1
     do ic=1,nc
        id = idfield(ic,ifield)
        do jc=1,nc
           jd = idfield(jc,jfield)
           fieldmat_loc(id,jd) = fieldmat_loc(id,jd) &
                - jvec_c(1,ic,iv_loc) * akmat(ic,jc) * jvec_c(1,jc,iv_loc) &
                * rfac
        enddo
     enddo

     if(n_field > 1) then

        ! int d^3 v vpar (z G) akmat (z/T f0 G)
        ifield=1
        jfield=2
        do ic=1,nc
           id = idfield(ic,ifield)
           do jc=1,nc
              jd = idfield(jc,jfield)
              fieldmat_loc(id,jd) = fieldmat_loc(id,jd) &
                   - jvec_c(1,ic,iv_loc) * akmat(ic,jc) *jvec_c(2,jc,iv_loc) &
                   * rfac
           enddo
        enddo

        ! int d^3 v vpar^2 (z G) akmat (z/T f0 G)
        ifield=2
        jfield=2
        do ic=1,nc
           id = idfield(ic,ifield)
           do jc=1,nc
              jd = idfield(jc,jfield)
              fieldmat_loc(id,jd) = fieldmat_loc(id,jd) &
                   + jvec_c(2,ic,iv_loc) * akmat(ic,jc) *jvec_c(2,jc,iv_loc) &
                   * rfac
           enddo
        enddo

        if(n_field > 2) then

           ! int d^3 v (z G) akmat (vperp/vt)^2 Bunit/B G1
           ifield=1
           jfield=3
           do ic=1,nc
              id = idfield(ic,ifield)
              do jc=1,nc
                 jt = it_c(jc)
                 jd = idfield(jc,jfield)
                 fieldmat_loc(id,jd) = fieldmat_loc(id,jd) &
                      - jvec_c(1,ic,iv_loc) * akmat(ic,jc) * jvec_c(3,jc,iv_loc) &
                      * rfac 
              enddo
           enddo

           ! int d^3 v (z G) vpar akmat (vperp/vt)^2 Bunit/B G1
           ifield=2
           jfield=3
           do ic=1,nc
              id = idfield(ic,ifield)
              do jc=1,nc
                 jt = it_c(jc)
                 jd = idfield(jc,jfield)
                 fieldmat_loc(id,jd) = fieldmat_loc(id,jd) &
                      + jvec_c(2,ic,iv_loc) * akmat(ic,jc) * jvec_c(3,jc,iv_loc) &
                      * rfac
              enddo
           enddo

           ! (bfac) int d^3 v (z G1 (vperp/vt)^2 Bunit/B) akmat G0
           ifield=3
           jfield=1
           do ic=1,nc
              it = it_c(ic)
              id = idfield(ic,ifield)
              do jc=1,nc
                 jt = it_c(jc)
                 jd = idfield(jc,jfield)
                 fieldmat_loc(id,jd) = fieldmat_loc(id,jd) &
                      - (-0.5*betae_unit)/(dens_ele*temp_ele) &
                      * jvec_c(3,ic,iv_loc) * akmat(ic,jc) * jvec_c(1,jc,iv_loc) &
                      * rfac 
              enddo
           enddo

           ! (bfac) int d^3 v (z G1 (vperp/vt)^2 Bunit/B) akmat G0 vpar
           ifield=3
           jfield=2
           do ic=1,nc
              it = it_c(ic)
              id = idfield(ic,ifield)
              do jc=1,nc
                 jt = it_c(jc)
                 jd = idfield(jc,jfield)
                 fieldmat_loc(id,jd) = fieldmat_loc(id,jd) &
                      - (-0.5*betae_unit)/(dens_ele*temp_ele) &
                      * jvec_c(3,ic,iv_loc) * akmat(ic,jc) * jvec_c(2,jc,iv_loc) &
                      * rfac 
              enddo
           enddo

           ! (bfac) int d^3 v (T G1 (vperp/vt)^2 Bunit/B) akmat 
           ! G1 (vperp/vt)^2 Bunit/B
           ifield=3
           jfield=3
           do ic=1,nc
              it = it_c(ic)
              id = idfield(ic,ifield)
              do jc=1,nc
                 jt = it_c(jc)
                 jd = idfield(jc,jfield)
                 fieldmat_loc(id,jd) = fieldmat_loc(id,jd) &
                      - (-0.5*betae_unit)/(dens_ele*temp_ele) &
                      * jvec_c(3,ic,iv_loc) * akmat(ic,jc) * jvec_c(3,jc,iv_loc) &
                      * rfac 
              enddo
           enddo

        endif

     endif

  enddo

  if(n_field > 1) then
     do ic=1,nc
        id  = idfield(ic,2)
        idp = idfield(ic,1)
        do jc=1,nc
           jd  = idfield(jc,1)
           jdp = idfield(jc,2)
           fieldmat_loc(id,jd) = -fieldmat_loc(idp,jdp)
        enddo
     enddo
  endif

  ! Special case for n=0,p=0
  if(n == 0) then
     do ic=1,nc
        ir = ir_c(ic)
        it = it_c(ic)
        if(px(ir) == 0) then

           do jc=1,nc
              do ifield=1,n_field
                 do jfield=1,n_field
                    id = idfield(ic,ifield)
                    jd = idfield(jc,jfield)
                    fieldmat_loc(id,jd) = 0.0
                 enddo
              enddo
           enddo

           do jc=1,nc
              jr = ir_c(jc)
              jt = it_c(jc)
              if(jt==it .and. px(jr)==0) then
                 id = idfield(ic,1)
                 jd = idfield(jc,1)
                 fieldmat_loc(id,jd) = 1.0
                 if(n_field > 1) then
                    id = idfield(ic,2)
                    jd = idfield(jc,2)
                    fieldmat_loc(id,jd) = 1.0
                    if(n_field > 2) then
                       id = idfield(ic,3)
                       jd = idfield(jc,3)
                       fieldmat_loc(id,jd) = 1.0
                    endif
                 endif
              endif
           enddo

        endif
     enddo
  endif

  deallocate(akmat)
  deallocate(i_piv_field)
  deallocate(work_field)

  allocate(fieldmat(nc*n_field,nc*n_field))

  ! Velocity space reduce

  call MPI_ALLREDUCE(fieldmat_loc(:,:),&
       fieldmat(:,:),&
       size(fieldmat(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  deallocate(fieldmat_loc)

  ! Add Poisson/Ampere factors 

  ifield=1
  do ic=1,nc
     id = idfield(ic,ifield)
     ir = ir_c(ic) 
     it = it_c(ic)
     fieldmat(id,id) = fieldmat(id,id) &
          + (k_perp(ic)**2*lambda_debye**2* &
          dens_ele/temp_ele+sum_den_h(it))
  enddo

  if(n_field > 1) then

     ifield=2
     jfield=2
     do ic=1,nc
        id = idfield(ic,ifield)
        ir = ir_c(ic) 
        it = it_c(ic)
        fieldmat(id,id) = fieldmat(id,id) &
             + (2.0*k_perp(ic)**2*rho**2 &
             /betae_unit*dens_ele*temp_ele)
     enddo

     if(n_field > 2) then
        ifield=3
        jfield=3
        do ic=1,nc
           id = idfield(ic,ifield)
           ir = ir_c(ic) 
           it = it_c(ic)
           fieldmat(id,id) = fieldmat(id,id) &
                + 1.0
        enddo
     endif

  endif

  ! Lapack factorization
  allocate(i_piv_field(nc*n_field))
  call ZGETRF(nc*n_field,nc*n_field,fieldmat(:,:),nc*n_field,i_piv_field,info)

  ! Extra allocations needed for solve
  allocate(gkvec(nc,nv_loc))
  allocate(gksvec(nc))
  allocate(gkwvec(2*nc))
  allocate(fieldvec(nc*n_field))
  allocate(fieldvec_loc(nc*n_field))

end subroutine cgyro_init_implicit_gk

subroutine cgyro_clean_implicit_gk

  use cgyro_globals

  implicit none

  if (implicit_flag /= 1) return

  if(allocated(gkvec))        deallocate(gkvec)
  if(allocated(fieldmat))     deallocate(fieldmat)
  if(allocated(idfield))      deallocate(idfield)
  if(allocated(i_piv_field))  deallocate(i_piv_field)
  if(allocated(fieldvec))     deallocate(fieldvec)
  if(allocated(fieldvec_loc)) deallocate(fieldvec_loc)

  if(allocated(gksp_mat))     deallocate(gksp_mat)
  if(allocated(gksp_indx))    deallocate(gksp_indx)
  if(allocated(gksp_cntl))    deallocate(gksp_cntl)
  if(allocated(gksp_icntl))   deallocate(gksp_icntl)
  if(allocated(gksp_keep))    deallocate(gksp_keep)
  if(allocated(gksvec))       deallocate(gksvec)
  if(allocated(gkwvec))       deallocate(gkwvec)

end subroutine cgyro_clean_implicit_gk

subroutine cgyro_step_implicit_gk

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none

  integer :: is, ir, it, ie, ix
  integer :: id, jc, ifield
  real    :: rval
  complex :: rfac(nc)

  if (implicit_flag /= 1) return

  call timer_lib_in('stream')

  ! Solve the gk eqn for the part of RHS depending on old H,fields
  ! RHS = (1 - delta_t/2 * stream)*H_old - (Z f0/T)G field_old

  ! form the rhs

  ! Special case for n=0, p=0
  if (n == 0) then
     do ic=1,nc
        ir = ir_c(ic)
        it = it_c(ic)
        if (px(ir) == 0) then
           cap_h_c(ic,:) = 0.0
           field(:,ic)  = 0.0
        endif
     enddo
  endif

  gkvec(:,:) = cap_h_c(:,:)
  iv_loc = 0
  do iv=nv1,nv2
     iv_loc = iv_loc+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc
        ir = ir_c(ic)
        it = it_c(ic)

        rval = omega_stream(it,is)*vel(ie)*xi(ix) 

        do id=-nup_theta,nup_theta
           jc = icd_c(ic,id)
           gkvec(ic,iv_loc) = gkvec(ic,iv_loc) &
                + (rval*dtheta(ic,id) &
                + abs(rval)*dtheta_up(ic,id)) &
                * (-0.5 * delta_t) * cap_h_c(jc,iv_loc)
        enddo

        gkvec(ic,iv_loc) = gkvec(ic,iv_loc) &
             - z(is)/temp(is)*sum(jvec_c(:,ic,iv_loc)*field(:,ic))

     enddo

     ! matrix solve
     call UMZ2SO(nc,0,.false.,gksp_nmax,2*gksp_nmax,&
          gksp_mat(:,iv_loc),gksp_indx(:,iv_loc),gksp_keep(:,iv_loc),&
          gkvec(:,iv_loc),gksvec,gkwvec,&
          gksp_cntl(:,iv_loc),gksp_icntl(:,iv_loc),&
          gksp_uinfo,gksp_rinfo)
     gkvec(:,iv_loc) = gksvec(:)

  enddo

  ! Field solve -- form the RHS

  fieldvec_loc(:) = (0.0,0.0)
  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     rfac(:) = z(is) * w_xi(ix)*w_e(ie)*dens(is)*gkvec(:,iv_loc)
     ifield = 1

     do ic=1,nc
        id = idfield(ic,ifield)
        fieldvec_loc(id) = fieldvec_loc(id) + rfac(ic)*jvec_c(1,ic,iv_loc) 
     enddo

     if(n_field > 1) then
        ifield = 2
        do ic=1,nc
           id = idfield(ic,ifield)
           fieldvec_loc(id) = fieldvec_loc(id)-rfac(ic)*jvec_c(2,ic,iv_loc)
        enddo

        if(n_field > 2) then
           ifield = 3
           do ic=1,nc
              id = idfield(ic,ifield)
              fieldvec_loc(id) = fieldvec_loc(id)+rfac(ic)*jvec_c(3,ic,iv_loc) &
                   *(-0.5*betae_unit)/(dens_ele*temp_ele)           
           enddo
        endif

     endif

  enddo

  call MPI_ALLREDUCE(fieldvec_loc(:),&
       fieldvec(:),&
       size(fieldvec(:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  ! matrix solve
  call ZGETRS('N',nc*n_field,1,fieldmat(:,:),nc*n_field,i_piv_field(:),&
       fieldvec(:),nc*n_field,info)

  ! map back into field(ic)
  do ic=1,nc
     do ifield=1,n_field
        id = idfield(ic,ifield)
        field(ifield,ic) = fieldvec(id)
     enddo
  enddo

  ! Solve the gk eqn for the part of RHS depending on new fields
  ! RHS = (Z f0/T)G field_new

  cap_h_c(:,:) = gkvec(:,:)

  gkvec(:,:) = (0.0,0.0)
  iv_loc = 0
  do iv=nv1,nv2
     iv_loc = iv_loc+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc
        gkvec(ic,iv_loc) = gkvec(ic,iv_loc) &
             +z(is)/temp(is)*sum(jvec_c(:,ic,iv_loc)*field(:,ic))
     enddo

     ! matrix solve
     call UMZ2SO(nc,0,.false.,gksp_nmax,2*gksp_nmax,&
          gksp_mat(:,iv_loc),gksp_indx(:,iv_loc),gksp_keep(:,iv_loc),&
          gkvec(:,iv_loc),gksvec,gkwvec,&
          gksp_cntl(:,iv_loc),gksp_icntl(:,iv_loc),&
          gksp_uinfo,gksp_rinfo)
     gkvec(:,iv_loc) = gksvec(:)

  enddo

  cap_h_c(:,:) = cap_h_c(:,:) + gkvec(:,:)

  ! Reform little h and psi

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc
        psi(ic,iv_loc) = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        h_x(ic,iv_loc) = cap_h_c(ic,iv_loc)-z(is)*psi(ic,iv_loc)/temp(is)
     enddo
  enddo

  call timer_lib_out('stream')

end subroutine cgyro_step_implicit_gk
