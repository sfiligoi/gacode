subroutine cgyro_init_implicit_gk

  use mpi
  use parallel_lib
  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer :: is, ir, it, ie, ix
  integer :: id, jd, idp, jdp, jt, jr, jc, ifield, jfield
  real :: rval
  real :: rfac, vfac, vfac2
  complex, dimension(:,:), allocatable   :: akmat, fieldmat_loc
  complex, dimension(:), allocatable     :: work_field

  if(implicit_flag == 0) return

  ! Kinetic eqn solve matrix(ic,ic,nv_loc) 
  ! gkmat = (1 + delta_t/2 * stream)
  ! EAB: TO DO: MAKE THIS SPARSE

  allocate(gkmat(nc,nc,nv_loc))
  allocate(i_piv_gk(nc,nv_loc))
  gkmat(:,:,:) = (0.0,0.0)

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
        
        do id=-2,2
           jt = thcyc(it+id)
           jr = rcyc(ir,it,id)
           jc = ic_c(jr,jt)

           ! EAB: TO DO: MAKE THIS SPARSE
           gkmat(ic,jc,iv_loc) = gkmat(ic,jc,iv_loc) &
                + (rval*dtheta(ir,it,id) &
                + abs(rval)*dtheta_up(ir,it,id)) * 0.5 * delta_t
           
        enddo

        gkmat(ic,ic,iv_loc) = gkmat(ic,ic,iv_loc) + 1.0

     enddo

     ! Lapack factorization 
     ! EAB: TO DO: MAKE THIS SPARSE
     call ZGETRF(nc,nc,gkmat(:,:,iv_loc),nc,i_piv_gk(:,iv_loc),info)
  
  enddo

  ! Field eqn matrix (ic*ifield,ic*ifield) (dense)
  allocate(idfield(nc,n_field))
  allocate(fieldmat_loc(nc*n_field,nc*n_field))
  allocate(akmat(nc,nc))
  allocate(i_piv_field(nc))
  allocate(work_field(nc))

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

     rfac  = 0.5*w_xi(ix)*w_e(ie)*dens(is) * z(is)**2/temp(is)
     vfac  = xi(ix) * sqrt(2.0*energy(ie)) * vth(is)
     vfac2 = xi(ix)**2 * 2.0*energy(ie) * vth(is)**2

     akmat(:,:) = (0.0,0.0)
     do ic=1,nc
        
        ir = ir_c(ic) 
        it = it_c(ic)

        rval = omega_stream(it,is)*sqrt(energy(ie))*xi(ix) 
        
        do id=-2,2
           jt = thcyc(it+id)
           jr = rcyc(ir,it,id)
           jc = ic_c(jr,jt)

           akmat(ic,jc) = akmat(ic,jc) &
                + (rval*dtheta(ir,it,id) &
                + abs(rval)*dtheta_up(ir,it,id)) * 0.5 * delta_t
           
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
                - j0_c(ic,iv_loc) * akmat(ic,jc) * j0_c(jc,iv_loc) &
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
                   + j0_c(ic,iv_loc) * akmat(ic,jc) *j0_c(jc,iv_loc) &
                   * rfac * vfac
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
                   + j0_c(ic,iv_loc) * akmat(ic,jc) *j0_c(jc,iv_loc) &
                   * rfac * vfac2
           enddo
        enddo

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
          + (k_perp(it,ir)**2*lambda_debye**2* &
          dens_ele/temp_ele+sum_den_h)
  enddo

  if(n_field > 1) then

     ifield=2
     jfield=2
     do ic=1,nc
        id = idfield(ic,ifield)
        ir = ir_c(ic) 
        it = it_c(ic)
        fieldmat(id,id) = fieldmat(id,id) &
             + (2.0*k_perp(it,ir)**2*rho**2 &
             /betae_unit*dens_ele*temp_ele)
     enddo

  endif

  ! Lapack factorization
  allocate(i_piv_field(nc*n_field))
  call ZGETRF(nc*n_field,nc*n_field,fieldmat(:,:),nc*n_field,i_piv_field,info)

  ! Extra allocations needed for solve
  allocate(gkvec(nc,nv_loc))
  allocate(fieldvec(nc*n_field))
  allocate(fieldvec_loc(nc*n_field))

end subroutine cgyro_init_implicit_gk

subroutine cgyro_clean_implicit_gk
  use cgyro_globals

  implicit none

  if(implicit_flag == 0) return

  if(allocated(gkmat))        deallocate(gkmat)
  if(allocated(i_piv_gk))     deallocate(i_piv_gk)
  if(allocated(gkvec))        deallocate(gkvec)
  if(allocated(fieldmat))     deallocate(fieldmat)
  if(allocated(idfield))      deallocate(idfield)
  if(allocated(i_piv_field))  deallocate(i_piv_field)
  if(allocated(fieldvec))     deallocate(fieldvec)
  if(allocated(fieldvec_loc)) deallocate(fieldvec_loc)

end subroutine cgyro_clean_implicit_gk

subroutine cgyro_step_implicit_gk

  use mpi
  use timer_lib
  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer :: is, ir, it, ie, ix
  integer :: id, jt, jr, jc, ifield
  complex :: efac(n_field)
  real    :: rval, rfac, vfac

  if(implicit_flag == 0) return

  call timer_lib_in('rhs_imp')

  ! Solve the gk eqn for the part of RHS depending on old H,fields
  ! RHS = (1 - delta_t/2 * stream)*H_old - (Z f0/T)G field_old

  ! form the rhs

  gkvec(:,:) = cap_h_c(:,:)
  iv_loc = 0
  do iv=nv1,nv2
     iv_loc = iv_loc+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)
  
     efac(1) = 1.0
     if (n_field > 1) then
        efac(2) = -xi(ix)*sqrt(2.0*energy(ie))*vth(is)
     endif

     do ic=1,nc
        ir = ir_c(ic)
        it = it_c(ic)

        rval = omega_stream(it,is)*sqrt(energy(ie))*xi(ix) 
        
        do id=-2,2
           jt = thcyc(it+id)
           jr = rcyc(ir,it,id)
           jc = ic_c(jr,jt)

           gkvec(ic,iv_loc) = gkvec(ic,iv_loc) &
                + (rval*dtheta(ir,it,id) &
                + abs(rval)*dtheta_up(ir,it,id)) &
                * (-0.5 * delta_t) * cap_h_c(jc,iv_loc)
           
        enddo

        gkvec(ic,iv_loc) = gkvec(ic,iv_loc) &
             - z(is)/temp(is)*j0_c(ic,iv_loc)*sum(efac(:)*field(ir,it,:))

     enddo

     ! matrix solve
     ! EAB: TO DO: MAKE THIS SPARSE
     call ZGETRS('N',nc,1,gkmat(:,:,iv_loc),nc,i_piv_gk(:,iv_loc),&
          gkvec(:,iv_loc),nc,info)

  enddo

  ! Field solve

  ! form the rhs

  fieldvec_loc(:) = (0.0,0.0)
  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1
     
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     rfac = z(is) * 0.5*w_xi(ix)*w_e(ie)*dens(is)
     vfac = xi(ix) * sqrt(2.0*energy(ie)) * vth(is)

     ifield = 1
     do ic=1,nc
        id = idfield(ic,ifield)
        ir = ir_c(ic) 
        it = it_c(ic)
        fieldvec_loc(id) = fieldvec_loc(id) + gkvec(ic,iv_loc) &
             * rfac * j0_c(ic,iv_loc)
     enddo

     if(n_field > 1) then
        ifield = 2
        do ic=1,nc
           id = idfield(ic,ifield)
           ir = ir_c(ic) 
           it = it_c(ic)
           fieldvec_loc(id) = fieldvec_loc(id) + gkvec(ic,iv_loc) &
                * rfac * j0_c(ic,iv_loc) * vfac
        enddo
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

  ! map back into field(ir,it)
  do ic=1,nc
     do ifield=1,n_field

        id = idfield(ic,ifield)
        ir = ir_c(ic) 
        it = it_c(ic)

        field(ir,it,ifield) = fieldvec(id)
        
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
  
     efac(1) = 1.0
     if (n_field > 1) then
        efac(2) = -xi(ix)*sqrt(2.0*energy(ie))*vth(is)
     endif

     do ic=1,nc
        ir = ir_c(ic)
        it = it_c(ic)

        gkvec(ic,iv_loc) = gkvec(ic,iv_loc) &
             + z(is)/temp(is)*j0_c(ic,iv_loc)*sum(efac(:)*field(ir,it,:))

     enddo

     ! matrix solve
     call ZGETRS('N',nc,1,gkmat(:,:,iv_loc),nc,i_piv_gk(:,iv_loc),&
          gkvec(:,iv_loc),nc,info)

  enddo

  cap_h_c(:,:) = cap_h_c(:,:) + gkvec(:,:)

  ! Reform little h and psi

  iv_loc = 0
  do iv=nv1,nv2
     
     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     efac(1) = 1.0
     if (n_field > 1) then
        efac(2) = -xi(ix)*sqrt(2.0*energy(ie))*vth(is)
     endif
     do ic=1,nc

        ir = ir_c(ic)
        it = it_c(ic)

        psi(ic,iv_loc) = j0_c(ic,iv_loc)*sum(efac(:)*field(ir,it,:))

        h_x(ic,iv_loc) = cap_h_c(ic,iv_loc) - z(is)*psi(ic,iv_loc)/temp(is)

     enddo
  enddo

  call timer_lib_out('rhs_imp')

end subroutine cgyro_step_implicit_gk
