subroutine cgyro_init_implicit_gk

  use mpi
  use parallel_lib
  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer :: is, ir, it, ie, ix
  integer :: id, jd, jt, jr, jc, ifield, jfield
  real :: rval
  complex, dimension(:,:), allocatable   :: gkmat_loc, akmat_temp
  complex, dimension(:), allocatable :: zwork
  complex, parameter    :: znum1 = (1.0,0.0)
  complex, parameter    :: znum0 = (0.0,0.0)

  if(implicit_flag == 0) return

  ! Set-up arrays for implicit time-stepping of streaming terms
  allocate(akmat(nc,nc,nv_loc))
  allocate(akmat_temp(nc,nc))
  allocate(gkhmat(nc,nc,nv_loc))
  allocate(gkmat_loc(nc,nc))
  allocate(gk11mat(nc,nc))
  allocate(gk12mat(nc,nc))
  allocate(gk22mat(nc,nc))
  allocate(gkmat(nc*n_field,nc*n_field))
  allocate(i_piv_gk(nc*n_field))
  allocate(gkrhsvec(nc*n_field))
  allocate(gkhvec_loc(nc))
  allocate(gkhvec(nc,n_field))

  ! matrix solve parameters
  allocate(zwork(nc))
  allocate(i_piv(nc))

  akmat(:,:,:)    = (0.0,0.0)
  akmat_temp(:,:) = (0.0,0.0)
  gkhmat(:,:,:)   = (0.0,0.0)
  gkmat_loc(:,:)  = (0.0,0.0)
  gk11mat(:,:)    = (0.0,0.0)
  gk12mat(:,:)    = (0.0,0.0)
  gk22mat(:,:)    = (0.0,0.0)
  gkmat(:,:)      = (0.0,0.0)

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

           akmat(ic,jc,iv_loc) = akmat(ic,jc,iv_loc) &
                + (rval*dtheta(ir,it,id) &
                + abs(rval)*dtheta_up(ir,it,id)) * 0.5 * delta_t
           
        enddo

     enddo

     ! (1 - delta_t/2 * stream)
     akmat_temp(:,:)   = -akmat(:,:,iv_loc) 
     do ic=1,nc
        akmat_temp(ic,ic) = akmat_temp(ic,ic) + 1.0
     enddo

     ! (1 + delta_t/2 * stream)
     do ic=1,nc
        akmat(ic,ic,iv_loc) = akmat(ic,ic,iv_loc) + 1.0
     enddo

     ! Lapack factorization and inverse of LHS
     ! (1 + delta_t/2 * stream)^-1
     call ZGETRF(nc,nc,akmat(:,:,iv_loc),nc,i_piv,info)
     call ZGETRI(nc,akmat(:,:,iv_loc),nc,i_piv,zwork,nc,info)

     ! Matrix multiply
     ! gkhmat = 1 + delta_t/2 * stream)^-1 * (1 - delta_t/2 * stream)
     call ZGEMM('N','N',nc,nc,nc,znum1,akmat(:,:,iv_loc),&
          nc,akmat_temp(:,:),nc,znum0,gkhmat(:,:,iv_loc),nc)
           

     ! akmat = (1 + delta_t/2 * stream)^-1 * (Ze/T) G
     do ic=1,nc
        do jc=1,nc
           akmat(ic,jc,iv_loc) = akmat(ic,jc,iv_loc) &
                * z(is)/temp(is)*j0_c(jc,iv_loc)
        enddo
     enddo

     ! (Ze) G and integration factors 
     do ic=1,nc
        do jc=1,nc
           gkmat_loc(ic,jc) = gkmat_loc(ic,jc) &
                + akmat(ic,jc,iv_loc) * z(is)*j0_c(ic,iv_loc) &
                * 0.5*w_xi(ix)*w_e(ie)*dens(is)
        enddo
     enddo
  enddo

  deallocate(zwork)
  deallocate(i_piv)
        
  ! int d^3v (Ze) G  akmat 
  call MPI_ALLREDUCE(gkmat_loc(:,:),&
       gk11mat(:,:),&
       size(gk11mat(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  ! int d^3v (Ze) G akmat vpar
  gkmat_loc(:,:) = (0.0,0.0)
  iv_loc = 0
  do iv=nv1,nv2
     iv_loc = iv_loc+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)
     do ic=1,nc
        do jc=1,nc
           gkmat_loc(ic,jc) = gkmat_loc(ic,jc) &
                + akmat(ic,jc,iv_loc) * z(is)*j0_c(ic,iv_loc) &
                * 0.5*w_xi(ix)*w_e(ie)*dens(is) &
                * xi(ix) * sqrt(2.0*energy(ie)) * vth(is)
        enddo
     enddo
  enddo
  call MPI_ALLREDUCE(gkmat_loc(:,:),&
       gk12mat(:,:),&
       size(gk12mat(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

   ! int d^3v (Ze) G akmat vpar^2
  gkmat_loc(:,:) = (0.0,0.0)
  iv_loc = 0
  do iv=nv1,nv2
     iv_loc = iv_loc+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)
     do ic=1,nc
        do jc=1,nc
           gkmat_loc(ic,jc) = gkmat_loc(ic,jc) &
                + akmat(ic,jc,iv_loc) * z(is)*j0_c(ic,iv_loc) &
                * 0.5*w_xi(ix)*w_e(ie)*dens(is) &
                * xi(ix)**2 * 2.0*energy(ie) * vth(is)**2
        enddo
     enddo
  enddo
  call MPI_ALLREDUCE(gkmat_loc(:,:),&
       gk22mat(:,:),&
       size(gk22mat(:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  ! form field implicit matrix
  id=1
  do ifield=1,n_field
     do ic=1,nc
        ir = ir_c(ic) 
        it = it_c(ic)
        jd=1
        do jfield=1,n_field
           do jc=1,nc
              
              if(ifield==1 .and. jfield ==1) then
                 gkmat(id,jd) = -gk11mat(ic,jc)
                 if(ic == jc) then
                    gkmat(id,jd) = gkmat(id,jd) &
                         + (k_perp(it,ir)**2*lambda_debye**2* &
                         dens_ele/temp_ele+sum_den_h)
                 endif

              else if(ifield==1 .and. jfield ==2) then
                 gkmat(id,jd) = gk12mat(ic,jc)
              
              
              else if(ifield==2 .and. jfield ==1) then
                 gkmat(id,jd) = -gk12mat(ic,jc)
              
              else if(ifield==2 .and. jfield ==2) then
                 gkmat(id,jd) = gk22mat(ic,jc)
                 if(ic == jc) then
                    gkmat(id,jd) = gkmat(id,jd) &
                         + (2.0*k_perp(it,ir)**2*rho**2 &
                         /betae_unit*dens_ele*temp_ele)
                 endif

              endif
              jd=jd+1
           enddo
        enddo
        id=id+1
     enddo
  enddo
  
  ! Lapack factorization of field LHS
  call ZGETRF(nc*n_field,nc*n_field,gkmat(:,:),nc,i_piv_gk,info)

  ! clean-up
  deallocate(akmat_temp)
  deallocate(gkmat_loc)

end subroutine cgyro_init_implicit_gk

subroutine cgyro_clean_implicit_gk
  use cgyro_globals

  implicit none

  if(implicit_flag == 0) return

  if(allocated(gkhmat))     deallocate(gkhmat)
  if(allocated(akmat))      deallocate(akmat)
  if(allocated(gk11mat))    deallocate(gk11mat)
  if(allocated(gk12mat))    deallocate(gk12mat)
  if(allocated(gk22mat))    deallocate(gk22mat)
  if(allocated(gkmat))      deallocate(gkmat)
  if(allocated(i_piv_gk))   deallocate(i_piv_gk)
  if(allocated(gkrhsvec))   deallocate(gkrhsvec)
  if(allocated(gkhvec))     deallocate(gkhvec)
  if(allocated(gkhvec_loc)) deallocate(gkhvec_loc)

end subroutine cgyro_clean_implicit_gk

subroutine cgyro_step_implicit_gk

  use mpi
  use cgyro_globals

  implicit none

  integer :: is, ir, it, ie, ix
  integer :: id, jt, jr, jc, ifield
  complex :: efac(n_field)

  if(implicit_flag == 0) return

  ! Store the old H and fields
  h0_x(:,:)        = cap_h_c(:,:)
  field_loc(:,:,:) = field(:,:,:)

  ! form the rhs
  gkhvec(:,:) = (0.0,0.0)
  
  ! first integrate the old H
  gkhvec_loc(:) = (0.0,0.0)
  do ic=1,nc
     do jc=1,nc
        iv_loc = 0
        do iv=nv1,nv2
           iv_loc = iv_loc+1
           is = is_v(iv)
           ix = ix_v(iv)
           ie = ie_v(iv)
           gkhvec_loc(ic) = gkhvec_loc(ic) &
                + gkhmat(ic,jc,iv_loc) * z(is)*j0_c(ic,iv_loc) &
                * 0.5*w_xi(ix)*w_e(ie)*dens(is) * cap_h_c(jc,iv_loc)
        enddo
     enddo
        
     call MPI_ALLREDUCE(gkhvec_loc(:),&
          gkhvec(:,1),&
          size(gkhvec(:,1)),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)
  enddo

  ! integrate old H * vpar
  if(n_field > 1) then
     gkhvec_loc(:) = (0.0,0.0)
     do ic=1,nc
        do jc=1,nc
           iv_loc = 0
           do iv=nv1,nv2
              iv_loc = iv_loc+1
              is = is_v(iv)
              ix = ix_v(iv)
              ie = ie_v(iv)
              gkhvec_loc(ic) = gkhvec_loc(ic) &
                   + gkhmat(ic,jc,iv_loc) * cap_h_c(jc,iv_loc) &
                   * 0.5*w_xi(ix)*w_e(ie)*dens(is) &
                   * xi(ix) * sqrt(2.0*energy(ie)) * vth(is)
           enddo
        enddo
        
        call MPI_ALLREDUCE(gkhvec_loc(:),&
             gkhvec(:,2),&
             size(gkhvec(:,2)),&
             MPI_DOUBLE_COMPLEX,&
             MPI_SUM,&
             NEW_COMM_1,&
             i_err)
     enddo
  endif
  

  gkrhsvec(:) = (0.0,0.0)
  
  id=1
  do ifield=1,n_field
     do ic=1,nc
        gkrhsvec(id) = gkrhsvec(id) + gkhvec(ic,ifield)
        id = id + 1
     enddo
  enddo

  id=1
  do ifield=1,n_field
     do ic=1,nc
        ir = ir_c(ic) 
        it = it_c(ic)

        do jc=1,nc
           jr = ir_c(jc) 
           jt = it_c(jc)
           
           if(ifield == 1) then
              gkrhsvec(id) = gkrhsvec(id) - gk11mat(ic,jc)*field(jr,jt,1)
                   
              if(n_field > 1) then
                 gkrhsvec(id) = gkrhsvec(id) + gk12mat(ic,jc)*field(jr,jt,2)
              endif

           else if(ifield ==2) then
              gkrhsvec(id) = gkrhsvec(id) - gk12mat(ic,jc)*field(jr,jt,1)
              if(n_field > 1) then
                 gkrhsvec(id) = gkrhsvec(id) + gk22mat(ic,jc)*field(jr,jt,2)
              endif

           endif

        enddo

        id = id + 1

     enddo
  enddo

  ! Solve the matrix system for the new fields

  call ZGETRS('N',nc*n_field,1,gkmat(:,:),nc*n_field,i_piv_gk,&
       gkrhsvec(:),nc*n_field,info)

  id=1
  do ifield=1,n_field
     do ic=1,nc
        ir = ir_c(ic) 
        it = it_c(ic)
        if(ifield == 1) then
           field(ir,it,1) = gkrhsvec(id)
        else if(ifield == 2) then
           field(ir,it,2) = gkrhsvec(id)
        endif
        id = id + 1
     enddo
  enddo

  ! Now solve for the new capital H
  cap_h_c(:,:) = 0.0
  iv_loc = 0
  do iv=nv1,nv2
     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        do jc=1,nc

           jr = ir_c(jc)
           jt = it_c(jc)

           cap_h_c(ic,iv_loc) = cap_h_c(ic,iv_loc) &
                + gkhmat(ic,jc,iv_loc)*h0_x(jc,iv_loc) &
                + akmat(ic,jc,iv_loc)*(field(jr,jt,1)-field_loc(jr,jt,1))

           if(n_field > 1) then
              cap_h_c(ic,iv_loc) = cap_h_c(ic,iv_loc) &
                   + akmat(ic,jc,iv_loc)*(field(jr,jt,2)-field_loc(jr,jt,2)) &
                   * (-xi(ix)*sqrt(2.0*energy(ie))*vth(is))
           endif

        enddo
     enddo
  enddo


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

end subroutine cgyro_step_implicit_gk
