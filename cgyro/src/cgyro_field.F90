!-----------------------------------------------------------------
! cgyro_field.f90
!
! PURPOSE:
!  Perform field solves for data distributed in both the velocity 
!  and configuration indices:
!  
!  - cgyro_field_v_notae() (configuration distributed)
!  - cgyro_field_c() (velocity distributed)
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Velocity (configuration-distributed) field solve
!-----------------------------------------------------------------
subroutine cgyro_field_v

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none
  
  integer :: itor

  call timer_lib_in('field')

  field_loc(:,:,:) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of H

!$omp parallel do collapse(2) private(ic_loc,iv)
  do itor=nt1,nt2
   do ic=nc1,nc2
     ic_loc = ic-nc1+1
     do iv=1,nv
        field_loc(:,ic,itor) = field_loc(:,ic,itor)+dvjvec_v(:,ic_loc,itor,iv)*cap_h_v(ic_loc,itor,iv)
     enddo
   enddo
  enddo

  call timer_lib_out('field')

  call timer_lib_in('field_com')

  call MPI_ALLREDUCE(field_loc(:,:,:),&
       field(:,:,:),&
       size(field(:,:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

   call timer_lib_out('field_com')
  
  call timer_lib_in('field')

  ! Poisson LHS factors
!$omp parallel do
  do itor=nt1,nt2
   if (itor == 0 .and. ae_flag == 1) then
     call cgyro_field_ae('v')
   else
     field(:,:,itor) = fcoef(:,:,itor)*field(:,:,itor)
   endif
  enddo

  call timer_lib_out('field')

end subroutine cgyro_field_v

! like cgyro_field_v_notae, but with parametrized start_t
subroutine cgyro_field_v_notae_s(start_t)

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none
  ! ------------------ 
  integer, intent(in) :: start_t
  !
  integer :: itor

  call timer_lib_in('field')

  field_loc(:,:,start_t:nt2) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of H

!$omp parallel do collapse(2) private(ic_loc,iv)
  do itor=start_t,nt2
   do ic=nc1,nc2
     ic_loc = ic-nc1+1
     do iv=1,nv
        field_loc(:,ic,itor) = field_loc(:,ic,itor)+dvjvec_v(:,ic_loc,itor,iv)*cap_h_v(ic_loc,itor,iv)
     enddo
   enddo
  enddo

  call timer_lib_out('field')

  call timer_lib_in('field_com')

  call MPI_ALLREDUCE(field_loc(:,:,start_t:nt2),&
       field(:,:,start_t:nt2),&
       size(field(:,:,start_t:nt2)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

   call timer_lib_out('field_com')
  
  call timer_lib_in('field')

  ! Poisson LHS factors
!$omp parallel do
  do itor=start_t,nt2
     ! assuming  (.not.(itor == 0 .and. ae_flag == 1))
     field(:,:,itor) = fcoef(:,:,itor)*field(:,:,itor)
  enddo

  call timer_lib_out('field')
end subroutine cgyro_field_v_notae_s

! like cgyro_field_v, but skip (itor == 0 .and. ae_flag == 1)
subroutine cgyro_field_v_notae

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none

  if (nt1 == 0 .and. ae_flag == 1) then
     if (nt2>0) then
        call cgyro_field_v_notae_s(1)
     endif
     ! else no-op
  else
     ! don't have to worry about ae_flag, just use all the elements
     call cgyro_field_v_notae_s(nt1)
  endif

end subroutine cgyro_field_v_notae

!
! GPU versions
!

! Note: Not supporting the ae-version of cgyro_field_v_gpu

subroutine cgyro_field_v_notae_s_gpu(start_t)
  use mpi
  use timer_lib
  use cgyro_globals

  implicit none
  ! ------------------ 
  integer, intent(in) :: start_t
  !
  integer :: i_f,itor
  complex :: field_loc_l 

  call timer_lib_in('field')
!$acc data present(cap_h_v)
!$acc data present(field,field_loc)


  ! Poisson and Ampere RHS integrals of H

!$acc parallel loop collapse(3) independent copyin(start_t) &
!$acc&         present(nt2,nc,n_field) default(none)
  do itor=start_t,nt2
   do ic=1,nc
       do i_f=1,n_field
        field_loc(i_f,ic,itor) = (0.0,0.0)
       enddo
   enddo
  enddo

!$acc parallel loop collapse(3) gang private(ic_loc,field_loc_l) copyin(start_t) &
!$acc&         present(dvjvec_v,cap_h_v,field_loc) &
!$acc&         present(nt2,nc1,nc2,n_field,nv) default(none)
  do itor=start_t,nt2
   do ic=nc1,nc2
    do i_f=1,n_field
      ic_loc = ic-nc1+1
      field_loc_l = (0.0,0.0)
!$acc loop vector reduction(+:field_loc_l)
      do iv=1,nv
        field_loc_l = field_loc_l+dvjvec_v(i_f,ic_loc,itor,iv)*cap_h_v(ic_loc,itor,iv)
     enddo
     field_loc(i_f,ic,itor) = field_loc_l
    enddo
   enddo
  enddo

  call timer_lib_out('field')

  call timer_lib_in('field_com')

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update host(field_loc(:,:,start_t:nt2))
#else
!$acc host_data use_device(field_loc,field)
#endif

  call MPI_ALLREDUCE(field_loc(:,:,start_t:nt2),&
       field(:,:,start_t:nt2),&
       size(field(:,:,start_t:nt2)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update device(field(:,:,start_t:nt2))
#else
!$acc end host_data
#endif

  call timer_lib_out('field_com')

  call timer_lib_in('field')
  ! Poisson LHS factors
!$acc parallel loop collapse(3) independent present(fcoef) copyin(start_t) &
!$acc&         present(nt2,nc,n_field) default(none)
  do itor=start_t,nt2
     ! assuming  (.not.(itor == 0 .and. ae_flag == 1))
     do ic=1,nc
       do i_f=1,n_field
        field(i_f,ic,itor) = fcoef(i_f,ic,itor)*field(i_f,ic,itor)
       enddo
     enddo
  enddo

!$acc end data
!$acc end data

  call timer_lib_out('field')

end subroutine cgyro_field_v_notae_s_gpu

! like cgyro_field_v, but skip (itor == 0 .and. ae_flag == 1)
subroutine cgyro_field_v_notae_gpu

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none

  if (nt1 == 0 .and. ae_flag == 1) then
     if (nt2>0) then
        call cgyro_field_v_notae_s_gpu(1)
     endif
     ! else no-op
  else
     ! don't have to worry about ae_flag, just use all the elements
     call cgyro_field_v_notae_s_gpu(nt1)
  endif

end subroutine cgyro_field_v_notae_gpu

!-----------------------------------------------------------------
! Configuration (velocity-distributed) field solve
!-----------------------------------------------------------------
subroutine cgyro_field_c_cpu

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none

  integer :: is,itor
  complex :: my_psi
  
  complex, dimension(nc) :: tmp
  
  call timer_lib_in('field')

  field_loc(:,:,:) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of h

!$omp parallel private(iv_loc,ic)
!$omp do collapse(2) reduction(+:field_loc)
  do itor=nt1,nt2
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     do ic=1,nc
        field_loc(:,ic,itor) = field_loc(:,ic,itor)+dvjvec_c(:,ic,iv_loc,itor)*h_x(ic,iv_loc,itor)
     enddo
   enddo
  enddo
!$omp end do
!$omp end parallel

  call timer_lib_out('field')

  call timer_lib_in('field_com')

  call MPI_ALLREDUCE(field_loc(:,:,:),&
       field(:,:,:),&
       size(field(:,:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  call timer_lib_out('field_com')

  call timer_lib_in('field')

!$omp parallel do private(tmp) shared(field)
  do itor=nt1,nt2
   if (n_field > 2) then
     field(3,:,itor) = field(3,:,itor)*fcoef(3,:,itor)
   endif

   ! Poisson LHS factors
   if (itor == 0 .and. ae_flag == 1) then
    call cgyro_field_ae('c')
   else
     if (n_field > 2) then
        tmp(:) = field(1,:,itor)
        field(1,:,itor) = gcoef(1,:,itor)*field(1,:,itor) + &
                gcoef(4,:,itor)*field(3,:,itor)
        field(2,:,itor) = gcoef(2,:,itor)*field(2,:,itor)
        field(3,:,itor) = gcoef(3,:,itor)*field(3,:,itor) + &
                gcoef(5,:,itor)*tmp(:)
     else
        field(:,:,itor) = gcoef(:,:,itor)*field(:,:,itor)
     endif
   endif
  enddo

!$omp parallel do collapse(2) private(iv_loc,is,ic,my_psi)
  do itor=nt1,nt2
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        my_psi = sum( jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
        cap_h_c(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)+my_psi*z(is)/temp(is)
     enddo
   enddo
  enddo

  call timer_lib_out('field')

end subroutine cgyro_field_c_cpu

! like cgyro_field_c, but assume (my_toroidal == 0 .and. ae_flag == 1)
subroutine cgyro_field_c_ae_cpu

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none

  integer :: is,itor
  complex :: my_psi
  
  call timer_lib_in('field')

  field_loc(:,:,0:0) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of h

!$omp parallel private(iv_loc,ic)
!$omp do collapse(2) reduction(+:field_loc)
  do itor=0,0
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     do ic=1,nc
        field_loc(:,ic,itor) = field_loc(:,ic,itor)+dvjvec_c(:,ic,iv_loc,itor)*h_x(ic,iv_loc,itor)
     enddo
   enddo
  enddo
!$omp end do
!$omp end parallel

  call timer_lib_out('field')

  call timer_lib_in('field_com')

  call MPI_ALLREDUCE(field_loc(:,:,0:0),&
       field(:,:,0:0),&
       size(field(:,:,0:0)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  call timer_lib_out('field_com')

  call timer_lib_in('field')

  do itor=0,0
   if (n_field > 2) then
     field(3,:,itor) = field(3,:,itor)*fcoef(3,:,itor)
   endif

   ! Poisson LHS factors
   call cgyro_field_ae('c')
  enddo

!$omp parallel do collapse(2) private(iv_loc,is,ic,my_psi)
  do itor=0,0
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        my_psi = sum( jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
        cap_h_c(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)+my_psi*z(is)/temp(is)
     enddo
   enddo
  enddo

  call timer_lib_out('field')

end subroutine cgyro_field_c_ae_cpu

#ifdef _OPENACC
subroutine cgyro_field_c_gpu
  use mpi
  use timer_lib
  use cgyro_globals
  implicit none
  integer :: is,i_f,itor
  integer :: itor1,itor2
  complex :: tmp,field_loc_l
  complex :: my_psi

  call timer_lib_in('field')
!$acc data present(h_x,cap_h_c)

!$acc data present(field,field_loc)

  ! Poisson and Ampere RHS integrals of h

!$acc parallel loop collapse(3) independent private(field_loc_l) &
!$acc&         present(dvjvec_c) present(nt1,nt2,nc,n_field,nv1,nv2) default(none)
  do itor=nt1,nt2
   do ic=1,nc
    do i_f=1,n_field
      field_loc_l = (0.0,0.0)    
!$acc loop seq private(iv_loc)
      do iv=nv1,nv2
         iv_loc = iv-nv1+1
         field_loc_l = field_loc_l+dvjvec_c(i_f,ic,iv_loc,itor)*h_x(ic,iv_loc,itor)
      enddo
      field_loc(i_f,ic,itor) = field_loc_l
    enddo
   enddo
  enddo
  call timer_lib_out('field')
  call timer_lib_in('field_com')

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update host(field_loc)
#else
!$acc host_data use_device(field_loc,field)
#endif

  call MPI_ALLREDUCE(field_loc(:,:,:),&
       field(:,:,:),&
       size(field(:,:,:)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update device(field)
#else
!$acc end host_data
#endif

  call timer_lib_out('field_com')
  call timer_lib_in('field')
  if (n_field > 2) then
!$acc parallel loop collapse(2) independent present(fcoef) &
!$acc&         present(nt1,nt2,nc) default(none)
    do itor=nt1,nt2
      do ic=1,nc
       field(3,ic,itor) = field(3,ic,itor)*fcoef(3,ic,itor)
      enddo
     enddo
  endif

  ! Poisson LHS factors
  itor1=nt1
  itor2=nt2

  if (nt1 == 0 .and. ae_flag == 1) then
    ! Note: Called rarely, use the CPu version
!$acc update host(field)
    call cgyro_field_ae('c')
!$acc update device(field)
    ! mark we already processed ==0, nothing else to do there
    itor1=1
  endif

  if (itor1<=itor2) then
     if (n_field > 2) then
!$acc parallel loop collapse(2) independent private(tmp) present(gcoef) &
!$acc&         copyin(itor1,itor2) present(nc) default(none)
        do itor=itor1,itor2
         do ic=1,nc
          tmp = field(1,ic,itor)
          field(1,ic,itor) = gcoef(1,ic,itor)*field(1,ic,itor)+ &
                  gcoef(4,ic,itor)*field(3,ic,itor)
          field(2,ic,itor) = gcoef(2,ic,itor)*field(2,ic,itor)
          field(3,ic,itor) = gcoef(3,ic,itor)*field(3,ic,itor)+ &
                  gcoef(5,ic,itor)*tmp
         enddo
        enddo
     else
!$acc parallel loop collapse(3) independent present(gcoef) &
!$acc&         copyin(itor1,itor2) present(nc,n_field) default(none)
        do itor=itor1,itor2
         do ic=1,nc
          do i_f=1,n_field
            field(i_f,ic,itor) = gcoef(i_f,ic,itor)*field(i_f,ic,itor)
          enddo
         enddo
        enddo
     endif
  endif

!$acc parallel loop collapse(3) gang vector private(iv_loc,is,my_psi) &
!$acc&         present(jvec_c,z,temp,is_v) present(nt1,nt2,nv1,nv2,nc) default(none)
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        my_psi = sum( jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
        cap_h_c(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)+my_psi*z(is)/temp(is)
     enddo
   enddo
  enddo

!$acc end data

!$acc end data

  call timer_lib_out('field')
end subroutine cgyro_field_c_gpu

! like cgyro_field_c, but assume (my_toroidal == 0 .and. ae_flag == 1)
subroutine cgyro_field_c_ae_gpu
  use mpi
  use timer_lib
  use cgyro_globals
  implicit none
  integer :: is,i_f,itor
  complex :: tmp,field_loc_l
  complex :: my_psi

  call timer_lib_in('field')
!$acc data present(h_x,cap_h_c)

!$acc data present(field,field_loc)

  ! Poisson and Ampere RHS integrals of h

!$acc parallel loop collapse(3) independent private(field_loc_l) &
!$acc&         present(dvjvec_c) present(nc,n_field,nv1,nv2) default(none)
  do itor=0,0
   do ic=1,nc
    do i_f=1,n_field
      field_loc_l = (0.0,0.0)    
!$acc loop seq private(iv_loc)
      do iv=nv1,nv2
         iv_loc = iv-nv1+1
         field_loc_l = field_loc_l+dvjvec_c(i_f,ic,iv_loc,itor)*h_x(ic,iv_loc,itor)
      enddo
      field_loc(i_f,ic,itor) = field_loc_l
    enddo
   enddo
  enddo
  call timer_lib_out('field')
  call timer_lib_in('field_com')

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update host(field_loc(:,:,0:0))
#else
!$acc host_data use_device(field_loc,field)
#endif

  call MPI_ALLREDUCE(field_loc(:,:,0:0),&
       field(:,:,0:0),&
       size(field(:,:,0:0)),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

#ifdef DISABLE_GPUDIRECT_MPI
!$acc update device(field(:,:,0:0))
#else
!$acc end host_data
#endif

  call timer_lib_out('field_com')
  call timer_lib_in('field')
  if (n_field > 2) then
!$acc parallel loop collapse(2) independent present(fcoef) present(nc) default(none)
    do itor=0,0
      do ic=1,nc
       field(3,ic,itor) = field(3,ic,itor)*fcoef(3,ic,itor)
      enddo
     enddo
  endif
  ! Poisson LHS factors
  ! Note: Called rarely, use the CPU version
!$acc update host(field(:,:,0:0))
    call cgyro_field_ae('c')
!$acc update device(field(:,:,0:0))

!$acc parallel loop collapse(3) gang vector private(iv_loc,is,my_psi) &
!$acc&         present(jvec_c,z,temp,is_v) present(nv1,nv2,nc) default(none)
  do itor=0,0
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        my_psi = sum( jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
        cap_h_c(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)+my_psi*z(is)/temp(is)
     enddo
   enddo
  enddo

!$acc end data

!$acc end data

  call timer_lib_out('field')
end subroutine cgyro_field_c_ae_gpu

#endif


subroutine cgyro_field_c
  implicit none
#ifdef _OPENACC
   call cgyro_field_c_gpu
#else
   call cgyro_field_c_cpu
#endif
end subroutine cgyro_field_c

! like cgyro_field_c, but only for (itor == 0 .and. ae_flag == 1)
subroutine cgyro_field_c_ae
  implicit none
#ifdef _OPENACC
   call cgyro_field_c_ae_gpu
#else
   call cgyro_field_c_ae_cpu
#endif
end subroutine cgyro_field_c_ae

!-----------------------------------------------------------------
! Adiabatic electron field solves for n=0
! Can only be called if itor==0
!-----------------------------------------------------------------
subroutine cgyro_field_ae(space)

  use cgyro_globals
  use timer_lib  

  implicit none

  character(len=1), intent(in) :: space
  integer :: ir,i,j
  complex, dimension(n_theta) :: pvec_in,pvec_out

  if (space == 'c') then
     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_mode == 0) then
           field(1,ic_c(ir,:),0) = 0.0
        else
           do i=1,n_theta
              pvec_out(i) = 0.0
              pvec_in(i)  = field(1,ic_c(ir,i),0)
           enddo
           do j=1,n_theta
              do i=1,n_theta
                 pvec_out(i) = pvec_out(i)+xzf(ir,i,j)*pvec_in(j)
              enddo
           enddo
           do i=1,n_theta
              field(1,ic_c(ir,i),0) = pvec_out(i)
           enddo
        endif
     enddo
  else
     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_mode == 0) then
           field(1,ic_c(ir,:),0) = 0.0
        else
           do i=1,n_theta
              pvec_out(i) = 0.0
              pvec_in(i)  =  field(1,ic_c(ir,i),0)
           enddo
           do j=1,n_theta
              do i=1,n_theta
                 pvec_out(i) = pvec_out(i)+hzf(ir,i,j)*pvec_in(j)
              enddo
           enddo
           do i=1,n_theta
              field(1,ic_c(ir,i),0) = pvec_out(i)
           enddo
        endif
     enddo
  endif

end subroutine cgyro_field_ae
