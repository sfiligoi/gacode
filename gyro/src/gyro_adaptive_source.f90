!-----------------------------------------------------------
! gyro_adaptive_source.f90 [caller: gyro_fulladvance]
!
! PURPOSE:
!  Compute long-wavelength adaptive source term,
!  as used in gyro_rhs_total.  Also compute diagnostics.
!-----------------------------------------------------------

subroutine gyro_adaptive_source

  use mpi
  use gyro_globals
  use gyro_pointers

  !----------------------------------------------------
  implicit none
  !
  integer :: p
  !
  real, dimension(n_kinetic,n_energy,n_x) :: h0_loc,h0_mod
  real, dimension(:,:,:), allocatable :: src_3
  real, dimension(:,:,:,:), allocatable :: src_4
  !----------------------------------------------------

  !------------------
  ! Total RHS source
  !------------------

  h0_loc(:,:,:) = 0.0

  if (n_1(in_1) == 0) then
     p_nek_loc = 0

     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   

        ck = class(k)

        do i=1,n_x
           do m=1,n_stack
              h0_loc(:,ie,i) = h0_loc(:,ie,i)+real(h(m,i,p_nek_loc,:))* &
                   w_lambda(i,k)*d_tau(ck)
           enddo ! m
        enddo ! i

     enddo ! p_nek
  endif

  call MPI_ALLREDUCE(h0_loc,&
       h0_mod,&
       size(h0_mod),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  call MPI_BCAST(h0_mod,&
       size(h0_mod),&
       MPI_DOUBLE_PRECISION,&
       0,&
       NEW_COMM_2,&
       i_err)

  !----------------------------------------------------------------
  ! Compute 'simple' density and energy moments (diagnostic) from
  ! fluctuating h:
  !
  h0_n(:,:) = 0.0
  h0_e(:,:) = 0.0
  !
  do is=1,n_kinetic
     do ie=1,n_energy
        h0_n(is,:) = h0_n(is,:)+h0_mod(is,ie,:)*w_energy(ie)
        h0_e(is,:) = h0_e(is,:)+h0_mod(is,ie,:)*w_energy(ie)*energy(ie)
     enddo
  enddo
  !----------------------------------------------------------------

  h0_eq(:,:,:)  = 0.0
  source_n(:,:) = 0.0
  source_e(:,:) = 0.0

  if (source_method == 2) then

     !----------------------------------------------------------------
     ! Determine the (long-wavelength) finite-element amplitudes:
     !
     allocate(src_3(n_kinetic,n_energy,n_lump))
     src_3(:,:,:) = 0.0

     do p=1,n_lump
        do i=1,n_x
           src_3(:,:,p) = src_3(:,:,p)+b_src(i,p)*h0_mod(:,:,i)
        enddo
     enddo
     !----------------------------------------------------------------

     !----------------------------------------------------------------
     ! Reconstruct h0_eq, the long-wavelength part of h0_mod.
     !
     do is=1,n_kinetic 
        do ie=1,n_energy

           call DGETRS('N',&
                n_lump,&
                1,&
                m_src,&
                n_lump,&
                src_piv,&
                src_3(is,ie,:),&
                n_lump,&
                info)

        enddo
     enddo

     do p=1,n_lump
        do i=1,n_x
           h0_eq(:,:,i) = h0_eq(:,:,i)+src_3(:,:,p)*b_src(i,p)
        enddo
     enddo
     deallocate(src_3)
     !----------------------------------------------------------------

  endif

  if (source_method == 3) then

     !----------------------------------------------------------------
     ! Determine the (long-wavelength) finite-element amplitudes:
     !
     allocate(src_4(n_stack,n_lump,n_nek_loc_1,n_kinetic))
     src_4 = 0.0

     if (n_1(in_1) == 0) then
        do p=1,n_lump
           do i=1,n_x
              src_4(:,p,:,:) = src_4(:,p,:,:)+b_src(i,p)*h(:,i,:,:)
           enddo
        enddo
     endif
     !----------------------------------------------------------------

     !----------------------------------------------------------------
     ! Reconstruct h_source, the long-wavelength part of h.
     !
     do is=1,n_kinetic 
        do m=1,n_stack
           do p_nek_loc=1,n_nek_loc_1
              call DGETRS('N',&
                   n_lump,&
                   1,&
                   m_src,&
                   n_lump,&
                   src_piv,&
                   src_4(m,:,p_nek_loc,is),&
                   n_lump,&
                   info)
           enddo
        enddo
     enddo

     h_source = 0.0
     do p=1,n_lump
        do i=1,n_x
           h_source(:,i,:,:) = h_source(:,i,:,:)+src_4(:,p,:,:)*b_src(i,p)
        enddo
     enddo
     deallocate(src_4)

     p_nek_loc = 0

     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   
        ck = class(k)

        if (n_1(in_1) == 0) then
           do i=1,n_x
              do m=1,n_stack
                 h0_loc(:,ie,i) = h0_loc(:,ie,i)+real(h_source(m,i,p_nek_loc,:))* &
                      w_lambda(i,k)*d_tau(ck)
              enddo ! m
           enddo ! i
        endif

     enddo ! p_nek

     call MPI_ALLREDUCE(h0_loc,&
          h0_eq,&
          size(h0_eq),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     call MPI_BCAST(h0_eq,&
          size(h0_eq),&
          MPI_DOUBLE_PRECISION,&
          0,&
          NEW_COMM_2,&
          i_err)

     !----------------------------------------------------------------

  endif

  !----------------------------------------------------------------
  ! Diagnostics (not used in time evolution)
  !
  do is=1,n_kinetic
     do ie=1,n_energy
        source_n(is,:) = source_n(is,:)-nu_source*h0_eq(is,ie,:)* &
             w_energy(ie)
        source_e(is,:) = source_e(is,:)-nu_source*h0_eq(is,ie,:)* &
             w_energy(ie)*energy(ie)
     enddo
  enddo
  !----------------------------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_adaptive_source done]'
  endif

end subroutine gyro_adaptive_source
