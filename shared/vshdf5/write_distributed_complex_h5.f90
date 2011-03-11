!------------------------------------------------------
! write_distributed_complex.f90
!
! PURPOSE:
!  Control merged output of complex distributed array.
!------------------------------------------------------

subroutine write_distributed_complex_h5(vname,rGid,r3Did,&
                     n_fn,n1,n2,n3,fn,plot3d,h5in,h5err)

  use math_constants
  use hdf5_mod
  use hdf5_api
  use gyro_globals, only : &
       q, &
       n0,&
       d_n,&
       n_n,&
       n_n_1,&
       n_proc_1,&
       n_theta_plot,&
       omega_exp,&
       t_current,&
       debug_flag,&
       recv_status,&
       data_step,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err, &
       electron_method,&
       n_alpha_fine,&
       n_alpha_plot

  !------------------------------------------------------
  !   mom: n1,n2,n3=n_theta_plot,n_x,n_kinetic
  !   phi: n1,n2,n3=n_theta_plot,n_x,n_field
  !------------------------------------------------------
  implicit none
  !
  character*(*), intent(in) :: vname
  integer(HID_T), intent(in) :: rGid,r3Did
  integer, intent(in) :: n_fn,n1,n2,n3
  complex, intent(in) :: fn(n_fn)
  logical, intent(in) :: plot3d
  character(128) :: tempVarName
  character(3) :: n_name
  character(1) :: ikin_name
  type(hdf5InOpts), intent(inout) :: h5in
  type(hdf5ErrorType), intent(inout) :: h5err
  !
  integer :: data_loop
  integer :: i_group_send
  integer :: i_send, iphi, istart,nn,i,ikin,in, ix,nphi
  !
  complex :: fn_recv(n_fn), c_i
  complex, dimension(:,:,:,:), allocatable :: buffn
  real, dimension(:,:,:,:), allocatable:: real_buff
  real, dimension(:,:), allocatable:: alpha_loc
  logical :: iscoarse

  !------------------------------------------------------
  include 'mpif.h'
  c_i=(0,1)

  if(n1==n_theta_plot) then
     iscoarse=.true.
  else
     iscoarse=.false.
  endif

  allocate(buffn(0:n1,n2,n3,n_n)); buffn=0.

     do in=1,n_n
  !WRITE(*,*) "in ", in, i_proc
        !-----------------------------------------
        ! Subgroup collector:
        !
        i_group_send = (in-1)/n_n_1
        if (i_group_send /= 0) then
           i_send = i_group_send*n_proc_1
           if (i_proc == 0) then
              call MPI_RECV(fn_recv,&
                   n_fn, MPI_DOUBLE_COMPLEX, i_send, in,&
                   GYRO_COMM_WORLD, recv_status, i_err)
           else if (i_proc == i_send) then
              call MPI_SEND(fn,&
                   n_fn, MPI_DOUBLE_COMPLEX, 0, in,&
                   GYRO_COMM_WORLD, i_err)
           endif
        else
           fn_recv(:) = fn(:)
        endif
        
        if (i_proc == 0) then
           buffn(0:n1-1,:,:,in)=reshape(fn_recv,(/n1,n2,n3/))
        endif
            

     enddo ! in
     !-----------------------------------------
     if (i_proc /= 0) return
     !-----------------------------------------
     !-----------------------------------------
     ! Dump each mode in the same format that gyro does
     !-----------------------------------------
     if (iscoarse) then
       do in=1,n_n
         WRITE(n_name,fmt='(i3.3)') in
         tempVarName=trim(vname)//"_real"//n_name
         call dump_h5(rGid,trim(tempVarName),real(buffn(0:n1-1,:,:,in)),h5in,h5err)
         tempVarName=trim(vname)//"_imag"//n_name
         call dump_h5(rGid,trim(tempVarName),aimag(buffn(0:n1-1,:,:,in)),h5in,h5err)
       enddo ! in
     else
!      if (debug_flag == 1) then
        do in=1,n_n
         WRITE(n_name,fmt='(i3.3)') in
         tempVarName=trim(vname)//"fine_real"//n_name
         call dump_h5(rGid,trim(tempVarName),real(buffn(0:n1-1,:,:,in)),h5in,h5err)
         tempVarName=trim(vname)//"fine_imag"//n_name
         call dump_h5(rGid,trim(tempVarName),aimag(buffn(0:n1-1,:,:,in)),h5in,h5err)
       enddo ! in
!      endif
     endif
     if(.not.plot3d) then
       deallocate(buffn)
       return
     endif
     !-----------------------------------------
     ! Apply boundary conditions
     !-----------------------------------------
     do in=1,n_n
       nn=n0+(in-1)*d_n
       !apply theta BC: z_n(r,,2*pi) = z_n(r,0)exp(I*n*(nu(r,2*pi)-nu(r,0)))
       !with nu(r,2*pi) - nu(r,0) = -2*pi*q by definition
       ! phase[*] = EXP(-2*!PI*C_I*n[i_n]*profile_data.q[*])
       do ix=1,n2
         buffn(n1,ix,:,in)=buffn(0,ix,:,in)*exp(-2.*pi*c_i*nn*q(ix))
       enddo
     enddo ! in

     !-----------------------------------------
     ! Tranform into real space
     !   mom: n1,n2,n3=n_theta_plot,n_x,n_kinetic
     !   phi: n1,n2,n3=n_theta_plot,n_x,n_field
     !-----------------------------------------
     if (iscoarse) then
       nphi=n_alpha_plot
     else
       nphi=n_alpha_fine
     endif
     allocate(real_buff(0:n1,n2,n3,nphi))
     allocate(alpha_loc(0:n1,n2))
!sv     Default for n0=30, because of k_rho_s scaling.
!       I think we always want to count over all toroidal modes calculated
!       which means below, going from istart=1 to n_n (toroidal grid).
!     if (n0==0) then
         istart=1
         do iphi=1,nphi
           real_buff(:,:,:,iphi)=real(buffn(:,:,:,1))
         enddo
!     else
!         istart=2
!         real_buff(:,:,:,:)=0.
!         !where does real_buff come from in thise case??/
!     endif
     do iphi=1,nphi
       !Get alpha coordinate on either the coarse or fine mesh.
       ! Include doppler shift here
       if (iscoarse) then
         alpha_loc=alpha_phi(:,:,iphi)+omega_exp*t_current
       else
         alpha_loc=alpha_phi_fine(:,:,iphi)+omega_exp*t_current
       endif
       do in=istart,n_n
          nn=n0+(in-1)*d_n
          do ikin=1,n3
            real_buff(:,:,ikin,iphi)=real_buff(:,:,ikin,iphi)&
                  +2.*real(buffn(:,:,ikin,in)*exp(-c_i*nn*alpha_loc(:,:)))
          enddo
       enddo
     enddo
        
     

     deallocate(buffn,alpha_loc)


     ! Mapping of the variable names to array indices depends on input types
      if (iscoarse) then
       do ikin=1,n3
         if (trim(vname) /= "phi") then 
           ! See gyro_select_methods for understanding this logic
           if(electron_method==2 .and. ikin==n3) THEN
              tempVarName=trim(vname)//"_electron"
           elseif(electron_method==3) THEN
              tempVarName=trim(vname)//"_electron"
           else
              write(ikin_name,fmt='(i1.1)') ikin
              tempVarName=trim(vname)//"_ion"//ikin_name
           endif
         else
           if(ikin==1) tempVarName="phi"
           if(ikin==2) tempVarName="A_par"
           if(ikin==3) tempVarName="B_par"
         endif
         call dump_h5(r3Did,trim(tempVarName),real_buff(:,:,ikin,:),h5in,h5err)
       enddo
     else
       ! Dump each phi slice as a separate variable
       do iphi=1,nphi
         do ikin=1,n3
           if (trim(vname) /= "phi") then 
             if(electron_method==2 .and. ikin==n3) THEN
                tempVarName=trim(vname)//"_electron"
             elseif(electron_method==3) THEN
                tempVarName=trim(vname)//"_electron"
             else
                write(ikin_name,fmt='(i1.1)') ikin-1
                tempVarName=trim(vname)//"_ion"//ikin_name
             endif
           else
             if(ikin==1) tempVarName="phi"
             if(ikin==2) tempVarName="A_par"
             if(ikin==3) tempVarName="B_par"
           endif
           write(n_name,fmt='(i2.2)') iphi
           tempVarName=trim(tempVarName)//"_phi"//TRIM(n_name)
           call dump_h5(r3Did,trim(tempVarName),real_buff(:,:,ikin,iphi),h5in,h5err)
         enddo
       enddo
     endif

     deallocate(real_buff)

return
end subroutine write_distributed_complex_h5
