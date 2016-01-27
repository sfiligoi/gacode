subroutine cgyro_check_memory(datafile)

  use cgyro_globals

  implicit none

  character (len=*), intent(in) :: datafile
  integer :: mult
  
  if (i_proc == 0) then

     open(unit=io,file=datafile,status='replace')

     total_memory = 0
     write(io,*) ' == Non-Distributed =='

     write(io,*)
     write(io,*) 'Pointers'
     write(io,*)    
     call alloc_add(io,nv*4.0,'ie_v')
     call alloc_add(io,nv*4.0,'ix_v')
     call alloc_add(io,nc*4.0,'ir_c')
     call alloc_add(io,nc*4.0,'it_c')
     call alloc_add(io,n_radial*n_theta*4.0,'ic_c')
     call alloc_add(io,n_energy*n_xi*n_species*4.0,'iv_v')
     call alloc_add(io,nc*(2*nup_theta+1)*4.0,'icd_c')
     call alloc_add(io,nc*(2*nup_theta+1)*16.0,'dtheta')
     call alloc_add(io,nc*(2*nup_theta+1)*16.0,'dtheta_up')

     write(io,*)
     write(io,*) 'Fields and field solve'
     write(io,*)
     call alloc_add(io,n_field*nc*16.0,'field')
     call alloc_add(io,n_field*nc*16.0,'field_loc')
     call alloc_add(io,n_field*nc*16.0,'field_old')
     call alloc_add(io,n_field*nc*16.0,'field_old2')
     call alloc_add(io,n_field*nc*16.0,'field_old3')
     call alloc_add(io,n_field*nc*8.0,'fcoef')
     if (n_field < 3) then
        call alloc_add(io,n_field*nc*8.0,'gcoef')
     else
        call alloc_add(io,5*nc*8.0,'gcoef')
     endif
     call alloc_add(io,nc*8.0,'sum_den_x')
     call alloc_add(io,nc*8.0,'sum_cur_x')

     if (nonlinear_flag == 1) then
        write(io,*)
        write(io,*) 'Nonlinear'
        write(io,*)
        ! nsplit * n_toroidal = nv_loc * n_theta
        if (nonlinear_method /= 1) then
           nx0 = n_radial
           ny0 = 2*n_toroidal-1
           nx = (3*nx0)/2
           ny = (3*ny0)/2
           call alloc_add(io,(ny/2+1)*nx*16.0,'fx')
           call alloc_add(io,(ny/2+1)*nx*16.0,'gx')
           call alloc_add(io,(ny/2+1)*nx*16.0,'fy')
           call alloc_add(io,(ny/2+1)*nx*16.0,'gy')
           call alloc_add(io,ny*nx*8.0,'ux')
           call alloc_add(io,ny*nx*8.0,'vx')
           call alloc_add(io,ny*nx*8.0,'uy')
           call alloc_add(io,ny*nx*8.0,'vy')
           call alloc_add(io,ny*nx*8.0,'uv')
        endif
     endif

     write(io,*)
     write(io,*) 'TOTAL'
     write(io,10) total_memory/1e6,' MB [persistent memory per core]'
     write(io,*)

     total_memory = 0
     write(io,*) ' == Distributed =='

     if (implicit_flag == 1) then
        write(io,*)
        write(io,*) 'Implicit time advance'
        write(io,*)
        call alloc_add(io,nc*nv_loc*16,'gkvec')
        call alloc_add(io,(2*nup_theta+1)*nc*10*nv_loc*16.0,'gksp_mat')
        call alloc_add(io,(2*nup_theta+1)*nc*10*nv_loc*2*4.0,'gksp_indx')
        call alloc_add(io,10*nv_loc*8.0,'gksp_cntl')
        call alloc_add(io,20*nv_loc*4.0,'gksp_icntl')
        call alloc_add(io,20*nv_loc*4.0,'gksp_keep')
     endif

     write(io,*)
     write(io,*) 'Distribution-like arrays'
     write(io,*)
     call alloc_add(io,4*nc*nv_loc*16.0,'rhs')
     call alloc_add(io,n_field*nc*nv_loc*16.0,'omega_s')
     call alloc_add(io,nc*nv_loc*16.0,'omega_cap_h')
     call alloc_add(io,nc*nv_loc*16.0,'omega_h')
     call alloc_add(io,nc*nv_loc*16.0,'h_x')
     call alloc_add(io,nc*nv_loc*16.0,'g_x')
     call alloc_add(io,nc*nv_loc*16.0,'psi')
     call alloc_add(io,nc*nv_loc*16.0,'h0_x')
     call alloc_add(io,nc*nv_loc*16.0,'cap_h_c')
     call alloc_add(io,nc*nv_loc*16.0,'cap_h_ct')
     call alloc_add(io,n_field*nc*nv_loc*8.0,'jvec_c')
     call alloc_add(io,n_field*nc_loc*nv*8.0,'jvec_v')
     call alloc_add(io,nc_loc*nv*16.0,'cap_h_v')
     call alloc_add(io,nc_loc*nv*16.0,'cap_h_v_prime')

     if (nonlinear_flag == 1) then
        write(io,*)
        write(io,*) 'Nonlinear bracket'
        write(io,*)
        ! nsplit * n_toroidal = nv_loc * n_theta
        if (nonlinear_method == 1) then
           call alloc_add(io,nc*nsplit*n_toroidal*16.0,'f_nl')
           call alloc_add(io,nc*nsplit*n_toroidal*16.0,'g_nl')
        else
           call alloc_add(io,n_radial*nsplit*n_toroidal*16.0,'f_nl')
           call alloc_add(io,n_radial*nsplit*n_toroidal*16.0,'g_nl')
        endif
     endif

     write(io,*)
     write(io,*) 'Collision operator'
     write(io,*)

     call alloc_add(io,(8.0*nv)*nv*nc_loc,'cmat')

     if (implicit_flag == 1) then
        write(io,*)
        write(io,*) 'implicit gk'
        write(io,*)
        call alloc_add(io,(nc*n_field)**2*16.0,'fieldmat')
        call alloc_add(io,nc*n_field*4.0,'idfield')
        call alloc_add(io,nc*n_field*4.0,'i_piv_field')
        call alloc_add(io,nc*n_field*16.0,'fieldvec')
        call alloc_add(io,nc*n_field*16.0,'fieldvec_loc')
        call alloc_add(io,nc*16.0,'gksvec')
        call alloc_add(io,2*nc*16.0,'gkwvec')
     endif

     write(io,*)
     write(io,*) 'TOTAL'
     if (test_flag == 1) then
        write(io,10) total_memory/1e9,&
             ' GB [per toroidal mode ; halved with every doubling of cores]'
        write(io,*) ' '
        do mult=2,16,2
           write(io,20) total_memory/1e9/mult,' GB ',mult*n_toroidal,' cores'
        enddo           
     else
        write(io,10) total_memory/1e9,' GB [per core]'
     endif    
     close(io)

  end if

10 format(t2,f8.3,a,3x,a)
20 format(t2,f8.3,a,i5,a)

end subroutine cgyro_check_memory

!------------------------------------------------
! alloc_add.f90
!
! PURPOSE:
!  Primitive allocation addition routine.
!------------------------------------------------

subroutine alloc_add(my_io,bytes,name)

  use cgyro_globals

  implicit none
  !
  integer, intent(in) :: my_io
  real, intent(in) :: bytes
  character (len=*), intent(in) :: name
 
  if (i_proc == 0) then
     
     total_memory = total_memory+bytes

     if (bytes < 1e8) then 
        write(my_io,10) bytes/1e6,' MB',name
     else
        write(my_io,10) bytes/1e9,' GB',name
     endif
     
  endif

10 format(t2,f8.3,a,3x,a)
  
end subroutine alloc_add
