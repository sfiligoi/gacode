subroutine cgyro_check_memory(data_file)

  use cgyro_globals

  implicit none
  
  character (len=*), intent(in) :: data_file

  if(i_proc == 0) then

     open(unit=io,file=data_file,status='replace')

     write(io,*) 'grids'
     write(io,*)
     
     call alloc_add(io,n_energy,8,'energy')
     call alloc_add(io,n_energy,8,'w_e')
     call alloc_add(io,n_energy**2,8,'e_deriv1_mat')
     call alloc_add(io,n_energy**2,8,'e_deriv2_mat')

     call alloc_add(io,n_xi,8,'xi')
     call alloc_add(io,n_xi,8,'w_xi')
     call alloc_add(io,n_xi**2,8,'xi_lor_mat')
     call alloc_add(io,n_xi**2,8,'xi_deriv_mat')
     call alloc_add(io,n_xi**2,8,'xi_upderiv_mat')
     call alloc_add(io,n_xi,4,'indx_xi')

     call alloc_add(io,n_radial,4,'px')

     call alloc_add(io,n_theta,8,'theta')
     call alloc_add(io,3*n_theta,4,'thcyc')
     call alloc_add(io,nc*(2*nup_theta+1),4,'icd_c')
     call alloc_add(io,nc*(2*nup_theta+1),16,'dtheta')
     call alloc_add(io,nc*(2*nup_theta+1),16,'dtheta_up')
     call alloc_add(io,2*nup_theta+1,8,'cderiv')
     call alloc_add(io,2*nup_theta+1,8,'uderiv')
     call alloc_add(io,n_radial,8,'spec_uderiv')

     write(io,*)
     write(io,*) 'mpi grids'
     write(io,*)

     call alloc_add(io,nv,4,'ie_v')
     call alloc_add(io,nv,4,'ix_v')
     call alloc_add(io,nc,4,'ir_c')
     call alloc_add(io,nc,4,'it_c')
     call alloc_add(io,n_radial*n_theta,4,'ic_c')
     call alloc_add(io,n_energy*n_xi*n_species,4,'iv_v')

     write(io,*)
     write(io,*) 'gk rhs'
     write(io,*)

     call alloc_add(io,4*nc*nv_loc,16,'rhs')
     call alloc_add(io,n_theta*n_species,8,'omega_stream')
     call alloc_add(io,n_theta*n_species,8,'omega_trap')
     call alloc_add(io,n_theta*n_species,8,'omega_rdrift')
     call alloc_add(io,n_theta*n_species,8,'omega_adrift')
     call alloc_add(io,n_theta*n_species,8,'omega_aprdrift')
     call alloc_add(io,n_theta*n_species,8,'omega_cdrift')
     call alloc_add(io,n_theta,8,'omega_gammap')
     call alloc_add(io,nc*nv_loc,16,'omega_cap_h')
     call alloc_add(io,nc*nv_loc,16,'omega_h')
     call alloc_add(io,n_field*nc*nv_loc,16,'omega_s')

     write(io,*)
     write(io,*) 'geometry'
     write(io,*)

     call alloc_add(io,n_theta,8,'w_theta')
     call alloc_add(io,n_theta,8,'b_mag')
     call alloc_add(io,nc,8,'k_perp')
     if(equilibrium_model == 1) then
        call alloc_add(io,8,8,'geo_yin')
     else if(equilibrium_model == 2) then
     else
        call alloc_add(io,8,8,'geo_yin')
     endif
        call alloc_add(io,8*(geo_ny+1),8,'geo_yin')

     write(io,*)
     write(io,*) 'fields and field solve'
     write(io,*)

     call alloc_add(io,n_field*nc,16,'field')
     call alloc_add(io,n_field*nc,16,'field_loc')
     call alloc_add(io,n_field*nc,16,'field_old')
     call alloc_add(io,n_field*nc,16,'field_old2')
     call alloc_add(io,n_field*nc,16,'field_old3')
     call alloc_add(io,n_field*nc,8,'fcoef')
     if(n_field < 3) then
        call alloc_add(io,n_field*nc,8,'gcoef')
     else
        call alloc_add(io,5*nc,8,'gcoef')
     endif
     call alloc_add(io,nc,8,'sum_den_x')
     call alloc_add(io,nc,8,'sum_cur_x')

     if (n == 0 .and. ae_flag == 1) then
        call alloc_add(io,n_radial*n_theta**2,8,'hzf')
        call alloc_add(io,n_radial*n_theta**2,8,'xzf')
        call alloc_add(io,n_theta,8,'pvec_in')
        call alloc_add(io,n_theta,8,'pvec_outr')
        call alloc_add(io,n_theta,8,'pvec_outi')
     endif

     write(io,*)
     write(io,*) 'fluxes and diagnostics'
     write(io,*)

     call alloc_add(io,n_radial*n_species,16,'moment_loc')
     call alloc_add(io,n_radial*n_species,16,'moment')
     call alloc_add(io,(n_radial/box_size)*n_theta,8,'thetab')
     call alloc_add(io,(n_radial/box_size)*n_theta,16,'f_balloon')
     call alloc_add(io,n_radial*n_species,8,'flux')
     call alloc_add(io,n_radial*n_species,8,'flux_loc')

     write(io,*)
     write(io,*) 'distributions'
     write(io,*)

     call alloc_add(io,nc*nv_loc,16,'h_x')
     call alloc_add(io,nc*nv_loc,16,'g_x')
     call alloc_add(io,nc*nv_loc,16,'psi')
     call alloc_add(io,nc*nv_loc,16,'h0_x')
     call alloc_add(io,nc*nv_loc,16,'cap_h_c')
     call alloc_add(io,nc*nv_loc,16,'cap_h_ct')

     write(io,*)
     write(io,*) 'bessel functions'
     write(io,*)

     call alloc_add(io,n_field*nc*nv_loc,8,'jvec_c')
     call alloc_add(io,n_field*nc_loc*nv,8,'jvec_v')

     call alloc_add(io,nc_loc*nv,16,'cap_h_v')
     call alloc_add(io,nc_loc*nv,16,'cap_h_v_prime')

     if(nonlinear_flag == 1) then
        write(io,*)
        write(io,*) 'nonlinear'
        write(io,*)
        if(nonlinear_method == 1) then
           call alloc_add(io,nc*nsplit*n_toroidal,16,'f_nl')
           call alloc_add(io,nc*nsplit*n_toroidal,16,'g_nl')
        else
           call alloc_add(io,n_radial*nsplit*n_toroidal,16,'f_nl')
           call alloc_add(io,n_radial*nsplit*n_toroidal,16,'g_nl')
           nx0 = n_radial
           ny0 = 2*n_toroidal-1
           nx = (3*nx0)/2
           ny = (3*ny0)/2
           call alloc_add(io,(ny/2_1)*nx,16,'fx')
           call alloc_add(io,(ny/2_1)*nx,16,'gx')
           call alloc_add(io,(ny/2_1)*nx,16,'fy')
           call alloc_add(io,(ny/2_1)*nx,16,'gy')
           call alloc_add(io,ny*nx,8,'ux')
           call alloc_add(io,ny*nx,8,'vx')
           call alloc_add(io,ny*nx,8,'uy')
           call alloc_add(io,ny*nx,8,'vy')
           call alloc_add(io,ny*nx,8,'uv')
        endif
     endif

     write(io,*)
     write(io,*) 'collisions'
     write(io,*)
     
     call alloc_add(io,nv*nv*nc_loc,8,'cmat')
     call alloc_add(io,nv,16,'bvec')
     call alloc_add(io,nv,16,'cvec')

     if(implicit_flag == 1) then
        write(io,*)
        write(io,*) 'implicit gk'
        write(io,*)
        call alloc_add(io,nc*nv_loc,16,'gkvec')
        call alloc_add(io,(nc*n_field)**2,16,'fieldmat')
        call alloc_add(io,nc*n_field,4,'idfield')
        call alloc_add(io,nc*n_field,4,'i_piv_field')
        call alloc_add(io,nc*n_field,16,'fieldvec')
        call alloc_add(io,nc*n_field,16,'fieldvec_loc')
        
        call alloc_add(io,(2*nup_theta+1)*nc*10*nv_loc,16,'gksp_mat')
        call alloc_add(io,(2*nup_theta+1)*nc*10*nv_loc*2,4,'gksp_indx')
        call alloc_add(io,10*nv_loc,8,'gksp_cntl')
        call alloc_add(io,20*nv_loc,4,'gksp_icntl')
        call alloc_add(io,20*nv_loc,4,'gksp_keep')
        call alloc_add(io,nc,16,'gksvec')
        call alloc_add(io,2*nc,16,'gkwvec')

     endif

     write(io,*)
     write(io,*) '---------------------'
     write(io,'(f7.3,a,3x,a)') total_memory/1048576.0,' MB'

     close(io)
     
  end if

end subroutine cgyro_check_memory

!------------------------------------------------
! alloc_add.f90
!
! PURPOSE:
!  Primitive allocation addition routine.
!------------------------------------------------

subroutine alloc_add(my_io,n_size,bytes,name)

  use cgyro_globals

  implicit none
  !
  integer, intent(in) :: my_io
  integer, intent(in) :: n_size
  integer, intent(in) :: bytes
  character (len=*), intent(in) :: name
  !
  real :: this_memory

  if (i_proc == 0) then
     
     this_memory  = 1.0*n_size*bytes
     total_memory = total_memory+this_memory 
     
     write(my_io,10) this_memory/1048576.0,' MB',name
     
  endif

10 format(f7.3,a,3x,a)

end subroutine alloc_add
