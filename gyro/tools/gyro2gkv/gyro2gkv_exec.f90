!--------------------------------------------------------------
! main.f90
!
! PURPOSE:
!  Driver for the gyro2gkv converter.
!
! NOTES:
!  Use "gyro2gkv -h" to see summary of output.
!   
! REVISIONS
! 17 Mar 04: jc
!  Documented.
! 10 May 06: jc
!  Rewritten to read profiles_vugyro.out.
! 13 Feb 07: ew
!  Added flag -f to output full potentials
!---------------------------------------------------------------

program main

  use gyro2gkv_globals

  implicit none

  include 'netcdf.inc'

  !------------------------------------------------------------------
  ! netCDF variable IDs
  !
  integer :: n_x_NC
  integer :: n_theta_NC
  integer :: n_time_NC
  integer :: n_field_NC
  integer :: n_kinetic_NC
  integer :: n_n_NC
  integer :: n_ri_NC
  !
  integer :: q_NC
  integer :: r_NC
  integer :: theta_NC
  integer :: t_NC
  integer :: n_NC
  integer :: diffp_i_NC
  integer :: difft_i_NC
  integer :: diffp_i_em_NC
  integer :: difft_i_em_NC
  integer :: phi_NC
  integer :: phi_full_NC
  integer :: a_parallel_NC
  integer :: a_parallel_full_NC
  integer :: moment_n_NC
  integer :: moment_n_full_NC
  integer :: moment_e_NC
  integer :: moment_e_full_NC
  integer :: phi_r0_NC
  !
  integer :: err
  integer :: fileid
  !
  integer :: dd(5)
  !------------------------------------------------------------------

  call read_all

  !------------------------------------------------------------------
  ! Open the netCDF file
  !
  err = NF_CREATE(trim(tag(3)),NF_CLOBBER,fileid)
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Dimensions
  !
  err = NF_DEF_DIM(fileid,'r',n_x,n_x_NC)
  err = NF_DEF_DIM(fileid,'theta',n_theta_plot,n_theta_NC)
  err = NF_DEF_DIM(fileid,'t',n_time,n_time_NC)
  err = NF_DEF_DIM(fileid,'species',n_kinetic,n_kinetic_NC)
  err = NF_DEF_DIM(fileid,'n',n_n,n_n_NC)
  err = NF_DEF_DIM(fileid,'ri',2,n_ri_NC)
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Global attributes
  !
  err = NF_PUT_ATT_TEXT(fileid,NF_GLOBAL,"CodeID",4,"GYRO")
  err = NF_PUT_ATT_TEXT(fileid,NF_GLOBAL,"CodePI",20,"candy@fusion.gat.com")
  err = NF_PUT_ATT_TEXT(fileid,NF_GLOBAL,"RunID-version",lt(1),trim(tag(1)))
  err = NF_PUT_ATT_TEXT(fileid,NF_GLOBAL,"RunID-date",lt(2),trim(tag(2)))
  err = NF_PUT_ATT_TEXT(fileid,NF_GLOBAL,"RunID-case",lt(3),trim(tag(3)))
  err = NF_PUT_ATT_TEXT(fileid,NF_GLOBAL,"RunID-boundary",lt(4),trim(tag(4)))
  err = NF_PUT_ATT_TEXT(fileid,NF_GLOBAL,"RunID-model",lt(5),trim(tag(5)))
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! 1D variables and their attributes
  !
  dd(1) = n_x_NC
  err = NF_DEF_VAR(fileid,'q',NF_REAL,1,dd,q_NC)
  err = NF_PUT_ATT_TEXT(fileid,q_NC,"long_name",13,"safety factor")
  err = NF_PUT_ATT_TEXT(fileid,q_NC,"units",13,"dimensionless")

  dd(1) = n_x_NC
  err = NF_DEF_VAR(fileid,'r',NF_REAL,1,dd,r_NC)
  err = NF_PUT_ATT_TEXT(fileid,r_NC,"long_name",1,"r")
  err = NF_PUT_ATT_TEXT(fileid,r_NC,"units",1,"a")

  dd(1) = n_theta_NC
  err = NF_DEF_VAR(fileid,'theta',NF_REAL,1,dd,theta_NC)
  err = NF_PUT_ATT_TEXT(fileid,theta_NC,"long_name",5,"theta")
  err = NF_PUT_ATT_TEXT(fileid,theta_NC,"units",7,"radians")

  dd(1) = n_time_NC
  err = NF_DEF_VAR(fileid,'t',NF_REAL,1,dd,t_NC)
  err = NF_PUT_ATT_TEXT(fileid,t_NC,"long_name",1,"t")
  err = NF_PUT_ATT_TEXT(fileid,t_NC,"units",8,"a/c!ds!n")

  dd(1) = n_n_NC
  err = NF_DEF_VAR(fileid,'n',NF_INT,1,dd,n_NC)
  err = NF_PUT_ATT_TEXT(fileid,n_NC,"long_name",1,"n")
  err = NF_PUT_ATT_TEXT(fileid,n_NC,"units",14,"dimensionless")
  !
  ! 2D variables and their attributes
  !
  if (exists_diff_i == 1) then
     dd(3) = n_time_NC
     dd(2) = n_kinetic_NC
     dd(1) = n_x_NC
     err = NF_DEF_VAR(fileid,'diff_p',NF_REAL,3,dd,diffp_i_NC)
     err = NF_PUT_ATT_TEXT(fileid,diffp_i_NC,"long_name",1,"D")
     err = NF_PUT_ATT_TEXT(fileid,diffp_i_NC,"units",18,"!4q!3!u2!nc!ds!n/a")

     dd(3) = n_time_NC
     dd(2) = n_kinetic_NC
     dd(1) = n_x_NC
     err = NF_DEF_VAR(fileid,'diff_t',NF_REAL,3,dd,difft_i_NC)
     err = NF_PUT_ATT_TEXT(fileid,difft_i_NC,"long_name",5,"!4v!3")
     err = NF_PUT_ATT_TEXT(fileid,difft_i_NC,"units",18,"!4q!3!u2!nc!ds!n/a")

     if (n_field ==2) then
        dd(3) = n_time_NC
        dd(2) = n_kinetic_NC
        dd(1) = n_x_NC
        err = NF_DEF_VAR(fileid,'diff_p_em',NF_REAL,3,dd,diffp_i_em_NC)
        err = NF_PUT_ATT_TEXT(fileid,diffp_i_em_NC,"long_name",1,"D_em")
        err = NF_PUT_ATT_TEXT(fileid,diffp_i_em_NC,"units",18,"!4q!3!u2!nc!ds!n/a")
        
        dd(3) = n_time_NC
        dd(2) = n_kinetic_NC
        dd(1) = n_x_NC
        err = NF_DEF_VAR(fileid,'diff_t_em',NF_REAL,3,dd,difft_i_em_NC)
        err = NF_PUT_ATT_TEXT(fileid,difft_i_em_NC,"long_name",5,"!4v!3")
        err = NF_PUT_ATT_TEXT(fileid,difft_i_em_NC,"units",18,"!4q!3!u2!nc!ds!n/a")
     endif
  endif

  if (exists_u == 1) then
     if (f_flag == 0) then
     dd(4) = n_time_NC
     dd(3) = n_n_NC
     dd(2) = n_x_NC
     dd(1) = n_ri_NC
     err = NF_DEF_VAR(fileid,'potential_r',NF_REAL,4,dd,phi_NC)
     err = NF_PUT_ATT_TEXT(fileid,phi_NC,"long_name",8,"!4u!3(r)")
     err = NF_PUT_ATT_TEXT(fileid,phi_NC,"units",3,"T/e")
     endif
  endif

  if (n_field == 2) then
     if (f_flag == 0) then
     dd(4) = n_time_NC
     dd(3) = n_n_NC
     dd(2) = n_x_NC
     dd(1) = n_ri_NC
     err = NF_DEF_VAR(fileid,'a_parallel',NF_REAL,4,dd,a_parallel_NC)
     err = NF_PUT_ATT_TEXT(fileid,phi_NC,"long_name",8,"!4a!3(r)")
     err = NF_PUT_ATT_TEXT(fileid,phi_NC,"units",3,"Tc/(e c_s)")
     endif
  endif


  if (field_r0_flag == 1) then
     dd(4) = n_time_NC
     dd(3) = n_n_NC
     dd(2) = n_theta_NC
     dd(1) = n_ri_NC
     err = NF_DEF_VAR(fileid,'potential_theta',NF_REAL,4,dd,phi_r0_NC)
     err = NF_PUT_ATT_TEXT(fileid,phi_r0_NC,"long_name",8,"!4u(h)!3")
     err = NF_PUT_ATT_TEXT(fileid,phi_r0_NC,"units",3,"T/e")
  endif

  if (exists_n == 1) then
     if (f_flag == 0) then
     dd(4) = n_time_NC
     dd(3) = n_n_NC
     dd(2) = n_x_NC
     dd(1) = n_ri_NC
     err = NF_DEF_VAR(fileid,'moment_n',NF_REAL,4,dd,moment_n_NC)
     err = NF_PUT_ATT_TEXT(fileid,moment_n_NC,"long_name",6,"!3n(r)")
     err = NF_PUT_ATT_TEXT(fileid,moment_n_NC,"units",14,"n!de!n(r!d0!n)")
     endif
  endif

  if (exists_e == 1) then
     if (f_flag == 0) then
     dd(4) = n_time_NC
     dd(3) = n_n_NC
     dd(2) = n_x_NC
     dd(1) = n_ri_NC
     err = NF_DEF_VAR(fileid,'moment_e',NF_REAL,4,dd,moment_e_NC)
     err = NF_PUT_ATT_TEXT(fileid,moment_e_NC,"long_name",6,"!3E(r)")
     err = NF_PUT_ATT_TEXT(fileid,moment_e_NC,"units",14,"T!de!n(r!d0!n)")
     endif
  endif

  !
  ! 3D variables and their atributes
  !
  if (f_flag == 1) then
     if (exists_u == 1) then
        dd(5) = n_time_NC
        dd(4) = n_n_NC
        dd(3) = n_x_NC
        dd(2) = n_theta_NC
        dd(1) = n_ri_NC
        err = NF_DEF_VAR(fileid,'potential_full',NF_REAL,5,dd,phi_full_NC)
        err = NF_PUT_ATT_TEXT(fileid,phi_NC,"long_name",10,"!4u!3(r.h)")
        err = NF_PUT_ATT_TEXT(fileid,phi_NC,"units",3,"T/e")
        if (n_field == 2) then
           dd(5) = n_time_NC
           dd(4) = n_n_NC
           dd(3) = n_x_NC
           dd(2) = n_theta_NC
           dd(1) = n_ri_NC
           err = NF_DEF_VAR(fileid,'a_parallel_full',NF_REAL,5,dd,a_parallel_full_NC)
           err = NF_PUT_ATT_TEXT(fileid,phi_NC,"long_name",10,"!4a!3(r.h)")
           err = NF_PUT_ATT_TEXT(fileid,phi_NC,"units",3,"T/e")
        endif
     endif
  endif
  
  if (f_flag == 1) then
     if (exists_n ==1) then
        dd(5) = n_time_NC
        dd(4) = n_n_NC
        dd(3) = n_x_NC
        dd(2) = n_theta_NC
        dd(1) = n_ri_NC
        err = NF_DEF_VAR(fileid,'moment_n_full',NF_REAL,5,dd,moment_n_full_NC)
        err = NF_PUT_ATT_TEXT(fileid,moment_n_NC,"long_name",6,"!3n(r)")
        err = NF_PUT_ATT_TEXT(fileid,moment_n_NC,"units",14,"n!de!n(r!d0!n)")
     endif
     if (exists_e ==1) then
        dd(5) = n_time_NC
        dd(4) = n_n_NC
        dd(3) = n_x_NC
        dd(2) = n_theta_NC
        dd(1) = n_ri_NC
        err = NF_DEF_VAR(fileid,'moment_e_full',NF_REAL,5,dd,moment_e_full_NC)
        err = NF_PUT_ATT_TEXT(fileid,moment_n_NC,"long_name",6,"!3E(r)")
        err = NF_PUT_ATT_TEXT(fileid,moment_n_NC,"units",14,"n!de!n(r!d0!n)")
        
     endif
  endif

  err = NF_ENDDEF(fileid)
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Write the data.
  !
  err = NF_PUT_VAR_REAL(fileid,q_NC,q(:))
  err = NF_PUT_VAR_REAL(fileid,r_NC,r(:))
  err = NF_PUT_VAR_REAL(fileid,theta_NC,theta_plot(:))
  err = NF_PUT_VAR_REAL(fileid,t_NC,t(:))
  err = NF_PUT_VAR_INT(fileid,n_NC,n(:))

  if (exists_diff_i == 1) then
     err = NF_PUT_VAR_REAL(fileid,diffp_i_NC,diff_density(:,:,:))
     err = NF_PUT_VAR_REAL(fileid,difft_i_NC,diff_energy(:,:,:))
     if (n_field == 2) then
        err = NF_PUT_VAR_REAL(fileid,diffp_i_em_NC,diff_density_em(:,:,:))
        err = NF_PUT_VAR_REAL(fileid,difft_i_em_NC,diff_energy_em(:,:,:))
     endif
  endif

  if (exists_u == 1) then     
     if (f_flag == 0) then
     err = NF_PUT_VAR_REAL(fileid,phi_NC,phi0(:,:,:,:))     
     endif
  endif

  if (n_field == 2) then
     if (f_flag == 0) then 
     err = NF_PUT_VAR_REAL(fileid,a_parallel_NC,a_parallel0(:,:,:,:))
     endif
  endif

  if (field_r0_flag == 1) then
     err = NF_PUT_VAR_REAL(fileid,phi_r0_NC,phi_r0(:,:,:,:))
  endif

  if (exists_n ==1 ) then
     if (f_flag == 0) then 
     err = NF_PUT_VAR_REAL(fileid,moment_n_NC,moment_n(:,:,:,:))
     endif
  endif

  if (exists_e == 1) then
     if (f_flag == 0)then 
        err = NF_PUT_VAR_REAL(fileid,moment_e_NC,moment_e(:,:,:,:))
     endif
  endif

  if (f_flag == 1) then ! saving full data
     err = NF_PUT_VAR_REAL(fileid,phi_full_NC,phi0_full(:,:,:,:,:))
     if (n_field == 2 )then ! a_parallel
        err = NF_PUT_VAR_REAL(fileid,a_parallel_full_NC,a_parallel0_full(:,:,:,:,:))
     endif
     if (exists_n==1) then ! ion density
        err = NF_PUT_VAR_REAL(fileid,moment_n_full_NC,moment_n_full(:,:,:,:,:))
     endif
     if (exists_e==1) then ! electron density
        err = NF_PUT_VAR_REAL(fileid,moment_e_full_NC,moment_e_full(:,:,:,:,:))
     endif
  endif
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Close the netCDF file   
  !
  err = NF_CLOSE(fileid)
  !------------------------------------------------------------------

end program main
