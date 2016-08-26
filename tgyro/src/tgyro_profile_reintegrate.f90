subroutine tgyro_profile_reintegrate

  use tgyro_globals
  use tgyro_ped

  use mpi
  use tgyro_iteration_variables
  use EXPRO_interface

  implicit none

  integer :: i_ion
  integer :: i_star  

  CHARACTER(LEN=6) :: NUM

  if (tgyro_ped_model > 1) then

     ! Map data over r(n_r) < r < a

     call tgyro_pedestal_map(dlnnedr(n_r),zn_top,n_top(1),nn_vec(:,2),i_star,exp_ne)
     call tgyro_pedestal_map(dlntedr(n_r),zt_top,t_top(1),t_vec(:),i_star,exp_te)
     call tgyro_pedestal_map(dlntidr(1,n_r),zt_top,t_top(1),t_vec(:),i_star,exp_ti(1,:))

     ! Set ion densities
     exp_ni(1,i_star:n_exp) = exp_ne(i_star:n_exp)
     if (loc_n_ion > 1) then
        exp_ni(1,i_star:n_exp) = exp_ni(1,i_star:n_exp)-zi_vec(2)*exp_ni(2,i_star:n_exp)
     endif
     if (loc_n_ion > 2) then
        exp_ni(1,i_star:n_exp) = exp_ni(1,i_star:n_exp)-zi_vec(3)*exp_ni(3,i_star:n_exp)
     endif
     ! Set thermal ion temperatures
     do i_ion=2,loc_n_ion
        if (therm_flag(i_ion) == 1) exp_ti(i_ion,i_star:n_exp) = exp_ti(1,i_star:n_exp)
     enddo
  endif

  ! Map data inside r < r(n_r)

  call tgyro_expro_map(r,dlnnedr,n_r,rmin_exp,exp_ne,n_exp)
  call tgyro_expro_map(r,dlntedr,n_r,rmin_exp,exp_te,n_exp)
  do i_ion=1,loc_n_ion
     if (therm_flag(i_ion) == 1) then
        call tgyro_expro_map(r,dlnnidr(i_ion,:),n_r,rmin_exp,exp_ni(i_ion,:),n_exp)
        call tgyro_expro_map(r,dlntidr(i_ion,:),n_r,rmin_exp,exp_ti(i_ion,:),n_exp)
     endif
  enddo
  
  ptot_exp = exp_ne*exp_te
  do i_ion=1,loc_n_ion
     ptot_exp = ptot_exp + exp_ni(i_ion,:)*exp_ti(i_ion,:)
  enddo

  ! Convert to Pa: n[1/cm^3]*(kT[ev])/10  
  ptot_exp = ptot_exp*k/10.0

  if ((tgyro_write_profiles_flag==-1) .and. (i_tran.ge.2) .and. (i_tran_old .ne. i_tran)) then
     call EXPRO_palloc(MPI_COMM_WORLD,'./',1)
     call EXPRO_pread

     EXPRO_ptot = ptot_exp
     EXPRO_ne   = exp_ne*1e-13
     EXPRO_te   = exp_te*1e-3
     EXPRO_ni(1:loc_n_ion,:) = exp_ni(1:loc_n_ion,:)*1e-13
     EXPRO_ti(1:loc_n_ion,:) = exp_ti(1:loc_n_ion,:)*1e-3
     EXPRO_ptot = ptot_exp ! already in Pa

     if (i_proc_global == 0) then
        ! Write data to file
        write(NUM,'(i0)')i_tran-1
        write(*,*)'write input.profiles.'//trim(NUM)
        call EXPRO_write_original(&
             1,'input.profiles',&
             2,'input.profiles.'//trim(NUM),&
             'Profiles modified by TGYRO')
     endif

     call EXPRO_palloc(MPI_COMM_WORLD,'./',0)
     i_tran_old=i_tran
  endif

end subroutine tgyro_profile_reintegrate
