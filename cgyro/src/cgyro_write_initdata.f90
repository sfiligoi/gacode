!-----------------------------------------------------------------
! cgyro_write_timedata.f90
!
! PURPOSE:
!  Output of time-independent (initial) data.
!-----------------------------------------------------------------

subroutine cgyro_write_initdata

  use mpi
  use cgyro_globals
  use cgyro_experimental_globals
  
  implicit none

  integer :: in,is
  real :: kymax,z_eff

  !----------------------------------------------------------------------------
  ! Runfile to give complete summary to user
  ! 
  ! NOTE: Some of this data is reproduced in other data files.
  !
  if (silent_flag == 0 .and. i_proc == 0) then

     open(unit=io,file=trim(path)//runfile_info,status='old',position='append')

  
     write(io,*)
     write(io,'(a)') ' n_theta | n_species | n_energy | n_xi '
     write(io,'(t4,i3,t16,i1,t26,i2,t36,i2)') n_theta,n_species,n_energy,n_xi
     if (test_flag == 0) then
        write(io,*) 
        write(io,'(a)') ' nc_loc | nv_loc | nsplit | n_MPI | n_OMP'
        write(io,'(t3,i4,t12,i4,t21,i4,t29,i4,t36,i4)') nc_loc,nv_loc,nsplit,n_proc,n_omp
     endif

     if (zf_test_flag == 0) then

        ! Compute kymax
        if (n_toroidal == 1) then
           kymax = q/rmin*rho
        else
           kymax = q/rmin*(n_toroidal-1)*rho
        endif

        write(io,*)
        write(io,*) '          n    Delta      Max     L/rho'
        write(io,'(a,i4,2x,2(f7.3,2x),2x,f6.2)') ' kx*rho:',&
             n_radial,2*pi*rho/length,2*pi*rho*(n_radial/2-1)/length,length/rho
        write(io,'(a,i4,2x,2(f7.3,2x),2x,f6.2)') ' ky*rho:',&
             n_toroidal,q/rmin*rho,kymax,2*pi/ky

     else

        write(io,*) ' kx*rho:',2*pi*rho/length

     endif

     write(io,*) 
     write(io,20) '    r/a:',rmin
     write(io,20) '    R/a:',rmaj, '  shift:',shift
     write(io,20) '      q:',q,    '      s:',s
     write(io,20) '  kappa:',kappa,'s_kappa:',s_kappa
     write(io,20) '  delta:',delta,'s_delta:',s_delta
     write(io,20) '   zeta:',zeta, ' s_zeta:',s_zeta
     write(io,20) '   zmag:',zmag, '  dzmag:',dzmag
     write(io,*)
     write(io,20) '  betae:',betae_unit, ' beta_*:',beta_star,'lambda_*:',lambda_star
     write(io,*)
     write(io,20) 'gamma_e:', gamma_e,   'gamma_p:', gamma_p, '    mach:', mach

     z_eff = 0.0
     do is=1,n_species
        if (z(is) > 0.0) then 
           z_eff = z_eff+dens(is)*z(is)**2/dens_ele
        endif
     enddo
     write(io,*)
     write(io,20) '    rho:',rho,'  z_eff:',z_eff

     write(io,*)
     write(io,'(a)') &
          'indx  z     n/n_norm     T/T_norm     m/m_norm      a/Ln         a/Lt         nu'
     do is=1,n_species
        write(io,'(t2,i2,2x,i2,2x,6(1pe11.4,2x))') &
             is,z(is),dens(is),temp(is),mass(is),dlnndr(is),dlntdr(is),nu(is)
     enddo

     if (profile_model == 2) then
        write(io,*)
        write(io,20) ' a_meters:',a_meters, '   b_unit:',b_unit
        write(io,20) 'dens_norm:',dens_norm,'temp_norm:',temp_norm
        write(io,20) ' vth_norm:',vth_norm
     endif

     write(io,*)

     close(io)

  endif
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  ! Write the initial equilibrium data
  !
  if (silent_flag == 0 .and. i_proc == 0) then

     open(unit=io,file=trim(path)//'out.cgyro.equil',status='replace')
     write (io,fmtstr,advance='no') rmin
     write (io,fmtstr,advance='no') rmaj
     write (io,fmtstr,advance='no') q
     write (io,fmtstr,advance='no') s
     write (io,fmtstr,advance='no') rho
     write (io,fmtstr,advance='no') ky
     do is=1,n_species
        write (io,fmtstr,advance='no') dens(is)
        write (io,fmtstr,advance='no') temp(is)
        write (io,fmtstr,advance='no') dlnndr(is)
        write (io,fmtstr,advance='no') dlntdr(is)
        write (io,fmtstr,advance='no') nu(is)
     enddo
     write (io,*)
     close(io)

  endif

  if (silent_flag == 0 .and. i_proc == 0 .and. profile_model == 2) then

     open(unit=io,file=trim(path)//'out.cgyro.expnorm',status='replace')
     write (io,fmtstr,advance='no') a_meters
     write (io,fmtstr,advance='no') b_unit
     write (io,fmtstr,advance='no') dens_norm
     write (io,fmtstr,advance='no') temp_norm
     write (io,fmtstr,advance='no') vth_norm
     write (io,*)
     close(io)

  endif

  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Write the initial equilibrium geometry data
  !
  if (silent_flag == 0 .and. i_proc == 0) then

     open(unit=io,file=trim(path)//'out.cgyro.geo',status='replace')
     write(io,fmtstr) theta(:)
     write(io,fmtstr) w_theta(:)
     write(io,fmtstr) bmag(:)
     write(io,fmtstr) omega_stream(:,1)
     write(io,fmtstr) omega_trap(:,1)
     write(io,fmtstr) omega_rdrift(:,1)
     write(io,fmtstr) omega_adrift(:,1)
     write(io,fmtstr) omega_aprdrift(:,1)
     write(io,fmtstr) omega_cdrift(:,1)
     write(io,fmtstr) omega_gammap(:)
     write(io,fmtstr) k_perp(ic_c(n_radial/2+1,:))
     close(io)

  endif
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Write the initial grid data 
  !
  if (silent_flag == 0 .and. i_proc == 0) then

     open(unit=io,file=trim(path)//runfile_grids,status='replace')
     write(io,'(i4)') n_toroidal
     write(io,'(i4)') n_species
     write(io,'(i4)') n_field
     write(io,'(i4)') n_radial
     write(io,'(i4)') n_theta
     write(io,'(i4)') n_energy
     write(io,'(i4)') n_xi
     write(io,'(i4)') box_size
     write(io,'(1pe12.5)') length/rho
     write(io,'(i4)') px(:)
     write(io,'(1pe12.5)') theta(:)
     write(io,'(1pe12.5)') energy(:)
     write(io,'(1pe12.5)') xi(:)
     write(io,'(1pe12.5)') transpose(thetab(:,:))
     write(io,'(1pe12.5)') (rho*q/rmin*in,in=0,n_toroidal-1)
     close(io)

  endif
  !----------------------------------------------------------------------------

20 format(t2,3(a,1x,1pe11.4,4x)) 

end subroutine cgyro_write_initdata
