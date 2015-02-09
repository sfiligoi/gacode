!-----------------------------------------------------------------
! cgyro_write_timedata.f90
!
! PURPOSE:
!  Output of time-independent (initial) data.
!-----------------------------------------------------------------

subroutine cgyro_write_initdata

  use mpi
  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer :: in,is, it

  !----------------------------------------------------------------------------
  ! Runfile to give complete summary to user
  ! 
  ! NOTE: Some of this data is reproduced in other data files.
  !
  if (silent_flag == 0 .and. i_proc == 0) then

     open(unit=io,file=trim(path)//runfile_info,status='old',position='append')

     write(io,*)
     write(io,'(a)') 'n_radial  n_theta   n_species n_energy  n_xi'
     write(io,'(t1,5(i4,6x))') n_radial,n_theta,n_species,n_energy,n_xi

     write(io,*) 
     write(io,'(a,f7.2,a,f7.2)') '(Lx,Ly)/rho',length/rho,',',2*pi/ky
     write(io,*) 

     if (zf_test_flag == 1) then
        write(io,20) 'ky*rho',0.0
     else 
        write(io,20) 'ky*rho'
        write(io,'(8f6.3,1x)') (in*q/rmin*rho,in=0,n_toroidal-1)
     endif

     write(io,*) 

     if (zf_test_flag == 1) then
        write(io,20) 'kx*rho',2*pi*rho/length
     else
        write(io,20) 'min(kx*rho)',2*pi*rho/length
        write(io,20) 'max(kx*rho)',2*pi*rho*(n_radial/2-1)/length
     endif

     write(io,*) 
     write(io,20) 'rho',rho
     write(io,20) 'r/a',rmin
     write(io,20) 'R/a',rmaj
     write(io,20) 'q',q
     write(io,20) 's',s
     write(io,20) 'shift',shift
     write(io,20) 'kappa',kappa
     write(io,20) 's_kappa',s_kappa
     write(io,20) 'delta',delta
     write(io,20) 's_delta',s_delta
     write(io,20) 'zeta',zeta
     write(io,20) 's_zeta',s_zeta
     write(io,20) 'zmag',zmag
     write(io,20) 's_zmag',s_zmag

     write(io,*)
     write(io,'(a)') &
          'indx  z    n/n_norm    T/T_norm    m/m_norm     a/Ln        a/Lt        nu'
     do is=1,n_species
        write(io,'(t2,i2,2x,i2,2x,6(1pe10.4,2x))') &
             is,z(is),dens(is),temp(is),mass(is),dlnndr(is),dlntdr(is),nu(is)
     enddo

     write(io,*)

     write(io,*) 'nc',nc
     write(io,*) 'nv',nv

     close(io)

  endif
  !----------------------------------------------------------------------------

  if (test_flag == 1) return

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
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Write the initial equilibrium geometry data
  !
  if (silent_flag == 0 .and. i_proc == 0) then

     open(unit=io,file=trim(path)//'out.cgyro.geo',status='replace')
     do it=1,n_theta
        write(io,fmtstr) theta(it)
     enddo
     do it=1,n_theta
        write(io,fmtstr) w_theta(it)
     enddo
     do it=1,n_theta
        write(io,fmtstr) Bmag(it)
     enddo
     do it=1,n_theta
        write(io,fmtstr) omega_stream(it,1)
     enddo
     do it=1,n_theta
        write(io,fmtstr) omega_trap(it,1)
     enddo
     do it=1,n_theta
        write(io,fmtstr) omega_rdrift(it,1)
     enddo
     do it=1,n_theta
        write(io,fmtstr) omega_adrift(it,1)
     enddo
     do it=1,n_theta
        write(io,fmtstr) omega_aprdrift(it,1)
     enddo
     do it=1,n_theta
        write(io,fmtstr) k_perp(it,n_radial/2+1)
     enddo
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
     write(io,'(i4)') px(:)
     write(io,'(1pe12.5)') theta(:)
     write(io,'(1pe12.5)') energy(:)
     write(io,'(1pe12.5)') xi(:)
     write(io,'(1pe12.5)') transpose(thetab(:,:))
     write(io,'(1pe12.5)') (rho*q/rmin*in,in=0,n_toroidal-1)
     close(io)

  endif
  !----------------------------------------------------------------------------

20 format(t2,a,':',t13,1pe12.4) 

end subroutine cgyro_write_initdata
