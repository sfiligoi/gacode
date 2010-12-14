subroutine read_input(dir)

  use fieldline_input_data

  implicit none

  character (len=*) :: dir

  integer :: i 
  integer :: j
  integer :: k ! WG
  integer :: i_n
  integer :: i_t
  integer :: i_field

  real :: dum
  real :: x_real,x_imag
  real, dimension(:), allocatable :: vec
  real, dimension(:,:), allocatable :: svec

  open(unit=1,file=trim(dir)//'profile_vugyro.out')
  read(1,*) n_x
  read(1,*) n_theta_section
  read(1,*) n_pass
  read(1,*) n_trap
  read(1,*) n_energy
  read(1,*) n_theta_plot
  read(1,*) n0
  read(1,*) n_n
  read(1,*) d_n
  read(1,*) n_explicit_damp
  read(1,*) nonlinear_flag
  read(1,*) electron_method
  read(1,*) n_field
  read(1,*) n_ion
  read(1,*) n_kinetic
  read(1,*) n_spec
  read(1,*) field_r0_flag
  read(1,*) field_r0_grid
  read(1,*) n_grid_exp
  read(1,*) boundary_method

  allocate(r(n_x))
  allocate(r_s(n_x))
  allocate(q(n_x))
  allocate(shat_s(n_x))
  allocate(aspect_s(n_x))
  allocate(vec(n_x))
  allocate(svec(n_spec,n_x))
  allocate(ky(0:n_n-1))

  read(1,*) r(:)
  read(1,*) q(:)
  read(1,*) r_s(:)
  read(1,*) vec(:)    ! q_s
  read(1,*) svec(:,:) ! dlntdr_s
  read(1,*) svec(:,:) ! dlnndr_s
  read(1,*) svec(:,:) ! tem_s
  read(1,*) svec(:,:) ! den_s
  read(1,*) aspect_s(:) 
  read(1,*) vec(:)    ! delta_s
  read(1,*) vec(:)    ! zeta_s
  read(1,*) vec(:)    ! kappa_s
  read(1,*) vec(:)    ! drmaj_s
  read(1,*) shat_s(:) 
  read(1,*) vec(:)    ! s_delta_s
  read(1,*) vec(:)    ! s_zeta_s
  read(1,*) vec(:)    ! s_kappa_s
  read(1,*) vec(:)    ! zmag_s
  read(1,*) vec(:)    ! dzmag_s
  read(1,*) vec(:)    ! beta_unit_s
  read(1,*) vec(:)    ! gamma_e_s 
  read(1,*) vec(:)    ! gamma_p_s 
  read(1,*) vec(:)    ! mach_s 
  read(1,*) vec(:)    ! b_unit_s
  read(1,*) vec(:)    ! dr_eodr
  read(1,*) vec(:)    ! z_eff_s
  read(1,*) vec(:)    ! nu_s
  read(1,*) vec(:)    ! w0_s

  read(1,*) dum ! box_multiplier
  do i=1,n_pass+n_trap+n_energy+1
     read(1,*) dum
  enddo
  read(1,*) ky
  read(1,*) rhos_norm
  close(1)

  deallocate(vec)

  !=====================================================
  ! Read G_theta from geometry_arrays.out
  open(unit=1,file=trim(dir)//'geometry_arrays.out') ! WG

  allocate(g_theta(n_theta_plot+1,-n_x/2:n_x/2-1))

  do i=-n_x/2,n_x/2-1
     do j=1,n_theta_plot
        do k=1,11
           read(1,*) dum
           if (k == 8) g_theta(j,i) = dum
        enddo
     enddo
  enddo
  close(1)

  !=====================================================

  open(unit=1,file=trim(dir)//'u.out')

  allocate(a(n_theta_plot+1,-n_x/2:n_x/2-1,0:n_n-1))

  do i_t=0,n_time_slice
     do i_n=0,n_n-1
        do i_field=1,n_field
           do i=-n_x/2,n_x/2-1
              do j=1,n_theta_plot
                 read(1,*,IOSTAT=ierr) x_real,x_imag
                 if (ierr /= 0) then
                    n_time_slice = i_t
                    goto 100
                 endif
                 if (i_field == 2) a(j,i,i_n) = cmplx(x_real,x_imag)
              enddo
           enddo
        enddo
     enddo
  enddo

100 return

end subroutine read_input
