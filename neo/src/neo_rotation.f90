module neo_rotation

  implicit none

  public :: ROT_alloc, ROT_solve_phi, ROT_write
  real, dimension(:), allocatable   :: phi_rot        ! phi(r)-phi(r,theta0)
  real, dimension(:), allocatable   :: phi_rot_deriv  ! theta derivative
  real, dimension(:), allocatable   :: phi_rot_rderiv ! radial derivative
  real                              :: phi_rot_avg    ! flux surface avg
  real, dimension(:,:), allocatable :: dens_fac       ! (ns,nth):
                                                      ! n=n(th0)*dens_fac
  real, dimension(:), allocatable   :: dens_avg       ! (ns): <n>/n0
  real, dimension(:), allocatable   :: dens_avg_cos
  real, dimension(:), allocatable   :: rotavg_e0, rotavg_e1, &
       rotavg_e2, rotavg_e3, rotavg_e4

  logical, private :: initialized = .false.
  integer, private :: io_rot=50
  character(len=80),private :: runfile = 'out.neo.rotation'

  contains

    subroutine ROT_alloc(flag)
      use neo_globals
      implicit none
      integer, intent (in) :: flag  ! flag=1: allocate; else deallocate

      if(flag == 1) then
         if(initialized) return
         allocate(phi_rot(n_theta))
         allocate(phi_rot_deriv(n_theta))
         allocate(phi_rot_rderiv(n_theta))
         allocate(dens_fac(n_species,n_theta))
         allocate(dens_avg(n_species))
         allocate(dens_avg_cos(n_species))
         allocate(rotavg_e0(n_species))
         allocate(rotavg_e1(n_species))
         allocate(rotavg_e2(n_species))
         allocate(rotavg_e3(n_species))
         allocate(rotavg_e4(n_species))
         if(silent_flag == 0 .and. i_proc == 0 .and. rotation_model == 2) then
            open(unit=io_rot,file=trim(path)//runfile,status='replace')
            close(io_rot)
         end if
         initialized = .true.

      else
         if(.NOT. initialized) return
         deallocate(phi_rot)
         deallocate(phi_rot_deriv)
         deallocate(phi_rot_rderiv)
         deallocate(dens_fac)
         deallocate(dens_avg)
         deallocate(dens_avg_cos)
         deallocate(rotavg_e0)
         deallocate(rotavg_e1)
         deallocate(rotavg_e2)
         deallocate(rotavg_e3)
         deallocate(rotavg_e4)
         initialized = .false.

      endif

    end subroutine ROT_alloc

    ! solve the quasi-neutrality relation for poloidal part of phi
    subroutine ROT_solve_phi(ir)
      use neo_globals
      use neo_equilibrium
      implicit none
      integer, intent (in) :: ir
      integer, parameter :: nmax = 200
      integer :: it, is, jt, id, n
      real :: x, x0, sum_zn, dsum_zn, fac
      
      if(rotation_model == 1 .or. spitzer_model==1) then
         do it=1, n_theta
            phi_rot(it) = 0.0
            phi_rot_deriv(it) = 0.0
            phi_rot_rderiv(it) = 0.0
            do is=1, n_species
               dens_fac(is,it) = 1.0
            enddo
         enddo
         phi_rot_avg = 0.0
         omega_rot(ir) = 0.0
         omega_rot_deriv(ir) = 0.0
         do is=1, n_species
            dens_avg(is) = dens(is,ir)
         enddo
         rotavg_e0(:) = 0.0
         rotavg_e1(:) = 0.0
         rotavg_e2(:) = 0.0
         rotavg_e3(:) = 0.0
         rotavg_e4(:) = 0.0

      else 
         
         phi_rot_avg = 0.0
         dens_avg    = 0.0
         dens_avg_cos = 0.0
         x = 0.05          ! initial guess for phi_rot(1)
         do it=1,n_theta
            
            n=1

            do 
               ! use Newton's method to solve the quasi-neutrality relation
               
               sum_zn  = 0.0
               dsum_zn = 0.0
               do is=1, n_species
                  fac = z(is) * dens(is,ir) * &
                       exp(  omega_rot(ir)**2 * 0.5 / vth(is,ir)**2 &
                       * (bigR(it)**2 - bigR_th0**2) &
                       - z(is) / temp(is,ir) * x)
                  sum_zn  = sum_zn  + fac
                  dsum_zn = dsum_zn - z(is) / temp(is,ir) * fac
               enddo
               
               if(adiabatic_ele_model == 1) then
                  fac = -ne_ade(ir) * exp(1.0/te_ade(ir) * x)
                  sum_zn  = sum_zn  + fac
                  dsum_zn = dsum_zn + 1.0 / te_ade(ir) * fac
               endif

               x0 = x
               x  = x0 - sum_zn / dsum_zn
               
               if (abs(x-x0) < 1.0e-12) exit
               
               n = n + 1
               if(n > nmax) exit
               
            enddo
            
            if(n > nmax) then
               call neo_error('ERROR: (NEO) Rotation density computation failed to converge')
               return
            endif
            
            ! phi_rot is total phi(r) - total phi(r,theta0)
            phi_rot(it) = x  
            
            ! flux surface avg
            phi_rot_avg = phi_rot_avg &
                 + w_theta(it) * phi_rot(it)  
            
            do is=1, n_species
               dens_fac(is,it) = exp(omega_rot(ir)**2 * 0.5/vth(is,ir)**2 &
                       * (bigR(it)**2 - bigR_th0**2) &
                       - z(is) / temp(is,ir) * phi_rot(it))
               dens_avg(is) = dens_avg(is) + &
                    w_theta(it) * dens(is,ir) * dens_fac(is,it)
               dens_avg_cos(is) = dens_avg_cos(is) + cos(theta(it)) * &
                    w_theta(it) * dens(is,ir) * dens_fac(is,it)
            enddo
            
         enddo

         ! solve the radial derivative of QN for d phi_rot/dr
         do it=1,n_theta
            phi_rot_rderiv(it) = 0.0
            sum_zn  = 0.0
            do is=1, n_species
               fac = z(is) * dens(is,ir) &
                    * exp(omega_rot(ir)**2 * 0.5 / vth(is,ir)**2 &
                    * (bigR(it)**2 - bigR_th0**2) &
                    - z(is) / temp(is,ir) * phi_rot(it))

               sum_zn = sum_zn + z(is) / temp(is,ir) * fac
            enddo
            if(adiabatic_ele_model == 1) then
               fac = -ne_ade(ir) * exp(1.0/te_ade(ir) * phi_rot(it))
               sum_zn = sum_zn - 1.0/te_ade(ir) * fac
            endif

            do is=1,n_species
               fac = z(is) * dens(is,ir) &
                    * exp(omega_rot(ir)**2 * 0.5 / vth(is,ir)**2 &
                    * (bigR(it)**2 - bigR_th0**2) &
                    - z(is) / temp(is,ir) * phi_rot(it))
               phi_rot_rderiv(it) = phi_rot_rderiv(it) &
                    + fac * (-dlnndr(is,ir) + omega_rot_deriv(ir) &
                    * omega_rot(ir)/vth(is,ir)**2 &
                    * (bigR(it)**2 - bigR_th0**2) &
                    + omega_rot(ir)**2 / vth(is,ir)**2 &
                    * (bigR(it)*bigR_rderiv(it) &
                    - bigR_th0*bigR_th0_rderiv) &
                    - phi_rot(it) * z(is)/temp(is,ir) * dlntdr(is,ir) &
                    + 0.5 * omega_rot(ir)**2 / vth(is,ir)**2 * dlntdr(is,ir) &
                    * (bigR(it)**2 - bigR_th0**2) ) 
            enddo
            if(adiabatic_ele_model == 1) then
               fac = -ne_ade(ir) * exp(1.0/te_ade(ir) * phi_rot(it))
               phi_rot_rderiv(it) = phi_rot_rderiv(it) &
                    + fac * (-dlnndre_ade(ir) &
                    + phi_rot(it)/te_ade(ir) * dlntdre_ade(ir))
            endif
            
            phi_rot_rderiv(it) = phi_rot_rderiv(it) / sum_zn
         enddo

         do it=1,n_theta
            phi_rot_deriv(it) = 0.0
            do id=-2,2
               if (id /= 0) then
                  jt = thcyc(it+id)
                  phi_rot_deriv(it) = phi_rot_deriv(it) &
                       + phi_rot(jt) * cderiv(id) / (12.0*d_theta)
               endif
            enddo
         enddo

         do is=1,n_species
            rotavg_e0(is) = 0.0
            rotavg_e1(is) = 0.0
            rotavg_e2(is) = 0.0
            rotavg_e3(is) = 0.0
            rotavg_e4(is) = 0.0
            do it=1,n_theta
               fac = exp(omega_rot(ir)**2 * 0.5 / vth(is,ir)**2 &
                    * (bigR(it)**2 - bigR_th0**2) &
                    - z(is) / temp(is,ir) * phi_rot(it))

               rotavg_e0(is) = rotavg_e0(is) + w_theta(it) * fac
               rotavg_e1(is) = rotavg_e1(is) + w_theta(it) * fac &
                    * z(is)/temp(is,ir) * phi_rot(it)
               rotavg_e2(is) = rotavg_e2(is) + w_theta(it) * fac &
                    * z(is)/temp(is,ir) * phi_rot_rderiv(it)
               rotavg_e3(is) = rotavg_e3(is) + w_theta(it) * fac &
                    * (bigR(it)**2 - bigR_th0**2)
               rotavg_e4(is) = rotavg_e4(is) + w_theta(it) * fac &
                    * 2.0*(bigR(it)*bigR_rderiv(it) &
                    - bigR_th0*bigR_th0_rderiv)
            enddo
         enddo

         !!!!!!!!!!!!!!!!!!!!!!!!!
         ! checks for pure plasma
         !!!!!!!!!!!!!!!!!!!!!!!!!
         !do it=1,n_theta
         !   if(adiabatic_ele_model == 1) then
         !      print *, theta(it), phi_rot(it), &
         !           0.5*omega_rot(ir)**2 * (bigR(it)**2 - bigR_th0**2) &
         !           / (1.0/temp(1,ir) + 1.0/te_ade(ir))
         !   else
         !      print *, theta(it), phi_rot(it), &
         !           0.5*omega_rot(ir)**2 * (bigR(it)**2 - bigR_th0**2) &
         !           / (1.0/temp(1,ir) + 1.0/temp(2,ir))
         !   endif
         !enddo
         !do it=1,n_theta
         !   if(adiabatic_ele_model == 1) then
         !      print *, theta(it), &
         !           phi_rot_rderiv(it)*(1.0/temp(1,ir) + 1.0/te_ade(ir)), &
         !           -phi_rot(it)*(1.0/temp(1,ir)*dlntdr(1,ir) &
         !           + 1.0/te_ade(ir)*dlntdre_ade(ir)) &
         !           + omega_rot(ir) &
         !           * omega_rot_deriv(ir)*(bigR(it)**2 - bigR_th0**2) &
         !           + 0.5*omega_rot(ir)**2  * dlntdr(1,ir) &
         !           *(bigR(it)**2 - bigR_th0**2) &
         !           + omega_rot(ir)**2 *(bigR(it)*bigR_rderiv(it) &
         !           - bigR_th0*bigR_th0_rderiv)
         !   else
         !      print *, theta(it), &
         !           phi_rot_rderiv(it)*(1.0/temp(1,ir) + 1.0/temp(2,ir)), &
         !           -phi_rot(it)*(1.0/temp(1,ir)*dlntdr(1,ir) &
         !           + 1.0/temp(2,ir)*dlntdr(2,ir)) &
         !           + omega_rot(ir) &
         !           * omega_rot_deriv(ir)*(bigR(it)**2 - bigR_th0**2) &
         !           + 0.5*omega_rot(ir)**2  * dlntdr(1,ir) &
         !           *(bigR(it)**2 - bigR_th0**2) &
         !           + omega_rot(ir)**2 *(bigR(it)*bigR_rderiv(it) &
         !           - bigR_th0*bigR_th0_rderiv)
         !   endif
         !enddo

      endif   

    end subroutine ROT_solve_phi

    subroutine ROT_write(ir)
      use neo_globals
      implicit none
      integer, intent (in) :: ir
      integer :: is, it
      real :: gamma_theory

      if(silent_flag == 0 .and. i_proc == 0 .and. rotation_model == 2) then
         open(io_rot,file=trim(path)//runfile,status='old',position='append')
         write (io_rot,'(e16.8)',advance='no') r(ir)
         write (io_rot,'(e16.8)',advance='no') phi_rot_avg
         do is=1, n_species
            write (io_rot,'(e16.8)',advance='no') 1.0/rotavg_e0(is)
            write (io_rot,'(e16.8)',advance='no') ( -rotavg_e2(is) &
                 + rotavg_e3(is) * omega_rot(ir) / vth(is,ir)**2 &
                 * omega_rot_deriv(ir) &
                 + rotavg_e4(is) * 0.5 * omega_rot(ir)**2 / vth(is,ir)**2 &
                 + rotavg_e1(is) * (-dlntdr(is,ir)) &
                 - rotavg_e3(is) * 0.5 * omega_rot(ir)**2 / vth(is,ir)**2 &
                 * (-dlntdr(is,ir)) ) / rotavg_e0(is)
         enddo
         do it=1, n_theta
            write (io_rot,'(e16.8)',advance='no') phi_rot(it)
         enddo
         write(io_rot,*)
         close(io_rot)
      endif

      !call ROT_theory(ir,gamma_theory)
      !open(io_rot,file='out.neo.rot_theory',status='replace')
      !write (io_rot,'(e16.8)',advance='no') gamma_theory
      !write(io_rot,*)
      !close(io_rot)

    end subroutine ROT_write

    subroutine ROT_theory(ir,gamma_theory)
      use neo_globals
      use neo_equilibrium
      implicit none
      integer, intent (in) :: ir
      real, intent (inout) :: gamma_theory
      integer :: is_ion, is_imp
      real    :: geofac
      integer :: it
      
      ! assume primary ion species is is=1
      is_ion = 1
      is_imp = 2

      geofac = 0.0
      do it=1,n_theta
         geofac = geofac + w_theta(it)*dens_fac(is_imp,it)/bmag(it)**2 &
              * (1.0-bmag(it)**2/Bmag2_avg*(1.0-ftrap))
         print *, dens_fac(is_imp,it)/bmag(it)**2
      enddo

      gamma_theory = dens(is_ion,ir)*temp(is_ion,ir) * rho(ir)*rho(ir) &
           * (4.0/3.0/sqrt(pi)) * (dlnndr(is_ion,ir) - 0.5*dlntdr(is_ion,ir)) &
           * nu(is_ion,ir)*dens(is_imp,ir)*z(is_imp)**2/dens(is_ion,ir) &
           * geofac * I_div_psip**2

    end subroutine ROT_theory


end module neo_rotation
