module neo_rotation

  implicit none

  public :: ROT_alloc, ROT_solve_phi, ROT_write
  real, dimension(:), allocatable   :: phi_rot        ! phi(r)-phi(r,theta0)
  real, dimension(:), allocatable   :: phi_rot_deriv  ! theta derivative
  real, dimension(:), allocatable   :: phi_rot_rderiv ! radial derivative
  real                              :: phi_rot_avg    ! flux surface avg
  real, dimension(:,:), allocatable :: dens_fac       ! (ns,nth):
                                                      ! n=n(th0)*dens_fac
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
      integer :: it, is, jt, id, n, is_ele
      real :: x, x0, sum_zn, dsum_zn, fac, fac_add
      real, dimension(:,:), allocatable :: fac_aniso, dlnndr_fac, dens_fac_new
      real, dimension(:), allocatable :: dens_avg, dlnndr_avg, dens_avg_new
      real :: dens_new
      integer :: test_flag

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
         rotavg_e0(:) = 0.0
         rotavg_e1(:) = 0.0
         rotavg_e2(:) = 0.0
         rotavg_e3(:) = 0.0
         rotavg_e4(:) = 0.0

      else 
         
         ! check for equilibrium-scale QN at theta=0
         sum_zn  = 0.0
         dsum_zn = 0.0
         do is=1, n_species
            sum_zn  = sum_zn  + z(is) * dens(is,ir)
            dsum_zn = dsum_zn + z(is) * dens(is,ir) * dlnndr(is,ir)
         enddo
         if(adiabatic_ele_model == 1) then
            sum_zn  = sum_zn - ne_ade(ir)
            sum_zn  = sum_zn / ne_ade(ir)
            dsum_zn = dsum_zn - ne_ade(ir) * dlnndre_ade(ir)
            dsum_zn = dsum_zn / ne_ade(ir)
         else
            do is=1, n_species
               if(Z(is) == -1) then
                  is_ele = is
                  exit
               endif
            enddo
            sum_zn   = sum_zn / dens(is_ele,ir)
            dsum_zn  = dsum_zn / dens(is_ele,ir)
         endif

         if(abs(sum_zn) > 1.0e-3) then
            if(silent_flag == 0 .and. i_proc == 0) then
               open(unit=io_neoout,file=trim(path)//runfile_neoout,&
                    status='old',position='append')
               write(io_neoout,*) &
                    'WARNING: (NEO) ROTATION IS BEING RUN WITHOUT QN DENSITIES'
               close(io_neoout)
            endif
         endif

         allocate(fac_aniso(n_species,n_theta))
         do is=1,n_species
            if(aniso_model(is) == 2) then
               do it=1,n_theta
                  fac_aniso(is,it) = (Bmag(it)/Bmag_th0)**(-eta_perp(is,ir))
               enddo
            else
               fac_aniso(is,:) = 1.0
            endif
         enddo

         allocate(dens_avg(n_species))
         allocate(dlnndr_fac(n_species,n_theta))
         allocate(dens_avg_new(n_species))
         allocate(dens_fac_new(n_species,n_theta))

         ! partial component of n/n(theta0) 
         ! -- will add the phi_rot component after the QN solve
         do is=1, n_species
            do it=1,n_theta
               dens_fac(is,it) = fac_aniso(is,it) &
                    * exp(omega_rot(ir)**2 * 0.5/vth_para(is,ir)**2 &
                    * (bigR(it)**2 - bigR_th0**2)) 
            enddo
         enddo

         phi_rot_avg = 0.0
         x = 0.05          ! initial guess for phi_rot(1)
         do it=1,n_theta
            
            n=1

            do 
               ! use Newton's method to solve the quasi-neutrality relation
               
               sum_zn  = 0.0
               dsum_zn = 0.0
               do is=1, n_species
                  fac = z(is) * dens(is,ir) * dens_fac(is,it) &
                       * exp(-z(is) / temp_para(is,ir) * x)
                  sum_zn  = sum_zn  + fac
                  dsum_zn = dsum_zn - z(is) / temp_para(is,ir) * fac
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
            
         enddo

         ! n(theta)/n(0)
         do is=1, n_species
            do it=1,n_theta
               dens_fac(is,it) = dens_fac(is,it) &
                    * exp(-z(is) / temp_para(is,ir) * phi_rot(it))
               
            enddo
         enddo

         ! <n>
         dens_avg(:)   = 0.0
         do is=1, n_species
            do it=1,n_theta
               dens_avg(is) = dens_avg(is) + &
                    w_theta(it) * dens(is,ir) * dens_fac(is,it)
            enddo
         enddo

         ! dln n(theta)/dr
         do is=1, n_species
            do it=1,n_theta
               if(aniso_model(is) == 2) then
                  fac_add   = -log(Bmag(it)/Bmag_th0) &
                       * eta_perp_rderiv(is,ir) &
                       - eta_perp(is,ir)/Bmag(it) * Bmag_rderiv(it) &
                       + eta_perp(is,ir)/Bmag_th0 * Bmag_th0_rderiv
                  
               else
                  fac_add   = 0.0
               endif
               dlnndr_fac(is,it) = -(-dlnndr(is,ir) + omega_rot_deriv(ir) &
                    * omega_rot(ir)/vth_para(is,ir)**2 &
                    * (bigR(it)**2 - bigR_th0**2) &
                    + omega_rot(ir)**2 / vth_para(is,ir)**2 &
                    * (bigR(it)*bigR_rderiv(it) &
                    - bigR_th0*bigR_th0_rderiv) &
                    - phi_rot(it) * z(is)/temp_para(is,ir) &
                    * dlntdr_para(is,ir) &
                    + 0.5 * omega_rot(ir)**2 / vth_para(is,ir)**2 &
                    * dlntdr_para(is,ir) * (bigR(it)**2 - bigR_th0**2) &
                    + fac_add)
            enddo
         enddo

         ! solve the radial derivative of QN for dphi*/dr
         do it=1,n_theta
            phi_rot_rderiv(it) = 0.0
            sum_zn  = 0.0
            do is=1, n_species
               fac = z(is) * dens(is,ir) * dens_fac(is,it)
               sum_zn = sum_zn + z(is) / temp_para(is,ir) * fac
            enddo
            if(adiabatic_ele_model == 1) then
               fac = -ne_ade(ir) * exp(1.0/te_ade(ir) * phi_rot(it))
               sum_zn = sum_zn - 1.0/te_ade(ir) * fac
            endif
            
            do is=1,n_species
               fac = z(is) * dens(is,ir) * dens_fac(is,it)
               phi_rot_rderiv(it) = phi_rot_rderiv(it) &
                    + fac * (-dlnndr_fac(is,it))
            enddo
            if(adiabatic_ele_model == 1) then
               fac = -ne_ade(ir) * exp(1.0/te_ade(ir) * phi_rot(it))
               phi_rot_rderiv(it) = phi_rot_rderiv(it) &
                    + fac * (-dlnndre_ade(ir) &
                    + phi_rot(it)/te_ade(ir) * dlntdre_ade(ir))
            endif
            phi_rot_rderiv(it) = phi_rot_rderiv(it) / sum_zn
         enddo

         ! d phi*/d theta
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
               fac = fac_aniso(is,it) &
                    * exp(omega_rot(ir)**2 * 0.5 / vth_para(is,ir)**2 &
                    * (bigR(it)**2 - bigR_th0**2) &
                    - z(is) / temp_para(is,ir) * phi_rot(it))
               rotavg_e0(is) = rotavg_e0(is) + w_theta(it) * fac
               rotavg_e1(is) = rotavg_e1(is) + w_theta(it) * fac &
                    * z(is)/temp_para(is,ir) * phi_rot(it)
               rotavg_e2(is) = rotavg_e2(is) + w_theta(it) * fac &
                    * z(is)/temp_para(is,ir) * phi_rot_rderiv(it)
               rotavg_e3(is) = rotavg_e3(is) + w_theta(it) * fac &
                    * (bigR(it)**2 - bigR_th0**2)
               rotavg_e4(is) = rotavg_e4(is) + w_theta(it) * fac &
                    * 2.0*(bigR(it)*bigR_rderiv(it) &
                    - bigR_th0*bigR_th0_rderiv)
            enddo
            rotavg_e3(is) = rotavg_e3(is) &
                 * (omega_rot(ir)/vth_para(is,ir)**2 &
                 * omega_rot_deriv(ir) - 0.5 * omega_rot(ir)**2 &
                 / vth_para(is,ir)**2 * (-dlntdr_para(is,ir)) )
            rotavg_e4(is) = rotavg_e4(is) * 0.5 * omega_rot(ir)**2 &
                 / vth_para(is,ir)**2
            rotavg_e1(is) = rotavg_e1(is) * (-dlntdr_para(is,ir))
         enddo
         
         ! redefine nm wrt Teff and keeping <nm> fixed
         do is=1,n_species
            if(aniso_model(is) == 2) then
               dens_new = dens(is,ir)
               dens_fac_new(is,:) = fac_aniso(is,:) &
                    * exp(omega_rot(ir)**2 * 0.5/vth(is,ir)**2 &
                    * (bigR(:)**2 - bigR_th0**2) &
                    - z(is) / temp(is,ir) * phi_rot(:))
               dens_avg_new(is)   = 0.0
               do it=1,n_theta
                  dens_avg_new(is) = dens_avg_new(is) &
                       + w_theta(it) * dens_fac_new(is,it)
               enddo
               dens_new   = dens_avg(is) / dens_avg_new(is)
               open(unit=1,file=trim(path)//'out.neo.diagnostic_rot',&
                    status='replace')
               do it=1,n_theta
                  write (1,'(e16.8)',advance='no') theta(it)
                  write (1,'(e16.8)',advance='no') dens(is,ir) &
                       * dens_fac(is,it)
                  write (1,'(e16.8)',advance='no') dens_new &
                       * dens_fac_new(is,it)
                  write(1,*)
               enddo
               close(1)
               nu(is,ir)  = nu(is,ir) / dens(is,ir)
               dens(is,ir) = dens_new
               dens_fac(is,:) = dens_fac_new(is,:)
               nu(is,ir)  = nu(is,ir) * dens(is,ir)
            endif
         enddo


         ! EAB test cases
         ! Teff=2/3Tperp + 1/3 Tpara
         test_flag = 0
         if(test_flag /= 0) then
            do is=1,n_species
               if(aniso_model(is) == 2) then
                  nu(is,ir)  = nu(is,ir) * temp(is,ir)**1.5 
                  temp(is,ir)   = 1.0/3.0*temp_para(is,ir) &
                       + 2.0/3.0*temp_perp(is,ir)
                  nu(is,ir)  = nu(is,ir) / temp(is,ir)**1.5
                  dlntdr(is,ir) = 1.0/3.0*temp_para(is,ir)/temp(is,ir) &
                       *dlntdr_para(is,ir) + 2.0/3.0*temp_perp(is,ir) &
                       /temp(is,ir) *dlntdr_perp(is,ir)
                  vth(is,ir) = sqrt(temp(is,ir)/mass(is))

                  ! test case b and d -- keep eta
                  if(test_flag == 1 .or. test_flag == 3) then
                     dens_new = dens(is,ir)
                     dens_fac_new(is,:) = fac_aniso(is,:) &
                          * exp(omega_rot(ir)**2 * 0.5/vth(is,ir)**2 &
                          * (bigR(:)**2 - bigR_th0**2) &
                          - z(is) / temp(is,ir) * phi_rot(:))
                  endif

                  ! test case c and e -- no eta terms in n or DKE
                  if(test_flag == 2 .or. test_flag == 4) then
                     eta_perp(is,ir)    = 0.0
                     eta_perp_rderiv(is,ir) = 0.0
                     dens_new = dens(is,ir)
                     dens_fac_new(is,:) = &
                          exp(omega_rot(ir)**2 * 0.5/vth(is,ir)**2 &
                          * (bigR(:)**2 - bigR_th0**2) &
                          - z(is) / temp(is,ir) * phi_rot(:))
                  endif
                  
                  ! test cases d and e redefine n 
                  if(test_flag == 3 .or. test_flag == 4) then
                     dens_avg_new(is)   = 0.0
                     do it=1,n_theta
                        dens_avg_new(is) = dens_avg_new(is) &
                             + w_theta(it) * dens_fac_new(is,it)
                     enddo
                     dens_new   = dens_avg(is) / dens_avg_new(is)
                  endif

                  open(unit=1,file=trim(path)//'out.neo.diagnostic_rot',&
                       status='replace')
                  do it=1,n_theta
                     write (1,'(e16.8)',advance='no') theta(it)
                     write (1,'(e16.8)',advance='no') dens(is,ir) &
                          * dens_fac(is,it)
                     write (1,'(e16.8)',advance='no') dens_new &
                          * dens_fac_new(is,it)
                     write(1,*)
                  enddo
                  close(1)

                  nu(is,ir)  = nu(is,ir) / dens(is,ir)
                  dens(is,ir) = dens_new
                  dens_fac(is,:) = dens_fac_new(is,:)
                  nu(is,ir)  = nu(is,ir) * dens(is,ir)
                  
               endif
            enddo
         endif
         
         deallocate(fac_aniso)
         deallocate(dens_avg)
         deallocate(dlnndr_fac)
         deallocate(dens_avg_new)
         deallocate(dens_fac_new)

      endif

    end subroutine ROT_solve_phi

    subroutine ROT_write(ir)
      use neo_globals
      implicit none
      integer, intent (in) :: ir
      integer :: is, it

      if(silent_flag == 0 .and. i_proc == 0 .and. rotation_model == 2) then
         open(io_rot,file=trim(path)//runfile,status='old',position='append')
         write (io_rot,'(e16.8)',advance='no') r(ir)
         write (io_rot,'(e16.8)',advance='no') phi_rot_avg
         do is=1, n_species
            write (io_rot,'(e16.8)',advance='no') 1.0/rotavg_e0(is)
            write (io_rot,'(e16.8)',advance='no') ( -rotavg_e2(is) &
                 + rotavg_e3(is) + rotavg_e4(is) + rotavg_e1(is)) &
                 / rotavg_e0(is) 
         enddo
         do it=1, n_theta
            write (io_rot,'(e16.8)',advance='no') phi_rot(it)
         enddo
         do is=1, n_species
            do it=1, n_theta
                write (io_rot,'(e16.8)',advance='no') dens_fac(is,it)
             enddo
          enddo
         write(io_rot,*)
         close(io_rot)
      endif

    end subroutine ROT_write

end module neo_rotation
