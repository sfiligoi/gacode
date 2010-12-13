module neo_rotation

  implicit none

  public :: ROT_alloc, ROT_solve_phi, ROT_write
  real, dimension(:), allocatable   :: phi_rot        ! phi(r)-phi(r,theta0)
  real, dimension(:), allocatable   :: phi_rot_deriv  ! theta derivative
  real                              :: phi_rot_avg    ! flux surface avg
  real, dimension(:,:), allocatable :: dens_fac       ! (ns,nth):
                                                      ! n=n(th0)*dens_fac
  real, dimension(:), allocatable   :: dens_avg       ! (ns): <n>/n0
  real, dimension(:), allocatable   :: dens_avg_cos

  logical, private :: initialized = .false.
  integer, private :: io_rot=50

  contains

    subroutine ROT_alloc(flag)
      use neo_globals, only : n_species, n_theta, write_out_mode, rotation_model
      implicit none
      integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
      
      if(flag == 1) then
         if(initialized) return
         allocate(phi_rot(n_theta))
         allocate(phi_rot_deriv(n_theta))
         allocate(dens_fac(n_species,n_theta))
         allocate(dens_avg(n_species))
         allocate(dens_avg_cos(n_species))
         if(write_out_mode > 0) then
            if(rotation_model == 2) then
               open(unit=io_rot,file='rotation.out',status='replace')
            endif
         end if
         initialized = .true.
         
      else
         if(.NOT. initialized) return
         deallocate(phi_rot)
         deallocate(phi_rot_deriv)
         deallocate(dens_fac)
         deallocate(dens_avg)
         deallocate(dens_avg_cos)
         if(write_out_mode > 0) then
            if(rotation_model == 2) then
               close(io_rot)
            endif
         endif
         initialized = .false.

      endif

    end subroutine ROT_alloc

    ! solve the quasi-neutrality relation for poloidal part of phi
    subroutine ROT_solve_phi(ir)
      use neo_globals
      use neo_equilibrium
      implicit none
      integer, intent (in) :: ir
      integer, parameter :: nmax = 100
      integer :: it, is, jt, id, n
      real :: x, x0, sum_zn, dsum_zn, fac
      
      if(rotation_model == 1 .or. case_spitzer .or. zf_model == 1) then
         do it=1, n_theta
            phi_rot(it) = 0.0
            phi_rot_deriv(it) = 0.0
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
               print *, 'Rotation density computation failed to converge'
               stop
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

         !do is=1, n_species
         !   print *, dens_avg(is), dens_avg_cos(is)
         !enddo
         !stop

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

         x = 0.0
         do it=1, n_theta
            x = x + bigR(it)**2 * w_theta(it)
         enddo

      endif
      
      ! call ROT_theory(ir)
      
    end subroutine ROT_solve_phi

    subroutine ROT_write(ir)
      use neo_globals
      implicit none
      integer, intent (in) :: ir
      integer :: is, it

      if(write_out_mode == 0 .or. rotation_model == 1) return

      write (io_rot,'(e16.8,$)') r(ir)
      write (io_rot,'(e16.8,$)') phi_rot_avg
      do is=1, n_species
         write (io_rot,'(e16.8,$)') dens_avg(is)
      enddo
      do it=1, n_theta
         write (io_rot,'(e16.8,$)') phi_rot(it)
      enddo
      write(io_rot,*)

    end subroutine ROT_write

    subroutine ROT_theory(ir)
      use neo_globals
      use neo_equilibrium
      implicit none
      integer, intent (in) :: ir
      integer :: is, it, Zefftype, is_ion1, is_ion2
      real :: ft, nis1, nis2, nis2_inv, L31, L32, L33, Pi_flux, bigR2_avg
      
      Zefftype = 0
      do is=1, n_species
         if(Z(is) > 0) then
            if(Zefftype == 0) then
               is_ion1 = is
               Zefftype = 1
            else if(Zefftype == 1) then
               is_ion2 = is
               Zefftype = 2
               exit
            endif
         endif
      enddo

      if(Zefftype == 0 .or. Zefftype == 1) return
      

      call compute_fractrap_Zeff(ir, Zefftype, is_ion1, is_ion2, ft)

      if(Zefftype == 1) then
         nis2 = 1.0
         nis2_inv = 1.0
      else
         nis1 = 0.0
         nis2 = 0.0
         nis2_inv = 0.0
         do it=1, n_theta
            nis1 = nis1 + w_theta(it) &
                 * (dens(is_ion1,ir) * dens_fac(is_ion1,it))
            nis2 = nis2 + w_theta(it) &
                 * (dens(is_ion2,ir) * dens_fac(is_ion2,it))
            nis2_inv = nis2_inv + w_theta(it) &
                 / (dens(is_ion2,ir) * dens_fac(is_ion2,it))
         enddo
         nis2_inv = nis2_inv * nis2
      end if

      bigR2_avg = 0.0
      do it=1, n_theta
         bigR2_avg = bigR2_avg + w_theta(it) * bigR(it)**2
      enddo

      L31 = ft * (1.0 - ft/(nis2_inv - (1.0 - ft)))
      L32 = -1.5 * L31
      L33 = L31 * omega_rot(ir) * omega_rot_deriv(ir) * bigR2_avg &
           * (mass(is_ion1)/temp(is_ion1,ir) &
           - mass(is_ion2)/(temp(is_ion2,ir) * Z(is_ion2)))

      ! EAB NOTE: L31 term dN/dpsi does not contain the d<phi_star>/dpsi term
      ! -- ok as long as er0 is zero
      Pi_flux = -nis1 * temp(is_ion1,ir) * I_div_psip**2 &
           * rho(ir)**2 / (Z(is_ion1)**2 * Bmag2_avg) &
           * nu(is_ion1,ir) * Z(is_ion2)**2 / Z(is_ion1)**2 &
           * nis2/ dens(is_ion1,ir) &
           * (4.0 / (3.0 * sqrt(pi))) &
           * mass(is_ion1) * omega_rot(ir) * bigR2_avg &
           * (L31 * (-dlnndr(is_ion1,ir) &
           - dlntdr(is_ion1,ir) * (1.0 &
           + (omega_rot(ir) / vth(is,ir))**2 * 0.5 * bigR_th0**2 &
           + Z(is) / temp(is,ir) * phi_rot(it)) &
           - omega_rot(ir) * bigR_th0**2 &
           / vth(is,ir)**2 * omega_rot_deriv(ir)) &
           - L32 * dlntdr(is_ion1,ir) + L33)

      if(write_out_mode > 1) then
         print *, 'ROT_theory Pi: ', Pi_flux
         print *, 'ROT_theory coeffs: ', L31, L32, L33
      endif

    end subroutine ROT_theory

    ! Computes the fraction of trapped particles
    subroutine compute_fractrap_Zeff(ir, Zefftype, is_ion1, is_ion2, ft)
      use neo_globals
      use neo_equilibrium, only : w_theta, Bmag, Bmag2_avg
      implicit none
      integer, intent (in) :: ir, Zefftype, is_ion1, is_ion2
      real, intent (inout) :: ft
      integer :: nlambda=500
      real, dimension(:), allocatable :: lambda
      real :: dlambda, fac_lambda, sum_th
      integer i,it
      real :: Bmax, Zeff_avg, Zeff, eps
      
      eps = r(ir)/rmaj(ir)

      Bmax = Bmag(1)
      do it=2, n_theta
         if(Bmag(it) > Bmax) then
            Bmax = Bmag(it)
         endif
      enddo
      
      allocate(lambda(nlambda))
      dlambda = 1.0/(nlambda-1)
      do i=1, nlambda
         lambda(i) = (i-1)*dlambda
      enddo
      
      ft = 0.0
      
      do i=2, nlambda-1
         
         ! open integration for lambda
         if(i==2 .or. i == nlambda-1)  then
            fac_lambda = 23.0/12.0
         else if(i==3 .or. i == nlambda-2) then
            fac_lambda = 7.0/12.0
         else
            fac_lambda = 1.0
         endif
         
         ! closed integration for th: <sqrt(1-lambda*B/Bmax)>
         sum_th = 0.0
         do it=1, n_theta
            if(Zefftype == 1) then
               Zeff = 1.0
            else
               Zeff = 1.0 + Z(is_ion2)*Z(is_ion2)&
                    * dens(is_ion2,ir)/dens(is_ion1,ir) &
                    * dens_fac(is_ion2,it) / dens_fac(is_ion1,it)
            end if
            sum_th = sum_th + w_theta(it) &
                 * Zeff * sqrt(1.0 - lambda(i) * Bmag(it) / Bmax)
         enddo
         
         ft = ft + fac_lambda * lambda(i) / sum_th;
         
      end do
      
      if(Zefftype == 1) then
         Zeff_avg = 1.0
      else
         Zeff_avg = 1.0
         do it=1, n_theta
             Zeff_avg =  Zeff_avg + w_theta(it)* Z(is_ion2)*Z(is_ion2)&
                  * dens(is_ion2,ir)/dens(is_ion1,ir) &
                  * dens_fac(is_ion2,it) / dens_fac(is_ion1,it)
         enddo
      endif

      ft = Zeff_avg * ft * dlambda * 0.75 * Bmag2_avg / Bmax**2
      ft = 1.0 - ft
      
      deallocate(lambda)
      
    end subroutine compute_fractrap_Zeff

end module neo_rotation
