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
         allocate(dens_fac(n_species,n_theta))
         allocate(dens_avg(n_species))
         allocate(dens_avg_cos(n_species))
         if(silent_flag == 0 .and. i_proc == 0 .and. rotation_model == 2) then
            open(unit=io_rot,file=trim(path)//runfile,status='replace')
            close(io_rot)
         end if
         initialized = .true.

      else
         if(.NOT. initialized) return
         deallocate(phi_rot)
         deallocate(phi_rot_deriv)
         deallocate(dens_fac)
         deallocate(dens_avg)
         deallocate(dens_avg_cos)
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
      
      if(rotation_model == 1 .or. spitzer_model==1) then
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
     
      
    end subroutine ROT_solve_phi

    subroutine ROT_write(ir)
      use neo_globals
      implicit none
      integer, intent (in) :: ir
      integer :: is, it

      if(silent_flag == 0 .and. i_proc == 0 .and. rotation_model == 2) then
         open(io_rot,file=trim(path)//runfile,status='old',position='append')
         write (io_rot,'(e16.8,$)') r(ir)
         write (io_rot,'(e16.8,$)') phi_rot_avg
         do is=1, n_species
            write (io_rot,'(e16.8,$)') dens_avg(is)
         enddo
         do it=1, n_theta
            write (io_rot,'(e16.8,$)') phi_rot(it)
         enddo
         write(io_rot,*)
         close(io_rot)
      endif

    end subroutine ROT_write


end module neo_rotation
