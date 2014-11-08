module cgyro_field

  implicit none

  public :: FIELD_alloc, FIELDh_do, FIELDx_do
  logical, private :: initialized = .false.
  real, private :: sum_den_h
  real, dimension(:,:), allocatable, private :: sum_den_x, sum_cur_x
  real, dimension(:,:,:), allocatable, private :: hzf, xzf ! for n=0 test
  real, dimension(:), allocatable, private :: pvec_in, pvec_outr, pvec_outi

contains
  
  subroutine FIELD_alloc(flag)
    use cgyro_globals
    use cgyro_gyro
    use cgyro_equilibrium
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: is, ie, ix, ir, it, jt
    integer :: info
    integer, dimension(:), allocatable :: i_piv
    real, dimension(:), allocatable :: work

    if (flag == 1) then
       if (initialized) return
    
       ! Pre-factors for Poisson eqn

       sum_den_h = 0.0
       do is=1,n_species
          do ie=1,n_energy
             do ix=1,n_xi
                sum_den_h = sum_den_h &
                     + 0.5 * w_xi(ix) &
                     * z(is)**2/temp(is) *dens(is) * w_e(ie)
             enddo
          enddo
       enddo
       if(adiabatic_ele_model == 1) then
          sum_den_h = sum_den_h + dens_ele / temp_ele
       endif

       allocate(sum_den_x(n_radial,n_theta))
       do ir=1,n_radial
          do it=1,n_theta
             sum_den_x(ir,it) = 0.0
             do is=1,n_species
                do ie=1,n_energy
                   do ix=1,n_xi
                      sum_den_x(ir,it) = sum_den_x(ir,it) &
                           + 0.5 * w_xi(ix) &
                           * (1.0 - gyrox_J0(is,ir,it,ie,ix)**2) &
                           * z(is)**2/temp(is) *dens(is) * w_e(ie)
                   enddo
                enddo
             enddo
             if(adiabatic_ele_model == 1) then
                sum_den_x(ir,it) = sum_den_x(ir,it) + dens_ele / temp_ele
             endif
          enddo
       enddo

       ! Special factors in zf n=0 Poisson eqn if adiabatic electrons

       if(toroidal_model == 2 .and. adiabatic_ele_model == 1) then

          allocate(hzf(n_radial,n_theta,n_theta))
          hzf(:,:,:) = 0.0      
          do ir=1,n_radial
             do it=1,n_theta
                hzf(ir,it,it) = -k_perp(it,ir)**2 * lambda_debye**2 &
                     * dens_ele / temp_ele + sum_den_h
                do jt=1,n_theta
                   hzf(ir,it,jt) = hzf(ir,it,jt) &
                        - dens_ele / temp_ele * w_theta(jt)
                enddo
             enddo
          enddo
          allocate(work(n_theta))
          allocate(i_piv(n_theta))
          do ir=1,n_radial
             call DGETRF(n_theta,n_theta,hzf(ir,:,:),n_theta,i_piv,info)
             call DGETRI(n_theta,hzf(ir,:,:),n_theta,i_piv,work,n_theta,info)
          enddo
          deallocate(i_piv)
          deallocate(work)

          allocate(xzf(n_radial,n_theta,n_theta))
          xzf(:,:,:) = 0.0     
          do ir=1,n_radial
             do it=1,n_theta
                xzf(ir,it,it) = -k_perp(it,ir)**2 * lambda_debye**2 &
                     * dens_ele / temp_ele + sum_den_x(ir,it)
                do jt=1,n_theta
                   xzf(ir,it,jt) = xzf(ir,it,jt) &
                        - dens_ele / temp_ele * w_theta(jt)
                enddo
             enddo
          enddo
          allocate(work(n_theta))
          allocate(i_piv(n_theta))
          do ir=1,n_radial
             call DGETRF(n_theta,n_theta,xzf(ir,:,:),n_theta,i_piv,info)
             call DGETRI(n_theta,xzf(ir,:,:),n_theta,i_piv,work,n_theta,info)
          enddo
          deallocate(i_piv)
          deallocate(work)

          allocate(pvec_in(n_theta))
          allocate(pvec_outr(n_theta))
          allocate(pvec_outi(n_theta))

       endif

       ! Pre-factors for Ampere eqn

       if(n_field > 1) then
          allocate(sum_cur_x(n_radial,n_theta))
          do ir=1,n_radial
             do it=1,n_theta
                sum_cur_x(ir,it) = 0.0
                do is=1,n_species
                   do ie=1,n_energy
                      do ix=1,n_xi
                         sum_cur_x(ir,it) = sum_cur_x(ir,it) &
                              + 0.5 * w_xi(ix) * w_e(ie) &
                              * gyrox_J0(is,ir,it,ie,ix)**2 &
                              * xi(ix)**2 * 2.0 * energy(ie) &
                              * z(is)**2/temp(is) *dens(is) 
                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif

       initialized = .true.
       
    else
       if(.NOT. initialized) return

       deallocate(sum_den_x)

       if(toroidal_model == 2 .and. adiabatic_ele_model == 1) then
          deallocate(hzf)
          deallocate(xzf)
          deallocate(pvec_in)
          deallocate(pvec_outr)
          deallocate(pvec_outi)
       endif

       if(n_field > 1) then
          deallocate(sum_cur_x)
       endif

       initialized = .false.
    endif

  end subroutine FIELD_alloc

  subroutine FIELDh_do
     use mpi
     use timer_lib
     use cgyro_globals
     use cgyro_gyro
     use cgyro_equilibrium
     implicit none
     integer :: is, ie, ix, ir, it
     complex :: fac

     field_loc(:,:,:) = (0.0,0.0)

     ! Poisson and Ampere RHS integrals of H

     ic_loc = 0
     do ic=nc1,nc2
        ic_loc = ic_loc+1
        
        it = it_c(ic)
        ir = ir_c(ic)
        
        do iv=1,nv
           
           is = is_v(iv)
           ix = ix_v(iv)
           ie = ie_v(iv)
           
           fac = gyrox_J0(is,ir,it,ie,ix) &
                * z(is)*dens(is) * w_e(ie) * 0.5 * w_xi(ix)*cap_h_v(ic_loc,iv)

           field_loc(ir,it,1) = field_loc(ir,it,1) + fac 

           if(n_field > 1) then
              field_loc(ir,it,2) = field_loc(ir,it,2) + fac &
                   * xi(ix) * sqrt(2.0 * energy(ie))
           endif
           
        enddo
     enddo

     call MPI_ALLREDUCE(field_loc(:,:,:),&
          field(:,:,:),&
          size(field(:,:,:)),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     ! Poisson LHS factors

     if (toroidal_model == 2 .and. adiabatic_ele_model == 1) then
        
        do ir=1,n_radial
           pvec_in(:) = real(field(ir,:,1))
           call DGEMV('N',n_theta,n_theta,num1,hzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outr(:),1)
           pvec_in(:) = imag(field(ir,:,1))
           call DGEMV('N',n_theta,n_theta,num1,hzf(ir,:,:),&
                n_theta,pvec_in(:),1,num0,pvec_outi(:),1)
           field(ir,:,1) = pvec_outr(:) + i_c * pvec_outi(:)
        enddo
        
     else
        
        do ir=1,n_radial
           do it=1,n_theta
              field(ir,it,1) = field(ir,it,1) &
                   / (-k_perp(it,ir)**2 * lambda_debye**2 &
                   * dens_ele / temp_ele + sum_den_h)
           enddo
        enddo
        
     endif

     ! Ampere LHS factors

     print *,dens_ele, temp_ele

     if(n_field > 1) then
        do ir=1,n_radial
           do it=1,n_theta
              field(ir,it,2) = field(ir,it,2) &
                   / (-2.0*k_perp(it,ir)**2 * rho**2 / betae_unit & 
                   * dens_ele * temp_ele)
           enddo
        enddo
     endif

  end subroutine FIELDh_do

  subroutine FIELDx_do
    use mpi
    use timer_lib
    use cgyro_globals
    use cgyro_gyro
    use cgyro_equilibrium
    implicit none
    integer :: is, ie, ix, ir, it
    complex :: fac

    call timer_lib_in('fieldx')

    field_loc(:,:,:) = (0.0,0.0)

    ! Poisson and Ampere RHS integrals of h
    
    iv_loc = 0
    do iv=nv1,nv2
       
       iv_loc = iv_loc+1
       
       is = is_v(iv)
       ix = ix_v(iv)
       ie = ie_v(iv)

       do ic=1,nc

          ir = ir_c(ic)
          it = it_c(ic)

          
          fac = gyrox_J0(is,ir,it,ie,ix) &
               * z(is)*dens(is) * w_e(ie) * 0.5 * w_xi(ix)*h_x(ic,iv_loc)
          
          field_loc(ir,it,1) = field_loc(ir,it,1) + fac 

          if(n_field > 1) then
             field_loc(ir,it,2) = field_loc(ir,it,2) + fac &
                  * xi(ix) * sqrt(2.0 * energy(ie))
          endif
          
       enddo
    enddo

     call MPI_ALLREDUCE(field_loc(:,:,:),&
          field(:,:,:),&
          size(field(:,:,:)),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

    ! Poisson LHS factors

    if (toroidal_model == 2 .and. adiabatic_ele_model == 1) then
       do ir=1,n_radial
          pvec_in(:) = real(field(ir,:,1))
          call DGEMV('N',n_theta,n_theta,num1,xzf(ir,:,:),&
               n_theta,pvec_in(:),1,num0,pvec_outr(:),1)
          pvec_in(:) = imag(field(ir,:,1))
          call DGEMV('N',n_theta,n_theta,num1,xzf(ir,:,:),&
               n_theta,pvec_in(:),1,num0,pvec_outi(:),1)
          field(ir,:,1) = pvec_outr(:) + i_c * pvec_outi(:)
       enddo
    else
       do ir=1,n_radial
          do it=1,n_theta
             field(ir,it,1) = field(ir,it,1) &
                  / (-k_perp(it,ir)**2 * lambda_debye**2 &
                  * dens_ele / temp_ele + sum_den_x(ir,it))
          enddo
       enddo
    endif

    ! Ampere LHS factors

    if(n_field > 1) then
       do ir=1,n_radial
          do it=1,n_theta
             field(ir,it,2) = field(ir,it,2) &
                  / (-2.0*k_perp(it,ir)**2 * rho**2 / betae_unit & 
                  * dens_ele * temp_ele + sum_cur_x(ir,it))
          enddo
       enddo
    endif

    ! Compute H given h and phi(h) apar(h)

    iv_loc = 0
    do iv=nv1,nv2
 
       iv_loc = iv_loc+1

       is = is_v(iv)
       ix = ix_v(iv)
       ie = ie_v(iv)
       
       do ic=1,nc
          
          ir = ir_c(ic)
          it = it_c(ic)
          
          cap_h_c(ic,iv_loc) = h_x(ic,iv_loc) &
               + z(is)/temp(is) * gyrox_J0(is,ir,it,ie,ix) &
               * field(ir,it,1)
          
          if(n_field > 1) then
             cap_h_c(ic,iv_loc) = cap_h_c(ic,iv_loc) &
             - z(is)/temp(is) * gyrox_J0(is,ir,it,ie,ix) &
                  * xi(ix)*sqrt(2.0*energy(ie))*field(ir,it,2)
          endif

       enddo
    enddo


  call timer_lib_out('fieldx')

  end subroutine FIELDx_do

end module cgyro_field
