module neo_3d_transport

  implicit none

  public :: ThreeD_TRANSP_alloc, ThreeD_TRANSP_do, ThreeD_TRANSP_write

  real, dimension(:,:), allocatable :: pflux        ! (ns): gamma/(n0*vt0)
  real, dimension(:,:), allocatable :: eflux        ! (ns): Q/(n0*vt0*T0)
  real, dimension(:,:), allocatable :: uparB        ! (ns): upar*B/(vto*B0)
  real, dimension(:),   allocatable :: jpar         ! sum(Z*upar*B*n) 

  integer, parameter, private :: io2=1, io3=2
  character(len=80),private :: runfile_transp2 = 'out.neo.transport_2d'
  character(len=80),private :: runfile_transp3 = 'out.neo.transport_3d'
  logical, private :: initialized = .false.

contains

  subroutine ThreeD_TRANSP_alloc(flag)
    use neo_globals
    use neo_3d_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate

    if(flag == 1) then
       if(initialized) return

       allocate(pflux(n_species,3))
       allocate(eflux(n_species,3))
       allocate(uparB(n_species,3))
       allocate(jpar(3))

       if(silent_flag == 0 .and. i_proc == 0) then
          open(unit=io2,file=trim(path)//runfile_transp2,status='replace')
          close(io2)
          open(unit=io3,file=trim(path)//runfile_transp3,status='replace')
          close(io3)
       endif

       initialized = .true.

    else
       if(.NOT. initialized) return
       
       deallocate(pflux)
       deallocate(eflux)
       deallocate(uparB)
       deallocate(jpar)

       initialized = .false.

    endif
    
  end subroutine ThreeD_TRANSP_alloc

  subroutine ThreeD_TRANSP_do(ir,nonaxisym_flag)
    use neo_globals
    use neo_energy_grid
    use neo_3d_globals
    use neo_3d_equilibrium
    implicit none
    integer :: i, is, ie, ix, it, ip, j1, j2, j
    integer, intent (in) :: ir
    integer, intent (in) :: nonaxisym_flag

    ! Compute the neoclassical transport coefficients

    if(nonaxisym_flag == 1) then
       j1 = 3
       j2 = 3
    else
       j1 = 1
       j2 = 2
    endif

    do j=j1,j2
       pflux(:,j)  = 0.0
       eflux(:,j)  = 0.0
       uparB(:,j)  = 0.0
       jpar(j)     = 0.0
    enddo

    do i=1,n_row

       is = is_indx(i)
       ie = ie_indx(i)
       ix = ix_indx(i)
       it = it_indx(i)
       ip = ip_indx(i)


       if (ix == 0) then

          pflux(is,j1) = pflux(is,j1) + g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) / sum_w_theta_0 &
               * (4.0/3.0) * evec_e1(ie,ix) &
               * driftx_2d(is,it) * w_theta_0(it) 
          
          eflux(is,j1) = eflux(is,j1) + g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) / sum_w_theta_0 &
               * temp(is,ir) * (4.0/3.0) * evec_e2(ie,ix) &
               * driftx_2d(is,it) * w_theta_0(it)

          if(nonaxisym_flag /= 0) then
             pflux(is,j2) = pflux(is,j2) + g(i) &
                  * 4.0/sqrt(pi) * dens(is,ir) / sum_w_theta_0 &
                  * (4.0/3.0) * evec_e1(ie,ix) &
                  * (-driftx_2d(is,it) * w_theta_0(it) &
                  * sum_w_theta_1/sum_w_theta_0 &
                  + driftx_2d(is,it) * w_theta_1(it,ip) &
                  + driftx_3d(is,it,ip) * w_theta_0(it))

             eflux(is,j2) = eflux(is,j2) + g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) / sum_w_theta_0 &
               * temp(is,ir) * (4.0/3.0) * evec_e2(ie,ix) &
               * (-driftx_2d(is,it) * w_theta_0(it) &
               * sum_w_theta_1/sum_w_theta_0 &
               + driftx_2d(is,it) * w_theta_1(it,ip) &
               + driftx_3d(is,it,ip) * w_theta_0(it))
          endif
          
       else if(ix == 1) then
          ! uparB = < B * 1/n * int vpar * (F0 g)>
          uparB(is,j1) = uparB(is,j1) + g(i) &
                * 4.0/sqrt(pi) / sum_w_theta_0 &
                * (1.0/3.0) * evec_e05(ie,ix) * sqrt(2.0) * vth(is,ir) &
                * Bmag_0(it) * w_theta_0(it)

          if(nonaxisym_flag /= 0) then
             uparB(is,j2) = uparB(is,j2) + g(i) &
                  * 4.0/sqrt(pi) / sum_w_theta_0 &
                  * (1.0/3.0) * evec_e05(ie,ix) * sqrt(2.0) * vth(is,ir) &
                  * (-Bmag_0(it) * w_theta_0(it) &
                  * sum_w_theta_1/sum_w_theta_0 &
                  + Bmag_0(it) * w_theta_1(it,ip) &
                  + Bmag_1(it,ip) * w_theta_0(it))
          endif

       else if(ix == 2) then
          pflux(is,j1) = pflux(is,j1) + g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) / sum_w_theta_0 &
               * (2.0/15.0) * evec_e1(ie,ix) &
               * driftx_2d(is,it) * w_theta_0(it) 
          
          eflux(is,j1) = eflux(is,j1) + g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) / sum_w_theta_0 &
               * temp(is,ir) * (2.0/15.0) * evec_e2(ie,ix) &
               * driftx_2d(is,it) * w_theta_0(it)

          if(nonaxisym_flag /= 0) then
             pflux(is,j2) = pflux(is,j2) + g(i) &
                  * 4.0/sqrt(pi) * dens(is,ir) / sum_w_theta_0 &
                  * (2.0/15.0) * evec_e1(ie,ix) &
                  * (-driftx_2d(is,it) * w_theta_0(it) &
                  * sum_w_theta_1/sum_w_theta_0 &
                  + driftx_2d(is,it) * w_theta_1(it,ip) &
                  + driftx_3d(is,it,ip) * w_theta_0(it))

             eflux(is,j2) = eflux(is,j2) + g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) / sum_w_theta_0 &
               * temp(is,ir) * (2.0/15.0) * evec_e2(ie,ix) &
               * (-driftx_2d(is,it) * w_theta_0(it) &
               * sum_w_theta_1/sum_w_theta_0 &
               + driftx_2d(is,it) * w_theta_1(it,ip) &
               + driftx_3d(is,it,ip) * w_theta_0(it))
          endif
          
       end if
       
    end do

    ! Bootstrap current = sum <Z*n*upar B>
    do j=j1,j2
       do is=1, n_species
          jpar(j) = jpar(j) + Z(is) * dens(is,ir) * uparB(is,j)
       enddo
    enddo
    
  end subroutine ThreeD_TRANSP_do

  subroutine ThreeD_TRANSP_write(ir)
    use neo_globals
    implicit none
    integer, intent (in) :: ir
    integer :: is

    if(silent_flag > 0 .or. i_proc > 0) return

    ! 2D transport coefficients (normalized)
    open(io2,file=trim(path)//runfile_transp2,status='old',position='append')
    write (io2,'(e16.8)',advance='no') r(ir)
    write (io2,'(e16.8)',advance='no') jpar(1)
    do is=1, n_species
       write (io2,'(e16.8)',advance='no') pflux(is,1)
       write (io2,'(e16.8)',advance='no') eflux(is,1)
       write (io2,'(e16.8)',advance='no') uparB(is,1)
    enddo
    write (io2,*)
    close(io2)

    ! 3D transport coefficients (normalized)
    open(io3,file=trim(path)//runfile_transp3,status='old',position='append')
    write (io3,'(e16.8)',advance='no') r(ir)
    write (io3,'(e16.8)',advance='no') jpar(2)+jpar(3)
    do is=1, n_species
       write (io3,'(e16.8)',advance='no') pflux(is,2)+pflux(is,3)
       write (io3,'(e16.8)',advance='no') eflux(is,2)+eflux(is,3)
       write (io3,'(e16.8)',advance='no') uparB(is,2)+uparB(is,3)
    enddo
    write (io3,*)
    close(io3)

  end subroutine ThreeD_TRANSP_write

end module neo_3d_transport
