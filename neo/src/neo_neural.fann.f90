module neo_neural

  implicit none

  public :: NEURAL_alloc, NEURAL_do, NEURAL_write
  real :: jpar_nn_neo, jtor_nn_neo
  character(len=80),private :: runfile_nn = 'out.neo.transport_nn'
  integer, parameter, private :: io=1
  logical, private :: initialized = .false.

contains

  subroutine NEURAL_alloc(flag)
    use neo_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate

    if(flag == 1) then
       if(initialized) return

       if(silent_flag == 0 .and. i_proc == 0) then
          open(unit=io,file=trim(path)//runfile_nn,status='replace')
          close(io)
       endif

       initialized = .true.

    else
       if(.NOT. initialized) return

       initialized = .false.

    endif

  end subroutine NEURAL_alloc
       
  subroutine NEURAL_do(ir)
    use neo_globals
    implicit none
    integer, intent (in) :: ir
    
    call compute_nn_jpar(ir)
    if(error_status > 0) return
    
  end subroutine NEURAL_do

  ! The NN for the bootstrap current presently assumes 2 ion species
  ! (D + C) and kinetic electrons.  Strong rotation effects are not included.
  subroutine compute_nn_jpar(ir)
    use mpi
    use neo_globals
    use neo_equilibrium
    implicit none
    
    integer, intent (in) :: ir
    integer :: is, k, ierr
    integer :: is_i1, is_i2
    real    :: d_max
    real(4), dimension(6)  :: nn_in
    real(4), dimension(6) :: nn_out
    real :: CTi2_neo, CTi1_neo, CTe_neo, CNi2_neo, CNi1_neo, CNe_neo
    ! min NN training inputs
    real, dimension(5) :: xmin = (/ 0.05,1.0,-2.0,0.6,1.0 /)
    ! max NN training outputs
    real, dimension(5) :: xmax = (/ 0.35,10.0,1.0,0.99,3.0 /)
    real :: cscale_nu
    character(len=218) :: root
    character(len=255) :: data
    
    include 'brainfuse_lib.inc'
    
    ! The nn assumes 2 ion species + electrons
    if(n_species /= 3) then
       call neo_error('ERROR: (NEO) NN requires 2 ion species + electrons')
       return
    endif
    if(adiabatic_ele_model == 1) then
       call neo_error('ERROR: (NEO) NN requires 2 ion species + electrons')
       return
    endif
    
    ! Identify the main ion species
    is_i1 = -1
    d_max = -1.0
    do is=1, n_species
       if(Z(is) > 0 .and. dens(is,ir) > d_max) then
          is_i1 = is
          d_max  = dens(is,ir) 
       endif
    enddo
    if(abs(Z(is_i1) - 1.0) > epsilon(0.)) then
       if(silent_flag==0 .and. i_proc==0) then
          open(unit=io_neoout,file=trim(path)//runfile_neoout,&
               status='old',position='append')
          write(io_neoout,*) 'WARNING: (NEO) NN assumes Z(ion1)=1.0'
          close(io_neoout)
       endif
    endif
    
    ! Identify the secondary ion species
    do is=1, n_species
       if(is /= is_ele .and. is /= is_i1) then
          is_i2 = is
       endif
    enddo
    if(abs(Z(is_i2) - 6.0) > epsilon(0.)) then
       if(silent_flag==0 .and. i_proc==0) then
          open(unit=io_neoout,file=trim(path)//runfile_neoout,&
               status='old',position='append')
          write(io_neoout,*) 'WARNING: (NEO) NN assumes Z(ion2)=6.0'
          close(io_neoout)
       endif
    endif

    if(rotation_model == 2) then
       if(silent_flag==0 .and. i_proc==0) then
          open(unit=io_neoout,file=trim(path)//runfile_neoout,&
               status='old',position='append')
          write(io_neoout,*) 'WARNING: (NEO) NN assumes weak rotation'
          close(io_neoout)
       endif
    endif
    
    ! Set the input parameters for the NN
    nn_in(1) = r(ir)/rmaj(ir)                      ! r/R
    nn_in(2) = abs(q(ir))                          ! q
    ! log(nuee/cs/R)
    nn_in(3) = log10(nu(is_ele,ir)*rmaj(ir)/sqrt(temp(is_ele,ir)))
    nn_in(4) = dens(is_i1,ir)/dens(is_ele,ir)      ! ni1/ne
    nn_in(5) = temp(is_i1,ir)/temp(is_ele,ir)      ! Ti1/Te
    nn_in(6) = geo_param(ir,2)                     ! geo param: ftrap

    ! Re-scale nn_in if out-of-range of training data
    cscale_nu = 1.0
    do k=2,5
       if(nn_in(k) > xmax(k))  then
          nn_in(k) = xmax(k)
          if(k == 3) then
             cscale_nu = 10.0**xmax(k)/10.0**nn_in(k)
          endif
       endif
       if(nn_in(k) < xmin(k))  then
          nn_in(k) = xmin(k)
       endif
    enddo

    ! Run the NN
    call get_environment_variable('GACODE_ROOT',root)
    data = trim(root)//'/../neural/neonn/jbsnn/'
    ierr=load_anns(0, trim(data)//char(0),'brainfuse'//char(0))
    if(ierr == 0) then
       ! returns number of brainfuse files
       call neo_error('ERROR: (NEO) Neural network loading failed.')
       return
    endif
    ierr=load_anns_inputs(nn_in)
    ierr=run_anns()
    ierr=get_anns_avg_array(nn_out)
    if(ierr > 0) then
       call neo_error('ERROR: (NEO) Neural network failed.')
       return
    endif
    
    ! Get coeffcients computed by the NN
    CNe_neo  = nn_out(1)*cscale_nu
    CNi1_neo = nn_out(2)*cscale_nu
    CNi2_neo = nn_out(3)*cscale_nu
    CTe_neo  = nn_out(4)*cscale_nu
    CTi1_neo = nn_out(5)*cscale_nu
    CTi2_neo = nn_out(6)*cscale_nu
    
    ! Reconstruct jpar from NEO NN 
    
    jpar_nn_neo = geo_param(ir,1) * rho(ir) * temp(is_ele,ir) &
         * (abs(Z(is_ele))*dens(is_ele,ir) &
         * (CTe_neo*dlntdr(is_ele,ir) + Cne_neo*dlnndr(is_ele,ir)) &
         + abs(Z(is_i1))*dens(is_i1,ir) &
         * (CTi1_neo*dlntdr(is_i1,ir) + Cni1_neo*dlnndr(is_i1,ir)) &
         +  abs(Z(is_i2))*dens(is_i2,ir) &
            * (CTi2_neo*dlntdr(is_i2,ir) + Cni2_neo*dlnndr(is_i2,ir)))
    
    ! Toroidal component of Bootstrap current = sum <Z*n*utor/R>/<1/R>
    
    jtor_nn_neo = jpar_nn_neo*Btor2_avg/Bmag2_avg
    do is=1, n_species
       jtor_nn_neo = jtor_nn_neo + rho(ir)*geo_param(ir,1) &
            * dens(is,ir)*temp(is,ir) &
            * (dlnndr(is,ir) + dlntdr(is,ir))*(1.0-Btor2_avg/Bmag2_avg)
    enddo
    jtor_nn_neo = jtor_nn_neo/(Btor_th0*bigR_th0*bigRinv_avg)
    
  end subroutine compute_nn_jpar
  
  subroutine NEURAL_write(ir)
    use neo_globals
    implicit none
    integer, intent (in) :: ir
    
    if(silent_flag == 0 .and. i_proc == 0) then
       open(io,file=trim(path)//runfile_nn,status='old',position='append')
       write(io,'(e16.8)',advance='no') r(ir)
       write(io,'(e16.8)',advance='no') jpar_nn_neo
    endif
    
  end subroutine NEURAL_write

end module neo_neural
