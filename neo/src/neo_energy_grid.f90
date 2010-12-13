module neo_energy_grid
  
  implicit none

  public :: ENERGY_basis_ints_alloc, ENERGY_basis_ints, &
       ENERGY_coll_ints_alloc, ENERGY_coll_ints, get_coll_freqs

  ! private variables for integration
  integer, private :: iebasis, jebasis, ietype  ! global params for integration
  integer, private :: ir, is, js                ! for collision ints
  integer, parameter :: integ_order = 64
  integer, parameter, private :: Nx_start=5, Nx_max=200000
  logical, private :: initialized_basis = .false.
  logical, private :: initialized_coll  = .false.
  real, private :: itol                         ! integ convergence tol (input)
  real, private :: energy_min                   ! min ene for integration

  ! energy integral values
  !
  ! matrix index
  !
  integer, dimension(:,:), allocatable :: eindx !indx for matrices (symmetric)
  !
  ! Kinetic vars
  !
  ! vectors (ie)
  real, dimension(:), allocatable :: evec_e05,evec_e105 ! ene^0.5; ene^1.5
  real, dimension(:), allocatable :: evec_e0, evec_e1, evec_e2, evec_e3
  real, dimension(:), allocatable :: evec_e1de, evec_e2de, evec_e3de
  ! matrices (ie/je)
  real, dimension(:), allocatable :: emat_e05,emat_e1, emat_e2, emat_en05, &
       emat_e0
  ! non-symmetric matrices
  real, dimension(:,:), allocatable :: ematij_e2de, ematij_e05de, ematij_e1de
                                       ! ene^2*dBn/dene, ene^0.5*dBn/dene 
  !
  ! Collision vars
  !
  ! matrices (is,js,ie/je)
  real, dimension(:,:,:), allocatable :: emat_coll_lorentz, emat_coll_ru, &
       emat_coll_rpi, emat_coll_rn
  ! vector (is,js,ie)
  real, dimension(:,:,:), allocatable :: evec_coll_rs,evec_coll_rh, &
       evec_coll_rk, evec_coll_rp, evec_coll_rq,evec_coll_rqd, &
       evec_coll_rn_krook
  ! constant (is,js)
  real, dimension(:,:), allocatable :: econ_coll_rs, econ_coll_rh, &
       econ_coll_rk, econ_coll_rp, econ_coll_rq, econ_coll_rn_krook
  
contains
  
  subroutine ENERGY_basis_ints_alloc(flag)
    use neo_globals, only : n_energy, energy_max, &
         energy_tol, energy_min_connor, collision_model
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: ie, je, ke
    
    if(flag == 1) then
       if(initialized_basis) return
       
       allocate(eindx(n_energy,n_energy))
       ke=1
       do ie=1, n_energy
          do je=ie, n_energy
             eindx(ie,je) = ke
             if(je .ne. ie) then
                eindx(je,ie) = ke
             endif
             ke = ke + 1
          enddo
       enddo
       
       allocate(evec_e05(n_energy));   allocate(evec_e105(n_energy))
       allocate(evec_e0(n_energy));    allocate(evec_e2(n_energy))
       allocate(evec_e1(n_energy));    allocate(evec_e3(n_energy))
       allocate(evec_e2de(n_energy));  allocate(evec_e3de(n_energy))
       allocate(evec_e1de(n_energy))

       ! matrices are symmetric: (n^2+n)/2 elements
       ke = (n_energy**2 + n_energy)/2
       
       allocate(emat_en05(ke))
       allocate(emat_e05(ke))
       allocate(emat_e0(ke))
       allocate(emat_e1(ke))
       allocate(emat_e2(ke))
       ! non-symmetric matrix
       allocate(ematij_e2de(n_energy, n_energy))
       allocate(ematij_e05de(n_energy, n_energy))
       allocate(ematij_e1de(n_energy, n_energy))

       itol = energy_tol
       if(collision_model == 1) then
          energy_min = energy_min_connor
       else
          energy_min = 0.0
       endif

       initialized_basis = .true.
       
    else
       if(.NOT. initialized_basis) return
       
       deallocate(eindx)

       deallocate(evec_e05);   deallocate(evec_e105)
       deallocate(evec_e0);    deallocate(evec_e2)
       deallocate(evec_e1);    deallocate(evec_e3)
       deallocate(evec_e2de);  deallocate(evec_e3de)
       deallocate(evec_e1de)

       deallocate(emat_e05);   deallocate(emat_e1)
       deallocate(emat_e2);    deallocate(emat_e0)
       deallocate(ematij_e2de)
       deallocate(emat_en05);  deallocate(ematij_e05de)
       deallocate(ematij_e1de)

       initialized_basis = .false.

    endif
    
  end subroutine ENERGY_basis_ints_alloc

  subroutine ENERGY_coll_ints_alloc(flag)
    use neo_globals, only : n_energy, n_species
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: ke
    
    if(flag == 1) then
       if(initialized_coll) return
       
       ! matrices are symmetric: (n^2+n)/2 elements
       ke = (n_energy**2 + n_energy)/2
       
       allocate(emat_coll_lorentz(n_species,n_species,ke))
       allocate(emat_coll_ru(n_species,n_species,ke))
       allocate(emat_coll_rpi(n_species,n_species,ke))
       allocate(emat_coll_rn(n_species,n_species,ke))
       allocate(evec_coll_rs(n_species,n_species,n_energy))
       allocate(econ_coll_rs(n_species,n_species))
       allocate(evec_coll_rh(n_species,n_species,n_energy))
       allocate(econ_coll_rh(n_species,n_species))
       allocate(evec_coll_rk(n_species,n_species,n_energy))
       allocate(econ_coll_rk(n_species,n_species))
       allocate(evec_coll_rp(n_species,n_species,n_energy))
       allocate(econ_coll_rp(n_species,n_species))
       allocate(evec_coll_rq(n_species,n_species,n_energy))
       allocate(econ_coll_rq(n_species,n_species))
       allocate(evec_coll_rqd(n_species,n_species,n_energy))
       allocate(econ_coll_rn_krook(n_species,n_species))
       allocate(evec_coll_rn_krook(n_species,n_species,n_energy))

       initialized_coll = .true.
       
    else
       if(.NOT. initialized_coll) return

       deallocate(emat_coll_lorentz); deallocate(emat_coll_ru)
       deallocate(emat_coll_rpi);     deallocate(emat_coll_rn)
       deallocate(evec_coll_rs);      deallocate(econ_coll_rs)
       deallocate(evec_coll_rh);      deallocate(econ_coll_rh)
       deallocate(evec_coll_rk);      deallocate(econ_coll_rk)
       deallocate(evec_coll_rp);      deallocate(econ_coll_rp)
       deallocate(evec_coll_rq);      deallocate(econ_coll_rq)
       deallocate(evec_coll_rqd)
       deallocate(evec_coll_rn_krook)
       deallocate(econ_coll_rn_krook)

       initialized_coll = .false.

    endif
    
  end subroutine ENERGY_coll_ints_alloc

  ! Compute the basis function value for B_ie(x_arg) -> x_val
  subroutine  get_basis_val(ie,x_arg,x_val,dtype)
    use neo_globals, only : pi
    implicit none
    integer, intent(in) :: ie
    real, intent(in) ::  x_arg
    real, intent(out) :: x_val
    integer, intent(in) :: dtype ! type of deriv (0->B; 1->B_prime)
    real :: b_arg
    
    ! Chebyshev polynomials for z
    ! T_n(x) = cos(n*arcos(x))
    ! T_n(x)_prime = n / sqrt(1-x^2) * sin(n*arcos(x))
    
    if(x_arg <= -1.0) then
       b_arg = pi  ! use x = -1
    else if(x_arg >= 1.0) then
       b_arg = 0.0 ! use x = 1
    else
       b_arg = acos(x_arg)
    endif

    if(dtype == 0) then
       x_val = cos((ie-1) * b_arg)
    else if(dtype == 1) then
       if(x_arg >= 1.0) then
          x_val = 1.0 * (ie-1) * (ie-1)
       else if(x_arg <= -1.0) then
          if(modulo(ie-1,2) == 0) then
             x_val = -1.0 * (ie-1) * (ie-1)
          else
             x_val = 1.0 * (ie-1) * (ie-1)
          endif
       else
          x_val = ( (ie-1) * 1.0 / sqrt(1-x_arg**2)) * sin((ie-1) * b_arg)
       endif
    else
       print *, 'Invalid dtype in get_basis_val'
       stop
    endif

  end subroutine get_basis_val

  real*8 function myenefunc(x)
    ! Uses global iebasis, jebasis, ietype 
    use neo_globals, only : energy_max
    implicit none
    real :: x, xa, xb
    real :: ene, de, bi, bj, val
    real :: nu_d, nu_s, nu_par, nu_e, nu_h, nu_k, nu_p
    
    ! x grid: -1 -> 1 (equally spaced) for int 0..emax
    ! x = 2 (v/vth)/sqrt(2*emax)-1

    xa = 2.0 / (1.0 - sqrt(energy_min / energy_max))
    xb = - (1.0 + sqrt(energy_min / energy_max)) &
         / (1.0 - sqrt(energy_min / energy_max))
    ene = energy_max * ( (x-xb) / xa )**2 
    de = 2.0 * sqrt(energy_max) / xa

    if(iebasis <= 0) then
       bi = 1.0                                ! no B_i
    else
       if(ietype == 31) then
          call get_basis_val(iebasis,x,bi,1)   ! deriv of B_i
          bi = bi / (de * sqrt(ene))
       else
          call get_basis_val(iebasis,x,bi,0)   ! B_i
       endif
    endif
    
    if(jebasis <= 0) then
       bj = 1.0                                ! no B_j
    else
       if(ietype == 31 .or. ietype == 7 .or. ietype == 8 .or. ietype == 9 &
            .or. ietype == 10) then
          call get_basis_val(jebasis,x,bj,1)   ! deriv of B_j
          bj = bj / (de * sqrt(ene))
       else
          call get_basis_val(jebasis,x,bj,0)   ! B_j
       endif
    endif
    
    val = de * exp(-ene) * bi * bj  ! weight w/o ene factor
    
    ! Integral types
    if(ietype == 0) then
       val = val * ene                   ! ene^0
    else if(ietype == 1) then
       val = val * ene * ene             ! ene
    else if(ietype == 2) then
       val = val * ene * ene**2          ! ene^2
    else if(ietype == 3) then
       val = val * ene * ene**3          ! ene^3
    else if(ietype == 4) then
       val = val * ene * sqrt(ene)       ! sqrt(ene)
    else if(ietype == 5) then 
       val = val * ene * ene**1.5        ! ene^(3/2)
    else if(ietype == 6) then
       val = val * sqrt(ene)             ! 1.0/sqrt(ene)   
    else if(ietype == 7) then 
       val = val * ene * ene**2          ! ene^2 (*de)
    else if(ietype == 8) then 
       val = val * ene * ene**3          ! ene^3 (*de)
    else if(ietype == 9) then 
       val = val * ene * sqrt(ene)       ! sqrt(ene) (*de)   
    else if(ietype == 10) then 
       val = val * ene * ene             ! ene (*de)

    ! Collision Integral types
    else
       call get_coll_freqs(1,ir,is,js,ene,&
            nu_d, nu_s, nu_par, nu_e, nu_h, nu_k, nu_p)
       if(ietype == 20) then
          val = val * nu_d
       else if(ietype == 21) then
          val = val * nu_s
       else if(ietype == 22) then
          val = val * nu_s * sqrt(ene) 
       else if(ietype == 23) then
          val = val * nu_s * ene
       else if(ietype == 24) then
          val = val * nu_h * ene**1.5
       else if(ietype == 25) then
          val = val * nu_h * ene**3
       else if(ietype == 26) then
          val = val * nu_k * ene**1.5
       else if(ietype == 27) then
          val = val * nu_k * ene**3 
       else if(ietype == 28) then
          val = val * nu_p * ene
       else if(ietype == 29) then
          val = val * nu_p * ene**2
       else if(ietype == 30) then
          val = val * nu_e
       else if(ietype == 31) then
          val = val * nu_par * ene**2
       else if(ietype == 32) then
          val = val * (2.0*nu_d + nu_par) * ene
       else if(ietype == 33) then
          val = val * (2.0*nu_d + nu_par) * ene**2   
       else if(ietype == 34) then
          val = val * nu_s * ene**2 
       else
          print *, 'invalid itype in energy_grid module'
          stop
       end if
    end if
    
    myenefunc = val 

  end function myenefunc

  subroutine ENERGY_basis_ints

    use neo_globals, only : n_energy, write_out_mode

    implicit none
    integer :: ie, je, ke, Nx
    logical :: first
    real :: eii_val, peii_val

    if(write_out_mode > 1) then
       print *, 'start basis_ints'
    endif
    
    ! symmetric matrices
    do ietype = 0,6
       ke = 1
       do ie=1,n_energy
          do je=ie,n_energy

             ! Get the num int points from the B_ie * B_je integral convergence
             
             if(ie == je) then
                
                iebasis = ie
                jebasis = je
                Nx = Nx_start
                first  = .true.          
                int_basis_loop: do 
                   ! do the integral (uses iebasis, jebasis, ietype)
                   call gauss_integ(-1.0,1.0,myenefunc,&
                        integ_order,Nx,eii_val)
                   
                   if (first) then
                      peii_val  = eii_val
                      first = .false.
                      
                   else
                      ! check if converged
                      if(abs((eii_val-peii_val)/peii_val) < itol) then
                         exit int_basis_loop
                      endif
                   endif
                   
                   ! if not connverged, refine the grid and repeat
                   peii_val = eii_val
                   Nx = Nx * 2
                   if(Nx > Nx_max) then
                      print *, 'Energy int not converging for ietype', &
                           ietype
                      print *, 'ie = ', ie
                      stop
                   endif
                   
                end do int_basis_loop
                
             endif ! end if ie==je
             
             ! use the converged Nx to compute the vectors and matrics
             
             ! Matrices
             if (ietype==0 .or. ietype == 1 .or. ietype == 2 .or. &
                  ietype == 4 .or. ietype == 6) then
                iebasis = ie; jebasis = je
                call gauss_integ(-1.0,1.0,myenefunc,integ_order,Nx,eii_val)
                if(ietype == 0) then
                   emat_e0(ke)   = eii_val
                else if(ietype == 1) then
                   emat_e1(ke)   = eii_val
                else if(ietype == 2) then
                   emat_e2(ke)   = eii_val
                else if(ietype == 4) then
                   emat_e05(ke)  = eii_val
                else if(ietype == 6) then
                   emat_en05(ke)  = eii_val
                endif
             endif
             
             ! Vectors 
             if (ie == je) then
                iebasis = ie; jebasis = -1
                call gauss_integ(-1.0,1.0,myenefunc,integ_order,Nx,eii_val)
                if(ietype == 0) then
                   evec_e0(ie)   = eii_val
                else if(ietype == 1) then
                   evec_e1(ie)   = eii_val
                else if(ietype == 2) then
                   evec_e2(ie)   = eii_val
                else if(ietype == 3) then
                   evec_e3(ie)   = eii_val
                else if(ietype == 4) then
                   evec_e05(ie)  = eii_val
                else if(ietype == 5) then
                   evec_e105(ie) = eii_val
                endif
             endif
             
             ke = ke + 1
          enddo 
       enddo  
    enddo
    
    ! non-symmetric
    ! base convergence on ie=je and then loop over ie
    do ietype = 7,10
       do je=1,n_energy
          do ie=1,n_energy
             
             if(je .ne. 1) then
                iebasis = ie
                jebasis = je
                Nx = Nx_start
                first  = .true.          
                int_basis_loop2: do 
                   ! do the integral (uses iebasis, jebasis, ietype)
                   call gauss_integ(-1.0,1.0,myenefunc,&
                        integ_order,Nx,eii_val)
                   
                   if (first) then
                      peii_val  = eii_val
                      first = .false.
                      
                   else
                      ! check if converged
                      if(abs((eii_val-peii_val)/peii_val) < itol) then
                         exit int_basis_loop2
                      endif
                   endif
                   
                   ! if not connverged, refine the grid and repeat
                   peii_val = eii_val
                   Nx = Nx * 2
                   if(Nx > Nx_max) then
                      print *, 'Energy int not converging for ietype', &
                           ietype
                      print *, 'ie = ', ie
                      stop
                   endif
                   
                end do int_basis_loop2
                
             endif

             if(ietype == 7) then
                ! Matrix
                if(je == 1) then
                   ematij_e2de(ie,je) = 0.0
                else 
                   iebasis = ie; jebasis = je
                   call gauss_integ(-1.0,1.0,myenefunc,integ_order,&
                        Nx,eii_val)
                   ematij_e2de(ie,je) = eii_val
                endif
                ! Vector
                if(ie == je) then
                   if(je == 1) then
                      evec_e2de(je) = 0.0
                   else
                      iebasis = -1; jebasis = je
                      call gauss_integ(-1.0,1.0,myenefunc,integ_order,&
                           Nx,eii_val)
                      evec_e2de(je) = eii_val
                   endif
                endif
                
             else if(ietype == 8) then
                ! Vector
                if (ie == je) then
                   if(je == 1) then
                      evec_e3de(je) = 0.0
                   else
                      iebasis = -1; jebasis = je
                      call gauss_integ(-1.0,1.0,myenefunc,integ_order,&
                           Nx,eii_val)
                      evec_e3de(je) = eii_val
                   endif
                endif
                
             else if(ietype == 9) then
                ! Matrix
                if(je == 1) then
                   ematij_e05de(ie,je) = 0.0
                else
                   iebasis = ie; jebasis = je
                   call gauss_integ(-1.0,1.0,myenefunc,integ_order,&
                        Nx,eii_val)
                   ematij_e05de(ie,je) = eii_val
                endif

             else if(ietype == 10) then
                ! Matrix
                if(je == 1) then
                   ematij_e1de(ie,je) = 0.0
                else
                   iebasis = ie; jebasis = je
                   call gauss_integ(-1.0,1.0,myenefunc,integ_order,&
                        Nx,eii_val)
                   ematij_e1de(ie,je) = eii_val
                endif
                ! Vector
                if(ie == je) then
                   if(je == 1) then
                      evec_e1de(je) = 0.0
                   else
                      iebasis = -1; jebasis = je
                      call gauss_integ(-1.0,1.0,myenefunc,integ_order,&
                           Nx,eii_val)
                      evec_e1de(je) = eii_val
                   endif
                endif
                
             endif
             
          enddo
       enddo
    enddo
    
    if(write_out_mode > 1) then
       print *, 'done basis_ints'
    endif

  end subroutine ENERGY_basis_ints


  ! basis integrals for collision terms
  subroutine ENERGY_coll_ints(ir_in)

    use neo_globals, only : n_energy, n_species, case_spitzer, &
         collision_model, write_out_mode
    implicit none
    integer, intent(in) :: ir_in
    integer :: ie, je, ke, Nx, ietype_end, case_nofullhs
    logical :: first
    real :: eii_val, peii_val
    
    if(write_out_mode > 1) then
       print *, 'start coll_ints'
    end if

    ir = ir_in

    do is=1, n_species
       do js=1, n_species
          
          if (collision_model == 1 .or. collision_model == 2 &
               .or. ( (case_spitzer) .and. js .ne. is)) then
             ietype_end = 23
             emat_coll_rpi(is,js,:) = 0.0
             emat_coll_rn(is,js,:) = 0.0
             evec_coll_rh(is,js,:)  = 0.0
             econ_coll_rh(is,js)     = 0.0
             evec_coll_rk(is,js,:)  = 0.0
             econ_coll_rk(is,js)     = 0.0
             evec_coll_rp(is,js,:)  = 0.0
             econ_coll_rp(is,js)     = 0.0
             evec_coll_rq(is,js,:)  = 0.0
             econ_coll_rq(is,js)     = 0.0
             evec_coll_rqd(is,js,:) = 0.0
             evec_coll_rn_krook(is,js,:) = 0.0
             econ_coll_rn_krook(is,js) = 0.0
             case_nofullhs = 1
          else
             ietype_end = 34
             case_nofullhs = 0
          end if

          do ietype = 20, ietype_end

             ke = 1
             do ie=1, n_energy
                do je=ie, n_energy
                   
                   ! Get the num int points from the 
                   ! B_ie * B_je integral convergence
                   
                   ! Special cases: ietype=31
                   ! B_prime(ie=1) = 0 so all ints are zero
                   ! do not do convergence test and manually set
                   ! matrices and vectors to zero below
                   
                   if(ie == je .and. (.not.(ietype==31 .and. ie==1)) ) then
                      iebasis = ie
                      jebasis = je
                      
                      Nx = Nx_start
                      first  = .true.
                      int_coll_loop: do 
                         ! do the integral (uses iebasis, jebasis, ietype)
                         ! also uses (ir,is,js) for coll freqs
                         call gauss_integ(-1.0,1.0,myenefunc,&
                              integ_order,Nx,eii_val)
                         
                         if(first) then
                            if(abs(eii_val) < epsilon(0.)) then
                               print *, 'Coll convergence test failing'
                               print *, ietype, ir, is, js, ie, eii_val
                               stop
                            endif
                            peii_val  = eii_val
                            first = .false.
                            
                         else
                            ! check if converged
                            if(abs((eii_val-peii_val)/peii_val) < itol) then
                               exit int_coll_loop
                            endif
                         endif
                         ! if not converged, refine the grid and repeat
                         peii_val = eii_val
                         Nx = Nx * 2
                         if(Nx > Nx_max) then
                            print *, 'Energy int not converging'
                            print *, ietype, ir, is, js
                            stop
                         endif
                      end do int_coll_loop
                      
                   endif ! end if ie==je
                   
                   ! use the converged Nx to compute the vectors and matrics
                   
                   ! Matrices
                   if(ietype == 20 .or. ietype == 21 .or. &
                        ietype == 30 .or. ietype == 31) then
                      iebasis = ie; jebasis = je
                      call gauss_integ(-1.0,1.0,myenefunc,&
                           integ_order,Nx,eii_val)
                      if(ietype == 20) then
                         emat_coll_lorentz(is,js,ke) = eii_val
                      else if(ietype == 21) then
                         ! nud - nus
                         emat_coll_ru(is,js,ke) = &
                              emat_coll_lorentz(is,js,ke) - eii_val
                      else if(ietype == 30) then
                         ! -nue
                         emat_coll_rpi(is,js,ke) = -eii_val
                      else if(ietype == 31) then
                         if(ie == 1) then
                            ! EAB: B_prime(ie=1) = 0 (exactly)
                            emat_coll_rn(is,js,ke) = 0.0
                         else
                            emat_coll_rn(is,js,ke) = eii_val
                         endif
                      endif
                   endif
 
                   ! Vectors and constants
                   if(ie == je) then
                      
                      ! Vectors
                      if(ietype == 22 .or. ietype == 24 .or. &
                           ietype == 26 .or. ietype == 28 .or. &
                           ietype == 31 .or. ietype == 32 .or. &
                           ietype == 23 .or. ietype == 30) then
                         iebasis = ie; jebasis = -1
                         call gauss_integ(-1.0,1.0,myenefunc,&
                              integ_order,Nx,eii_val)
                         if(ietype == 22) then
                            evec_coll_rs(is,js,ie) = eii_val
                         else if(ietype == 24) then
                            evec_coll_rh(is,js,ie) = eii_val
                         else if(ietype == 26) then
                            evec_coll_rk(is,js,ie) = eii_val
                         else if(ietype == 28) then
                            evec_coll_rp(is,js,ie) = eii_val
                         else if(ietype == 31) then
                            if(ie == 1) then
                               ! EAB: B_prime(ie=1) = 0 (exactly)
                               evec_coll_rqd(is,js,ie) = 0.0
                            else
                               evec_coll_rqd(is,js,ie) = eii_val
                            endif
                         else if(ietype == 23) then
                            if(ie == 1 .or. case_nofullhs == 1) then
                               evec_coll_rq(is,js,ie) = 0.0
                            else 
                               evec_coll_rq(is,js,ie) = 2.0 * eii_val
                            endif
                         else if(ietype == 32) then
                            if(ie == 1 .or. case_nofullhs == 1) then
                               evec_coll_rq(is,js,ie) = 0.0
                            else 
                               evec_coll_rq(is,js,ie) = &
                                    evec_coll_rq(is,js,ie) - eii_val
                            endif
                         else if(ietype == 30) then
                            evec_coll_rn_krook(is,js,ie) = eii_val
                         endif
                      endif
                      
                      ! Constants
                      if(ie== 1 .and. (ietype == 23 .or. ietype == 25 &
                           .or. ietype == 27 .or. ietype == 29 &
                           .or. ietype == 33 .or. ietype == 34 &
                           .or. ietype == 30)) then
                         iebasis = -1; jebasis = -1
                         call gauss_integ(-1.0,1.0,myenefunc,&
                              integ_order,Nx,eii_val)
                         if(ietype == 23) then
                            econ_coll_rs(is,js) = eii_val
                         else if(ietype == 25) then
                            econ_coll_rh(is,js) = eii_val
                         else if(ietype == 27) then
                            econ_coll_rk(is,js) = eii_val
                         else if(ietype == 29) then
                            econ_coll_rp(is,js) = eii_val 
                         else if(ietype == 33) then
                            econ_coll_rq(is,js) = -eii_val
                         else if(ietype == 34) then
                            econ_coll_rq(is,js) = &
                                 econ_coll_rq(is,js) + 2.0*eii_val   
                         else if(ietype == 30) then
                            econ_coll_rn_krook(is,js) = eii_val
                         endif
                      endif
                         
                   endif
                   
                   ke = ke + 1
                enddo
             enddo 
          enddo    
          
       enddo       ! js
    enddo          ! is
    
    if(write_out_mode > 1) then
       print *, 'done coll_ints'
    end if
    
  end subroutine ENERGY_coll_ints
  
  
  ! Compute the energy-dependent coll freqs for a given ir, is, js
  ! returns nu -> (nu * ene)
  subroutine get_coll_freqs(flag, ir_loc,is_loc,js_loc,ene, &
       nu_d, nu_s, nu_par, nu_e, nu_h, nu_k, nu_p)
    use neo_globals
    implicit none
    real, external :: derf
    integer, intent(in)  :: flag
    ! 0->use HS0 model; else->use collision_model
    integer, intent(in) :: ir_loc, is_loc, js_loc
    real, intent(in)    :: ene
    real, intent(inout) :: nu_d, nu_s, nu_par, nu_e, nu_h, nu_k, nu_p
    real :: xa, xb, Hd_coll, Xd_coll, Hs_coll, Xs_coll, fac

    ! (Note: nu(is,ir) and pol part of dens from rotation will be added 
    !  in coll term in kinetic equation)
    fac = (1.0*Z(js_loc))**2 / (1.0*Z(is_loc))**2 & 
         * dens(js_loc,ir_loc)/dens(is_loc,ir_loc)
    xa = sqrt(ene)
    xb = xa * (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))

    if (flag .ne. 0 .and. &
         (collision_model == 1 .or. &
         (case_spitzer .and. js_loc .ne. is_loc)) ) then
       ! Connor model: nu_d = nu_s, form dependent on (is,js)
       if ((is_loc == js_loc) .or. &  
            (abs(mass(is_loc) - mass(js_loc)) < epsilon(0.))) then
          ! case 1: like-species/same mass collisions
          if(xb < 1e-4) then
             ! special case limit ene->0 
             Hd_coll = (1.0/sqrt(pi)) * &
                  (4.0/3.0    *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc)) &
                  - 4.0/15.0  *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**3 &
                  * ene & 
                  + 2.0/35.0  *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**5 &
                  * ene**2 &
                  - 2.0/189.0 *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**7 &
                  * ene**3)
             Xd_coll = 1.0
          else
             Hd_coll = exp(-xb*xb) / (xb*sqrt(pi)) &
                  + (1.0-1.0/(2.0*xb*xb)) * DERF(xb)
             Xd_coll = 1.0 / (xa)
          endif
       else if(mass(is_loc) < mass(js_loc)) then
          ! case 2: ele-ion and ion-imp(heavy) collisions
          Hd_coll = 1.0
          Xd_coll = 1.0 / (xa)
       else
          ! case 3: ion-ele and imp(heavy)-ion collisions
          Hd_coll = 1.0
          Xd_coll = 1.0 * (xa**2)
          fac = fac * 4.0/(3.0*sqrt(pi)) &
               * sqrt(mass(js_loc) / mass(is_loc)) &
               * (temp(is_loc,ir_loc) / temp(js_loc,ir_loc))**1.5
       endif
       nu_d = fac * Hd_coll * Xd_coll
       nu_s = nu_d
       
    else if(xb < 1e-4) then
       ! Hirshman-Sigmar model: nu_d != nu_s
       ! special case limit ene->0 
       ! return nu_d * ene (finite); nu_s * ene (zero)

       nu_d = fac * (1.0/sqrt(pi)) * &
            (4.0/3.0    * (vth(is_loc,ir_loc)/vth(js_loc,ir_loc)) &
            - 4.0/15.0  * (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**3 * ene & 
            + 2.0/35.0  * (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**5 * ene**2 &
            - 2.0/189.0 * (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**7 * ene**3)

       nu_par = fac * 2.0 * (1.0/sqrt(pi)) * &
            (2.0/3.0   *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc)) & 
            - 2.0/5.0  *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**3 * ene &
            + 1.0/7.0  *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**5 * ene**2 &
            - 1.0/27.0 *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**7 * ene**3)

       nu_s =  (2.0*temp(is_loc,ir_loc)/temp(js_loc,ir_loc)) &
            * (1.0 + mass(js_loc) / mass(is_loc)) &
            * fac * (1.0/sqrt(pi)) * &
            (2.0/3.0   *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc)) * ene & 
            - 2.0/5.0  *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**3 * ene**2 &
            + 1.0/7.0  *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**5 * ene**3)
       
    else
       ! Hirshman-Sigmar model: nu_d != nu_s
       Hd_coll = exp(-xb*xb) / (xb*sqrt(pi)) & 
            + (1.0-1.0/(2.0*xb*xb)) * DERF(xb)
       Xd_coll = 1.0 / (xa)
       Hs_coll = -exp(-xb*xb) / (xb*sqrt(pi)) & 
            + (1.0/(2.0*xb*xb)) * DERF(xb)
       Xs_coll = (xa) * (2.0 * temp(is_loc,ir_loc)/temp(js_loc,ir_loc)) &
            * (1.0 + mass(js_loc) / mass(is_loc))
       nu_d = fac * Hd_coll * Xd_coll
       nu_s = fac * Hs_coll * Xs_coll
       nu_par = fac * Hs_coll * Xd_coll * 2.0

    endif

    ! heating friction/energy diffusion frequencies
    ! (only for the full Hirshman-Sigmar model)
    if (flag == 0 &
         .or. collision_model == 1 .or. collision_model == 2 &
         .or. ( (case_spitzer) .and. js_loc .ne. is_loc)) then
       nu_par = 0.0
       nu_e   = 0.0
       nu_h   = 0.0
       nu_k   = 0.0
       nu_p   = 0.0
       
    else
       nu_e   = 2.0 * nu_s - 2.0 * nu_d - nu_par
       nu_h   = 1.5 * (nu_par &
            * (2.0*mass(is_loc)/mass(js_loc)-1) + 2.0 * nu_s &
            - 2.0 * (1.0+mass(is_loc)/mass(js_loc)) * nu_d )
       nu_k   = 2.0 * nu_par - nu_e
       nu_p   = 2.0 * nu_s &
            - 2.0 * (1.0+1.5*mass(is_loc)/mass(js_loc)) &
            * (nu_d - nu_par)
    endif
    
  end subroutine get_coll_freqs
  
end module neo_energy_grid
