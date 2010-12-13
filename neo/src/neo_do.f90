!-----------------------------------------------------------------
! neo_do.f90
!
! PURPOSE:
!  Subroutinized main neo program.  
!
! NOTES:
!  This can be called directly using the driver routine neo 
!  (in which case input data will read from input.dat) or called 
!  as a subroutine using neo_sub.
!-----------------------------------------------------------------

subroutine neo_do

  use neo_globals
  use neo_energy_grid
  use neo_sparse_solve
  use neo_equilibrium
  use neo_rotation
  use neo_transport
  use neo_theory
  use neo_g_velocitygrids
  use neo_allocate_profile
  implicit none

  integer :: i_iter, ir, is, ie, ix, it, jr, js, je, jx, jt, ks, ke
  integer, dimension(:,:,:,:), allocatable :: mindx  ! (ns,ne,nxi+1,nth)
  integer :: n_elem
  integer :: i, j, k, id
  logical :: temps_equal
  integer :: ierr

  ! higher-order terms:
  real, dimension(:,:),   allocatable :: g_global, g_old_global
  real, dimension(:,:,:), allocatable :: flux_drift
  real, dimension(:,:),   allocatable :: pflux_global, eflux_global
  real, dimension(:,:),   allocatable :: dpflux, deflux
  real, dimension(:),     allocatable :: spflux, seflux

  ! kinetic equation terms: 
  real :: stream, trap, rotkin
  
  integer, parameter :: io_neo=10, io_f=11, io_src=12

  call neo_make_profiles
  call neo_check

  ! cyclic index (for theta-periodicity)
  allocate(thcyc(1-n_theta:2*n_theta))
  do it=1,n_theta
     thcyc(it-n_theta) = it
     thcyc(it) = it
     thcyc(it+n_theta) = it
  enddo
  ! coefficients for 4th order centered derivative
  cderiv(-2) =  1
  cderiv(-1) = -8
  cderiv(0)  =  0
  cderiv(1)  =  8
  cderiv(2)  = -1
  
  if(sim_model == 0) then
     ! Theory calculation only -- no numerical kinetic calculation
     call EQUIL_alloc(1)
     call THEORY_alloc(1)
     do ir=1, n_radial
        call EQUIL_do(ir)
        call THEORY_do(ir)
        ! Store the local neo transport values at ir=1 in neo_xtheory_out
        ! (n_species_max, transport coeff)
        ! transport coeff: 1-> gamma, 2-> Q, 3->Pi, 4-> upar
        if(ir == 1) then
           neo_dke_out(:,:) = 0.0
           neo_dke_1d_out   = 0.0
           neo_gv_out(:,:)  = 0.0
           neo_th_out(:)    = 0.0
           neo_th_out(1)    = pflux_HH
           neo_th_out(2)    = efluxi_HH
           do is=1, n_species
              if(Z(is) == -1) then
                 neo_th_out(3) = efluxe_HH
              endif
           enddo
           neo_th_out(4)     = efluxi_CH1
           neo_th_out(5)     = jpar_HH
           neo_th_out(6)     = jpar_S
           neo_thHS_out(:,:) = 0.0
           do is=1, n_species
              neo_thHS_out(is,1) = pflux_multi_HS2(is) 
              neo_thHS_out(is,2) = eflux_multi_HS2(is) 
           enddo
        end if
     enddo
     call THEORY_alloc(0)
     call EQUIL_alloc(0)
     deallocate(thcyc)
     return
  end if

  ! Set-up the energy grids for basis ints (indep of r)
  call ENERGY_basis_ints_alloc(1)
  call ENERGY_basis_ints
  call ENERGY_coll_ints_alloc(1)

  ! Matrix solve allocations
  n_row = n_species*n_energy*(n_xi+1)*n_theta
  !n_max = n_species*n_energy*n_energy*(n_xi+1)*n_theta*350
  n_max = n_species*n_energy*n_energy*(n_xi+1)*n_theta*450
  allocate(a(n_max),stat=ierr)
  if(ierr /= 0) then
     print *, 'allocation failed'
     stop
  end if
  allocate(a_indx(2*n_max),stat=ierr)
  if(ierr /= 0) then
     print *, 'allocation failed'
     stop
  end if
  allocate(g(n_row))
  allocate(mindx(n_species,n_energy,0:n_xi,n_theta))
  allocate(is_indx(n_row))
  allocate(ie_indx(n_row))
  allocate(ix_indx(n_row))
  allocate(it_indx(n_row))

  ! matrix indices
  i = 0
  do is=1,n_species
     do ie=1,n_energy
        do ix=0,n_xi
           do it=1,n_theta
              i = i+1
              mindx(is,ie,ix,it) = i
              is_indx(i) = is
              ie_indx(i) = ie
              ix_indx(i) = ix
              it_indx(i) = it
           enddo
        enddo
     enddo
  enddo

  allocate(driftx(n_species,n_theta))
  allocate(driftth(n_species,n_theta))
  allocate(driftx_rot0(n_species,n_theta))
  allocate(driftx_rot1(n_species,n_theta))
  allocate(driftx_rot2(n_species,n_theta))

  call EQUIL_alloc(1)
  call ROT_alloc(1)
  call TRANSP_alloc(1)

  ! Set-up the matrix equation: LHS radially local matrix
  if (case_spitzer .or. zf_model == 1) then
     n_order = 1  ! only do local problem
  endif
  
  if(write_out_mode > 0) then
     open(unit=io_neo,file='grid.out',status='replace')
     write(io_neo,*) n_species
     write(io_neo,*) n_energy
     write(io_neo,*) n_xi
     write(io_neo,*) n_theta
     do it=1, n_theta
        write(io_neo,*) theta(it)
     enddo
     write(io_neo,*) n_radial
     do ir=1, n_radial
        write(io_neo,*) r(ir)
     enddo
     close(io_neo)
     
     open(unit=io_f,file='f.out',status='replace')
  end if

  if (n_order > 1) then
     if(write_out_mode > 0) then
        open(unit=io_src,file='source.out',status='replace')
     end if
     allocate(g_global(n_radial,n_row))
     allocate(g_old_global(n_radial,n_row))
     allocate(flux_drift(n_species,n_energy,0:n_xi))
     allocate(pflux_global(n_species,n_radial))
     allocate(eflux_global(n_species,n_radial))
     allocate(dpflux(n_species,n_radial))
     allocate(deflux(n_species,n_radial))
     allocate(spflux(n_species))
     allocate(seflux(n_species))
  endif

  do i_iter = 1, n_order
     if(write_out_mode > 1) print *, 'i_iter = ', i_iter

     do ir=1, n_radial
 
        if(write_out_mode > 1) print *, 'ir = ', ir 
        if(write_out_mode > 1) print *, 'Begin matrix set-up'
        
        temps_equal = .true.
        do is=2, n_species
           if(abs(temp(is,ir)-temp(1,ir)) > 1e-3) then
              temps_equal = .false.
              exit
           endif
        enddo

        ! Get the equilibrium parameters (th)
        call EQUIL_do(ir)
        
        ! Get the rotation phi
        call ROT_solve_phi(ir)
        
        ! Get the energy grids for coll ints
        call ENERGY_coll_ints(ir)
        
        ! Set the LHS (same form for each i_iter order)

        a_indx(:) = 0.0
        a(:) = 0.0
        k = 0
        
        do is=1,n_species
           do ie=1,n_energy
              do ix=0,n_xi
                 do it=1,n_theta
                    
                    ! Set the kinetic terms (funcs of is,ir,it)
                    if (case_spitzer) then
                       ! Spitzer problem: 
                       ! C_ee + nu_ei * L g_e 
                       ! = -v_par*(Z_e e E_par/T_e - grad_par ln p_e 
                       !           + (ene-5/2) grad_par ln T_e
                       stream   = 0.0
                       trap     = 0.0
                       rotkin   = 0.0
                       driftx(is,it)  = 0.0
                       driftth(is,it) = 0.0
                       driftx_rot0(is,it)  = 0.0
                       driftx_rot1(is,it)  = 0.0
                       driftx_rot2(is,it)  = 0.0
                       
                    else
                       ! vpar bhat dot grad 
                       ! -> stream * xi * sqrt(ene)
                       stream = sqrt(2.0) * vth(is,ir) * k_par(it) &
                            / (12*d_theta)
                       ! mu bdot grad B d/dvpar 
                       ! -> trap * (1-xi^2)*d/dxi * sqrt(ene)
                       trap   = ( gradpar_Bmag(it) / Bmag(it) ) &
                            * sqrt(0.5) * vth(is,ir)
                       ! vdrift dot grad r 
                       ! -> driftx * (1+xi^2) * d/d(r/a) * ene(ie)
                       driftx(is,it) = v_drift_x(it) &
                            * mass(is)/(1.0*Z(is)) * (vth(is,ir))**2 
                       ! r * vdrift dot grad th 
                       ! -> driftth * (1+xi^2) * d/d(th)* ene(ie)
                       driftth(is,it)  = v_drift_th(it) &
                            * mass(is)/(1.0*Z(is))* (vth(is,ir))**2 &
                            / (12*d_theta)
                       ! rotation terms
                       ! no rotation if zf problem
                       if(zf_model == 1) then
                          rotkin = 0.0
                          driftx_rot0(is,it) = 0.0
                          driftx_rot1(is,it) = 0.0
                          driftx_rot2(is,it) = 0.0
                       else
                          ! kinetic eqn rotation term
                          rotkin = sqrt(2.0) * vth(is,ir) *  k_par(it) & 
                               * (-Z(is)/temp(is,ir)*phi_rot_deriv(it) &
                               + omega_rot(ir)**2 * bigR(it) &
                               / vth(is,ir)**2 &
                               * bigR_tderiv(it)) 
                          ! vdrift dot grad r -- ix=0 rotation part
                          driftx_rot0(is,it) = I_div_psip* k_par(it) &
                               * mass(is)/(1.0*Z(is)) * rho(ir) / Bmag(it) &
                               * (vth(is,ir))**2 &
                               * (-Z(is)/temp(is,ir)*phi_rot_deriv(it) &
                               + omega_rot(ir)**2 * bigR(it)/ vth(is,ir)**2 &
                               * bigR_tderiv(it))
                          driftx_rot1(is,it) = I_div_psip* k_par(it) &
                               / Btor(it) &
                               * mass(is)/(1.0*Z(is)) * rho(ir) &
                               * vth(is,ir) * 2.0 * sqrt(2.0) &
                               * bigR_tderiv(it) * omega_rot(ir)
                          driftx_rot2(is,it) = 1.0/sqrt(2.0) &
                               * vth(is,ir)**2 * mass(is)/(1.0*Z(is)) &
                               * rho(ir) / Bmag(it) &
                               * Btor(it)/(Bmag(it)* I_div_psip) &
                               * (2.0*k_par(it) * gradr(it) &
                               * gradr_tderiv(it) &
                               - gradpar_Bmag(it) / Bmag(it) &
                               * gradr(it)**2)
                       endif
                    endif
                    
                    i = mindx(is,ie,ix,it)
                    
                    ! Impose constant for g0
                    ! local  -> at each species and ene 
                    
                    if (case_spitzer .and. ix==0) then
                       js = is; je = ie; jt = it
                       j = mindx(js,je,ix,jt)
                       k = k+1
                       a(k) = 1.0
                       a_indx(k) = i
                       a_indx(k+n_max) = j
                       
                    else if ((.not. case_spitzer) .and. ix==0 .and. it==1 &
                         .and. (collision_model==1 .or. collision_model==2))&
                         then
                       if(ie==1) then
                          ! <int f>=0
                          js = is
                          do je=1, n_energy
                             do jt=1, n_theta
                                j = mindx(js,je,ix,jt)
                                k = k+1
                                a(k) = w_theta(jt) * evec_e0(je)
                                a_indx(k) = i
                                a_indx(k+n_max) = j
                             enddo
                          enddo
                       else if(ie==2) then
                          ! <int e*f>=0
                          js = is
                          do je=1, n_energy
                             do jt=1, n_theta
                                j = mindx(js,je,ix,jt)
                                k = k+1
                                a(k) = w_theta(jt) * evec_e1(je)
                                a_indx(k) = i
                                a_indx(k+n_max) = j
                             enddo
                          enddo
                       else
                          ! <f_ie> = 0
                          js = is; je = ie
                          do jt=1, n_theta
                             j = mindx(js,je,ix,jt)
                             k = k+1
                             a(k) = w_theta(jt)
                             a_indx(k) = i
                             a_indx(k+n_max) = j
                          enddo
                       endif
                       
                    else if ((.not. case_spitzer) .and. ix==0 .and. it==1 &
                         .and. (ie==1 .or. &
                         (is==1 .and. ie==2 .and. temps_equal))) then
                       if(ie==1) then
                          ! <int f>=0
                          js = is
                          do je=1, n_energy
                             do jt=1, n_theta
                                j = mindx(js,je,ix,jt)
                                k = k+1
                                a(k) = w_theta(jt) * evec_e0(je)
                                a_indx(k) = i
                                a_indx(k+n_max) = j
                             enddo
                          enddo
                       else
                          ! sum_s <int e*f>=0
                          do js=1, n_species
                             do je=1, n_energy
                                do jt=1, n_theta
                                   j = mindx(js,je,ix,jt)
                                   k = k+1
                                   a(k) = w_theta(jt) * evec_e1(je)
                                   a_indx(k) = i
                                   a_indx(k+n_max) = j
                                enddo
                             enddo
                          enddo
                       endif
                       
                    else
                       
                       ! Collisions 
                       jt = it; jx = ix
                       if(ix == 0) then
                          ! only n_a1, q_ba energy diffusion terms
                          ! (no Lorentz term)
                          do je=1,n_energy
                             do js=1,n_species
                                j = mindx(js,je,jx,jt)
                                k = k+1
                                
                                ! q_ba term
                                if(abs(econ_coll_rq(is,js)) > epsilon(0.)) then
                                   if(collision_model == 4) then
                                      ! approx krook-like HS diffusion term
                                      a(k) = - (temp(js,ir)/temp(is,ir)) &
                                           * (dens(js,ir)/dens(is,ir)) &
                                           / econ_coll_rq(is,js) &
                                           * evec_coll_rq(is,js,ie) &
                                           * evec_coll_rq(js,is,je) &
                                           * dens_fac(js,it) * nu(js,ir)
                                   else
                                      ! full HS diffusion term
                                      a(k) = -2.0 * (temp(js,ir)/temp(is,ir)) &
                                           * (dens(js,ir)/dens(is,ir)) &
                                           / econ_coll_rq(is,js) &
                                           * evec_coll_rqd(is,js,ie) &
                                           * evec_coll_rq(js,is,je) &
                                           * dens_fac(js,it) * nu(js,ir)
                                   endif
                                else
                                   a(k) = 0.0
                                endif
                                
                                if(js == is) then
                                   if(collision_model == 4) then
                                      ! approx krook-like HS diffusion term
                                      do ks=1, n_species
                                         a(k) = a(k) &
                                              -emat_coll_rpi(is,ks,eindx(ie,je)) &
                                              * dens_fac(ks,it) * nu(is,ir)
                                         if(abs(econ_coll_rn_krook(is,ks)) &
                                              > epsilon(0.)) then
                                            a(k) = a(k) &
                                                 -1.0 &
                                                 *evec_coll_rn_krook(is,ks,ie) &
                                                 *evec_coll_rn_krook(is,ks,je) &
                                                 /econ_coll_rn_krook(is,ks) &
                                                 * dens_fac(ks,it) * nu(is,ir)
                                         endif
                                      enddo
                                   else
                                      ! full HS diffusion term
                                      do ks=1, n_species
                                         a(k) = a(k) + 2.0 * &
                                              emat_coll_rn(is,ks,eindx(ie,je)) &
                                              * dens_fac(ks,it) * nu(is,ir)
                                      enddo
                                   endif
                                endif
                                
                                a_indx(k) = i
                                a_indx(k+n_max) = j
                             enddo
                          enddo
                          
                       else if(ix == 1) then
                          ! Lorentz + slowing down + h,k heating terms
                          do je=1,n_energy
                             do js=1,n_species
                                j = mindx(js,je,jx,jt)
                                k = k+1 
                                
                                ! slowing-down term
                                a(k) = -(mass(js)/mass(is)) &
                                     * (vth(js,ir)/vth(is,ir)) & 
                                     * (dens(js,ir)/dens(is,ir)) &
                                     / econ_coll_rs(is,js) &
                                     * evec_coll_rs(is,js,ie) &
                                     * evec_coll_rs(js,is,je) &
                                     * dens_fac(js,it) * nu(js,ir)
                                
                                ! h heating term
                                if(abs(econ_coll_rh(is,js)) > epsilon(0.)) then
                                   a(k) = a(k) - (mass(js)/mass(is)) &
                                        * (vth(js,ir)/vth(is,ir))**3 &
                                        * (dens(js,ir)/dens(is,ir)) &
                                        / econ_coll_rh(is,js) &
                                        * evec_coll_rh(is,js,ie) &
                                        * evec_coll_rh(js,is,je) &
                                        * dens_fac(js,it) * nu(js,ir)
                                endif
                                
                                if (js == is) then  
                                   do ks=1, n_species
                                      ! Lorentz and u_a1 terms
                                      a(k) = a(k) &
                                           + 0.5*ix*(ix+1) &
                                           * emat_coll_lorentz(is,ks,eindx(ie,je)) &
                                           * dens_fac(ks,it) * nu(is,ir) &
                                           - emat_coll_ru(is,ks,eindx(ie,je)) &
                                           * dens_fac(ks,it) * nu(is,ir)
                                      ! k heating term
                                      if(abs(econ_coll_rk(is,ks)) &
                                           > epsilon(0.)) then
                                         a(k) = a(k) &
                                              - 1.0 / econ_coll_rk(is,ks) &
                                              * evec_coll_rk(is,ks,ie) &
                                              * evec_coll_rk(is,ks,je) &
                                              * dens_fac(ks,it) * nu(is,ir)
                                      endif
                                   enddo
                                endif
                                a_indx(k) = i
                                a_indx(k+n_max) = j
                             enddo
                          enddo
                          
                       else if (ix == 2) then
                          ! Lorentz + pi energy diffusion terms
                          do je=1,n_energy
                             do js=1,n_species
                                j = mindx(js,je,jx,jt)
                                k = k+1
                                
                                ! pi_ba term
                                if(abs(econ_coll_rp(is,js)) > epsilon(0.)) then
                                   a(k) = -(temp(js,ir)/temp(is,ir)) &
                                        * (dens(js,ir)/dens(is,ir)) &
                                        / econ_coll_rp(is,js) &
                                        * evec_coll_rp(is,js,ie) &
                                        * evec_coll_rp(js,is,je) &
                                        * dens_fac(js,it) * nu(js,ir)
                                else
                                   a(k) = 0.0
                                endif
                                
                                if (js == is) then  
                                   ! Lorentz and pi_a1 terms
                                   do ks=1, n_species
                                      a(k) = a(k) &
                                           + 0.5*ix*(ix+1) &
                                           * emat_coll_lorentz(is,ks,eindx(ie,je)) &
                                           * dens_fac(ks,it) * nu(is,ir) &
                                           - emat_coll_rpi(is,ks,eindx(ie,je)) &
                                           * dens_fac(ks,it) * nu(is,ir)
                                   enddo
                                endif
                                
                                a_indx(k) = i
                                a_indx(k+n_max) = j
                             enddo
                          enddo
                          
                       else
                          ! only Lorentz term
                          js = is
                          do je=1, n_energy
                             j = mindx(js,je,jx,jt)
                             k = k+1
                             a(k) = 0.0
                             do ks=1, n_species
                                a(k) = a(k) + 0.5*ix*(ix+1) &
                                     * emat_coll_lorentz(is,ks,eindx(ie,je)) &
                                     * dens_fac(ks,it) * nu(is,ir)
                             enddo
                             a_indx(k) = i
                             a_indx(k+n_max) = j
                          enddo
                       endif
                       
                       ! Streaming
                       js = is
                       do je=1, n_energy
                          do id=-2,2
                             if (id /= 0) then
                                jt = thcyc(it+id)
                                jx = ix-1
                                if (jx >= 0) then
                                   j = mindx(js,je,jx,jt)
                                   k = k+1
                                   a(k) = stream &
                                        * ix/(2*ix-1.0) * cderiv(id) &
                                        * emat_e05(eindx(ie,je))
                                   a_indx(k) = i
                                   a_indx(k+n_max) = j
                                endif
                                jx = ix+1
                                if (jx <= n_xi) then
                                   j = mindx(js,je,jx,jt)
                                   k = k+1
                                   a(k) = stream &
                                        * (ix+1.0)/(2*ix+3.0) * cderiv(id) &
                                        * emat_e05(eindx(ie,je))
                                   a_indx(k) = i
                                   a_indx(k+n_max) = j
                                endif
                             endif
                          enddo
                       enddo
                       
                       ! Trapping and Rotation
                       js = is; jt = it
                       do je=1, n_energy
                          jx = ix-1
                          if (jx >= 0) then
                             j = mindx(js,je,jx,jt)
                             k = k+1
                             a(k) = trap &
                                  * ix*(ix-1.0)/(2*ix-1.0) &
                                  * emat_e05(eindx(ie,je)) &
                                  + rotkin &
                                  * (ix/(2*ix-1.0) * ematij_e05de(ie,je) &
                                  - ix*(ix-1.0)/(2*ix-1.0) * 0.5 &
                                  * emat_en05(eindx(ie,je)))
                             a_indx(k) = i
                             a_indx(k+n_max) = j                      
                          endif
                          jx = ix+1
                          if (jx <= n_xi) then
                             j = mindx(js,je,jx,jt)
                             k = k+1
                             a(k) = -trap &
                                  * (ix+1.0)*(ix+2.0)/(2*ix+3.0) &
                                  * emat_e05(eindx(ie,je)) &
                                  + rotkin &
                                  * ((ix+1.0)/(2*ix+3.0) &
                                  * ematij_e05de(ie,je) &
                                  + (ix+1.0)*(ix+2.0)/(2*ix+3.0) * 0.5 &
                                  * emat_en05(eindx(ie,je)))
                             a_indx(k) = i
                             a_indx(k+n_max) = j
                          endif
                       enddo
                       
                    endif ! bc if/else
                 enddo ! it
              enddo ! ix
           enddo ! ie
        enddo ! is
        
        n_elem = k
        do k=1,n_elem
           a_indx(n_elem+k) = a_indx(n_max+k)
        enddo
        
        ! Factor the Matrix -- uses a(:) and a_indx(:)
        if(write_out_mode > 1) print *, 'Begin matrix factor'
        call SOLVE_factor(n_elem)
        if(write_out_mode > 1) print *, 'Done matrix factor'
        
        ! Set the RHS source (form depends on i_iter order)
        call set_RHS_source
        
        ! Matrix solve -- uses g(:), a(:), and a_indx(:)
        if(write_out_mode > 1) print *, 'Begin matrix solve'
        call SOLVE_do
        if(write_out_mode > 1) print *, 'Done matrix solve'

        ! Compute the neo transport coefficients
        call TRANSP_do(ir)
        call TRANSP_write(ir)

        ! re-construct the energy dependence
        call g_energy(i_iter,ir)

        if(i_iter == 1) then
           ! Write the rotation parameters
           call ROT_write(ir)

           ! Compute the theory transport coefficients
           if(ir == 1) then
              call THEORY_alloc(1)
           end if
           call  THEORY_do(ir)

           ! Store the local neo transport values at ir=1 in neo_x_out
           ! (n_species_max, transport coeff)
           ! transport coeff: 1-> gamma, 2-> Q, 3->Pi, 4-> upar
           if(ir == 1) then
              neo_dke_out(:,:) = 0.0
              neo_dke_1d_out   = 0.0
              neo_gv_out(:,:)  = 0.0
              do is=1, n_species
                 neo_dke_out(is,1) = pflux(is)
                 neo_dke_out(is,2) = eflux(is)
                 neo_dke_out(is,3) = mflux(is)
                 neo_dke_out(is,4) = eflux(is) - omega_rot(ir)*mflux(is)
                 neo_dke_out(is,5) = vpol_th0(is)
                 neo_dke_out(is,6) = vtor_th0(is) + vtor_0order_th0
                 neo_gv_out(is,1)  = pflux_gv(is)
                 neo_gv_out(is,2)  = eflux_gv(is)
                 neo_gv_out(is,3)  = mflux_gv(is)
                 neo_gv_out(is,4)  = eflux_gv(is) - omega_rot(ir)*mflux_gv(is)
              enddo
              neo_dke_1d_out    = jpar
              neo_th_out(:) = 0.0
              neo_th_out(1) = pflux_HH
              neo_th_out(2) = efluxi_HH
              do is=1, n_species
                 if(Z(is) == -1) then
                    neo_th_out(3) = efluxe_HH
                 endif
              enddo
              neo_th_out(4) = efluxi_CH1
              neo_th_out(5)     = jpar_HH
              neo_th_out(6)     = jpar_S
              neo_thHS_out(:,:) = 0.0
              do is=1, n_species
                 neo_thHS_out(is,1) = pflux_multi_HS2(is) 
                 neo_thHS_out(is,2) = eflux_multi_HS2(is) 
              enddo
           end if

           if(ir == n_radial) then
              call THEORY_alloc(0)
           end if

        endif

        if(write_out_mode > 0) then
           write(io_f,*) g(:)
        endif


        if(n_order > 1) then
           g_global(ir,:) = g(:)
           do is=1, n_species
              pflux_global(is,ir) = pflux(is)
              eflux_global(is,ir) = eflux(is)
           enddo
        endif

     end do ! ir

     if(n_order > 1) then
        g_old_global(:,:) = g_global(:,:)
        
        ! Diagnostics: 
        ! compute 1/vprime d/dr (vprime * flux) for previous iter
        do is=1, n_species
           do ir=1, n_radial
              dpflux(is,ir) = 0.0
              deflux(is,ir) = 0.0
              do jr=1, n_radial
                 dpflux(is,ir) =  dpflux(is,ir) & 
                      + wd_rad(ir,jr) * pflux_global(is,jr) &
                      * v_prime_g(jr) * dens_norm(jr) * vth_norm(jr) &
                      / (v_prime_g(ir) * dens_norm(ir) * vth_norm(ir))
                 deflux(is,ir) =  deflux(is,ir) & 
                      + wd_rad(ir,jr) * eflux_global(is,jr) &
                      * v_prime_g(jr) * dens_norm(jr) * vth_norm(jr) &
                      * temp_norm(jr) &
                      / (v_prime_g(ir) * dens_norm(ir) * vth_norm(ir) &
                      * temp_norm(ir))
              enddo
           enddo
        enddo
     endif

  end do    ! i_iter

  ! Clean-up
  if(write_out_mode > 0) then
     close(io_f)
  end if
  if(n_order > 1) then
     if(write_out_mode > 0) then
        close(io_src)
     end if
     deallocate(g_old_global); deallocate(g_global)
     deallocate(flux_drift)
     deallocate(pflux_global); deallocate(eflux_global)
     deallocate(spflux);       deallocate(seflux)
     deallocate(dpflux);       deallocate(deflux)
  endif
  call ENERGY_basis_ints_alloc(0)
  call ENERGY_coll_ints_alloc(0)
  call EQUIL_alloc(0)
  call ROT_alloc(0)
  call TRANSP_alloc(0)
  call PROFILE_SIM_alloc(0)
  deallocate(a)
  deallocate(a_indx)
  deallocate(g)
  deallocate(mindx)
  deallocate(is_indx)
  deallocate(ie_indx)
  deallocate(ix_indx)
  deallocate(it_indx)
  deallocate(thcyc)
  deallocate(driftx)
  deallocate(driftth)
  deallocate(driftx_rot0)
  deallocate(driftx_rot1)
  deallocate(driftx_rot2)

contains

  subroutine set_RHS_source
    implicit none
    real :: src_F0_Ln, src_F0_Lt, src_P0, src_Rot1, src_Rot2

    g(:) = 0.0

    if(i_iter == 1) then
       ! First-Order Source Term

       do is=1,n_species
          
          ! src = (1/F0) * dF0/dr + Ze/T dPhi0/dr
          src_F0_Ln = -(dlnndr(is,ir) - 1.5*dlntdr(is,ir))  ! ene^0 part
          src_F0_Lt = -dlntdr(is,ir)                        ! ene^1 part
          src_P0    = (1.0*Z(is))/temp(is,ir) * dphi0dr(ir) ! ene^0 Er part
          
          do ie=1,n_energy
             do ix=0,n_xi
                do it=1,n_theta
                   
                   i = mindx(is,ie,ix,it)
                   
                   src_Rot1   = -omega_rot(ir) * bigR_th0**2 &
                        / vth(is,ir)**2 * omega_rot_deriv(ir) &
                        - dlntdr(is,ir) * Z(is) / temp(is,ir) * phi_rot(it) &
                        + dlntdr(is,ir) * (omega_rot(ir)/vth(is,ir))**2 * 0.5 &
                        * (bigR(it)**2 - bigR_th0**2) &
                        - omega_rot(ir)**2 * bigR_th0 &
                        / vth(is,ir)**2 * bigR_th0_rderiv
                   
                   src_Rot2 = omega_rot_deriv(ir) * bigR(it) / vth(is,ir) 

                   if(case_spitzer) then
                      
                      if(ix == 0) then
                         ! Impose constant for g0
                         g(i) = 0.0                      
                      else if (ix == 1) then
                         ! Equilibrium source: 
                         ! -(vpar dot grad (F0 + Ze/T Phi0))
                         g(i) = -sqrt(2.0) * vth(is,ir) &
                              * ( (src_F0_Ln + src_P0) * evec_e05(ie) &
                              + src_F0_Lt * evec_e105(ie) )
                      endif
                      
                   else
                      
                      ! Impose constant for g0
                      
                      if(ix==0 .and. it==1 .and. &
                           (collision_model==1 .or. collision_model==2)) then
                         g(i) = 0.0
                         
                      else if(ix==0 .and. it==1 &
                           .and. (ie==1 .or. &
                           (is==1 .and. ie==2 .and. temps_equal))) then
                         g(i) = 0.0
                         
                      else 
                         ! Equilibrium source:
                         ! Ze/T (I/Omega d^2phi/drdt) vpar F0
                         ! - vpar b_dot_grad (Ze/T (I/Omega dphi/dr) vpar F0)
                       
                         if(zf_model == 1) then
                            ! Equilibrium source             
                            if(ix == 0) then
                               g(i) = 4.0/3.0 &
                                    * gradpar_Bmag(it) /  Bmag(it) &
                                    * I_div_psip * vth(is,ir)**2 &
                                    * rho(ir) &
                                    * mass(is)/(1.0*Z(is)) / Bmag(it) &
                                    * src_P0 * evec_e1(ie)
                            else if(ix == 1) then
                               g(i) = I_div_psip * sqrt(2.0) &
                                    * vth(is,ir) * rho(ir) &
                                    * mass(is)/(1.0*Z(is)) / Bmag(it) &
                                    * zf_time * src_P0 * evec_e05(ie)
                            else if(ix == 2) then
                               g(i) = 2.0/3.0 &
                                    * gradpar_Bmag(it) /  Bmag(it) &
                                    * I_div_psip * vth(is,ir)**2 &
                                    * rho(ir) &
                                    * mass(is)/(1.0*Z(is)) / Bmag(it) &
                                    * src_P0 * evec_e1(ie) 
                            endif
                          
                         else 
                            ! Equilibrium source:
                            ! -(vdrift dot grad (F0 + Ze/T Phi0))
                            
                            if (ix == 0) then
                               g(i) = -(4.0/3.0) * driftx(is,it) &
                                    * ( (src_F0_Ln + src_P0 + src_Rot1) &
                                    * evec_e1(ie) &
                                    + src_F0_Lt * evec_e2(ie) ) &
                                    - driftx_rot0(is,it) &
                                    * ( (src_F0_Ln + src_P0 + src_Rot1) &
                                    * evec_e0(ie) &
                                    + src_F0_Lt * evec_e1(ie) ) &
                                    - src_Rot2 * 1.0/3.0* driftx_rot1(is,it) &
                                    * sqrt(2.0) * Btor(it)/Bmag(it) &
                                    * evec_e1(ie) &
                                    - src_Rot2 * 4.0/3.0 * driftx(is,it) &
                                    * omega_rot(ir) * bigR(it)/vth(is,ir) &
                                    * evec_e1(ie) &
                                    - src_Rot2 * driftx_rot0(is,it) &
                                    * omega_rot(ir) * bigR(it)/vth(is,ir) &
                                    * evec_e0(ie)
                               
                            else if(ix == 1) then
                               g(i) =  - driftx_rot1(is,it) &
                                    * ( (src_F0_Ln + src_P0 + src_Rot1) &
                                    * evec_e05(ie) &
                                    + src_F0_Lt * evec_e105(ie) ) &
                                    - src_Rot2 * 8.0/5.0 * driftx(is,it) &
                                    * sqrt(2.0) * Btor(it)/Bmag(it) &
                                    * evec_e105(ie) &
                                    - src_Rot2 * driftx_rot0(is,it) &
                                    * sqrt(2.0) * Btor(it)/Bmag(it) &
                                    * evec_e05(ie) &
                                    - src_Rot2 * driftx_rot1(is,it) &
                                    * omega_rot(ir) * bigR(it)/vth(is,ir) &
                                    * evec_e05(ie) &
                                    - src_Rot2 * 2.0/5.0* driftx_rot2(is,it) &
                                    * evec_e105(ie)
                               g(i) = g(i) + sqrt(2.0) * vth(is,ir) &
                                    * (1.0*Z(is))/temp(is,ir) &
                                    * epar0(ir) * evec_e05(ie)

                            else if (ix == 2) then
                               g(i) = -(2.0/3.0) * driftx(is,it) &
                                    * ( (src_F0_Ln + src_P0 + src_Rot1) &
                                    * evec_e1(ie) &
                                    + src_F0_Lt * evec_e2(ie) ) &
                                    - src_Rot2 * 2.0/3.0* driftx_rot1(is,it) &
                                    * sqrt(2.0) * Btor(it)/Bmag(it) &
                                    * evec_e1(ie) &
                                    - src_Rot2 * 2.0/3.0 * driftx(is,it) &
                                    * omega_rot(ir) * bigR(it)/vth(is,ir) &
                                    * evec_e1(ie)
                            else if(ix == 3) then
                               g(i) = -src_Rot2 * 2.0/5.0 * driftx(is,it) &
                                    * sqrt(2.0) * Btor(it)/Bmag(it) &
                                    * evec_e105(ie) &
                                    + src_Rot2 * 2.0/5.0* driftx_rot2(is,it) &
                                    * evec_e105(ie)
                            endif
                            
                         endif
                         
                      endif

                   endif

                enddo
             enddo
          enddo
       enddo
  
    else
       ! higher-order RHS (i_iter >= 2)

       ! Compute the drift term and the source=<drift>
       flux_drift(:,:,:) = 0.0
       spflux(:) = 0.0
       seflux(:) = 0.0
       do is=1,n_species
          do ie=1,n_energy
             do ix=0,n_xi
                do it=1,n_theta
                   i = mindx(is,ie,ix,it)
                   call get_drift
                   flux_drift(is,ie,ix) = flux_drift(is,ie,ix) &
                        + w_theta(it) * g(i)
                enddo
             enddo
          enddo
       enddo

       ! RHS = - (drift - source)
       do is=1,n_species
          do ie=1,n_energy
             do ix=0,n_xi
                do it=1,n_theta
                   i = mindx(is,ie,ix,it)
                   ! Impose constant for g0
                   if(ix==0 .and. it==1 .and. &
                        (collision_model==1 .or. collision_model==2)) then
                      g(i) = 0.0
                   else if(ix==0 .and. it==1 &
                        .and. (ie==1 .or. &
                        (is==1 .and. ie==2 .and. temps_equal))) then
                      g(i) = 0.0
                   else if(ix==0) then
                      g(i) = -(g(i) - flux_drift(is,ie,ix))
                   else
                      g(i) = -(g(i))
                   endif
                enddo
             enddo
          enddo
       enddo
       
       ! Diagnostics: compute vel integral of source term
       do is=1, n_species
          spflux(is) = spflux(is) * dens(is,ir) * 2.0/sqrt(pi)
          seflux(is) = seflux(is) * dens(is,ir) * temp(is,ir) * 2.0/sqrt(pi)
       enddo
       
       
       if(write_out_mode > 0) then
          write (io_src,'(e16.8,$)') r(ir)
          do is=1, n_species
             write (io_src,'(e16.8,$)') spflux(is)
             write (io_src,'(e16.8,$)') dpflux(is,ir)
             write (io_src,'(e16.8,$)') seflux(is)
             write (io_src,'(e16.8,$)') deflux(is,ir)
          enddo
          write (io_src,*)
       end if

    endif
    
  end subroutine set_RHS_source

  ! Computes the drift term for a given (is,ie,ix,it,ir)
  ! -> stores in g(i); uses g_old_global(:,:)
  subroutine get_drift
    implicit none
    real, parameter :: driftx_const    =  1.0
    real, parameter :: drifte_const    =  1.0
    real, parameter :: driftth_const   =  1.0
    real, parameter :: drifter1_const  =  1.0
    real, parameter :: drifter2_const  =  1.0
    real, parameter :: drifterth_const =  1.0
    integer :: jr, kr, je, kt
    real :: rfac, rfac_spflux, rfac_seflux, val, src_F0_Ln, src_F0_Lt, src_P0
    real :: val1, val2

    ! Radial drift (from B dot grad B)
    
    js = is; jt = it
    
    do kr=1,n_radial+1
       
       do je=1, n_energy
          
          ! d(g*F0) = F0*dg + g * dF0 
          if(kr == n_radial+1) then
             ! rad and ene part of g * 1/F0 vdrift dot grad F0 
             jr = ir
             src_F0_Ln = -(dlnndr(is,ir) - 1.5*dlntdr(is,ir)) ! ene^0 part
             src_F0_Lt = -dlntdr(is,ir)                       ! ene^1 part
             rfac = src_F0_Ln * emat_e1(eindx(ie,je)) &
                  + src_F0_Lt * emat_e2(eindx(ie,je))
             if(ix==0 .and. je == ie) then
                rfac_spflux = src_F0_Ln * evec_e1(ie) &
                     + src_F0_Lt * evec_e2(ie)
                rfac_seflux = src_F0_Ln * evec_e2(ie) &
                     + src_F0_Lt * evec_e3(ie)
             endif
          else
             jr = kr
             rfac = wd_rad(ir,kr) * emat_e1(eindx(ie,je))
             if(ix==0 .and. je == ie) then
                rfac_spflux = wd_rad(ir,kr) * evec_e1(ie)
                rfac_seflux = wd_rad(ir,kr) * evec_e2(ie)
             endif
          end if

          jx = ix
          j = mindx(js,je,jx,jt)
          val = driftx(is,it) &
               * (1.0 + (ix+1.0)*(ix+1.0)/((2*ix+1.0)*(2*ix+3.0)) &
               + (ix*ix)/((2*ix+1.0)*(2*ix-1.0))) &
               * driftx_const * g_old_global(jr,j)
          g(i) = g(i) +  rfac * val    
          if(ix==0 .and. je==ie) then
             spflux(is) = spflux(is) + rfac_spflux * val &
                  * w_theta(it)
             seflux(is) = seflux(is) + rfac_seflux * val &
                  * w_theta(it)
          endif

          jx = ix-2
          if (jx >= 0) then
             j = mindx(js,je,jx,jt)
             val = driftx(is,it) &
                  * (ix)*(ix-1.0)/((2*ix-3.0)*(2*ix-1.0)) &
                  * driftx_const * g_old_global(jr,j)
             g(i) = g(i) +  rfac * val    
             if(ix==0 .and. je==ie) then
                spflux(is) = spflux(is) + rfac_spflux * val &
                     * w_theta(it)
                seflux(is) = seflux(is) + rfac_seflux * val &
                     * w_theta(it)
             endif
          endif

          jx = ix+2
          if (jx <= n_xi) then
             j = mindx(js,je,jx,jt)
             val = driftx(is,it) &
                  * (ix+1.0)*(ix+2.0)/((2*ix+3.0)*(2*ix+5.0)) &
                  * driftx_const * g_old_global(jr,j)
             g(i) = g(i) + rfac * val
             if(ix==0 .and. je==ie) then
                spflux(is) = spflux(is) + rfac_spflux * val &
                     * w_theta(it)
                seflux(is) = seflux(is) + rfac_seflux * val &
                     * w_theta(it)
             endif
          endif

       enddo
    enddo

    ! Energy deriv part of Radial drift (from B dot grad B)

    js = is; jt = it

    do je=1,n_energy

       jx = ix
       j = mindx(js,je,jx,jt)
       val = dlntdr(is,ir) * driftx(is,it) &
            * (1.0 + (ix+1.0)*(ix+1.0)/((2*ix+1.0)*(2*ix+3.0)) &
            + (ix*ix)/((2*ix+1.0)*(2*ix-1.0))) &
            * drifte_const * g_old_global(ir,j) 
       g(i) = g(i) + ematij_e2de(ie,je) * val
       if(ix==0 .and. je==ie) then
          spflux(is) = spflux(is) + evec_e2de(ie) * val &
               * w_theta(it)
          seflux(is) = seflux(is) + evec_e3de(ie) * val &
               * w_theta(it)
       endif

       jx = ix-2
       if (jx >= 0) then
          j = mindx(js,je,jx,jt)
          val = dlntdr(is,ir) * driftx(is,it) &
               * (ix)*(ix-1.0)/((2*ix-3.0)*(2*ix-1.0)) &
               * drifte_const * g_old_global(ir,j)
          g(i) = g(i) + ematij_e2de(ie,je) * val
          if(ix==0 .and. je==ie) then
             spflux(is) = spflux(is) + evec_e2de(ie) * val &
                  * w_theta(it)
             seflux(is) = seflux(is) + evec_e3de(ie) * val &
                  * w_theta(it)
          endif
       endif

       jx = ix+2
       if (jx <= n_xi) then
          j = mindx(js,je,jx,jt)
          val = dlntdr(is,ir) * driftx(is,it) &
               * (ix+1.0)*(ix+2.0)/((2*ix+3.0)*(2*ix+5.0)) &
               * drifte_const * g_old_global(ir,j)
          g(i) = g(i) +  ematij_e2de(ie,je) * val
          if(ix==0 .and. je==ie) then
             spflux(is) = spflux(is) + evec_e2de(ie) * val &
                  * w_theta(it)
             seflux(is) = seflux(is) + evec_e3de(ie) * val &
                  * w_theta(it)
          endif
       endif

    enddo

    ! Theta drift (from B dot grad B and grad cross B)

    js = is

    do je = 1,n_energy
       do id=-2,2
          if (id /= 0) then
             jt = thcyc(it+id)

             jx = ix
             j = mindx(js,je,jx,jt)
             rfac = (1.0/r(ir)) * g_old_global(ir,j)
             val = cderiv(id) * rfac &
                  * driftth_const * driftth(is,it) &
                  * (1.0 + (ix+1.0)*(ix+1.0)/((2*ix+1.0)*(2*ix+3.0)) &
                  + (ix*ix)/((2*ix+1.0)*(2*ix-1.0)))
             g(i) = g(i) + emat_e1(eindx(ie,je)) * val
             if(ix==0 .and. je==ie) then
                spflux(is) = spflux(is) + evec_e1(ie) * val &
                     * w_theta(it)
                seflux(is) = seflux(is) + evec_e2(ie) * val &
                     * w_theta(it)
             endif

             jx = ix-2
             if (jx >= 0) then
                j = mindx(js,je,jx,jt)
                rfac = (1.0/r(ir)) * g_old_global(ir,j)
                val = cderiv(id) * rfac &
                     * driftth_const * driftth(is,it) &
                     * (ix)*(ix-1.0)/((2*ix-3.0)*(2*ix-1.0))
                g(i) = g(i) +  emat_e1(eindx(ie,je)) * val
                if(ix==0 .and. je==ie) then
                   spflux(is) = spflux(is) + evec_e1(ie) * val &
                        * w_theta(it)
                   seflux(is) = seflux(is) + evec_e2(ie) * val &
                        * w_theta(it)
                endif
             endif

             jx = ix+2
             if (jx <= n_xi) then
                j = mindx(js,je,jx,jt)
                rfac = (1.0/r(ir)) * g_old_global(ir,j)
                val = cderiv(id) * rfac &
                     * driftth_const * driftth(is,it) &
                     * (ix+1.0)*(ix+2.0)/((2*ix+3.0)*(2*ix+5.0))
                g(i) = g(i) + emat_e1(eindx(ie,je)) * val
                if(ix==0 .and. je==ie) then
                   spflux(is) = spflux(is) + evec_e1(ie) * val &
                        * w_theta(it)
                   seflux(is) = seflux(is) + evec_e2(ie) * val &
                        * w_theta(it)
                endif
             endif
          endif
       enddo
    enddo

    ! Er parts of the drift

    src_P0    = (1.0*Z(is))/temp(is,ir) * dphi0dr(ir)
    js = is; jt = it
    
    do je=1,n_energy
       
       jx = ix
       j = mindx(js,je,jx,jt)
       val1 = -driftx(is,it) * src_P0 * g_old_global(ir,j) &
            * (1.0 + (ix+1.0)*(ix+1.0)/((2*ix+1.0)*(2*ix+3.0)) &
            + (ix*ix)/((2*ix+1.0)*(2*ix-1.0))) * drifter1_const
       val2 = -driftx(is,it) * src_P0 * g_old_global(ir,j) &
            * 0.5 * (ix*ix/(2*ix-1.0) &
            - ix*(ix+1.0)*(ix+1.0)/((2*ix+1.0)*(2*ix+3.0)) &
            - ix*ix*ix/((2*ix+1.0)*(2*ix-1.0))) * drifter2_const
       g(i) = g(i)  &
            + val1 * (ematij_e1de(ie,je) - emat_e1(eindx(ie,je))) &
            + val2 * emat_e0(eindx(ie,je)) 
       if(ix==0 .and. je==ie) then
          spflux(is) = spflux(is) &
               + val1 * (evec_e1de(ie) - evec_e1(ie)) * w_theta(it) &
               + val2 * evec_e0(ie) * w_theta(it)
          seflux(is) = seflux(is) &
               + val1 * (evec_e2de(ie) - evec_e2(ie)) * w_theta(it) &
               + val2 * evec_e1(ie) * w_theta(it)
       endif
       

       ! theta drivative part -- only for ix=jx
       do id=-2,2
          if (id /= 0) then
             kt = thcyc(it+id)
             j = mindx(js,je,jx,kt)
             val = cderiv(id) / (12*d_theta) * g_old_global(ir,j) &
                  * I_div_psip * k_par(it) / Bmag(it) * rho(ir) * dphi0dr(ir) &
                  * drifterth_const
             g(i) = g(i) + val * emat_e0(eindx(ie,je)) 
             if(ix==0 .and. je==ie) then
                spflux(is) = spflux(is) + val * evec_e0(ie) * w_theta(it)
                seflux(is) = seflux(is) + val * evec_e1(ie) * w_theta(it)
             endif
          endif
       enddo

       jx = ix-2
       if (jx >= 0) then
          j = mindx(js,je,jx,jt)
          val1 = -driftx(is,it) * src_P0 * g_old_global(ir,j) &
               * (ix)*(ix-1.0)/((2*ix-3.0)*(2*ix-1.0)) * drifter1_const
          val2 = -driftx(is,it) * src_P0 * g_old_global(ir,j) &
               * 0.5 * (-ix*(ix-2.0)*(ix-1.0)/((2*ix-3.0)*(2*ix-1.0))) &
               * drifter2_const
          g(i) = g(i) + val1 * (ematij_e1de(ie,je) - emat_e1(eindx(ie,je))) &
               + val2 * emat_e0(eindx(ie,je))
          if(ix==0 .and. je==ie) then
             spflux(is) = spflux(is) &
                  + val1 * (evec_e1de(ie) - evec_e1(ie)) * w_theta(it) &
                  + val2 * evec_e0(ie) * w_theta(it)
             seflux(is) = seflux(is) &
                  + val1 * (evec_e2de(ie) - evec_e2(ie)) * w_theta(it) &
                  + val2 * evec_e1(ie) * w_theta(it)
          endif
       endif
       
       jx = ix+2
       if (jx <= n_xi) then
          j = mindx(js,je,jx,jt)
          val1 = -driftx(is,it) * src_P0 * g_old_global(ir,j) &
               * (ix+1.0)*(ix+2.0)/((2*ix+3.0)*(2*ix+5.0)) * drifter1_const
          val2 = -driftx(is,it) * src_P0 * g_old_global(ir,j) &
               * 0.5 * ((ix+2.0)*(ix+1.0)/(2*ix+3.0) &
               - (ix+1.0)*(ix+2.0)*(ix+2.0)/((2*ix+3.0)*(2*ix+5.0))) &
               * drifter2_const
          g(i) = g(i) + val1 * (ematij_e1de(ie,je) - emat_e1(eindx(ie,je))) &
               + val2 * emat_e0(eindx(ie,je))
          if(ix==0 .and. je==ie) then
             spflux(is) = spflux(is) &
                  + val1 * (evec_e1de(ie) - evec_e1(ie)) * w_theta(it) &
                  + val2 * evec_e0(ie) * w_theta(it)
             seflux(is) = seflux(is) &
                  + val1 * (evec_e2de(ie) - evec_e2(ie)) * w_theta(it) &
                  + val2 * evec_e1(ie) * w_theta(it)
          endif
       endif
       
    enddo
    
  end subroutine get_drift

end subroutine neo_do
