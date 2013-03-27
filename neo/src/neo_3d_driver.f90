module neo_3d_driver
  implicit none

  public :: ThreeD_do

  integer, parameter,private :: io_neo=10, io_f=11
  character(len=80),private  :: runfile_f = 'out.neo.f_3d'
  integer, dimension(:,:,:,:,:), allocatable, private :: mindx 
  real, dimension(:), allocatable, private :: g0
  integer, dimension(:,:,:,:,:), allocatable, private :: mindx0

contains

  subroutine ThreeD_do
    use neo_globals
    use neo_energy_grid
    use neo_allocate_profile
    use neo_3d_globals
    use neo_3d_equilibrium
    use neo_3d_transport
    implicit none
    integer :: ir, it, ip
    integer :: status

    ! cyclic index (for theta-periodicity)
    allocate(thcyc(1-n_theta:2*n_theta))
    do it=1,n_theta
       thcyc(it-n_theta) = it
       thcyc(it) = it
       thcyc(it+n_theta) = it
    enddo
    allocate(vpcyc(1-n_varphi:2*n_varphi))
    do ip=1,n_varphi
       vpcyc(ip-n_varphi) = ip
       vpcyc(ip) = ip
       vpcyc(ip+n_varphi) = ip
    enddo
    ! coefficients for 4th order centered derivative
    cderiv(-2) =  1
    cderiv(-1) = -8
    cderiv(0)  =  0
    cderiv(1)  =  8
    cderiv(2)  = -1

    ! Set-up the energy grids for basis ints (indep of r)
    call ENERGY_basis_ints_alloc(1)
    call ENERGY_basis_ints
    call ENERGY_coll_ints_alloc(1)

    call ThreeD_EQUIL_alloc(1)
    call ThreeD_TRANSP_alloc(1)

    if (silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_neo,file=trim(path)//'out.neo.grid_3d',status='replace')
       write(io_neo,*) n_species
       write(io_neo,*) n_energy
       write(io_neo,*) n_xi
       write(io_neo,*) n_theta
       do it=1, n_theta
          write(io_neo,*) theta(it)
       enddo
       write(io_neo,*) n_varphi
       do ip=1, n_varphi
          write(io_neo,*) varphi(ip)
       enddo
       write(io_neo,*) n_radial
       do ir=1, n_radial
          write(io_neo,*) r(ir)
       enddo
       close(io_neo)
       
       open(unit=io_f,file=trim(path)//runfile_f,status='replace')
       close(io_f)
    end if

    do ir=1, n_radial
       call threed_matrix_solve(ir,0,status)
       if(status == 1) goto 200
       call threed_matrix_solve(ir,1,status)
       if(status == 1) goto 200
       call ThreeD_TRANSP_write(ir)
    enddo

200 call ENERGY_basis_ints_alloc(0)
    call ENERGY_coll_ints_alloc(0)
    call ThreeD_EQUIL_alloc(0)
    call ThreeD_TRANSP_alloc(0)
    call PROFILE_SIM_alloc(0)
    deallocate(thcyc)
    deallocate(vpcyc)

  end subroutine ThreeD_do


  subroutine threed_matrix_solve(ir,nonaxisym_flag,status)
    use neo_globals
    use neo_energy_grid
    use neo_sparse_solve
    use neo_3d_globals
    use neo_3d_equilibrium
    use neo_3d_transport
    implicit none
    integer, intent(in)  :: ir
    integer, intent(in)  :: nonaxisym_flag ! 1 for 3d, otherwise 2d
    integer, intent(out) :: status
    integer :: np_loc
    integer :: is, ie, ix, it, ip, js, je, jx, jt, jp, ks
    integer :: i, j, k, id
    integer :: ierr
    integer :: n_elem
    real, dimension(:), allocatable :: a
    integer, dimension(:), allocatable :: a_iindx
    integer, dimension(:), allocatable :: a_jindx
    integer :: ifac, matfac_err, max_ifac=5
    real :: stream_t, stream_p, trap

    status = 0

    if(nonaxisym_flag == 1) then
       np_loc = n_varphi
    else
       np_loc = 1
    endif

    ! Matrix solve allocations
    n_row = n_species*(n_energy+1)*(n_xi+1)*n_theta*np_loc
    ! Estimate number of elements
    i = 0
    ! constraint
    if(collision_model == 1 .or. collision_model == 2) then
       i = i + n_species*n_theta*np_loc * (n_energy+1)
    else
       i = i + n_species*n_theta*np_loc * 2
    endif
    ! collisions
    i = i + n_species*(n_xi+1)*n_theta*np_loc*(n_energy+1)**2 &
         * (1 + n_species)
    ! streaming (theta and varphi derivs)
    i = i + n_species*n_theta*np_loc*(n_energy+1)**2*((n_xi-1)*4*2 + 2*4)
    if(nonaxisym_flag == 1) then
       i = i + n_species*n_theta*np_loc*(n_energy+1)**2*((n_xi-1)*4*2 + 2*4)
    endif
    ! trapping
    i = i + n_species*n_theta*np_loc*(n_energy+1)**2*((n_xi-1)*2 + 2)
    n_max = i
    allocate(a(n_max),stat=ierr)
    if(ierr /= 0) then
       call neo_error('ERROR: (NEO) Array allocation failed')
       status=1
       goto 100
    end if
    allocate(a_iindx(n_max),stat=ierr)
    if(ierr /= 0) then
       call neo_error('ERROR: (NEO) Array allocation failed')
       status=1
       goto 100
    end if
    allocate(a_jindx(n_max),stat=ierr)
    if(ierr /= 0) then
       call neo_error('ERROR: (NEO) Array allocation failed')
       status=1
       goto 100
    end if
    allocate(g(n_row))
    allocate(mindx(n_species,0:n_energy,0:n_xi,n_theta,np_loc))
    allocate(is_indx(n_row))
    allocate(ie_indx(n_row))
    allocate(ix_indx(n_row))
    allocate(it_indx(n_row))
    allocate(ip_indx(n_row))

    ! matrix indices
    i = 0
    do is=1,n_species
       do ie=0,n_energy
          do ix=0,n_xi
             do it=1,n_theta
                do ip=1,np_loc
                   i = i+1
                   mindx(is,ie,ix,it,ip) = i
                   is_indx(i) = is
                   ie_indx(i) = ie
                   ix_indx(i) = ix
                   it_indx(i) = it
                   ip_indx(i) = ip
                enddo
             enddo
          enddo
       enddo
    enddo
    
    allocate(driftx_2d(n_species,n_theta))
    allocate(driftx_3d(n_species,n_theta,n_varphi))
    
    ! Set-up the matrix equation: LHS radially local matrix
     
    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_neoout,file=trim(path)//runfile_neoout,status='old',position='append')
       write(io_neoout,*) 'ir = ', ir
       write(io_neoout,*) 'Begin matrix set-up'
       close(io_neoout)
    endif
    
    if(collision_model == 3) then
       do is=2, n_species
          if(abs(temp(is,ir)-temp(1,ir)) > 1e-3) then
             if(silent_flag == 0 .and. i_proc == 0) then
                open(unit=io_neoout,file=trim(path)//runfile_neoout,status='old',position='append')
                write(io_neoout,*) 'WARNING: (NEO) Full HS collisions with unequal temps'
                close(io_neoout)
             endif
          endif
       enddo
    endif
    
    ! Get the equilibrium parameters (th)
    call ThreeD_EQUIL_do(ir)
    
    ! Get the energy grids for coll ints
    call ENERGY_coll_ints(ir)
    
    ! Set the LHS
    a_iindx(:) = 0
    a_jindx(:) = 0
    a(:) = 0.0
    k = 0
    
    do is=1,n_species
       do it=1,n_theta
          ! vdrift dot grad r 
          ! -> driftx * (1+xi^2) * d/d(r/a) * ene(ie)
          driftx_2d(is,it) = v_drift_x_overB2_0(it) &
               * mass(is)/(1.0*Z(is)) * (vth(is,ir))**2
          do ip=1,n_varphi
             driftx_3d(is,it,ip) = v_drift_x_overB2_1(it,ip) &
                  * mass(is)/(1.0*Z(is)) * (vth(is,ir))**2
          enddo
       enddo
    enddo
    
    do is=1,n_species
       do ie=0,n_energy
          do ix=0,n_xi
             do it=1,n_theta
                do ip=1,np_loc
                   
                   ! vpar bhat dot grad 
                   ! -> stream * xi * sqrt(ene)
                   stream_t = sqrt(2.0) * vth(is,ir) * k_par_t_0(it) &
                        / (12*d_theta)
                   stream_p = sqrt(2.0) * vth(is,ir) * k_par_p_0(it) &
                        / (12*d_varphi)
                   ! mu bdot grad B d/dvpar 
                   ! -> trap * (1-xi^2)*d/dxi * sqrt(ene)
                   trap   = gradpar_Bmag_overB_0(it) &
                        * sqrt(0.5) * vth(is,ir)
                   
                   i = mindx(is,ie,ix,it,ip)
                   
                   ! Impose constant for g0
                   ! local  -> at each species and ene 
                   
                   if (ix==0 .and. it == 1 .and. ip == 1 &
                        .and. (collision_model==1 .or. collision_model==2))&
                        then
                      ! <f_ie> = 0
                      js = is; je = ie; jx=0
                      do jt=1, n_theta
                         do jp=1, np_loc
                            j = mindx(js,je,jx,jt,jp)
                            k = k+1
                            a(k) = w_theta_0(jt)/sum_w_theta_0
                            a_iindx(k) = i
                            a_jindx(k) = j
                         enddo
                      enddo
                      
                   else if(ix==0 .and. it == 1 .and. ip == 1 &
                        .and. (ie==0 .or. ie==1) &
                        .and. (collision_model==3 .or. &
                        collision_model==4 .or. collision_model==5)) then
                      ! <f_ie> = 0
                      js = is; je = ie; jx=0
                      do jt=1, n_theta
                         do jp=1, np_loc
                            j = mindx(js,je,jx,jt,jp)
                            k = k+1
                            a(k) = w_theta_0(jt)/sum_w_theta_0
                            a_iindx(k) = i
                            a_jindx(k) = j
                         enddo
                      enddo
                      
                   else
                      
                      ! Collisions 
                      jt = it; jp = ip; jx = ix
                      
                      ! test particle
                      js = is
                      do je=0, n_energy
                         j = mindx(js,je,jx,jt,jp)
                         k = k+1
                         a(k) = 0.0
                         do ks=1, n_species
                            a(k) = a(k) &
                                 - emat_coll_test(is,ks,ie,je,ix)
                         enddo
                         a_iindx(k) = i
                         a_jindx(k) = j
                      enddo
                      
                      ! field particle                     
                      do je=0, n_energy
                         do js=1,n_species
                            j = mindx(js,je,jx,jt,jp)
                            k = k+1
                            a(k) = -emat_coll_field(is,js,ie,je,ix)
                            a_iindx(k) = i
                            a_jindx(k) = j
                         enddo
                      enddo
                      
                      
                      ! Streaming -- d/dtheta
                      js = is; jp = ip
                      do je=0, n_energy
                         do id=-2,2
                            if (id /= 0) then
                               jt = thcyc(it+id)
                               jx = ix-1
                               if (jx >= 0) then
                                  j = mindx(js,je,jx,jt,jp)
                                  k = k+1
                                  a(k) = stream_t &
                                       * ix/(2*ix-1.0) * cderiv(id) &
                                       * emat_e05(ie,je,ix,1) 
                                  a_iindx(k) = i
                                  a_jindx(k) = j
                               endif
                               jx = ix+1
                               if (jx <= n_xi) then
                                  j = mindx(js,je,jx,jt,jp)
                                  k = k+1
                                  a(k) = stream_t &
                                       * (ix+1.0)/(2*ix+3.0) * cderiv(id) &
                                       * emat_e05(ie,je,ix,2) 
                                  a_iindx(k) = i
                                  a_jindx(k) = j
                               endif
                            endif
                         enddo
                      enddo
                      
                      ! Streaming -- d/dvarphi
                      if(nonaxisym_flag == 1) then
                         js = is; jt = it
                         do je=0, n_energy
                            do id=-2,2
                               if (id /= 0) then
                                  jp = vpcyc(ip+id)
                                  jx = ix-1
                                  if (jx >= 0) then
                                     j = mindx(js,je,jx,jt,jp)
                                     k = k+1
                                     a(k) = stream_p &
                                          * ix/(2*ix-1.0) * cderiv(id) &
                                          * emat_e05(ie,je,ix,1) 
                                     a_iindx(k) = i
                                     a_jindx(k) = j
                                  endif
                                  jx = ix+1
                                  if (jx <= n_xi) then
                                     j = mindx(js,je,jx,jt,jp)
                                     k = k+1
                                     a(k) = stream_p &
                                          * (ix+1.0)/(2*ix+3.0) &
                                          * cderiv(id) &
                                          * emat_e05(ie,je,ix,2) 
                                     a_iindx(k) = i
                                     a_jindx(k) = j
                                  endif
                               endif
                            enddo
                         enddo
                      endif
                      
                      ! Trapping 
                      js = is; jt = it; jp = ip
                      do je=0, n_energy
                         jx = ix-1
                         if (jx >= 0) then
                            j = mindx(js,je,jx,jt,jp)
                            k = k+1
                            a(k) = trap &
                                 * ix*(ix-1.0)/(2*ix-1.0) &
                                 * emat_e05(ie,je,ix,1) 
                            a_iindx(k) = i
                            a_jindx(k) = j                      
                         endif
                         jx = ix+1
                         if (jx <= n_xi) then
                            j = mindx(js,je,jx,jt,jp)
                            k = k+1
                            a(k) = -trap &
                                 * (ix+1.0)*(ix+2.0)/(2*ix+3.0) &
                                 * emat_e05(ie,je,ix,2)
                            a_iindx(k) = i
                            a_jindx(k) = j
                         endif
                      enddo
                      
                   endif ! bc if/else
                enddo ! ip
             enddo ! it
          enddo ! ix
       enddo ! ie
    enddo ! is
    
    ! Factor the Matrix -- uses a(:) and a_indx(:)
    
    n_elem = k
    n_max = n_elem*matsz_scalefac
    matfac_err = 0
    
    do ifac = 1, max_ifac
       
       if(allocated(amat))       deallocate(amat)
       allocate(amat(n_max),stat=ierr)
       if(ierr /= 0) then
          call neo_error('ERROR: (NEO) Array allocation failed')
          status=1
          goto 100
       end if
       if(allocated(amat_indx))  deallocate(amat_indx)
       allocate(amat_indx(2*n_max),stat=ierr)
       if(ierr /= 0) then
          call neo_error('ERROR: (NEO) Array allocation failed')
          status=1
          goto 100
       end if
       
       amat(:) = 0.0
       amat_indx(:) = 0
       do k=1,n_elem
          amat(k) = a(k)
          amat_indx(k) = a_iindx(k)
          amat_indx(n_elem+k) = a_jindx(k)
       enddo
       
       if(silent_flag == 0 .and. i_proc == 0) then
          open(unit=io_neoout,file=trim(path)//runfile_neoout,&
               status='old',position='append')
          if(ifac == 1) then
             write(io_neoout,*) 'Begin matrix factor'
          else
             write(io_neoout,*) 'Re-trying matrix factorization'
          endif
          close(io_neoout)
       endif
       call SOLVE_factor(n_elem)
       if(error_status > 0) then
          error_status = 0
          n_max = n_max * 2
       else
          matfac_err = 1
          exit
       endif
    enddo
    if(matfac_err == 0) then
       call neo_error('ERROR: (NEO) Matrix factorization failed. Try increasing MATSZ_SCALEFAC.')
       status=1
       goto 100
    endif
    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_neoout,file=trim(path)//runfile_neoout,&
            status='old',position='append')
       write(io_neoout,*) 'Done matrix factor'
       close(io_neoout)
    endif
    
    ! Set the RHS source 
    if(nonaxisym_flag == 1) then
       call set_threed_source(ir)
    else
       call set_twod_source(ir)
    endif
    
    ! Matrix solve -- uses g(:), a(:), and a_indx(:)
    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_neoout,file=trim(path)//runfile_neoout,&
            status='old',position='append')
       write(io_neoout,*) 'Begin matrix solve'
       close(io_neoout)
    endif
    call SOLVE_do
    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_neoout,file=trim(path)//runfile_neoout,&
            status='old',position='append')
       write(io_neoout,*) 'Done matrix solve'
       close(io_neoout)
    endif
    
    ! Compute the neo transport coefficients
    call ThreeD_TRANSP_do(ir,nonaxisym_flag)
    
    if(error_status > 0) then
       status=1
       goto 100
    endif
    
    if(silent_flag == 0 .and. i_proc == 0) then
       open(io_f,file=trim(path)//runfile_f,status='old',position='append')
       write(io_f,*) g(:)
       close(io_f)
    endif
    
    if(nonaxisym_flag /= 1) then
       allocate(g0(n_row))
       allocate(mindx0(n_species,0:n_energy,0:n_xi,n_theta,np_loc))
       g0(:) = g(:)
       mindx0(:,:,:,:,:) = mindx(:,:,:,:,:)
    endif
    
100 continue
    if(allocated(a))          deallocate(a)
    if(allocated(a_iindx))    deallocate(a_iindx)
    if(allocated(a_jindx))    deallocate(a_jindx)
    if(allocated(amat))       deallocate(amat)
    if(allocated(amat_indx))  deallocate(amat_indx)
    if(allocated(g))          deallocate(g)
    if(allocated(mindx))      deallocate(mindx)
    if(allocated(is_indx))    deallocate(is_indx)
    if(allocated(ie_indx))    deallocate(ie_indx)
    if(allocated(ix_indx))    deallocate(ix_indx)
    if(allocated(it_indx))    deallocate(it_indx)
    if(allocated(ip_indx))    deallocate(ip_indx)
    if(allocated(driftx_2d))  deallocate(driftx_2d)
    if(allocated(driftx_3d))  deallocate(driftx_3d)
    if(nonaxisym_flag == 1) then
       if(allocated(g0))          deallocate(g0)
       if(allocated(mindx0))      deallocate(mindx0)
    endif

  end subroutine threed_matrix_solve

  subroutine set_twod_source(ir)
    use neo_globals
    use neo_energy_grid
    use neo_3d_globals
    implicit none
    integer, intent(in) :: ir
    real :: src_F0_Ln, src_F0_Lt, src_P0
    integer :: i, is, ie, ix, it, ip

    ! First-Order Source Term
    g(:) = 0.0

    ip = 1
    do is=1,n_species
       
       ! src = (1/F0) * dF0/dr + Ze/T dPhi0/dr
       src_F0_Ln = -(dlnndr(is,ir) - 1.5*dlntdr(is,ir))  ! ene^0 part
       src_F0_Lt = -dlntdr(is,ir)                        ! ene^1 part
       src_P0    = (1.0*Z(is))/temp(is,ir) * dphi0dr(ir) ! ene^0 Er part
       
       do ie=0,n_energy
          do ix=0,n_xi
             do it=1,n_theta
                   
                i = mindx(is,ie,ix,it,ip)
                
                ! Impose constant for g0
                
                if(ix==0 .and. it==1 .and. ip==1 .and. &
                     (collision_model==1 .or. collision_model==2)) then
                   g(i) = 0.0
                   
                else if(ix==0 .and. it==1 .and. ip==1 &
                     .and. (ie==0 .or. ie==1) &
                     .and. (collision_model==3 .or. &
                     collision_model==4 .or. collision_model==5)) then
                   g(i) = 0.0
                   
                else 

                   ! Equilibrium source:
                   ! -(vdrift dot grad (F0 + Ze/T Phi0))
                   
                   if (ix == 0) then
                      g(i) = -(4.0/3.0) * driftx_2d(is,it) &
                           * ( (src_F0_Ln + src_P0) &
                           * evec_e1(ie,ix) &
                           + src_F0_Lt * evec_e2(ie,ix) )
                      
                   else if (ix == 2) then
                      g(i) = -(2.0/3.0) * driftx_2d(is,it) &
                           * ( (src_F0_Ln + src_P0) &
                           * evec_e1(ie,ix) &
                           + src_F0_Lt * evec_e2(ie,ix) )
                      
                   endif
                   
                endif
                
             enddo
          enddo
       enddo
    enddo

  end subroutine set_twod_source
    

  subroutine set_threed_source(ir)
    use neo_globals
    use neo_energy_grid
    use neo_3d_globals
    use neo_3d_equilibrium
    implicit none
    integer, intent(in) :: ir
    real :: src_F0_Ln, src_F0_Lt, src_P0
    integer :: i, is, ie, ix, it, ip, j, js, je, jx, jt, jp, ks, id
    real :: stream_t, trap

    ! First-Order Source Term = -L1 g0 + S1 f0

    !!!!!!!!!!!!!!!!!!!!!!!!!
    ! -L1 * g0
    !!!!!!!!!!!!!!!!!!!!!!!!!

    g(:) = 0.0
    do is=1,n_species
       do ie=0,n_energy
          do ix=0,n_xi
             do it=1,n_theta
                do ip=1,n_varphi
                
                   i = mindx(is,ie,ix,it,ip)

                   ! vpar bhat dot grad 
                   ! -> stream * xi * sqrt(ene)
                   stream_t = sqrt(2.0) * vth(is,ir) * k_par_t_1(it,ip) &
                        / (12*d_theta)
                   ! mu bdot grad B d/dvpar 
                   ! -> trap * (1-xi^2)*d/dxi * sqrt(ene)
                   trap   = gradpar_Bmag_overB_1(it,ip) &
                        * sqrt(0.5) * vth(is,ir)
       
                   ! Streaming -- d/dtheta
                   js = is; jp = 1
                   do je=0, n_energy
                      do id=-2,2
                         if (id /= 0) then
                            jt = thcyc(it+id)
                            jx = ix-1
                            if (jx >= 0) then
                               j = mindx0(js,je,jx,jt,jp)
                               g(i) = g(i) - stream_t &
                                    * ix/(2*ix-1.0) * cderiv(id) &
                                    * emat_e05(ie,je,ix,1) * g0(j)
                            endif
                            jx = ix+1
                            if (jx <= n_xi) then
                               j = mindx0(js,je,jx,jt,jp)
                               g(i) = g(i) - stream_t &
                                    * (ix+1.0)/(2*ix+3.0) * cderiv(id) &
                                    * emat_e05(ie,je,ix,2) * g0(j) 
                            endif
                         endif
                      enddo
                   enddo
                   
                                            
                   ! Trapping 
                   js = is; jt = it; jp = 1
                   do je=0, n_energy
                      jx = ix-1
                      if (jx >= 0) then
                         j = mindx0(js,je,jx,jt,jp)
                         g(i) = g(i) - trap &
                              * ix*(ix-1.0)/(2*ix-1.0) &
                              * emat_e05(ie,je,ix,1) * g0(j)            
                      endif
                      jx = ix+1
                      if (jx <= n_xi) then
                         j = mindx0(js,je,jx,jt,jp)
                         g(i) = g(i)  + trap &
                              * (ix+1.0)*(ix+2.0)/(2*ix+3.0) &
                              * emat_e05(ie,je,ix,2) * g0(j)
                      endif
                   enddo
                   
                enddo
             enddo
          enddo
       enddo
    enddo
                   
    !!!!!!!!!!!!!!!!!!!!!!!!!
    ! S1 * f0 and constraint
    !!!!!!!!!!!!!!!!!!!!!!!!!

    do is=1,n_species
       
       ! src = (1/F0) * dF0/dr + Ze/T dPhi0/dr
       src_F0_Ln = -(dlnndr(is,ir) - 1.5*dlntdr(is,ir))  ! ene^0 part
       src_F0_Lt = -dlntdr(is,ir)                        ! ene^1 part
       src_P0    = (1.0*Z(is))/temp(is,ir) * dphi0dr(ir) ! ene^0 Er part
       
       do ie=0,n_energy
          do ix=0,n_xi
             do it=1,n_theta
                do ip=1, n_varphi
                   
                   i = mindx(is,ie,ix,it,ip)
                   
                   ! Impose constant for g0
                   
                   if(ix==0 .and. it==1 .and. ip==1 .and. &
                        (collision_model==1 .or. collision_model==2)) then
                      g(i) = 0.0
                      do jt=1,n_theta
                         do jp=1,n_varphi
                            g(i) = g(i) - w_theta_1(jt,jp)/sum_w_theta_0
                         enddo
                      enddo
                      
                   else if(ix==0 .and. it==1 .and. ip==1 &
                        .and. (ie==0 .or. ie==1) &
                        .and. (collision_model==3 .or. &
                        collision_model==4 .or. collision_model==5)) then
                      g(i) = 0.0
                      do jt=1,n_theta
                         do jp=1,n_varphi
                            g(i) = g(i) - w_theta_1(jt,jp)/sum_w_theta_0
                         enddo
                      enddo
                      
                   else 
                      
                      ! Equilibrium source:
                      ! -(vdrift dot grad (F0 + Ze/T Phi0))
                      
                      if (ix == 0) then
                         g(i) = g(i) - (4.0/3.0) * driftx_3d(is,it,ip) &
                              * ( (src_F0_Ln + src_P0) &
                              * evec_e1(ie,ix) &
                              + src_F0_Lt * evec_e2(ie,ix) )
                         
                      else if (ix == 2) then
                         g(i) = g(i) - (2.0/3.0) * driftx_3d(is,it,ip) &
                              * ( (src_F0_Ln + src_P0) &
                              * evec_e1(ie,ix) &
                              + src_F0_Lt * evec_e2(ie,ix) )
                         
                      endif
                      
                   endif
                   
                enddo
             enddo
          enddo
       enddo
    enddo
    
  end subroutine set_threed_source

end module neo_3d_driver
