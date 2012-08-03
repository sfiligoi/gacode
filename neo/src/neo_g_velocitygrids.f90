module neo_g_velocitygrids
  
  implicit none
  
  public :: g_energy, g_xi

contains  

  subroutine g_energy(ir)
    use neo_globals
    use neo_equilibrium, only : w_theta
    use neo_energy_grid
    implicit none
    integer, intent(in) :: ir
    real, dimension(:), allocatable :: ene
    real, dimension(:,:,:), allocatable :: xval
    real :: xvall
    real,parameter :: emin=0.0, emax=10.0
    integer,parameter :: ne=100
    integer,parameter :: io=52
    integer :: i, is, ie, ix, it, je
    real, dimension(:,:,:,:), allocatable :: gall

    allocate(ene(ne))
    do je=1, ne
       ene(je) = emin + (je-1) * (emax-emin) / (ne-1)
    enddo

    allocate(xval(0:n_energy,ne,0:n_xi))
    do ix=0, n_xi
       do ie=0,n_energy
          do je=1,ne
             call compute_laguerre(ie,e_lag(ix)*0.5,sqrt(ene(je))**e_alpha,xvall)
             xval(ie,je,ix) = xvall * sqrt(ene(je))**xi_beta_l(ix)
          enddo
       enddo
    enddo

    allocate(gall(n_species,0:n_xi,n_theta,ne))
    gall(:,:,:,:) = 0.0

    do je=1, ne
       do i=1,n_row      
          is = is_indx(i)
          ie = ie_indx(i)
          ix = ix_indx(i)
          it = it_indx(i)
          gall(is,ix,it,je) =  gall(is,ix,it,je) + xval(ie,je,ix) * g(i)
       enddo
    enddo

    if(silent_flag == 0 .and. i_proc == 0) then
       if(ir == 1) then
          open(unit=io,file=trim(path)//'out.neo.g_ene_x',status='replace')
          do je=1,ne
             write (io,'(e16.8)',advance='no') ene(je)
          enddo
          close(io)
          open(unit=io,file=trim(path)//'out.neo.g_ene_y',status='replace')
       else
          open(unit=io,file=trim(path)//'out.neo.g_ene_y',status='old',position='append')
       endif
       do is=1, n_species
          do ix=0, n_xi
             do it=1,n_theta
                do je=1, ne
                   write (io,'(e16.8)',advance='no') gall(is,ix,it,je)
                enddo
             enddo
          enddo
       enddo
       close(io)
    endif

    deallocate(ene)
    deallocate(xval)
    deallocate(gall)

  end subroutine g_energy

  
  subroutine g_xi
    use neo_globals
    implicit none
    real, dimension(:,:), allocatable :: xval
    real :: xval1
    integer :: i, is, ie, ix, it, jx
    real, dimension(:), allocatable :: xi
    real, dimension(:,:,:), allocatable :: g0, g1, gall
    integer, parameter :: nxi=100
    integer, parameter :: io=51
    
    allocate(xi(nxi))
    do jx=1, nxi
       xi(jx) = -1.0 + (jx-1)*2.0/(nxi-1)
    enddo
    
    allocate(xval(0:n_xi,nxi))

    do ix=0,n_xi
       do jx=1,nxi
          call compute_legendre(ix,xi(jx),xval1)
          xval(ix,jx) = xval1
       enddo
    enddo

    allocate(gall(n_species,n_theta,nxi))
    gall(:,:,:) = 0.0
    allocate(g0(n_species,n_theta,nxi))
    g0(:,:,:) = 0.0
    allocate(g1(n_species,n_theta,nxi))
    g1(:,:,:) = 0.0

    do jx=1,nxi
       do i=1,n_row      
          is = is_indx(i)
          ie = ie_indx(i)
          ix = ix_indx(i)
          it = it_indx(i)
          if(ie == 1) then
             gall(is,it,jx) = gall(is,it,jx) + xval(ix,jx) * g(i)
             if(mod(ix,2) == 0) then
                g0(is,it,jx) = g0(is,it,jx) + xval(ix,jx) * g(i)
             else
                g1(is,it,jx) = g1(is,it,jx) + xval(ix,jx) * g(i)
             endif
          endif
       enddo
    enddo

    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=io,file=trim(path)//'out.neo.g_xi',status='replace')
       do is=1, n_species
          do it=1, n_theta
             do jx=1, nxi
                write (io,'(e16.8)',advance='no') xi(jx)
                write (io,'(e16.8)',advance='no') gall(is,it,jx)
                write (io,'(e16.8)',advance='no') g0(is,it,jx)
                write (io,'(e16.8)',advance='no') g1(is,it,jx)
                write (io,*)
             enddo
          enddo
       enddo
       close(io)
    end if
    
    deallocate(xi)
    deallocate(xval)
    deallocate(gall)
    deallocate(g0)
    deallocate(g1)
    
  end subroutine g_xi

  subroutine compute_laguerre(n,k,arg,val)
    integer, intent (in) :: n
    real, intent (in) :: k, arg
    real, intent(out) :: val
    integer :: i
    real :: L0, L1

    val = 1.0
    L1  = 0.0
    
    do i=1,n
       L0  = L1
       L1  = val
       val = ((2*i-1+k-arg) * L1 - (i-1+k)*L0) / (1.0*i)
    enddo

  end subroutine compute_laguerre

  subroutine compute_legendre(n,arg,val)
    integer, intent (in) :: n
    real, intent (in) :: arg
    real, intent(out) :: val
    real :: pmm, pmmp1, pnn
    integer :: k
    
    pmm=1.0

    if(n==0) then
       val = pmm
    else
       pmmp1 = arg*pmm;
       if(n==1) then
          val = pmmp1
       else
          do k=2, n
             pnn = (arg*(2*k-1)*pmmp1 - (k-1)*pmm)/(1.0*k)
             pmm=pmmp1
             pmmp1=pnn
          enddo
          val = pnn
       end if
    end if

  end subroutine compute_legendre
  
end module neo_g_velocitygrids
