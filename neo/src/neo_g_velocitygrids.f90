module neo_g_velocitygrids
  
  implicit none
  
  public :: g_energy, g_xi

contains  

  subroutine g_energy(i_iter,ir)
    use neo_globals
    use neo_equilibrium, only : w_theta
    implicit none
    integer, intent(in) :: i_iter, ir
    real, dimension(:), allocatable :: ene, zene
    real :: energy_min, b_arg, xval
    integer :: i, is, ie, ix, it, je
    real, dimension(:,:,:,:), allocatable :: gall
    integer,parameter :: ne=100
    integer,parameter :: io=52

    allocate(ene(ne))
    allocate(zene(ne))
    if(collision_model == 1) then
       energy_min = energy_min_connor
    else
       energy_min = 0.0
    endif
    do je=1, ne
       ene(je) = energy_min + (je-1) * (energy_max-energy_min) / (ne-1)
       zene(je) = 2.0/(1.0-sqrt(energy_min/energy_max)) &
            * sqrt(ene(je)/energy_max) &
            - (1.0+sqrt(energy_min/energy_max)) &
            / (1.0-sqrt(energy_min/energy_max))
    enddo
    
    allocate(gall(n_species,0:n_xi,n_theta,ne))
    gall(:,:,:,:) = 0.0

    do je=1, ne
       if(zene(je) <= -1.0) then
          b_arg = pi  ! use x = -1
       else if(zene(je) >= 1.0) then
          b_arg = 0.0 ! use x = 1
       else
          b_arg = acos(zene(je))
       endif
       do i=1,n_row      
          is = is_indx(i)
          ie = ie_indx(i)
          ix = ix_indx(i)
          it = it_indx(i)
          xval = cos((ie-1) * b_arg)
          gall(is,ix,it,je) =  gall(is,ix,it,je) + xval * g(i)
       enddo
    enddo

    if(write_out_mode > 0) then
       if(i_iter == 1 .and. ir == 1) then
          open(unit=io,file='g_ene_x.out',status='replace')
          do je=1,ne
             write (io,'(e16.8,$)') ene(je)
          enddo
          close(io)
          open(unit=io,file='g_ene_y.out',status='replace')
       else
          open(unit=io,file='g_ene_y.out',status='old',position='append')
       endif
       do is=1, n_species
          do ix=0, n_xi
             do it=1,n_theta
                do je=1, ne
                   write (io,'(e16.8,$)') gall(is,ix,it,je)
                enddo
             enddo
          enddo
       enddo
       close(io)
    endif

    deallocate(ene)
    deallocate(zene)
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

    if(write_out_mode > 0) then
       open(unit=io,file='g_xi.out',status='replace')
       do is=1, n_species
          do it=1, n_theta
             do jx=1, nxi
                write (io,'(e16.8,$)') xi(jx)
                write (io,'(e16.8,$)') gall(is,it,jx)
                write (io,'(e16.8,$)') g0(is,it,jx)
                write (io,'(e16.8,$)') g1(is,it,jx)
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
