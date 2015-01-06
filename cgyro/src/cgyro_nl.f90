!--------------------------------------------------------------
! Evaluate nonlinear bracket with direct dealiased convolution
!--------------------------------------------------------------

subroutine cgyro_nl_direct(ij)

  use timer_lib
  use parallel_lib

  use cgyro_globals
  use cgyro_equilibrium

  implicit none

  integer, intent(in) :: ij
  integer :: nx,ny,nx0,ny0
  integer :: ix,ixp,iy,iyp
  integer :: ixpp,iypp
  integer :: ir,it,j,in
  complex, dimension(:,:), allocatable :: f
  complex, dimension(:,:), allocatable :: g
  complex, dimension(:,:), allocatable :: fg


  ny0 = n_toroidal-1
  nx0 = n_radial/2
  ny = int(1.5*ny0)+1
  nx = int(1.5*nx0)+1

  if (.not.allocated(pcyc)) then
     allocate(pcyc(-3*nx:3*nx-1))
     allocate(ncyc(-3*ny:3*ny))
     do ix=-nx,nx-1
        pcyc(ix-2*nx) = ix
        pcyc(ix) = ix
        pcyc(ix+2*nx) = ix
     enddo
     do iy=-ny,ny
        ncyc(iy-2*ny) = iy
        ncyc(iy) = iy
        ncyc(iy+2*ny) = iy
     enddo
  endif

  allocate( f(-nx:nx,-ny:ny) )
  allocate( g(-nx:nx,-ny:ny) )
  allocate(fg(-nx:nx,-ny:ny) )

  call timer_lib_in('comm_nl')
  call parallel_slib_f(h_x,f_nl)
  call parallel_slib_f(psi,g_nl)
  call timer_lib_out('comm_nl')

  call timer_lib_in('rhs_nl')
  do j=1,nsplit
     do it=1,n_theta 

        f = 0.0
        g = 0.0

        ! Array mapping
        do in=1,n_toroidal
           iy = in-1
           do ir=1,n_radial
              ic = ic_c(ir,it) 
              ix = ir-1-nx0
              f(ix,iy) = f_nl(ic,j,in)
              g(ix,iy) = g_nl(ic,j,in)
              f(-ix,-iy) = conjg(f(ix,iy))
              g(-ix,-iy) = conjg(g(ix,iy)) 
           enddo
        enddo

        ! Reality
        do ix=-nx0,nx0
           do iy=1,ny0
              f(-ix,-iy) = conjg(f(ix,iy))
              g(-ix,-iy) = conjg(g(ix,iy))
           enddo
        enddo

        fg = (0.0,0.0)
        do ix=-nx0,nx0-1
           do ixp=-nx,nx-1
              ixpp = pcyc(ix-ixp)
              do iy=0,ny0
                 do iyp=-ny+iy,ny
                    iypp = ncyc(iy-iyp)
                    fg(ix,iy) = fg(ix,iy)+f(ixpp,iypp)*g(ixp,iyp)*(iypp*ixp-iyp*ixpp)
                 enddo
              enddo
           enddo
        enddo

        do ir=1,n_radial
           ic = ic_c(ir,it) 
           ix = ir-1-nx0
           do in=1,n_toroidal
              iy = in-1
              g_nl(ic,j,in) = fg(ix,iy)
           enddo
        enddo

     enddo ! it
  enddo ! j
  call timer_lib_out('rhs_nl')

  call timer_lib_in('comm_nl')
  call parallel_slib_r(g_nl,psi)
  call timer_lib_out('comm_nl')

  deallocate( f)
  deallocate( g)
  deallocate(fg)

  ! RHS -> -[f,g] = (n'' p' - n' p'') f'' g' 

  rhs(ij,:,:) = rhs(ij,:,:)+(q*rho/rmin)*(2*pi/length)*psi(:,:)

end subroutine cgyro_nl_direct

!--------------------------------------------------------------
! Evaluate nonlinear bracket with dealiased FFT 
!--------------------------------------------------------------

subroutine cgyro_nl_fftw(ij)

  use timer_lib
  use parallel_lib

  use cgyro_globals

  ! Wrapper for FFTW++
  use fftwpp 


  implicit none

  integer, intent(in) :: ij

end subroutine cgyro_nl_fftw
