      subroutine 2d_lin_interp(tin,rhon,nj,tarray_in, rarray_in,
     .                          array_in,nx,nt,array_out)
c----------------------------------------------------------------------
c        given array_in(nx,nt) nx = space,nt= time,get array_out(nx) at
c        time tin using linear interpolation in space and time.
c INPUT:
c tin            time at which vales are to be found
c rhon(nj)       normalized [0,1] rho grid.Results will be obtained for
c                values of rho in rhon.
c nj             size of rhon
c tarray_in      input list of times at which second index of array_in is given
c                tarray_in must span the time over which the transport
c                simulation is done.Out of bounds interpolation is fatal.
c rarray_in      input list of normalized rho values at which first index of
c                array_in is given. Must cover [0,1]
c array_in(i,j)  i = 1 is assumed to corespond to rhon(1)  = magnetic axis
c                i = nx                           rhon(nj) = plasma edge
c                intermediate values, 1 < i < nj are NOT assumed
c                to match the values in rhon(i) since typically nx < nj
c nx,nt          actual space and time sizes of arrays
c     
c OUTPUT:
c array_out(nj) spatial values at time tin            
c------------------------------------------------------------HSJ--------
      implicit none
      integer nx,nt,nj,j,i1,i2,k,k1,k2
      real *8 tin,tarray_in(nt),array_out(nx),rarray_in(nx),
    .           rhon(nj),a,array_in(nx,nt),a1,a2,b,rhow


      !first get the time interval that brackets tin:
      do j=1,nt-1
         dt = (tin-tarray_in(j))*(tin-tarray_in(j+1))
         if(dt .le. 0.0)then
            i1 = j
            i2 = j+1
            go to 10
         endif
      enddo
      call STOP('Sub transp_2d_interp time out of bounds',0)


      !now get the space interval and do interpolation
 10   do j=2,nj-1
         rhow = rho(j)
         do k = 1,nx-1
           drho  = (rhow-rarray_in(k))*(rhow-rarray_in(k+1))
           if(drho .le. 0.0)then
             k1 = k
             k2 = k+1
           endif
           !first TIME INTERPOLATION:
           a = (array_in(k1,i2)-array_in(k1,i1))/(tarray_in(i2)
    .                                          -tarray_in(i1))
           b = array_in(k1,i1)-a*tarray_in(i1)
           a1 = a*tin+b  ! value of array_in at rho = rarray_in(k1),time =tin
           a = (array_in(k2,i2)-array_in(k2,i1))/(tarray_in(i2)
    .                                          -tarray_in(i1))
           b = array_in(k2,i1)-a*tarray_in(i1)
           a2 = a*tin+b  ! value of array_in at rho = rarray_in(k2) ,time =tin
           !next  SPACE INTERPOLATION:
           a = (a2-a1)/(rarray_in(k2) -rarray_in(k1))
           b = a1 -a*rarray_in(1)
           array_out(j) = a*rhow + b  !interpolated value at rhow and tin
           go to 20
         enddo
         call STOP('Sub transp_2d_interp rho  out of bounds',0)
 20      continue
      enddo

      !value at the magnetic axis:
      a = (array_in(1,i2)-array_in(1,i1))/(tarray_in(i2)-tarray_in(i1))
      b = array_in(1,i1) - a*tarray_in(i1)
      array_out(1) = a*tin + b

      !value at the plasma edge
      a = (array_in(nx,i2)-array_in(nx,i1))/(tarray_in(i2)-tarray_in(i1))
      b = array_in(nx,i1) - a*tarray_in(i1)
      array_out(nj) = a*tin + b

      return
      end
      
