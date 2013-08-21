c---------------------
      subroutine eqi_binvec(id_axis,ivec,zval,indx,ierr)
c
c  use eqi_binget (below) to get the bin index -- **vector**
c
      implicit NONE
c
c  input:
      integer id_axis                   ! axis id code number
c                       to suppress warning message:  give -(axis id code)
      integer ivec                      ! vector size
      real*8 zval(ivec)                 ! vector of values to look up
c
c  output:
      integer indx(ivec)                ! indices of zones containing zval's
      integer ierr                      ! completion code, 0=OK
c
c$r8real_input: zval
c---------------------
c
      real*8 zdum1(ivec),zdum2(ivec),zdum3(ivec)
c
c---------------------
      call eqi_binget(id_axis,ivec,zval,indx,zdum1,zdum2,zdum3,ierr)
      return
      end
c---------------------
      subroutine eqi_binget(id_axis,ivec,zval,indx,zparam,zh,zhi,ierr)
c
c  **vectorized**
c  find the bin (zone) in which `zval' belongs, in axis `id_axis'
c  return the index, relative displacement, and zone width info.
c  set error flag if zval is out of range.
c
      use eq_xyz_cmn
c
      implicit NONE
c
c  input:
      integer id_axis                   ! axis id code number
c                       to suppress warning message:  give -(axis id code)
      integer ivec                      ! vector dimensioning
      real*8 zval(ivec)                 ! values to look up
c
c  output:
      integer indx(ivec)                ! index of zone containing zval
c
c  if zval is less that the first val, indx=1, ierr=1 will be returned
c  if zval is greater thatn the last val, indx=N-1, ierr=1 will be returned
c  zval out of range error message will be written
c
      real*8 zparam(ivec)               ! displacement parameter w/in zone
c
c  zparam=0:  zval at low edge of zone
c  zparam=1:  zval at high edge of zone
c
      real*8 zh(ivec)                   ! width of zone, h
      real*8 zhi(ivec)                  ! 1/h
c
      integer ierr                      ! completion code, 0=OK
c
c$r8real_input: zval
c$r8real_output: zparam zh zhi
c-------------------------------
      integer iaddr,isize,ilast,i,id_axis_use
      real*8 zndx,zdeli,zzh,zzhi,zvalu,zmargin,zuse,ztola
c-------------------------------
c
      ierr=0
      id_axis_use=id_axis
      if(id_axis_use.lt.0) id_axis_use=-id_axis_use
c
      if((id_axis_use.le.0).or.(id_axis_use.gt.naxes)) then
         write(lunerr,*) ' ?eqi_binget:  invalid axis id code:  ',
     >      id_axis_use
         ierr=1
         return
      endif
c
      iaddr=axes(id_axis_use)%address
      isize=axes(id_axis_use)%size
c
      if(axes(id_axis_use)%lin.eq.0) then
         call r8xlookup(ivec,zval,isize,eqbuf(iaddr),2,
     >      indx,zparam,zh,zhi,ierr)
      else if(axes(id_axis_use)%per.eq.0) then
c  linear, w/bounds check
         ilast=iaddr+isize-1
         do i=1,ivec
            zuse=max(eqbuf(iaddr),min(eqbuf(ilast),zval(i)))
            if(zuse.ne.zval(i)) then
               ztola=max(ceps10,min(ceps5,axes(id_axis_use)%tol))*
     >            max(abs(eqbuf(iaddr)),abs(eqbuf(ilast)))
               if(abs(zuse-zval(i)).gt.ztola) ierr=ierr+1
            endif
            zndx=cone+(isize-1)*(zuse-eqbuf(iaddr))*
     >         axes(id_axis_use)%deli
            indx(i)=min((isize-1),int(zndx))
            zparam(i)=min(cone,(zndx-indx(i)))
            zh(i)=axes(id_axis_use)%hav
            zhi(i)=axes(id_axis_use)%havi
         enddo
      else
c  linear, periodic
         ilast=iaddr+isize-1
         do i=1,ivec
            if((zval(i).lt.eqbuf(iaddr)).or.
     >         (zval(i).gt.eqbuf(ilast))) then
               zuse=mod((zval(i)-eqbuf(iaddr)),
     >            (eqbuf(ilast)-eqbuf(iaddr)))
               if(zuse.lt.czero) zuse=zuse+(eqbuf(ilast)-eqbuf(iaddr))
               zuse=zuse+eqbuf(iaddr)
               zuse=max(eqbuf(iaddr),min(eqbuf(ilast),zuse))
            else
               zuse=zval(i)
            endif
            zndx=cone+(isize-1)*(zuse-eqbuf(iaddr))*
     >         axes(id_axis_use)%deli
            indx(i)=min((isize-1),int(zndx))
            zparam(i)=min(cone,(zndx-indx(i)))
            zh(i)=axes(id_axis_use)%hav
            zhi(i)=axes(id_axis_use)%havi
         enddo
      endif
c
      if((ierr.gt.0).and.(id_axis_use.gt.0)) then
         write(lunerr,*) ' %eqi_binget:  ',ierr,' points out of range.'
      print *,'eqbuf(iaddr) ',eqbuf(iaddr:ilast)
      print *,'zval =',zval
      endif
c
      return
      end
