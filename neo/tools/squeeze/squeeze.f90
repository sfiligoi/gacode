program squeeze

  implicit none

  integer, parameter :: deriv_flag = 2

  integer :: i, j
  integer :: nspec 
  integer :: ne
  integer :: nxi
  integer :: ntheta
  integer :: nr
  real, dimension(:), allocatable :: r
  real, dimension(:), allocatable :: dphi0
  real, dimension(:), allocatable :: q
  real, dimension(:), allocatable :: rho
  real, dimension(:), allocatable :: rmaj
  real, dimension(:), allocatable :: den
  real, dimension(:), allocatable :: theta
  real, dimension(:), allocatable :: df
  real, dimension(:), allocatable :: ipsi
  real, dimension(:), allocatable :: aln
  real, dimension(:), allocatable :: a0
  real, dimension(:), allocatable :: capf
  real, dimension(:), allocatable :: capg
  real, dimension(:), allocatable :: caph
  real, dimension(:), allocatable :: capf_p
  real, dimension(:), allocatable :: v1
  real, dimension(:), allocatable :: v3
  real, dimension(:), allocatable :: v2

  real, dimension(:,:), allocatable :: vec
  real, dimension(:,:), allocatable :: wd

  open(unit=1,file='grid.out',status='old')
  read(1,*) nspec
  read(1,*) ne
  read(1,*) nxi
  read(1,*) ntheta
  allocate(theta(ntheta))
  read(1,*) theta(:)
  read(1,*) nr
  allocate(r(nr))
  read(1,*) r(:)
  close(1)

  allocate(vec(10+5*nspec,nr))

  open(unit=1,file='equil.out',status='old')
  read(1,*) vec(:,:)
  close(1)

  allocate(dphi0(nr))
  allocate(q(nr))
  allocate(rho(nr))
  allocate(rmaj(nr))
  allocate(den(nr))
  allocate(df(nr))
  allocate(ipsi(nr))
  allocate(capf(nr))
  allocate(capg(nr))
  allocate(caph(nr))
  allocate(capf_p(nr))
  allocate(aln(nr))
  allocate(a0(nr))
  allocate(v1(nr))
  allocate(v3(nr))
  allocate(v2(nr))

  dphi0(:) = vec(2,:)
  q(:)     = vec(3,:)
  rho(:)   = vec(4,:)
  rmaj(:)  = vec(5,:)
  den(:)   = vec(11,:) * vec(8,:)

  !q = 3.0
  !den = exp(-r*r)
  !print *, den(:)-exp(-r**2)

  ! I/psi'
  ipsi(:) = q*rmaj/r

  if(deriv_flag == 2) then
     allocate(wd(nr,nr))
     open(unit=1,file='curtis.out',status='old')
     do i=1,nr
        do j=1,nr
           read(1,*) wd(i,j)
        enddo
     enddo
     close(1)
  endif

  if(deriv_flag == 2) then
     df(:) = 0.0
     do i=1,nr
        do j=1,nr
           df(i) =  df(i) + wd(i,j) * den(j)
        enddo
     enddo
  else if(deriv_flag == 3) then
     call cub_spline_deriv2(r,den,nr,df)
  else 
     call bound_deriv(df,den,r,nr)
  endif

  aln = -df/den

  a0 = -dphi0(:)+aln

  capf = a0*den*ipsi

  if(deriv_flag == 2) then
     df(:) = 0.0
     do i=1,nr
        do j=1,nr
           df(i) =  df(i) + wd(i,j) * capf(j)
        enddo
     enddo
  else if(deriv_flag == 3) then
     call cub_spline_deriv2(r,capf,nr,df)
  else
     call bound_deriv(df,capf,r,nr)
  endif

  capf_p = df

  capg = (df+dphi0*capf)*ipsi

  if(deriv_flag == 2) then
     df(:) = 0.0
     do i=1,nr
        do j=1,nr
           df(i) =  df(i) + wd(i,j) * capg(j)
        enddo
     enddo
  else if(deriv_flag == 3) then
     call cub_spline_deriv2(r,capg,nr,df)
  else
     call bound_deriv(df,capg,r,nr)
  endif

  caph = (df+dphi0*capg)*ipsi

  capf = capf/den
  capg = capg/den
  caph = caph/den
  capf_p = capf_p/den

  v1 = capf
  v2 = capg
  v3 = caph

  print '(a)','   r/a         F/rho_*     G/rho_*^2  H/rho_*^3'  
  print '(a)',' -------------------------------------'
  do i=1,nr
     print '(5(1pe12.5,1x))',r(i),v1(i),v2(i), v3(i)
  enddo

  open(unit=1,file='squeeze.out',status='replace')
  do i=1,nr
     write(1,'(5(1pe12.5,1x))') r(i),v1(i),v2(i), v3(i)
  enddo
  close(1)

  

end program squeeze
