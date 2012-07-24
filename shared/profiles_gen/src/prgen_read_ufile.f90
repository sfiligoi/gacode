!--------------------------------------------------------------
! prgen_read_ufile.f90
!
! PURPOSE:
!  Extract data from UFILE using external library (read
!--------------------------------------------------------------

subroutine prgen_read_ufile

  use prgen_read_globals
  use legacyread_interface

  implicit none

  integer :: i

  character (len=6) :: phase
  character (len=50) :: cudir
  real :: xp_time
  real :: endtime
  integer :: time_series
  integer :: itorque
  integer :: iptot
  integer :: nsmooth
  integer :: idatzero
  integer :: iproc_d

  open(unit=1,file=raw_data_file,status='old')
  read(1,*) ufile_tok
  read(1,*) ufile_shot
  read(1,*) xp_time ! 3.09
  read(1,*) endtime ! 1000.0
  read(1,*) itorque
  read(1,*) iptot
  read(1,*) ufile_nj
  close(1)

  phase = ''
  cudir = './'
  time_series = 0
  nsmooth     = 0
  idatzero    = 0
  iproc_d     = 0

  call readufiles(&
       ufile_tok,&
       ufile_shot,&
       phase,&
       cudir,&
       xp_time,&
       endtime,&
       time_series,&
       itorque,&
       iptot,&
       ufile_nj,&
       nsmooth,&
       idatzero,&
       iproc_d)

  ufile_nj = ufile_nj+1
  nx       = ufile_nj

  ufile_nion  = nion_d
  ufile_nprim = nprim_d
  ufile_nimp  = nimp_d
  ufile_ibion = ibion_d

  call allocate_internals
  call allocate_ufile_vars

  do i=1,nx
     ! Normalized rho
     rho(i)   = rho_d(i)
     rmin(i)  = rminavnpsi_d(i)
     rmaj(i)  = rmajavnpsi_d(i)
     q(i)     = q_d(i)
     kappa(i) = elongx_d(i)  
     delta(i) = deltax_d(i)
     dpsi(i)  = psir_d(i)-psir_d(1)

     ! Unnormalized rho in m.
     ufile_rho(i)  = r_d(i)
     ufile_te(i)   = te_d(i)
     ufile_ti(i)   = ti_d(i)
     ufile_zeff(i) = zeff_d(i)
     ufile_ene(i)  = ene_d(i)
     ufile_en(i,1:ufile_nion) = en_d(i,1:ufile_nion)
     ufile_enbeam(i) = enbeam_d(i)
     ufile_pfast(i) = pfast_d(i)
     ufile_ptot(i) = ptot_d(i)
     ufile_angrot(i) = angrot_d(i)
  enddo

  ufile_bref = -btor_d

end subroutine prgen_read_ufile

