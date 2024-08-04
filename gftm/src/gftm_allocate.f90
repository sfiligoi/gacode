      SUBROUTINE gftm_allocate
!
! allocate internal modules
!
!-------------------------------------------------
      USE gftm_dimensions
      USE gftm_species
      USE gftm_coeff
      USE gftm_GFS
      USE gftm_eigen
      IMPLICIT NONE
!
!      write(*,*)"allocate",ns,nbasis,nx,nky
!
!     MODULE gftm_species
! 
! species parameters
!
      ALLOCATE(ei_exch(ns,ns))
      ALLOCATE(resist(ns,ns))
      ALLOCATE(zs(ns))
      ALLOCATE(mass(ns))
      ALLOCATE(vs(ns))
      ALLOCATE(fts(ns))
      ALLOCATE(rlts(ns))
      ALLOCATE(rlns(ns))
      ALLOCATE(vpar_shear_s(ns))
      ALLOCATE(as(ns))
      ALLOCATE(taus(ns))
      ALLOCATE(vpar_s(ns))         
!
!
!-------------------------------------------------
!
!      MODULE gftm_coeff
!
! store the hermite basis matrix coefficients
!
!  ave_theta
      ALLOCATE(ave_kx(nbasis,nbasis))
      ALLOCATE(ave_c_tor_par(nbasis,nbasis))
      ALLOCATE(ave_c_tor_per(nbasis,nbasis))
      ALLOCATE(ave_c_par_par(nbasis,nbasis))
      ALLOCATE(ave_wdpar(nbasis,nbasis))
      ALLOCATE(ave_wdper(nbasis,nbasis))
      ALLOCATE(ave_gradB(nbasis,nbasis))
      ALLOCATE(ave_lnB(nbasis,nbasis))
      ALLOCATE(ave_b0(nbasis,nbasis))
      ALLOCATE(ave_b0inv(nbasis,nbasis))
      ALLOCATE(ave_kpar(nbasis,nbasis))
      ALLOCATE(ave_modkpar(nbasis,nbasis))
      ALLOCATE(ave_p0(nbasis,nbasis))
      ALLOCATE(ave_p0inv(nbasis,nbasis))
      ALLOCATE(ave_bp(nbasis,nbasis))
      ALLOCATE(ave_bpinv(nbasis,nbasis))
!-------------------------------------------------
!
!      MODULE gftm_GFS
!
! allocate gftm_GFS
      ALLOCATE(zomega(ntot))
      ALLOCATE(matmirror(nune,nune))
      ALLOCATE(matu(nu,nu))
      ALLOCATE(matuc(nu,nu))
      ALLOCATE(matdu(nu,nu))
      ALLOCATE(matuu(nu,nu))
      ALLOCATE(mate(ne,ne))
      ALLOCATE(matde(ne,ne))
      ALLOCATE(phib(nbasis,ntot))
      ALLOCATE(psib(nbasis,ntot))
      ALLOCATE(sigb(nbasis,ntot))
      ALLOCATE(pe1j0phib(ns,ne,nbasis,ntot))
      ALLOCATE(pe1j0psib(ns,ne,nbasis,ntot))
      ALLOCATE(pe1j1sigb(ns,ne,nbasis,ntot))
      ALLOCATE(pe2j0phib(ns,ne,nbasis,ntot))
      ALLOCATE(pe2j0psib(ns,ne,nbasis,ntot))
      ALLOCATE(pe2j1sigb(ns,ne,nbasis,ntot))
      ALLOCATE(pe1j0phi(ns,ne,nbasis))
      ALLOCATE(pe1j0psi(ns,ne,nbasis))
      ALLOCATE(pe1j1sig(ns,ne,nbasis))
      ALLOCATE(pe2j0phi(ns,ne,nbasis))
      ALLOCATE(pe2j0psi(ns,ne,nbasis))
      ALLOCATE(pe2j1sig(ns,ne,nbasis))
      ALLOCATE(pu1pe1PSI(ntot,ntot))
      ALLOCATE(pu1pe2PSI(ntot,ntot))
      ALLOCATE(pu3pe1PSI(ntot,ntot))
      ALLOCATE(pu2pe1PSI(ntot,ntot))
      ALLOCATE(p0invpe1j0s(ns,ne,nbasis,nbasis))
      ALLOCATE(b0invpe1j0s(ns,ne,nbasis,nbasis))
      ALLOCATE(mateq(ntot,ntot))
      ALLOCATE(mats(ntot,ntot))
      ALLOCATE(matb(ntot,ntot))
!
!   gftm_LS local
      ALLOCATE(di(ntot))
      ALLOCATE(de(ntot))
      ALLOCATE(ipive(nphase))
      ALLOCATE(ipiv(ntot))
      ALLOCATE(emat(ntot,ntot))
      ALLOCATE(zmat(nphase,nphase))
      ALLOCATE(gmat(nphase,nphase))
      ALLOCATE(gmatinv(ns,nphase,nphase))
!   gftm_eigen module
      ALLOCATE(grow_index(ntot))
      ALLOCATE(fv1(ntot))
      ALLOCATE(fv2(ntot))
      ALLOCATE(fv3(ntot))
      ALLOCATE(rr(ntot))
      ALLOCATE(ri(ntot))
      ALLOCATE(ar(ntot,ntot))
      ALLOCATE(ai(ntot,ntot))
      ALLOCATE(vr(ntot,ntot))
      ALLOCATE(vi(ntot,ntot))
      ALLOCATE(amat(ntot,ntot))
      ALLOCATE(bmat(ntot,ntot))
      ALLOCATE(he(ntot))
      ALLOCATE(hetot(ntot,ntot))
      ALLOCATE(alpha(ntot))
      ALLOCATE(beta(ntot))
!
      ALLOCATE(phi(nbasis))
      ALLOCATE(psi(nbasis))
      ALLOCATE(sig(nbasis))
      ALLOCATE(pwns(ns,nphase))
      ALLOCATE(pwps(ns,nphase))
      ALLOCATE(pwtpars(ns,nphase))
      ALLOCATE(pwtpers(ns,nphase))
      ALLOCATE(pwes(ns,nphase))
      ALLOCATE(hes(ns,ntot))
      ALLOCATE(hwns(ns,ntot))
      ALLOCATE(hwps(ns,ntot))
      ALLOCATE(hwtpars(ns,ntot))
      ALLOCATE(hwtpers(ns,ntot))
      ALLOCATE(hwes(ns,ntot))
      ALLOCATE(PSIpu1pe1(ns,nphase))
      ALLOCATE(PSIpu2pe1(ns,nphase))
      ALLOCATE(PSIpu3pe1(ns,nphase))
      ALLOCATE(PSIpu1pe2(ns,nphase))
!
      ALLOCATE(rwn(ns,5))
      ALLOCATE(rwp(ns,5))
      ALLOCATE(rwtpar(ns,5))
      ALLOCATE(rwtper(ns,5))
      ALLOCATE(rwe(ns,5))
      ALLOCATE(fluxe(ns,5))
!
      END SUBROUTINE gftm_allocate
!-----------------------------------------------------
!
