      SUBROUTINE gftm_deallocate
!
! deallocate internal modules
!
!-------------------------------------------------
      USE gftm_species
      USE gftm_coeff
      USE gftm_GFS
      USE gftm_eigen
      IMPLICIT NONE
!
!
!     MODULE gftm_species
! 
! species parameters
!
      IF(ALLOCATED(ei_exch))DEALLOCATE(ei_exch)
      IF(ALLOCATED(resist))DEALLOCATE(resist)
      IF(ALLOCATED(zs))DEALLOCATE(zs)
      IF(ALLOCATED(mass))DEALLOCATE(mass)
      IF(ALLOCATED(vs))DEALLOCATE(vs)
      IF(ALLOCATED(fts))DEALLOCATE(fts)
      IF(ALLOCATED(rlts))DEALLOCATE(rlts)
      IF(ALLOCATED(rlns))DEALLOCATE(rlns)
      IF(ALLOCATED(vpar_shear_s))DEALLOCATE(vpar_shear_s)
      IF(ALLOCATED(as))DEALLOCATE(as)
      IF(ALLOCATED(taus))DEALLOCATE(taus)
      IF(ALLOCATED(vpar_s))DEALLOCATE(vpar_s)

!      write(*,*)"deallocated species"
!
!
!-------------------------------------------------
!
!      MODULE gftm_coeff
!
!  ave_theta
!
      IF(ALLOCATED(ave_kx))DEALLOCATE(ave_kx)
      IF(ALLOCATED(ave_c_tor_par))DEALLOCATE(ave_c_tor_par)
      IF(ALLOCATED(ave_c_tor_per))DEALLOCATE(ave_c_tor_per)
      IF(ALLOCATED(ave_c_par_par))DEALLOCATE(ave_c_par_par)
      IF(ALLOCATED(ave_wdpar))DEALLOCATE(ave_wdpar)
      IF(ALLOCATED(ave_wdper))DEALLOCATE(ave_wdper)
      IF(ALLOCATED(ave_gradB))DEALLOCATE(ave_gradB)
      IF(ALLOCATED(ave_lnB))DEALLOCATE(ave_lnB)
      IF(ALLOCATED(ave_b0))DEALLOCATE(ave_b0)
      IF(ALLOCATED(ave_b0inv))DEALLOCATE(ave_b0inv)
      IF(ALLOCATED(ave_kpar))DEALLOCATE(ave_kpar)
      IF(ALLOCATED(ave_modkpar))DEALLOCATE(ave_modkpar)
      IF(ALLOCATED(ave_p0))DEALLOCATE(ave_p0)
      IF(ALLOCATED(ave_p0inv))DEALLOCATE(ave_p0inv)
      IF(ALLOCATED(ave_bp))DEALLOCATE(ave_bp)
      IF(ALLOCATED(ave_bpinv))DEALLOCATE(ave_bpinv)
!      write(*,*)"deallocated ave_theta"
!
!
!-------------------------------------------------
!
!      MODULE gftm_GFS
!
! allocate gftm_GFS
      IF(ALLOCATED(zomega))DEALLOCATE(zomega)
      IF(ALLOCATED(matmirror))DEALLOCATE(matmirror)
      IF(ALLOCATED(matu))DEALLOCATE(matu)
      IF(ALLOCATED(matuc))DEALLOCATE(matuc)
      IF(ALLOCATED(matdu))DEALLOCATE(matdu)
      IF(ALLOCATED(matuu))DEALLOCATE(matuu)
      IF(ALLOCATED(mate))DEALLOCATE(mate)
      IF(ALLOCATED(matde))DEALLOCATE(matde)
      IF(ALLOCATED(phib))DEALLOCATE(phib)
      IF(ALLOCATED(psib))DEALLOCATE(psib)
      IF(ALLOCATED(sigb))DEALLOCATE(sigb)
      IF(ALLOCATED(pe1j0phi))DEALLOCATE(pe1j0phi)
      IF(ALLOCATED(pe1j0psi))DEALLOCATE(pe1j0psi)
      IF(ALLOCATED(pe1j1sig))DEALLOCATE(pe1j1sig)
      IF(ALLOCATED(pe2j0phi))DEALLOCATE(pe2j0phi)
      IF(ALLOCATED(pe2j0psi))DEALLOCATE(pe2j0psi)
      IF(ALLOCATED(pe2j1sig))DEALLOCATE(pe2j1sig)
      IF(ALLOCATED(pe1j0phib))DEALLOCATE(pe1j0phib)
      IF(ALLOCATED(pe1j0psib))DEALLOCATE(pe1j0psib)
      IF(ALLOCATED(pe1j1sigb))DEALLOCATE(pe1j1sigb)
      IF(ALLOCATED(pe2j0phib))DEALLOCATE(pe2j0phib)
      IF(ALLOCATED(pe2j0psib))DEALLOCATE(pe2j0psib)
      IF(ALLOCATED(pe2j1sigb))DEALLOCATE(pe2j1sigb)
      IF(ALLOCATED(pu1pe1PSI))DEALLOCATE(pu1pe1PSI)
      IF(ALLOCATED(pu1pe2PSI))DEALLOCATE(pu1pe2PSI)
      IF(ALLOCATED(pu3pe1PSI))DEALLOCATE(pu3pe1PSI)
      IF(ALLOCATED(pu2pe1PSI))DEALLOCATE(pu2pe1PSI)
      IF(ALLOCATED(p0invpe1j0s))DEALLOCATE(p0invpe1j0s)
      IF(ALLOCATED(b0invpe1j0s))DEALLOCATE(b0invpe1j0s)
      IF(ALLOCATED(mateq))DEALLOCATE(mateq)
      IF(ALLOCATED(mats))DEALLOCATE(mats)
      IF(ALLOCATED(matb))DEALLOCATE(matb)
!
! dealocate eigenvalues, eigenvectors
!
!    write(*,*)"ready to deallocate gftm_LS"
!  gftm_LS local
      if(ALLOCATED(di))DEALLOCATE(di)
      if(ALLOCATED(de))DEALLOCATE(de)
      if(ALLOCATED(ipive))DEALLOCATE(ipive)
      if(ALLOCATED(ipiv))DEALLOCATE(ipiv)
      if(ALLOCATED(emat))DEALLOCATE(emat)
      if(ALLOCATED(zmat))DEALLOCATE(zmat)
      if(ALLOCATED(gmat))DEALLOCATE(gmat)
      if(ALLOCATED(gmatinv))DEALLOCATE(gmatinv)
!      write(*,*)"deallocated gftm_LS local"
!  gftm_eigen module
      if(ALLOCATED(grow_index))DEALLOCATE(grow_index)
      if(ALLOCATED(fv1))DEALLOCATE(fv1)
      if(ALLOCATED(fv2))DEALLOCATE(fv2)
      if(ALLOCATED(fv3))DEALLOCATE(fv3)
      if(ALLOCATED(rr))DEALLOCATE(rr)
      if(ALLOCATED(ri))DEALLOCATE(ri)
      if(ALLOCATED(ar))DEALLOCATE(ar)
      if(ALLOCATED(ai))DEALLOCATE(ai)
      if(ALLOCATED(vr))DEALLOCATE(vr)
      if(ALLOCATED(vi))DEALLOCATE(vi)
      if(ALLOCATED(amat))DEALLOCATE(amat)
      if(ALLOCATED(bmat))DEALLOCATE(bmat)
      if(ALLOCATED(he))DEALLOCATE(he)
      if(ALLOCATED(hetot))DEALLOCATE(hetot)
      if(ALLOCATED(alpha))DEALLOCATE(alpha)
      if(ALLOCATED(beta))DEALLOCATE(beta)
!
      if(ALLOCATED(phi))DEALLOCATE(phi)
      if(ALLOCATED(psi))DEALLOCATE(psi)
      if(ALLOCATED(sig))DEALLOCATE(sig)
      if(ALLOCATED(pwns))DEALLOCATE(pwns)
      if(ALLOCATED(pwps))DEALLOCATE(pwps)
      if(ALLOCATED(pwtpars))DEALLOCATE(pwtpars)
      if(ALLOCATED(pwtpers))DEALLOCATE(pwtpers)
      if(ALLOCATED(pwes))DEALLOCATE(pwes)
      if(ALLOCATED(hes))DEALLOCATE(hes)
      if(ALLOCATED(hwns))DEALLOCATE(hwns)
      if(ALLOCATED(hwps))DEALLOCATE(hwps)
      if(ALLOCATED(hwtpars))DEALLOCATE(hwtpars)
      if(ALLOCATED(hwtpers))DEALLOCATE(hwtpers)
      if(ALLOCATED(hwes))DEALLOCATE(hwes)
      if(ALLOCATED(PSIpu1pe1))DEALLOCATE(PSIpu1pe1)
      if(ALLOCATED(PSIpu2pe1))DEALLOCATE(PSIpu2pe1)
      if(ALLOCATED(PSIpu3pe1))DEALLOCATE(PSIpu3pe1)
      if(ALLOCATED(PSIpu1pe2))DEALLOCATE(PSIpu1pe2)
      if(ALLOCATED(rwn))DEALLOCATE(rwn)
      if(ALLOCATED(rwp))DEALLOCATE(rwp)
      if(ALLOCATED(rwtpar))DEALLOCATE(rwtpar)
      if(ALLOCATED(rwtper))DEALLOCATE(rwtper)
      if(ALLOCATED(rwe))DEALLOCATE(rwe)
      if(ALLOCATED(fluxe))DEALLOCATE(fluxe)
!
!      write(*,*)"deallocated eigen"
!
      END SUBROUTINE gftm_deallocate
!-----------------------------------------------------
!
