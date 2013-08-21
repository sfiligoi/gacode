      MODULE psig
        USE nrtype,                         ONLY : I4B,DP
        USE mhdpar,                         ONLY : kpsi,nw
!
!--- quantities defined on the MHD grid (i.e., kpsi):
!--- note that these vectors are defined in reverse order from the
!--- corresponding quantitites defined on the rho grid
!--- (see INCLUDE file rhog.i)
!--- that is for example, psival(1) = psi on plasma boundary,
!--- psival(kpsi) = psi on magnetic axis!
!--- ratave = <1./R**2>, ratavei = <R**2>
!--- vprime = dv/dpsi
!--- qpsi = safety factor
!--- fpsi = f(psi)  (poloidal current FUNCTION)
!--- ravg = <R>
!--- ravgi= <1/R>   added for kinetic efit output
!--- bsqinvavg =<(Bt0/B)**2>
!--- bsq_avg =<B**2/Bt0**2>
!--- h_factr = <SQRT(1-B/Bmax)>
!--- elongx elongation of flux surfaces
!--- ffppsival = ffprime on psival grid,
!                computed in routine fofpsi (D.F.)
!
      REAL(DP)      psival(kpsi), rho(kpsi), ratave(kpsi),               &
                   vprime(kpsi), psivolp(kpsi), ffppsi(kpsi),           &
                   drhdt(kpsi), qpsi(kpsi), fpsi(kpsi), presspsi(kpsi), &
                   ratavei(kpsi), ravg(kpsi), pppsi(kpsi),elongx(kpsi), &
                   prseqdsk(nw), psieqdsk(nw), ffppsival(kpsi),         &
                   qpsival(kpsi),bsqinvavg(kpsi),bsq_avg(kpsi),         &
                   ravgi(kpsi),b_avg(kpsi),h_factr(kpsi) 
  
      INTEGER(I4B)  nxeqd
!
!--- the following is for the ITER database.
!--- It is included here rather than in iterdb.i because it is always
!--- calculated.
!
      REAL(DP) rmajavnpsi (kpsi), rminavnpsi(kpsi),                      &
                         triangnpsi (kpsi), sfareanpsi(kpsi),           &
                         grho1npsi  (kpsi), grho2npsi (kpsi),           &
                         pindentnpsi(kpsi), cxareanpsi(kpsi),           &
                         torfluxnpsi(kpsi),triangnpsi_l(kpsi)
!
!   rmajavnpsi   : avg major radius (average of min and max major radii
!                  of each flux surface at elevation of magnetic axis.)
!   rminavnpsi   : avg minor radius at elevation of magnetic axis
!   triangnpsi   : avg upper triangularity
!   triangnpsi_l : avg lower triangularity
!   pindentnpsi  : avg indentation
!   sfareanpsi   : surface area
!   torfluxnpsi  : toroidal flux
!   grho1npsi    : <ABS(grad rho)   >
!   grho2npsi    : <   (grad rho)**2>
!   b_avg(kpsi)  ; <B/Bt0>
!   h_factr(kpsi);  = <SQRT(1-B/Bmax)>
!
     END MODULE psig
