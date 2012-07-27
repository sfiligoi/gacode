!---------------------------------------------------------
! Dep_mainsub.f90
!
!  Called by NUBEAM to get D_EP_starOchi_i_TGLF(ie,k,isig) (and cA_EP) 
!      at each local r_hat
!
!            NUBEAM supplies local T_EP_hat
!
!--------------------------------------------------------
! PURPOSE:
!  Entry from kinetic EP transport code like NUBEAM to compute the 
!   EP diffusivity D_EP(r_hat,e_hat,lambda,sigma) from TGLF quasilinear ratio 
!     M_EP  is an 
!                rr  re
!                er  ee
!             positive definite matrix
!
!   D_EP  = D_EP_star*M_EP 
!           
!           M_EP_rr = 1         M_EP_re = -c A_EP
!           M_EP_er = c -A_EP     M_EP_ee= c A_EP**2
!   "c" reminds there is an n-mode (c)onvolution involved
!
!   Main OUTPUT to NUBEAM:       D_EP_starOchi_i_TGLF ( and A_EP  later)
!   
!          typical call  D_EP = D_EP_starOchi_i_TGLF*M_EP*chi_i_TGLF 
!                or
!                    D_EP = D_EP_starOchi_i_TGLF*M_EP*chi_exp
! 
!                    until NUBEAM can resolve  -d f_EP/d e_hat
!                    use  only M_EP_rr part
!
!                    D_EP_starOchi_i_TGLF=D_EP_star/chi_i_TGLF
! 
!                    A_EP ~ 2*Ln_EP/R depends on -d lnf_EP/dr_hat
!                                     and set to 0 until later.
! 
!   INPUTS:
!     local plasma:   r_hat = r/a   minor radius
!                         q         safety factor 
!                     R_hat =R0/a   major radius flux surface center
!                     Ti_hat = T0_i/T0_e
!                     aoLT_i = a/LT_i  = -d ln T_i/dr_hat
!                     aoLn_i = a/Ln_i  = -d ln n_i/dr_hat
!
!     local EP:       T_EP_hat = 2/3*ave_EP_energy/T0_e
!                       where
!                      ave_EP_energy = [Int dv3  m*v**2/2 * f_EP ]
!                                      ---------------------------
!                                      [Int dv3  1* f_EP ]
! 
!                     aoLn_EP  = -d ln n_EP/d r_hat
!
!    Some control parameters:   ??
!
!    Some grid parameters:  hardwired in Dep_globla
!                           ie_max = 8 k_max = 8 (4 pass & 4 trapped) 
!                           isig_max = 2 sign(1) = 1. & sign(2) = -1. 
!  
!
!                                                   n_max = 16  or TGLF std
!                          
!  OUTPUT: 
!    grid   EP_hat(ie) = T_EP_hat*e_hat(ie)  with e_hat(ie) hardwired
!                                             to cover thermal ion Maxwellian
!                 lambda(k)    hardwired from formula
!
!  Dep_kernal and TGLF INPUTS:
!
!           D_EP_starOchi_i_TGLF=Sum_n Sun_nb 
!                              D_EP_starOchi_i_kernal(n,nb)*chi_tglf_wt(n,nb)
!             
!                 chi_i_tglf_wt(n,nb) = chi_i_tglf(n,nb)/chi_i_tglf
!
!
!                   D_EP_starOchi_i_kernal(n,nb) provided by 
!                                        call to Dep_kernal.f90
! 
!                      chi_i_tglf_wt(n,nb) provide by call to TGLF
!
!        n is krho grid lable & nb is branch lable
!                 
!---------------------------------------------------------

subroutine Dep_mainsub

    use from_tglf_to_Dep     !krho(n), omega_tglf(n,nb)
                             !chi_i_tglf_wt(n,nb), chi_i_tglf
                             !r_hat,rmaj_hat,q_saf, Ti_hat, aoLT_i, aoLn_i

    use Dep_global           !e_hat(ie),lambda(k),sig(isig)
                             !e_wts(ie),lambda_wts(k),Fmax(ie)
                             !ie_max,ie,k_max,n_max,nb_max
                             !D_EP_starOchi_i_kernal(ie,k,isig,n,nb)
                             !A_EP(ie,k,isig,n)

    use from_nubeam_to_Dep   !T_EP_hat, aoLf_EP
    
    use from_Dep_to_nubeam   !D_EP_starOchi_i_TGLF(ie,k,isig)

!
  implicit none

    integer :: n
    integer :: nb
    integer :: ie
    integer :: k
    integer :: isig
 
    real :: D_EPOchi_i_TGLF

   
    call Dep_kernel !calls Dep_grid_wts

    D_EP_starOchi_i_TGLF(:,:,:) = 0.
 
    do n=1,n_max
     do nb=1,nb_max
      D_EP_starOchi_i_TGLF(:,:,:)=D_EP_starOchi_i_TGLF(:,:,:)+ &
           D_EP_starOchi_i_kernal(:,:,:,n,nb)*chi_i_tglf_wt(n,nb)
     enddo
    enddo

! check on Maxwellian  D_EPOchi_i_TGLF  assuming  aoLT_EP = 0.

   D_EPOchi_i_TGLF=0.0

   do n = 1,n_max
    do nb = 1,nb_max
     do ie = 1,ie_max
      do k = 1,k_max
       do isig = 1,2
         D_EPOchi_i_TGLF=D_EPOchi_i_TGLF+ &
           D_EP_starOchi_i_kernal(ie,k,isig,n,nb)*(1. - A_EP(ie,k,isig,n))*&
              Fmax(ie)*chi_i_tglf_wt(n,nb)
       enddo !isig
      enddo !k
     enddo !ie
    enddo !nb
   enddo !n

    print *, 'r_hat=',r_hat,'  D_EPOchi_i_TGLF=',D_EPOchi_i_TGLF
     
    


end subroutine Dep_mainsub
