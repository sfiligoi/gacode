!---------------------------------------------------------
! Dep_grid_wts.f90
!
! Called by Dep_kernal.f90  to compute  e_hat(ie), lambda(k), sig(isig),
!  and velocity space weights e_wts(ie), lambda_wts(k)
!  and Fmax(ie)
! 
!-------------------------------------------------------------

subroutine  Dep_grid_wts

    use Dep_global           !e_hat(ie),lambda(k),sig(isig)
                             !e_wts(ie),lambda_wts(k),Fmax(ie)
                             !ie_max,ie,k_max,n_max,nb_max
                             !D_EP_starOchi_i_kernal(ie,k,isig,n,nb)
                             !A_EP(ie,k,isig,n,nb)

    use Dep_from_tglf        !krho(n), omega_tglf(n,nb)
                             !chi_i_tglf_wt(n,nb), chi_i_tglf
                             !r_hat,rmaj_hat,q_saf, Ti_hat, aoLT_i, aoLn_i

!
  implicit none

    integer :: n
    integer :: ie
    integer :: k
    integer :: isig

    integer :: i_chk_print

   
   real :: lambda_m
   real :: sum_energy_wts
   real :: sum_lambda_wts
   real :: pi


   i_chk_print =1
  !i_chk_print =0

!!!   print *, 'start inside Dep_grid_wts'
!!!   print *,n_max,nb_max,ie_max,k_max

   pi=3.14159265

   ! make energy grids
     e_hat(1) = 0.1086724547
     e_hat(2) = 0.3531529368  
     e_hat(3) = 0.7008934061 
     e_hat(4) = 1.1731964060  
     e_hat(5) = 1.8231909594  
     e_hat(6) = 2.7751284427  
     e_hat(7) = 4.3805639357  
     e_hat(8) = 6.0000000000  

   ! make energy wts
     e_wts(1) = 0.0642644789
     e_wts(2) = 0.1388201408
     e_wts(3) = 0.1895054690
     e_wts(4) = 0.2074366620
     e_wts(5) = 0.1895054690
     e_wts(6) = 0.1388201408
     e_wts(7) = 0.0642644789
     e_wts(8) = 0.0073831605

   ! make Fmax
     do ie=1,8
      Fmax(ie)=1./(pi*sqrt(pi))*exp(-e_hat(ie))
     enddo
!  Fmax is actually not used anywhere except Fmax(ie)/Fmax(ie)

   !test sum
    sum_energy_wts=0.
    do ie=1,8
     print *,'ie=',ie,'e_hat=',e_hat(ie),'e_wts=',e_wts(ie) 
     sum_energy_wts=sum_energy_wts+e_wts(ie)
    enddo
    if(i_chk_print .eq. 1) &
    print *, 'sum_energy_wts=',sum_energy_wts

   !note:  e_wts(ie) has Fmax(ie) built-in, so delta_f --->delta_f/Fmax                      


   !make lambda grids and wts
    lambda_m = (1.-r_hat/rmaj_hat)/(1.+r_hat/rmaj_hat)
        if(i_chk_print .eq. 1) then
     print *, 'lambda_m =',lambda_m 
     print *, 'passing wts=',1.0-sqrt(1-lambda_m)
     print *, 'trapped wts=', sqrt(1-lambda_m)
        endif

    !passing
    do k=1,4
     lambda(k)=lambda_m/8.*(2.*REAL(k)-1.0)
     lambda_wts(k) = sqrt(1.-(lambda_m/4.)* &
        (REAL(k)-1.0))- &
                     sqrt(1.-(lambda_m/4.)*REAL(k))
    if(i_chk_print .eq. 1) &
      print *, 'lambda=',lambda(k),'  lambda_wts=',lambda_wts(k)
    enddo

    !trapped
    do k=5,8
     lambda(k)=(1.-lambda_m)/8.*(2.*(REAL(k)-4.)-1.0) +lambda_m
     lambda_wts(k) = sqrt(1.-((1.-lambda_m)/4.)* &
        ((REAL(k)-4.0)-1.0)-lambda_m)- &
                     sqrt(1.-((1.-lambda_m)/4.)*(REAL(k)-4.)-lambda_m)
    if(i_chk_print .eq. 1) &
     print *, 'lambda=',lambda(k),'  lambda_wts=',lambda_wts(k)
    enddo
    !test sum
    sum_lambda_wts=0.
    do k=1,4
     sum_lambda_wts=sum_lambda_wts+lambda_wts(k)
    enddo
    if(i_chk_print .eq. 1) &
     print *, 'passing sum_lambda_wts=',sum_lambda_wts
    do k=5,8
     sum_lambda_wts=sum_lambda_wts+lambda_wts(k)
    enddo

    if(i_chk_print .eq. 1) then 
     print *, 'total sum_lambda_wts=',sum_lambda_wts
     lambda_wts(:)=lambda_wts(:)/2.
     print *, 'lambda_wts divided by 2. to compensate isig + and - sum'
     print *, 'Sum_isig Sum_ie Sum_k e_wts(ie)*lambda_wts(k) = 1.0 = Int dv**3*Fmax'
    endif 

     print *, 'Dep_grid_wts done'

end  subroutine Dep_grid_wts
