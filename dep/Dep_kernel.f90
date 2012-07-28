!---------------------------------------------------------
! Dep_kernel.f90
!
! Called by Dep_mainsub.f90  to compute  D_EP_starOchi_i_kernal(ie,k,isig,n,nb)
! from TGLF output omega_tglf(n,nb)
! 
! or
!
! Called by Dep_driver.f90 for testing
!
!-------------------------------------------------------------
! INPUT from TGLF (or Dep_driver):
!  r_hat, rmaj_hat, q_saf, T_EP_hat, Ti_hat, aoLT_i, aoLn_i
!  krho(n), omega_tglf(n,nb)
!
! OUTPUT:   D_EP_starOchi_i_kernal(ie,k,isig,n,nb)
!           A_EP(ie,k,isig,n) 
!                 
!---------------------------------------------------------

subroutine Dep_kernel

    use from_tglf_to_Dep     !krho(n), omega_tglf(n,nb)
                             !chi_i_tglf_wt(n,nb), chi_i_tglf
                             !r_hat,rmaj_hat,q_saf, Ti_hat, aoLT_i, aoLn_i

    use Dep_global           !e_hat(ie),lambda(k),sig(isig)
                             !e_wts(ie),lambda_wts(k),Fmax(ie)
                             !ie_max,ie,k_max,n_max,nb_max
                             !D_EP_starOchi_i_kernal(ie,k,isig,n,nb)
                             !A_EP(ie,k,isig,n)

    use from_nubeam_to_Dep   !T_EP_hat, aoLf_EP
                             
!
  implicit none

    integer :: n
    integer :: nb
    integer :: ie
    integer :: k
    integer :: isig

    real :: Z_EP
    real :: mu_EP
    real :: aoLn_EP

    real :: krho_i
    real :: j0_i
    real :: krho_EP
    real :: j0_EP

    real, dimension(k_max) :: k_par
    real :: k_par_fit
    real :: Mysign(2)

    real :: omega_d_fit
    real, dimension(ie_max) :: omega_star_i
    real, dimension(ie_max,k_max,2,n_max) :: omega_db_i
    real, dimension(ie_max,k_max,2,n_max) :: omega_db_EP

    real, dimension(ie_max,k_max,2,n_max,nb_max) :: D_chi_i_den_part

    real :: pi
    complex :: xi

    real, external :: BESJ0

    pi = 3.14159265
    xi = (0.,1.0)
 
!    j0_i = BESJ0(0.1)
!    print *, 'j0 start=',j0_i

    mu_EP = 1.   !  mu_EP = sqrt(m_i/m_EP)  must reset if neede
    Z_EP  = 1.   !  Z_EP  = z_Ep/z_i    must reset if neede 


    aoLn_EP =1.  !default


    omega_d_fit = 1.0  ! to be adjusted
   
    k_par_fit = 1.0  ! to be adjusted 


    k_par(:) = 0.0
    do k=1,4
     k_par(k) = k_par_fit/rmaj_hat/q_saf 
    enddo
    ! trapped k_par = 0 from bounce averaging
    Mysign(1) = 1.
    Mysign(2) = -1.
  
    call Dep_grid_wts

!   krho(n) from TGLF
!   omega_tglf(n) from TGLF

    D_EP_starOchi_i_kernal(:,:,:,:,:) = 0.0
!         D_EP_starOchi_i_kernal(ie,k,isig,n)
    A_EP(:,:,:,:) = 0.0
!         A_EP(ie,k,isig,n)  does NOT depend on nb 
!                            since no omega_tglf dependence, but
!                            does depend on n via k_par/krho part of omega_db
!   which means the A_EP d/de parts must be convolved with chi_i_tglf_wt 

   D_EP_starOchi_i_kernal_den=0.0

    do n=1,n_max
     do nb=1,nb_max
     
      do ie=1,ie_max
       do k=1,k_max
        do isig=1,2
         

     krho_i = krho(n)*sqrt(Ti_hat)*sqrt(2.*e_hat(ie)*lambda(k))
     !!!call CALJY0(krho_i,j0_i,0)
     j0_i = BESJ0(krho_i)

     krho_EP = sqrt(T_EP_hat/Ti_hat)/Z_EP/mu_EP*krho_i
     !!!call CALJY0(krho_EP,j0_EP,0)
     j0_EP = BESJ0(krho_EP)
   
     !test print
     !     print *, ie,k,isig,krho_i,krho_EP,j0_i,j0_EP
     !test print

     omega_star_i(ie) = krho(n)*(aoLn_i + aoLT_i*(e_hat(ie)-1.5))

     omega_db_i(ie,k,isig,n)=omega_d_fit*krho(n)*Ti_hat*2./Rmaj_hat*&
              (e_hat(ie)*(1.-lambda(k))+1./2.*e_hat(ie)*lambda(k)) & 
                  +k_par(k)*Mysign(isig)*sqrt(Ti_hat)*sqrt(2.*e_hat(ie)*(1.0-lambda(k)))

     omega_db_EP(ie,k,isig,n)=omega_d_fit*krho(n)*T_EP_hat/Z_EP*2./Rmaj_hat*&
              (e_hat(ie)*(1.-lambda(k))+1./2.*e_hat(ie)*lambda(k)) &
          +k_par(k)*Mysign(isig)*mu_EP*sqrt(T_EP_hat)*sqrt(2.*e_hat(ie)*(1.0-lambda(k)))

! EP numerator
!-----------------------------------------------------------------------------
    D_EP_starOchi_i_kernal(ie,k,isig,n,nb) = e_wts(ie)*lambda_wts(k)/Fmax(ie)* &
                j0_EP**2*Real(xi*krho(n)/ &
                (omega_tglf(n,nb)+omega_db_EP(ie,k,isig,n)))
!1/Fmax  since the v-space wts have a Maxwellian factor built-in & must remove
!
!
!    energy derivative off diagonal
!-----------------------------------------------------------------------------
!     aoLf_EP(ie,k,isig) = -d(ln f_Ep(ie,k,isig))/dr_hat

     if(aoLf_EP(ie,k,isig) .eq. 0.)  aoLf_EP(ie,k,isig) = aoLn_EP

       A_EP(ie,k,isig,n)=Z_EP/T_EP_hat*(omega_db_EP(ie,k,isig,n)/krho(n))/ &
                aoLf_EP(ie,k,isig)
!-----------------------------------------------------------------------------

    

! ion chi denominator
!-----------------------------------------------------------------------------
     D_chi_i_den_part(ie,k,isig,n,nb) =  e_hat(ie)* &
                                           e_wts(ie)*lambda_wts(k)* &
                j0_i**2*Real(xi*(omega_star_i(ie)-omega_db_i(ie,k,isig,n))/ &
                (omega_tglf(n,nb)+omega_db_i(ie,k,isig,n))/ &
                                                           aoLT_i)

     D_EP_starOchi_i_kernal_den = D_EP_starOchi_i_kernal_den + D_chi_i_den_part(ie,k,isig,n,nb)
!  D_EP_starOchi_i_kernal_den = chi_i(n,nb)  GYRO defined
!-----------------------------------------------------------------------------

        enddo ! isig
       enddo ! k
      enddo ! ie


     enddo !nb
    enddo !n
    
! divide numerator by denominator
!-----------------------------------------------------------------------------        

   print *, 'D_EP_starOchi_i_kernal_den=',D_EP_starOchi_i_kernal_den

   D_EP_starOchi_i_kernal(:,:,:,:,:) = D_EP_starOchi_i_kernal(:,:,:,:,:)/ &
                                             D_EP_starOchi_i_kernal_den 

   do n = 1,n_max
    do nb = 1,nb_max
     do k=1,k_max
  if(k .eq. 3 .or. k .eq. 6) then
      print *, '----------------------------------------------------------------'
    print *, ' lambda=',lambda(k)
      do ie = 1,ie_max
      print *, '----------------------------------------------------------------'
       do isig = 1,2
      print *, 'e_hat=',e_hat(ie),' sign=',Mysign(isig)

      print *, 'D_EP_star=',D_EP_starOchi_i_kernal(ie,k,isig,n,nb), ' A_EP=',A_EP(ie,k,isig,n)
      print *, 'D_chi_i_den_part=',D_chi_i_den_part(ie,k,isig,n,nb)
       enddo ! isig
      enddo ! ie
   endif
     enddo !k
    enddo !nb
   enddo !n

  print *, 'Dep_kernel Done'

end  subroutine Dep_kernel
