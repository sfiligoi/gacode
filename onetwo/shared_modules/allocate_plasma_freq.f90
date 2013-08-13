    MODULE pl_freq
       USE nrtype,                    ONLY : DP,I4B,I2B

       USE ions_gcnmp,                ONLY : nprimp1

       USE Plasma_properties ,        ONLY : plasma_frequencies

       USE vector_class,              ONLY : zero_vector

       USE grid_class,                ONLY : nj

       INTEGER(I4B),ALLOCATABLE,DIMENSION(:),SAVE ::  id_plas_freq
       INTEGER(I4B),ALLOCATABLE,DIMENSION(:),SAVE ::  id_ci_freq,id_lh_freq,id_uh_freq

    CONTAINS

    SUBROUTINE allocate_plasma_freq     
     !------------------------------------------------------------
     ! define arrays for common frequencies:
     ! define them and set to zero
     ! frequencies that are not flux functions are on the rho 
     !grid at zma=0  on the outboard side
     !-----------------------------------------------------------------

       IMPLICIT NONE
       INTEGER(I4B) i


     ! nprimp1 to separate 'dt' to 'd' and't' if necessary
     IF(.NOT. ASSOCIATED(plasma_frequencies%omega_pi))THEN 
          ALLOCATE (id_plas_freq(nprimp1),id_ci_freq(nprimp1), &
                    id_lh_freq(nprimp1),id_uh_freq(nprimp1))

          ALLOCATE (plasma_frequencies%omega_pi(nprimp1))
          ALLOCATE (plasma_frequencies%omega_ci(nprimp1))
          ALLOCATE (plasma_frequencies%omega_lh(nprimp1))
          ALLOCATE (plasma_frequencies%omega_uh(nprimp1))

          DO i=1,nprimp1
             plasma_frequencies%omega_pi(i) = zero_vector(nj)
             plasma_frequencies%omega_ci(i) = zero_vector(nj)
             plasma_frequencies%omega_lh(i) = zero_vector(nj)
             plasma_frequencies%omega_uh(i) = zero_vector(nj)
          ENDDO
          plasma_frequencies%omega_ce       = zero_vector(nj)
          plasma_frequencies%omega_pe       = zero_vector(nj)  
     ENDIF
        
     RETURN

    END SUBROUTINE allocate_plasma_freq 

    SUBROUTINE common_frequencies
!--------------------------------------------------------------------------
! determine rf and plasma fequencies
! NOTE:
!    Non flux surface average quantites are given on the outboard side
!    on the transport (rho) grid. 
! INPUT:
!     btotrmaj total B field on rmajor(rho) at z= mag axis, outborad side
!     rho_grid,ene,en,z,Atw
!     
!--------------------------------------------------------------------------
    USE grid_class,                                    ONLY : nj

    USE plasma_properties,                             ONLY : plasma_frequencies, &
                                                              mhd_dat,profile

    USE ions_gcnmp,                                    ONLY : namep,namei,nion,nprim, &
                                                              atw,z,nprimp1


    USE Vector_class,                                  ONLY : delete_vector,zero_vector


    INTEGER(I4B) j,i,k,lden
    REAL(DP) ne,B,ni,Ai,zi



   
    ! we use nprimp1 in case one in species is 'dt'
    IF(.NOT. ASSOCIATED(plasma_frequencies%omega_pi))THEN
         CALL allocate_plasma_freq
    ELSE
       DO i=1,nprimp1
          CALL delete_vector(plasma_frequencies%omega_pi(i))
          CALL delete_vector(plasma_frequencies%omega_ci(i))
          CALL delete_vector(plasma_frequencies%omega_lh(i))
          CALL delete_vector(plasma_frequencies%omega_uh(i))
       ENDDO

       DO i=1,nprimp1
          plasma_frequencies%omega_pi(i) = zero_vector(nj)
          plasma_frequencies%omega_ci(i) = zero_vector(nj)
          plasma_frequencies%omega_lh(i) = zero_vector(nj)
          plasma_frequencies%omega_uh(i) = zero_vector(nj)
       ENDDO
       plasma_frequencies%omega_ce       = zero_vector(nj)
       plasma_frequencies%omega_pe       = zero_vector(nj)
    ENDIF

    DO j=1,nj
       ne = profile%ene%data(j)
       B = mhd_dat%btotrmaj%data(j)
       k = 0
       DO i = 1,nprim
          k = k+1
          ni = profile%en(i)%data(j)
          Zi = z(j,i) 
          Ai = atw(i) 
          IF(namep(i) == 'dt')THEN ! 
             ni = ni*0.5_DP  ! assumes dt 50/50
             Ai  = 2._DP
             plasma_frequencies%omega_pi(k)%data(j) = w_pi(ni,Ai)
             plasma_frequencies%omega_ci(k)%data(j) = w_ci(Zi,Ai,B)
             plasma_frequencies%omega_lh(k)%data(j) = w_ci(Zi,Ai,B)
             plasma_frequencies%omega_uh(k)%data(j) = w_ci(Zi,Ai,B)
             k = k +1
             Ai  = 3._DP
             plasma_frequencies%omega_pi(k)%data(j) = w_pi(ni,Ai)
             plasma_frequencies%omega_ci(k)%data(j) = w_ci(Zi,Ai,B)
             plasma_frequencies%omega_lh(k)%data(j) = w_ci(Zi,Ai,B)
             plasma_frequencies%omega_uh(k)%data(j) = w_uh(B,ne)
          ELSE
             plasma_frequencies%omega_pi(k)%data(j) = w_pi(ni,Ai)
             plasma_frequencies%omega_ci(k)%data(j) = w_ci(Zi,Ai,B)
             plasma_frequencies%omega_lh(k)%data(j) = w_ci(Zi,Ai,B)
             plasma_frequencies%omega_uh(k)%data(j) = w_uh(B,ne)
          ENDIF
       ENDDO
       plasma_frequencies%omega_pe%data(j)     = w_pe(ne)
       plasma_frequencies%omega_ce%data(j)     = ABS(w_ce(B))

    ENDDO

     
    END SUBROUTINE common_frequencies



    REAL(DP) FUNCTION  w_ce(B) 
!--------------------------------------------------------------------------
!    electron cyclotron freq, rad/sec (negative) , B in Tesla
!--------------------------------------------------------------------------
     REAL(DP) B
       w_ce =  -0.176e12*B    ! - sign because of qe, B could also be negative
     RETURN
     END FUNCTION w_ce


    REAL(DP) FUNCTION  w_pe(ne) 
!--------------------------------------------------------------------------
! electron plasma freq,ne in m^-3, output in rad/sec
!--------------------------------------------------------------------------
     REAL(DP) ne 
      w_pe = 56.4_DP*SQRT(ne)
     RETURN
     END FUNCTION w_pe


     REAL(DP) FUNCTION w_ci(Zi,Ai,B) 
!----------------------------------------------------------------------------
!    ion cyclotron freq, rad/sec, B in Tesla,A Mass no, z charge
!----------------------------------------------------------------------------
    REAL(DP) zi,Ai,B
      w_ci =  95.5e6_DP*Zi*B/Ai

    RETURN
    END FUNCTION w_ci



    REAL(DP) FUNCTION  w_pi(ni,Ai)
!--------------------------------------------------------------------------
! ion plasma freq. ni in m^-3,A mass no." output in rad/sec
!--------------------------------------------------------------------------
     REAL(DP) Ai,ni
       w_pi = 1.32_DP*SQRT(ni/Ai)

     RETURN
     END FUNCTION w_pi


    REAL(DP) FUNCTION  w_lh(Zi,Ai,B,ne,ni) 
!--------------------------------------------------------------------------
! lower hybrid freq, rad/sec"
!--------------------------------------------------------------------------
     REAL(DP) Ai,ni,Zi,B,ne,wci,wce,wpe,wpi,w1,w2,w3
       wci = w_ci(Zi,Ai,B) ; wce = ABS(w_ce(B))
       wpi = w_pi(ni,Ai)   ; wpe = w_pe(ne)
       w1 = 1._DP/(wci*wci + wpi*wpi)
       w2 = 1._DP/(wce*wci)
       w3 = 1._DP/(w1 + w2 )
       w_lh = SQRT( w3 )
     RETURN

     END FUNCTION w_lh


    REAL(DP) FUNCTION  w_uh(B,ne) 
!--------------------------------------------------------------------------
! upper  hybrid freq, rad/sec"
!--------------------------------------------------------------------------
     REAL(DP) B,ne,wce,wpe
       wce = w_ce(B)
       wpe = w_pe(ne)
       w_uh = SQRT(wce*wce + wpe*wpe)
     RETURN

     END FUNCTION w_uh




    END MODULE pl_freq
