 MODULE dep_var
    USE nrtype, ONLY : DP,I4B
    IMPLICIT NONE
      REAL(DP),DIMENSION(:),ALLOCATABLE,SAVE       :: te,ti,angrot,rbp,ene,             &
                                                      dangrotdt,dnedt,enesav,qhsj,      &
                                                      etor,tot_press,ion_press,         &
                                                      vphi,vpol,vdia,vparm,vexb,        &
                                                      dtedt,dtidt,dpedt,drbpdt,         &
                                                      dangrot_drho,vpol_nclass,         &
                                                      vpar_nclass,er_tot_nclass,        &
                                                      er_vtor_nclass,er_vpol_nclass,    &
                                                      er_dia_nclass,vpinch_d1_nclass,   &
                                                      vpinch_h1_nclass,vpinch_t1_nclass,&
                                                      vpinch_imp_nclass,vpinch_e1_nclass,&
                                                      vpinch_dt1_nclass,vpinch_he4_nclass

      REAL(DP),DIMENSION(:,:),ALLOCATABLE,SAVE     :: en,dnidt,dpidt

      REAL(DP),ALLOCATABLE,DIMENSION(:),SAVE       :: vminus, vplus

      REAL(DP),SAVE                                :: vminuspr,vpluspr

      INTEGER(I4b),PARAMETER                       :: dp4 = 4       !stands for 4 dep var :te ,ti,rbp,angrot


      INTEGER(I4B),DIMENSION(:),ALLOCATABLE,SAVE ::   itran_local

      INTEGER(I4B),DIMENSION(:),ALLOCATABLE,SAVE ::   itran_save,itenp_save,iteni_save

      REAL(DP), DIMENSION(:), ALLOCATABLE,  SAVE ::   change,xval_global

      REAL(DP), DIMENSION(:), ALLOCATABLE,  SAVE ::   xscale_global,typf_global

      REAL(DP), DIMENSION(:), ALLOCATABLE,  SAVE ::   typx_global

      REAL(DP), DIMENSION(:), ALLOCATABLE,  SAVE ::   xval,xval_in

      REAL(DP), DIMENSION(:), ALLOCATABLE,  SAVE ::   xval_old
 

 CONTAINS
   SUBROUTINE allocate_dep_var(nj,nion,ntot)
         USE nrtype,                                  ONLY : I4B,DP
         USE common_constants,                        ONLY : zeroc,izero

         IMPLICIT NONE
         INTEGER(I4B) nj,nion,ntot


         IF(.NOT. ALLOCATED(rbp))         ALLOCATE(rbp(nj))
         rbp(:)      =    zeroc
         IF(.NOT. ALLOCATED(te))          ALLOCATE(te(nj))
         te(:)      =    zeroc
         IF(.NOT. ALLOCATED(ti))          ALLOCATE(ti(nj))
         ti(:)      =    zeroc
         IF(.NOT. ALLOCATED(ene))         ALLOCATE(ene(nj))
         ene(:)    =    zeroc
         IF(.NOT. ALLOCATED(angrot))      ALLOCATE(angrot(nj))
         angrot(:)  =    zeroc
         IF(.NOT. ALLOCATED(dangrot_drho))ALLOCATE(dangrot_drho(1:nj-1))
         dangrot_drho(:)  =    zeroc
         IF(.NOT. ALLOCATED(en))          ALLOCATE(en(nj,nion))
         en(:,:)    =    zeroc
         IF(.NOT. ALLOCATED(dpidt))       ALLOCATE(dpidt(nj,nion))
         dpidt(:,:)    =    zeroc
         IF(.NOT. ALLOCATED(dnidt))       ALLOCATE(dnidt(nj,nion))
         dnidt(:,:)    =    zeroc
         IF(.NOT. ALLOCATED(enesav))      ALLOCATE(enesav(nj))
         enesav(:)  =    zeroc
         IF(.NOT. ALLOCATED(vminus))      ALLOCATE(vminus(ntot),vplus(ntot))
         vminus =zeroc ; vplus = zeroc
         IF( .NOT. ALLOCATED(etor))       ALLOCATE(etor(nj))
         etor(1:nj) = zeroc
         IF(.NOT. ALLOCATED(tot_press))   ALLOCATE(tot_press(nj))
         tot_press(:)    =    zeroc
         IF(.NOT. ALLOCATED(ion_press))   ALLOCATE(ion_press(nj))
         ion_press(:)    =    zeroc
         IF( .NOT. ALLOCATED(vphi))       ALLOCATE(vphi(nj))
         vphi(1:nj) = zeroc
         IF( .NOT. ALLOCATED(vdia))       ALLOCATE(vdia(nj))
         vdia(1:nj) = zeroc
         IF( .NOT. ALLOCATED(vpol))       ALLOCATE(vpol(nj))
         vpol(1:nj) = zeroc


         IF( .NOT. ALLOCATED(vpol_nclass))ALLOCATE(vpol_nclass(nj))
         vpol_nclass(1:nj) = zeroc
         IF( .NOT. ALLOCATED(vpar_nclass))ALLOCATE(vpar_nclass(nj))
         vpar_nclass(1:nj) = zeroc
         IF( .NOT. ALLOCATED(er_tot_nclass))ALLOCATE(er_tot_nclass(nj))
         er_tot_nclass(1:nj) = zeroc
         IF( .NOT. ALLOCATED(er_vtor_nclass))ALLOCATE(er_vtor_nclass(nj))
         er_vtor_nclass(1:nj) = zeroc
         IF( .NOT. ALLOCATED(er_vpol_nclass))ALLOCATE(er_vpol_nclass(nj))
         er_vpol_nclass(1:nj) = zeroc
         IF( .NOT. ALLOCATED(er_dia_nclass))ALLOCATE(er_dia_nclass(nj))
         er_dia_nclass(1:nj) = zeroc
         IF( .NOT. ALLOCATED(vpinch_e1_nclass))ALLOCATE(vpinch_e1_nclass(nj))
         vpinch_e1_nclass(1:nj) = zeroc
         IF( .NOT. ALLOCATED(vpinch_d1_nclass))ALLOCATE(vpinch_d1_nclass(nj))
         vpinch_d1_nclass(1:nj) = zeroc
         IF( .NOT. ALLOCATED(vpinch_h1_nclass))ALLOCATE(vpinch_h1_nclass(nj))
         vpinch_h1_nclass(1:nj) = zeroc
         IF( .NOT. ALLOCATED(vpinch_t1_nclass))ALLOCATE(vpinch_t1_nclass(nj))
         vpinch_t1_nclass(1:nj) = zeroc
         IF( .NOT. ALLOCATED(vpinch_dt1_nclass))ALLOCATE(vpinch_dt1_nclass(nj))
         vpinch_dt1_nclass(1:nj) = zeroc
        IF( .NOT. ALLOCATED(vpinch_he4_nclass))ALLOCATE(vpinch_he4_nclass(nj))
         vpinch_he4_nclass(1:nj) = zeroc
         IF( .NOT. ALLOCATED(vpinch_imp_nclass))ALLOCATE(vpinch_imp_nclass(nj))
         vpinch_imp_nclass(1:nj) = zeroc


         IF( .NOT. ALLOCATED(vparm))      ALLOCATE(vparm(nj))
         vparm(1:nj) = zeroc
         IF( .NOT. ALLOCATED(vexb))       ALLOCATE(vexb(nj))
         vexb(1:nj) = zeroc
         IF( .NOT. ALLOCATED(dtedt))      ALLOCATE(dtedt(nj))
         dtedt(1:nj) = zeroc
         IF( .NOT. ALLOCATED(dtidt))      ALLOCATE(dtidt(nj))
         dtidt(1:nj) = zeroc
         IF( .NOT. ALLOCATED(dpedt))      ALLOCATE(dpedt(nj))
         dpedt(1:nj) = zeroc
         IF( .NOT. ALLOCATED(dnedt))      ALLOCATE(dnedt(nj))
         dnedt(1:nj) = zeroc
         IF( .NOT. ALLOCATED(dangrotdt))  ALLOCATE(dangrotdt(nj))
         dangrotdt(1:nj) = zeroc
         IF( .NOT. ALLOCATED(drbpdt))     ALLOCATE(drbpdt(nj))
         drbpdt(1:nj) = zeroc

         vminuspr = zeroc ; vpluspr = zeroc

        RETURN

   END SUBROUTINE allocate_dep_var

 END MODULE dep_var
