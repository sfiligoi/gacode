   MODULE reqdsk
!------------------------------------------------------------------
! ---  f90 module for eqdsk related items
!------------------------------------------------------------------
  USE nrtype,    ONLY: DP,I4B
 

  IMPLICIT NONE
  INTEGER,PARAMETER               :: nbr_profiles = 15
  TYPE o12_data
      INTEGER(I4B) nj, nprim, nimp, nti,nion,npsi,kj,kk,ke,kb,     &
                   p12_start,p12_end
      REAL(DP) time0,time,dt,realn,eqtim0,dtt,realibeam
      REAL(DP),DIMENSION(:),POINTER   ::  atw,zeff,r,psir,curden,  &
                                          te,ti,ene,spbb,pressb,   &
                                          fcap,gcap,hcap,r2cap,bp, &
                                          pprim,ffprim,rho,        &
                                          q,r2capi,rcap,enasav,    &
                                          wasav  => NULL()
      REAL(DP),DIMENSION(:,:),POINTER :: z,en,u => NULL()
      REAL(DP),DIMENSION(:,:,:),POINTER ::enbsav,wbsav => NULL()
  END TYPE o12_data


  TYPE ptr_array
     REAL(DP),DIMENSION(:),  POINTER :: ptr_1d_profile => NULL()
  END TYPE ptr_array

  TYPE eqdsk
     INTEGER(I4b)         nreqd,nzeqd,ncontr,ipestg,nlimiter,       &
                          nptr,peq_start,peq_end,ffpindex,ppindex,  &
                          ffp_nsfit,pp_nsfit
     CHARACTER(len =8)    ntitle(5),dat
     CHARACTER(len=13)    eqdsrce 
     CHARACTER(LEN=256)   name
     LOGICAL eqtype
     REAL(DP)              rdimeqd, zdimeqd, reqd, redeqd, zmideqd, &
                           rma,zma,psimag,psilim,beqd,toteqd,       &
                           psimx1, psimx2, rax1, rax2 ,zax1,zax2,   &
                           psisep,rsep,zsep,eqdskvol,eqdskarea,     &
                           ffpweight,ppweight

     REAL(DP),DIMENSION(:),POINTER   :: fpsi,presspsi,ffpsi,ppsi,   &
                                        qpsi,rcontr,zcontr,rlimiter,&
                                        zlimiter,rmhdgrid,zmhdgrid, &
                                        psival,norm_psi_grid,       &
                                        smth_ffpsi,smth_ppsi,       &
                                        smth_ffpsi_fit,             &
                                        smth_ppsi_fit,ffpsi_knots,  &
                                        ppsi_knots,                 &
                                        ffpsi_spline_eval,          &
                                        ppsi_spline_eval => NULL()
    
     REAL(DP),DIMENSION(:,:),POINTER :: psi => NULL()
     REAL(DP),DIMENSION(:,:),POINTER :: pressure => NULL()
     TYPE(o12_data) eq12
     !to define an array of pointers as opposed to pointer arrays,
     !which are different,we have to do the following:
     TYPE(ptr_array),DIMENSION(nbr_profiles)    :: profile_ptr_1d
   END TYPE eqdsk

   TYPE  eqdsk_ptr_array
      TYPE(eqdsk), POINTER :: eqdsk_ptr => NULL()
   END TYPE eqdsk_ptr_array
     
   TYPE eq_list
      type(eqdsk_ptr_array),DIMENSION(2) :: eq_list_ptr 
   END TYPE eq_list

     CONTAINS 

        SUBROUTINE read_eqdsk(eqdsk_in,nout,ncrt,nin,eqdskname)
!
!
! ----------------------------------------------------------------------
! this subroutine reads an eqdsk as a starting guess for an equilibrium.
! (quantities indexed by nreqd are read from the eqdsk in reverse
! order so that for example qpsi(1) = edge,qpsi(nreqd) = axis value.)
!  presspsi(1..nreqd)    pressure (nt/m**2)
!  fpsi(1..nreqd)        f(psi),tesla/m
!  ffpsi(1..nreqd)      f(psi)*d(f)/dpsi
!  pppsi(1..nreqd)       d(presspsi)/dpsi
!  qpsi(1..nreqd)        safety factor
!  pressb(1...nreqd)     beam pressure (if known)
!
!  on exit from this subroutine all quantities are in guassian units!
!  the switch between gaussian and MKS (and MKS back to gaussian) is
!  done in subroutine MHD.
!  The actual MHD calculations are done in MKS units.
! ------------------------------------------------------------------ HSJ
!

      USE nrtype,              ONLY : DP,I4B,I2B


       IMPLICIT NONE 

       TYPE(eqdsk)  eqdsk_in

       INTEGER(I4b) limread,nout,ncrt,contrpts_set,nin,size,i,j,         &
                    k,uniform_grid
       CHARACTER(len=*) eqdskname


!
      eqdsk_in%name(:) =' '
      eqdsk_in%name(1:LEN_TRIM(eqdskname))  = eqdskname(1:LEN_TRIM(eqdskname))
      eqdsk_in%nptr =0
      eqdsk_in%ncontr  = 0 
      READ (nin, 8190) (eqdsk_in%ntitle(i), i=1,5), &
           eqdsk_in%dat, eqdsk_in%ipestg, eqdsk_in%nreqd,       &
           eqdsk_in%nzeqd,eqdsk_in%eqdsrce                       !line1

      CALL associate_eqdsk(eqdsk_in)



      READ (nin, 8200)  eqdsk_in%rdimeqd, eqdsk_in%zdimeqd,            &
                   eqdsk_in%reqd, eqdsk_in%redeqd, eqdsk_in%zmideqd 
      IF (ABS (eqdsk_in%zmideqd) .GT. 1.0e-5) THEN
          WRITE  (nout, 8251)  eqdsk_in%zmideqd
          WRITE  (ncrt, 8251)  eqdsk_in%zmideqd
 8251     FORMAT (' subroutine REQDSK detects an error:' / &
                  ' zmideqd =', 1pe12.6                  / &
                  ' we expect zmideqd = 0.0'        / &
                  ' code  must stop'                     )
          CALL STOP ('subroutine REQDSK: problem #2', 32)
      END IF
!
      eqdsk_in%eqtype = .FALSE.
      IF (eqdsk_in%eqdsrce .EQ.  'ONETWO EQDSK '  .OR. &
          eqdsk_in%eqdsrce .EQ. 'FIXBDRY EQDSK')  eqdsk_in%eqtype = .TRUE.
!
      contrpts_set = 0    ! set to 1 below if boundary pts are available
      READ (nin, 8200)  eqdsk_in%rma, eqdsk_in%zma, eqdsk_in%psimag,  &
                            eqdsk_in%psilim, eqdsk_in%beqd               !line3
      READ (nin, 8200)  eqdsk_in%toteqd, eqdsk_in%psimx1, eqdsk_in%psimx2,&
                            eqdsk_in%rax1, eqdsk_in%rax2                  !line4
      READ (nin, 8200)  eqdsk_in%zax1, eqdsk_in%zax2, eqdsk_in%psisep,     &
                        eqdsk_in%rsep, eqdsk_in%zsep
      READ (nin, 8200)  ( eqdsk_in%fpsi(i), i=eqdsk_in%nreqd,1,-1)
!                                               fpsi(nreqd) = axis value


      READ (nin, 8200)  (eqdsk_in%presspsi(i), i=eqdsk_in%nreqd,1,-1)
!                                           presspsi(nreqd) = axis value
      READ (nin, 8200)  (eqdsk_in%ffpsi(i), i=eqdsk_in%nreqd,1,-1)
      READ (nin, 8200)  (eqdsk_in%ppsi(i), i=eqdsk_in%nreqd,1,-1)
      READ (nin, 8200)  ((eqdsk_in%psi(i,j), i=1,eqdsk_in%nreqd), j=1,eqdsk_in%nzeqd)
      READ (nin, 8200)  (eqdsk_in%qpsi(i), i=eqdsk_in%nreqd,1,-1)
      READ (nin, 8210, err=4000, END=4000)   eqdsk_in%ncontr, eqdsk_in%nlimiter
      

      CALL associate_eqdsk_1(eqdsk_in)

      READ (nin, 8200, err=4000, END=4000)  (eqdsk_in%rcontr(i),       &
                                        eqdsk_in%zcontr(i),i = 1,eqdsk_in%ncontr)


      IF (eqdsk_in%nlimiter .GT. 0) THEN
          CALL volcalc (eqdsk_in%rcontr,eqdsk_in%zcontr,eqdsk_in%ncontr,     &
                        eqdsk_in%rma,eqdsk_in%zma,eqdsk_in%eqdskvol,         &
                        eqdsk_in%eqdskarea)
          WRITE (ncrt, '(1x, a, f12.6)') &
            'volume determined from eqdsk boundary values = ', eqdsk_in%eqdskvol
          WRITE (ncrt, '(1x, a, f12.6)') &
            'area   determined from eqdsk boundary values = ', eqdsk_in%eqdskarea
      END IF

 

!
! --- pick up the limiter
!

        limread  = 1 
        READ (nin, 8200, err = 8220, END = 8230) &
             (eqdsk_in%rlimiter(i), eqdsk_in%zlimiter(i), i=1,eqdsk_in%nlimiter)
        print *,'eqdsk ',eqdsk_in%name,' has ',eqdsk_in%nlimiter,' limiter points'

        go to 4000
        size=LEN_TRIM(eqdsk_in%name)
 8220   WRITE  (nout, 8240)  eqdsk_in%name(1:size)
        WRITE  (ncrt, 8240)  eqdsk_in%name(1:size)
 8240   FORMAT (/ ' ERROR: eqdsk file "', a, '"'                    / &
                      8x, 'has generated an error in while reading' / &
                      8x, 'the limiter points in subroutine REQDSK' /)
        eqdsk_in%nlimiter = -1
        go to 4000
 8230   WRITE  (nout, 8250)  eqdsk_in%name(1:size)
        WRITE  (ncrt, 8250)  eqdsk_in%name(1:size)
 8250   FORMAT (/ ' ERROR: eqdsk file "', a, '"'            / &
                      8x, 'does not contain limiter points' /)
        eqdsk_in%nlimiter = -1
    
     
!
 8190 FORMAT (6a8, 3i4, t73, a)
 8200 FORMAT (5e16.9)
 8210 FORMAT (2i5)
 8235 FORMAT (4(2x, i5))
!
! ----------------------------------------------------------------------
! --- the following output is not part of the "official" eqdsk. It is
! --- what distinguishes eqdsk created in ONETWO from those created by EFIT:
! --- we put out enough information so that the eqdsk could be used as
! --- a restart file.
! ----------------------------------------------------------------------
!
 

 4000 CONTINUE
       
      eqdsk_in%peq_end   = eqdsk_in%nptr !may be redefined in associate_eqdsk_eq12
      IF (eqdsk_in%eqdsrce .EQ. 'ONETWO EQDSK') THEN
 
           READ (nin, 8235)   eqdsk_in%eq12%nj, eqdsk_in%eq12%nprim,         &
              eqdsk_in%eq12%nimp, eqdsk_in%eq12%nti,eqdsk_in%eq12%npsi
           eqdsk_in%eq12%nion = eqdsk_in%eq12%nprim + eqdsk_in%eq12%nimp
      
            call associate_eqdsk_eq12(eqdsk_in)
            READ (nin, 8200)  (eqdsk_in%eq12%atw(k), k=1,eqdsk_in%eq12%nion)
            DO  k=1,eqdsk_in%eq12%nion
               READ (nin, 8200)  (eqdsk_in%eq12%z     (j,k), j=1,eqdsk_in%eq12%nj)
            ENDDO
            READ (nin, 8200)  (eqdsk_in%eq12%zeff  (j  ), j=1,eqdsk_in%eq12%nj)
            READ (nin, 8200)  (eqdsk_in%eq12%r     (j  ), j=1,eqdsk_in%eq12%nj)
            READ (nin, 8200)  (eqdsk_in%eq12%psir  (j  ), j=1,eqdsk_in%eq12%nj)
            READ (nin, 8200)  (eqdsk_in%eq12%curden(j  ), j=1,eqdsk_in%eq12%nj)
            READ (nin, 8200)  (eqdsk_in%eq12%te    (j  ), j=1,eqdsk_in%eq12%nj)

            DO  k=1,eqdsk_in%eq12%nti
               READ (nin, 8200)  (eqdsk_in%eq12%ti    (j  ), j=1,eqdsk_in%eq12%nj)
            ENDDo
            READ (nin, 8200)  (eqdsk_in%eq12%ene   (j  ), j=1,eqdsk_in%eq12%nj)
            DO  k=1,eqdsk_in%eq12%nprim
                  READ (nin, 8200)  (eqdsk_in%eq12%en    (j,k), j=1,eqdsk_in%eq12%nj)
            ENDDO

            IF (eqdsk_in%eq12%nimp .NE. 0)   THEN
                     DO 740 k=eqdsk_in%eq12%nprim+1,eqdsk_in%eq12%nion
740                     READ (nin,8200) (eqdsk_in%eq12%en    (j,k), j=1,eqdsk_in%eq12%nj)
            END IF
                     READ (nin, 8200) (eqdsk_in%eq12%spbb(i), i=eqdsk_in%nreqd,1,-1)
                     !
                     READ (nin, 8200) (eqdsk_in%eq12%pressb(i), i=1,eqdsk_in%eq12%nj)
print *,'line 245'
                     READ (nin, 8200) ((eqdsk_in%eq12%u(i,j), i=1,eqdsk_in%eq12%nion+4), j=1,eqdsk_in%eq12%nj)
print *,'line 247'
                     READ (nin, 8200) eqdsk_in%eq12%time0,eqdsk_in%eq12%time,&
                          eqdsk_in%eq12%dt,eqdsk_in%eq12%realn,  &
                          eqdsk_in%eq12%eqtim0,eqdsk_in%eq12%dtt,&
                          eqdsk_in%eq12%realibeam
 print *,'line 251'
                     READ (nin, 8200) (eqdsk_in%eq12%fcap(j),j=1,eqdsk_in%eq12%nj)
                     READ (nin, 8200) (eqdsk_in%eq12%gcap(j),j=1,eqdsk_in%eq12%nj)
                     READ (nin, 8200) (eqdsk_in%eq12%hcap(j),j=1,eqdsk_in%eq12%nj)
                     READ (nin, 8200) (eqdsk_in%eq12%r2cap(j),j=1,eqdsk_in%eq12%nj)

                     READ (nin, 8200) (eqdsk_in%eq12%bp(j),j=1,eqdsk_in%eq12%nj)
                     READ (nin, 8200) (eqdsk_in%eq12%pprim(j),j=1,eqdsk_in%eq12%nj)
                     READ (nin, 8200) (eqdsk_in%eq12%ffprim(j),j=1,eqdsk_in%eq12%nj)
!                     READ (nin, 8200) (eqdsk_in%eq12%psival(j),j=1,eqdsk_in%eq12%npsi)
 
                     READ (nin, 8200) (eqdsk_in%eq12%rho(j),j=1,eqdsk_in%eq12%npsi)
                     READ (nin, 8200) (eqdsk_in%eq12%q(j),j=1,eqdsk_in%eq12%nj)
                     READ (nin, 8200) (eqdsk_in%eq12%r2capi(j),j=1,eqdsk_in%eq12%nj)
                     READ (nin, 8200) (eqdsk_in%eq12%rcap(j),j=1,eqdsk_in%eq12%nj)
                     ! rcapi is not currently included

                     READ (nin, 8200) (((eqdsk_in%eq12%enbsav(i,j,k),i=1,eqdsk_in%eq12%nj),               &
                                                           j=1,eqdsk_in%eq12%ke),k=1,eqdsk_in%eq12%kb)
                     READ (nin, 8200) (((eqdsk_in%eq12%wbsav(i,j,k),i=1,eqdsk_in%eq12%nj),                &
                                                           j=1,eqdsk_in%eq12%ke),k=1,eqdsk_in%eq12%kb)
                     READ (nin, 8200) (eqdsk_in%eq12%enasav(i),i=1,eqdsk_in%eq12%nj)
                     READ (nin, 8200) (eqdsk_in%eq12%wasav(i),i=1,eqdsk_in%eq12%nj)

           END IF



!
! --- set rmhdgrid and zmhdgrid for eqdsk:
      CALL set_rz_grid(eqdsk_in)

!
! --- set up the psival grid: psival(1) = edge, psival(nreqd) = axis
!
      uniform_grid = 1
      CALL psiset (eqdsk_in,uniform_grid)

!
! ----------------------------------------------------------------------
! close eqdsk  file
! ----------------------------------------------------------------------
!
      CLOSE (unit = nin)


      RETURN



!
! --- fatal errors, stop code
!
 5000 WRITE  (ncrt, 5100)  eqdsk_in%eq12%npsi,eqdsk_in%nreqd,eqdsk_in%ncontr,    &
                                                         eqdsk_in%nlimiter
      WRITE  (nout, 5100)  eqdsk_in%eq12%npsi,eqdsk_in%nreqd,eqdsk_in%ncontr,    &
                                                         eqdsk_in%nlimiter
 5100 FORMAT (/ ' ERROR detected in subroutine REQDSK'          / &
                '   parameter settings inconsistent with eqdsk' / &
                '   kpsi  = ', 3(2x, i6)     / &
                '   nreqd, ncontr , nlimiter   = ', 3(2x, i6)     / &
                '   ONETWO cannot continue')
      CALL STOP ('subroutine REQDSK: eqdsk consistency problem', 33)
!

!
      END SUBROUTINE read_eqdsk



      SUBROUTINE associate_eqdsk(eqdsk_in)
        USE nrtype, ONLY                          : I4B,DP
        IMPLICIT NONE
        TYPE(eqdsk) eqdsk_in


      eqdsk_in%nptr =0
      eqdsk_in%ncontr  = 0 

      IF(.NOT. ASSOCIATED(eqdsk_in%fpsi))                               &
                      ALLOCATE(eqdsk_in%fpsi(eqdsk_in%nreqd))
      eqdsk_in%nptr = eqdsk_in%nptr +1  ; eqdsk_in%peq_start = eqdsk_in%nptr 
      eqdsk_in%profile_ptr_1d(eqdsk_in%nptr)%ptr_1d_profile => eqdsk_in%fpsi

      IF(.NOT. ASSOCIATED(eqdsk_in%presspsi))                           &
                      ALLOCATE(eqdsk_in%presspsi(eqdsk_in%nreqd))
      eqdsk_in%nptr = eqdsk_in%nptr +1 
      eqdsk_in%profile_ptr_1d(eqdsk_in%nptr)%ptr_1d_profile => eqdsk_in%presspsi

      IF(.NOT. ASSOCIATED(eqdsk_in%ffpsi))                              &
                      ALLOCATE(eqdsk_in%ffpsi(eqdsk_in%nreqd))
      eqdsk_in%nptr = eqdsk_in%nptr +1
      eqdsk_in%profile_ptr_1d(eqdsk_in%nptr)%ptr_1d_profile => eqdsk_in%ffpsi
      eqdsk_in%ffpindex =  eqdsk_in%nptr 

      IF(.NOT. ASSOCIATED(eqdsk_in%ppsi))                               &
                      ALLOCATE(eqdsk_in%ppsi(eqdsk_in%nreqd)) 
      eqdsk_in%nptr = eqdsk_in%nptr +1
      eqdsk_in%profile_ptr_1d(eqdsk_in%nptr)%ptr_1d_profile => eqdsk_in%ppsi
      eqdsk_in%ppindex = eqdsk_in%nptr


      IF(.NOT. ASSOCIATED(eqdsk_in%qpsi))                               &
                      ALLOCATE(eqdsk_in%qpsi(eqdsk_in%nreqd)) 
      eqdsk_in%nptr = eqdsk_in%nptr +1
      eqdsk_in%profile_ptr_1d(eqdsk_in%nptr)%ptr_1d_profile => eqdsk_in%qpsi

      IF(.NOT. ASSOCIATED(eqdsk_in%rmhdgrid))                                &
                      ALLOCATE(eqdsk_in%rmhdgrid(eqdsk_in%nreqd))

      IF(.NOT. ASSOCIATED(eqdsk_in%zmhdgrid))                                &
                      ALLOCATE(eqdsk_in%zmhdgrid(eqdsk_in%nzeqd))


      IF(.NOT. ASSOCIATED(eqdsk_in%psi))                                &
                      ALLOCATE(eqdsk_in%psi(eqdsk_in%nreqd,eqdsk_in%nzeqd)) 

      IF(.NOT. ASSOCIATED(eqdsk_in%psival))               &
                    ALLOCATE(eqdsk_in%psival(eqdsk_in%nreqd))


      END SUBROUTINE associate_eqdsk


      SUBROUTINE associate_eqdsk_1(eqdsk_in)
        USE nrtype, ONLY                          : I4B,DP
        IMPLICIT NONE
        TYPE(eqdsk) eqdsk_in

        IF(.NOT. ASSOCIATED(eqdsk_in%rcontr))                                &
                      ALLOCATE(eqdsk_in%rcontr(eqdsk_in%ncontr)) 

        IF(.NOT. ASSOCIATED(eqdsk_in%zcontr))                                &
                      ALLOCATE(eqdsk_in%zcontr(eqdsk_in%ncontr)) 

        IF(.NOT. ASSOCIATED(eqdsk_in%rlimiter))                           &
                        ALLOCATE(eqdsk_in%rlimiter(eqdsk_in%nlimiter))
        IF(.NOT. ASSOCIATED(eqdsk_in%zlimiter))                           &
                        ALLOCATE(eqdsk_in%zlimiter(eqdsk_in%nlimiter))
      END SUBROUTINE associate_eqdsk_1



      SUBROUTINE associate_eqdsk_eq12(eqdsk_in)



        USE nrtype, ONLY                          : I4B,DP
        IMPLICIT NONE
        TYPE(eqdsk) eqdsk_in




           !we have to assume that kpsi = nreqd  since kpsi is not in the eqdsk file:
           !we have to assume that kj =nj since kj is not in the eqdsk file:
           !we have to assume that kb = 8 since kb is not in the eqdsk file:
           !we have to assume that kk = 9 since kk is not in the eqdsk file:
           !we have to assume that ke = 3 since ke is not in the eqdsk file:





           eqdsk_in%peq_end   = eqdsk_in%nptr
!           eqdsk_in%eq12%kpsi =  eqdsk_in%nreqd 
!           eqdsk_in%nptr = eqdsk_in%nptr +1
           eqdsk_in%eq12%p12_start = eqdsk_in%nptr+1


 
           eqdsk_in%eq12%kj  = eqdsk_in%eq12%nj
           eqdsk_in%eq12%kj = 201

           eqdsk_in%eq12%kk = 9

           eqdsk_in%eq12%kb = 8

           eqdsk_in%eq12%ke = 3

           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%atw))               &
                    ALLOCATE(eqdsk_in%eq12%atw(eqdsk_in%eq12%nion))

           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%zeff))               &
                    ALLOCATE(eqdsk_in%eq12%zeff(eqdsk_in%eq12%nj))
           eqdsk_in%nptr = eqdsk_in%nptr +1
           eqdsk_in%profile_ptr_1d(eqdsk_in%nptr)%ptr_1d_profile => eqdsk_in%eq12%zeff
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%r))               &
                    ALLOCATE(eqdsk_in%eq12%r(eqdsk_in%eq12%nj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%psir))               &
                    ALLOCATE(eqdsk_in%eq12%psir(eqdsk_in%eq12%nj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%curden))               &
                    ALLOCATE(eqdsk_in%eq12%curden(eqdsk_in%eq12%nj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%te))               &
                    ALLOCATE(eqdsk_in%eq12%te(eqdsk_in%eq12%nj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%ti))               &
                    ALLOCATE(eqdsk_in%eq12%ti(eqdsk_in%eq12%nj))

           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%ene))               &
                    ALLOCATE(eqdsk_in%eq12%ene(eqdsk_in%eq12%nj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%spbb))               &
                    ALLOCATE(eqdsk_in%eq12%spbb(eqdsk_in%nreqd))

           eqdsk_in%nptr = eqdsk_in%nptr +1
           eqdsk_in%profile_ptr_1d(eqdsk_in%nptr)%ptr_1d_profile => eqdsk_in%eq12%spbb

           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%pressb))               &
                    ALLOCATE(eqdsk_in%eq12%pressb(eqdsk_in%eq12%nj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%fcap))               &
                    ALLOCATE(eqdsk_in%eq12%fcap(eqdsk_in%eq12%kj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%gcap))               &
                    ALLOCATE(eqdsk_in%eq12%gcap(eqdsk_in%eq12%kj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%hcap))               &
                    ALLOCATE(eqdsk_in%eq12%hcap(eqdsk_in%eq12%kj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%rcap))               &
                    ALLOCATE(eqdsk_in%eq12%rcap(eqdsk_in%eq12%kj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%r2cap))               &
                    ALLOCATE(eqdsk_in%eq12%r2cap(eqdsk_in%eq12%kj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%r2capi))               &
                    ALLOCATE(eqdsk_in%eq12%r2capi(eqdsk_in%eq12%kj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%bp))               &
                    ALLOCATE(eqdsk_in%eq12%bp(eqdsk_in%eq12%kj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%pprim))               &
                    ALLOCATE(eqdsk_in%eq12%pprim(eqdsk_in%eq12%kj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%ffprim))               &
                    ALLOCATE(eqdsk_in%eq12%ffprim(eqdsk_in%eq12%kj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%rho))               &
                    ALLOCATE(eqdsk_in%eq12%rho(eqdsk_in%eq12%npsi))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%q))               &
                    ALLOCATE(eqdsk_in%eq12%q(eqdsk_in%eq12%npsi))

           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%z))                 &
                    ALLOCATE(eqdsk_in%eq12%z(eqdsk_in%eq12%nj,eqdsk_in%eq12%nion))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%en))                &
                    ALLOCATE(eqdsk_in%eq12%en(eqdsk_in%eq12%nj,eqdsk_in%eq12%nion))
 
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%u))                 &
                    ALLOCATE(eqdsk_in%eq12%u(eqdsk_in%eq12%nion+4,eqdsk_in%eq12%nj))

           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%enbsav))                 &
                    ALLOCATE(eqdsk_in%eq12%enbsav(eqdsk_in%eq12%nj,eqdsk_in%eq12%ke,eqdsk_in%eq12%kb))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%wbsav))                 &
                    ALLOCATE(eqdsk_in%eq12%wbsav(eqdsk_in%eq12%nj,eqdsk_in%eq12%ke,eqdsk_in%eq12%kb))

           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%enasav))               &
                    ALLOCATE(eqdsk_in%eq12%enasav(eqdsk_in%eq12%kj))
           IF(.NOT. ASSOCIATED(eqdsk_in%eq12%wasav))               &
                    ALLOCATE(eqdsk_in%eq12%wasav(eqdsk_in%eq12%kj))


           eqdsk_in%eq12%p12_end = eqdsk_in%nptr

      END SUBROUTINE associate_eqdsk_eq12





      SUBROUTINE define_eqdsk(eqdsk_in,eqdsk_out)
        USE nrtype, ONLY                          : I4B,DP
        IMPLICIT NONE
        INTEGER(i4B) leq
        TYPE(eqdsk) eqdsk_in,eqdsk_out

!       define some quantities for eqdsk_out, based on values from eqdsk_in:
        eqdsk_out%name(:) = ' '
        leq =LEN_TRIM(ADJUSTL(eqdsk_in%name))
        eqdsk_out%name(1:leq) = eqdsk_in%name(1:leq)
        eqdsk_out%name(leq+1:leq+13) = '_interpolated'
        eqdsk_out%ntitle(:) = eqdsk_in%ntitle(:)
        eqdsk_out%dat      = eqdsk_in%dat  
        eqdsk_out%ipestg   = eqdsk_in%ipestg
        eqdsk_out%eqdsrce  = eqdsk_in%eqdsrce    
        eqdsk_out%rdimeqd  = eqdsk_in%rdimeqd
        eqdsk_out%zdimeqd  = eqdsk_in%zdimeqd
        eqdsk_out%redeqd   = eqdsk_in%redeqd
        eqdsk_out%psimag   = eqdsk_in%psimag
        eqdsk_out%psilim   = eqdsk_in%psilim
        eqdsk_out%reqd     = eqdsk_in%reqd
        eqdsk_out%zmideqd  = eqdsk_in%zmideqd
        eqdsk_out%rma      = eqdsk_in%rma
        eqdsk_out%zma      = eqdsk_in%zma
        eqdsk_out%beqd     = eqdsk_in%beqd
        eqdsk_out%toteqd   = eqdsk_in%toteqd 
        eqdsk_out%psimx1   = eqdsk_in%psimx1    
        eqdsk_out%psimx2   = eqdsk_in%psimx2
        eqdsk_out%rax1     = eqdsk_in%rax1   
        eqdsk_out%rax2     = eqdsk_in%rax2
        eqdsk_out%zax1     = eqdsk_in%zax1 
        eqdsk_out%zax2     = eqdsk_in%zax2
        eqdsk_out%psisep   = eqdsk_in%psisep
        eqdsk_out%rsep     = eqdsk_in%rsep
        eqdsk_out%zsep     = eqdsk_in%zsep
 


        CALL associate_eqdsk(eqdsk_out)
        CALL associate_eqdsk_1(eqdsk_out)
        eqdsk_out%ncontr      =  eqdsk_in%ncontr
        eqdsk_out%nlimiter    =  eqdsk_in%nlimiter
        eqdsk_out%rlimiter(:) = eqdsk_in%rlimiter(:) 
        eqdsk_out%zlimiter(:) = eqdsk_in%zlimiter(:)
        eqdsk_out%rcontr(:)   = eqdsk_in%rcontr
        eqdsk_out%zcontr(:)   = eqdsk_in%zcontr
 
        CALL set_rz_grid(eqdsk_out)
        CALL associate_eqdsk_eq12(eqdsk_out)


      END SUBROUTINE define_eqdsk




       SUBROUTINE write_eqdsk(eqdsk_in,nout,ncrt)
!
!
! ----------------------------------------------------------------------
! this subroutine writes an eqdsk .
! 
! ------------------------------------------------------------------ HSJ
!

      USE nrtype,              ONLY : DP,I4B,I2B


       IMPLICIT NONE 

       TYPE(eqdsk)  eqdsk_in

       INTEGER(I4b) limread,nout,iostat,ncrt,i,j,k
       
       OPEN (unit = nout, file = eqdsk_in%name , status = 'UNKNOWN', iostat = iostat)
       PRINT *,'Creating eqdsk : ', eqdsk_in%name(1:LEN_TRIM( eqdsk_in%name))
       PRINT *,' with nr,nz =',eqdsk_in%nreqd,eqdsk_in%nzeqd
      WRITE(nout, 8190) (eqdsk_in%ntitle(i), i=1,5), &
           eqdsk_in%dat, eqdsk_in%ipestg, eqdsk_in%nreqd,       &
           eqdsk_in%nzeqd,eqdsk_in%eqdsrce(1:LEN_TRIM(eqdsk_in%eqdsrce))
      WRITE (nout, 8200)  eqdsk_in%rdimeqd, eqdsk_in%zdimeqd,            &
                   eqdsk_in%reqd, eqdsk_in%redeqd, eqdsk_in%zmideqd 

      WRITE (nout, 8200)  eqdsk_in%rma, eqdsk_in%zma, eqdsk_in%psimag,  &
                            eqdsk_in%psilim, eqdsk_in%beqd               !line3
      WRITE (nout, 8200)  eqdsk_in%toteqd, eqdsk_in%psimx1, eqdsk_in%psimx2,&
                            eqdsk_in%rax1, eqdsk_in%rax2                  !line4
      WRITE (nout, 8200)  eqdsk_in%zax1, eqdsk_in%zax2, eqdsk_in%psisep,     &
                        eqdsk_in%rsep, eqdsk_in%zsep
      WRITE (nout, 8200)  ( eqdsk_in%fpsi(i), i=eqdsk_in%nreqd,1,-1)
!                                               fpsi(nreqd) = axis value


      WRITE (nout, 8200)  (eqdsk_in%presspsi(i), i=eqdsk_in%nreqd,1,-1)
!                                           presspsi(nreqd) = axis value
      WRITE (nout, 8200)  (eqdsk_in%ffpsi(i), i=eqdsk_in%nreqd,1,-1)
      WRITE (nout, 8200)  (eqdsk_in%ppsi(i), i=eqdsk_in%nreqd,1,-1)
      WRITE (nout, 8200)  ((eqdsk_in%psi(i,j), i=1,eqdsk_in%nreqd), j=1,eqdsk_in%nzeqd)
      WRITE (nout, 8200)  (eqdsk_in%qpsi(i), i=eqdsk_in%nreqd,1,-1)
      WRITE (nout, 8210)   eqdsk_in%ncontr, eqdsk_in%nlimiter
      


      WRITE (nout, 8200)  (eqdsk_in%rcontr(i),       &
                                        eqdsk_in%zcontr(i),i = 1,eqdsk_in%ncontr)
        WRITE (nout, 8200) &
             (eqdsk_in%rlimiter(i), eqdsk_in%zlimiter(i), i=1,eqdsk_in%nlimiter)

!
 8190 FORMAT (6a8, 3i4, t73, a)
 8200 FORMAT (5e16.9)
 8210 FORMAT (2i5)
 8235 FORMAT (4(2x, i5))
!
! ----------------------------------------------------------------------
! --- the following output is not part of the "official" eqdsk. It is
! --- what distinguishes eqdsk created in ONETWO from those created by EFIT:
! --- we put out enough information so that the eqdsk could be used as
! --- a restart file.
! ----------------------------------------------------------------------
!

      IF (eqdsk_in%eqdsrce .EQ. 'ONETWO EQDSK') THEN
           WRITE (nout, 8235)   eqdsk_in%eq12%nj, eqdsk_in%eq12%nprim,         &
              eqdsk_in%eq12%nimp, eqdsk_in%eq12%nti, eqdsk_in%eq12%npsi
            WRITE (nout, 8200)  (eqdsk_in%eq12%atw(k), k=1,eqdsk_in%eq12%nion)
            DO  k=1,eqdsk_in%eq12%nion
               WRITE (nout, 8200)  (eqdsk_in%eq12%z     (j,k), j=1,eqdsk_in%eq12%nj)
            ENDDO
            WRITE (nout, 8200)  (eqdsk_in%eq12%zeff  (j  ), j=1,eqdsk_in%eq12%nj)
            WRITE (nout, 8200)  (eqdsk_in%eq12%r     (j  ), j=1,eqdsk_in%eq12%nj)
            WRITE (nout, 8200)  (eqdsk_in%eq12%psir  (j  ), j=1,eqdsk_in%eq12%nj)
            WRITE (nout, 8200)  (eqdsk_in%eq12%curden(j  ), j=1,eqdsk_in%eq12%nj)
            WRITE (nout, 8200)  (eqdsk_in%eq12%te    (j  ), j=1,eqdsk_in%eq12%nj)
            DO  k=1,eqdsk_in%eq12%nti
               WRITE (nout, 8200)  (eqdsk_in%eq12%ti    (j  ), j=1,eqdsk_in%eq12%nj)
            ENDDo
            WRITE (nout, 8200)  (eqdsk_in%eq12%ene   (j  ), j=1,eqdsk_in%eq12%nj)
            DO  k=1,eqdsk_in%eq12%nprim
                  WRITE (nout, 8200)  (eqdsk_in%eq12%en    (j,k), j=1,eqdsk_in%eq12%nj)
            ENDDO
            IF (eqdsk_in%eq12%nimp .NE. 0)   THEN
                     DO 740 k=eqdsk_in%eq12%nprim+1,eqdsk_in%eq12%nion
740                     WRITE (nout,8200) (eqdsk_in%eq12%en    (j,k), j=1,eqdsk_in%eq12%nj)
            END IF


                     WRITE (nout, 8200) (eqdsk_in%eq12%spbb(i), i=eqdsk_in%nreqd,1,-1)

                     WRITE (nout, 8200) (eqdsk_in%eq12%pressb(i), i=1,eqdsk_in%eq12%kj)
                     WRITE (nout, 8200) ((eqdsk_in%eq12%u(i,j), i=1,eqdsk_in%eq12%kk), j=1,eqdsk_in%eq12%kj)

                     WRITE (nout, 8200) eqdsk_in%eq12%time0,eqdsk_in%eq12%time,&
                          eqdsk_in%eq12%dt,eqdsk_in%eq12%realn,  &
                          eqdsk_in%eq12%eqtim0,eqdsk_in%eq12%dtt,&
                          eqdsk_in%eq12%realibeam

                     WRITE (nout, 8200) (eqdsk_in%eq12%fcap(j),j=1,eqdsk_in%eq12%nj)
                     WRITE (nout, 8200) (eqdsk_in%eq12%gcap(j),j=1,eqdsk_in%eq12%nj)
                     WRITE (nout, 8200) (eqdsk_in%eq12%hcap(j),j=1,eqdsk_in%eq12%nj)
                     WRITE (nout, 8200) (eqdsk_in%eq12%r2cap(j),j=1,eqdsk_in%eq12%nj)
                     WRITE (nout, 8200) (eqdsk_in%eq12%bp(j),j=1,eqdsk_in%eq12%nj)
                     WRITE (nout, 8200) (eqdsk_in%eq12%pprim(j),j=1,eqdsk_in%eq12%nj)
                     WRITE (nout, 8200) (eqdsk_in%eq12%ffprim(j),j=1,eqdsk_in%eq12%nj)

!                     WRITE (nout, 8200) (eqdsk_in%eq12%psival(j),j=1,eqdsk_in%eq12%npsi)
                     WRITE (nout, 8200) (eqdsk_in%eq12%rho(j),j=1,eqdsk_in%eq12%npsi)
                     WRITE (nout, 8200) (eqdsk_in%eq12%q(j),j=1,eqdsk_in%eq12%nj)
                     WRITE (nout, 8200) (eqdsk_in%eq12%r2capi(j),j=1,eqdsk_in%eq12%nj)
                     WRITE (nout, 8200) (eqdsk_in%eq12%rcap(j),j=1,eqdsk_in%eq12%nj)
                     ! rcapi is not currently included
                     WRITE (nout, 8200) (((eqdsk_in%eq12%enbsav(i,j,k),i=1,eqdsk_in%eq12%nj),               &
                                                           j=1,eqdsk_in%eq12%ke),k=1,eqdsk_in%eq12%kb)
                     WRITE (nout, 8200) (((eqdsk_in%eq12%wbsav(i,j,k),i=1,eqdsk_in%eq12%nj),                &
                                                           j=1,eqdsk_in%eq12%ke),k=1,eqdsk_in%eq12%kb)
                     WRITE (nout, 8200) (eqdsk_in%eq12%enasav(i),i=1,eqdsk_in%eq12%nj)
                     WRITE (nout, 8200) (eqdsk_in%eq12%wasav(i),i=1,eqdsk_in%eq12%nj)



           END IF






      RETURN




!
      END SUBROUTINE write_eqdsk

 END MODULE reqdsk
