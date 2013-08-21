! -----------------------------------------------SPS 24 Jan 2011-----------------------
MODULE eval_paleo

    USE curden_terms,       ONLY : q,eta
    
    USE dep_var,            ONLY : ene, te, en, ti
    
    USE grid_class,         ONLY : r, hcap, nj
    
    USE ions_gcnmp,         ONLY : zeff
    
    USE paleo,              ONLY : fsa_div_Gamma_pc, eV2J, fsa_div_qe_pc, fsa_div_qi_pc
    
    USE plasma_properties,  ONLY : dischg
    
    USE source_terms_gcnmp, ONLY : spaleo, qepaleo, qipaleo, qbeame
    
    USE solcon_gcnmp,       ONLY : use_paleo
    
    implicit none
     
    
CONTAINS

SUBROUTINE eval_paleo_sources
    USE nrtype,             ONLY : I4B, I2B
    USE MPI_data,           ONLY : myid,master
    IMPLICIT NONE
    INTEGER(I4b) :: i
    INTEGER(I2b) :: iout,get_next_io_unit

    IF(myid == master) iout   = get_next_io_unit ()
    !IF(myid == master) write (*,*) 'Calling eval_paleo_sources'

! Add paleoclassical particle flux as negative source term
    IF (use_paleo(1)==1) THEN
        !write (*,*) 'Setting spaleo'
        spaleo(2:) = fsa_div_gamma_pc( ene(2:), r(2:), dischg%elongxnj%data(2:), hcap(2:), eta(1:) )
    ENDIF
! Add paleoclassical electron energy flux as negative to sources
    IF (use_paleo(2)==1) THEN 
        !write (*,*) 'Setting qepaleo'
        !Convert temperature to eV then flux from W/m^3 to keV/m^3
        qepaleo(2:) = fsa_div_qe_pc(ene(2:), te(2:)*1.D3, zeff(2:), q(2:), &
            & r(2:), dischg%elongxnj%data(2:), dischg%rmajor, hcap(2:), eta(1:) )/(eV2J*1D3)
    ENDIF
! Add paleoclassical ion energy flux as negative to sources
    IF (use_paleo(3)==1) THEN 
        !write (*,*) 'Setting qipaleo'
        !Convert temperature to eV then flux from W/m^3 to keV/m^3
        qipaleo(2:) = fsa_div_qi_pc(en(2:,1), ti(2:)*1.D3, r(2:), &
            & dischg%elongxnj%data(2:), hcap(2:), eta(1:) )/(eV2J*1D3)
    ENDIF
!    IF(myid == master)THEN
    IF(myid == -1 )THEN ! disable the writes
       OPEN(unit=iout,file='paleo_sources.dat',status='unknown')
       write(iout,'(12a20)') '#r','spaleo','qepaleo','qipaleo', &
            & 'ni','hcap','kappa','eta','ne',&
            & 'te','ti','qbeame'
       write(iout,'(12a20)') '#m','1/m^3/s','keV/m^3/s','keV/m^3/s',&
            & '1/m^3','1','1','ohm-m','1/m^3',&
            & 'keV','keV','?'
       DO i=1,nj
          write(iout,'(12g20.8)') r(i), spaleo(i), qepaleo(i), qipaleo(i),&
               & en(i,1),hcap(i),dischg%elongxnj%data(i),eta(min(i,nj-1)),ene(i),&
               & te(i), ti(i),qbeame(i)
       ENDDO
       close(iout)
    ENDIF
ENDSUBROUTINE   

    
ENDMODULE
