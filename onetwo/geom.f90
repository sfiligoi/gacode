!
      MODULE geom
      USE param,only: kj
      integer       comp_methd_eqdsk, itorfluse
      real *8                                                         &
            timcap, rhoa0, drhoadt_geom, rhoa, fcap0(kj),             &
            dfdt(kj), rhoa_save,                                      &
            fcap(kj), gcap0(kj), dgdt(kj), gcap(kj), hcap0(kj),       &
            dhdt(kj), hcap(kj), hcapctr(kj), fcapctr(kj),             &
            dhdtsv(kj), dfdtsv(kj), dhdtcap(kj), dfdtcap(kj),         &
            r2cap0(kj), dr2dt(kj), r2cap(kj), rdif, ddifdt,           &
            rhosp0, xi(kj), sfarea, cxarea, volfac, r2capi(kj),       &
            r2capi0(kj), dr2idt(kj), rcap(kj),bsqncap(kj),rcapi(kj),  &
            ali, rcap0(kj),bsq_avg_cap(kj),kappa_nj(kj),              &
            drcapdt(kj), rcap0i(kj),drcapidt(kj),b_avg_cap(kj),       &
            h_factr_cap(kj), psivloop, psivlop0, dpsivlop, pvbar,     &
            bpsqsurf, betapmhd, betatmhd
      character*8   codeid, machinei, eqgrdsze
!     in fluxav we calculate
!                    rm2 = <1/R**2>,rm2i=<R**2>,ravg=<R>,ravgi=<1/R>
!                    bsqinvavg=<Bt0**2/B**2>
!     based on npsi psi grid points with (1) corresponding to the
!     plasma edge and (npsi) corresponding to the magnetica axis.
!     These quantitties are interpolated onto the psir(1..nj) psi grid
!     rm2       ==> <R0**2/R**2>       ==> r2cap(1..nj)
!     rm2i      ==> <R**2>             ==> r2capi(1..nj)
!     ravg      ==> <R>                ==> rcap(1..nj)
!     ravgi     ==> <1./R>             ==> rcapi(1..nj)
!     bsqinvavg ==> <Bt0**2/B**2>      ==> bsqncap
!     bsq_avg   ==> <B**2/Bto**2>      ==> bsq_avg_cap
!     b_avg     ==> <B/Bt0>            ==> b_avg_cap
!     h_factr   ==> <SQRT(1. -B/Bmax)> ==> h_factr_cap
!     etc
     END MODULE geom
