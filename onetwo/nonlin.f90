     MODULE nonlin
      USE param,only : kj,kk
      IMPLICIT none
      integer,parameter  :: ku = kj*kk
      integer                                                      &
              nrmax, nequations,                                   &
              iteration_method,non_lin_method,                     &
              bandwidth,maxfev,fdigits,switch_method,              &
              get_typf,jacobian_type,tot_iters_max,                &
              ihalve_nl,freeze_nl,freeze_type,                     &
              iters_freeze       
      real*8                                                       &
              xscale(ku),resid(ku),diag(ku),                       &
             fvec(ku),tresid,ssqrmin,steptol,                      &
             gradtol,fvectol,typf(ku)                              
!
     END MODULE nonlin 
