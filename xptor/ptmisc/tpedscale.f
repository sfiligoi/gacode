      subroutine tpedscale(iproc,istep,time,pow_ped,
     >           iohm,palpha,pohm,pbrem,psync)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c Use power dependent pedestal scalings for time-dependent Tped
c including Ohmic and alpha heating while subtracting off radiation
c Note: radiation must be self-consistently computed (irad=-1 or 1)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
      include '../inc/input.m'
      include '../inc/data.m'
      include '../inc/tport.m'
      include '../inc/model.m'
c
       integer iproc, istep, iohm
       real*8 time, tocur, eps_exp
       real*8 palpha, pohm, pbrem, psync
       real*8 zpi, q95, q_cyl, q_sh, volume, vol_ped, convert,
     >        n_ped, powe_ped, powi_ped, pow_ped, w_ped(8)
       real*8 r_tar(1), q_tar(1)
c
c...setup
c
       zpi = atan2(0.D0,-1.D0)
       convert = 1.6022D-22
c
       tocur=dabs(tocur_d)*1.D-6  ! Ip (MA)
       eps_exp = amin_d / rmajor_exp
       r_tar(1)=0.95D0
       call w_lin_interp_r8(nj_d,rho_d,q_d,1,r_tar,q_tar)
       q95=q_tar(1)
       q_cyl = 5.D0*elonga_exp*amin_d**2.D0*bt_exp/rmajor_exp/tocur
       q_sh = q95/q_cyl
       volume = 2.D0*elonga_exp*
     >          zpi**2.D0*rmajor_exp*amin_d**2.D0
       volume = vol_exp(nj_d-1) ! above formula not accurate
       vol_ped = rho(jout_m)*volume
       n_ped = ne_exp(jout_m)
c
c...power thru the pedestal including alpha power
c
       if(iohm.eq.1) then
         powe_ped = powe_beam_exp(nj_d-1)+powe_rf_exp(nj_d-1)+
     >              powe_lh_exp(nj_d-1)+pohm
       else
         powe_ped = powe_beam_exp(nj_d-1)+powe_rf_exp(nj_d-1)+
     >              powe_lh_exp(nj_d-1)+powe_oh_exp(nj_d-1)
       endif
       powe_ped = powe_ped + palpha - pbrem - psync
       powi_ped = powi_beam_exp(nj_d-1)+powi_rf_exp(nj_d-1)
       
       pow_ped = powe_ped + powi_ped
c       write(*,*) 'Palpha = ',palpha
c       write(*,*) 'Pohmic = ',pohm
c       write(*,*) 'Pbrem = ',pbrem
c       write(*,*) 'Psync = ',psync
c
c Power dependent scaling - Eqn 2 from 2-term Cordey IAEA02 paper
c Auxiliary heating + Ohmic only, use line-averaged ne
c
       w_ped(3) = 7.32D-4*tocur**1.45D0*rmajor_exp**1.54D0*
     >            pow_ped**0.16D0*neline_exp**(0.08D0)*
     >            bt_exp**0.32D0*elonga_exp**1.74D0*
     >            (amin_d/rmajor_exp)**(-1.74D0)*
     >            amassgas_exp**0.2D0*q_sh**2.3D0
       t_ped(3) = w_ped(3)/(convert*3.D0*n_ped*1.D19*vol_ped)
c
c Power dependent scaling - Eqn 3 from 2-term Cordey IAEA02 paper
c w/o type III ELMy data
c Auxiliary heating + Ohmic only, use line-averaged ne
c
       w_ped(4) = 0.014D0*tocur**1.60D0*rmajor_exp**1.36D0*
     >            pow_ped**0.29D0*neline_exp**(-0.07D0)*
     >            bt_exp**0.31D0*amassgas_exp**0.2D0*q_sh**2.34D0
       t_ped(4) = w_ped(4)/(convert*3.D0*n_ped*1.D19*vol_ped)
c
c Power dependent scaling - Eqn 2 from 2-term Cordey IAEA02 paper
c w/ DeBoo correction included in 10/10/02 poster, RMSE=25.4%
c Auxiliary heating + Ohmic only, use line-averaged ne
c
       w_ped(5) = 5.095D-4*tocur**1.57D0*rmajor_exp**1.37D0*
     >            pow_ped**0.20D0*neline_exp**(0.03D0)*
     >            bt_exp**0.21D0*elonga_exp**1.84D0*
     >            (amin_d/rmajor_exp)**(-2.16D0)*
     >            amassgas_exp**0.2D0*q_sh**2.37D0
       t_ped(5) = w_ped(5)/(convert*3.D0*n_ped*1.D19*vol_ped)
c
       if(ibound.eq.5 .and. iproc.eq.0) 
     >    write(*,50) istep, time, pow_ped, t_ped(3)
       if(ibound.eq.6 .and. iproc.eq.0) 
     >    write(*,50) istep, time, pow_ped, t_ped(4)
c
 50    format(' step=',i4,', time=',f8.4,', Pped=',f8.3,
     >        ', Tped=',f8.4)
c
      return
      end
