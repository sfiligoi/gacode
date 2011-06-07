c&gks_out.m  02-June-02 G. Staebler
c
c
      integer jmx
      parameter (jmx=200)
      character(15) version
      real *8 phi_bar_k,ne_bar_k
      real *8 te_bar_k,ti_bar_k
      real *8 ne_te_phase_k(4)
      real *8 mhd_DR_k
      real *8 anrate_m(0:jmx),dnrate_m(0:jmx),anfreq_m(0:jmx)
      real *8 dnfreq_m(0:jmx),anrate_sum(0:jmx),dnrate_sum(0:jmx)
      real *8 gamma_ion(0:jmx),gamma_electron(0:jmx)
      real *8 freq_ion(0:jmx),freq_electron(0:jmx)
      real *8 ky_m(0:jmx),aky_m(0:jmx)
      real *8 peflx_m(0:jmx),qeflx_m(0:jmx),qiflx_m(0:jmx)
      real *8 gamma_mks(0:jmx),dgamma_mks(0:jmx)
      real *8 freq_mks(0:jmx),ky_mks(0:jmx)
      real *8 zpte_crit_m(0:jmx),zpti_crit_m(0:jmx)
      real *8 phi_bar_m(0:jmx),ne_bar_m(0:jmx)
      real *8 te_bar_m(0:jmx),ti_bar_m(0:jmx)
      real *8 ne_te_phase_ion(0:jmx),ne_te_phase_electron(0:jmx)
      real *8 mhd_DR_m(0:jmx)
c
      common /gks_out/ phi_bar_k,ne_bar_k
     & ,te_bar_k,ti_bar_k,ne_te_phase_k
     & ,mhd_DR_k,ky_m,aky_m
     & ,anrate_m,dnrate_m,anfreq_m,dnfreq_m
     & ,zpte_crit_m,zpti_crit_m
     & ,peflx_m,qeflx_m,qiflx_m
     & ,gamma_mks,dgamma_mks,freq_mks,ky_mks
     & ,gamma_ion,freq_ion,gamma_electron,freq_electron
     & ,phi_bar_m,ne_bar_m,te_bar_m,ti_bar_m
     & ,ne_te_phase_ion,ne_te_phase_electron,mhd_DR_m
     & ,version
c
