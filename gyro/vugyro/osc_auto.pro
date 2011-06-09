; Makes both .ps and .idlout:
; osc-density-elec.ps
; osc-ddensity-elec.ps
; osc-temperature-elec.ps
; osc-dtemperature-elec.ps
; osc-density-ion1.ps
; osc-ddensity-ion1.ps
; osc-temperature-ion1.ps
; osc-dtemperature-ion1.ps
; ... for all ions

pro osc_auto

  common GLOBAL
  common ZMOMENT_DATA
  common PRIVATE_OSC
  
 ; Initialize variables that control plots
i_f=0       ; up to n_field
i_zero=1    ; 0 or 1, 0 forces Ymin=0.
;  modified osc_plot attempts to include entire range of fluct variable
zoom=1.
i_units=0   ; 0, 1, 2
n_ss=0
i_div=0 ; ???


; i_gradient: 1 for gradient, 0 for function
equil_flag=1  ; multiplies background profile
i_ss_mn_den_plot=0
n_ss_bin=0
i_ss_mn=0

;  Do particle fluxes first, then energy fluxes:
i_moment=0  ; up to p_moment-1; density, energy, momentum, exchange
i_spec=0    ; up to n_kinetic-1, electrons are last

for i_gradient=0,1 do begin
  for i_moment=0,1 do begin
    for i_spec=0,n_kinetic-1 do begin 
;     Stack IDL commands with &
; plot_mode: anything but 1 produces a .ps file
; Write .idlout file with plot_export=1; can also make .ps file
      plot_mode=2 & plot_export=1 & osc_plot
    endfor ; loop over species
  endfor ; loop over moments
endfor ; loop over f/grad-f

end
