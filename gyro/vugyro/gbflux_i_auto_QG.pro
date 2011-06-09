; Makes both .ps and .idlout
; gbflux_i-elec-density-phi.ps
; gbflux_i-elec-energy-phi.ps
; gbflux_i-ion1-density-phi.ps
; gbflux_i-ion1-energy-phi.ps
; ... for all ions

pro gbflux_i_auto_QG

  common GLOBAL
  common PRIVATE_GBFLUX_I  

; Initialize variables that control plots
i_f=0       ; up to n_field
i_zero=0    ; 0 or 1
zoom=1.
i_units=0   ; 0, 1, 2
n_ss=0


;  Do particle fluxes first, then energy fluxes:
; i_moment=0  ; up to p_moment-1; density, energy, momentum, exchange
; i_spec=0    ; up to n_kinetic-1, electrons are last

for i_moment=0,1 do begin
  i_zero=1-i_moment ; works only for particle and energy fluxes

  for i_spec=0,n_kinetic-1 do begin 
;   Stack IDL commands with &
; plot_mode: anything but 1 produces a .ps file
; Write .idlout file with plot_export=1; can also make .ps file
     plot_mode=2 & plot_export=1 & gbflux_i_plot
  endfor ; loop over species
endfor ; loop over moments

end
