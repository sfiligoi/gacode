; Make .ps files
; pwr_n_potential.ps
; pwr_n_density_elec.ps
; pwr_n_energy_elec.ps
; pwr_n_tem_elec.ps
; pwr_n_density_ion1.ps
; pwr_n_energy_ion1.ps
; pwr_n_tem_ion1.ps
; ... for all ions

pro midplane_fluc_auto

  common GLOBAL
  common MIDPLANE_DATA
  common POLOIDAL_DATA

t_c=n_time1    ; the end of the simulation

mUt = 4 ; fixed "Refine" value

; i_pwr runs 0 to n_pwr-1,   the index for what is plotted
for i_pwr=0,n_pwr-1 do begin
; plot_mode: anything but 1 produces a .ps file
  plot_mode=2 & midplane_fluc_plot
endfor

end
