; "cd" to the sim location so "spawn,'pwd'" below works
; in IDL, do "@comp_auto.pro" first
; then @batch_auto.pro

pro main_auto,sim_dir,tavg_s,tavg_e,noload=noload,save_it=save_it
; setting noload skips the loadsim call
; need to add variables to save info for making sav_* directories
;  use keyword save_it with the last entry in batch_auto IF it is /no_load

  common SAV_DIR, sim_base,sav_dir_path, sav_dir_name, sav_dir_flag

  common GLOBAL
; the next two are needed by spectrum_n_event
  common PLOT_VARIABLES
  common PRIVATE_SPECTRUM_N

in_flag=0
if (n_elements(sim_dir) le 0) then in_flag=1
if (n_elements(tavg_s) le 0) then in_flag=1
if (n_elements(tavg_e) le 0) then in_flag=1
if (in_flag ne 0 ) then begin
  print," Quitting, need four arguments: sim_dir,tavg_s,tavg_e"
   return
endif

; initialization
if(n_elements(sim_base) le 0) then begin
   spawn,'pwd',sim_base 
   print, " sim_base is ", sim_base
endif

; Should existing .ps and .idlout be moved to a sav_ directory?
if(n_elements(sav_dir_flag) le 0) then sav_dir_flag=0 ; initialization
sav_now=0
sav_end=0
sav_and_quit=0
if (n_elements(noload) ne 0 ) then begin
   sav_dir_flag=1 ; remember that this pass has /noload
; save the files made when the data files were read before:
   sav_now=1
;  use keyword save_it with the last entry in batch_auto IF it is /no_load
   if (n_elements(save_it) ne 0 ) then sav_end=1
endif else begin
; save the files made when the data files were last made with /noload:
   if(sav_dir_flag ne 0) then sav_now=1
; remember that this pass does not have /noload
   sav_dir_flag=0 
endelse

save_files_here:
if (sav_now ne 0) then begin
;  make sure sav_dir_name exists 
; look for "sav" in sav_dir_name to test if it is available?
   if(n_elements(sav_dir_name) le 0) then begin
      print," Quit: sav_dir_name is not set, but need to make sav_dir"
      print," problem with ",sav_dir_name
      sav_dir_flag=0 ; so the next sim_dir will run
      return
   endif
   if(n_elements(sav_dir_path) le 0) then begin
      print," Quit: sav_dir_path is not set, but need to make sav_dir"
      print," problem with ",sav_dir_path
      sav_dir_flag=0 ; so the next sim_dir will run
      return
   endif
   SIM_LOC=sim_base+"/"+sav_dir_path
   cd,SIM_LOC,current=curr_dir
   print,' Now do   "mkdir '+sav_dir_name+'"  and save files'
   spawn,'mkdir ' + sav_dir_name,exit_status=exstat
   if(exstat ne 0) then begin
      print," Quit: problem with   mkdir "+sav_dir_name
      print," problem with ",sav_dir_path+'  '+sav_dir_name
      sav_dir_flag=0 ; so the next sim_dir will run
      return
   endif
; Not for general use, next block with packtsumm:
   spawn,'mv *.ptar ' + sav_dir_name,exit_status=exstat
   if(exstat ne 0) then begin
      print," Quit: problem with   mv *.ptar "+sav_dir_name
      return
   endif
   spawn,'mv *.ps *.idlout ' + sav_dir_name,exit_status=exstat
   if(exstat ne 0) then begin
      print," Quit: problem with mv *.ps *.idlout "+sav_dir_name
      return
   endif
   cd,curr_dir
   if (sav_and_quit ne 0) then return ; don't need to remake .ps, ...
endif

;  Begin working on the current directory
SIM_PATH=sim_base+"/"+sim_dir
if (n_elements(noload) eq 0 ) then loadsim,SIM_PATH

; set up variables for quick testing:
;n_time1=800
;t=indgen(n_time1+1)

print,"t_min,t_max,n_time1= ",t_min,t_max,n_time1

t_min=tavg_s
; use tavg_e=-1. to not reset the averaging end time:
if (tavg_e gt tavg_s) then t_max=tavg_e else t_max=t[n_time1]
; reset time-interval string and indices:

t_indices,t_min,t_max,it1,it2
get_t_string

; Finished using old sav_dir_name; save new sav_dir_name for future use:
sav_dir_name='sav_' + strtrim(string(tavg_s,format='(I)'),2) + '_'+ $
strtrim(string(t_max,format='(I)'),2)
;debug print, " sav_dir_name is ",sav_dir_name
sav_dir_path=sim_dir ; to be used in next call to main_auto
;debug print, " sav_dir_path is ",sav_dir_path

cd,SIM_PATH,current=curr_dir

; for quick debugging, make two files:
; spawn,'touch test.ps'
; spawn,'touch test.idlout'

; Now make .ps and .idlout files automatically:

; spectrum_np_ave.ps
;  spectrum_np_ave_plot writes .ps automatically, but the next line
;  insures the same averaging window is used for this plot, too.
plot_mode=2 & i_log=1 & spectrum_np_ave_plot

; t_error.ps
plot_mode=2
t_error_event

; spectrum_n_all.ps
plot_mode=2
active_spectrum_n = 2
title  = '!3RMS Zonal and finite-n potentials'
ytitle = '!4<u>!3(t)'
pname0  = 'spectrum_n_all'
zoom=1.
spectrum_n_event

; New 'auto' procedures loop over species and moment:
; gbflux-elec-density-phi.ps; etc.
gbflux_auto
; gbflux_i-elec-energy-phi.ps; etc.
;   gb_i_overlay flag controls Nun_pl=2 in gbflux_i_plot.pro
gb_i_overlay=1 ; non-zero makes overlays
gbflux_i_auto_QG
; elec-energy-phi.ps; etc.
gbflux_n_auto
; pwr_n_potential.ps; etc.
midplane_fluc_auto
; osc-dtemperature-elec.ps
osc_auto

; packtsumm belongs at the end:
; Not for general use next line:
spawn,'packtsumm -q '+sim_dir+' d',exit_status=exstat
if(exstat ne 0) then begin
   print," Quit: problem with packtsumm "+sim_dir
   return
endif
;  end of packtsumm block

if((sav_now ne 0) and (sav_end ne 0)) then begin
; this should occur the end of batch_auto, so forget there is a /noload
; so that the next execution will start from scratch:
   sav_dir_flag=0 
   sav_and_quit=1 ; prevents repeat creation of .ps, .idlout, .ptar
   sav_end=0 ; prevents endless loop with next goto:
   goto,save_files_here
endif

cd,curr_dir

end
