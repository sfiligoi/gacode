; "cd" to the sim location so "spawn,'pwd'" in main_auto.pro works
; in IDL, do "@comp_auto.pro" first
; then @batch_auto.pro

; run vugyro to set up common blocks before compiling routines that
; use the COMMONs:
;  do this inside comp_auto.pro:   vugyro,/quick_exit

; define the home for sim dirs:
auto_home_def,'/global/scratch/sd/u4146/gacode/sim'

; main_auto,'',tavg_s,tavg_e,noload=noload
;  comment out lines with "noload" to save other files
;  execute them by hand after saving other files


main_auto,'TP27DHB_2k15/P27D_2_gT10',200.,-1.
main_auto,'TP27DHB_2k15/P27DH_2',200.,-1.
main_auto,'TP27DHB_2k15/P27DHBe_2',200.,-1.
main_auto,  'TP27DHB_2k15/P27DHBe_2',200.,350.,/noload

main_auto,'TP27DHB_k15/P27D_2_gT10',200.,-1.
main_auto,  'TP27DHB_k15/P27D_2_gT10',800.,-1.,/noload
main_auto,  'TP27DHB_k15/P27D_2_gT10',200.,600.,/noload
main_auto,'TP27DHB_k15/P27DH_2',200.,-1.
main_auto,  'TP27DHB_k15/P27DH_2',500.,900.,/noload
main_auto,'TP27DHB_k15/P27DHBe_2',200.,-1.
main_auto,  'TP27DHB_k15/P27DHBe_2',500.,-1.,/noload,/save_it
;  use keyword /save_it with the last entry ONLY IF it is /no_load

; NOT in a non-procedure:  end
