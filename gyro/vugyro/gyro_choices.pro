;;--------------------------------------
;; Choices for toplevel menu
;;--------------------------------------

pro gyro_choices, event

  common GLOBAL

  widget_control, get_uvalue=control, event.id

  ;;--------------------------------------------
  ;; Test if control is a string beginning with 
  ;; a digit.  If so, load directory:
  ;;
  if (control eq 'sim_pick_dialog') then begin

     if (remotedir_flag eq 0) then begin
        simroot = getenv('PWD')
     endif else begin
        simroot = remotedir
     endelse

     chooser: picked_dir = DIALOG_PICKFILE(/DIRECTORY,PATH=simroot, $
                                           TITLE='Select a valid GYRO directory')

     ;; Exit chooser if cancel is selected 
     if picked_dir eq '' then return

     ;; Test to see if simdir is valid 
     spawn,getenv("GACODE_ROOT")+'/gyro/bin/'+'simdir_test '+picked_dir,exit_status=status
     if status eq 1 then begin
        ;; Load if valid
        loadsim, picked_dir 
     endif else begin
        ;; Go back to chooser if directory not valid
        simroot=picked_dir
        goto, chooser
     endelse

     return

  endif else if (control le '9') then begin

     loadsim, dir(fix(control))
     return

  endif
  ;;--------------------------------------------

  case control of

     'exit': widget_control, event.top,$
       /destroy

     'view_run': view_run

     'view_efficiency': view_efficiency

     else: $

     returnvalue = execute(control+',group=event.top')

  endcase

end


