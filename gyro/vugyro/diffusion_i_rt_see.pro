;
; C. Holland, UCSD
; v1.0: July 19, 2007: plots flux-surface avg chi/D as function of
; (r,t)
;
PRO diffusion_i_rt_event, the_event

  common GLOBAL

  WIDGET_CONTROL, the_event.top, GET_UVALUE = state, /NO_COPY
  WIDGET_CONTROL, the_event.id, get_uvalue=uvalue

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  smf_tag,title,pname,ytitle
  pname = 'diff_i-'+pname+'-rt'

  WSET, state.winID
  CASE (uvalue) OF 

     'PLOT': begin
        diffusion_i_rt_plot
     end   

     'spec_btn': begin
         i_spec = i_spec + 1
         if (i_spec gt n_kinetic-1) then i_spec = 0
         diffusion_i_rt_plot
     end

     'moment_btn': begin
         i_moment = 1 - i_moment
         diffusion_i_rt_plot
     end

     'field_btn': begin
         i_f = i_f+1
         if (i_f ge n_field) then i_f = -1
         diffusion_i_rt_plot
     end

     'units_btn': begin
         i_moment = 1 - i_moment
         diffusion_i_rt_plot
     end
 
     'ssplus_btn': begin
         n_ss = n_ss+dn_ss
         if (n_ss gt (n_n-1)*dn_ss) then n_ss=0
         diffusion_i_rt_plot
     end
    
     'ssminus_btn': begin
         n_ss = n_ss-dn_ss
         if (n_ss lt 0) then n_ss=(n_n-1)*dn_ss
         diffusion_i_rt_plot
     end

     'PS_DUMP': begin
         plot_mode = 2
         diffusion_i_rt_plot
     end

     'SAV': BEGIN
          savefilename = DIALOG_PICKFILE(FILE = pname+'.sav',/WRITE,$
                                      /OVERWRITE_PROMPT)
         IF (STRLEN(savefilename) NE 0) THEN BEGIN
             diffusion_i_rt_plot, SAVEFILE = savefilename
         ENDIF ELSE MESSAGE, '.sav save cancelled', /INFO
     END

     'ASCII': BEGIN
          asciifilename = DIALOG_PICKFILE(FILE = pname+'.dat',/WRITE,$
                                      /OVERWRITE_PROMPT)
         IF (STRLEN(asciifilename) NE 0) THEN BEGIN
             diffusion_i_rt_plot, ASCIIFILE=asciifilename
         ENDIF ELSE MESSAGE, 'text save cancelled', /INFO
     END

     'DONE': WIDGET_CONTROL, the_event.top, /DESTROY

  ENDCASE

  IF (uvalue NE 'DONE') THEN $
    WIDGET_CONTROL, the_event.top, SET_UVALUE = state, /NO_COPY

END ;diffusion_i_rt_event

PRO diffusion_i_rt_see, group=group

  common GLOBAL

  ;;---------------------------------------------------------
  ;; Return conditions
  if (exists_diff_i EQ 0) then return ;no diff_i.out file read
  if xregistered('diffusion_i_rt_see') then return
  ;;
  ;;----------------------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1, $
                    value='plot', $
                    uvalue='PLOT')

  x = widget_button(row1, $
                    value='i/e', $
                    uvalue='spec_btn')

  x = widget_button(row1, $
                    value='n/E', $
                    uvalue='moment_btn')

  x = widget_button(row1, $
                    value='es/em', $
                    uvalue='field_btn')

  x = widget_button(row1, $
                    value='units', $
                    uvalue='units_btn')

  x = widget_button(row1, $
                    value='ss+', $
                    uvalue='ssplus_btn')

  x = widget_button(row1, $
                    value='ss-', $
                    uvalue='ssminus_btn')

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue='PS_DUMP')

  x = widget_button(row1, $
                    value='Write to .sav file', $
                    uvalue='SAV')

  x = widget_button(row1, $
                    value='Write to text file', $
                    uvalue='ASCII')

  x = widget_button(row1, $
                    value='Done', $
                    uvalue='DONE')

  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base,     $
                     xsize=sx, $
                     ysize=sy)

  widget_control, base, $
                  set_uvalue=state,$
                  /no_copy, $
                  /realize


  widget_control, draw, $
                  get_value=winID

  ;associate state variable with window ID so we don't need another common 
  ;block
  state = {winID:winID}
  WIDGET_CONTROL, base, SET_UVALUE = state, /NO_COPY

  xmanager,'diffusion_i_rt_see', $
           base, $
           event='diffusion_i_rt_event',$
           group_leader=group

END ;diffusion_i_rt_see
