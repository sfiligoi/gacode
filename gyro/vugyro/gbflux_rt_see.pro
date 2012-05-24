;
; C. Holland, UCSD
; v1.0: May 10, 2012: plots flux-surface avg fluxes as function of
; (r,t)
;
PRO gbflux_rt_event, the_event

  common GLOBAL

  WIDGET_CONTROL, the_event.top, GET_UVALUE = state, /NO_COPY
  WIDGET_CONTROL, the_event.id, get_uvalue=uvalue

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  gbflux_tag,title,pname,ytitle,units1,units2,units3
  pname = 'gbflux_rt-'+pname

  WSET, state.winID
  CASE (uvalue) OF 

     'PLOT': begin
        gbflux_rt_plot
     end   

     'spec_btn': begin
         i_spec = i_spec + 1
         if (i_spec ge n_kinetic) then i_spec = 0
         gbflux_rt_plot
     end

     'moment_btn': begin
         i_moment = i_moment+1
         if (i_moment ge p_moment) then i_moment = 0
         gbflux_rt_plot
     end

     'field_btn': begin
         i_f = i_f+1
         if (i_f ge n_field) then i_f = -1
         gbflux_rt_plot
     end

     'units_btn': begin
         if (exists_exp_derived eq 1) then begin
             i_units = i_units+1
             if (i_units eq 3) then i_units = 0
         endif
         gbflux_rt_plot
     end
 
     'ssplus_btn': begin
         n_ss = n_ss+dn_ss
         if (n_ss gt (n_n-1)*dn_ss) then n_ss=0
         gbflux_rt_plot
     end
    
     'ssminus_btn': begin
         n_ss = n_ss-dn_ss
         if (n_ss lt 0) then n_ss=(n_n-1)*dn_ss
         gbflux_rt_plot
     end

     'PS_DUMP': begin
         plot_mode = 2
         gbflux_rt_plot
     end

     'SAV': BEGIN
          savefilename = DIALOG_PICKFILE(FILE = pname+'.sav',/WRITE,$
                                      /OVERWRITE_PROMPT)
         IF (STRLEN(savefilename) NE 0) THEN BEGIN
             gbflux_rt_plot, SAVEFILE = savefilename
         ENDIF ELSE MESSAGE, '.sav save cancelled', /INFO
     END

     'DONE': WIDGET_CONTROL, the_event.top, /DESTROY

  ENDCASE

  IF (uvalue NE 'DONE') THEN $
    WIDGET_CONTROL, the_event.top, SET_UVALUE = state, /NO_COPY

END ;gbflux_rt_event

PRO gbflux_rt_see, group=group

  common GLOBAL

  ;;---------------------------------------------------------
  ;; Return conditions
  if (exists_gbflux_i EQ 0) then return ;no gbflux_i.out file read
  if xregistered('gbflux_rt_see') then return
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
                    value='Plot', $
                    uvalue='PLOT')

  x = widget_button(row1, $
                    value='+species', $
                    uvalue='spec_btn')

  x = widget_button(row1, $
                    value='+moment', $
                    uvalue='moment_btn')

  x = widget_button(row1, $
                    value='+field', $
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

  xmanager,'gbflux_rt_see', $
           base, $
           event='gbflux_rt_event',$
           group_leader=group

END ;gbflux_rt_see
