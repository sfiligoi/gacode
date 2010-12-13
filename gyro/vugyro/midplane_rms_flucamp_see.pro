PRO midplane_rms_flucamp_event, the_event

  common GLOBAL
  common MIDPLANE_DATA

  WIDGET_CONTROL, the_event.top, GET_UVALUE = state, /NO_COPY
  WIDGET_CONTROL, the_event.id, get_uvalue=uvalue

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  pname = 'midplane_rms_flucamp_'+tag[i_pwr]
  CASE (uvalue) OF 

     'PLOT': begin
        WSET, state.winID
        midplane_rms_flucamp_plot, LOCAL_NORM = state.lnorm, $
                           FINITE_N = state.finite_n
     end   

     'idx_up': begin
         counter_up, i_pwr, n_pwr-1,1
         midplane_rms_flucamp_plot, LOCAL_NORM = state.lnorm, $
                           FINITE_N = state.finite_n
     end
     
     'idx_dn': begin
         counter_dn, i_pwr, 0, 1
         midplane_rms_flucamp_plot, LOCAL_NORM = state.lnorm, $
                           FINITE_N = state.finite_n
     end

     'lnorm': begin
         state.lnorm = 1 - state.lnorm
         midplane_rms_flucamp_plot, LOCAL_NORM = state.lnorm, $
                           FINITE_N = state.finite_n
     end

     'finite_n': begin
         state.finite_n = 1 - state.finite_n
         midplane_rms_flucamp_plot, LOCAL_NORM = state.lnorm, $
                           FINITE_N = state.finite_n
     end

     'ssplus_btn': begin
         n_ss = n_ss+dn_ss
         if (n_ss gt (n_n-1)*dn_ss) then n_ss=0
         midplane_rms_flucamp_plot, LOCAL_NORM = state.lnorm, $
                           FINITE_N = state.finite_n
    end
     
     'ssminus_btn': begin
         n_ss = n_ss-dn_ss
         if (n_ss lt 0) then n_ss=(n_n-1)*dn_ss
         midplane_rms_flucamp_plot, LOCAL_NORM = state.lnorm, $
                           FINITE_N = state.finite_n
     end

     'finte_n_btn': begin
         state.finite_n = 1 - state.finite_n
         midplane_rms_flucamp_plot, LOCAL_NORM = state.lnorm, $
                           FINITE_N = state.finite_n
     end

     'PS_DUMP': begin
         plot_mode = 2
         midplane_rms_flucamp_plot, LOCAL_NORM = state.lnorm, $
                           FINITE_N = state.finite_n
     end

     'SAV': BEGIN
          savefilename = DIALOG_PICKFILE(FILE = pname+'.sav',/WRITE,$
                                      /OVERWRITE_PROMPT)
         IF (STRLEN(savefilename) NE 0) THEN BEGIN
             midplane_rms_flucamp_plot, LOCAL_NORM = state.lnorm, $
                               FINITE_N = state.finite_n, $
                               SAVEFILE = savefilename
         ENDIF ELSE MESSAGE, '.sav save cancelled', /INFO
     END

     'ASCII': BEGIN
          asciifilename = DIALOG_PICKFILE(FILE = pname+'.dat',/WRITE,$
                                      /OVERWRITE_PROMPT)
         IF (STRLEN(asciifilename) NE 0) THEN BEGIN
             midplane_rms_flucamp_plot, LOCAL_NORM = state.lnorm, $
                               FINITE_N = state.finite_n, $
                               ASCIIFILE=asciifilename
         ENDIF ELSE MESSAGE, 'text save cancelled', /INFO
     END

     'DONE': WIDGET_CONTROL, the_event.top, /DESTROY

  ENDCASE

  IF (uvalue NE 'DONE') THEN $
    WIDGET_CONTROL, the_event.top, SET_UVALUE = state, /NO_COPY

END ;midplane_rms_flucamp_event

PRO midplane_rms_flucamp_see, group=group

  common GLOBAL

  ;;---------------------------------------------------------
  ;; Return conditions
;  if (N_ELEMENTS(exists_u) EQ 0) then return ;no data loaded
;  if (exists_u EQ 0) then return ;no u.out file read
  if xregistered('midplane_rms_flucamp_see') then return
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
                    value='index up', $
                    uvalue='idx_up')

  x = widget_button(row1, $
                    value='index dn', $
                    uvalue='idx_dn')

  x = widget_button(row1, $
                    value='Local norm', $
                    uvalue='lnorm')

  x = widget_button(row1, $
                    value='Finite-n only', $
                    uvalue='finite_n')

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
  state = {winID:winID, lnorm:0, finite_n: 0}
  WIDGET_CONTROL, base, SET_UVALUE = state, /NO_COPY

  xmanager,'midplane_rms_flucamp_see', $
           base, $
           event='midplane_rms_flucamp_event',$
           group_leader=group

END ;midplane_rms_flucamp_see
