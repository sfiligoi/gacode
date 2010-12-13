PRO midplane_avgkspect_event, the_event
;
; C. Holland, UCSD
; Plotting routine for improved midplane <S(kx,ky)>
;
; v1.0: Jan 26, 2007- debugged (I think), ready for release into the
;                     wild
;
; v1.1: Jul 19 2007- updated to use w/ vuGyro plotting options, no_kr option
;
; v1.2: Aug 16, 2007- added use_n feature to allow n or r-dep k_theta plots
;
  common GLOBAL
  common MIDPLANE_DATA

  WIDGET_CONTROL, the_event.top, GET_UVALUE = state, /NO_COPY
  WIDGET_CONTROL, the_event.id, get_uvalue=uvalue

  ;add field name in here
  pname = 'midplane_' + tag[state.i_pwr] + '_avgkspect'

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  WSET, state.winID
  CASE (uvalue) OF 

     'PLOT': BEGIN
        midplane_avgkspect_plot, state.i_pwr, state.krmax, $
          PLOTMIN = state.plotmin, PLOTMAX = state.plotmax, $
          NO_KR = state.xaxis_style, USE_N = state.use_n
     END

     'i_pwr++': BEGIN
         IF (state.i_pwr LT n_pwr-1) THEN BEGIN
             state.i_pwr++
             midplane_avgkspect_plot, state.i_pwr, state.krmax, $
               PLOTMIN = state.plotmin, PLOTMAX = state.plotmax, $
               NO_KR = state.xaxis_style, USE_N = state.use_n
         ENDIF
     END

     'i_pwr--': BEGIN
         IF (state.i_pwr GT 0) THEN BEGIN
             state.i_pwr--
             midplane_avgkspect_plot, state.i_pwr, state.krmax, $
               PLOTMIN = state.plotmin, PLOTMAX = state.plotmax, $
               NO_KR = state.xaxis_style, USE_N = state.use_n
         ENDIF
     END

     'xaxis_style': BEGIN
        state.xaxis_style = 1 - state.xaxis_style
        midplane_avgkspect_plot, state.i_pwr, state.krmax, $
          PLOTMIN = state.plotmin, PLOTMAX = state.plotmax, $
          NO_KR = state.xaxis_style, USE_N = state.use_n
     END

     'use_n': BEGIN
        state.use_n = 1 - state.use_n
        midplane_avgkspect_plot, state.i_pwr, state.krmax, $
          PLOTMIN = state.plotmin, PLOTMAX = state.plotmax, $
          NO_KR = state.xaxis_style, USE_N = state.use_n
     END

     'MAXKR': BEGIN
         newkr = the_event.value
         IF (newkr LE 0.0) THEN BEGIN
             MESSAGE, 'new kr max must be greater than zero,' + $
                      'restoring to default value', /INFO
             WIDGET_CONTROL, the_event.id, SET_VALUE = MAX(kt_rho)
         ENDIF ELSE IF (newkr GE MAX(kr_rho)) THEN BEGIN
             MESSAGE, 'new kr max must be less than ' + $
                      NUMTOSTRING(MAX(kr_rho)) + $
                      ', restoring to default value', /INFO
            WIDGET_CONTROL, the_event.id, SET_VALUE = MAX(kt_rho)
         ENDIF ELSE BEGIN
             state.krmax = newkr
             midplane_avgkspect_plot, state.i_pwr, state.krmax, $
               PLOTMIN = state.plotmin, PLOTMAX = state.plotmax, $
               NO_KR = state.xaxis_style, USE_N = state.use_n
         ENDELSE
     END

     'PLOTMIN': BEGIN
         state.plotmin = the_event.value
         midplane_avgkspect_plot, state.i_pwr, state.krmax, $
           PLOTMIN = state.plotmin, PLOTMAX = state.plotmax, $
           NO_KR = state.xaxis_style, USE_N = state.use_n
     END

     'PLOTMAX': BEGIN
         state.plotmax = the_event.value
         midplane_avgkspect_plot, state.i_pwr, state.krmax, $
           PLOTMIN = state.plotmin, PLOTMAX = state.plotmax, $
           NO_KR = state.xaxis_style, USE_N = state.use_n
    END


     'PS_DUMP': begin
         plot_mode = 2
         midplane_avgkspect_plot, state.i_pwr, state.krmax, $
           PLOTMIN = state.plotmin, PLOTMAX = state.plotmax, $
           NO_KR = state.xaxis_style, USE_N = state.use_n
     end

     'SAV': BEGIN
         savefilename = DIALOG_PICKFILE(FILE = pname+'.sav',/WRITE,$
                                      /OVERWRITE_PROMPT)
         IF (STRLEN(savefilename) NE 0) THEN BEGIN
             midplane_avgkspect_plot, state.i_pwr, state.krmax, $
               PLOTMIN = state.plotmin, PLOTMAX = state.plotmax, $
               SAVEFILE = savefilename, NO_KR = state.xaxis_style, $
               USE_N = state.use_n
         ENDIF ELSE MESSAGE, '.sav save cancelled', /INFO
     END

     'ASCII': BEGIN
         asciifilename = DIALOG_PICKFILE(FILE = pname+'.dat',/WRITE,$
                                      /OVERWRITE_PROMPT)
         IF (STRLEN(asciifilename) NE 0) THEN BEGIN
             midplane_avgkspect_plot, state.i_pwr, state.krmax, $
               PLOTMIN = state.plotmin, PLOTMAX = state.plotmax, $
               ASCIIFILE=asciifilename, NO_KR = state.xaxis_style, $
               USE_N = state.use_n
         ENDIF ELSE MESSAGE, 'text save cancelled', /INFO
     END

     'DONE': WIDGET_CONTROL, the_event.top, /DESTROY

  ENDCASE

  IF (uvalue NE 'DONE') THEN $
    WIDGET_CONTROL, the_event.top, SET_UVALUE = state, /NO_COPY

END ;midplane_avg_kspect_event

PRO midplane_avgkspect_see, group=group
;
; C. Holland, UCSD
; Plotting routine for improved midplane <S(kx,ky)>
;
; v1.0: Jan 26, 2007- debugged (I think), ready for release into the
;                     wild
;

  common GLOBAL

  ;;---------------------------------------------------------
  ;; Return conditions
  if (N_ELEMENTS(exists_u) EQ 0) then return ;no data loaded
  if (exists_u EQ 0) then return ;no u.out file read
  if xregistered('midplane_avgkspect_see') then return
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
                    value='field up', $
                    uvalue='i_pwr++')

  x = widget_button(row1, $
                    value='field dn', $
                    uvalue='i_pwr--')

  x = widget_button(row1, $
                    value='kr/r', $
                    uvalue='xaxis_style')

  x = widget_button(row1, $
                    value='kt/n', $
                    uvalue='use_n')

  x = cw_field(row1, /RETURN_EVENTS, TITLE = 'Max kr', $
               /FLOATING, VALUE = MAX(kt_rho), XSIZE=8, $
               UVALUE = 'MAXKR')

  x = cw_field(row1, /RETURN_EVENTS, TITLE = 'Min power', $
               /FLOATING, VALUE = -10., XSIZE=8, $
               UVALUE = 'PLOTMIN')

  x = cw_field(row1, /RETURN_EVENTS, TITLE = 'Max power', $
               /FLOATING, VALUE = 0., XSIZE=8, $
               UVALUE = 'PLOTMAX')

  row2 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row2, $
                    value='PS dump', $
                    uvalue='PS_DUMP')

  x = widget_button(row2, $
                    value='Write to .sav file', $
                    uvalue='SAV')

  x = widget_button(row2, $
                    value='Write to text file', $
                    uvalue='ASCII')

  x = widget_button(row2, $
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
  state = {winID:winID, i_pwr:0, xaxis_style: 0, krmax:MAX(kt_rho), $
           plotmin:-10., plotmax:0., use_n:0}
  WIDGET_CONTROL, base, SET_UVALUE = state, /NO_COPY

  xmanager,'midplane_avgkspect_see', $
           base, $
           event='midplane_avgkspect_event',$
           group_leader=group

END ;midplane_avgkspect_see
