;
; C.Holland, UCSD
; v1.0: July 19, 2007
;
; Plots midplane (r,alpha) correlation functions, centered at
; i_r = N_R/2
;
;
; 11/9/2007: added finite-n option
;
PRO midplane_corr_event, the_event

  common GLOBAL
  common MIDPLANE_DATA

  WIDGET_CONTROL, the_event.top, GET_UVALUE = state, /NO_COPY
  WIDGET_CONTROL, the_event.id, GET_UVALUE=uvalue

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  pname = 'midplane_corr'
  WSET, state.winID

  CASE (uvalue) OF 

     'PLOT': begin
         midplane_corr_plot, state.i_r0, state.which_plot, $
                             LOCAL_CORR = state.local, $
                             FINITE_N = state.finite_n
     end   

     'which_plot': begin
         state.which_plot = the_event.index
         midplane_corr_plot, state.i_r0, state.which_plot, $
                             LOCAL_CORR = state.local, $
                             FINITE_N = state.finite_n
     end

     'local': begin
         state.local = 1 - state.local
         midplane_corr_plot, state.i_r0, state.which_plot, $
                             LOCAL_CORR = state.local, $
                             FINITE_N = state.finite_n
     end

     'finite_n': begin
         state.finite_n = 1-state.finite_n
         midplane_corr_plot, state.i_r0, state.which_plot, $
                             LOCAL_CORR = state.local, $
                             FINITE_N = state.finite_n
     end

     'idx_up': begin
         counter_up,i_pwr,n_pwr-1,1
         midplane_corr_plot, state.i_r0, state.which_plot, $
                             LOCAL_CORR = state.local, $
                             FINITE_N = state.finite_n
     end   

     'idx_dn': begin
         counter_dn,i_pwr,0,1
         midplane_corr_plot, state.i_r0, state.which_plot, $
                             LOCAL_CORR = state.local, $
                             FINITE_N = state.finite_n
     end   

     'ir0++': begin
         i_r0 = state.i_r0
         counter_up,i_r0,N_R-1,1
         state.i_r0 = i_r0
         midplane_corr_plot, state.i_r0, state.which_plot, $
                             LOCAL_CORR = state.local, $
                             FINITE_N = state.finite_n
     end   

     'ir0--': begin
         i_r0 = state.i_r0
         counter_dn,i_r0,0,1
         state.i_r0 = i_r0
         midplane_corr_plot, state.i_r0, state.which_plot, $
                             LOCAL_CORR = state.local, $
                             FINITE_N = state.finite_n
     end   

     'PS_DUMP': begin
         plot_mode = 2
         midplane_corr_plot, state.i_r0, state.which_plot, $
                             LOCAL_CORR = state.local, $
                             FINITE_N = state.finite_n
     end

     'SAV': BEGIN
         savefilename = DIALOG_PICKFILE(FILE = pname+'.sav',/WRITE,$
                                      /OVERWRITE_PROMPT)
         IF (STRLEN(savefilename) NE 0) THEN BEGIN
             midplane_corr_plot, state.i_r0, state.which_plot, $
                                 LOCAL_CORR = state.local, $
                                 FINITE_N = state.finite_n, $
                                 SAVEFILE = savefilename
         ENDIF ELSE MESSAGE, '.sav save cancelled', /INFO
     end

     'ASCII': BEGIN
         asciifilename = DIALOG_PICKFILE(FILE = pname+'.dat',/WRITE,$
                                      /OVERWRITE_PROMPT)
         IF (STRLEN(asciifilename) NE 0) THEN BEGIN
             midplane_corr_plot, state.i_r0, state.which_plot, $
                                 LOCAL_CORR = state.local, $
                                 FINITE_N = state.finite_n, $
                                 ASCIIFILE=asciifilename
         ENDIF ELSE MESSAGE, 'text save cancelled', /INFO
     end

     'DONE': WIDGET_CONTROL, the_event.top, /DESTROY

  ENDCASE

  IF (uvalue NE 'DONE') THEN $
    WIDGET_CONTROL, the_event.top, SET_UVALUE = state, /NO_COPY

END ;midplane_corr_event

PRO midplane_corr_see, group=group

  common GLOBAL

  ;;---------------------------------------------------------
  ;; Return conditions
  if (N_ELEMENTS(exists_u) EQ 0) then return ;no data loaded
  if (exists_u EQ 0) then return ;no u.out file read
  if (n0 NE 0) then return ;weed out linear or strange NL runs
  if xregistered('midplane_corr_see') then return
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

  x = widget_droplist(row1, $
                      value = ['C(dr)', $
                               'C(da)',$
                               'C(dr,da)'], $
                      uvalue = 'which_plot')
                      
  x = widget_button(row1, $
                    value='Local/box-avg', $
                    uvalue='local')

  x = widget_button(row1, $
                    value='Finite-n only', $
                    uvalue='finite_n')

  x = widget_button(row1, $
                    value='index up', $
                    uvalue='idx_up')

  x = widget_button(row1, $
                    value='index dn', $
                    uvalue='idx_dn')

  x = widget_button(row1, $
                    value='i_r0++', $
                    uvalue='ir0++')

  x = widget_button(row1, $
                    value='i_r0--', $
                    uvalue='ir0--')

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
  state = {winID:winID, i_r0:N_R/2, which_plot:0, local: 1, finite_n: 1}
  WIDGET_CONTROL, base, SET_UVALUE = state, /NO_COPY

  xmanager,'midplane_corr_see', $
           base, $
           event='midplane_corr_event',$
           group_leader=group

END ;midplane_corr_see
