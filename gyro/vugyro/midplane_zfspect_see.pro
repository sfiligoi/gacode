PRO midplane_zfspect_event, the_event
;
; C. Holland, UCSD
; v1.0: 5/20/2008
;
; Event handler for midplane_zfspect_see
;

  common GLOBAL
  common MIDPLANE_DATA

  WIDGET_CONTROL, the_event.top, GET_UVALUE = state, /NO_COPY
  WIDGET_CONTROL, the_event.id, get_uvalue=uvalue

  IF (state.plot_shear EQ 1) THEN pname = 'midplane_zfshearspect_'+tag[i_pwr] $
  ELSE pname = 'midplane_zfspect_'+tag[i_pwr]
  WSET, state.winID
       
  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  CASE (uvalue) OF 

     'PLOT': begin
         midplane_zfspect_plot, state.plot_log, state.krmax, state.NFFT, $
           state.plot_shear
     end   

     'linlog': begin
         state.plot_log = 1 - state.plot_log
         midplane_zfspect_plot, state.plot_log, state.krmax, state.NFFT, $
           state.plot_shear
      end   

     'shear_toggle': begin
         state.plot_shear = 1 - state.plot_shear
         midplane_zfspect_plot, state.plot_log, state.krmax, state.NFFT, $
           state.plot_shear
      end   

     'idx_up': begin
         counter_up,i_pwr,n_pwr-1,1
         midplane_zfspect_plot, state.plot_log, state.krmax, state.NFFT, $
           state.plot_shear
      end   

     'idx_dn': begin
         counter_dn,i_pwr,0,1
         midplane_zfspect_plot, state.plot_log, state.krmax, state.NFFT, $
           state.plot_shear
      end   

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
             midplane_zfspect_plot, state.plot_log, state.krmax, state.NFFT, $
           state.plot_shear
         ENDELSE
     END

     'NFFT': BEGIN
         state.NFFT = the_event.value
         midplane_zfspect_plot, state.plot_log, state.krmax, state.NFFT, $
           state.plot_shear
     END

     'PS_DUMP': begin
         plot_mode = 2
         midplane_zfspect_plot, state.plot_log, state.krmax, state.NFFT, $
           state.plot_shear
     end

     'SAV': BEGIN
         savefilename = DIALOG_PICKFILE(FILE = pname+'.sav',/WRITE,$
                                      /OVERWRITE_PROMPT)
         IF (STRLEN(savefilename) NE 0) THEN BEGIN
             midplane_zfspect_plot, state.plot_log, state.krmax, $
               state.NFFT, state.plot_shear, SAVEFILE = savefilename
         ENDIF ELSE MESSAGE, '.sav save cancelled', /INFO
     END

     'ASCII': BEGIN
         asciifilename = DIALOG_PICKFILE(FILE = pname+'.dat',/WRITE,$
                                      /OVERWRITE_PROMPT)
         IF (STRLEN(asciifilename) NE 0) THEN BEGIN
             midplane_zfspect_plot, state.plot_log, state.krmax, $
               state.NFFT, state.plot_shear, ASCIIFILE=asciifilename
           
         ENDIF ELSE MESSAGE, 'text save cancelled', /INFO
     END

     'DONE': WIDGET_CONTROL, the_event.top, /DESTROY

  ENDCASE

  IF (uvalue NE 'DONE') THEN $
    WIDGET_CONTROL, the_event.top, SET_UVALUE = state, /NO_COPY

END ;midplane_zfspect_event

PRO midplane_zfspect_see, group=group
;
; C. Holland, UCSD
; v1.0: 5/20/2008
;
; Plots the spectrum of a n=0 fluctuation field at the outboard
; midplane, as a function of k_r and omega (plot_shear = 0).  
;
; Setting plot_shear = 1 plots the shear spectrum 
; |d^2 field/dr^2|^2 = k_r^4 |field|^2
;
; krmax specified max kr value to be plotted.  Routine assumes results
; are from a local sim with periodic radial BC's.
;
; NFFT_in specfies # pts used in frequency FFT, set less than
; (it2-it1+1)/2 to average multiple realizations
;

  common GLOBAL
  common MIDPLANE_DATA

  ;;---------------------------------------------------------
  ;; Return conditions
  if (n_pwr EQ 0) then return 
  if (n_n EQ 1) then return
  if xregistered('midplane_zfspect_see') then return
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
                    value='lin/log', $
                    uvalue='linlog')

  x = widget_button(row1, $
                    value='plot shear on/off', $
                    uvalue='shear_toggle')

  x = widget_button(row1, $
                    value='idx++', $
                    uvalue='idx_up')

  x = widget_button(row1, $
                    value='idx--', $
                    uvalue='idx_dn')

  x = cw_field(row1, /RETURN_EVENTS, TITLE = 'Max kr', $
               /FLOATING, VALUE = MAX(kt_rho), XSIZE=8, $
               UVALUE = 'MAXKR')

  x = cw_field(row1, /RETURN_EVENTS, TITLE = 'NFFT', $
               /FLOATING, VALUE = it2-it1+1, XSIZE=8, $
               UVALUE = 'NFFT')

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
  state = {winID:winID, krmax:MAX(kt_rho), NFFT:it2-it1+1, plot_log:0, $
          plot_shear:0}
  WIDGET_CONTROL, base, SET_UVALUE = state, /NO_COPY

  xmanager,'midplane_zfspect_see', $
           base, $
           event='midplane_zfspect_event',$
           group_leader=group

END ;midplane_zfspect_see
