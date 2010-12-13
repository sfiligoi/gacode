; $Id: xdisplayfile.pro,v 1.1.1.1 2010/03/04 00:31:39 candy Exp $
;
; Copyright (c) 1991-2003, Research Systems, Inc.  All rights reserved.
;	Unauthorized reproduction prohibited.

PRO XdispFile_write, wText, filename

  COMPILE_OPT hidden

  WIDGET_CONTROL, /HOURGLASS
  OPENW, unit, FILENAME, /GET_LUN, ERROR=i		;open the file and then
  if i lt 0 then begin		;OK?
	a = [ !error_state.msg, filename + ' could not be opened for writing.']  ;No
	void = DIALOG_MESSAGE(a, /ERROR, DIALOG_PARENT=wText)
  endif else begin
	WIDGET_CONTROL, wText, get_value=txtArray
	ON_IOERROR, done_writing
	; print out each line separately in order to get desired line breaks
	for j = 0, N_ELEMENTS(txtArray)-1 do PRINTF, unit, txtArray[j]
done_writing:
	ON_IOERROR, null
	FREE_LUN, unit				;free the file unit.
  endelse
END

PRO XDispFile_evt, event

  COMPILE_OPT hidden

  WIDGET_CONTROL, event.top, GET_UVALUE=state
  CASE TAG_NAMES(event, /STRUCTURE_NAME) OF
    'WIDGET_BASE': BEGIN
      WIDGET_CONTROL, event.top, /MAP, /UPDATE
      RETURN
      END
    'WIDGET_KILL_REQUEST': retval = "EXIT"
    ELSE: WIDGET_CONTROL, event.id, GET_UVALUE = retval
  ENDCASE

  CASE retval OF
  	"SAVE": BEGIN
               if (LMGR(/DEMO)) then begin
                  tmp = DIALOG_MESSAGE( /ERROR, $
                        'Save: Feature disabled for demo mode.')
                  return
                endif
		IF (STRLEN(state.filename) EQ 0) THEN BEGIN
			state.filename = DIALOG_PICKFILE(/WRITE)
		ENDIF
		IF (STRLEN(state.filename) GT 0) THEN BEGIN
			XdispFile_write, state.filetext, state.filename
			WIDGET_CONTROL, event.top, SET_UVALUE=state
			IF state.notitle THEN WIDGET_CONTROL, event.top, $
				TLB_SET_TITLE=state.filename
		ENDIF

		RETURN
  	END
	"SAVE_AS": BEGIN
               if (LMGR(/DEMO)) then begin
                  tmp = DIALOG_MESSAGE( /ERROR, $
                        'Save As: Feature disabled for demo mode.')
                  return
                endif
		state.filename = DIALOG_PICKFILE(/WRITE)
		IF (STRLEN(state.filename) GT 0) THEN BEGIN
			XdispFile_write, state.filetext, state.filename
			WIDGET_CONTROL, event.top, SET_UVALUE=state
			IF state.notitle THEN WIDGET_CONTROL, event.top, $
				TLB_SET_TITLE=state.filename
		ENDIF
		RETURN
	END
	"EXIT": BEGIN
		WIDGET_CONTROL, event.top, /DESTROY
		IF (WIDGET_INFO(state.ourGroup, /VALID)) THEN $
			WIDGET_CONTROL, state.ourGroup, /DESTROY
	END
	ELSE:
  ENDCASE
END


PRO XDisplayFile, FILENAME, TITLE = TITLE, GROUP = GROUP, WIDTH = WIDTH, $
                  HEIGHT = HEIGHT, TEXT = TEXT, FONT = font, $
                  DONE_BUTTON=done_button, MODAL=MODAL, $
                  EDITABLE=editable, $
                  WTEXT=filetext, BLOCK=block
;+
; NAME:
;	XDISPLAYFILE
;
; PURPOSE:
;	Display an ASCII text file using widgets and the widget manager.
;
; CATEGORY:
;	Widgets.
;
; CALLING SEQUENCE:
;	XDISPLAYFILE, Filename
;
; INPUTS:
;     Filename:	A scalar string that contains the filename of the file
;		to display.  The filename can include a path to that file.
;
; KEYWORD PARAMETERS:
;	BLOCK:  Set this keyword to have XMANAGER block when this
;		application is registered.  By default the Xmanager
;               keyword NO_BLOCK is set to 1 to provide access to the
;               command line if active command 	line processing is available.
;               Note that setting BLOCK for this application will cause
;		all widget applications to block, not only this
;		application.  For more information see the NO_BLOCK keyword
;		to XMANAGER.
;
;	DONE_BUTTON: the text to use for the Done button.  If omitted,
;		the text "Done with <filename>" is used.
;
;	EDITABLE: Set this keyword to allow modifications to the text
;		displayed in XDISPLAYFILE.  Setting this keyword also
;		adds a "Save" button in addition to the Done button.
;
;	FONT:   The name of the font to use.  If omitted use the default
;		font.
;	GROUP:	The widget ID of the group leader of the widget.  If this
;		keyword is specified, the death of the group leader results in
;		the death of XDISPLAYFILE.
;
;	HEIGHT:	The number of text lines that the widget should display at one
;		time.  If this keyword is not specified, 24 lines is the
;		default.
;
;	TEXT:	A string or string array to be displayed in the widget
;		instead of the contents of a file.  This keyword supercedes
;		the FILENAME input parameter.
;
;	TITLE:	A string to use as the widget title rather than the file name
;		or "XDisplayFile".
;
;	WIDTH:	The number of characters wide the widget should be.  If this
;		keyword is not specified, 80 characters is the default.
;
;	WTEXT:	Output parameter, the id of the text widget.  This allows
;		setting text selections and cursor positions programmatically.
;
; OUTPUTS:
;	No explicit outputs.  A file viewing widget is created.
;
; SIDE EFFECTS:
;	Triggers the XMANAGER if it is not already in use.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	Open a file and create a widget to display its contents.
;
; MODIFICATION HISTORY:
;	Written By Steve Richards, December 1990
;	Graceful error recovery, DMS, Feb, 1992.
;       12 Jan. 1994  - KDB
;               If file was empty, program would crash. Fixed.
;       4 Oct. 1994     MLR Fixed bug if /TEXT was present and /TITLE was not.
;	2 jan 1997	DMS Added DONE_BUTTON keyword, made Done
;			button align on left, removed padding.
;-

IF(NOT(KEYWORD_SET(EDITABLE))) THEN editable = 0          ;use the defaults if
IF(NOT(KEYWORD_SET(HEIGHT))) THEN HEIGHT = 24		;the keywords were not
IF(NOT(KEYWORD_SET(WIDTH))) THEN WIDTH = 80		;passed in
IF N_ELEMENTS(block) EQ 0 THEN block=0
noTitle = N_ELEMENTS(title) EQ 0

IF(NOT(KEYWORD_SET(TEXT))) THEN BEGIN
  IF noTitle THEN TITLE = FILENAME
  OPENR, unit, FILENAME, /GET_LUN, ERROR=i		;open the file and then
  if i lt 0 then begin		;OK?
	a = [ !error_state.msg, ' Unable to display ' + filename]  ;No
  endif else begin
	  maxlines = 10000		;Maximum # of lines in file
	  a = strarr(maxlines)
	  on_ioerror, done_reading
	  readf, unit, a
done_reading: s = fstat(unit)		;Get # of lines actually read
	  a = a[0: (s.transfer_count-1) > 0]
	  on_ioerror, null
	  FREE_LUN, unit				;free the file unit.
  endelse
ENDIF ELSE BEGIN
    IF(N_ELEMENTS(FILENAME) EQ 0) THEN FILENAME=''
    IF noTitle THEN TITLE = 'XDisplayFile'
    a = TEXT
ENDELSE

ourGroup = 0L
if KEYWORD_SET(MODAL) then begin
    if N_ELEMENTS(GROUP) GT 0 then begin
        filebase = WIDGET_BASE(TITLE = TITLE, $
        			/TLB_KILL_REQUEST_EVENTS, TLB_FRAME_ATTR=1, $
                               /BASE_ALIGN_LEFT, /COLUMN, $
                               /MODAL, GROUP_LEADER = GROUP)
    endif else begin
	; modal requires a group leader
        ourGroup = WIDGET_BASE()
        filebase = WIDGET_BASE(TITLE = TITLE, $
        			/TLB_KILL_REQUEST_EVENTS, TLB_FRAME_ATTR=1, $
                               /BASE_ALIGN_LEFT, /COLUMN, $
                               /MODAL, GROUP_LEADER = ourGroup)
    endelse
    menu_bar = filebase
endif else begin
    filebase = WIDGET_BASE(TITLE = TITLE, $
        			/TLB_KILL_REQUEST_EVENTS, /TLB_SIZE_EVENTS, $
                           /BASE_ALIGN_LEFT, /COLUMN, MBAR=menu_bar)
endelse


extra = ''
IF (menu_bar NE filebase) THEN BEGIN
	IF (!VERSION.OS_FAMILY EQ 'Windows') THEN extra = '&'
	menu_bar = WIDGET_BUTTON(menu_bar, VALUE=extra+'File', /MENU)
ENDIF ELSE $
	menu_bar = WIDGET_BASE(filebase, /ROW)


IF (editable) THEN BEGIN
	; add 'Save', 'Save as...' buttons here
	saveButton = WIDGET_BUTTON(menu_bar, $
		VALUE = extra+'Save', UVALUE = "SAVE")
	saveAsButton = WIDGET_BUTTON(menu_bar, $
		VALUE = 'Save '+extra+'As...', UVALUE = "SAVE_AS")

ENDIF

; Default Done button name:
if n_elements(done_button) eq 0 then done_button = "Done with " + TITLE
filequit = WIDGET_BUTTON(menu_bar, $			;create a Done Button
		SEPARATOR=editable, $
		VALUE = extra+done_button, UVALUE = "EXIT")

IF n_elements(font) gt 0 then $
 filetext = WIDGET_TEXT(filebase, $			;create a text widget
		XSIZE = WIDTH, $			;to display the file's
		YSIZE = HEIGHT, $			;contents
		EDITABLE = editable, $
		UVALUE='TEXT', $
		/SCROLL, FONT = font, $
		VALUE = a) $
ELSE filetext = WIDGET_TEXT(filebase, $			;create a text widget
		XSIZE = WIDTH, $			;to display the file's
		YSIZE = HEIGHT, $			;contents
		EDITABLE = editable, $
		UVALUE='TEXT', $
		/SCROLL, $
		VALUE = a)


WIDGET_CONTROL, filebase, /REALIZE			;instantiate the widget

state={ourGroup:ourGroup, $
       filename: filename, $
       filetext:filetext, $
       notitle:noTitle}
WIDGET_CONTROL, filebase, SET_UVALUE = state


Xmanager, "XDisplayFile", $				;register it with the
		filebase, $				;widget manager
		GROUP_LEADER = GROUP, $
		EVENT_HANDLER = "XDispFile_evt", NO_BLOCK=(NOT(FLOAT(block)))

END  ;--------------------- procedure XDisplayFile ----------------------------



