; DEVSETUP                    C.M. Greenfield            May 5, 1999
;                             Last modification    November 19, 2001
; NAME:
;       DEVSETUP
;
; PURPOSE:
;       Set up the plot device (PostScript or X)
;
; CATEGORY:
;       IDL utility
;
; CALLING SEQUENCE:
;       DS_STRUCT=devsetup([ORIENTATION=ORIENTATION][,CTABLE=CTABLE]
;                          [,WINDOW=WINDOW][,NONEW=NONEW]
;                          [,WTITLE=WTITLE],[,GEOMETRY=GEOMETRY])
;
; INPUTS:
;       NONE
;
; KEYWORD PARAMETERS:
;       ORIENTATION
;               Sets plot orientation in PS mode. The possibilities are:
;                       'landscape' (the default)
;                       'portrait'
;                       'square'    (square area on portrait page)
;               No effect on X-window device
;       CTABLE  Load the color table specified, with color 0 as
;               foreground and last color as background. Special
;               case: CTABLE='CMG' forces load of my personal color
;               table (CMGCOLORS). CTABLE='CMGSMALL' loads my small
;               color table.
;       WINDOW  Window ID to use. Only valid for X-windows device
;       NONEW   Do not create new windows.
;       WTITLE  Window title
;       GEOMETRY[xsize,ysize,x0,y0]
;
; RETURNS:
;       DS_STRUCT
;               Structure containing previous graphics state and
;               strings defining symbols for use in plots.
;
; NOTES:
;
;
; OTHER REQUIRED ROUTINES:
;
;  NONE
;
;-----------------------------------------------------------------------------
function devsetup,orientation=orientation,ctable=ctable,window=window, $
                  nonew=nonew,wtitle=wtitle,geometry=geometry
common XXX_DEVSETUP_INFO_XXX,oldwtitle
if (n_elements(oldwtitle) eq 0) then oldwtitle='$NULL$'
psave=!p & xsave=!x & ysave=!y
tvlct,rsave,gsave,bsave,/get
case !d.name of
    'X':begin
        device,decomposed=0
        !p.font=-1
        Lrho='!4q!X'    & Urho='!4Q!X'
        Lchi='!4v!X'    & Uchi='!4V!X'
        Lomega='!4x!X'  & Uomega='!4X!X'
        Lgamma='!4c!X'  & Ugamma='!4C!X'
        Lphi='!4u!X'    & Uphi='!4U!X'
        Lpsi='!4w!X'    & Upsi='!4W!X'
        Lalpha='!4a!X'  & Ualpha='!4A!X'
        Lbeta='!4b!X'   & Ubeta='!4B!X'
        Lkappa='!4j!X'  & Ukappa='!4J!X'
        Ltheta='!4h!X'  & Utheta='!4H!X'
        Ldelta='!4d!X'  & Udelta='!4D!X'
        times='!9X!X'   & div='!9/!X'
        del='!9G!X'     & pder='!9D!X'
        ntilde=string(241b)
        if (n_elements(window) ne 0) then wset,window
        device,set_character_size=[6,10]
        !p.thick=2
    end
    'PS':begin
        device,/color,/helvetica,/narrow,/bold,isolatin1=1
        if (n_elements(orientation) eq 1) then case strupcase(orientation) of
            'LANDSCAPE':device,/landscape
            'PORTRAIT':device,/portrait,ysize=24.9,yoffset=1.5
            'SQUARE':device,/portrait,xsize=20,ysize=20, $
              xoffset=0.795,yoffset=3.97
            default:device,/landscape
        endcase else device,/landscape
        !p.font=0
        Lrho='!Mr!X'   & Urho='!MR!X'
        Lchi='!Mc!X'   & Uchi='!MC!X'
        Lomega='!Mw!X' & Uomega='!MW!X'
        Lgamma='!Mg!X' & Ugamma='!MG!X'
        Lphi='!Mf!X'   & Uphi='!MF!X'
        Lpsi='!My!X'   & Upsi='!MY!X'
        Lalpha='!Ma!X' & Ualpha='!MA!X'
        Lbeta='!Mb!X'   & Ubeta='!MB!X'
        Lkappa='!Mk!X'  & Ukappa='!MK!X'
        Ltheta='!Mq!X'  & Utheta='!MQ!X'
        Ldelta='!Md!X'  & Udelta='!MD!X'
        times='!M'+string(180b)+'!X'
        div='!M'+string(184b)+'!X'
        del='!M'+string(209b)+'!X'
        pder='!M'+string(182b)+'!X'
        ntilde='!X'+string(241b)
        !p.thick=4
    end
endcase
ds_struct={state:{p:psave,x:xsave,y:ysave,r:rsave,g:gsave,b:bsave}, $
           lc:{rho:Lrho,chi:Lchi,omega:Lomega,gamma:Lgamma, $
               psi:Lpsi,phi:Lphi, $
               alpha:Lalpha,beta:Lbeta,kappa:Lkappa,theta:Ltheta, $
               delta:Ldelta}, $
           uc:{rho:Urho,chi:Uchi,omega:Uomega,gamma:Ugamma, $
               psi:Upsi,phi:Uphi, $
               alpha:Ualpha,beta:Ubeta,kappa:Ukappa,theta:Utheta, $
               delta:Udelta}, $
           sym:{times:times,div:div,del:del,pder:pder,ntilde:ntilde}}
if keyword_set(ctable) then begin
    if ((reverse(size(ctable)))[1] eq 7) then begin
        case strupcase(ctable) of
            'CMG':cnames=cmgcolors()
            'CMGSMALL':cnames=cmgcolors(/small)
            ELSE:cnames=cmgcolors()
        endcase
        !p.color=cnames.black
        !p.background=cnames.white
        ds_struct=create_struct(ds_struct,'cnames',cnames)
    endif else begin
        loadct,ctable,ncolors=!d.table_size,/silent
        !p.color=0
        !p.background=!d.table_size-1
        ds_struct=create_struct(ds_struct,'ctable',ctable)
    endelse
endif
;Create window if X-windows and either title or geometry has changed
if (!d.name eq 'X') then begin
    if (n_elements(wtitle) eq 0) then begin
        help,calls=a
        wtitle=(strsplit(a[1],' ',/extract))[0]
    endif
    if keyword_set(nonew) then newwin=0b $
    else begin
        device,get_window_position=pos
        newwin=(wtitle ne oldwtitle)
        case n_elements(geometry) of
            2:BEGIN
                newwin=newwin or $
                  (geometry[0] ne !d.x_size) or (geometry[1] ne !d.y_size)
                geometry=[geometry,pos]
                usegeom=1b
            END
            4:BEGIN
                newwin=newwin or $
                  (geometry[0] ne !d.x_size) or $
                  (geometry[1] ne !d.y_ysize) or $
                  (geometry[2] ne pos[0]) or (geometry[3] ne pos[1])
                usegeom=1b
            END
            else:usegeom=0b
        endcase
    endelse
    if (newwin) then begin
        if (usegeom) then window,title=wtitle,xsize=geometry[0], $
          ysize=geometry[1],xpos=geometry[2],ypos=geometry[3] $
        else window,title=wtitle
        oldwtitle=wtitle
    endif
endif
return,ds_struct
end
