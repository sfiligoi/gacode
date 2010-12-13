pro geometry_event, geometry_obj

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_GEOMETRY
  common PROFILE_SIM_DATA

  widget_control, geometry_obj.id, $
    get_uvalue=uvalue

  wset,widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  if uvalue lt 30 then begin

     case (uvalue) of

        1: begin
           counter_up,i_c,n_r-1,1 
           goto, plot_it
        end   

        2: begin
           counter_dn,i_c,0,1 
           goto, plot_it
        end

        3: begin
           plot_mode = 2
           goto, plot_it
        end

        4: begin
           plot_export = 1
           goto, plot_it
        end

        5: widget_control, geometry_obj.top, /destroy

     endcase

     return

  endif else begin

     i_geometry = uvalue-30

  endelse

  plot_it:
  
  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  pname = geometry_label[i_geometry]

  plot_def_new,pname

  xtitle = '!4h/p!3'

  y = fltarr(n_fine)
  y0 = fltarr(n_fine)

  x = -1.0+2.0*findgen(n_fine)/n_fine

  case (i_geometry) of

     0: y0 = -q_i(i_c)*(x*!pi) ; nu
     1: y0 = sin(x*!pi) ; gsin
     2: y0 = cos(x*!pi) ; gcos1
     3: y0 = 0.0*x ; gcos2
     4: y0 = sin(x*!pi) ; usin
     5: y0 = cos(x*!pi) ; ucos
     6: y0 = 1.0/(1.0+cos(x*!pi)/aspect_s(i_c)) ; b
     7: y0 = 1.0+0.0*x ; g_theta
     8: y0 = 1.0+0.0*x ; grad_r
     9: y0 = 1.0+0.0*x ; g_q
    10: y0 = shat_s(i_c)*(x*!pi) ; captheta

  endcase

  y = fltarr(n_fine)
  y = geometry[i_geometry,*,i_c]

  dy = max(y)-min(y)
  if (dy lt 0.5) then dy = 0.5
  ymax = max(y)+0.1*dy
  ymin = min(y)-0.1*dy
 
  plot,x,y0,linestyle=1,$
    ytitle=pname,$
    xstyle=1,$
    xminor=0,$
    xrange=[-1,1],$
    ystyle=1,$
    yrange=[ymin,ymax],$
    xtitle=xtitle,$
    color=line,$
    charsize=csize

  oplot,x,y,color=color_vec[1]

  if (plot_export eq 1) then begin

     openw,1,pname+'.idlout'
     for j=0,n_fine-1 do begin
        printf,1,x[j],y[j]
     endfor
     close,1

     print,'Exported data to ',pname+'.idlout'
     plot_export = 0

  endif

  plot_finish

  return

end
