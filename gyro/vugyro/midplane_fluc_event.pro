pro midplane_fluc_event, midplane_fluc

  common GLOBAL
  common MIDPLANE_DATA
  common POLOIDAL_DATA

  widget_control, midplane_fluc.id, $
    get_uvalue=uvalue

  wset, midplane_fluc_wid

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  ;; Refinement
  if (uvalue ge 40) then begin
     mUt = v_mUt[uvalue-40]
     goto, plot_it
  endif

  case (uvalue) of 

     1: begin
        goto, plot_it
     end   

     2: begin
        counter_up,i_pwr,n_pwr-1,1
        goto, plot_it
     end   

     3: begin
        counter_dn,i_pwr,0,1
        goto, plot_it
     end   

     4: begin
        counter_up,t_c,n_time1,1
        goto, plot_it
     end

     5: begin
        counter_dn,t_c,0,1
        goto, plot_it
     end

     6: begin
        counter_up,t_c,n_time1,10 
        goto, plot_it
     end

     7: begin 
        counter_dn,t_c,0,10
        goto, plot_it
     end

     8: begin
        plot_mode = 2
        goto, plot_it
     end

     9: begin
        movie_flag = 1
        goto, plot_it
     end

     10: widget_control, midplane_fluc.top, /destroy

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------


  if movie_flag eq 0 then begin  

     midplane_fluc_plot

  endif else begin

     t_0 = t_c

     for t_c=it1,it2 do begin

        midplane_fluc_plot

        ind  = strmid(string(long(t_c)+10000,format='(i5)'),1)
        file = strcompress('frame'+ind+'.jpg',/remove_all)

        tmp_image = tvrd(true=3,/order) 
        write_jpeg,file,tmp_image,true=3,/order
        print,'Wrote '+file

     endfor
     
     t_c = t_0

     close,1
     movie_flag = 0

  endelse

  return
  
end
