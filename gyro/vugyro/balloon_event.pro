pro balloon_event, balloon

  common GLOBAL
  common PRIVATE_BALLOON

  widget_control, balloon.id, $
  get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: goto, plot_it

     2: begin
        counter_up,t_c,n_time1,1
        goto, plot_it
     end

     3: begin
        counter_dn,t_c,0,1
        goto, plot_it
     end

     4: begin
        counter_up,t_c,n_time1,10 
        goto, plot_it
     end

     5: begin 
        counter_dn,t_c,0,10
        goto, plot_it
     end

     6: begin
        balloon_index = balloon_index+1
        if (balloon_index ge n_balloon) then balloon_index = 0 
        goto, plot_it
     end

     7: begin
        balloon_l0 = balloon_l0+1
        if (balloon_l0 ge balloon_m) then balloon_l0 = 0 
        goto, plot_it
     end
   
     8: begin
        counter_up,balloon_norm,n_balloon,1 
        goto, plot_it
     end

     9: begin
        counter_dn,balloon_norm,0,1 
        goto, plot_it
     end

     10: begin
        counter_up,i_p,n_r/2,1 
        goto, plot_it
     end

     11: begin
        counter_dn,i_p,0,1 
        goto, plot_it
     end

     12: begin
        dy_balloon = dy_balloon*2^0.5
        goto, plot_it
     end

     13: begin
        dy_balloon = dy_balloon/2^0.5
        goto, plot_it
     end

     14: begin
        plot_mode = 2
        goto, plot_it
     end

     15: begin
        plot_export = 1
        goto, plot_it
     end

     16: begin
        plot_export = 2
        goto, plot_it
     end

     17: widget_control, balloon.top, /destroy
     
  endcase

  return

plot_it:

;;-------------------------------------------------------
;; PLOTTING
;;-------------------------------------------------------

  if (plot_export eq 2) then begin

     ;; MPEG GENERATION

     plot_export = 0

     t_store = t_c

     openw,1,'incl'

     for tt=0,n_time1 do begin
        
        if (t[tt] ge t_min) and (t[tt] le t_max) then begin

           t_str = 't = '+strtrim(string(t[t_c]),2)
           t_c = tt

           balloon_label_update,balloon,t_str
           balloon_plot,balloon_index

           file = balloon_tag[balloon_index]+strtrim(string(tt),2)+'.jpg'
           tmp_image = tvrd(true=3,/order) 
           write_jpeg,file,tmp_image,true=3,/order
           print,'Wrote '+file
           printf,1,file

        endif
        
     endfor

     close,1

     t_c = t_store

  endif else begin

     ;; NORMAL PLOT

     t_str = 't = '+strtrim(string(t[t_c]),2)
     if balloon_norm eq 0 then begin
        n_str = 'normalization: standard'
     endif else begin
        n_str = 'normalization: '+balloon_tag[balloon_norm-1]
     endelse

     balloon_label_update,balloon,t_str,n_str
     balloon_plot,balloon_index

  endelse

  return


end
