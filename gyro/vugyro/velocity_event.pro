pro velocity_event, velocity_obj

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_VELOCITY

  widget_control, velocity_obj.id, $
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
        counter_up,in_c,n_n-1,1 
        goto, plot_it
     end

     7: begin
        counter_dn,in_c,0,1  
        goto, plot_it
     end

     8: begin
        i_spec = i_spec+1
        if (i_spec ge n_kinetic) then i_spec = 0
        goto, plot_it
     end   

     9: begin
        i_moment = 1+i_moment
        if (i_moment ge 2) then i_moment = 0
        goto, plot_it
     end   

     10: begin
        if (n_field gt 1) then begin
           i_f = i_f+1
           if (i_f ge n_field) then i_f = 0
        endif
        goto, plot_it
     end

     11: begin
        plot_mode = 2
        goto, plot_it
     end

     12: widget_control, velocity_obj.top, /destroy
     
     20: begin
        i_plot = 0
        goto, plot_it
     end
     21: begin
        i_plot = 1
        goto, plot_it
     end
     22: begin
        i_plot = 2
        goto, plot_it
     end

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  smf_tag,title,pname,ytitle

  title_e = '!4e!3=E/T'
  title_l = '!4k=l!3/E'

  pname = 'diff_v'+pname

  plot_def_new,pname+strtrim(string(in_c),2)

  title = '!3'+title+' !3flux at !3k!4!dh!nq!d!3s!n!4 = '+$
    strmid(strtrim(string(kt_rho[in_c]),2),0,5)

  z = fltarr(n_energy,n_lambda)
  z[*,*] = f_nek[*,*,i_spec,i_f,i_moment,in_c,t_c]
  ;;z = z/max(z)

  case (i_plot) of 

     0: begin
        
        loadct,CTab,/SILENT

        contour,z,energy,lambda,$
          nlevels=nlevels,/fill,$
          xtitle=title_e,$
          ytitle=title_l,$
          title=title

        oplot,[0,energy[n_energy-1]],lambda_tp*[1,1],linestyle=2

        for ie=0,n_energy-1 do begin
           oplot,0.0*lambda+energy[ie],lambda,psym=8
        endfor

        set_line_colors

     end
     
     1: begin

        z0 = fltarr(n_energy)
        for i=0,n_lambda-1 do begin
           z0[*] = z0[*]+z[*,i]
        endfor

        plot,[0],[0],$
          /nodata,$
          title=title,$
          xstyle=1,$
          xrange=[0,energy[n_energy-1]],$
          xtitle=title_e,$
          ystyle=1,$
          yrange=max(abs(z0))*[-1,1],$
          color=line

        oplot,energy,z0,color=color_vec[1]
        oplot,energy,z0,color=line,psym=8
        oplot,[0,energy[n_energy-1]],[0,0],linestyle=1

     end

     2: begin

        z0 = fltarr(n_lambda)
        for i=0,n_energy-1 do begin
           z0[*] = z0[*]+z[i,*]
        endfor

        plot,[0],[0],/nodata,$
          title=title,$
          xstyle=1,$
          xminor=0,$
          xrange=[0,lambda[n_lambda-1]],$
          xtitle=title_l,$
          ystyle=1,$
          yminor=0,$
          yrange=max(abs(z0))*[-1,1],$
          color=line
        
        oplot,lambda[0:n_pass-1],z0[0:n_pass-1],$
          color=color_vec[1]

        oplot,lambda[0:n_pass-1],z0[0:n_pass-1],$
          color=line,psym=8

        oplot,lambda[n_pass:n_lambda-1],z0[n_pass:n_lambda-1],$
          color=color_vec[2]

        oplot,lambda[n_pass:n_lambda-1],z0[n_pass:n_lambda-1],$
          color=line,psym=8

        oplot,[0,lambda[n_lambda-1]],[0,0],linestyle=1

        oplot,lambda_tp*[1,1],max(abs(z0))*[-1,1],linestyle=2

     end

  endcase

  plot_finish

end
