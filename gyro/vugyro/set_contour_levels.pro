pro set_contour_levels

  common GLOBAL
  common PLOT_VARIABLES
  common POLOIDAL_DATA

  make_fine_grid

  if c_table_max lt 0.0 then begin

    extract_contour,n_time/2
    clevels = min(a)+findgen(nlevels)*(max(a)-min(a))/(nlevels-1.0)

  endif else begin

    clevels = c_table_min+$
        findgen(nlevels)*(c_table_max-c_table_min)/(nlevels-1.0)

  endelse

end
