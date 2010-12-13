pro balloon_label_update,bm,t_str,n_str

  common GLOBAL

  ;-------------------------------------------------
  ; widget control for label functions
  ;
  widget_control, bm.top, $
                  get_uvalue=state, $
                  /no_copy

  widget_control, state.t_label, set_value=t_str      
  widget_control, state.n_label, set_value=n_str      

  widget_control, bm.top, $
                  set_uvalue=state, $
                  /no_copy
  ;-------------------------------------------------

end
