pro balloon_common

  common PRIVATE_BALLOON,$
    n_balloon, $
    balloon_tag, $
    balloon_plot, $
    dy_balloon, $
    balloon_index, $ 
    balloon_l0, $
    balloon_norm, $ 
    balloon_m, $
    balloon_np, $
    balloon_next, $
    widget

  ;; Initialize

  dy_balloon    = 1.0
  balloon_index = 0
  balloon_l0    = 0
  balloon_norm  = 0
  balloon_m     = 0
  balloon_np    = 0
  balloon_next  = 0

end
