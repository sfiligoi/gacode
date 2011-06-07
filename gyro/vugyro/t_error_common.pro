pro t_error_common

  common T_ERROR_DATA, $
    t_error_min, $
    t_error_min_fac, $
    t_error_wid, $
    exists_t_error, $
    t_error

  ;; Initialize

  t_error_min = 1.e-8
  t_error_min_fac = 1.

end
