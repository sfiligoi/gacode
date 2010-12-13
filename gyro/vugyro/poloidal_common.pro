pro poloidal_common

  common POLOIDAL_DATA, $
    mUt, $
    n_mUt, $
    s_mUt, $
    v_mUt, $
    a, $
    th, $
    qt, $
    i_field,$
    active_poloidal,$
    one_contour

  ;;--------------------------------------------
  ;; Theta gridding
  ;;
  n_mUt = 8
  ;;
  s_mUt = ['x1','x2','x4','x8','x16','x24','x32','x64']
  v_mUt = [1,2,4,8,16,24,32,64]
  ;;--------------------------------------------

  ;; Initializations

  mUt = 1

  ;; field selector
  i_field = 1

  active_poloidal = 1
  one_contour = 0

end
