pro color_table_setup

  common GLOBAL

  ;--------------------------------------------
  ; Contour table and levels 
  ;
  n_nLevels = 8
  ;
  s_nLevels = ['4','8','16','24','32','64','96','128']
  v_nLevels = [4,8,16,24,32,64,96,128]
  ;--------------------------------------------

  ;--------------------------------------------
  ; Color tables
  ;
  n_CTab = 37
  ;
  s_CTab = ['0 B-W linear',$
            '1 Blue/White',$
            '2 Green/Red/Blue/White',$
            '3 Red Temp',$
            '4 Blue/Green/Red/Yellow',$
            '5 Gamma-II',$
            '6 Prism',$
            '7 Red-Purple',$
            '8 Green/White linear',$
            '9 Green/White exp',$
            '10 Green-Pink',$
            '11 Blue-Red', $
            '12 16 level',$
            '13 Rainbow',$
            '14 Steps',$
            '15 Stern Special',$
            '16 Haze',$
            '17 Blue-Pastel-Red',$
            '18 Pastels',$
            '19 HSL 1',$
            '20 HSL 2',$
            '21 HSV 1',$
            '22 HSV 2',$
            '23 Purple-Red',$
            '24 Beach',$
            '25 Mac Style',$
            '26 Eos A',$
            '27 Eos B',$
            '28 Hardcandy',$
            '29 Nature',$
            '30 Ocean',$
            '31 Peppermint',$
            '32 Plasma',$
            '33 Blue-Red',$
            '34 Rainbow',$
            '35 Blue waves',$
            '36 Volcano']

  v_CTab = indgen(n_CTab)
  ;--------------------------------------------

end 
