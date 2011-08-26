; xunits[0]  : m_ref (kg)
; xunits[1]  : b_unit (Tesla)
; xunits[2]  : a (m)
; xunits[3]  : csD/a (1/s)
; xunits[4]  : csD (m/s)
; xunits[5]  : Te (keV)
; xunits[6]  : ne (10^19/m^3)
; xunits[7]  : rho_sD (m)
; xunits[8]  : chi_gBD (m^2/s)
; xunits[9]  : Gamma_gBD (0.624e22/m^2/s) = (MW/keV/m^2)
; xunits[10] : Q_gBD (MJ/m^2/s)           = (MW/m^2)
; xunits[11] : Pi_gBD (J/m^2)             = (Nm/m^2)
; xunits[12] : S_gBD (MJ/m^3/s)           = (MW/m^3)

pro read_units

  common GLOBAL
  common COLLISION_DATA

  file = 'out.gyro.units'

  openr,1,file,err=i_err
  if (i_err eq 0) then begin

     xunits = fltarr(13)
     for i=0,12 do begin
        readf,1,dummy
        xunits[i] = dummy
     endfor

     print_found,file,1
     exists_units = 1

  endif else begin

     print_found,file,0
     exists_units = 0

  endelse
  close,1

  return

end
