pro itgcrit

; Procedure to the critical ITG profile

id = ncdf_open('gksout.nc')
ncdf_varget,id,ncdf_varid(id,'shot'),shot
ncdf_varget,id,ncdf_varid(id,'jout_m'),jout
ncdf_varget,id,ncdf_varid(id,'rho'),x
ncdf_varget,id,ncdf_varid(id,'ky_m'),ky
ncdf_varget,id,ncdf_varid(id,'zpti_exp'),zptix
ncdf_varget,id,ncdf_varid(id,'zpti_crit_m'),zptim

print,shot
print,'         rho            zpti_exp       zpti_crit          ky'
for m=1,jout do print, x(m), zptix(m), zptim(m), ky(m)

end
