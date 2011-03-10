pro gamread

; Procedure to read growth rate profile

id = ncdf_open('xptor_out.nc')
ncdf_varget,id,ncdf_varid(id,'rho'),x
ncdf_varget,id,ncdf_varid(id,'anrate_m'),gamma
ncdf_varget,id,ncdf_varid(id,'egamma_m'),egamma
;ncdf_varget,id,ncdf_varid(id,'anfreq_m'),omega

print,'         rho           gamma           egamma           omega'
for m=1,51 do print, x(m), gamma(m),egamma(m)

end
