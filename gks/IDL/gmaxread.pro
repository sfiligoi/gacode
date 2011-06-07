pro gmaxread

; Procedure to read profiles vs time

id = ncdf_open('gksout.nc')
ncdf_varget,id,ncdf_varid(id,'rho'),x
ncdf_varget,id,ncdf_varid(id,'anrate_m'),gamma
ncdf_varget,id,ncdf_varid(id,'dnrate_m'),dgamma
ncdf_varget,id,ncdf_varid(id,'anfreq_m'),omega
ncdf_varget,id,ncdf_varid(id,'ky_m'),ky
print,'         rho           gamma_max        dgamma          omega            ky'
for m=1,50 do print, x(m),gamma(m),dgamma(m),omega(m),ky(m)

end
