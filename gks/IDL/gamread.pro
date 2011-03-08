pro gamread

; Procedure to read gamma,omega profile

id = ncdf_open('gksout.nc')
ncdf_varget,id,ncdf_varid(id,'shot'),shot
ncdf_varget,id,ncdf_varid(id,'rho_k'),x
ncdf_varget,id,ncdf_varid(id,'rho'),ky
ncdf_varget,id,ncdf_varid(id,'gamma_k'),gamma
ncdf_varget,id,ncdf_varid(id,'dgamma_k'),dgamma
ncdf_varget,id,ncdf_varid(id,'freq_k'),freq
ncdf_varget,id,ncdf_varid(id,'gamma_mks'),gammamks
ncdf_varget,id,ncdf_varid(id,'ky_mks'),kymks

print,shot
print,'         rho            ky              gamma         dgamma            freq           gamma_mks       ky_mks'
for m=1,21 do print, x, ky(m), gamma(m), dgamma(m), freq(m),gammamks(m),kymks(m)

end
