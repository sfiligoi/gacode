pro datread

; Procedure to read profiles vs time

id = ncdf_open('xptor_out.nc')
ncdf_varget,id,ncdf_varid(id,'nsteps'),n
ncdf_varget,id,ncdf_varid(id,'nrho'),nr
ncdf_varget,id,ncdf_varid(id,'rho'),x
ncdf_varget,id,ncdf_varid(id,'time'),t
;ncdf_varget,id,ncdf_varid(id,'ti_t'),ti
;ncdf_varget,id,ncdf_varid(id,'te_t'),te
;ncdf_varget,id,ncdf_varid(id,'vphi_t'),vphi
ncdf_varget,id,ncdf_varid(id,'ti_m'),ti
ncdf_varget,id,ncdf_varid(id,'te_m'),te
ncdf_varget,id,ncdf_varid(id,'vphi_m'),vphi
ncdf_varget,id,ncdf_varid(id,'ti_exp'),tix
ncdf_varget,id,ncdf_varid(id,'te_exp'),tex
; Here, ti=ti(x,t) where t=no. timesteps-1
print,n
print,'         rho           Te-Model         Ti-Model'
for m=1,nr do print, x(m), te(m), ti(m)
print,'         rho           Te-Data          Ti-Data'
for m=1,nr do print, x(m), tex(m), tix(m)

end
