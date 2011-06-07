pro chiplt,names

; Procedure to plot chi-profiles

id = ncdf_open('xptor_out.nc')
ncdf_varget,id,ncdf_varid(id,'nrho'),nr
ncdf_varget,id,ncdf_varid(id,'rho'),x
ncdf_varget,id,ncdf_varid(id,'chie_m'),chiem
ncdf_varget,id,ncdf_varid(id,'chii_m'),chiim
ncdf_varget,id,ncdf_varid(id,'chii_neo_m'),chiineo
ncdf_varget,id,ncdf_varid(id,'chii_t'),chiit
; printf,-1,f

names=strarr(8)
names[0]='chii-GLF23'
names[1]='chie-GLF23'
names[2]='DIII-D #104276'
names[3]='t=5.71 secs'
names[4]='GLF23 v1.61 model'
names[5]='100 pt grid, BC@rho=0.90'
names[6]='Case DV #11, Te,Ti,vpol evolved'
names[7]='chii-neo'

set_plot,'x'
plotwhite
!p.charsize=1.25
a=findgen(17) * (!pi*2/16.)
usersym,cos(a),sin(a),/fill
dev=devsetup(ctable='cmg')
plot, x(*),chiem(*),title='Diffusivity vs rho',xtitle='rho',ytitle='chi (m^2/s)',line=0,thick=2,yrange=[0,6],ystyle=1,color=dev.cnames.red
usersym,cos(a),sin(a)
oplot, x(*),chiim(*),line=0,thick=2,color=dev.cnames.blue
oplot, x(*),chiineo(*),line=0,thick=2,color=dev.cnames.green
oplot, x(*),chiit(1,*),line=1,thick=2,color=dev.cnames.blue
xyouts, .70,  5.65, names(0),charsize=1.25,color=dev.cnames.blue
xyouts, .70,  5.35, names(1),charsize=1.25,color=dev.cnames.red
xyouts, .70,  5.05, names(7),charsize=1.25,color=dev.cnames.green
xyouts, .04,  5.65, names(2),charsize=1.25,color=dev.cnames.black
xyouts, .04,  5.35, names(3),charsize=1.25,color=dev.cnames.black
xyouts, .04,  5.05, names(4),charsize=1.25,color=dev.cnames.black
xyouts, .04,  4.75, names(5),charsize=1.00,color=dev.cnames.black
xyouts, .04,  4.45, names(6),charsize=1.00,color=dev.cnames.black
for m=1,nr do print, x(m), chiim(m), chiit(m,1)

;B & W stuff: (comment out "dev=..." line above)
;oplot,x(*),tix(*),thick=2,psym=8,symsize=.65
;oplot, x(*),tim(*),line=0,thick=2
;oplot, x(*),tem(*),line=2,thick=2


set_plot,'PS'
device, file='chiprof.ps', bits=8,/color
plot, x(*),chiem(*),title='Diffusivity vs rho',xtitle='rho',ytitle='chi (m^2/s)',line=0,thick=2,yrange=[0,6],ystyle=1,color=dev.cnames.red
usersym,cos(a),sin(a)
oplot, x(*),chiim(*),line=0,thick=2,color=dev.cnames.blue
oplot, x(*),chiineo(*),line=0,thick=2,color=dev.cnames.green
xyouts, .70,  5.65, names(0),charsize=1.25,color=dev.cnames.blue
xyouts, .70,  5.35, names(1),charsize=1.25,color=dev.cnames.red
xyouts, .70,  5.05, names(7),charsize=1.25,color=dev.cnames.green
xyouts, .04,  5.65, names(2),charsize=1.25,color=dev.cnames.black
xyouts, .04,  5.35, names(3),charsize=1.25,color=dev.cnames.black
xyouts, .04,  5.05, names(4),charsize=1.25,color=dev.cnames.black
xyouts, .04,  4.75, names(5),charsize=1.00,color=dev.cnames.black
xyouts, .04,  4.45, names(6),charsize=1.00,color=dev.cnames.black
for m=1,nr do print, x(m), chiim(m), chiit(m,1)

device,/close
end
