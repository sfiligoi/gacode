pro powiplt,names

; Procedure to plot powi_m, powi_exp profiles

; openr,1,'out'
; readf,1,ny
; f=fltarr(8,ny)
; readf,1,f
; close,1

id = ncdf_open('xptor_out.nc')
ncdf_varget,id,ncdf_varid(id,'nrho_bc'),nr
ncdf_varget,id,ncdf_varid(id,'rho'),x
ncdf_varget,id,ncdf_varid(id,'powi_exp'),powix
ncdf_varget,id,ncdf_varid(id,'powi_m'),powim
; printf,-1,f

names=strarr(7)
names[0]='Ti'
names[1]='Te'
names[2]='DIII-D #101381'
names[3]='t=2.63 secs'
names[4]='GLF23 v1.61 model'
names[5]='Case DV #7, ne,Te,Ti,vphi,vpol evolved'
names[6]='w/ NCLASS, BC@rho=0.90'

set_plot,'x'
plotwhite
!p.charsize=1.25
a=findgen(17) * (!pi*2/16.)
usersym,cos(a),sin(a),/fill
dev=devsetup(ctable='cmg')
plot, x(*),powix(*),title='Ion Power',xtitle='rho',ytitle='T (kev)',thick=2,yrange=[0,1],color=dev.cnames.red
usersym,cos(a),sin(a)
oplot,x(*),powim(1:nr),thick=2,color=dev.cnames.blue
xyouts, .85, 0.95, names(0),charsize=1.50,color=dev.cnames.blue
xyouts, .85, 0.90, names(1),charsize=1.50,color=dev.cnames.red
print,'         rho           Powi_exp         Powi_m'
for m=1,nr do print, x(m), powix(m), powim(m)

;B & W stuff: (comment out "dev=..." line above)
;oplot,x(*),tix(*),thick=2,psym=8,symsize=.65
;oplot, x(*),tim(*),line=0,thick=2
;oplot, x(*),tem(*),line=2,thick=2


set_plot,'PS'
device, file='powi.ps', /color, bits=8
plotwhite
!p.charsize=1.25
a=findgen(17) * (!pi*2/16.)
usersym,cos(a),sin(a),/fill
plot, x(*),powix(*),title='Ion Power vs rho',xtitle='rho',ytitle='T (kev)',thick=2,yrange=[0,1],color=dev.cnames.red
usersym,cos(a),sin(a)
oplot,x(*),powim(1:nr),thick=2,color=dev.cnames.blue
xyouts, .85, 0.95, names(0),charsize=1.50,color=dev.cnames.blue
xyouts, .85, 0.90, names(1),charsize=1.50,color=dev.cnames.red

device,/close
end
