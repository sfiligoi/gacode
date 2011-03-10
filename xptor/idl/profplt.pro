pro profplt,names

; Procedure to plot profiles

; openr,1,'out'
; readf,1,ny
; f=fltarr(8,ny)
; readf,1,f
; close,1

id = ncdf_open('xptor_out.nc')
ncdf_varget,id,ncdf_varid(id,'discharge'),shot
ncdf_varget,id,ncdf_varid(id,'xp_time'),time
ncdf_varget,id,ncdf_varid(id,'nrho_bc'),nr
ncdf_varget,id,ncdf_varid(id,'rho'),x
ncdf_varget,id,ncdf_varid(id,'ti_exp'),tix
ncdf_varget,id,ncdf_varid(id,'te_exp'),tex
ncdf_varget,id,ncdf_varid(id,'ti_m'),tim
ncdf_varget,id,ncdf_varid(id,'te_m'),tem
; printf,-1,f

names=strarr(7)
names[0]='Ti'
names[1]='Te'
names[2]=string(shot)
names[3]=time
names[4]='TGLF model'
names[5]='Case #1, Ti,Te evolved'
names[6]='w/ Chang-Hinton, BC@rho=0.90'

tmax=strarr(4)
tixmax=max(tix)
timmax=max(tim)
texmax=max(tex)
temmax=max(tem)
tmax[0]=tixmax
tmax[1]=timmax
tmax[2]=texmax
tmax[3]=temmax
tempmax=max(tmax)
shotime=strarr(3)
shotime[0]=string(shot)
shotime[1]=time
;read,'tokamak = ',tok
;read,'shot = ',shot
print,'        Tix-max          Tim-max          Tex-max          Tem-max'
print,tmax
print,'        T-max =  ',tempmax

set_plot,'x'
plotwhite
!p.charsize=1.25
a=findgen(17) * (!pi*2/16.)
usersym,cos(a),sin(a),/fill
dev=devsetup(ctable='cmg')
plot, x(*),tex(*),title='Temperature',xtitle='rho',ytitle='T (kev)',thick=2,psym=8,symsize=.65,yrange=[0,tempmax+2],color=dev.cnames.red
usersym,cos(a),sin(a)
oplot,x(*),tix(*),thick=2,psym=8,symsize=.65,color=dev.cnames.blue
oplot, x(*),tim(1:nr),line=0,thick=2,color=dev.cnames.blue
oplot, x(*),tem(1:nr),line=0,thick=2,color=dev.cnames.red
xyouts, .85, tempmax+2-0.5, names(0),charsize=1.50,color=dev.cnames.blue
xyouts, .85, tempmax+2-0.8, names(1),charsize=1.50,color=dev.cnames.red
xyouts, .03, tempmax+2-0.4, shotime,charsize=1.25,color=dev.cnames.black
; xyouts, .03, tempmax+2-0.6, names(3),charsize=1.25,color=dev.cnames.black,alignment=0.4
xyouts, .03, tempmax+2-0.6, names(4),charsize=1.25,color=dev.cnames.black
xyouts, .03, tempmax+2-0.8, names(5),charsize=1.25,color=dev.cnames.black
;xyouts, .03, 9.25, names(6),charsize=1.00,color=dev.cnames.black


;B & W stuff: (comment out "dev=..." line above)
;oplot,x(*),tix(*),thick=2,psym=8,symsize=.65
;oplot, x(*),tim(*),line=0,thick=2
;oplot, x(*),tem(*),line=2,thick=2


set_plot,'PS'
device, file='profs.ps', /color, bits=8
plotwhite
!p.charsize=1.25
a=findgen(17) * (!pi*2/16.)
usersym,cos(a),sin(a),/fill
plot, x(*),tex(*),title='Temperature vs rho',xtitle='rho',ytitle='T (kev)',thick=2,psym=8,symsize=.65,yrange=[0,tempmax+2],color=dev.cnames.red
usersym,cos(a),sin(a)
oplot,x(*),tix(*),thick=2,psym=8,symsize=.65,color=dev.cnames.blue
oplot, x(*),tim(1:nr),line=0,thick=2,color=dev.cnames.blue
oplot, x(*),tem(1:nr),line=0,thick=2,color=dev.cnames.red
xyouts, .85, tempmax+2-0.5, names(0),charsize=1.50,color=dev.cnames.blue
xyouts, .85, tempmax+2-0.8, names(1),charsize=1.50,color=dev.cnames.red
xyouts, .03, tempmax+2-0.4, shotime,charsize=1.25,color=dev.cnames.black
; xyouts, .03, tempmax+2-0.6, names(3),charsize=1.25,color=dev.cnames.black,alignment=0.4
xyouts, .03, tempmax+2-0.6, names(4),charsize=1.25,color=dev.cnames.black
xyouts, .03, tempmax+2-0.8, names(5),charsize=1.25,color=dev.cnames.black
;xyouts, .03, 9.25, names(6),charsize=1.00,color=dev.cnames.black

device,/close
end
