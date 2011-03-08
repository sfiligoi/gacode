pro nprofplt

; Procedure to plot ne profile

;openr,1,'out2'
;readf,1,ny
;f=fltarr(8,ny)
;readf,1,f
;close,1

id = ncdf_open('xptor_out.nc')
ncdf_varget,id,ncdf_varid(id,'rho'),x
ncdf_varget,id,ncdf_varid(id,'ne_exp'),nex
ncdf_varget,id,ncdf_varid(id,'ni_exp'),nix
ncdf_varget,id,ncdf_varid(id,'ne_m'),nem
ncdf_varget,id,ncdf_varid(id,'ni_m'),nim

label=strarr(6)
label[0]='ne'
label[1]='ni'
label[2]='DIII-D #82205'
label[3]='time=2.50 secs'
label[4]='alpha_E=1.35'
label[5]='Case DV #3'

set_plot,'x'
plotwhite
!p.charsize=1.25
a=findgen(17) * (!pi*2/16.)
usersym,cos(a),sin(a),/fill
;plot, f(1,*),f(2,*),title='Electron Density',xtitle='rho',ytitle='ne (10^19 m^-3)',thick=2,line=0,symsize=.65,yrange=[0,3]
dev=devsetup(ctable='cmg')
plot, x(*),nex(*),title='Density Profile',xtitle='rho',ytitle='n (10^19 m^-3)',thick=2,psym=8,symsize=0.65,yrange=[0,12],ystyle=1,color=dev.cnames.red
usersym,cos(a),sin(a)
oplot, x(*),nix(*),thick=2,psym=8,symsize=.65,color=dev.cnames.blue
oplot, x(*),nim(*),line=0,thick=2,color=dev.cnames.blue
oplot, x(*),nem(*),line=0,thick=2,color=dev.cnames.red
;oplot, f(1,*),f(4,*),line=0,thick=2
;oplot, f(1,*),f(5,*),line=2,thick=2
xyouts, .85, 11.00, label(0),charsize=1.25,color=dev.cnames.red
xyouts, .85, 10.25, label(1),charsize=1.25,color=dev.cnames.blue
xyouts, .03, 11.00, label(2),charsize=1.25,color=dev.cnames.black
xyouts, .03, 10.35, label(3),charsize=1.25,color=dev.cnames.black
xyouts, .03, 9.80, label(4),charsize=1.00,color=dev.cnames.black
xyouts, .03, 9.30, label(5),charsize=1.00,color=dev.cnames.black


set_plot,'PS'
device, file='nprof.ps', /color, bits=8
plotwhite
!p.charsize=1.25
a=findgen(17) * (!pi*2/16.)
usersym,cos(a),sin(a),/fill
plot, x(*),nex(*),title='Density Profile',xtitle='rho',ytitle='n (10^19 m^-3)',thick=2,psym=8,symsize=0.65,yrange=[0,12],ystyle=1,color=dev.cnames.red
usersym,cos(a),sin(a)
oplot, x(*),nix(*),thick=2,psym=8,symsize=.65,color=dev.cnames.blue
oplot, x(*),nim(*),line=0,thick=2,color=dev.cnames.blue
oplot, x(*),nem(*),line=0,thick=2,color=dev.cnames.red
;oplot, f(1,*),f(4,*),line=0,thick=2
;oplot, f(1,*),f(5,*),line=2,thick=2
xyouts, .85, 11.00, label(0),charsize=1.25,color=dev.cnames.red
xyouts, .85, 10.25, label(1),charsize=1.25,color=dev.cnames.blue
xyouts, .03, 11.00, label(2),charsize=1.25,color=dev.cnames.black
xyouts, .03, 10.35, label(3),charsize=1.25,color=dev.cnames.black
xyouts, .03, 9.80, label(4),charsize=1.00,color=dev.cnames.black
xyouts, .03, 9.30, label(5),charsize=1.00,color=dev.cnames.black

device,/close
end
