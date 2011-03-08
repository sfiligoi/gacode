pro qprofplt

; Procedure to plot q profile

;openr,1,'out2'
;readf,1,ny
;f=fltarr(6,ny)
;readf,1,f
;close,1

id = ncdf_open('xptor_out.nc')
ncdf_varget,id,ncdf_varid(id,'rho'),x
ncdf_varget,id,ncdf_varid(id,'q_exp'),q
ncdf_varget,id,ncdf_varid(id,'alpha_exp'),alpx

names=strarr(8)
names[0]='DIII-D #104276'
names[1]='t=5.71 secs'
names[2]='q-Exp'
names[3]='alpha-Exp'

set_plot,'x'
plotwhite
!p.charsize=1.25
a=findgen(17) * (!pi*2/16.)
usersym,cos(a),sin(a),/fill
dev=devsetup(ctable='cmg')
;plot, f(1,*),f(5,*),title='Safety Factor',xtitle='rho',ytitle='q',line=2,thick=2,yrange=[0,5]
plot, x(*),q(*),title='Safety Factor',xtitle='rho',ytitle='q',line=0,thick=2,yrange=[0,6]
usersym,cos(a),sin(a)
;oplot, f(1,*),f(7,*),line=0,thick=2
;oplot, f(1,*),f(4,*),line=0,thick=2
;oplot, f(1,*),f(5,*),line=2,thick=2
oplot, x(*),alpx(*),line=0,thick=2,color=dev.cnames.blue
xyouts, .03, 5.50, names(0),charsize=1.25,color=dev.cnames.black
xyouts, .03, 5.15, names(1),charsize=1.25,color=dev.cnames.black
xyouts, .75, 5.50, names(2),charsize=1.25,color=dev.cnames.black
xyouts, .75, 5.15, names(3),charsize=1.25,color=dev.cnames.blue

set_plot,'PS'
device, file='qprof.ps', bits=8
plotwhite
!p.charsize=1.25
a=findgen(17) * (!pi*2/16.)
usersym,cos(a),sin(a),/fill
plot, x(*),q(*),title='Safety Factor',xtitle='rho',ytitle='q',thick=2,psym=8,symsize=.65,yrange=[0,6]
usersym,cos(a),sin(a)
;oplot, f(1,*),f(7,*),line=0,thick=2
;oplot, f(1,*),f(4,*),line=0,thick=2
;oplot, f(1,*),f(5,*),line=2,thick=2
xyouts, .03, 5.50, names(0),charsize=1.25,color=dev.cnames.black
xyouts, .03, 5.15, names(1),charsize=1.25,color=dev.cnames.black
xyouts, .75, 5.50, names(2),charsize=1.25,color=dev.cnames.black
xyouts, .75, 5.15, names(3),charsize=1.25,color=dev.cnames.blue

device,/close
end
