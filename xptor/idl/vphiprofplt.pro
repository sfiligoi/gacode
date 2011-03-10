pro vphiprofplt

; Procedure to plot experimental and model vphi profiles (in 2nd group of 'out')

;openr,1,'out'
;readf,1,ny
;f=fltarr(8,ny)
;readf,1,f
;f=fltarr(7,ny)
;readf,1,f
;close,1

id = ncdf_open('xptor_out.nc')
ncdf_varget,id,ncdf_varid(id,'rho'),x
ncdf_varget,id,ncdf_varid(id,'nrho_bc'),nr
ncdf_varget,id,ncdf_varid(id,'vphip_exp'),vphi_exp
ncdf_varget,id,ncdf_varid(id,'vphi_m'),vphi_m
vphi_m=abs(vphi_m)

label=strarr(6)
label[0]='Data'
label[1]='Model'
label[2]='DIII-D #98549'
label[3]='time=1.50 secs'
label[4]='Case DV #5'
label[5]='chiphi=chieneo+x11*chie'

set_plot,'x'
plotwhite
!p.charsize=1.25
a=findgen(17) * (!pi*2/16.)
usersym,cos(a),sin(a),/fill
dev=devsetup(ctable='cmg')
plot, x(*),vphi_exp(*),title='Main Ion Toroidal Velocity',xtitle='rho',ytitle='vphi (m/s)',thick=2,psym=8,symsize=.65,yrange=[0.e5,5.e5],ystyle=1,color=dev.cnames.black
usersym,cos(a),sin(a)
oplot, x(*), vphi_m(1:nr+1),line=0,thick=2,color=dev.cnames.blue
;oplot, f(1,*),f(6,*),line=1,thick=2
;oplot, f(1,*),f(5,*),line=2,thick=2
xyouts, .85, 4.6e5, label(0),charsize=1.25,color=dev.cnames.black
xyouts, .85, 4.2e5, label(1),charsize=1.25,color=dev.cnames.blue
xyouts, .05, 4.6e5, label(2),charsize=1.25,color=dev.cnames.black
xyouts, .05, 4.35e5, label(3),charsize=1.25,color=dev.cnames.black
xyouts, .05, 4.10e5, label(4),charsize=1.25,color=dev.cnames.black
xyouts, .05, 3.85e5, label(5),charsize=1.0,color=dev.cnames.black

set_plot,'PS'
device, file='vphiprof.ps', bits=8
plotwhite
!p.charsize=1.25
a=findgen(17) * (!pi*2/16.)
usersym,cos(a),sin(a),/fill
plot, x(*),vphi_exp(*),title='Main Ion Toroidal Velocity',xtitle='rho',ytitle='vphi (m/s)',thick=2,psym=8,symsize=.65,yrange=[0.e5,5.e5],ystyle=1,color=dev.cnames.black
usersym,cos(a),sin(a)
oplot, x(*),vphi_m(1:nr+1),line=0,thick=2,color=dev.cnames.blue
xyouts, .85, 4.6e5, label(0),charsize=1.25,color=dev.cnames.black
xyouts, .85, 4.2e5, label(1),charsize=1.25,color=dev.cnames.blue
xyouts, .05, 4.6e5, label(2),charsize=1.25,color=dev.cnames.black
xyouts, .05, 4.35e5, label(3),charsize=1.25,color=dev.cnames.black
xyouts, .05, 4.10e5, label(4),charsize=1.25,color=dev.cnames.black

device,/close
end
