function cmgcolors,show=show,small=small
r=[255b,  0b,  0b,  0b,255b,140b,255b,  0b,127b,127b,255b,  0b,192b,255b]
g=[255b,  0b,  0b,192b,  0b, 23b,127b,192b, 30b,127b,230b,191b,192b,  0b]
b=[255b,  0b,255b,  0b,  0b,136b,  0b,192b,  2b,127b,  0b,127b,192b,127b]
if keyword_set(small) then begin
    r=[r[0:4],r[10]]
    g=[g[0:4],g[10]]
    b=[b[0:4],b[10]]
    colors={white:0,black:1,blue:2,green:3,red:4,yellow:5}
endif else colors={white:0,black:1,blue:2,green:3,red:4,purple:5, $
                   orange:6,cyan:7,brown:8,grey:9,yellow:10,bluegreen:11, $
                   lightgrey:12, magenta:13}
tvlct,r,g,b
!p.color=colors.black & !p.background=colors.white
if keyword_set(show) then begin
    ncolors=n_elements(r)
    wsave=!d.window
    psave=!p & xsave=!x & ysave=!y
    window,/free,title='Current color table'
    nrow=round(sqrt(ncolors))
    ncol=round(float(ncolors+1)/float(nrow))
    !p.multi=[0,nrow,ncol]
    !x.tickname=replicate(' ',60) & !y.tickname=!x.tickname
    !x.margin=[0,0] & !y.margin=[0,0]
    dcolor=ncolors/2 & dt1=ncolors-1
    for i=0,dt1 do begin
        plot,[0,1],[0,1],/nodata
        polyfill,[0,1,1,0,0],[0,0,1,1,0],color=i
        xyouts,0.5,0.5,strcompress(string(i),/remove_all),alignment=0.5, $
          charsize=3,charthick=3,color=((i+dcolor) mod dt1)
    endfor
    !p=psave & !x=xsave & !y=ysave
    wset,wsave
endif

return,create_struct(colors,'table_size',n_elements(r))
end
