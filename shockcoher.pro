resolve_routine,'shockid_fun'

;ti=70
;tiarr=[30,60,80,100]
tiarr=[30:100:1]
car=['k','b','r','g']

rad=3

cocount=dblarr(2,(2*rad)^2+1,n_elements(tiarr))

cocount2=dblarr(2,n_elements(tiarr))

for ti=0,n_elements(tiarr)-1 do begin

fname='/home/snow/datadrive/simdata/MHD_OZ_1024/Data' & tread=tiarr(ti)

shockid_fun,shocks,fname=fname,tread=tread

fastco=dblarr(n_elements(shocks.fastx))
slowco=dblarr(n_elements(shocks.slowx))

for i=0,n_elements(shocks.fastx)-1 do begin

    nearx=abs(shocks.fastx-shocks.fastx(i))
    neary=abs(shocks.fasty-shocks.fasty(i))
    nearl=where((nearx le rad) AND (neary le rad))
    fastco(i)=n_elements(nearl)-1

endfor

for i=0,n_elements(shocks.slowx)-1 do begin

    nearx=abs(shocks.slowx-shocks.slowx(i))
    neary=abs(shocks.slowy-shocks.slowy(i))
    nearl=where((nearx le rad) AND (neary le rad))
    slowco(i)=n_elements(nearl)-1

endfor

for i=0,(2*rad)^2 do begin

    cocount(0,i,ti)=0.0
    cocount(1,i,ti)=0.0

    a=where(fastco eq i)
    if (a(0) ne -1) then begin
    cocount(0,i,ti)=double(n_elements(a))/double(n_elements(shocks.fastx))
    endif
    a=where(slowco eq i)
    if (a(0) ne -1) then begin
    cocount(1,i,ti)=double(n_elements(a))/double(n_elements(shocks.slowx))
    endif
;print,a
endfor

cocount2(0,ti)=total(cocount(0,6:(2*rad)^2,ti))
cocount2(1,ti)=total(cocount(1,6:(2*rad)^2,ti))

endfor

;p=plot(cocount(0,*,0),color=car(0),xtitle='Number of similar shocks in radius', ytitle='Normalised occurance',thick=2,xr=[0,15])
;for ti=1,n_elements(tiarr)-1 do begin
;    p=plot(cocount(0,*,ti),/overplot,color=car(ti),thick=2)
;endfor
;for ti=0,n_elements(tiarr)-1 do begin
;    p=plot(cocount(1,*,ti),/overplot,color=car(ti),thick=2,linesty=2)
;endfor

save,cocount2,tiarr,filename='shockcoher_data.sav'

p=plot(0.01*double(tiarr),cocount2(0,*),'b',thick=2,dim=[600,600],xtitle='Time',ytitle='Proportion of coherent shocks',font_size=16)
p=plot(0.01*tiarr,cocount2(1,*),/overplot,'r',thick=2)

p.save,'shockcoher.pdf'

END
