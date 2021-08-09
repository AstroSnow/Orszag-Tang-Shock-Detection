;Plot the current density with shock ID's

;use t=70,MHD
resolve_routine,'shockid_fun'

tarr= dblarr(101)
nslow=dblarr(101)
nfast=dblarr(101)
nint1=dblarr(101)
nint2=dblarr(101)
nint3=dblarr(101)
nint4=dblarr(101)

for ti=0,100 do begin

fname='/home/snow/datadrive/simdata/MHD_OZ_1024/Data' & tread=ti

shockid_fun,shocks,fname=fname,tread=tread

tarr(ti)=ti
if (shocks.slowx(0) eq 0) and (shocks.slowy(0) eq 0) and (N_elements(shocks.slowx) eq 1) then nslow(ti)=0 else $ 
nslow(ti)=n_elements(shocks.slowx)
if (shocks.fastx(0) eq 0) and (shocks.fasty(0) eq 0) and (N_elements(shocks.fastx) eq 1) then nfast(ti)=0 else $ 
nfast(ti)=n_elements(shocks.fastx)
if (shocks.int1x(0) eq 0) and (shocks.int1y(0) eq 0) and (N_elements(shocks.int1x) eq 1) then nint1(ti)=0 else $ 
nint1(ti)=n_elements(shocks.int1x)
if (shocks.int2x(0) eq 0) and (shocks.int2y(0) eq 0) and (N_elements(shocks.int2x) eq 1) then nint2(ti)=0 else $ 
nint2(ti)=n_elements(shocks.int2x)
if (shocks.int3x(0) eq 0) and (shocks.int3y(0) eq 0) and (N_elements(shocks.int3x) eq 1) then nint3(ti)=0 else $ 
nint3(ti)=n_elements(shocks.int3x)
if (shocks.int4x(0) eq 0) and (shocks.int4y(0) eq 0) and (N_elements(shocks.int4x) eq 1) then nint4(ti)=0 else $ 
nint4(ti)=n_elements(shocks.int4x)

endfor

save,nslow,nfast,nint1,nint2,nint3,nint4,tarr,filename='shockcounts_MHD.sav'

tarr=tarr*0.01

ps=plot(tarr,nslow,/ylog,yr=[1,10000],'r',thick=3,dim=[970,800],xtitle='Time',ytitle='Counts', position=[0.15,0.1,0.95,0.8],font_size=16,name='slow (3-4)')
pf=plot(tarr,nfast,/overplot,'b',thick=3,name='fast (1-2)')
pi1=plot(tarr,nint1,/overplot,'y',thick=3,name='intermediate (1-3)')
pi2=plot(tarr,nint2,/overplot,'m',thick=3,name='intermediate (1-4)')
pi3=plot(tarr,nint3,/overplot,'c',thick=3,name='intermediate (2-3)')
pi4=plot(tarr,nint4,/overplot,'g',thick=3,name='intermediate (2-4)')
l1=legend(target=[ps,pf,pi1],position=[0.49,0.97],font_size=16)
l2=legend(target=[pi2,pi3,pi4],position=[0.95,0.97],font_size=16)


END
