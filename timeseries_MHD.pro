;Plot the current density with shock ID's

;use t=70,MHD
resolve_routine,'shockid_fun'

tarr=[30,50,60,80,100]

for ti=0,n_elements(tarr)-1 do begin

fname='Sims/MHD_test' & tread=tarr(ti)

shockid_fun,shocks,fname=fname,tread=tread

rdmpi,ds,datapath=fname,time_step=tread,current=1,var=['bx','by','bz']
i=image(ds.j_z,ds.x,ds.y,xr=[0.0,1.0],yr=[0.0,1.0],axis_sty=2,dim=[600,600],xtitle='x',ytitle='y',/buffer,font_size=14)
;v=streamline(ds.bx,ds.by,ds.x+5*dx,ds.y+5*dx,streamline_stepsize=0.01,x_stream=101,y_stream=101,arrow_size=0.1,/overplot,streamline_nstep=200)

xshift=0
yshift=0

sth=3

ps=plot(1.0*shocks.slowx/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.slowy/(1.0*n_elements(shocks.y)-yshift),/overplot,'.r',sym_thick=sth)

pf=plot(1.0*shocks.fastx/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.fasty/(1.0*n_elements(shocks.y)-yshift),/overplot,'.b',sym_thick=sth)

pi1=plot(1.0*shocks.int1x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int1y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.y',sym_thick=sth)

pi2=plot(1.0*shocks.int2x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int2y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.m',sym_thick=sth)

pi3=plot(1.0*shocks.int3x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int3y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.c',sym_thick=sth)

pi4=plot(1.0*shocks.int4x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int4y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.g',sym_thick=sth)

i.save,'timeseries_MHD_'+strtrim(tread,1)+'.pdf'

i.close

endfor

END
