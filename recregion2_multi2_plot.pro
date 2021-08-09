;Plot the reconnection region using data from recregion2_multi2.pro

w=window(dim=[1000,600],/buffer)

tread=[0,66,67,68,69,70,71,72,73,74,75]

for ti=1,8 do begin

restore,'MHD_shockdata_t_'+strtrim(tread(ti),1)+'.sav'

ypos=1.0
if (ti le 4) then ypos=0.0
tip=ti
if tip gt 4 then tip=tip-4


i=image(ds.j_z,ds.x,ds.y,xr=[0.7,1.0],yr=[0.2,0.5],axis_sty=2,dim=[970,600],xtitle='x'$
, ytitle='y',/current,position=[0.08+0.2*(tip-1.0),0.5-0.3*ypos,0.08+0.2*(tip),0.8-0.3*ypos],$
aspect_rat=0)
;i=contour(ds.j_z,ds.x,ds.y,xr=[0.7,1.0],yr=[0.2,0.5],axis_sty=2,dim=[970,600],xtitle='x'$
;, ytitle='y',/current,position=[0.08+0.2*(tip-1.0),0.5-0.3*ypos,0.08+0.2*(tip),0.8-0.3*ypos],$
;/fill,rgb_tab=0,n_levels=101)
;v=streamline(ds.bx,ds.by,ds.x+5*dx,ds.y+5*dx,streamline_stepsize=0.01,x_stream=101,y_stream=101,arrow_size=0.1,/overplot,streamline_nstep=200)

ia=i.axes
ia[1].tickv=[0.25,0.3,0.35,0.4,0.45]
ia[0].tickv=[0.75,0.8,0.85,0.9,0.95]
ia[1].show=0
if (ti eq 1) or (ti eq 5) then ia[1].show=1
if (ti lt 5) then ia[0].show=0

xshift=5
yshift=5

symth=3

ps=plot(1.0*shocks.slowx/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.slowy/(1.0*n_elements(shocks.y)-yshift),/overplot,'.r',sym_thick=symth)

pf=plot(1.0*shocks.fastx/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.fasty/(1.0*n_elements(shocks.y)-yshift),/overplot,'.b',sym_thick=symth)

pi1=plot(1.0*shocks.int1x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int1y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.y',sym_thick=symth)

pi2=plot(1.0*shocks.int2x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int2y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.m',sym_thick=symth)

pi3=plot(1.0*shocks.int3x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int3y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.c',sym_thick=symth)

pi4=plot(1.0*shocks.int4x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int4y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.g',sym_thick=symth)

endfor

i.save,'recregion2_multi2_plot.pdf'

END
