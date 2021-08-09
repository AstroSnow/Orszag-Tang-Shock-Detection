;Identify and clasify shocks based on the SHOCKFIND algorythm (https://arxiv.org/abs/1608.02050)
pro shockid_fun,shocks,fname=fname,tread=tread
  Compile_Opt DEFINT32  

;Compile a few necessary function
RESOLVE_ROUTINE, 'exist', /IS_FUNCTION
RESOLVE_ROUTINE, 'where_vector', /IS_FUNCTION
RESOLVE_ROUTINE, 'interp2d', /IS_FUNCTION
RESOLVE_ROUTINE, 'datatype', /IS_FUNCTION
RESOLVE_ROUTINE, 'shocknormvals', /IS_FUNCTION

;Input directory
;fname='../MHD_test'
;fname='Sims/MHD_test' & tread=70
;fname='Sims/PIP_test' & tread=70
;fname='Sims/PIP_test_ac100' & tread=100
;fname='../Shock_stab/PIP/Sims/Data_4000x200_MHD_test' & tread=50
;fname='../Shock_stab/PIP/Sims/MHD_sf_o4_hr' & tread=30
;fname='/media/snow/datadrive/OZdata/MHD_OZ_16384/Data' & tread=100

;Read in the data
print,'Reading data'
rdmpi,ds,datapath=fname,time_step=tread

;Parameters for shock limits
convl=0.01 ;Convergence threshold
avecyl=3 ;Cylinder to average over
smthfac=1 ;smoothing factor (1=no smoothing)
shocktol=0.05 ;tolerence for the shock transitions
bulkspeed=0 ; Use the bulk sound and alfven speeds or individual
species='plasma' ; plasma or neutral
;species='neutral' ; plasma or neutral

;Data handeling flags
data_subset=0 ;flag to subset the data
data_save=0


;Define simulation grid
print,'Removing ghost cells'
margin=6

egx=N_elements(ds.x)-1
egy=N_elements(ds.y)-1

x=ds.x(margin:egx-margin)
y=ds.y(margin:egy-margin)
;MHD variables
if (ds.fl_pip eq 0) and (ds.fl_mhd eq 1) then begin
    print,'Using MHD quantities'
    ro=smooth(ds.ro_p(margin:egx-margin,margin:egy-margin),smthfac)
    vx=smooth(ds.vx_p(margin:egx-margin,margin:egy-margin),smthfac)
    vy=smooth(ds.vy_p(margin:egx-margin,margin:egy-margin),smthfac)
    vz=smooth(ds.vz_p(margin:egx-margin,margin:egy-margin),smthfac)
    bx=smooth(ds.bx(margin:egx-margin,margin:egy-margin),smthfac)
    by=smooth(ds.by(margin:egx-margin,margin:egy-margin),smthfac)
    bz=smooth(ds.bz(margin:egx-margin,margin:egy-margin),smthfac)
    pr=smooth(ds.pr_p(margin:egx-margin,margin:egy-margin),smthfac)
endif
;PIP plasma only variables
if (ds.fl_pip eq 1) and (bulkspeed eq 0) and (species eq 'plasma') then begin
    print,'Using plasma only quantities'
    ro=smooth(ds.ro_p(margin:egx-margin,margin:egy-margin),smthfac)
    vx=smooth(ds.vx_p(margin:egx-margin,margin:egy-margin),smthfac)
    vy=smooth(ds.vy_p(margin:egx-margin,margin:egy-margin),smthfac)
    vz=smooth(ds.vz_p(margin:egx-margin,margin:egy-margin),smthfac)
    bx=smooth(ds.bx(margin:egx-margin,margin:egy-margin),smthfac)
    by=smooth(ds.by(margin:egx-margin,margin:egy-margin),smthfac)
    bz=smooth(ds.bz(margin:egx-margin,margin:egy-margin),smthfac)
    pr=smooth(ds.pr_p(margin:egx-margin,margin:egy-margin),smthfac)
endif
;PIP plasma bulk variables
if (ds.fl_pip eq 1) and (bulkspeed eq 1) and (species eq 'plasma') then begin
    print,'Using plasma bulk quantities'
    ro=smooth(ds.ro_p(margin:egx-margin,margin:egy-margin),smthfac)+$
        smooth(ds.ro_n(margin:egx-margin,margin:egy-margin),smthfac)
    vx=smooth(ds.vx_p(margin:egx-margin,margin:egy-margin),smthfac)
    vy=smooth(ds.vy_p(margin:egx-margin,margin:egy-margin),smthfac)
    vz=smooth(ds.vz_p(margin:egx-margin,margin:egy-margin),smthfac)
    bx=smooth(ds.bx(margin:egx-margin,margin:egy-margin),smthfac)
    by=smooth(ds.by(margin:egx-margin,margin:egy-margin),smthfac)
    bz=smooth(ds.bz(margin:egx-margin,margin:egy-margin),smthfac)
    pr=smooth(ds.pr_p(margin:egx-margin,margin:egy-margin),smthfac)+$
        smooth(ds.pr_n(margin:egx-margin,margin:egy-margin),smthfac)
endif
;PIP plasma only variables
if (ds.fl_pip eq 1) and (bulkspeed eq 0) and (species eq 'neutral') then begin
    print,'Using neutral only quantities'
    ro=smooth(ds.ro_n(margin:egx-margin,margin:egy-margin),smthfac)
    vx=smooth(ds.vx_n(margin:egx-margin,margin:egy-margin),smthfac)
    vy=smooth(ds.vy_n(margin:egx-margin,margin:egy-margin),smthfac)
    vz=smooth(ds.vz_n(margin:egx-margin,margin:egy-margin),smthfac)
    pr=smooth(ds.pr_n(margin:egx-margin,margin:egy-margin),smthfac)
endif

;2D divergence of velocity
print,'Calculating divergence'
divv=dblarr(n_elements(x),n_elements(y))
if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
divv=(smooth(ds.vx_p(margin+1:egx-margin+1,margin:egy-margin),smthfac)-smooth(ds.vx_p(margin-1:egx-margin-1,margin:egy-margin),smthfac))/(ds.x(10)-ds.x(8)) $
+(smooth(ds.vy_p(margin:egx-margin,margin+1:egy-margin+1),smthfac)-smooth(ds.vy_p(margin:egx-margin,margin-1:egy-margin-1),smthfac))/(ds.y(10)-ds.y(8))
endif
if (species eq 'neutral') and (ds.fl_pip eq 1) then begin
divv=(smooth(ds.vx_n(margin+1:egx-margin+1,margin:egy-margin),smthfac)-smooth(ds.vx_n(margin-1:egx-margin-1,margin:egy-margin),smthfac))/(ds.x(10)-ds.x(8)) $
+(smooth(ds.vy_n(margin:egx-margin,margin+1:egy-margin+1),smthfac)-smooth(ds.vy_n(margin:egx-margin,margin-1:egy-margin-1),smthfac))/(ds.y(10)-ds.y(8))
endif

;Maximum gradient of density
print,'Calculating density gradients'
gradx=dblarr(n_elements(x),n_elements(y))
grady=dblarr(n_elements(x),n_elements(y))
gradz=dblarr(n_elements(x),n_elements(y))
if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
gradx=(smooth(ds.ro_p(margin+1:egx-margin+1,margin:egy-margin),smthfac)-smooth(ds.ro_p(margin-1:egx-margin-1,margin:egy-margin),smthfac))/(ds.x(10)-ds.x(8))
grady=(smooth(ds.ro_p(margin:egx-margin,margin+1:egy-margin+1),smthfac)-smooth(ds.ro_p(margin:egx-margin,margin-1:egy-margin-1),smthfac))/(ds.y(10)-ds.y(8))
gradz=0.0*dblarr(n_elements(x),n_elements(y))
endif
if (species eq 'neutral') and (ds.fl_pip eq 1) then begin
gradx=(smooth(ds.ro_n(margin+1:egx-margin+1,margin:egy-margin),smthfac)-smooth(ds.ro_n(margin-1:egx-margin-1,margin:egy-margin),smthfac))/(ds.x(10)-ds.x(8))
grady=(smooth(ds.ro_n(margin:egx-margin,margin+1:egy-margin+1),smthfac)-smooth(ds.ro_n(margin:egx-margin,margin-1:egy-margin-1),smthfac))/(ds.y(10)-ds.y(8))
gradz=0.0*dblarr(n_elements(x),n_elements(y))
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Step i: identify candidate cells based on 
print,'Identify candidate cells'
gradmag=sqrt(gradx^2+grady^2+gradz^2)

;Find candidate cells
index=where(divv LT -convl)
s = SIZE(divv)
ncol = s(1)
col = index MOD ncol
row = index / ncol

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Remove cells that are not a local maximum
print,'Removing non-local maximum density gradients from candidate cells'

col2=[]
row2=[]

for i=0,n_elements(col)-1 do begin
;for i=0,225 do begin

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Step ii: find shock normal based on density gradient
    normx=gradx(col(i),row(i))/gradmag(col(i),row(i))
    normy=grady(col(i),row(i))/gradmag(col(i),row(i))

    ;Find perpendicular direction too
    perpx=(-normy/normx)/sqrt(normy^2/normx^2+1.0)
    perpy=1.0/sqrt(normy^2/normx^2+1.0)

    tempx=indgen(2*avecyl+1)-avecyl+col(i)
    tempy=indgen(2*avecyl+1)-avecyl+row(i)
    tempx2=indgen(2*avecyl+1)-avecyl
    tempy2=indgen(2*avecyl+1)-avecyl

    ;use periodic BC to fix negative values
    for ii=0,2*avecyl do begin
        if tempx(ii) lt 0 then tempx(ii)=N_elements(divv(*,0))+tempx(ii)
        if tempy(ii) lt 0 then tempy(ii)=N_elements(divv(0,*))+tempy(ii)
    endfor

    ;use periodic BC to fix outside grid values
    for ii=0,2*avecyl do begin
        if tempx(ii) gt n_elements(divv(*,0))-1 then tempx(ii)=tempx(ii)-N_elements(divv(*,0))
        if tempy(ii) gt n_elements(divv(0,*))-1 then tempy(ii)=tempy(ii)-N_elements(divv(0,*))
    endfor


    ;Interpolate the values normal to the shock NOTE: THES ARE NOT PARALLEL TO SHOCK YET!
    ronorm=shocknormvals(ro,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)

    a=max(abs(deriv(ronorm)),b)
    if b eq avecyl then begin
        row2=[row2,row(i)]
        col2=[col2,col(i)]
    endif
endfor
row=row2
col=col2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Filter the candidate cells for testing
if data_subset eq 1 then begin
print,'Subsetting data'
xminf=0.28
xmaxf=0.35
yminf=0.22
ymaxf=0.3

col2=[]
row2=[]
for i=0,n_elements(row)-1 do begin
    if (x(col(i)) LT xmaxf) and (x(col(i)) GT xminf) and (y(row(i)) LT ymaxf) and (y(row(i)) GT yminf) then begin
            row2=[row2,row(i)]
            col2=[col2,col(i)]
    endif
endfor

col=col2
row=row2
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Arrays for the shock locations
fastx=[]
fasty=[]
fastm=[]

slowx=[]
slowy=[]
slowm=[]

int1x=[]
int1y=[]

int2x=[]
int2y=[]

int3x=[]
int3y=[]

int4x=[]
int4y=[]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if data_save eq 1 then begin
;create arrays to store shock normal data
arrronorm=dblarr(N_elements(col),avecyl*2+1)
;arrvxnorm=dblarr(N_elements(col),avecyl*2+1)
;arrvynorm=dblarr(N_elements(col),avecyl*2+1)
;arrbxnorm=dblarr(N_elements(col),avecyl*2+1)
;arrbynorm=dblarr(N_elements(col),avecyl*2+1)
arrprnorm=dblarr(N_elements(col),avecyl*2+1)
arrvpar=dblarr(N_elements(col),avecyl*2+1)
arrvperp=dblarr(N_elements(col),avecyl*2+1)
arrbpar=dblarr(N_elements(col),avecyl*2+1)
arrbperp=dblarr(N_elements(col),avecyl*2+1)
arrvslow2=dblarr(N_elements(col),avecyl*2+1)
arrvalf2=dblarr(N_elements(col),avecyl*2+1)
arrvfast2=dblarr(N_elements(col),avecyl*2+1)
arrmf=dblarr(N_elements(col),avecyl*2+1)
arrma=dblarr(N_elements(col),avecyl*2+1)
arrms=dblarr(N_elements(col),avecyl*2+1)
arrtype=dblarr(N_elements(col))
arrvs=dblarr(N_elements(col))
;arrnx=dblarr(N_elements(col))
;arrny=dblarr(N_elements(col))
;arrpx=dblarr(N_elements(col))
;arrpy=dblarr(N_elements(col))
;arrrescode=dblarr(N_elements(col))
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print,'Begin loop'
;loop over candidates
for i=0,n_elements(col)-1 do begin
;for i=0,225 do begin

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Step ii: find shock normal based on density gradient
    normx=gradx(col(i),row(i))/gradmag(col(i),row(i))
    normy=grady(col(i),row(i))/gradmag(col(i),row(i))

    ;Find perpendicular direction too
    perpx=(-normy/normx)/sqrt(normy^2/normx^2+1.0)
    perpy=1.0/sqrt(normy^2/normx^2+1.0)

    ;print,normx,normy,gradx(col(i),row(i)),grady(col(i),row(i))

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Step iii) 

    ;First create a (2*avecyl)^2 box of values
;    tempa=dblarr(2*avecyl+1,2*avecyl+1)
    tempx=indgen(2*avecyl+1)-avecyl+col(i)
    tempy=indgen(2*avecyl+1)-avecyl+row(i)
    tempx2=indgen(2*avecyl+1)-avecyl
    tempy2=indgen(2*avecyl+1)-avecyl

    ;use periodic BC to fix negative values
    for ii=0,2*avecyl do begin
        if tempx(ii) lt 0 then tempx(ii)=N_elements(divv(*,0))+tempx(ii)
        if tempy(ii) lt 0 then tempy(ii)=N_elements(divv(0,*))+tempy(ii)
    endfor

    ;use periodic BC to fix outside grid values
    for ii=0,2*avecyl do begin
        if tempx(ii) gt n_elements(divv(*,0))-1 then tempx(ii)=tempx(ii)-N_elements(divv(*,0))
        if tempy(ii) gt n_elements(divv(0,*))-1 then tempy(ii)=tempy(ii)-N_elements(divv(0,*))
    endfor


    ;Interpolate the values normal to the shock NOTE: THES ARE NOT PARALLEL TO SHOCK YET!
    ronorm=shocknormvals(ro,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
    vxnorm=shocknormvals(vx,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
    vynorm=shocknormvals(vy,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
    prnorm=shocknormvals(pr,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)

    if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
    bxnorm=shocknormvals(bx,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
    bynorm=shocknormvals(by,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
    endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check for maximum density gradient along the normal
    ;rgnorm=dblarr(2*avecyl+1)
    ;for ii=1,2*avecyl-2 do begin
    ;    rgnorm(ii)=abs(ronorm(ii+1)-ronorm(ii-1))
    ;endfor

    ;if (max(rgnorm)-abs(rgnorm(avecyl))) lt 1.0e-6 then begin

    ;a=max(abs(deriv(ronorm)),b)
    ;if b eq avecyl then begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;I think these are the wrong way round
    ;Calculate the parallel and perpendicular velocity and magnetic field
;    vpar=normx*vxnorm+normy*vynorm
    ;vperp=sqrt(vxnorm^2+vynorm^2-vpar^2)
;    vperp=perpx*vxnorm+perpy*vynorm
;    bpar=normx*bxnorm+normy*bynorm
    ;bperp=sqrt(bxnorm^2+bynorm^2-bpar^2)
;    bperp=perpx*bxnorm+perpy*bynorm
;    ang=atan(bpar/bperp) ;DOUBLE CHECK THIS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;Calculate the parallel and perpendicular velocity and magnetic field
    vperp=normx*vxnorm+normy*vynorm
    ;vperp=sqrt(vxnorm^2+vynorm^2-vpar^2)
    vpar=perpx*vxnorm+perpy*vynorm

    if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
    bperp=normx*bxnorm+normy*bynorm
    ;bperp=sqrt(bxnorm^2+bynorm^2-bpar^2)
    bpar=perpx*bxnorm+perpy*bynorm
    ang=atan(bpar/bperp) ;DOUBLE CHECK THIS
    endif

    ;define pre and post shock indicies
    if ronorm(0) LT ronorm(2*avecyl) then begin 
;            ipre=0
;            ipos=2*avecyl
            rominpre=min(ronorm(0:avecyl-1),ipre)
            rominpos=max(ronorm(avecyl+1:avecyl*2),ipos)
            ipos=ipos+avecyl+1
    endif else begin
;            ipos=0
;            ipre=2*avecyl
            rominpos=max(ronorm(0:avecyl-1),ipos)
            rominpre=min(ronorm(avecyl+1:avecyl*2),ipre)
            ipre=ipre+avecyl+1
    endelse

    ;calculate conserved quantities
    con_mas=(ronorm(ipos)*vperp(ipos)-ronorm(ipre)*vperp(ipre))/(ronorm(ipre)-ronorm(ipos))


    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Estimate shock speed 
;    vs=(vpar(ipre)-vpar(ipos))/(1.0-ronorm(ipre)/ronorm(ipos))
    ;vs=(vperp(ipre)-vperp(ipos))/(1.0-ronorm(ipre)/ronorm(ipos))
    vs=con_mas


    ;Calculate Mach numbers
    cs2=5.0/3.0*prnorm/ronorm
    if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
    va2=(bxnorm^2+bynorm^2)/ronorm
    vap2=bperp^2/ronorm
    vslow2=0.5*(cs2+va2-sqrt((cs2+va2)^2 - 4.0*va2*cs2*(cos(ang))^2))
    vfast2=0.5*(cs2+va2+sqrt((cs2+va2)^2 - 4.0*va2*cs2*(cos(ang))^2))
    endif

;    p=plot((vpar-vs)^2)
;    p=plot(cs2,/overplot,'r')
;    p=plot(va2,/overplot,'b')
;    p=plot(vslow2,/overplot,'g')
;    p=plot(vfast2,/overplot,'m')

    vf=vperp

    vf=vperp+vs ;I think this is the correct one
    ;if vs LT 0.0 then vf=vperp-vs
    ;if vs GT 0.0 then vf=vperp+vs

    ;if abs(vperp(ipre)+vs) GT abs(vperp(ipre)-vs) then vf=vperp+vs
    ;if abs(vperp(ipre)-vs) GT abs(vperp(ipre)+vs) then vf=vperp-vs

    ;if ronorm(ipre) LT ronorm(ipos) then vf=vperp-vs
    ;if ronorm(ipos) LT ronorm(ipre) then vf=vperp+vs

    ;if min(abs(vperp)) LT vs then vf=vperp+vs
    ;if min(abs(vperp)) GT vs then vf=vperp-vs

    if data_save eq 1 then begin
    ;alocate data to perminant array
    arrronorm(i,*)=ronorm
    ;arrvxnorm(i,*)=vxnorm
    ;arrvynorm(i,*)=vynorm
    ;arrbxnorm(i,*)=bxnorm
    ;arrbynorm(i,*)=bynorm
    arrprnorm(i,*)=prnorm
    arrvpar(i,*)=vpar
    arrvperp(i,*)=vperp
    arrbpar(i,*)=bpar
    arrbperp(i,*)=bperp
    arrvslow2(i,*)=vslow2
    arrvfast2(i,*)=vfast2
    arrvalf2(i,*)=vap2
    arrvs(i)=vs
    ;arrnx(i)=normx
    ;arrny(i)=normy
    ;arrpx(i)=perpx
    ;arrpy(i)=perpy
    arrms(i,*)=abs(vf/sqrt(vslow2))
    arrma(i,*)=abs(vf/sqrt(vap2))
    arrmf(i,*)=abs(vf/sqrt(vfast2))
    endif
    ;Test if it is a shock
;    if (max(vf/sqrt(vslow2)) GT 1) and (min(vf/sqrt(vslow2)) LT 1) then begin
;        slowx=[slowx,col(i)]
;        slowy=[slowy,row(i)]
;    endif

;    if (max(vf/sqrt(vfast2)) GT 1) and (min(vf/sqrt(vfast2)) LT 1) then begin
;        fastx=[fastx,col(i)]
;        fasty=[fasty,row(i)]
;    endif 

    ;Test if it is a shock
;    if (vf(ipre)/sqrt(vslow2(ipre)) GT 1) and (vf(ipos)/sqrt(vslow2(ipos)) LT 1) then begin
;        slowx=[slowx,col(i)]
;        slowy=[slowy,row(i)]
;    endif

;    if (vf(ipre)/sqrt(vfast2(ipre)) GT 1) and (vf(ipos)/sqrt(vfast2(ipos)) LT 1) then begin
;        fastx=[fastx,col(i)]
;        fasty=[fasty,row(i)]
;    endif 
    
;    if (max(vf(ipre:avecyl)/sqrt(vslow2(ipre:avecyl))) GT 1) and (min(vf(avecyl+1:ipos)/sqrt(vslow2(avecyl+1:ipos))) LT 1) then begin
;        slowx=[slowx,col(i)]
;        slowy=[slowy,row(i)]
;    endif

;    if (max(vf(ipre:avecyl)/sqrt(vfast2(ipre:avecyl))) GT 1) and (min(vf(avecyl+1:ipos)/sqrt(vfast2(avecyl+1:ipos))) LT 1) then begin
;        fastx=[fastx,col(i)]
;        fasty=[fasty,row(i)]
;    endif 


;Clasification based on transitions
;    vslowpre =mean(abs(vf(ipre:avecyl-3)/sqrt(vslow2(ipre:avecyl-3))))
;    valfpre  =mean(abs(vf(ipre:avecyl-3)/sqrt(vap2(ipre:avecyl-3))))
;    vfastpre =mean(abs(vf(ipre:avecyl-3)/sqrt(vfast2(ipre:avecyl-3))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Plasma transisitons
if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
    vslowpre =abs(vf(ipre)/sqrt(vslow2(ipre)))
    valfpre  =abs(vf(ipre)/sqrt(vap2(ipre)))
    vfastpre =abs(vf(ipre)/sqrt(vfast2(ipre)))

    tpre=4
;    if vfastpre gt 1.0 then tpre=1 else $
;    if (valfpre gt 1.0) and (vfastpre lt 1.0) then tpre=2 else $
;    if (vslowpre gt 1.0) and (valfpre lt 1.0) then tpre=3 

    diffmin=1000.0
    if vfastpre gt (1.0-shocktol) then begin 
        tpre=1 
        if (vfastpre-1.0) lt shocktol then diffmin=abs(vfastpre-1.0)
    endif

    if (valfpre gt (1.0-shocktol)) and (vfastpre lt (1.0+shocktol)) then begin 
        if tpre eq 4 then begin
            tpre=2 
            if (valfpre-1.0) lt shocktol then diffmin=abs(valfpre-1.0)
        endif else begin
            alerr=abs(valfpre-1.0)
            if (alerr lt shocktol) and (alerr lt diffmin) then diffmin=alerr & tpre=2
        endelse
    endif

    if (vslowpre gt (1.0-shocktol)) and (valfpre lt (1.0+shocktol)) then begin 
        if tpre eq 4 then begin
            tpre=3 
            if (vslowpre-1.0) lt shocktol then diffmin=abs(vslowpre-1.0)
        endif else begin
            alerr=abs(vslowpre-1.0)
            if (alerr lt shocktol) and (alerr lt diffmin) then diffmin=alerr & tpre=3
        endelse
    endif

;    vslowpos =mean(abs(vf(avecyl+4:ipos)/sqrt(vslow2(avecyl+4:ipos))))
;    valfpos  =mean(abs(vf(avecyl+4:ipos)/sqrt(vap2(avecyl+4:ipos))))
;    vfastpos =mean(abs(vf(avecyl+4:ipos)/sqrt(vfast2(avecyl+4:ipos))))
    vslowpos =abs(vf(ipos)/sqrt(vslow2(ipos)))
    valfpos  =abs(vf(ipos)/sqrt(vap2(ipos)))
    vfastpos =abs(vf(ipos)/sqrt(vfast2(ipos)))

    tpos=4
;    if vfastpos gt 1.0 then tpos=1 else $
;    if (valfpos gt 1.0) and (vfastpos lt 1.0) then tpos=2 else $
;    if (vslowpos gt 1.0) and (valfpos lt 1.0) then tpos=3 

    diffmin=1000.0
    if vfastpos gt (1.0-shocktol) then begin 
        tpos=1 
        if (vfastpos-1.0) lt shocktol then diffmin=abs(vfastpos-1.0)
    endif

    if (valfpos gt (1.0-shocktol)) and (vfastpos lt (1.0+shocktol)) then begin 
        if tpos eq 4 then begin
            tpos=2 
            if (valfpos-1.0) lt shocktol then diffmin=abs(valfpos-1.0)
        endif else begin
            alerr=abs(valfpos-1.0)
            if (alerr lt shocktol) and (alerr lt diffmin) then diffmin=alerr & tpos=2
        endelse
    endif

    if (vslowpos gt (1.0-shocktol)) and (valfpos lt (1.0+shocktol)) then begin 
        if tpos eq 4 then begin
            tpos=3 
            if (vslowpos-1.0) lt shocktol then diffmin=abs(vslowpos-1.0)
        endif else begin
            alerr=abs(vslowpos-1.0)
            if (alerr lt shocktol) and (alerr lt diffmin) then diffmin=alerr & tpos=3
        endelse
    endif
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Neutral transitions
if (species eq 'neutral') and (ds.fl_pip eq 1) then begin

    sonicm=abs(vf/sqrt(cs2))

    tpos=2
    if sonicm(ipos) gt 1 then tpos = 1

    tpre=2
    if sonicm(ipre) gt 1 then tpre = 1

endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check for maximum density gradient along the normal
    rescode='Not local maximum density gradient'
    a=max(abs(deriv(ronorm)),b)
    if b eq avecyl then begin

    rescode='No transition satisfied'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Compressible transitions

;Plasma transisitons
if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
    if (tpre eq 1) and (tpos eq 2) then begin
        fastx=[fastx,col(i)]
        fasty=[fasty,row(i)]
        if data_save eq 1 then arrtype(i)=6
        rescode='fast'
    endif

    if (tpre eq 3) and (tpos eq 4) then begin
        slowx=[slowx,col(i)]
        slowy=[slowy,row(i)]
        if data_save eq 1 then arrtype(i)=1
        rescode='slow'
    endif 

    if (tpre eq 1) and (tpos eq 3) and (bpar(ipos)*bpar(ipre) LT 0) then begin
        int1x=[int1x,col(i)]
        int1y=[int1y,row(i)]
        if data_save eq 1 then arrtype(i)=4
        rescode='int'
    endif  

    if (tpre eq 1) and (tpos eq 4) and (bpar(ipos)*bpar(ipre) LT 0) then begin
        int2x=[int2x,col(i)]
        int2y=[int2y,row(i)]
        if data_save eq 1 then arrtype(i)=3
        rescode='int'
    endif  
 

    if (tpre eq 2) and (tpos eq 3) and (bpar(ipos)*bpar(ipre) LT 0) then begin
        int3x=[int3x,col(i)]
        int3y=[int3y,row(i)]
        if data_save eq 1 then arrtype(i)=5
        rescode='int'
    endif     

    if (tpre eq 2) and (tpos eq 4) and (bpar(ipos)*bpar(ipre) LT 0) then begin
        int4x=[int4x,col(i)]
        int4y=[int4y,row(i)]
        if data_save eq 1 then arrtype(i)=2
        rescode='int'
    endif 
endif

;Neutral transitions
if (species eq 'neutral') and (ds.fl_pip eq 1) then begin
    if (tpre eq 1) and (tpos eq 2) then begin
        slowx=[slowx,col(i)]
        slowy=[slowy,row(i)]
        if data_save eq 1 then arrtype(i)=1
        rescode='sonic'
    endif
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Expansive transitions (non-shock?)

;    if (tpre eq 2) and (tpos eq 1) then begin
;        fastx=[fastx,col(i)]
;        fasty=[fasty,row(i)]
;        arrtype(i)=6
;    endif

;    if (tpre eq 4) and (tpos eq 3) then begin
;        slowx=[slowx,col(i)]
;        slowy=[slowy,row(i)]
;        arrtype(i)=1
;    endif    

;    if (tpre eq 3) and (tpos eq 1) then begin
;        int1x=[int1x,col(i)]
;        int1y=[int1y,row(i)]
;        arrtype(i)=4
;    endif  

;    if (tpre eq 4) and (tpos eq 1) then begin
;        int2x=[int2x,col(i)]
;        int2y=[int2y,row(i)]
;        arrtype(i)=3
;    endif  

;    if (tpre eq 3) and (tpos eq 2) then begin
;        int3x=[int3x,col(i)]
;        int3y=[int3y,row(i)]
;        arrtype(i)=5
;    endif 

;    if (tpre eq 4) and (tpos eq 2) then begin
;        int4x=[int4x,col(i)]
;        int4y=[int4y,row(i)]
;        arrtype(i)=2
;    endif    

    endif

;    print,i,(n_elements(col)-1),rescode
endfor

if n_elements(slowx) eq 0 then slowx=[0]
if n_elements(slowy) eq 0 then slowy=[0]

if n_elements(fastx) eq 0 then fastx=[0]
if n_elements(fasty) eq 0 then fasty=[0]

if n_elements(int1x) eq 0 then int1x=[0]
if n_elements(int1y) eq 0 then int1y=[0]

if n_elements(int2x) eq 0 then int2x=[0]
if n_elements(int2y) eq 0 then int2y=[0]

if n_elements(int3x) eq 0 then int3x=[0]
if n_elements(int3y) eq 0 then int3y=[0]

if n_elements(int4x) eq 0 then int4x=[0]
if n_elements(int4y) eq 0 then int4y=[0]

shocks=create_struct(['x','y','slowx','slowy','fastx','fasty',$
'int1x','int1y','int2x','int2y','int3x','int3y','int4x','int4y'],$
x,y,slowx,slowy,fastx,fasty,int1x,int1y,int2x,int2y,int3x,int3y,int4x,int4y)

;c=image(divv,x,y,rgb_tab=0,xtitle='x',ytitle='y',dim=[970,600], position=[0.1,0.1,0.7,0.9],axis_style=2,xr=[0,1],yr=[0,1],/buffer)
;c=contour(divv,x,y,/fill,n_levels=101,rgb_tab=0,xtitle='x',ytitle='y',dim=[970,600], position=[0.1,0.1,0.7,0.9],xr=[0,1],yr=[0,1])
;diva=divv                     
;diva[where(divv gt 0)]=1.0/0.0
;diva[where(divv gt -convl)]=1.0/0.0
;c=image(alog10(-diva),x,y,rgb_tab=0,xtitle='x',ytitle='y',dim=[970,600], position=[0.1,0.1,0.7,0.9],axis_style=2,xr=[0,1],yr=[0,1])
;c=image(alog10(abs(gradmag)),x,y,rgb_tab=0,xtitle='x',ytitle='y',dim=[970,600], position=[0.1,0.1,0.8,0.9],axis_style=2,xr=[0,1],yr=[0,1])

;xshift=5.0
;yshift=5.0

;if N_elements(slowx) gt 0 then $
;ps=plot(1.0*slowx/(1.0*n_elements(x)-xshift),1.0*slowy/(1.0*n_elements(y)-yshift),/overplot,'.r',name='slow',thick=0.1,sym_thick=2)
;if N_elements(fastx) gt 0 then $
;pf=plot(1.0*fastx/(1.0*n_elements(x)-xshift),1.0*fasty/(1.0*n_elements(y)-yshift),/overplot,'.b',name='fast',thick=0.1,sym_thick=2)
;if N_elements(int1x) gt 0 then $
;pi1=plot(1.0*int1x/(1.0*n_elements(x)-xshift),1.0*int1y/(1.0*n_elements(y)-yshift),/overplot,'.y',name='Int (1-3)',thick=0.1,sym_thick=2)
;if N_elements(int2x) gt 0 then $
;pi2=plot(1.0*int2x/(1.0*n_elements(x)-xshift),1.0*int2y/(1.0*n_elements(y)-yshift),/overplot,'.m',name='Int (1-4)',thick=0.1,sym_thick=2)
;if N_elements(int3x) gt 0 then $
;pi3=plot(1.0*int3x/(1.0*n_elements(x)-xshift),1.0*int3y/(1.0*n_elements(y)-yshift),/overplot,'.c',name='Int (2-3)',thick=0.1,sym_thick=2)
;if N_elements(int4x) gt 0 then $
;pi4=plot(1.0*int4x/(1.0*n_elements(x)-xshift),1.0*int4y/(1.0*n_elements(y)-yshift),/overplot,'.g',name='Int (2-4)',thick=0.1,sym_thick=2)

;to get the labels right:
;psl=plot([-1,-1],[-1,-1],/overplot,'r',name='slow (3-4)',thick=2)
;pfl=plot([-1,-1],[-1,-1],/overplot,'b',name='fast (1-2)',thick=2)
;pi1l=plot([-1,-1],[-1,-1],/overplot,'y',name='Int (1-3)',thick=2)
;pi2l=plot([-1,-1],[-1,-1],/overplot,'m',name='Int (1-4)',thick=2)
;pi3l=plot([-1,-1],[-1,-1],/overplot,'c',name='Int (2-3)',thick=2)
;pi4l=plot([-1,-1],[-1,-1],/overplot,'g',name='Int (2-4)',thick=2)

;l=legend(target=[psl,pfl,pi1l,pi2l,pi3l,pi4l])
;l.position=[0.9,0.9]

;Annotate the number of shocks
;ts=text(0.75,0.6,'Slow = '+strtrim(n_elements(slowx),1),/overplot)
;ts=text(0.75,0.55,'Fast = '+strtrim(n_elements(fastx),1),/overplot)
;ts=text(0.75,0.5,'Int (1-3) = '+strtrim(n_elements(int1x),1),/overplot)
;ts=text(0.75,0.45,'Int (1-4) = '+strtrim(n_elements(int2x),1),/overplot)
;ts=text(0.75,0.4,'Int (2-3) = '+strtrim(n_elements(int3x),1),/overplot)
;ts=text(0.75,0.35,'Int (2-4) = '+strtrim(n_elements(int4x),1),/overplot)

;if data_subset eq 1 then begin
;c.xr=[xminf,xmaxf]
;c.yr=[yminf,ymaxf]
;pf.sym_thick=4
;ps.sym_thick=4
;endif

;c.save,'shockid_PIP_MHD_t'+strtrim(tread,1)+'.pdf'
;c.save,'shockid_PIP_ac_100_t'+strtrim(tread,1)+'_'+species+'_bulk.pdf'

if data_save eq 1 then begin
save,arrronorm,arrvxnorm,arrvynorm,arrbxnorm,arrbynorm,$
arrprnorm,arrvpar,arrvperp,arrvslow2,arrvalf2,arrvfast2,$
arrvs,arrtype,col,row,filename='shockid_v7_data.dat'
;Order the perminant arrays
arrronorms=arrronorm(sort(arrtype),*)
arrvxnorms=arrvxnorm(sort(arrtype),*)
arrvynorms=arrvynorm(sort(arrtype),*)
arrbxnorms=arrbxnorm(sort(arrtype),*)
arrbynorms=arrbynorm(sort(arrtype),*)
arrprnorms=arrprnorm(sort(arrtype),*)
arrvpars=arrvpar(sort(arrtype),*)
arrvperps=arrvperp(sort(arrtype),*)
arrbpars=arrbpar(sort(arrtype),*)
arrbperps=arrbperp(sort(arrtype),*)
arrvslow2s=arrvslow2(sort(arrtype),*)
arrvfast2s=arrvfast2(sort(arrtype),*)
arrvalf2s=arrvalf2(sort(arrtype),*)
arrvss=arrvs(sort(arrtype),*)
arrtypes=arrtype(sort(arrtype))
endif

;p=plot(1.0*col/(n_elements(x)-1.0),1.0*row/(n_elements(y)-1.0),/overplot,'.')

close,/all

END
