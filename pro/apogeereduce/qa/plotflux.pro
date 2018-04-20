pro plotflux,planfile

; procedure to make plots of relative fluxes of fibers in flat field exposures

; load planfile
APLOADPLAN,planfile,planstr,/verbose,error=planerror
dirs=getdir()

; get plugmap
plugfile = planstr.plugmap
if tag_exist(planstr,'fixfiberid') then fixfiberid=planstr.fixfiberid
if tag_exist(planstr,'badfiberid') then badfiberid=planstr.badfiberid
plugmap=getplatedata(planstr.plateid,string(planstr.mjd,format='(i5.5)'),plugid=planstr.plugmap,fixfiberid=fixfiberid,badfiberid=badfiberid,mapper_data=mapper_data)

; read apFlux file
flux=apread('Flux',num=planstr.fluxid)

; plot setup
A = FINDGEN(17) * (!PI*2/16.)
usersym,cos(A),sin(A),/fill
thick=3
!p.thick=thick
!P.CharThick =thick
!X.Thick =thick
!Y.Thick =thick
!Z.Thick =thick

; output directory
platedir=apogee_filename('Plate',plate=planstr.plateid,mjd=planstr.mjd,chip='a',/dir)
outdir=platedir+'/plots/'
file=file_basename(apogee_filename('Flux',num=planstr.fluxid,chip=['a','b','c'],/base),'.fits')
set_plot,'PS'
loadct,39
ypos=300-plugmap.fiberdata.fiberid
if dirs.telescope eq 'lco25m' then lim=[-1.2,1.2] else lim=[-1.6,1.6]

; make the plots
for ichip=0,2 do begin
  med=median(flux[ichip].flux,dim=1)
  device,file=outdir+file[ichip]+'.eps',/encap,xsize=20,ysize=24,/color
  plotc,plugmap.fiberdata.zeta,plugmap.fiberdata.eta,med[ypos],min=0.5,max=1.5,psym=8,xr=lim,yr=lim,xstyle=1,ystyle=1,xtit='Zeta',ytit='Eta'
  ;for i=0,n_elements(plugmap.fiberdata)-1 do begin
  ;  xyouts,plugmap.fiberdata[i].zeta,plugmap.fiberdata[i].eta,string(format='(i3.3)',plugmap.fiberdata[i].fiberid),alignment=0.0,charsize=0.5
  ;endfor
  device,/close
  ps2jpg,outdir+file[ichip]+'.eps',/eps,chmod='664'o,/delete
endfor

longlink=[5,6,7,8,2,4,11,10,12,13]
block=fix((plugmap.fiberdata.fiberid-1)/30)+1
blockfile=file[0].replace('-a-','-block-')
device,file=outdir+blockfile+'.eps',/encap,xsize=20,ysize=24,/color
plotc,plugmap.fiberdata.zeta,plugmap.fiberdata.eta,block,min=0,max=10,psym=8,xr=lim,yr=lim,xstyle=1,ystyle=1,xtit='Zeta',ytit='Eta'
device,/close
ps2jpg,outdir+blockfile+'.eps',/eps,chmod='664'o,/delete


end

