pro check_airglow_lines

;; Use the end-of-night cals and last science exposure to
;; check the wavelengths of the best and brightest airglow lines

apsetver,vers=‘t9’,telescope=‘apo25m’
dirs = getdir(apogee_dir,cal_dir,spectro_dir,apogee_vers,lib_dir)
linelist_dir = lib_dir+'skylines/'
wdir = '/uufs/chpc.utah.edu/common/home/sdss06/apogeework/apogee/spectro/redux/t9/cal/wave/'
expdir = '/uufs/chpc.utah.edu/common/home/sdss06/apogeework/apogee/spectro/redux/t9/exposures/apogee-n/'
visitdir = '/uufs/chpc.utah.edu/common/home/sdss06/apogeework/apogee/spectro/redux/t9/visit/apo25m/'

airglow_all = importascii(linelist_dir+'airglow.txt',/header)
;goodid = [6,8,17,24,42,61,62,63,86,94,118,181,191,205,239,285,291,303,350,363]
goodid = [6,8,17,24,42,61,62,63,71,86,94,99,118,181,191,205,239,270,275,285,291,303,350,363]
match,airglow_all.id,goodid,ind1,ind2,/sort
airglow = airglow_all[ind1]
jonair = importascii(linelist_dir+'airglow.new',/header)


daynum = [2231, 2239, 2244, 2245, 2247, 2255]
ndaynum = n_elements(daynum)
plates = [9684, 9684, 9050, 9050, 8443, 8309]
fields = ['330+60', '330+60', '040+43_MGA', '040+43_MGA', '029+77_MGA', '102+61_MGA']

;; Loop over days
undefine,lines
for i=0,ndaynum-1 do begin
  idaynum = daynum[i]
  cmjd = getcmjd(long(strtrim(idaynum,2)+'0001'),mjd=mjd)
  plate = plates[i]
  field = fields[i]

  ;; Get wavelength solution
  wfile = file_search(wdir+'apWave-'+strtrim(idaynum,2)+'????.dat',count=nwfile)
  wframe = strmid(file_basename(wfile[0],'.dat'),7)
  print,strtrim(i+1,2),' ',strtrim(idayname,2),' ',wframe
  fits_read,wdir+'apWave-a-'+wframe+'.fits',wcoef1,exten=1
  fits_read,wdir+'apWave-a-'+wframe+'.fits',wim1,exten=2
  fits_read,wdir+'apWave-b-'+wframe+'.fits',wcoef2,exten=1
  fits_read,wdir+'apWave-b-'+wframe+'.fits',wim2,exten=2
  fits_read,wdir+'apWave-c-'+wframe+'.fits',wcoef3,exten=1
  fits_read,wdir+'apWave-c-'+wframe+'.fits',wim3,exten=2

  ;; Load the plan and plugmap files
  planfile = visitdir+field+'/'+strtrim(plate,2)+'/'+cmjd+'/apPlan-'+strtrim(plate,2)+'-'+cmjd+'.par'
  aploadplan,planfile,planstr
  plugfile = planstr.plugmap
  plugmap = getplatedata(planstr.plateid,string(planstr.mjd,format='(i5.5)'),plugid=planstr.plugmap,$
                         fixfiberid=fixfiberid,badfiberid=badfiberid,mapper_data=mapper_data)
  ;; Load the 1D spectra
  nframes = n_elements(planstr.apexp)
  ind = first_el(maxloc(long(planstr.apexp.name)))
  aploadframe,expdir+cmjd+'/ap1D-'+planstr.apexp[ind].name,frame
  frame.(0).wcoef = wcoef1
  frame.(0).wave = wim1
  frame.(1).wcoef = wcoef2
  frame.(1).wave = wim2
  frame.(2).wcoef = wcoef3
  frame.(2).wave = wim3
  ;; Get the sky fibers
  skyind = where(plugmap.fiberdata.spectrographid eq 2 and plugmap.fiberdata.objtype eq 'SKY',nskyind)
  skyindex = 300-plugmap.fiberdata[skyind].fiberid
  appeakfit,frame.(0),linestr1,fibers=skyindex
  add_tag,linestr1,'chipnum',1,linestr1
  add_tag,linestr1,'wave',0.0d0,linestr1
  for j=0,n_elements(linestr1)-1 do linestr1[j].wave=pix2wave(linestr1[j].gaussx,reform(wcoef1[linestr1[j].fiber,*]))
  appeakfit,frame.(1),linestr2,fibers=skyindex
  add_tag,linestr2,'chipnum',2,linestr2
  add_tag,linestr2,'wave',0.0d0,linestr2
  for j=0,n_elements(linestr2)-1 do linestr2[j].wave=pix2wave(linestr2[j].gaussx,reform(wcoef2[linestr2[j].fiber,*]))
  appeakfit,frame.(2),linestr3,fibers=skyindex
  add_tag,linestr3,'chipnum',3,linestr3
  add_tag,linestr3,'wave',0.0d0,linestr3
  for j=0,n_elements(linestr3)-1 do linestr3[j].wave=pix2wave(linestr3[j].gaussx,reform(wcoef3[linestr3[j].fiber,*]))
  linestr = [linestr1,linestr2,linestr3]
  push,lines,linestr

  print,'Saving to ',planstr.apexp[ind].name+'_airglow.fits'
  MWRFITS,linestr,planstr.apexp[ind].name+'_airglow.fits',/create
  ;save,linestr,file=planstr.apexp[ind].name+'_airglow.dat'

  ;stop

endfor
MWRFITS,lines,'check_airglow_lines.fits',/create


;; Match to "good" airglow linelist
undefine,final
add_tag,airglow,'wmed',0.0d0,airglow
add_tag,airglow,'wsig',0.0,airglow
add_tag,airglow,'wrms',0.0,airglow
add_tag,airglow,'werr',0.0,airglow
print,'   ID  NLINES  WORIG       WMED      WDIFF   WSIG   WRMS    WERR    WDIFF2'
for i=0,n_elements(airglow)-1 do begin
  ind = where(abs(airglow[i].wave-lines.wave) lt 0.1,nind)
  if nind gt 0 then begin
    newlines = lines[ind]
    add_tag,newlines,'AIRID',airglow[i].id,newlines
    push,final,newlines
    ;; Calculate median and sigma
    airglow[i].wmed = median(newlines.wave)
    airglow[i].wsig = mad(newlines.wave)
    airglow[i].wrms = stddev(newlines.wave)
    airglow[i].werr = airglow[i].wsig/sqrt(nind)
    ;; Compare to Jon's airglow.new
    match,airglow[i].id,jonair.id,ind3,ind4,count=njmatch
    if njmatch gt 0 then begin
      jwave = jonair[ind4].wave
      jdiff = string(airglow[i].wave-jwave,format='(F8.4)')
    endif else begin
      jwave = 0.0d0
      jdiff = '  ------'
    endelse
    print,airglow[i].id,nind,airglow[i].wave,airglow[i].wmed,airglow[i].wave-airglow[i].wmed,airglow[i].wsig,airglow[i].wrms,airglow[i].werr,jdiff,format='(I5,I5,2F12.4,4F8.4,A8)'
  endif
endfor

MWRFITS,final,'check_airglow_lines_matched.fits',/create
MWRFITS,airglow,'check_airglow_linelist.fits',/create

stop


end
