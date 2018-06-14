pro speclib_bundle,planfile,noflux=noflux,nonorm=nonorm,vmicro=vmicro,savenorm=savenorm

; read the plan file to tell us what files to bundle
aploadplan,planfile,planstr
if tag_exist(planstr,'npart') then npart=planstr.npart else npart=30
if tag_exist(planstr,'dopart') then dopart=planstr.dopart else dopart=-1
if tag_exist(planstr,'npca') then npca=planstr.npca else npca=30
if tag_exist(planstr,'addnoise') then addnoise=planstr.addnoise else addnoise=0
if tag_exist(planstr,'noflux') then noflux=planstr.noflux else noflux=1

; are we just doing one "part" of the spectrum
if dopart ge 0 then begin
  ipart1=dopart
  ipart2=dopart
  if dopart gt 0 then noflux=1
endif else begin
  ipart1=0
  ipart2=npart
endelse

; get file names. Prefix is for input files
outfile=planstr.name
if tag_exist(planstr,'elem') then outfile=outfile+'_'+planstr.elem
if tag_exist(planstr,'prefix') then prefix=planstr.prefix else prefix=''

; if we have vmicro coefficients and nvt>1, then we will use ferre to do the interpolation
; if we have vmicro coefficients and nvt=1, then we will use previously interpolated spectra
if tag_exist(planstr,'vmicrofit') then vmicrofit=planstr.vmicrofit else $
   if tag_exist(planstr,'vmicro') then vmicrofit=1
if tag_exist(planstr,'vmicro') then vcoef=planstr.vmicro
if tag_exist(planstr,'vmsuffix') then vmsuffix='_'+planstr.vmsuffix else vmsuffix='vfit'
ppre='p'
if keyword_set(vmicrofit) then if planstr.nvt gt 1 then ppre='p6' else ppre='pv'

; determine total number of dimensions of output grid
ndim=0
if planstr.ncm gt 1 then ndim+=1
if planstr.nnm gt 1 then ndim+=1
if planstr.nam gt 1 then ndim+=1
if planstr.nmh gt 1 then ndim+=1
if planstr.nlogg gt 1 then ndim+=1
if planstr.nteff gt 1 then ndim+=1
if planstr.nvt gt 1 then ndim+=1
if tag_exist(planstr,'twod') then twod=planstr.twod else twod=0
if tag_exist(planstr,'elem') then ndim+=1
if tag_exist(planstr,'nrot') then nrot=planstr.nrot else nrot=1
if nrot gt 1 then ndim+=1

; get the dimension names and limits
dim=strarr(ndim)
llimits=fltarr(ndim)
steps=fltarr(ndim)
nsteps=intarr(ndim)
idim=0
if planstr.nvt gt 1 then begin
  nvt=planstr.nvt
  llimits[idim]=planstr.vt0
  steps[idim]=planstr.dvt
  nsteps[idim]=planstr.nvt
  dim[idim++]='LOG10VDOP'
endif else nvt=1
if planstr.ncm gt 1 then begin
  ncm=planstr.ncm
  llimits[idim]=planstr.cm0
  steps[idim]=planstr.dcm
  nsteps[idim]=planstr.ncm
  dim[idim++]='C'
endif else ncm=1
if planstr.nnm gt 1 then begin
  nnm=planstr.nnm
  llimits[idim]=planstr.nm0
  steps[idim]=planstr.dnm
  nsteps[idim]=planstr.nnm
  dim[idim++]='N'
endif else nnm=1
if planstr.nam gt 1 then begin
  nam=planstr.nam
  llimits[idim]=planstr.am0
  steps[idim]=planstr.dam
  nsteps[idim]=planstr.nam
  dim[idim++]='O Mg Si S Ca Ti'
endif else nam=1
if tag_exist(planstr,'nrot') then begin
  if planstr.nrot gt 1 then begin
    llimits[idim]=planstr.rot0
    steps[idim]=planstr.drot
    nsteps[idim]=planstr.nrot
    dim[idim++]='LGVSINI'
  endif
endif
if tag_exist(planstr,'elem') then begin
  ;llimits[idim]=-0.75
  ;steps[idim]=0.25
  llimits[idim]=-0.75
  steps[idim]=0.25
  nsteps[idim]=10
  dim[idim++]=planstr.elem
endif
if planstr.nmh gt 1 then begin
  nmh=planstr.nmh
  llimits[idim]=planstr.mh0
  steps[idim]=planstr.dmh
  nsteps[idim]=planstr.nmh
  dim[idim++]='METALS'
endif else nmh=1
if planstr.nlogg gt 1 then begin
  nlogg=planstr.nlogg
  llimits[idim]=planstr.logg0
  steps[idim]=planstr.dlogg
  nsteps[idim]=planstr.nlogg
  dim[idim++]='LOGG'
endif else nlogg=1
if planstr.nteff gt 1 then begin
  nteff=planstr.nteff
  llimits[idim]=planstr.teff0
  steps[idim]=planstr.dteff
  nsteps[idim]=planstr.nteff
  dim[idim++]='TEFF'
endif else nteff=1

; read the file that gives which models are filled holes
if tag_exist(planstr,'holefile') then begin
  readcol,getenv('APOGEE_REDUX')+'/speclib/'+planstr.holefile,$
    mhhole,cmhole,amhole,teffhole,logghole,modelhole,format='(f,f,f,f,f,i)'
  gridhole=bytarr(nvt,ncm,nnm,nam,nmh,nlogg,nteff)
endif

; open output files
dir=getenv('APOGEE_LOCALDIR')
if dir ne '' then begin
  print,'mktemp -d '+dir+'/XXXXXX'
  spawn,'mktemp -d '+dir+'/XXXXXX',temp
  dir=temp[0]
  dir=dir+'/'
endif

nlun=indgen(npart+1)+20
for ipart=ipart1,ipart2 do begin
  if ipart eq 0 then $
    openw,nlun[ipart],dir+'n_aps'+outfile+'.dat' $
  else begin
    openw,nlun[ipart],dir+'n_aps'+outfile+'_'+string(ipart,format='(i2.2)')+'.hdr'
    openw,nlun[ipart]+npart,dir+'n_aps'+outfile+'_'+string(ipart,format='(i2.2)')+'.dat'
  endelse
endfor
if ipart1 eq 0 then begin
  openw,plun,dir+ppre+'_aps'+outfile+'_w123.dat',/get_lun
  openw,flun,dir+'f_'+outfile+'.dat',/get_lun
endif

if twod then nmh=planstr.nmh else nmh=1

; loop over all of the input files and append the data to the output files
nmod=0L
for ivt=0,nvt-1 do begin
 vt=planstr.vt0+ivt*planstr.dvt
 for icm=0,ncm-1 do begin
  cm=planstr.cm0+icm*planstr.dcm  
  for inm=0,nnm-1 do begin
    nm=planstr.nm0+inm*planstr.dnm  
    for iam=0,nam-1 do begin
     am=planstr.am0+iam*planstr.dam  
     for imh=0,nmh-1 do begin
      mh=planstr.mh0+imh*planstr.dmh  

      ; construct the input file name given the parameters
      if twod then  $
      name='m'+cval(mh)+'a'+cval(am)+'c'+cval(cm)+'n'+cval(nm) else $
      name='a'+cval(am)+'c'+cval(cm)+'n'+cval(nm)
      if planstr.nvt eq 1 and keyword_set(vmicrofit) then name=name+vmsuffix else $
      if tag_exist(planstr,'nvt') then name=name+'v'+cval(10^vt)
      name=prefix+name
      if tag_exist(planstr,'elem') then name=name+planstr.elem
      print,name

      ; read in the file
      if file_test(name+'.fits.fz') then spawn,['funpack','-O',dir+name+'.fits',name+'.fits.fz'],/noshell
      suffix='.fits'
      hdrraw=headfits(dir+name+suffix)
      spec=mrdfits(dir+name+suffix,0,hdr)
      s0=mrdfits(dir+name+suffix,1,hdr0)
      s1=mrdfits(dir+name+suffix,2,hdr1)
      s2=mrdfits(dir+name+suffix,3,hdr2)
      if file_test(name+'.fits.fz') then file_delete,dir+name+suffix

      npix=[sxpar(hdr0,'NAXIS1'),sxpar(hdr1,'NAXIS1'),sxpar(hdr2,'NAXIS1')]
      w0=[sxpar(hdr0,'CRVAL1'),sxpar(hdr1,'CRVAL1'),sxpar(hdr2,'CRVAL1')]
      dw=[sxpar(hdr0,'CDELT1'),sxpar(hdr1,'CDELT1'),sxpar(hdr2,'CDELT1')]

      ; if we are doing an element, read in the data for the elemental variation      
      ; This will need to be combined with the data from the nominal abundance
      ;if tag_exist(planstr,'elem') then begin
      ;  espec=mrdfits(name+'_'+planstr.elem+'.fits',0,hdr)
      ;  e0=mrdfits(name+'_'+planstr.elem+'.fits',1,hdr0)
      ;  e1=mrdfits(name+'_'+planstr.elem+'.fits',2,hdr1)
      ;  e2=mrdfits(name+'_'+planstr.elem+'.fits',3,hdr2)
      ;endif 

      ; check that file dimensions match requested grid dimensions.
      ;  in princple, might want to be able to make subgrids ... but not yet.
      sz=size(spec,/dim)
      sz0=size(s0,/dim)
      sz1=size(s1,/dim)
      sz2=size(s2,/dim)
      if sz0[1] ne planstr.nteff or sz0[2] ne planstr.nlogg then $
        stop,'dimensions of file do not match grid dimensions'+name

      ; set the output dimensions
      ; always have output npix, Teff, logg, mh, elem, rot
      outdim=intarr(5)
      outdim[0]=sz0[1]
      outdim[1]=sz0[2]
      if twod then outdim[2]=1 else outdim[2]=sz0[3]
      if tag_exist(planstr,'elem') then outdim[3]=sz0[4] else outdim[3]=1
      if nrot gt 1 then outdim[4]=nrot else outdim[4]=1
   
      ; if this is an element grid, need to integrate the default 0. abundance
      ; into the -0.75, -0.5, -0.25. .... 0.25, 0.5, 0.75 of elem grid 
      ;if tag_exist(planstr,'elem') then begin
      ;  tmp=fltarr([sz[0],outdim])
      ;  tmp[*,*,*,*,0:2,*]=espec[*,*,*,*,0:2,0:nrot-1]
      ;  tmp[*,*,*,*,3,*]=spec[*,*,*,*,0:nrot-1]
      ;  tmp[*,*,*,*,4:6,*]=espec[*,*,*,*,3:5,0:nrot-1]
      ;  espec=tmp
;
;        tmp=fltarr([sz0[0],outdim])
;        tmp[*,*,*,*,0:2,*]=e0[*,*,*,*,0:2,0:nrot-1]
;        tmp[*,*,*,*,3,*]=s0[*,*,*,*,0:nrot-1]
;        tmp[*,*,*,*,4:6,*]=e0[*,*,*,*,3:5,0:nrot-1]
;        e0=tmp
;
;        tmp=fltarr([sz1[0],outdim])
;        tmp[*,*,*,*,0:2,*]=e1[*,*,*,*,0:2,0:nrot-1]
;        tmp[*,*,*,*,3,*]=s1[*,*,*,*,0:nrot-1]
;        tmp[*,*,*,*,4:6,*]=e1[*,*,*,*,3:5,0:nrot-1]
;        e1=tmp
;
;        tmp=fltarr([sz2[0],outdim])
;        tmp[*,*,*,*,0:2,*]=e2[*,*,*,*,0:2,0:nrot-1]
;        tmp[*,*,*,*,3,*]=s2[*,*,*,*,0:nrot-1]
;        tmp[*,*,*,*,4:6,*]=e2[*,*,*,*,3:5,0:nrot-1]
;        e2=tmp
;
;        spec=reform(espec[*,*,*,*,*],[sz[0],outdim[0:3]])
;        s0=reform(e0[*,*,*,*,*,0:nrot-1],[sz0[0],outdim])
;        s1=reform(e1[*,*,*,*,*,0:nrot-1],[sz1[0],outdim])
;        s2=reform(e2[*,*,*,*,*,0:nrot-1],[sz2[0],outdim])
;      endif else begin
        ; otherwise reform into correct dimensions
        spec=reform(spec[*,*,*,*,*],[sz[0],outdim[0:3]])
        s0=reform(s0[*,*,*,*,*,0:nrot-1],[sz0[0],outdim])
        s1=reform(s1[*,*,*,*,*,0:nrot-1],[sz1[0],outdim])
        s2=reform(s2[*,*,*,*,*,0:nrot-1],[sz2[0],outdim])
;      endelse
      out=[s0[*,0,0,0,0,0],s1[*,0,0,0,0,0],s2[*,0,0,0,0,0]]
      szout=size(out,/dim)
      nout=szout[0]/npart
      ipca=intarr(npart)
      for ipart=0,npart-2 do ipca[ipart]=nout
      ipca[npart-1]=szout[0]-total(ipca)

near=fltarr(8,3)
near[0,*]=[0,-1,0]
near[1,*]=[0, 1,0]
near[2,*]=[0,0,-1]
near[3,*]=[0,0,1]
near[4,*]=[0,-1,-1]
near[5,*]=[0,-1,1]
near[6,*]=[0,1,-1]
near[7,*]=[0,1,1]
      if tag_exist(planstr,'holefile') then curhole=where(cmhole eq cm and amhole eq am)

      ; loop over all of the dimensions, writing out the spectra!
      for irot=0,outdim[4]-1 do begin
       rot = planstr.rot0+irot*planstr.drot
       for ielem=0,outdim[3]-1 do begin
        for jmh=0,outdim[2]-1 do begin
         if ~twod then mh=planstr.mh0+jmh*planstr.dmh  
         for ilogg=0,outdim[1]-1 do begin
          logg=planstr.logg0+ilogg*planstr.dlogg  
          for iteff=0,outdim[0]-1 do begin
           teff=planstr.teff0+iteff*planstr.dteff  
            if tag_exist(planstr,'holefile') then begin
              jhole=where(teffhole[curhole] eq teff and logghole[curhole] eq logg and mhhole[curhole] eq mh,nhole)
              if nhole eq 0 then gridhole[ivt,icm,inm,iam,jmh,ilogg,iteff]=0 else $
                 gridhole[ivt,icm,inm,iam,jmh,ilogg,iteff]=modelhole[curhole[jhole]]
              if spec[0,iteff,ilogg,jmh,ielem] lt 0. then $
                 gridhole[ivt,icm,inm,iam,jmh,ilogg,iteff]=3
;
;              if max(spec[*,iteff,ilogg,jmh,ielem]) lt 0.001 then begin
;                 gridhole[ivt,icm,inm,iam,jmh,ilogg,iteff]=3
;                 in=0
;                 while max(spec[*,iteff,ilogg,jmh,ielem]) lt 0.001 and in lt 7 do begin
;                 if iteff+near[in,0] ge 0 and iteff+near[in,0] lt outdim[0] and $
;                    ilogg+near[in,1] ge 0 and ilogg+near[in,1] lt outdim[1] and $
;                    jmh+near[in,2] ge 0 and jmh+near[in,2] lt outdim[2] then $
;                   spec[*,iteff,ilogg,jmh,ielem]=$
;                     spec[*,iteff+near[in,0],ilogg+near[in,1],jmh+near[in,2],ielem]
;                   in+=1
;                 endwhile
;                 gridhole[ivt,icm,inm,iam,jmh,ilogg,iteff]=3+in-1
;              endif
            endif

            ; write the convolved spectrum
            if ipart1 eq 0 and irot eq 0 then begin
              ym=double(spec[*,iteff,ilogg,jmh,ielem])
              if nmod eq 0 then begin
                printf,flun," &SYNTH"
                printf,flun," SYNTHFILE_INTERNAL = '"+outfile+"'"
                printf,flun," ID = '"+outfile+"'"
                printf,flun," DATE = '"+systime(0)+"'"
                printf,flun," N_OF_DIM = ",string(n_elements(dim))
                printf,flun," N_P =  ",string(nsteps)
                for i=0,n_elements(dim)-1 do $
                 printf,flun,' LABEL('+string(i+1,format='(i1)')+') = '+$
                          "'"+dim[i]+"'"
                printf,flun,' LLIMITS =  '+string(llimits,format='(30f11.5)')
                printf,flun,' STEPS =  '+string(steps,format='(30f11.5)')
                printf,flun," NPIX = "+string(n_elements(sz[0]))
                printf,flun," WAVE =  "+string(sxpar(hdr,'CRVAL1'))+$
                             ' '+string(sxpar(hdr,'CDELT1'))
                printf,flun," LOGW =  2 "
                printf,flun," VACUUM =  1 "
                printf,flun," RESOLUTION =  "+string(sxpar(hdr,'RESOLUTI'))
                printf,flun," ORIGINAL_SAMPLING=  "+string(sxpar(hdrraw,'CDELT1'))
                if sxpar(hdr,'LINELIST') then $
                  printf,flun," FILE_DATA19= '"+sxpar(hdr,"LINELIST")+"'"
                printf,flun," COMMENTS1 = "
                printf,flun," /"
              endif
              ; write out the spectrum!
              if ~keyword_set(noflux) then $
                printf,flun,ym,format='(100000(e13.5,1x))'
            endif

            ; continuum normalized, multi-file with header for each chip
            ;if ipart1 eq 0 and nmod eq 0 then begin
            if nmod eq 0 then begin
             for ipart=ipart1,ipart2 do begin
              printf,nlun[ipart]," &SYNTH"
              printf,nlun[ipart]," ID = '"+ppre+'_aps'+outfile+'_w123.dat'+"'"
              printf,nlun[ipart]," DATE = '"+systime(0)+"'"
              printf,nlun[ipart]," MULTI = 3"
              printf,nlun[ipart]," N_OF_DIM = ",string(n_elements(dim))
              printf,nlun[ipart]," N_P =  ",string(nsteps)
              for i=0,n_elements(dim)-1 do $
               printf,nlun[ipart],' LABEL('+string(i+1,format='(i1)')+') = '+$
                          "'"+dim[i]+"'"
              printf,nlun[ipart],' LLIMITS =  '+string(llimits,format='(30f11.5)')
              printf,nlun[ipart],' STEPS =  '+string(steps,format='(30f11.5)')
              printf,nlun[ipart]," /"
              if ipart eq 0 then begin
               printf,plun," &SYNTH"
               printf,plun," MULTI = 3"
               printf,plun," ID = '"+ppre+'_aps'+outfile+'_w123.dat'+"'"
               printf,plun," DATE = '"+systime(0)+"'"
               printf,plun,' NPCA = '+strjoin(string(ipca,format='(i6)'),' ')
               printf,plun,' NPIX = '+string(npca*npart,format='(i6)')
               if keyword_set(vmicrofit) and planstr.nvt gt 1 then begin
                 printf,plun," N_OF_DIM = ",string(n_elements(dim)-1)
                 printf,plun," N_P =  ",string(nsteps[1:*])
                 for i=1,n_elements(dim)-1 do $
                  printf,plun,' LABEL('+string(i,format='(i1)')+') = '+$
                             "'"+dim[i]+"'"
                 printf,plun,' LLIMITS =  '+string(llimits[1:*],format='(30f11.5)')
                 printf,plun,' STEPS =  '+string(steps[1:*],format='(30f11.5)')
               endif else begin
                 printf,plun," N_OF_DIM = ",string(n_elements(dim))
                 printf,plun," N_P =  ",string(nsteps)
                 for i=0,n_elements(dim)-1 do $
                  printf,plun,' LABEL('+string(i+1,format='(i1)')+') = '+$
                             "'"+dim[i]+"'"
                 printf,plun,' LLIMITS =  '+string(llimits,format='(30f11.5)')
                 printf,plun,' STEPS =  '+string(steps,format='(30f11.5)')
               endelse
               printf,plun," /"
              endif
             endfor
            endif
            ; loop over each of the three chips and write the individual headers
            ;if ipart1 eq 0 and nmod eq 0 then begin
            if nmod eq 0 then begin
             for ipart=ipart1,ipart2 do begin
              if ipart eq 0 then nloop=1 else nloop=0
              for ilun=0,nloop do begin
               if ipart eq 0 then if ilun eq 0 then lun=nlun[0] else lun=plun
               if ipart ne 0 then lun=nlun[ipart]
  	       for j=0,n_elements(npix)-1 do begin
                printf,lun," &SYNTH"
                printf,lun," SYNTHFILE_INTERNAL = '"+outfile+"'"
                printf,lun," ID = '"+outfile+"'"
                printf,lun," DATE = '"+systime(0)+"'"
                printf,lun," SPECLIB_VERS = '"+speclib_version()+"'"
                if keyword_set(vmicrofit) and planstr.nvt gt 1 and ipart eq 0 then begin
                  printf,lun," N_OF_DIM = ",string(n_elements(dim)-1)
                  printf,lun," N_P =  ",string(nsteps[1:*])
                  for i=1,n_elements(dim)-1 do $
                    printf,lun,' LABEL('+string(i,format='(i1)')+') = '+$
                                "'"+dim[i]+"'"
                  printf,lun,' LLIMITS =  '+string(llimits[1:*],format='(30f11.5)')
                  printf,lun,' STEPS =  '+string(steps[1:*],format='(30f11.5)')
                endif else begin 
                  printf,lun," N_OF_DIM = ",string(n_elements(dim))
                  printf,lun," N_P =  ",string(nsteps)
                  for i=0,n_elements(dim)-1 do $
                    printf,lun,' LABEL('+string(i+1,format='(i1)')+') = '+$
                                "'"+dim[i]+"'"
                  printf,lun,' LLIMITS =  '+string(llimits,format='(30f11.5)')
                  printf,lun,' STEPS =  '+string(steps,format='(30f11.5)')
                endelse
                if ipart eq 0 then npixout=npix[j] else if j eq 0 then npixout=ipca[ipart-1] else npixout=0
                printf,lun," NPIX = "+string(npixout)
                printf,lun," WAVE =  "+string(w0[j])+' '+string(dw[j])
                printf,lun," LOGW =  1 "
                printf,lun," VACUUM =  1 "
                printf,lun," RESOLUTION =  "+string(sxpar(hdr0,'RESOLUTI'))
                printf,lun," CONTINUUM =  "+string(format='(4f10.4)',$
                 sxpar(hdr0,'ORDER'),sxpar(hdr0,'NITER'),$
                 sxpar(hdr0,'LOWREJ'),sxpar(hdr0,'HIGHREJ'))
                ;printf,lun," ORIGINAL_SAMPLING=  "+ string(sxpar(hdrraw,'CDELT1'))
                if sxpar(hdr,'LINELIST') then $
                  printf,lun," FILE_DATA19= '"+sxpar(hdr,"LINELIST")+"'"
                printf,lun," COMMENTS1 = "
                printf,lun," /"
               endfor
              endfor
              endfor
            endif
            ; write out the spectrum!
            if not keyword_set(nonorm) then  begin
              out=[s0[*,iteff,ilogg,jmh,ielem,irot],$
                   s1[*,iteff,ilogg,jmh,ielem,irot],s2[*,iteff,ilogg,jmh,ielem,irot]]
junk=where(finite(out) eq 0,nbad)
if nbad gt 0 then stop,'non finite pixels',iteff,ilogg,jmh,name
              is=0
              for ipart=ipart1,ipart2 do begin
               if ipart gt 0 then begin
                is=(ipart-1)*nout
                if ipart eq npart then ie=szout[0]-1 else ie=is+nout-1
                printf,nlun[ipart]+npart,out[is:ie],format='(100000(e13.5,1x))'
               endif
              endfor
            endif
            ;printf,nlun,s0[*,iteff,ilogg,jmh,ielem,irot],$
            ;          s1[*,iteff,ilogg,jmh,ielem,irot],$
            ;          s2[*,iteff,ilogg,jmh,ielem,irot],format='(100000(e13.5,1x))'

            ; create input file for vmicro reduction
            if keyword_set(vmicrofit) and planstr.nvt gt 1 then begin
              if ivt eq 0 then begin
               ; only setup for one part!
               for ipart=ipart1,ipart1 do begin
                if nmod eq 0 then openw,vlun,dir+outfile+'_'+string(format='(i2.2)',ipart)+'.inp',/get_lun
                vout=vcoef[0]+vcoef[1]*cm+vcoef[2]*nm+vcoef[3]*am+$
                              vcoef[4]*mh+vcoef[5]*logg+vcoef[6]*teff
                ;vout=2.25-0.3*logg
                cmeps=0. & nmeps=0. & ameps=0. & mheps=0. & loggeps=0. & teffeps=0. & roteps=0.
                if icm eq 0 then cmeps=0.001 else if icm eq ncm-1 then cmeps=-0.001
                if inm eq 0 then nmeps=0.001 else if inm eq nnm-1 then nmeps=-0.001
                if iam eq 0 then ameps=0.001 else if iam eq nam-1 then ameps=-0.001
                if jmh eq 0 then mheps=0.001 else if jmh eq outdim[2]-1 then mheps=-0.001
                if ilogg eq 0 then loggeps=0.001 else if ilogg eq nlogg-1 then loggeps=-0.001
                if iteff eq 0 then teffeps=0.001 else if iteff eq nteff-1 then teffeps=-0.001
                if irot eq 0 then roteps=0.001 else if irot eq nrot-1 then roteps=-0.001
                if nrot gt 1 then vlunout = [alog10(vout), cm+cmeps, nm+nmeps,am+ameps,rot+roteps,mh+mheps, logg+loggeps, teff+teffeps] $
                else vlunout = [alog10(vout),cm+cmeps, nm+nmeps,am+ameps,mh+mheps, logg+loggeps, teff+teffeps]
                printf,vlun,'dummy ',vlunout, format='(a8,'+strtrim(ndim,2)+'f11.5)'
                ;printf,vlun,'dummy ',alog10(vout),mh+mheps,logg+loggeps,teff+teffeps,format='(a8,7f10.3)'
               endfor
              endif
            endif
            nmod+=1
          endfor
         endfor
        endfor
       endfor
      endfor
     endfor
    endfor
  endfor
 endfor
 if keyword_set(vmicrofit) and planstr.nvt gt 1 and ivt eq 0 then free_lun,vlun
endfor
if ipart1 eq 0 then begin
  free_lun,flun
  free_lun,plun
endif
for ipart=ipart1,ipart2 do begin
  close,nlun[ipart]
  close,nlun[ipart]+npart
endfor

if tag_exist(planstr,'holefile') then begin
  mkhdr,hdr,gridhole
  mwrfits,gridhole,outfile+'hole.fits',hdr,/create
endif

; if we are doing a 7D->6D collapse, use ferre to do it!
if keyword_set(vmicrofit) and planstr.nvt gt 1 then begin
  for ipart=ipart1,ipart2 do begin
   print,'running ferre...',ipart
   if ipart ne 0 then begin
    spawn,'cat '+dir+'n_aps'+outfile+'_'+string(ipart,format='(i2.2)')+'.hdr '+$
                 dir+'n_aps'+outfile+'_'+string(ipart,format='(i2.2)')+'.dat >'+$
                 dir+'temp_'+string(ipart,format='(i2.2)')+'.dat'
    ;writeferre
    openw,nml,dir+'input.nml',/get_lun
    printf,nml,'&LISTA'
    ;printf,nml,'NDIM = 4'
    ;printf,nml,'NDIM = 7'
    printf,nml,'NDIM = '+strtrim(ndim,2)
    printf,nml,'NOV = 0'
    ;printf,nml,'INDV = 1 2 3 4'
    ;printf,nml,'INDV = 1 2 3 4 5 6 7'
    printf,nml,'INDV = '+strjoin(strtrim(indgen(ndim)+1,2), ' ')
    printf,nml,"SYNTHFILE(1) = '"+dir+"temp_"+string(ipart,format='(i2.2)')+".dat'"
    printf,nml,"PFILE ='"+dir+outfile+"_"+string(ipart,format='(i2.2)')+".inp'"
    printf,nml,"OFFILE = '"+dir+"n_aps"+outfile+'_'+string(ipart,format='(i2.2)')+".dat'"
    printf,nml,'INTER=3'
    printf,nml,'F_FORMAT=0'
    printf,nml,'F_ACCESS=0'
    printf,nml,'/'
    free_lun,nml
    if dir ne '' then cd,dir,current=cwd
    spawn,'ferre.x'
    if dir ne '' then begin
      cd,cwd
      file_move,dir+outfile+'_'+string(format='(i2.2)',ipart)+'.inp',outfile+'_'+string(format='(i2.2)',ipart)+'.inp',/over
      file_move,dir+'input.nml',outfile+'_'+string(format='(i2.2)',ipart)+'.nml',/over
    endif
    print,'nlines: ', file_lines(dir+'n_aps'+outfile+'_'+string(ipart,format='(i2.2)')+'.dat')
    print,'nlines inp: ', file_lines(outfile+'_'+string(ipart,format='(i2.2)')+'.inp')
   endif
  endfor
endif

dopca:

if npart le 0 then begin
  ; need to determine how to get this to work without parts ...
  print,'pca_synth...'
  ;pca_synth,'n_aps'+outfile+'_w123.dat',24,sections=10 
  pca_synth,dir+'n_aps'+outfile+'_w123.dat',npca,sections=30 ,outsynthfile=dir+ppre+'_aps'+outfile+'_w123.dat'
endif else begin
  for ipart=ipart1,ipart2 do begin
    ; copy n_*.hdr files
    if ipart eq 0 then $
      file='n_aps'+outfile+'.hdr' $
    else $
      file='n_aps'+outfile+'_'+string(ipart,format='(i2.2)')+'.hdr'
    if dir ne '' then if file_test(dir+file) then file_move,dir+file,file,/over

    ; do pca and copy over n_ and p_ .dat
    if ipart eq 0 then $
      file='n_aps'+outfile+'.dat' $
    else $
      file='n_aps'+outfile+'_'+string(ipart,format='(i2.2)')+'.dat'
    if ipart ge 1 then begin
      speclib_pca,dir+file,npca,addnoise=addnoise
      pfile='_aps'+outfile+'_'+string(ipart,format='(i2.2)')+'.dat'
      if dir ne '' then file_move,dir+'p'+pfile,ppre+pfile,/over
    endif
    if dir ne '' and keyword_set(savenorm) then begin
      if keyword_set(vmicrofit) and planstr.nvt gt 1 then file_move,dir+file,'n6_'+strmid(file,2),/over $
      else file_move,dir+file,file,/over
    endif else begin
      if keyword_set(vmicrofit) and planstr.nvt gt 1 then file_delete,dir+file,'n6_'+strmid(file,2),/allow $
      else file_delete,dir+file,file,/allow
    endelse

    ; copy the header files for combined pieces
    if ipart eq 0 and dir ne '' then begin
      file=ppre+'_aps'+outfile+'_w123.dat'
      file_move,dir+file,file,/over
      file='f_'+outfile+'.dat'
      file_move,dir+file,file,/over
    endif
  endfor
  print,' You now need to run the PCA compression on the individual parts, e.g.: '
  print,'   runspeclib_pca.pbs {vers} '+'n_aps'+outfile+'_??.dat'
  print,' When the individual PCA files are done, enter .con to continue'
  ;stop
  ;spawn,'paste '+ppre+'_aps'+outfile+'_??.dat >> '+ppre+'_aps'+outfile+'_w123.dat'

endelse

; now convert to binary
;openw,com,ppre+'_aps'+outfile+'_w123.inp',/get_lun
;printf,com,dir+ppre+'_aps'+outfile+'_w123.dat'
;printf,com,'unf'
;free_lun,com
;spawn,'ascii2bin < '+ppre+'_aps'+outfile+'_w123.inp'


end
