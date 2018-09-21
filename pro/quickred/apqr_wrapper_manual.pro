;+
;
; APQR_WRAPPER_MANUAL
;
; This is the wrapper to the APOGEE quickreduce program without the actor
;
; INPUTS:
;   inlist: a csv file containing a list of exposures to be processed.  The format is: mjd, exp, exp_pk, plugfile
;   no_dbinsert: if this is not set, then the Apogee Quicklook DB will be updated  
;   data_dir, spectro_dir, archive_dir and quickred_dir can be set manually, otherwise they are set from the system values
; OUTPUTS:
;   A DB insert if no_dbinsert != 1.
; USAGE:
;   e.g.,
;   apqr_wrapper_manual,'inlist.txt', /no_dbinsert
;-
;

pro apqr_wrapper_manual, inlist, no_dbinsert=no_dbinsert, outfile=outfile, data_dir=data_dir, spectro_dir=spectro_dir, $
  archive_dir=archive_dir, quickred_dir=quickred_dir

   ;read in data from inlist
   if FILE_TEST(inlist) eq 0 then begin
      msg = "Input file not found: "+inlist
      print,msg
      return
   endif

   data = read_csv(inlist)
   mjds = strtrim(data.field1,2)
   frameids = strtrim(data.field2,2)
   exp_pks = strtrim(data.field3,2)
   plugfiles = strtrim(data.field4,2)
   nframes = n_elements(frameids)

   ; Get APOGEE directories if the environment variables exists
   if not keyword_set(data_dir) then data_dir = APGETDIR('APQLDATA_DIR',/exists,error=direrr)
   if not keyword_set(spectro_dir) then spectro_dir = APGETDIR('APQLSPECTRO_DIR',/exists,error=direrr)
   if not keyword_set(archive_dir) then archive_dir = APGETDIR('APQLARCHIVE_DIR',/exists,error=direrr)
   if not keyword_set(quickred_dir) then quickred_dir = APGETDIR('APQLQUICKRED_DIR',/exists,error=direrr)

   ; we need data_dir and spectro_dir
   if strlen(data_dir) eq 0 then begin
      msg = 'Missing data_dir -> aborted'
      print,msg
      return
   endif
   if strlen(spectro_dir) eq 0 then begin
      msg = 'Missing spectro_dir -> aborted'
      print,msg
      return
   endif
   if strlen(archive_dir) eq 0 then begin
      msg = 'Missing archive_dir -> aborted'
      print,msg
      return
   endif
      if strlen(quickred_dir) eq 0 then begin
      msg = 'Missing quickred_dir -> aborted'
      print,msg
      return
   endif

   psfdir = spectro_dir+'/cal/psf/'
   bpmdir = spectro_dir+'/cal/bpm/'

   ; Get Psf file
   chiptag = ['a','b','c']
   tinfo = APFILEINFO(FILE_SEARCH(psfdir+'apPSF-a-*.fits'),/silent)
   gdpsf = where(tinfo.exists eq 1 and tinfo.allchips eq 1 and tinfo.suffix ne '',ngdpsf)
   if ngdpsf eq 0 then begin
      msg='No good psf file found in '+psfdir
      print,msg
      return
   endif else begin
      ; sort through the files and use the latest one
      tinfogd = tinfo[gdpsf]       ; the good ones
      si = reverse(sort(tinfogd.mtime))
      psfid = psfdir+tinfogd[si[0]].suffix
   endelse

   ; Get BPM file
   binfo = APFILEINFO(FILE_SEARCH(bpmdir+'apBPM-a-*.fits'),/silent)
   gdbpm = where(binfo.exists eq 1 and binfo.allchips eq 1 and binfo.suffix ne '',ngdbpm)
   if ngdbpm eq 0 then begin
      msg='No good BPM file found in '+bpmdir
      print,msg
      ;return
   endif else begin
      ; sort through the files and use the latest one
      binfogd = binfo[gdbpm]       ; the good ones
      si = reverse(sort(binfogd.mtime))
      bpmid = bpmdir+binfogd[si[0]].suffix
   endelse

   ;for each exposure/frame
   for i=0,nframes-1 do begin
      exp_pk = exp_pks[i]
      mjd5 = mjds[i]
      frameid = frameids[i]
      plugfile = plugfiles[i]

      rawdir = data_dir+strtrim(mjd5,2)+'/'
      bundledir = archive_dir+strtrim(mjd5,2)+'/'
      quickreddir = quickred_dir+strtrim(mjd5,2)+'/'

      ;run quickred
      apquickred,frameid,plugfile=plugfile,rawdir=rawdir,bundledir=bundledir, $
          quickreddir=quickreddir,bpmid=bpmid,psfid=psfid,exp_pk=exp_pk,mjd5=mjd5, $
          no_dbinsert=no_dbinsert,outfile=outfile;,/no_compress

    endfor
end
