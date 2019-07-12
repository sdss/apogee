pro apditherpairs,shiftstr,pairstr,error=error,verbose=verbose,stp=stp,snsort=snsort

;+
;
; APDITHERPAIRS
;
; This program pairs up a number of dithered APOGEE frames
;
; INPUTS:
;  shiftstr     A structure that gives the shifts for all of the frames
;                 relative to the first. A positive shift means that
;                 frame2 is to the RIGHT of frame1 (the reference frame).
;                 This is normally measured by APDITHERSHIFT.PRO
;                 The SHIFTSTR structure should have the following tag:
;                  INDEX - A running index
;                  FRAMENUM - The ID8 frame number (string)
;                  SHIFT - The shift relative to the first frame
;                  SHIFTERR - Error in SHIFT
;  /verbose     Print a lot to the screen
;  /stp         Stop at the end of the program
;
; OUTPUTS:
;  pairstr      A structure that pairs up the dithered frames,
;                 gives their relative shifts and other information
;                 that is needed to for APDITHERCOMB.PRO do combine
;                 the frames.  The 1st frame of the 1st pair is the
;                 new reference frame and has the most negative shift
;                 (original).
;                 The PAIRSTR structure has tags:
;                  FRAMENAME - 2-element string array of frame numbers
;                  FRAMENUM - Same as FRAMENAME but long type
;                  OLDSHIFT - 2-element array of the original SHIFT
;                  SHIFT - 2-element array of the shifts relative to
;                           the NEW reference image (1st frame of 1st pair)
;                  RELSHIFT - The shift of 2nd frame wrt 1st frame
;                  NUSED - 2-element array giving the number of times
;                            that the frame is used in the PAIRSTR
;                  INDEX - 2-element array giving the index in the
;                            original SHIFTSTR.
;  =error       The error message if one occurred.
;
; USAGE:
;  IDL>apditherpairs,shiftstr,pairstr
;
; By D.Nidever  May 2010
;    J.Holtzman, various mods
;-

; TESTING
;shiftstr = replicate({index:0L,framenum:'',shift:0.0,shifterr:0.0},6)
;shiftstr.index = indgen(6)
;shiftstr.framenum = string(lindgen(6)+101,format='(I08)')
;;shiftstr.shift = [0.0, 0.48,0.01,0.52,0.51,-0.01]
;shiftstr.shift = [0.0, 0.01,0.49,-0.02,0.51,-0.51]
;shiftstr.shifterr = 0.001
;
;shiftstr = replicate({index:0L,framenum:'',shift:0.0,shifterr:0.0},5)
;shiftstr.index = indgen(5)
;shiftstr.framenum = string(lindgen(5)+101,format='(I08)')
;shiftstr.shift = [0.0, 0.01,0.49,-0.02,0.51]
;shiftstr.shifterr = 0.001

apgundef,pairstr

; Not enough inputs
nshiftstr = n_elements(shiftstr)
if nshiftstr eq 0 then begin
  print,'Syntax - apditherpairs,shiftstr,pairstr,error=error,verbose=verbose,stp=stp'
  error = 'Not enough inputs'
  return
endif

nframes = nshiftstr

; Not enough frames
if nframes lt 2 then begin
  error = 'Only ONE frame input. Need at least TWO'
  print,error
  return
endif

; This is how the frames are paired:
; 1.) need a frame that has a relative dither shift of at least 0.2 pixels
; 2.) Closer in time is preferred
; 3.) Not already taken as a pair is also preferred.


; Frame structure for internal use to keep track of things
framestr = REPLICATE({index:0L,framename:'',framenum:0L,shift:0.0,sn:0.0,$
                      pairframename:'',pairframenum:0L,$
                      pairshift:0.0,pairsn:0.0,pairindex:-1L,nused:0L},nframes)
framestr.index = shiftstr.index
framestr.framename = shiftstr.framenum
framestr.framenum = long(shiftstr.framenum)
framestr.shift = shiftstr.shift
framestr.sn = shiftstr.sn

; Print out information
if keyword_set(verbose) then begin
  print,strtrim(nframes,2),' frames input'
  print,' NUM   FRAME      SHIFT    S/N'
  for i=0,nframes-1 do print,format='(I4,A10,2F10.5)',i+1,framestr[i].framename,framestr[i].shift,framestr[i].sn
endif

minshift = 0.2  ; minimum shift for a dither pair
maxshift = 0.8  ; maximum shift for a dither pair
minsn = 3  ; minimum required S/N

; Loop through the frames
if keyword_set(snsort) then $
isort=reverse(sort(framestr.sn)) $
else $
isort=indgen(nframes)
verbose=1
For j=0,nframes-1 do begin
  i=isort[j]
  ; Not paired up yet
  if framestr[i].nused eq 0 then begin

    relshift = framestr[isort].shift-framestr[i].shift ; shift relative to this one
;    relshift[isort[i]] = 999999.0  ; don't want to pair with self

    fracrelshift = relshift-fix(relshift)  ; only the fraction of the shift

    gdframes = where(abs(fracrelshift) ge minshift,ngdframes)

    ; No observation to pair this one with
    if ngdframes eq 0 then begin
      print,'No frame shifted enough for ',framestr[i].framename
      goto,BOMB
    endif

    if keyword_set(snsort) then  begin
      passind = where(abs(fracrelshift) ge minshift AND $
                    abs(fracrelshift) le maxshift AND $
                    framestr[isort].sn gt minsn AND $
                    framestr[isort].framenum ne framestr[i].framenum AND $
                    framestr[isort].nused eq 0,npassind)
    endif else begin
      ; Try following shifted frames that aren't paired yet
      passind = where(abs(fracrelshift) ge minshift AND $
                    framestr[isort].framenum gt framestr[i].framenum AND $
                    framestr[isort].nused eq 0,npassind)

      ; Try preceding shifted frames that aren't paired yet
      if npassind eq 0 then $
      passind = where(abs(fracrelshift) ge minshift AND $
                    framestr[isort].framenum lt framestr[i].framenum AND $
                    framestr[isort].nused eq 0,npassind)
    endelse


    ; Try following shifted frames that ARE paired
;    if npassind eq 0 then $
;    passind = where(abs(fracrelshift) ge minshift AND $
;                    framestr.framenum gt framestr[i].framenum AND $
;                    framestr.nused gt 0,npassind)
;
;    ; Try preceding shifted frames that ARE paired
;    if npassind eq 0 then $
;    passind = where(abs(fracrelshift) ge minshift AND $
;                    framestr.framenum lt framestr[i].framenum AND $
;                    framestr.nused gt 0,npassind)

    ; No observation to pair with
    if npassind eq 0 then begin
      print,'No frame to pair with ',framestr[i].framename
      goto,BOMB
    endif

    ; Do the pairing
    ipair=isort[passind[0]]
    framestr[i].pairframename = framestr[ipair].framename
    framestr[i].pairframenum = framestr[ipair].framenum
    framestr[i].pairshift = framestr[ipair].shift
    framestr[i].pairsn = framestr[ipair].sn
    framestr[i].pairindex = ipair
    framestr[i].nused++
    framestr[ipair].nused++
    if keyword_set(verbose) then $
      print,'Pairing ',framestr[i].framename,' with ',framestr[i].pairframename

    ; Add to the pair structure
    newpairstr = {framename:strarr(2),framenum:lonarr(2),$
                  oldshift:fltarr(2),shift:fltarr(2),sn:fltarr(2),$
                  refshift:0.0,relshift:0.0,nused:lonarr(2),index:lonarr(2)}
    newpairstr.framename = [framestr[i].framename, framestr[i].pairframename]
    newpairstr.framenum = [framestr[i].framenum, framestr[i].pairframenum]
    newpairstr.oldshift = [framestr[i].shift,framestr[i].pairshift]
    newpairstr.shift = [framestr[i].shift,framestr[i].pairshift]
    newpairstr.sn = [framestr[i].sn,framestr[i].pairsn]
    newpairstr.relshift = newpairstr.shift[1]-newpairstr.shift[0] ; relative to first frame
    newpairstr.index = [i,framestr[i].pairindex]
    ; we'll update the Nused values later
    PUSH,pairstr,newpairstr

  endif  ; not paired up yet

  BOMB:

Endfor

npairs = n_elements(pairstr)

if npairs eq 0 then begin
  error = 'NO PAIRS'
  print,error
  return
end

; Put the pair with the most POSITIVE shift first
;  it will have the lowest wavelength on the left
;pairshifts = MIN(pairstr.shift,dim=1)
;refind = first_el(minloc(pairshifts))
pairshifts = MAX(pairstr.shift,dim=1)
refind = first_el(maxloc(pairshifts))
if refind ne 0 then begin
  refpair = pairstr[refind]
  REMOVE,refind,pairstr
  pairstr = [refpair,pairstr]
endif
; Change SHIFT so it is relative to the reference image
;refshift = pairstr[0].shift[0]
;pairstr.shift -= refshift
refshift = max(pairstr[0].shift)
pairstr.shift = refshift - pairstr.shift
; THIS SHIFT now indicates where this frame BEGINS (e.g. wavelength)
;  relative to the reference frame.

if keyword_set(verbose) then begin
  print,strtrim(npairs,2),' PAIRS'
  print,' NUM   FRAME1    FRAME2    SHIFT1    SHIFT2 NUSED1 NUSED2 RELSHIFT'
endif

; Loop through the pair structure
For i=0,npairs-1 do begin

  ; Update NUSED in the pair index
  frame1 = pairstr[i].framename[0]
  frame2 = pairstr[i].framename[1]

  ; How often was the first frame used
  used1 = where(pairstr.framename[0] eq frame1 or $
                pairstr.framename[1] eq frame1,nused1)

  ; How often was the second frame used
  used2 = where(pairstr.framename[0] eq frame2 or $
                pairstr.framename[1] eq frame2,nused2)

  pairstr[i].nused = [nused1,nused2]

  ;; They are in the wrong order, FLIP
  ;if pairstr[i].relshift lt 0.0 then begin
  ;  thispair = pairstr[i]
  ;  pairstr[i].framename = reverse(thispair.framename)
  ;  pairstr[i].framenum = reverse(thispair.framenum)
  ;  pairstr[i].shift = reverse(thispair.shift)
  ;  pairstr[i].nused = reverse(thispair.nused)
  ;  pairstr[i].index = reverse(thispair.index)
  ;  pairstr[i].relshift = -thispair.relshift
  ;endif

  ; They are in the wrong order, FLIP
  ;  want the first spectrum to have the lowest wavelength
  if pairstr[i].relshift gt 0.0 then begin
    thispair = pairstr[i]
    pairstr[i].framename = reverse(thispair.framename)
    pairstr[i].framenum = reverse(thispair.framenum)
    pairstr[i].oldshift = reverse(thispair.oldshift)
    pairstr[i].shift = reverse(thispair.shift)
    pairstr[i].sn = reverse(thispair.sn)
    pairstr[i].nused = reverse(thispair.nused)
    pairstr[i].index = reverse(thispair.index)
    ;pairstr[i].relshift = -thispair.relshift
  ; right order
  endif else begin
    ; want RELSHIFT to be POSITIVE
    pairstr[i].relshift = -pairstr[i].relshift
  endelse

  ; REFSHIFT
  ;  This is the shift of the dither pair relative to the
  ;  "reference" frame.  Should just be equal to the
  ;  "shift" value of the first frame of the pair.
  pairstr[i].refshift = pairstr[i].shift[0]

  ; Printing out the pairs
  if keyword_set(verbose) then begin
    fmt = '(I4,2A10,2F10.5,2I5,F10.5)'
    print,format=fmt,i+1,pairstr[i].framename[0],pairstr[i].framename[1],$
          pairstr[i].shift[0],pairstr[i].shift[1],pairstr[i].nused[0],$
          pairstr[i].nused[1],pairstr[i].relshift
  endif

End

; The "0th" frame might not be first.

; Use the indices of the shiftstr to correctly identify the
; right frame in the "allframes" structure


if keyword_set(stp) then stop

end

