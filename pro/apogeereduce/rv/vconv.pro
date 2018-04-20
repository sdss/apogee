PRO VCONV,vhel,glon,glat,vlsr,vgsr,vcirc=vcirc

;+
;
; VCONV
;
; Given an input array called vhel with values of heliocentric
; radial velocities and input arrays of the same length with
; the Galactic coordinates of the object (gl,gb),
; this code calculates Vlsr and Vgsr.
;
;  code assumes gl & gb given in degrees.
;
; INPUTS:
;  vhel   Array of heliocentric velocities (in km/s).
;  glon   Galactic longitude (in degrees)
;  glat   Galactic latitude (in degrees)
;  =vcirc Circular rotation velocity at the sun (220 km/s by default)
;
; OUTPUTS:
;  vlsr   Array of Local Standard of Rest velocities
;  vgsr   Array of Galactic Standard of Rest velocities
;
;
; USAGE:
;  IDL>vconv,vhelio,glon,glat,vlsr,vgsr
;
; By ??
;  updated by D.NIdever
;-

; Are there enough inputs
if n_params() lt 3 then begin
  print,'Syntax - vconv,vhel,glon,glat,vlsr,vgsr,vcirc=vcirc'
  return
endif

if not keyword_set(vcirc) then vcirc=220.0d

gl=glon*(!dpi/180.0d0)
gb=glat*(!dpi/180.0d0)

cgl=cos(gl) & sgl=sin(gl)
cgb=cos(gb) & sgb=sin(gb)

;  This equation takes the solar motion w.r.t. the LSR as
;  (9,11,6) km/sec (Mihalas & Binney)
;   Using updated values from Dehnen & Binney 1998, MNRAS, 298, 387
;   U = 10.00 +/- 0.36 km/s
;   V = 5.25 +/- 0.62 km/s
;   W = 7.17 +/- 0.38 km/s
;   This is in a right-handed system

vlsr=vhel+((10.00d0*cgb*cgl)+(5.25d0*cgb*sgl)+(7.17d0*sgb))
;vlsr=vhel+((9.0d0*cgb*cgl)+(11.0d0*cgb*sgl)+(6.0d0*sgb))

;  This equation takes the rotation velocity of the LSR to
;  be 220 km/sec

vgsr=vlsr+(vcirc*cgb*sgl)

return

end

