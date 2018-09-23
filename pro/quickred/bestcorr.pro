function bestcorr,files,frameid

  ; given list of calibration filenames that have framenumbers encoded in them
  ; choose the one closest to the input frameid, and return that filename
  i=strpos(file_basename(files[0]),'.fits')
  frameids = strmid(file_basename(files,'.fits'),i-8,8)
  framenums = long(frameids)
  best = first_el(minloc( abs(framenums-long(frameid)) ))
  corr =frameids[best]
  return,corr
end
