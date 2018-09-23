import sys
import os
import glob
import pdb
import numpy as np
from astropy.io import fits

def mjdcube(mjd, darkid=None, write=False, apred='current', clobber=False) :

  """
  Make a cube for a given night with the CDS images of all frames
  Optionally, write out individual uncompressed data cubes
  """

  print('mjd: ', mjd)
  print('apred: ', apred)
  print('write: ', write)
  print('darkid: ', darkid)
  datadir=os.getenv('APOGEE_DATA')
  outdir=os.getenv('APOGEE_REDUX')+'/'+apred+'/exposures/apogee-n/'+str(mjd)+'/'

  for chip in ['a','b','c'] :
    files = sorted(glob.glob(datadir+'/'+str(mjd)+'/apR-'+chip+'-*.apz')+glob.glob(datadir+'1m/'+str(mjd)+'/apR-'+chip+'-*.apz'))

    # output file name for CDS cube
    outfile = outdir+'apHist-'+chip+'-'+str(mjd)+'.fits'

    # does output file already exist?
    if not clobber and os.path.exists(outfile) : return

    out=fits.HDUList(fits.PrimaryHDU())

    # get dark frame if requested
    if darkid is not None :
        dark=fits.open(os.environ['APOGEE_REDUX']+'/'+apred+'/cal/darkcorr/apDark-'+chip+'-'+darkid+'.fits')[1].data

    # loop over all files
    for file in files :
      print file
      print 'file: ',file
      if write :
        # output file name for individual uncompressed images
        outfile = os.path.basename(file.strip('apz')+'fits')
        hdout=fits.HDUList(fits.PrimaryHDU())

      # open file and confirm checksums
      hd=fits.open(file, do_not_scale_image_data = True, uint = True, checksum = True)

      # file has initial header, avg_dcounts, then nreads
      nreads = len(hd)-2
      try:
          avg_dcounts=hd[1].data
      except:
          # fix header if there is a problem (e.g., MJD=55728, 01660046)
          hd[1].verify('fix')
          avg_dcounts=hd[1].data

      # first read is in extension 2
      ext = 2
  
      # loop over reads, processing into raw reads, and appending
      for read in range(1,nreads+1) :
        header = hd[ext].header
        try:
          raw = hd[ext].data
        except:
          hd[ext].verify('fix')
          raw = hd[ext].data
        if read == 1 :
          data = np.copy(raw)
        else :
          data = np.add(data,raw,dtype=np.int16)
          data = np.add(data,avg_dcounts,dtype=np.int16)
          if read == 2 : first = data
  
        if write : hdout.append(fits.ImageHDU(data,header))
          
        ext += 1
 
      # compute and add the cdsframe, subtract dark if we have one
      cds = (data[0:2048,0:2048] - first[0:2048,0:2048] ).astype(float)
      print cds.shape
      if darkid is not None :
          print dark.shape,nreads
          # if we don't have enough reads in the dark, do nothing
          try :
              cds -= (dark[nreads-1,:,:] - dark[2,:,:])
          except:
              print 'not halting: not enough reads in dark, skipping dark subtraction for mjdcube'
              pass
      out.append(fits.ImageHDU(cds,header))
      if write :
        hdout.writeto(outfile,clobber=True, checksum = True, output_verify='fix')
        hd.close()

    # write out the CDS frame
    out.writeto(outfile,clobber=True, checksum = True, output_verify='fix')
    out.close()

if __name__ == "__main__" :
    mjdcube(sys.argv[1],sys.argv[2:]) 
