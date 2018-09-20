# Edit and move fits
# A module to move and edit the raw FITS files as the become available from the ICS
import numpy as np
import pyfits 
import os
import datetime
import re
def movefits(CurrentFileName,DayNumber,current_obs_ra,current_obs_dec,current_obs_az,current_obs_alt,current_epoch):
    
        NEWCARDS  = ['RA','DEC','EPOCH','AZ','ALT','TELESCOP']
        NEWValues = [current_obs_ra,current_obs_dec,current_epoch,current_obs_az,current_obs_alt,'NMSU 1m']
        timestamp = datetime.datetime.now()
        rawdirec  = '/data-ics/'+DayNumber # Directory for Raw FITS files
        pathraw   = os.path.abspath(rawdirec) # Path to raw files
        MJD5      = int(DayNumber)+55562 # MJD calculted with daynumber
        direc = '/data-ql/data/'+str(MJD5) # Directory where edited FITS will be saved
	if os.path.exists(direc)!=1:
		os.mkdir(direc)
        editdirec = '/data-ql/data/'+str(MJD5)+'/1m/' # Directory where edited FITS will be saved
	t = editdirec
	if os.path.exists(t)!=1:
		os.mkdir(t)
        pathedit  = os.path.abspath(editdirec) # Path for edited FITS
        time      = str(datetime.datetime.now())


        # first extract the value of the checksum from the fits header (pyfits.getval removes
        # any checksum or datasum keywords)
        f=open(pathraw+'/'+CurrentFileName,'rb')
        checksum = None
        # only read the first 72 lines (which should be the whole header plus padding)
        for p in range(72):
          line = f.read(80)
          if line[0:8] == 'END     ':
              break
          if line[0:8] == 'CHECKSUM':
              checksum = line[11:27]
              cs_comment = line[33:80]
        f.close()

        # open the image with pyfits
        #img = pyfits.open(pathraw+'/'+CurrentFileName,do_not_scale_image_data=True)
        img = pyfits.open(pathraw+'/'+CurrentFileName,do_not_scale_image_data=True,uint16=True)
        print 'checksum: ' +checksum
        if checksum != None:
          # validate the value of the checksum found (corresponding to DATASUM in pyfits)
          # calulate the datasum
          ds = img[0]._calculate_datasum('standard')

          # add a new CHECKSUM line to the header (pyfits.open removes it) with same comment
          print 'updating header CHECKSUM ' + '0'*16 + cs_comment
          img[0].header.update("CHECKSUM",'0'*16, cs_comment)

          # calulate a new checksum
          cs=img[0]._calculate_checksum(ds,'standard')
          img[0].header.update("CHECKSUM",cs, cs_comment)
          print 'checksum ', checksum
          print 'ds ', ds
          print 'cs ', cs
          if cs != checksum:
              print "CHECKSUM Failed for file " + CurrentFileName

        # force these to be ints:
        # As of August 2013, the ICS writes them both as floats, but the
        # FITS standard wants them to be ints.
        bscale = int(img[0].header.get('BSCALE',1))
        bzero = int(img[0].header.get('BZERO',32768))
        del img[0].header['BSCALE']
        del img[0].header['BZERO']
        img[0].header.update('BSCALE',bscale,after='GCOUNT')
        img[0].header.update('BZERO',bzero,after='BSCALE')

        strp = re.sub('.fits',"",CurrentFileName) # strip .fits of file name
        new  = strp + '.fits' # add edit.fits to file name
        print 'before'
        print img[0].header
        for i in range(len(NEWCARDS)):
		img[0].header.update(NEWCARDS[i],NEWValues[i],'Taken from 1-meter')
        print 'after'
        print img[0].header
	img[0].header.add_history('FITS file edited'+' '+time)
	img.writeto(pathedit+'/'+new, checksum=True)
        print 'Done editing',CurrentFileName
	
	return
