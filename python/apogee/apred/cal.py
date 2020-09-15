# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: October 2018
# @Filename: wave.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

# general routines for APOGEE calibration

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import glob
import pdb
import numpy as np
from astropy.io import fits

def mkpar(mjdstart,mjdend,out='wave',lco=False,yearout='multiwave', append=False,maxsky=41) :
    """ Make calibration file list for wavecals between input dates
    """
    # open output file
    print(mjdstart,mjdend,out)
    if append: f = open(out,'a+')
    else : f = open(out,'w')
    indiv = []
    psf = []
    sky = []
    bad=np.loadtxt(os.environ['APOGEE_DIR']+'/data/cal/bad',dtype=int)
    for mjd in range(mjdstart,mjdend) :
        # get the files for this dat
        if lco :
            files = sorted(glob.glob(os.environ['APOGEE_DATA_2S']+'/'+str(mjd)+'/*-a-*.apz'))
        else :
            files = sorted(glob.glob(os.environ['APOGEE_DATA']+'/'+str(mjd)+'/*-a-*.apz'))
        print(mjd)
        dome=0
        if len(files) > 3:
            # look for sequences of QUARTZ, THARNE, UNE or THARNE, UNE, QUARTZ all at the same dither position
            hdr1 = fits.open(files[0])[1].header
            hdr2 = fits.open(files[1])[1].header
            for ifile in range(2,len(files)-1) :
                try :
                    hdr3 = fits.open(files[ifile])[1].header
                    #print(hdr1['DITHPIX'],hdr2['DITHPIX'],hdr3['DITHPIX'],hdr1['LAMPQRTZ'],hdr2['LAMPTHAR'],hdr3['LAMPUNE'])
                    gd = False
                    if hdr2['DITHPIX'] == hdr3['DITHPIX']  and hdr1['LAMPQRTZ'] and hdr2['LAMPTHAR'] and hdr3['LAMPUNE'] :
                        thar = ifile-1
                        une = ifile
                        qrtz = ifile-2
                        gd = True
                    elif hdr1['DITHPIX'] == hdr2['DITHPIX']  and hdr3['LAMPQRTZ'] and hdr1['LAMPTHAR'] and hdr2['LAMPUNE'] :
                        thar = ifile-2
                        une = ifile-1
                        qrtz = ifile
                        gd = True
                    if gd :
                        # remove if indicated as bad in calibration data file
                        num=int(files[thar].split('-')[2].replace('.apz',''))
                        if len(np.where(bad == num)[0]) > 0 : gd = False
                        num0 = (num // 10000 ) * 10000
                        if len(np.where(bad == num0)[0]) > 0 : gd = False

                    # write out the wave information and accumulate multiwave information
                    if gd :
                        f.write('wave 99999 99999 {:8s} {:8s},{:8s} {:8s}\n'.format(
                             files[thar].split('-')[2].replace('.apz',''),
                             files[thar].split('-')[2].replace('.apz',''),
                             files[une].split('-')[2].replace('.apz',''),
                             files[qrtz].split('-')[2].replace('.apz','')))
                        indiv.append(files[thar].split('-')[2].replace('.apz','') )
                        indiv.append(files[une].split('-')[2].replace('.apz','') )
                        psf.append(files[qrtz].split('-')[2].replace('.apz','') )
                        psf.append(files[qrtz].split('-')[2].replace('.apz','') )
                    # look for sky frames and preceding domeflat, for LSF product (should be following domeflat for APO!)
                    if hdr3['IMAGETYP'] == 'DomeFlat' : dome = ifile
                    if hdr1['IMAGETYP'].strip() == 'Object' and hdr1['NFRAMES']>10 and hdr1['NFRAMES']<maxsky :
                       sky.append(files[ifile-2].split('-')[2].replace('.apz','') )
                       f.write('lsf 99999 99999 {:8s} {:8s} {:8s}\n'.format(
                             files[ifile-2].split('-')[2].replace('.apz',''),
                             files[ifile-2].split('-')[2].replace('.apz',''),
                             files[dome].split('-')[2].replace('.apz','')))
                except :
                    pass
                hdr1=hdr2
                hdr2=hdr3
    # now write the multiwave lines for groups of 20 individual wavecals
    for i in range(0,len(indiv),20) :
        name = indiv[i][0:4]+'0000'
        f.write('multiwave 99999 99999 {:8s} {:8s}'.format(name,indiv[i]))
        for j in range(1,20) : 
            try: f.write(',{:8s}'.format(indiv[i+j]))
            except : pass
        f.write('\n')
    f.close()

    # now write the multiwave lines for the full period spaced by 10 days
    f = open(yearout,'a+')
    last = 0
    for i in range(0,len(indiv),2) :
        name = str((mjdstart-55562)*10000)
        if i == 0 : 
            f.write('multiwave {:d} {:d} {:8s} {:8s},{:8s}'.format(mjdstart, mjdend, name, indiv[i],indiv[i+1] ))
        elif int(indiv[i][0:4]) > last+10 : 
            f.write(',{:8s},{:8s}'.format(indiv[i],indiv[i+1]))
            last = int(indiv[i+1][0:4])
    f.write('\n')
    f.close()

def mkallpars(apo=True,lco=True) :
    if lco :
        try : os.remove('apogee-s-wave.par')
        except : pass
        mjds= [57829, 57966, 58360, 58700, 58950]
        mjds= [58360, 58700, 58950]
        for i in range(len(mjds)-1) :
            mkpar(mjds[i],mjds[i+1],out='apogee-s-year{:d}.par'.format(i),lco=True,yearout='apogee-s-multiwave.par', append=True,maxsky=41)
    if apo :
        try : os.remove('apogee-n-wave.par')
        except : pass
        mjds= [55800, 56130, 56512, 56876, 57230, 57600, 57966, 58360, 58700, 58950]
        mjds= [58360, 58700, 58950]
        for i in range(len(mjds)-1) :
            mkpar(mjds[i],mjds[i+1],out='apogee-n-year{:d}.par'.format(i),yearout='apogee-n-multiwave.par', append=True,maxsky=20)
        
