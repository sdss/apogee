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
from tools import html
from tools import plots
from apogee.utils import apload

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
        

def darkplot(apred='r14',telescope='apo25m'):
    """ Make some plots of dark frames for a reduction version
    """
    chips=['a','b','c'] 
    if telescope == 'apo25m' :
        prefix='ap'
        xbin=np.arange(1,1200,0.2)
        inst='apogee-n'
    else :
        prefix='as'
        xbin=np.arange(1,300,0.2)
        inst='apogee-s'

    # get dark frames
    load=apload.ApLoad(apred=apred,telescope=telescope)
    darkdir=os.path.dirname(load.filename('Dark',chips='a',num=0))
    try: os.makedirs(darkdir)
    except: pass
    files=np.sort(glob.glob(darkdir+'/'+prefix+'Dark-a-*.fits'))

    ny=len(files)//2
    fig,ax=plots.multi(2,3,hspace=0.001,wspace=0.001,figsize=(14,ny*4))
    imfig=[]
    imax=[]
    if ny%2 == 1 : ny+=1
    for ichip in range(3) :
        tfig,tax=plots.multi(2,ny,hspace=0.001,wspace=0.001,figsize=(10,8))
        imfig.append(tfig)
        imax.append(tax)
    i=0
    for file in files :
        im = os.path.basename(file).replace('.fits','').split('-')[-1]
        print(im)
        for ichip,chip in enumerate(chips) :
            dark=fits.open('{:s}/{:s}DarkRate-{:s}-{:s}.fits'.format(
                           darkdir,prefix,chip,im))[0].data
            ax[ichip,0].hist(dark.flatten(),bins=np.arange(0,1,0.02),label='{:s}'.format(im),histtype='step')
            ax[ichip,1].hist(dark.flatten(),bins=xbin,label='{:s}'.format(im),histtype='step')
            imax[ichip][i%ny,i//ny].imshow(dark,vmin=0,vmax=0.25,cmap='Greys')
            imax[ichip][i%ny,i//ny].text(0.1,0.9,im,transform=imax[ichip][i%ny,i//ny].transAxes)
            imax[ichip][i%ny,i//ny].get_xaxis().set_visible(False)
            imax[ichip][i%ny,i//ny].get_yaxis().set_visible(False)
        i+=1
    while i<2*ny :
        for ichip in range(3) :
            imax[ichip][i%ny,i//ny].set_visible(False)
            imax[ichip][i%ny,i//ny].set_visible(False)
        i+=1

    for ichip,chip in enumerate(chips) :
        ax[ichip,0].set_title('chip: {:s}'.format(chip))
        ax[ichip,1].set_title('chip: {:s}'.format(chip))
        ax[ichip,0].legend(fontsize='x-small')
        ax[ichip,0].set_yscale('log')
        ax[ichip,1].legend(fontsize='x-small')
        ax[ichip,1].set_yscale('log')
        #imfig[ichip].suptitle('chip: {:s}'.format(chip))
        imfig[ichip].savefig('{:s}/plots/{:s}_{:s}_darkimage.png'.format(darkdir,inst,chip))
    ax[2,0].set_xlabel('Dark rate (cnts/read)')
    ax[2,1].set_xlabel('Dark rate (cnts/read)')
    fig.tight_layout()
    fig.savefig('{:s}/plots/{:s}_darkhist.png'.format(darkdir,inst))

    grid=[['{:s}_darkhist.png'.format(inst)]]
    for ichip,chip in enumerate(chips) :
        grid.append(['{:s}_{:s}_darkimage.png'.format(inst,chip)])

    html.htmltab(grid,file='{:s}/plots/{:s}_darks.html'.format(darkdir,inst))

def flatplot(apred='r14',telescope='apo25m'):
    """ Make flat plot of flat frames for a reduction version
    """
    chips=['a','b','c'] 
    if telescope == 'apo25m' :
        prefix='ap'
        xbin=np.arange(1,1200,0.2)
        inst='apogee-n'
    else :
        prefix='as'
        xbin=np.arange(1,300,0.2)
        inst='apogee-s'

    # get flat frames
    load=apload.ApLoad(apred=apred,telescope=telescope)
    flatdir=os.path.dirname(load.filename('Flat',chips='a',num=0))
    try: os.makedirs(flatdir)
    except: pass
    files=np.sort(glob.glob(flatdir+'/'+prefix+'Flat-a-*.fits'))

    ny=len(files)//2
    if ny%2 == 1 : ny+=1
    fig,ax=plots.multi(1,3,figsize=(8,ny*4))
    imfig=[]
    imax=[]
    for ichip in range(3) :
        tfig,tax=plots.multi(2,ny,hspace=0.001,wspace=0.001,figsize=(10,8))
        imfig.append(tfig)
        imax.append(tax)
    i=0
    for file in files :
        im = os.path.basename(file).replace('.fits','').split('-')[-1]
        print(im)
        for ichip,chip in enumerate(chips) :
            flat=fits.open('{:s}/{:s}Flat-{:s}-{:s}.fits'.format(
                           flatdir,prefix,chip,im))[1].data
            ax[ichip].hist(flat.flatten(),bins=np.arange(0,2,0.02),label='{:s}'.format(im),histtype='step')
            med=np.median(flat)
            imax[ichip][i%ny,i//ny].imshow(flat,vmin=0.8*med,vmax=1.2*med,cmap='Greys')
            imax[ichip][i%ny,i//ny].text(0.1,0.9,im,transform=imax[ichip][i%ny,i//ny].transAxes)
            imax[ichip][i%ny,i//ny].get_xaxis().set_visible(False)
            imax[ichip][i%ny,i//ny].get_yaxis().set_visible(False)
        i+=1
    while i<2*ny :
        for ichip in range(3) :
            imax[ichip][i%ny,i//ny].set_visible(False)
            imax[ichip][i%ny,i//ny].set_visible(False)
        i+=1

    for ichip,chip in enumerate(chips) :
        ax[ichip].set_title('chip: {:s}'.format(chip))
        ax[ichip].legend(fontsize='x-small')
        ax[ichip].set_yscale('log')
        #imfig[ichip].suptitle('chip: {:s}'.format(chip))
        imfig[ichip].savefig('{:s}/plots/{:s}_{:s}_flatimage.png'.format(flatdir,inst,chip))
    fig.tight_layout()
    out='{:s}/plots/{:s}_flathist.png'.format(flatdir,inst)
    fig.savefig(out)
    grid=[[os.path.basename(out)]]
    for ichip,chip in enumerate(chips) :
        grid.append(['{:s}_{:s}_flatimage.png'.format(inst,chip)])

    html.htmltab(grid,file='{:s}/plots/{:s}_flats.html'.format(flatdir,inst))

