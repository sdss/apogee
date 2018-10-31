# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: October 2018
# @Filename: wave.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

# routines for APOGEE wavelength calibration

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import copy
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import pdb
from functools import wraps
from astropy.io import ascii
from astropy.io import fits
from scipy import signal
from scipy.optimize import curve_fit
from scipy.special import erf, erfc
from scipy.signal import medfilt
from scipy import interpolate
from apogee.utils import apload
from tools import plots
from sdss import yanny

colors=['r','g','b']
xlim=[[16400,17000],[15900,16500],[15100,15800]]

def gauss(x,a,x0,sig) :
    """ Evaluate Gaussian function """
    return a/np.sqrt(2*np.pi)/sig*np.exp(-(x-x0)**2/2./sig**2)

def gaussbin(x,a,x0,sig) :
    """ Evaluate integrated Gaussian function """
    # bin width
    xbin=1.
    t1=(x-x0-xbin/2.)/np.sqrt(2.)/sig
    t2=(x-x0+xbin/2.)/np.sqrt(2.)/sig
    y=(myerf(t2)-myerf(t1))/xbin
    return a*y

def peakfit(spec,pix0,estsig=5,sigma=None,mask=None,plot=False,func=gaussbin) :
    """ Return integrated-Gaussian centers near input pixel center
    """
    x = np.arange(len(spec))
    cen = int(round(pix0))
    sig=estsig
    back=0.
    for iter in range(11) :
        xwid=int(round(5*sig))
        if xwid < 3 : xwid=3
        y=spec[cen-xwid:cen+xwid+1]
        yerr=sigma[cen-xwid:cen+xwid+1]
        x0 = y.argmax()+(cen-xwid)
        peak = y.max()
        sig = np.sqrt(y.sum()**2/peak**2/(2*np.pi))
        pars=curve_fit(func,x[cen-xwid:cen+xwid+1],y,p0=[peak/sig/np.sqrt(2*np.pi),x0,sig],sigma=yerr)[0]
        # iterate unless new array range is the same
        if int(round(5*pars[2])) == xwid and int(round(pars[1])) == cen : break
        cen=int(round(pars[1]))
        sig=pars[2]
    if plot :
        plt.clf()
        plt.plot(x,spec)
        plt.plot(x,func(x,pars[0],pars[1],pars[2]))
        plt.xlim((pars[1]-50,pars[1]+50))
        plt.draw()
        pdb.set_trace()
    return(pars)

def myerf(t) :
    """ Evaluate function that integrates Gaussian from -inf to t"""
    neg = np.where(t<0.)[0]
    pos = np.where(t>=0.)[0]
    out = t*0.
    out[neg] = erfc(abs(t[neg]))/2.
    out[pos] = 0.5+erf(abs(t[pos]))/2.
    return out

def test() :
    spec=np.zeros([200])
    specbin=np.zeros([200])
    spec[50:151]=gauss(np.arange(50,151),100.,99.5,0.78)
    specbin[50:151]=gaussbin(np.arange(50,151),100.,99.5,0.78)
    plt.plot(spec)
    plt.plot(specbin)
    plt.show()
    plt.draw()
    pdb.set_trace()
    peakfit(spec,[95,99,102,107])

def skylines(num,fibers,out=None,plot=None,waveid=None,verbose=False,skyfile='airglow_oct18a') :
    """ Determine positions of skylines in input frame for specified fibers
    """
    #ap=apload.ApLoad(apred='t10')
    #frame=ap.ap1D(num)
    apload.apred='t10'
    print(num)
    frame=apload.ap1D(num)
    if waveid is not None :
        waveframe=apload.apWave(waveid)
        for chip in ['a','b','c'] : frame[chip][4].data = waveframe[chip][2].data
    skylines=ascii.read(os.environ['APOGEE_DIR']+'/data/skylines/'+skyfile+'.txt')
    nlines=len(skylines)
    nfibers=len(fibers)
    linestr = np.zeros(nlines*nfibers,dtype=[
                       ('chip','i4'), ('fiberrow','i4'), ('wave','f4'), ('pixel','f4'),
                       ('dpixel','f4'), ('frameid','i4')
                       ])
    nline=0
    for fiber in fibers :
        for ichip,chip in enumerate(['a','b','c']) :
            medspec = frame[chip][1].data[fiber,:]-medfilt(frame[chip][1].data[fiber,:],101)
            j=np.where(skylines['CHIPNUM'] == ichip+1)[0]
            for iline in j :
                wave=skylines['WAVE'][iline]
                pix0=wave2pix(wave,frame[chip][4].data[fiber,:])
                try :
                    # find peak in median-filtered subtracted spectrum
                    pars=peakfit(medspec,pix0,estsig=2,
                                 sigma=frame[chip][2].data[fiber,:],mask=frame[chip][3].data[fiber,:])
                    linestr['chip'][nline] = ichip+1
                    linestr['fiberrow'][nline] = fiber
                    linestr['wave'][nline] = wave
                    linestr['pixel'][nline] = pars[1]
                    linestr['dpixel'][nline] = pars[1]-pix0
                    linestr['frameid'][nline] = num
                    nline+=1
                    if out is not None :
                        out.write('{:5d}{:5d}{:12.3f}{:12.3f}{:12.3f}{:12d}\n'.format(
                                  ichip+1,fiber,wave,pars[1],pars[1]-pix0,num))
                    elif verbose :
                        print('{:5d}{:5d}{:12.3f}{:12.3f}{:12.3f}{:12d}'.format(
                              ichip+1,fiber,wave,pars[1],pars[1]-pix0,num))
                except :
                    print('failed: ',num,fiber,chip)
    if plot is not None :
        # plot the pixel shift for each chip derived from the airglow lines
        fig,ax = plots.multi(1,1)
        wfig,wax = plots.multi(1,3)
        for ichip in range(3) :
            gd=np.where(linestr['chip'] == ichip+1)[0]
            med=np.median(linestr['dpixel'][gd])
            x = linestr['fiberrow'][gd]
            y = linestr['dpixel'][gd]
            plots.plotp(ax,x,y,color=colors[ichip],xr=[0,300],yr=[med-0.5,med+0.5],
                        size=12,xt='Row',yt='Pixel shift')
            plots.plotc(wax[ichip],linestr['wave'][gd],y,linestr['fiberrow'][gd],zr=[0,300],yr=[med-0.5,med+0.5],
                        xr=xlim[ichip],size=12,xt='Wavelength',yt='Pixel shift')
            gdfit=np.where(np.abs(y-med) < 0.5)[0]
            if len(gdfit) > 1 :
                p=np.polyfit(x[gdfit],y[gdfit],1)
                xx=np.arange(300)
                plots.plotl(ax,xx,p[0]*xx+p[1],color=colors[ichip])
            if waveid : label = 'Frame: {:8d}  Waveid: {:8d}'.format(num,waveid)
            else : label = 'Frame: {:8d}  Delta from ap1dwavecal'.format(num)
            ax.text(0.1,0.9,label,transform=ax.transAxes)
        if type(plot) is str or type(plot) is unicode: 
            wfig.tight_layout()
            wfig.savefig(plot+'_wave.jpg')
            fig.savefig(plot+'.jpg')
        else : pdb.set_trace()
        plt.close('all')
    return linestr

def visit(planfile,out=None,lco=False,waveid=True,skyfile='airglow_oct18a') :
    """ Determine positions of skylines for all frames in input planfile
    """
    if out is not None : f=open(out,'a')
    else : f=None
    p=yanny.yanny(planfile)

    # get the plugmap to get the sky fibers
    plugmjd=p['plugmap'].split('-')[1]
    if lco : 
        apload.instrument='apogee-s'
        plugmap=yanny.yanny(
                os.environ['MAPPER_DATA_2S']+'/'+plugmjd+'/plPlugMapM-'+p['plugmap'].strip("'")+'.par')
    else :
        apload.instrument='apogee-n'
        plugmap=yanny.yanny(
                os.environ['MAPPER_DATA']+'/'+plugmjd+'/plPlugMapM-'+p['plugmap'].strip("'")+'.par')
    skyind=np.where((np.array(plugmap['PLUGMAPOBJ']['objType']) == 'SKY') & 
                   (np.array(plugmap['PLUGMAPOBJ']['holeType']) == 'OBJECT') &
                   (np.array(plugmap['PLUGMAPOBJ']['spectrographId']) == 2) )[0]
    skyfibers = np.array(plugmap['PLUGMAPOBJ']['fiberId'])[skyind]
    skyrows = np.sort(300-skyfibers)
    if p['platetype'].strip("'") == 'sky' : skyrows = np.arange(300)

    # loop over all frames in the planfile and assess skylines in each
    for iframe,frame in enumerate(p['APEXP']['name']) :
        if waveid :
            # use wavelength solution from specified wavecal
            plot = os.path.dirname(planfile)+'/plots/skypixshift-'+frame+'-'+skyfile
            linestr = skylines(int(frame),skyrows,out=f,plot=plot,waveid=int(p['waveid']),skyfile=skyfile)
        else :
            # use existing wavelength solution from ap1D file after ap1dwavecal has been run
            plot = os.path.dirname(planfile)+'/plots/skydeltapixshift-'+frame+'-'+skyfile
            linestr = skylines(int(frame),skyrows,out=f,plot=plot,skyfile=skyfile)

        # Get shifts relative to first frame for each line/fiber
        if iframe == 0 : 
            linestr0 = copy.copy(linestr)
            refnum = int(frame)
        for line in linestr :
            ref = np.where((linestr0['chip'] == line['chip']) & 
                           (linestr0['fiberrow'] == line['fiberrow']) &
                           (linestr0['wave'] == line['wave']))[0]
            if len(ref) > 0 : line['pixel'] -= linestr0['pixel'][ref].mean()
            else : line['pixel'] = -999
        med = np.median(linestr['pixel'])

        # plot shifts relative to first frame, i.e. dithershift via sky lines
        fig,ax=plots.multi(1,1)
        for ichip in range(3) :
            gd = np.where(linestr['chip'] == ichip+1)[0]   
            x = linestr['fiberrow'][gd]
            y = linestr['pixel'][gd]
            plots.plotp(ax,x,y, size=12,xr=[0,300],yr=[med-0.1,med+0.1],
                        xt='Row',yt='Pixel Shift',color=colors[ichip])
            gdfit=np.where(np.abs(y-med) < 0.5)[0]
            if len(gdfit) > 1 :
                pfit=np.polyfit(x[gdfit],y[gdfit],1)
                xx=np.arange(300)
                plots.plotl(ax,xx,pfit[0]*xx+pfit[1],color=colors[ichip])
            label = 'Frame: {:8d}  Waveid: {:8d}'.format(int(frame),refnum)
            ax.text(0.1,0.9,label,transform=ax.transAxes)
        fig.savefig(os.path.dirname(planfile)+'/plots/skydithershift-'+frame+'.jpg')
        plt.close()
        

def scalarDecorator(func):
    """Decorator to return scalar outputs for wave2pix and pix2wave
    """
    @wraps(func)
    def scalar_wrapper(*args,**kwargs):
        if np.array(args[0]).shape == ():
            scalarOut= True
            newargs= (np.array([args[0]]),)
            for ii in range(1,len(args)):
                newargs= newargs+(args[ii],)
            args= newargs
        else:
            scalarOut= False
        result= func(*args,**kwargs)
        if scalarOut:
            return result[0]
        else:
            return result
    return scalar_wrapper

@scalarDecorator
def wave2pix(wave,wave0) :
    """ convert wavelength to pixel given wavelength array
    Args :
       wave : wavelength (\AA) to get pixel of
       wave0 : array with wavelength as a function of pixel number 
    Returns :
       pixel in the chip
    """
    pix0= np.arange(len(wave0))
    # Need to sort into ascending order
    sindx= np.argsort(wave0)
    wave0= wave0[sindx]
    pix0= pix0[sindx]
    # Start from a linear baseline
    baseline= np.polynomial.Polynomial.fit(wave0,pix0,1)
    ip= interpolate.InterpolatedUnivariateSpline(wave0,pix0/baseline(wave0),k=3)
    out= baseline(wave)*ip(wave)
    # NaN for out of bounds
    out[wave > wave0[-1]]= np.nan
    out[wave < wave0[0]]= np.nan
    return out

@scalarDecorator
def pix2wave(pix,wave0) :
    """ convert pixel to wavelength
    Args :
       pix : pixel to get wavelength at
       wave0 : array with wavelength as a function of pixel number 
    Returns :
       wavelength in \AA
    """
    pix0= np.arange(len(wave0))
    # Need to sort into ascending order
    sindx= np.argsort(pix0)
    wave0= wave0[sindx]
    pix0= pix0[sindx]
    # Start from a linear baseline
    baseline= np.polynomial.Polynomial.fit(pix0,wave0,1)
    ip= interpolate.InterpolatedUnivariateSpline(pix0,wave0/baseline(pix0), k=3)
    out= baseline(pix)*ip(pix)
    # NaN for out of bounds
    out[pix < 0]= np.nan
    out[pix > 2047]= np.nan
    return out

def mkpar(mjdstart,mjdend,out=None,lco=False) :
    """ Make calibration file for wavecals between input dates
    """
    # open output file
    if out is not None : f = open(out,'w')
    indiv = []
    sky = []
    for mjd in range(mjdstart,mjdend) :
        # get the files for this dat
        if lco :
            files = sorted(glob.glob(os.environ['APOGEE_DATA_2S']+'/'+str(mjd)+'/*-a-*.apz'))
        else :
            files = sorted(glob.glob(os.environ['APOGEE_DATA']+'/'+str(mjd)+'/*-a-*.apz'))
        print(mjd)
        if len(files) > 3:
            # look for sequences of QUARTZ, THARNE, UNE all at the same dither position
            hdr1 = fits.open(files[0])[1].header
            hdr2 = fits.open(files[1])[1].header
            for ifile in range(2,len(files)-1) :
                try :
                    hdr3 = fits.open(files[ifile])[1].header
                    if ( hdr1['DITHPIX'] == hdr2['DITHPIX'] and hdr1['DITHPIX'] == hdr3['DITHPIX'] and
                         hdr1['LAMPQRTZ'] == 1 and hdr2['LAMPTHAR'] == 1 and hdr3['LAMPUNE'] == 1 ):
                       if out is not None : f.write('wave 99999 99999 {:8s} {:8s},{:8s} {:8s}\n'.format(
                             files[ifile-1].split('-')[2].replace('.apz',''),
                             files[ifile-1].split('-')[2].replace('.apz',''),
                             files[ifile].split('-')[2].replace('.apz',''),
                             files[ifile-2].split('-')[2].replace('.apz','')))
                       indiv.append(files[ifile-1].split('-')[2].replace('.apz','') )
                    if hdr1['IMAGETYP'].strip() == 'Object' and hdr1['NFRAMES']>10 and hdr1['NFRAMES']<40 :
                       sky.append(files[ifile-2].split('-')[2].replace('.apz','') )
                       if out is not None : f.write('lsf 99999 99999 {:8s} {:8s}\n'.format(
                             files[ifile-2].split('-')[2].replace('.apz',''),
                             files[ifile-2].split('-')[2].replace('.apz','')))
                except :
                    pass
                hdr1=hdr2
                hdr2=hdr3
    # now write the multiwave lines for groups of 10 individual wavecals
    for i in range(0,len(indiv),10) :
        name = indiv[i][0:4]+'0000'
        if out is not None : f.write('multiwave 99999 99999 {:8s} {:8s}'.format(name,indiv[i]))
        for j in range(1,10) : 
            try:
                f.write(',{:8s}'.format(indiv[i+j]))
            except :
                pass
        f.write('\n')
    f.close()

def mkallpars(apo=True,lco=True) :
    if lco :
        mjds= [57600, 57966, 58335]
        for i in range(len(mjds)-1) :
            mkpar(mjds[i],mjds[i+1],out='lco-year{:1d}.par'.format(i+1),lco=True)
    if apo :
        mjds= [55800, 56130, 56512, 56876, 57230, 57600, 57966, 58335]
        for i in range(len(mjds)-1) :
            mkpar(mjds[i],mjds[i+1],out='apo-year{:1d}.par'.format(i+1))
        
