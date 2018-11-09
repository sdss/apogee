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
from tools import html
from sdss import yanny
from astropy.table import Table
from pyvista import tv

chips=['a','b','c']
colors=['r','g','b','c','m','y']
xlim=[[16400,17000],[15900,16500],[15100,15800]]

def gauss(x,a,x0,sig) :
    """ Evaluate Gaussian function 
    """
    return a/np.sqrt(2*np.pi)/sig*np.exp(-(x-x0)**2/2./sig**2)

def gaussbin(x,a,x0,sig) :
    """ Evaluate integrated Gaussian function 
    """
    # bin width
    xbin=1.
    t1=(x-x0-xbin/2.)/np.sqrt(2.)/sig
    t2=(x-x0+xbin/2.)/np.sqrt(2.)/sig
    y=(myerf(t2)-myerf(t1))/xbin
    return a*y

def myerf(t) :
    """ Evaluate function that integrates Gaussian from -inf to t
    """
    neg = np.where(t<0.)[0]
    pos = np.where(t>=0.)[0]
    out = t*0.
    out[neg] = erfc(abs(t[neg]))/2.
    out[pos] = 0.5+erf(abs(t[pos]))/2.
    return out


def peakfit(spec,pix0,estsig=5,sigma=None,mask=None,plot=False,func=gaussbin) :
    """ Return integrated-Gaussian centers near input pixel center
    
    Args:
        spec (float) : data spectrum arraty
        pix0 (float) : initial pixel guess
        estsig (float ) : initial guess for window width=5*estsig (default=5)
        sigma (float)  : uncertainty array (default=None)
        mask (float)  : mask array (default=None), NOT CURRENTLY IMPLEMENT
        plot (bool) : plot spectrum and fit in current plot window (default=False)
        func (function) : user-supplied function to use to fit (default=gaussbin)
    """
    x = np.arange(len(spec))
    cen = int(round(pix0))
    sig=estsig
    back=0.
    for iter in range(11) :
        # window width to search
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

def test() :
    """ test routine for peakfity
    """
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

def func_multi_poly(x,*pars) :
    """ Convert pixel to wavelength using wavecal parameters
          w = poly(x + offset(group,chip))
          pars = [npoly coefficients, ngroup*3 chip offsets]
        Args:  
            x (float) : [npts,3] array of (pixel,chip,group)
         pars (float) : input parameter array

        Returns :
         wave (float) : wavelength array for input pixel(s), parameters
    """
    wave=np.zeros(x.shape[1])
    ngroup = len(set(x[2,:]))
    nchip = 3
    npoly = len(pars)-ngroup*nchip
    coef = pars[0:npoly]
    # loop over all chip/group combinations
    for ichip in range(nchip) :
        for igroup in range(ngroup) :
            offset = pars[npoly+igroup*nchip+ichip]
            j=np.where((x[1,:] == ichip+1) & (np.round(x[2,:]).astype(int) == igroup))[0]
            xglobal = x[0,j] - 1023.5 + (ichip-1)*2048 + offset
            wave[j] = np.polyval(coef,xglobal)
    return wave

def findlines(frame,rows,waves,lines,out=None,verbose=False,estsig=2) :
    """ Determine positions of lines from input file in input frame for specified rows

    Args:
        frame (dict) : dictionary with ['a','b','c'] keys for each chip containing HDULists with flux, error, and mask
        rows (list) : list of rows to look for lines in
        waves (list)  : list of wavelength arrays to be used to get initial pixel guess for input lines
        lines :  table with desired lines, must have at least CHIPNUM and WAVE tags
        out= (str) : optional name of output ASCII file for lines (default=None)

    Returns :
        structure with identified lines, with tags chip, row, wave, peak, pixrel, dpixel, frameid
    """
    num=int(os.path.basename(frame['a'][0].header['FILENAME']).split('-')[1])
    nlines=len(lines)
    nrows=len(rows)
    linestr = np.zeros(nlines*nrows,dtype=[
                       ('chip','i4'), ('row','i4'), ('wave','f4'), ('peak','f4'), ('pixel','f4'),
                       ('dpixel','f4'), ('frameid','i4')
                       ])
    nline=0
    for row in rows :
        for ichip,chip in enumerate(['a','b','c']) :
            medspec = frame[chip][1].data[row,:]-medfilt(frame[chip][1].data[row,:],101)
            j=np.where(lines['CHIPNUM'] == ichip+1)[0]
            for iline in j :
                wave=lines['WAVE'][iline]
                pix0=wave2pix(wave,waves[chip][row,:])
                try :
                    # find peak in median-filtered subtracted spectrum
                    pars=peakfit(medspec,pix0,estsig=estsig,
                                 sigma=frame[chip][2].data[row,:],mask=frame[chip][3].data[row,:])
                    linestr['chip'][nline] = ichip+1
                    linestr['row'][nline] = row
                    linestr['wave'][nline] = wave
                    linestr['peak'][nline] = pars[0]
                    linestr['pixel'][nline] = pars[1]
                    linestr['dpixel'][nline] = pars[1]-pix0
                    linestr['frameid'][nline] = num
                    nline+=1
                    if out is not None :
                        out.write('{:5d}{:5d}{:12.3f}{:12.3f}{:12.3f}{:12.3f}{:12d}\n'.format(
                                  ichip+1,row,wave,pars[0],pars[1],pars[1]-pix0,num))
                    elif verbose :
                        print('{:5d}{:5d}{:12.3f}{:12.3f}{:12.3f}{:12.3f}{:12d}'.format(
                              ichip+1,row,wave,pars[0],pars[1],pars[1]-pix0,num))
                except :
                    if verbose : print('failed: ',num,row,chip)
    return linestr[0:nline]

def wavecal(nums=[2420038],name=None,vers='t9',inst='apogee-n',rows=[150],npoly=4,reject=3,plot=False,hard=None,verbose=False,clobber=False,init=False) :
    """ APOGEE wavelength calibration

    Solves for wavelength calibration given input frameid(s) allowing for a single polynomial
    wavelength solution with offsets for each chip and each group of wavecals

    Args:
        nums (list of ints): list of input frame ids
        rows (list of ints): list of rows to get solutions at
        plot verbose (bool)  plot fits (default=False)
        verbose (bool) :  plot fits (default=False)
    """
    load=apload.ApLoad(apred=vers,instrument=inst)

    if name is None : name = nums[0]
    # Initial guess for wavelengths, used to find lines
    coef0 = {}
    if inst == 'apogee-n' :
        coef0['a'] = np.flip([ 16955.45703, -0.2128979266, -1.117692409e-05])
        coef0['b'] = np.flip([ 16434.20508, -0.2613874376, -1.035568130e-05])
        coef0['c'] = np.flip([ 15809.69238, -0.3065520823, -9.610030247e-06])
        pars0 = [1.19112154e-10,-1.03229705e-05,-2.82422124e-01,1.61568801e+04,
                -1.44043125e+02,0.00000000e+00,1.54456264e+02]
    else :
        coef0['a'] = np.flip([  16957.7252,  -2.14859462e-01,  -1.09959211e-05])
        coef0['b'] = np.flip([  16432.4720,  -2.63317139e-01,  -1.03074667e-05])
        coef0['c'] = np.flip([  15802.3346,  -3.08933509e-01,  -9.45618858e-06])
        pars0 = [ 1.19138048e-10,-1.03101159e-05,-2.84129914e-01,1.61531282e+04,
                 -1.49566344e+02,-9.99666738e-11, 1.58164526e+02]
    waves = {}
    pixels = np.arange(2048)
    for chip in chips : waves[chip] = np.tile(np.polyval(coef0[chip],pixels),(300,1))

    # find lines
    print('finding lines with initial wavelength guess: ',coef0)
    ngroup=1
    for inum,num in enumerate(nums) :
        # load 1D frame
        frame = load.ap1D(num)
        out = load.filename('Wave',num=num,chips=True)
        if frame is not None and frame != 0 :
            # get correct arclines
            if frame['a'][0].header['LAMPUNE'] : lampfile = 'UNe.vac.apogee'
            if frame['a'][0].header['LAMPTHAR'] : lampfile = 'tharne.lines.vac.apogee'
            arclines=ascii.read(os.environ['APOGEE_DIR']+'/data/arclines/'+lampfile)
            j=np.where(arclines['USEWAVE'])[0]
            arclines=arclines[j]
            # find lines or use previous found lines
            linesfile=out.replace('Wave','Lines')
            if os.path.exists(linesfile) and not clobber :
                print('Reading existing Lines data',num)
                flinestr = fits.open(linesfile)[1].data
            elif os.path.exists(linesfile.replace('as','ap')) :
                os.rename(linesfile.replace('as','ap'),linesfile)
                flinestr = fits.open(linesfile)[1].data
            else :
                print('Finding lines: ', num)
                flinestr = findlines(frame,rows,waves,arclines,verbose=verbose,estsig=2)
                Table(flinestr).write(linesfile,overwrite=True)
            # replace frameid tag with group identification, which must start at 0 for func_multi_poly indexing
            if inum > 0 and abs(num-nums[inum-1]) > 1 : ngroup +=1
            flinestr['frameid'] = ngroup-1
            if inum == 0 : linestr = flinestr
            else : linestr = np.append(linestr,flinestr)
            print(' Frame: {:d}  Nlines: {:d}  '.format(num,len(flinestr)))
        else :
            print('Error reading frame: ', num)

    out=load.filename('Wave',num=name,chips=True).replace('Wave','PWave')
    # do the wavecal fit
    # initial parameter guess for first row, subsequent rows will use guess from previous row
    npars=npoly+3*ngroup
    pars = np.zeros(npars)
    # initial quadratic relation from chip b and initial chip offsets or better guess if we have it
    if init :pars[npoly-3:npoly] = coef0['b']
    else : 
        pars[npoly-4:npoly] = pars0[0:4]
        for igroup in range(ngroup): pars[npoly+igroup*3:npoly+(igroup+1)*3] = pars0[4:7]

    # set up output arrays for all 300 fibers
    allpars=np.zeros([npars,300])
    rms=np.zeros(300)
    sig=np.zeros(300)
    if plot : 
        fig,ax=plots.multi(1,3,hspace=0.001,wspace=0.001)
        fig2,ax2=plots.multi(1,4,hspace=0.001,wspace=0.001)
    # loop over requested rows
    for irow,row in enumerate(rows) :
        # position of green chip in first group fixed to 0., otherwise no bounds
        bounds = ( np.zeros(len(pars))-np.inf, np.zeros(len(pars))+np.inf)
        bounds[0][npoly+1] =  -1.e-10
        bounds[1][npoly+1] =  1.e-10
        # set up independent variable array with pixel, chip, groupid, and dependent variable (wavelength)
        thisrow = np.where(linestr['row'] == row)[0]
        x = np.zeros([3,len(thisrow)])
        x[0,:] = linestr['pixel'][thisrow]
        x[1,:] = linestr['chip'][thisrow]
        x[2,:] = linestr['frameid'][thisrow]
        y = linestr['wave'][thisrow]
        res = y-func_multi_poly(x,*pars)
        if irow == 0 : maxiter=15
        else : maxiter=7
        for niter in range(maxiter) :
            # reject lines that have bad residuals
            gd = np.where(abs(res) < np.median(res)+reject*np.median(np.abs(res)))[0]
            # use curve_fit to optimize paramtyers
            try :
                popt,pcov = curve_fit(func_multi_poly,x[:,gd],y[gd],p0=pars,bounds=bounds)
                res = y-func_multi_poly(x,*pars)
                rms[row] = res[gd].std()
                sig[row] = np.median(np.abs(res[gd]))
                if verbose: print(niter,row,len(gd),np.median(res),np.median(np.abs(res)),res[gd].std())
                pars = popt
            except :
                print('Solution failed for row: ', row)
                popt = pars*0.
        # plot individual line residuals if requested
        if plot :
            for ichip in range(3) :
                plt = np.where(x[1,gd] == ichip+1)[0]
                z=np.zeros(len(plt))+row
                plots.plotc(ax[ichip],x[0,gd[plt]],y[gd[plt]]-func_multi_poly(x[:,gd[plt]],*popt),z,zr=[0,300],
                            xt='Pixel',yt='obs-fit wavelength',size=10)
        if verbose : print(row,pars)
        allpars[:,row] = popt
    # plot rms/sig and chip locations if requested
    if plot :
        plots.plotp(ax2[0],np.arange(300),rms,color='r',size=10,xt='Row',yt='rms/sig')
        plots.plotp(ax2[0],np.arange(300),sig,color='g',size=10)
        for ichip in range(3) : 
            for igroup in range(ngroup) :
                y=allpars[npoly+igroup*3+ichip,:]
                plots.plotp(ax2[1+ichip],np.arange(300),y,yr=[np.median(y)-2,np.median(y)+2],color=colors[ichip],
                            size=10,yt='chip location')
        fig.suptitle(name)
        fig2.suptitle(name)
        if hard :
            fig.savefig(os.path.dirname(out)+'/plots/'+os.path.basename(out).replace('.fits','.jpg'))
            fig2.savefig(os.path.dirname(out)+'/plots/'+os.path.basename(out).replace('.fits','_sig.jpg'))

    # output the apWavecal files with correct format
    x = np.zeros([3,2048])
    for ichip,chip in enumerate(chips) :
        hdu=fits.HDUList()
        x[0,:] = pixels
        x[1,:] = ichip+1
        x[2,:] = 0
        chippars=np.zeros([14,300])
        chipwaves=np.zeros([300,2048])
        for row in rows :
            pow=npoly-1-np.arange(npoly)
            polypars=allpars[0:npoly,row]*3000**pow
            chippars[:,row] = np.append([ -1023.5+(ichip-1)*2048 + allpars[npoly+ichip,row],0., 0., 1., 0., 0.], 
                              np.flip( np.append(np.zeros(8-npoly),polypars)))
            chipwaves[row,:] = func_multi_poly(x,*allpars[:,row])
        hdu.append(fits.PrimaryHDU())
        hdu.append(fits.ImageHDU(chippars))
        hdu.append(fits.ImageHDU(chipwaves))
        hdu.append(fits.ImageHDU(allpars))
        hdu.writeto(out.replace('Wave','Wave-'+chip),overwrite=True)

def visit(planfile,out=None,lco=False,waveid=True,skyfile='airglow_oct18a',vers='t9') :
    """ Determine positions of skylines for all frames in input planfile
    """
    if out is not None : f=open(out,'a')
    else : f=None
    p=yanny.yanny(planfile)

    if lco : inst='apogee-s'
    else : inst='apogee-n'
    load=apload.ApLoad(apred=vers,instrument=inst)

    # get the plugmap to get the sky fibers
    plugmjd=p['plugmap'].split('-')[1]
    if lco : 
        plugmap=yanny.yanny(
                os.environ['MAPPER_DATA_2S']+'/'+plugmjd+'/plPlugMapM-'+p['plugmap'].strip("'")+'.par')
    else :
        plugmap=yanny.yanny(
                os.environ['MAPPER_DATA']+'/'+plugmjd+'/plPlugMapM-'+p['plugmap'].strip("'")+'.par')
    skyind=np.where((np.array(plugmap['PLUGMAPOBJ']['objType']) == 'SKY') & 
                   (np.array(plugmap['PLUGMAPOBJ']['holeType']) == 'OBJECT') &
                   (np.array(plugmap['PLUGMAPOBJ']['spectrographId']) == 2) )[0]
    skyfibers = np.array(plugmap['PLUGMAPOBJ']['fiberId'])[skyind]
    skyrows = np.sort(300-skyfibers)
    if p['platetype'].strip("'") == 'sky' : skyrows = np.arange(300)
    skylines=ascii.read(os.environ['APOGEE_DIR']+'/data/skylines/'+skyfile+'.txt')

    # loop over all frames in the planfile and assess skylines in each
    for iframe,name in enumerate(p['APEXP']['name']) :
        frame = load.ap1D(int(name))
        waves=[]
        if waveid :
            # use wavelength solution from specified wavecal
            waveframe=load.apWave(waveid)
            for chip in chips : waves[chip] = waveframe[chip][2].data
            plot = os.path.dirname(planfile)+'/plots/skypixshift-'+frame+'-'+skyfile
        else :
            # use existing wavelength solution from ap1D file after ap1dwavecal has been run
            for chip in chips : waves[chip] = frame[chip][4].data
            plot = os.path.dirname(planfile)+'/plots/skydeltapixshift-'+frame+'-'+skyfile
        linestr = findlines(frame,skyrows,waves,skylines,out=f,plot=plot)
        if plot is not None :
            # plot the pixel shift for each chip derived from the airglow lines
            fig,ax = plots.multi(1,1)
            wfig,wax = plots.multi(1,3)
            for ichip in range(3) :
                gd=np.where(linestr['chip'] == ichip+1)[0]
                med=np.median(linestr['dpixel'][gd])
                x = linestr['row'][gd]
                y = linestr['dpixel'][gd]
                plots.plotp(ax,x,y,color=colors[ichip],xr=[0,300],yr=[med-0.5,med+0.5],
                            size=12,xt='Row',yt='Pixel shift')
                plots.plotc(wax[ichip],linestr['wave'][gd],y,linestr['row'][gd],zr=[0,300],yr=[med-0.5,med+0.5],
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

        # Get shifts relative to first frame for each line/fiber
        if iframe == 0 : 
            linestr0 = copy.copy(linestr)
            refnum = int(frame)
        for line in linestr :
            ref = np.where((linestr0['chip'] == line['chip']) & 
                           (linestr0['row'] == line['row']) &
                           (linestr0['wave'] == line['wave']))[0]
            if len(ref) > 0 : line['pixel'] -= linestr0['pixel'][ref].mean()
            else : line['pixel'] = -999
        med = np.median(linestr['pixel'])

        # plot shifts relative to first frame, i.e. dithershift via sky lines
        fig,ax=plots.multi(1,1)
        for ichip in range(3) :
            gd = np.where(linestr['chip'] == ichip+1)[0]   
            x = linestr['row'][gd]
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
       wave(s) : wavelength(s) (\AA) to get pixel of
       wave0 : array with wavelength as a function of pixel number 
    Returns :
       pixel(s) in the chip
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
    """ convert pixel(s) to wavelength(s)
    Args :
       pix : pixel(s) to get wavelength at
       wave0 : array with wavelength as a function of pixel number 
    Returns :
       wavelength(s) in \AA
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

def compare(npoly=4) :

    files=glob.glob('apPWave-b-*.fits')
    files.sort()
    dates=[]
    for file in files :
        dates.append(int(file.split('-')[2].replace('.fits',''))/10000)
    dates=np.array(dates)
    files=np.array(files)

    w=np.arange(15160.,16900.,50.)
    x=np.arange(-1024-2048-150,1024+2048+150,25)
    grid=[]
    for year in range(0,7,7) :
        print('year: ', year)
        i1 = 55757+year*365-55562
        i2 = i1+365*7
        j = np.where((dates >=i1) & (dates<=i2) & ((dates<2430) | (dates>2450)) )[0]
        wave=np.zeros([len(x),len(j)])
        pix=np.zeros([len(w),len(j)])
        chipa=np.zeros([300,len(j)])
        chipc=np.zeros([300,len(j)])
        fig,ax=plots.multi(1,3)
        for ifile,file in enumerate(files[j]) :
            print(file)
            a=fits.open(file)[3].data
            wave[:,ifile]=np.polyval(a[0:npoly,150],x)
            pix[:,ifile]=wave2pix(w,wave[:,ifile])
            chipa[:,ifile]=a[npoly]
            chipc[:,ifile]=a[npoly+2]

        wmed = np.median(wave,axis=1)
        pmed = np.nanmedian(pix,axis=1)
        chipamed = np.median(chipa,axis=1)
        chipcmed = np.median(chipc,axis=1)
        for ifile,file in enumerate(files[j]) :
            pix[:,ifile]-=pmed
            wave[:,ifile]-=wmed
            wave[:,ifile]-=np.median(wave[:,ifile])
            pix[:,ifile]-=np.median(pix[:,ifile])
            chipa[:,ifile]-=chipamed
            chipc[:,ifile]-=chipcmed
            plots.plotl(ax[0],x,wave[:,ifile],yr=[-0.5,0.5],xt='global pixel',yt=r'$\lambda-\lambda_{med}$')
            plots.plotl(ax[1],np.arange(300),chipa[:,ifile],yr=[-0.5,0.5],xt='Fiber',yt='chipa-chipamed')
            plots.plotl(ax[2],np.arange(300),chipc[:,ifile],yr=[-0.5,0.5],xt='Fiber',yt='chipc-chipcmed')
        fig.suptitle('Year: {:d}'.format(year))
        fig.savefig('year{:1d}.jpg'.format(year))
        plt.close()
        t=tv.TV()
        t.cmap='viridis'
        t.tv(wave,min=-0.1,max=0.05)
        t.fig.suptitle('Year: {:d}'.format(year))
        t.fig.savefig('year{:1d}wave.jpg'.format(year))
        t.ax.cla()
        t.tv(pix,min=-0.1,max=0.05)
        t.fig.suptitle('Year: {:d}'.format(year))
        t.fig.savefig('year{:1d}pix.jpg'.format(year))
        t.ax.cla()
        t.tv(chipa,min=-0.1,max=0.05)
        t.fig.suptitle('Year: {:d}'.format(year))
        t.fig.savefig('year{:1d}chipa.jpg'.format(year))
        t.ax.cla()
        t.tv(chipc,min=-0.1,max=0.05)
        t.fig.suptitle('Year: {:d}'.format(year))
        t.fig.savefig('year{:1d}chipc.jpg'.format(year))
        plt.close()
        root='year{:1d}'.format(year)
        grid.append([root+'.jpg',root+'wave.jpg',root+'pix.jpg',root+'chipa.jpg',root+'chipc.jpg'])
    html.htmltab(grid,file='apogee-n.html')


