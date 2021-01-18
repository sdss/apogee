# encoding: utf-8
#
# @Author: Jon Holtzman, some routines based off of IDL routines of David Nidever
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
from scipy.signal import medfilt, convolve, boxcar
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

def wavecal(nums=[2420038],name=None,vers='current',inst='apogee-n',rows=[150],npoly=4,reject=3,pixellim=[3,2043],
            plot=False,hard=True,verbose=False,clobber=False,init=False,nofit=False,test=False) :
    """ APOGEE wavelength calibration

    Solves for wavelength calibration given input frameid(s) allowing for a single polynomial
    wavelength solution with offsets for each chip and each group of wavecals

    Args:
        nums (list of ints): list of input frame ids
        name (int) : name for output ID (default = first frame from list of nums)
        vers (str) : apred version (defaut='current')
        inst (str) : instrument (default='apogee-n')
        rows (list of ints): list of rows to get solutions at
        npoly (int) : order of polynomial fit (default=4)
        reject (float) : factor for line rejection
        plot verbose (bool)  plot fits (default=False)
        verbose (bool) :  plot fits (default=False)
        clobber (bool) : forces remeasuring lines
        init (bool) : if true, use simple quadratic estimate for first guess (default=False)
        nofit (bool) : only find lines, skip fit (default=False)
        test (bool) : if True, use polynomial from first set, only let centers shift for all other sets
    """
    load=apload.ApLoad(apred=vers,instrument=inst)

    if name is None : name = nums[0]
    if test : name = int(name/10000)*10000+9999
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

    # find lines
    print('finding lines with initial wavelength guess: ',coef0)
    waves = {}
    pixels = np.arange(2048)
    for chip in chips : waves[chip] = np.tile(np.polyval(coef0[chip],pixels),(300,1))
    maxgroup=1
    frames=[]
    for inum,num in enumerate(nums) :
        # load 1D frame
        frame = load.ap1D(num)
        out = load.filename('Wave',num=num,chips=True)
        print(num,frame)
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
            else :
                print('Finding lines: ', num)
                flinestr = findlines(frame,rows,waves,arclines,verbose=verbose,estsig=1)
                Table(flinestr).write(linesfile,overwrite=True)
            # replace frameid tag with group identification, which must start at 0 for func_multi_poly indexing
            if inum > 0 and abs(num-nums[inum-1]) > 1 : maxgroup +=1
            flinestr['frameid'] = maxgroup-1
            if inum == 0 : linestr = flinestr
            else : linestr = np.append(linestr,flinestr)
            print(' Frame: {:d}  Nlines: {:d}  '.format(num,len(flinestr)))
            frames.append(num)
        else :
            print('Error reading frame: ', num)

    if nofit : return

    # do the wavecal fit
    # initial parameter guess for first row, subsequent rows will use guess from previous row
    npars=npoly+3*maxgroup
    pars = np.zeros(npars)

    # initial quadratic relation from chip b and initial chip offsets or better guess if we have it
    if init :pars[npoly-3:npoly] = coef0['b']
    else : 
        pars[npoly-4:npoly] = pars0[0:4]
        for igroup in range(maxgroup): pars[npoly+igroup*3:npoly+(igroup+1)*3] = pars0[4:7]
    initpars=copy.copy(pars)
 
    # set up output arrays for all 300 fibers
    allpars=np.zeros([npars,300])
    chipa=np.zeros([300,maxgroup])
    chipb=np.zeros([300,maxgroup])
    chipc=np.zeros([300,maxgroup])
    rms=np.zeros([300,maxgroup])
    sig=np.zeros([300,maxgroup])
    if plot : 
        fig,ax=plots.multi(1,3,hspace=0.001,wspace=0.001)
        fig2,ax2=plots.multi(1,3,hspace=0.001,wspace=0.001)

    # loop over requested rows
    for irow,row in enumerate(rows) :
        # set up independent variable array with pixel, chip, groupid, and dependent variable (wavelength)
        thisrow = np.where((linestr['row'] == row) & (linestr['peak'] > 100) & 
                           (linestr['pixel']>pixellim[0]) & (linestr['pixel']<pixellim[1]) )[0]
        x = np.zeros([3,len(thisrow)])
        x[0,:] = linestr['pixel'][thisrow]
        x[1,:] = linestr['chip'][thisrow]
        # we may have missing groups for this row
        groupid,groups = getgroup(linestr['frameid'][thisrow])
        x[2,:] = groupid
        #x[2,:] = linestr['frameid'][thisrow]
        ngroup = len(groups)
        y = linestr['wave'][thisrow]
        # if we don't have any groups, skip this row
        if ngroup<= 0: continue

        npars=npoly+3*ngroup
        pars = np.zeros(npars)
        # if we have more than one group, get starting polynomial guess from first group, to help
        #   to avoid local minima
        if ngroup > 1 :
            if abs(nums[1]-nums[0]) > 1 : 
                print('for multiple groups, first two frames must be from same group!')
                pdb.set_trace()
            print('running wavecal for: ', nums[0:2],maxgroup,ngroup,' row: ',row)
            pars0 = wavecal(nums=nums[0:2],name=None,vers=vers,inst=inst,rows=[row],npoly=npoly,reject=reject,init=init,verbose=verbose)
            pars[npoly-4:npoly] = pars0[0:4]
            for igroup in range(ngroup): pars[npoly+igroup*3:npoly+(igroup+1)*3] = pars0[4:7]
        else :
            pars = copy.copy(initpars)
        # initial residuals
        res = y-func_multi_poly(x,*pars)
        # iterate to allow outlier rejection
        if irow == 0 : maxiter=7
        else : maxiter=7
        for niter in range(maxiter) :
            # initialize bounds (to no bounds)
            bounds = ( np.zeros(len(pars))-np.inf, np.zeros(len(pars))+np.inf)
            # lock the middle chip position if we have one group, else the central wavelength
            if ngroup==1 :
                bounds[0][npoly+1] = -1.e-7
                bounds[1][npoly+1] = 1.e-7
            else :
                bounds[0][npoly-1] = pars[npoly-1]-1.e-7*abs(pars[npoly-1])
                bounds[1][npoly-1] = pars[npoly-1]+1.e-7*abs(pars[npoly-1])
                # if we have multiple groups, only fit for chip locations during first iterations and every 3rd thereafter
                if test or niter<3 or niter%3 == 1 : 
                    for i in range(npoly) : 
                        bounds[0][i] =  pars[i]-1.e-7*abs(pars[i])
                        bounds[1][i] =  pars[i]+1.e-7*abs(pars[i])

            # reject lines with bad residuals
            gd = np.where(abs(res) < np.median(res)+reject*np.median(np.abs(res)))[0]

            # calculate rms in each group and make sure we have lines in all groups
            for igroup in range(ngroup) : 
                j=np.where(x[2,gd] == igroup)[0]
                if len(j) <= 0 :
                  # if any group is missing, things will fail in func_multi_poly to determine correct npoly
                  print('missing lines from group: ', igroup)
                  pdb.set_trace()
            
            # use curve_fit to optimize parameters
            try :
                if verbose: 
                    print('Niter: ', niter, 'row: ', row, 'ngroup: ', ngroup, 'nlines: ', len(thisrow), 'gd: ', len(gd))
                    print(pars)
                popt,pcov = curve_fit(func_multi_poly,x[:,gd],y[gd],p0=pars,bounds=bounds)
                pars = copy.copy(popt)
                res = y-func_multi_poly(x,*pars)
                if verbose: 
                    print('res: ',len(gd),np.median(res),np.median(np.abs(res)),res[gd].std())
                    print(pars)
            except :
                print('Solution failed for row: ', row)
                pdb.set_trace()
                popt = pars*0.
            # plot individual line residuals if requested. Get rid of maxiter-1 if you want to see plots at each iteration
            if plot and niter == maxiter-1 :
                for ichip in range(3) :
                    gdplt = np.where(x[1,gd] == ichip+1)[0]
                    if niter == maxiter-1 : 
                        z=np.zeros(len(y))+row
                        zr=[0,300]
                    else : 
                        z=x[2,:]
                        zr=[0,ngroup]
                    if niter < maxiter-1 : ax[ichip].cla()
                    plots.plotc(ax[ichip],x[0,gd[gdplt]],res[gd[gdplt]],z[gd[gdplt]],zr=zr,
                                xt='Pixel',yt='obs-fit wavelength',size=10)
                    plt.show()
                if not hard : pdb.set_trace()
        
        allrms = res[gd].std()
        # throw out bad solutions
        #if allrms > 0.1 : popt = pars*0.
        if allrms < 0.1 : initpars = copy.copy(pars)
        if verbose : print(row,pars)
        allpars[0:npoly,row] = popt[0:npoly]
        ref=allpars[npoly+1,row]
        # save final fits in allpars. For chip locations, transform to chip offsets
        for jgroup in range(ngroup) :
            igroup=groups[jgroup]
            for ichip in range(3) : allpars[npoly+igroup*3+ichip,row] = popt[npoly+jgroup*3+ichip]
            j=np.where(x[2,gd] == igroup)[0]
            rms[row,igroup] = res[gd[j]].std()
            sig[row,igroup] = np.median(np.abs(res[gd[j]]))

    # now refine the solution by averaging zeropoint across all groups and
    # by fitting across different rows to require a smooth solution
    if ngroup > 1 : newpars,newwaves = refine(allpars)
    else : newpars = allpars

    # save results in apWave fies
    out=load.filename('Wave',num=name,chips=True)   #.replace('Wave','PWave')
    save_apWave(newpars,out=out,npoly=npoly,rows=rows,frames=frames,rms=rms,sig=sig,allpars=allpars)

    if plot : 
        plot_apWave([name],apred=vers,inst=inst,hard=hard)
        # individual lines from last row
        #if hard :
        #    try : os.mkdir(os.path.dirname(root))
        #    except : pass
        #    fig.savefig(root+'.png')

    return pars

def plot_apWave(nums,apred='current',inst='apogee-n',out=None,hard=False) :

  """ Diagnostic plots of wavecal
  """
  mjd=[]
  wfit=[]
  grid=[]
  yt=[]
  for num in nums :
    load=apload.ApLoad(apred=apred,instrument=inst)
    wave=load.apWave(num)
    outname=load.filename('Wave',num=num,chips=True)
    allpars=wave['a'][6].data
    rms=wave['a'][4].data
    sig=wave['a'][5].data
    ngroup=int(wave['a'][0].header['NGROUP'])
    npoly=wave['a'][0].header['NPOLY']

    name=os.path.basename(outname)
    fig2,ax2=plots.multi(1,3,hspace=0.001,wspace=0.001)
    # diagnostic plots
    # for plots, transform absolute chip location to relative to middle chip
    for jgroup in range(ngroup) :
        #igroup=groups[jgroup]
        igroup=jgroup
        for ichip in [0,2] : 
          for row in np.arange(300) : 
              allpars[npoly+igroup*3+ichip,row] -= allpars[npoly+igroup*3+1,row]
    for ichip in range(3) : 
        for igroup in range(ngroup) :
            y=allpars[npoly+igroup*3+ichip,:]
            gdlim=np.where(np.abs(y) > 0.0001)[0] 
            plots.plotc(ax2[ichip],np.arange(300),y,np.zeros(300)+igroup,yr=[np.median(y[gdlim])-5,np.median(y[gdlim])+5],
                        zr=[0,ngroup], size=10,yt='chip location')
    fig2.suptitle(name)
    xtit=['Individual lines','chip locations','chip locations','rms  and sig']
    root = os.path.dirname(outname)+'/plots/'+os.path.basename(outname).replace('.fits','')
    rootname = os.path.basename(root)
    if hard :
        try : os.mkdir(os.path.dirname(root))
        except : pass
        fig2.savefig(root+'_chiploc.png')
        plt.close()

    # summary figure of chip locations
    fig,ax=plots.multi(1,4,hspace=0.001)
    cb_ax=fig.add_axes((0.9,0.72,0.03,0.15))
    cb_ax2=fig.add_axes((0.9,0.15,0.03,0.4))
    # get chip positions relative to median postion across all groups
    chipa=allpars[4:200:3,:]-np.median(allpars[4:200:3,:],axis=0)
    chipc=allpars[6:200:3,:]-np.median(allpars[6:200:3,:],axis=0)
    chipb=allpars[5:200:3,:]-np.median(allpars[5:200:3,:],axis=0)
    # image of chip b shifts
    aximage=ax[0].imshow(chipb,vmin=-2,vmax=2,cmap='viridis',interpolation='nearest',aspect='auto')
    ax[0].set_ylabel('chip loc')
    fig.colorbar(aximage,cax=cb_ax,orientation='vertical')

    # get chip b shift relative to median across all rows
    chipb=(chipb.T-np.median(chipb,axis=1)).T
    vmin=-0.07
    vmax=0.07
    ax[1].imshow(chipb,vmin=vmin,vmax=vmax,cmap='viridis',interpolation='nearest',aspect='auto')
    ax[1].set_ylabel('rel chip loc')
    # chip gaps 
    ax[2].imshow(chipa,vmin=vmin,vmax=vmax,cmap='viridis',interpolation='nearest',aspect='auto')
    ax[2].set_ylabel('g-r gap')
    aximage=ax[3].imshow(chipc,vmin=vmin,vmax=vmax,cmap='viridis',interpolation='nearest',aspect='auto')
    ax[3].set_xlabel('Row')
    ax[3].set_ylabel('b-g gap')
    fig.suptitle(rootname)
    fig.colorbar(aximage,cax=cb_ax2,orientation='vertical')
    if hard: 
        fig.savefig(root+'_sum.png')
        plt.close()

    if rms is not None :
        fig,ax=plots.multi(1,2,hspace=0.5)
        cb_ax=fig.add_axes((0.9,0.6,0.03,0.3))
        cb_ax2=fig.add_axes((0.9,0.1,0.03,0.3))
        aximage=ax[0].imshow(rms,vmin=0.,vmax=0.1,cmap='viridis',interpolation='nearest',aspect='auto')
        fig.colorbar(aximage,cax=cb_ax,orientation='vertical')
        aximage=ax[1].imshow(sig,vmin=0.,vmax=0.05,cmap='viridis',interpolation='nearest',aspect='auto')
        fig.colorbar(aximage,cax=cb_ax2,orientation='vertical')
        if hard : 
            fig.savefig(root+'_rms.png')
            plt.close()

    wgroup=[]
    for group in range(1,ngroup ) :
        for ichip in [0,2] : 
            allpars[npoly+group*3+ichip,:] += allpars[npoly+group*3+1,:]
        frame=wave['a'][0].header['FRAME{:d}'.format(group*2)]
        # solve for 4 parameter fit to dpixel, with linear trend with row, plus 2 chip offsets
        design=np.zeros([900,4])
        y=np.zeros(900)
        # offset of each chip
        for ichip in range(3) :
            # global slope with rows
            design[ichip*300+np.arange(300),0] = np.arange(300)
            design[ichip*300+np.arange(300),ichip+1] = 1.
            y[ichip*300+np.arange(300)] = allpars[npoly+group*3+ichip,:]-allpars[npoly+ichip,:]
        mjd.append(frame//10000+55562)
        # reject bad fibers
        gd=np.where(abs(y) > 1.e-5)[0]
        design=design[gd,:]
        y=y[gd]
        # solve
        try : 
            w = np.linalg.solve(np.dot(design.T,design), np.dot(design.T, y))
            wgroup.append(w)
        except : 
            print('fit failed ....')
            pdb.set_trace()
        print('fit parameters: ', w)
    # subtract median parameter value for this wavecal, so we can compare across different wavecals
    wgroup=np.array(wgroup)
    for i in range(4) : wgroup[:,i]-=np.median(wgroup[:,i])
    wfit.extend(wgroup)

    grid.append(['../plots/'+rootname+'.png','../plots/'+rootname+'_chiploc.png','../plots/'+rootname+'_sum.png','../plots/'+rootname+'_rms.png'])
    yt.append('{:08d}'.format(num))

  mjd=np.array(mjd)
  wfit=np.array(wfit)
  fig,ax=plots.multi(1,4,hspace=0.001)
  plots.plotp(ax[0],mjd,wfit[:,0],yr=[-1.e-3,1.e-3],yt='slope')
  plots.plotp(ax[1],mjd,wfit[:,2],yr=[-5,5],yt='g')
  plots.plotp(ax[2],mjd,wfit[:,1],yr=[-0.2,0.2],yt='r-g')
  plots.plotp(ax[3],mjd,wfit[:,3],yr=[-0.2,0.2],yt='b-g')
  if hard :
      plt.close()
      if out is None : root = os.path.dirname(outname)+'/plots/'+os.path.basename(outname).replace('.fits','')
      else : root = os.path.dirname(outname)+'/plots/'+out
      fig.savefig(root+'_history.png')
      if out is None : root = os.path.dirname(outname)+'/html/'+os.path.basename(outname).replace('.fits','')
      else : root = os.path.dirname(outname)+'/html/'+out
      try : os.mkdir(os.path.dirname(root))
      except : pass
      header=('<TABLE BORDER=2> <TR><TD> Parameters from wavecals <TD>Parameters from exposures (relative to wavecal)'
              '<TR><TD><IMG SRC=../plots/'+os.path.basename(root)+'_history.png>'+
              '<TD><IMG SRC=../plots/'+os.path.basename(root)+'_exposures.png></TABLE>')
      html.htmltab(grid,file=root+'.html',xtitle=xtit,ytitle=yt,header=header)
  else: 
      pdb.set_trace()
      for i in range(4) : plt.close()


def save_apWave(pars,out=None,group=0,rows=np.arange(300),npoly=4,frames=[],rms=None,sig=None,allpars=None) :
    """ Write the apWave files in standard format given the wavecal parameters
    """
    x = np.zeros([3,2048])
    allhdu=[]
    for ichip,chip in enumerate(chips) :
        hdu=fits.HDUList()
        x[0,:] = np.arange(2048)
        x[1,:] = ichip+1
        x[2,:] = group
        chippars=np.zeros([14,300])
        chipwaves=np.zeros([300,2048])
        for row in rows :
            pow=npoly-1-np.arange(npoly)
            polypars=pars[0:npoly,row]*3000**pow
            chippars[:,row] = np.append([ -1023.5+(ichip-1)*2048 + pars[npoly+ichip,row],0., 0., 1., 0., 0.], 
                              np.flip( np.append(np.zeros(8-npoly),polypars)))
            chipwaves[row,:] = func_multi_poly(x,*pars[:,row],npoly=npoly)
        hdu.append(fits.PrimaryHDU())
        hdu[0].header['NFRAMES']=(len(frames),'number of frames in fit')
        for i in range(len(frames)) : hdu[0].header['FRAME{:d}'.format(i)] = frames[i]
        hdu[0].header['NPOLY']=(npoly,'polynomial order of fit')
        if allpars is not None : 
            ngroup = int(round((allpars.shape[0]-npoly)/3))
            hdu[0].header['NGROUP']=(ngroup,'number of groups in fit')
        hdu[0].header['COMMENT']='HDU#1 : wavelength calibration parameters [14,300]'
        hdu[0].header['COMMENT']='HDU#2 : wavelength calibration array [300,2048]'
        hdu[0].header['COMMENT']='HDU#3 : wavecal fit parameter array [npoly+3*ngroup,300]'
        hdu[0].header['COMMENT']='HDU#4 : rms from fit [300,ngroup]'
        hdu[0].header['COMMENT']='HDU#5 : sig from fit [300,ngroup]'
        if allpars is not None : hdu[0].header['COMMENT']='HDU#3 : wavecal fit parameter array [npoly+3*ngroup,300]'
        hdu[0].header.add_comment('APOGEE_VER:'+os.environ['APOGEE_VER'])
        hdu.append(fits.ImageHDU(chippars))
        hdu.append(fits.ImageHDU(chipwaves))
        hdu.append(fits.ImageHDU(pars))
        if rms is not None :
            hdu[0].header['MEDRMS']=(np.nanmedian(rms),'median rms')
            hdu.append(fits.ImageHDU(rms))
        if sig is not None :
            hdu[0].header['MEDSIG']=(np.nanmedian(sig),'median sig')
            hdu.append(fits.ImageHDU(sig))
        if allpars is not None :
            hdu.append(fits.ImageHDU(allpars))
        if out is not None: hdu.writeto(out.replace('Wave','Wave-'+chip),overwrite=True)
        allhdu.append(hdu)
    return allhdu

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
                       ('dpixel','f4'), ('wave_found','f4'), ('frameid','i4')
                       ])
    nline=0
    for ichip,chip in enumerate(['a','b','c']) :
        # Use median offset of previous row for starting guess
        # Add a dummy first row to get starting guess offset for the first row
        dpixel_median = 0.
        for irow,row in enumerate(np.append([rows[0]],rows)) :
            # subtract off median-filtered spectrum to remove background
            medspec = frame[chip][1].data[row,:]-medfilt(frame[chip][1].data[row,:],101)
            j=np.where(lines['CHIPNUM'] == ichip+1)[0]
            dpixel=[]
            # for dummy row, open up the search window by a factor of two
            if irow == 0 : estsig0=2*estsig
            else : estsig0=estsig
            for iline in j :
                wave=lines['WAVE'][iline]
                try :
                    pix0=wave2pix(wave,waves[chip][row,:])+dpixel_median
                    # find peak in median-filtered subtracted spectrum
                    pars=peakfit(medspec,pix0,estsig=estsig0,
                                 sigma=frame[chip][2].data[row,:],mask=frame[chip][3].data[row,:])
                    if lines['USEWAVE'][iline] == 1 : dpixel.append(pars[1]-pix0)
                    if irow > 0 :
                        linestr['chip'][nline] = ichip+1
                        linestr['row'][nline] = row
                        linestr['wave'][nline] = wave
                        linestr['peak'][nline] = pars[0]
                        linestr['pixel'][nline] = pars[1]
                        linestr['dpixel'][nline] = pars[1]-pix0
                        linestr['wave_found'][nline] = pix2wave(pars[1],waves[chip][row,:])
                        linestr['frameid'][nline] = num
                        nline+=1
                    if out is not None :
                        out.write('{:5d}{:5d}{:12.3f}{:12.3f}{:12.3f}{:12.3f}{:12d}\n'.format(
                                  ichip+1,row,wave,pars[0],pars[1],pars[1]-pix0,num))
                    elif verbose :
                        print('{:5d}{:5d}{:12.3f}{:12.3f}{:12.3f}{:12.3f}{:12d}'.format(
                              ichip+1,row,wave,pars[0],pars[1],pars[1]-pix0,num))
                except :
                    if verbose : print('failed: ',num,row,chip,wave)
            if len(dpixel) > 10 : dpixel_median = np.median(np.array(dpixel))
            if verbose: print('median offset: ',row,chip,dpixel_median)

    return linestr[0:nline]

def gaussbin(x,a,x0,sig) :
    """ Evaluate integrated Gaussian function 
    """
    # bin width
    xbin=1.
    t1=(x-x0-xbin/2.)/np.sqrt(2.)/sig
    t2=(x-x0+xbin/2.)/np.sqrt(2.)/sig
    y=(myerf(t2)-myerf(t1))/xbin
    return a*y

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

def func_multi_poly(x,*pars, **kwargs) :
    """ Convert pixel to wavelength using wavecal parameters
          w = poly(x + offset(group,chip))
          pars = [npoly coefficients, ngroup*3 chip offsets]
        Args:  
            x (float) : [3,npts] array of (pixel,chip,group)
         pars (float) : input parameter array

        Returns :
         wave (float) : wavelength array for input pixel(s), parameters
    """
    wave=np.zeros(x.shape[1])
    ngroup = int(round(x[2,:].max()))+1
    nchip = 3
    npoly=kwargs.get('npoly',None)
    if npoly is None : npoly = len(pars)-ngroup*nchip
    coef = pars[0:npoly]
    # loop over all chip/group combinations
    for ichip in range(nchip) :
        for igroup in range(ngroup) :
            offset = pars[npoly+igroup*nchip+ichip]
            j=np.where((x[1,:] == ichip+1) & (np.round(x[2,:]).astype(int) == igroup))[0]
            xglobal = x[0,j] - 1023.5 + (ichip-1)*2048 + offset
            wave[j] = np.polyval(coef,xglobal)
    return wave

def getgroup(groups) :
    """ Given input list of group ids that may not be consecutive, return consecutive list
    """
    group=sorted(set(groups))
    out = np.zeros(len(groups))
    for i in range(len(group)) :
        j=np.where(groups == group[i])[0]
        out[j] = i
    return out,group


def skycal(planfile,out=None,inst=None,waveid=None,group=-1,skyfile='airglow',vers=None,nosky=False) :
    """ Determine positions of skylines for all frames in input planfile
    """
    # read planfile
    if type(planfile) is dict :
        p=planfile
        dirname='.'
    else :
        p=yanny.yanny(planfile)
        dirname=os.path.dirname(planfile)
    if dirname == '' : dirname = '.'
    if inst is None : inst = p['instrument'].strip("'") if p.get('instrument') else 'apogee-n'
    if vers is None : vers = p['apred_vers'].strip("'") if p.get('apred_vers') else 'current'
    if waveid is None : waveid = int(p['waveid'].strip("'")) if p.get('waveid') else None

    # open output line data?
    if out is not None : f=open(out,'a') 
    else : f=None

    # get the plugmap to get the sky fibers
    if p['platetype'].strip("'") == 'sky' : 
        skyrows = np.arange(300)
    elif p['platetype'].strip("'") == 'single' : 
        if p['telescope'].strip("'") == 'apo1m' : 
            if int(p['fixfiberid']) == 1 :
                fiberid=np.array([218,220,222,223,226,227,228,229,230,231])
            else :
                fiberid=np.array([218,219,221,223,226,228,230])
            skyfibers = fiberid[np.where(fiberid != int(p['APEXP']['single'][0]))[0]]
            skyrows = np.sort(300-skyfibers)
            print(skyrows)
    else :
        plugmjd=p['plugmap'].split('-')[1]
        if inst == 'apogee-s' : 
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
    if not nosky :
        skylines=ascii.read(os.environ['APOGEE_DIR']+'/data/skylines/'+skyfile+'.txt')

    # set up file reader
    load=apload.ApLoad(apred=vers,instrument=inst,verbose=False)

    # if we have a wavecal, get the wavelength array
    if waveid > 0 :
        #use wavelength solution from specified wavecal
        print('loading waveid: ', waveid)
        waveframe=load.apWave(waveid)
        npoly=waveframe['a'][0].header['NPOLY']
        allpars=waveframe['a'][3].data
        if len(waveframe['a']) == 6 : allpars,waves=refine(waveframe['a'][3].data)
        else : 
            waves={}
            for chip in chips : waves[chip]=waveframe[chip][2].data
    
    # loop over all frames in the planfile and assess skylines in each
    grid=[]
    ytit=[]
    for iframe,name in enumerate(p['APEXP']['name']) :
        print('frame: ', name)
        frame = load.ap1D(int(name))
        if waveid > 0 :
            if not nosky : plot = dirname+'/plots/skypixshift-'+name+'-'+skyfile
            if group >= 0 :
                allpars=waveframe['a'][3].data
                waves={}
                x = np.zeros([3,2048])
                for ichip,chip in enumerate(chips) :
                    x[0,:] = np.arange(2048)
                    x[1,:] = ichip+1
                    x[2,:] = group
                    waves[chip]=np.zeros([300,2048])
                    for row in np.arange(300) :
                        waves[chip][row,:] = func_multi_poly(x,*allpars[:,row],npoly=4)
        else :
            # use existing wavelength solution from ap1D file after ap1dwavecal has been run
            for chip in chips : waves[chip] = frame[chip][4].data
            if not nosky : plot = dirname+'/plots/skydeltapixshift-'+name+'-'+skyfile

        if nosky :
            w = np.zeros(4)
        else :

            # only use USEWAVE=1 lines for fit (can output others to derive wavelengths)
            gd = np.where(skylines['USEWAVE'] == 1)[0]
            nuse = 0
            fact = 1.
            niter = 0
            while (nuse < 0.9*len(gd)) & (niter < 0.75*len(skyrows)) :
                linestr = findlines(frame,skyrows,waves,skylines,out=f,estsig=1.*fact)
                use=[]
                nuse = 0
                for line in skylines['WAVE'][gd] :
                    j=np.where((linestr['wave'] == line) & (linestr['peak'] > 500.) )[0]
                    print(line,len(j),linestr['wave_found'][j].mean(),linestr['wave_found'][j].std())
                    use.extend(j)
                    # if we found this line in more than 5 fibers, count it
                    if len(j) > 5: nuse+=1
                use=np.array(use)
                print('fact : {:f} nuse: {:d}  ngd: {:d}'.format(fact,nuse,len(gd)))
                fact *= 1.15
                niter += 1

            # remove persistence affected fibers
            if int(p['mjd']) < 56860 and inst == 'apogee-n' :
                mask = np.ones(len(use), dtype=bool) # all elements included/True.
                bd = np.where((linestr['row'][use] > 200) & (linestr['chip'][use] == 3) )[0]
                mask[bd] = False              # Set unwanted elements to False
                use=use[mask]

            # solve for 4 parameter fit to dpixel, with linear trend with row, plus 2 chip offsets
            design=np.zeros([len(use),4])
            # global slope with rows
            design[:,0] = linestr['row'][use]
            # offset of each chip
            for ichip in range(3) :
                gd = np.where(linestr['chip'][use] == ichip+1)[0]
                design[gd,ichip+1] = 1.
            y=linestr['dpixel'][use]
            # reject outliers
            med=np.median(y)
            gd=np.where(np.abs(y-med) < 2.5)[0]
            gd=np.where(np.abs(y-med) < 0.5)[0]
            design=design[gd,:]
            y=y[gd]
            # if 1m, don't solve for a slope, not enough information
            if p['telescope'].strip("'") == 'apo1m' : design = design[:,1:4]
            # solve
            try : w = np.linalg.solve(np.dot(design.T,design), np.dot(design.T, y))
            except : 
                print('fit failed ....')
                pdb.set_trace()
            if p['telescope'].strip("'") == 'apo1m' : w=np.append([0.],w)
            print('fit parameters: ', w)

        if waveid > 0 :
            # adjust wavelength solution based on skyline fit
            newpars=copy.copy(allpars)
            for ichip in range(3) :
                newpars[npoly+ichip,:] -= (w[0]*np.arange(300) + w[ichip+1])
            allhdu = save_apWave(newpars,npoly=npoly)

            # rewrite out 1D file with adjusted wavelength information
            outname=load.filename('1D',num=int(name),mjd=load.cmjd(int(name)),chips=True)
            for ichip,chip in enumerate(chips) :
                hdu=fits.HDUList()
                frame[chip][0].header['HISTORY'] = 'Added wavelengths from SKYCAL, waveid: {:08d}'.format(waveid)
                frame[chip][0].header['HISTORY'] = 'Wavelength shift parameters {:12.5e} {:8.3f} {:8.3f} {:8.3f}'.format(w[0],w[1],w[2],w[3])
                frame[chip][0].header['WAVEFILE'] = load.allfile('Wave',num=waveid,chips=True)
                frame[chip][0].header['WAVEHDU'] = 5

                hdu.append(frame[chip][0])
                hdu.append(frame[chip][1])
                hdu.append(frame[chip][2])
                hdu.append(frame[chip][3])
                hdu.append(allhdu[ichip][2])
                hdu.append(allhdu[ichip][1])
                hdu.writeto(outname.replace('1D-','1D-'+chip+'-'),overwrite=True)

        # plots
        if plot is not None :
            try: os.mkdir(dirname+'/plots')
            except: pass
            # plot the pixel shift for each chip derived from the airglow lines
            fig,ax = plots.multi(1,1)
            wfig,wax = plots.multi(1,3)
            for ichip in range(3) :
                gd=np.where(linestr['chip'][use] == ichip+1)[0]
                med=np.median(linestr['dpixel'][use[gd]])
                x = linestr['row'][use[gd]]
                y = linestr['dpixel'][use[gd]]
                plots.plotp(ax,x,y,color=colors[ichip],xr=[0,300],yr=[med-0.2,med+0.2],
                            size=12,xt='Row',yt='Pixel shift')
                plots.plotc(wax[ichip],linestr['wave'][use[gd]],y,linestr['row'][use[gd]],zr=[0,300],yr=[med-0.5,med+0.5],
                            xr=xlim[ichip],size=12,xt='Wavelength',yt='Pixel shift')
                gdfit=np.where(np.abs(y-med) < 0.5)[0]
                xx=np.arange(300)
                #if len(gdfit) > 1 :
                #    p=np.polyfit(x[gdfit],y[gdfit],1)
                #    xx=np.arange(300)
                #    plots.plotl(ax,xx,p[0]*xx+p[1],color=colors[ichip])
                yy=w[0]*xx
                yy+=w[ichip+1]
                plots.plotl(ax,xx,yy,color=colors[ichip])
                if waveid > 0 : label = 'Frame: {:s}  Waveid: {:8d}'.format(name,waveid)
                else : label = 'Frame: {:s}  Delta from ap1dwavecal'.format(name)
                ax.text(0.1,0.9,label,transform=ax.transAxes)
            if type(plot) is str or type(plot) is unicode: 
                wfig.tight_layout()
                wfig.savefig(plot+'_wave.png')
                fig.savefig(plot+'.png')
                grid.append(['../plots/'+os.path.basename(plot)+'.png','../plots/'+os.path.basename(plot)+'_wave.png'])
                ytit.append(name)
            else : 
                plt.show()
                plt.draw()
                pdb.set_trace()
            plt.close('all')

        # Get shifts relative to first frame for each line/fiber
        if iframe == 0 : 
            linestr0 = copy.copy(linestr)
            use0 = copy.copy(use)
            refnum = int(name)
        print(iframe,len(linestr0),len(linestr))
        for line in range(len(use)) :
            ref = np.where((linestr0['chip'][use0] == linestr['chip'][use[line]]) & 
                           (linestr0['row'][use0] == linestr['row'][use[line]]) &
                           (linestr0['wave'][use0] == linestr['wave'][use[line]]))[0]
            if len(ref) > 0 : 
                linestr['pixel'][use[line]] -= linestr0['pixel'][use0[ref]].mean()
            else : linestr['pixel'][use[line]] = -999
        med = np.median(linestr['pixel'][use])

        # plot shifts relative to first frame, i.e. dithershift via sky lines
        if plot is not None:
            fig,ax=plots.multi(1,1)
            for ichip in range(3) :
                gd = np.where(linestr['chip'][use] == ichip+1)[0]   
                x = linestr['row'][use[gd]]
                y = linestr['pixel'][use[gd]]
                plots.plotp(ax,x,y, size=12,xr=[0,300],yr=[med-0.1,med+0.1],
                            xt='Row',yt='Pixel Shift',color=colors[ichip])
                gdfit=np.where(np.abs(y-med) < 0.5)[0]
                if len(gdfit) > 1 :
                    pfit=np.polyfit(x[gdfit],y[gdfit],1)
                    xx=np.arange(300)
                    plots.plotl(ax,xx,pfit[0]*xx+pfit[1],color=colors[ichip])
                label = 'Frame: {:8d}  Waveid: {:8d}'.format(int(name),refnum)
                ax.text(0.1,0.9,label,transform=ax.transAxes)
            fig.savefig(dirname+'/plots/skydithershift-'+name+'.png')
            plt.close()
    if plot is not None : 
        try: os.mkdir(dirname+'/html')
        except: pass
        html.htmltab(grid,file=dirname+'/html/skywavecal.html',ytitle=ytit)
    return linestr

def getskywave(frame,waveid,group=-1,vers='test',telescope='apo25m',plugmap=None) :
    """ Given input frame and waveid/group for frame taken without dither move to waveid,
        return skyline wavelengths
    """
    p={}
    p['APEXP']={}
    p['APEXP']['name']=[str(frame)]
    p['mjd'] = frame // 10000  + 55562
    p['waveid'] = str(waveid)
    if plugmap is None :
        p['platetype'] = 'sky'
    else :
        p['platetype'] = 'object'
        p['plugmap'] = plugmap
    p['telescope'] = telescope
    if telescope == 'lco25m' : inst = 'apogee-s'
    else : inst = 'apogee-n'
    p['apred_vers'] = vers
    p['instrument'] = inst
    return skycal(p,group=group)
    
def skywaves() :
    """ get skyline wavelengths from some particular frames, 16390029 and 22430033
    """
    linestr1 = getskywave(16930029,16680000,group=0,plugmap='8615-57255-02')
    pdb.set_trace()
    linestr2 = getskywave(22430033,20380000,group=18,plugmap='9050-57804-01')
    linestr=np.append(linestr1,linestr2)
    # derived wavelengths for sky lines (for adjusting airglow file initially)
    for line in set(linestr['wave']) :
        j=np.where(linestr['wave'] == line)[0]
        print(line,len(j),linestr['wave_found'][j].mean(),linestr['wave_found'][j].std())

def plotskywave(apred='r11',inst='apogee-n') :

    os.chdir(os.environ['APOGEE_REDUX']+'/'+apred+'/exposures/') 
    files = glob.glob(inst+'/*/a?1D-a*.fits')
    wfit=[]
    mjd=[]
    for file in files :
        print(file)
        a=fits.open(file)[0].header
        hist=a['HISTORY']
        for line in hist :
            if line.split()[0] == 'Wavelength' :
                wfit.append(line.split()[3:])
                mjd.append(a['JD-MID']-2400000.5)
    mjd=np.array(mjd)
    wfit=np.array(wfit).astype(float)
    fig,ax=plots.multi(1,4,hspace=0.001)
    plots.plotp(ax[0],mjd,wfit[:,0],yr=[-1.e-3,1.e-3],yt='slope')
    plots.plotp(ax[1],mjd,wfit[:,2],yr=[-5,5],yt='g')
    plots.plotp(ax[2],mjd,wfit[:,1]-wfit[:,2],yr=[-0.2,0.2],yt='r-g')
    plots.plotp(ax[3],mjd,wfit[:,3]-wfit[:,2],yr=[-0.2,0.2],yt='b-g')
    fig.savefig(inst+'.png')
        
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

def compare(npoly=4,lco=False) :

    if lco :
        files=glob.glob('asWave-b-*.fits')
        out='apogee-s'
        root='lco'
    else :
        files=glob.glob('apWave-b-*0000.fits')
        out='apogee-n'
        root='apo'
    files.sort()
    dates=[]
    for file in files :
        dates.append(int(file.split('-')[2].replace('.fits',''))/10000)
    dates=np.array(dates)
    files=np.array(files)

    # wavelengths to compare solutions at
    w=np.arange(15160.,16900.)
    rows=np.arange(300.)
    x=np.arange(-1024-2048-150,1024+2048+150,25)
    x=np.arange(-1024-2048-150,1024+2048+150)
    grid=[]
    ytit=[]
    for year in range(-1,7) :
      if year == -1 :
          i1=55757-55562
          i2=99999
      else :
          i1 = 55757+year*365-55562
          i2 = i1+365
      j = np.where((dates >=i1) & (dates<=i2) & ((dates<2430) | (dates>2450)) )[0]
      print('year: ', year,i1,i2,len(j))
      maxgroup=20
      if len(j) > 0 :
        wave=np.zeros([len(x),len(j)])
        # in pix, store the global pixel corresponding to the range of wavelengths
        # this is better than looking at the wavelength comparison of different solutions, 
        # because if the chips have moved (shifted dither position), this is a constant
        # global pixel offset, but not a constant wavelength offset (because dispersion varies)
        pix=np.zeros([len(w[::50]),len(j)])
        pixraw=np.zeros([len(w[::50]),len(j)])
        chipa=np.zeros([300,len(j)*maxgroup])
        chipc=np.zeros([300,len(j)*maxgroup])
        chipafit=np.zeros([300,len(j)*maxgroup])
        chipcfit=np.zeros([300,len(j)*maxgroup])
        fig,ax=plots.multi(1,3)
        nfile=0
        noffset=0
        gdfiles=[]
        for ifile,file in enumerate(files[j]) :
            print(file)
            try :
              a=fits.open(file)[3].data
              wave[:,ifile]=np.polyval(a[0:npoly,150],x)
              pix[:,ifile]=wave2pix(w,wave[:,ifile])[::50]
              ngroup=fits.open(file)[0].header['NGROUP']
              for igroup in range(ngroup) :
                  chipa[:,noffset]=a[npoly+igroup*3]
                  # fit a lit to the chip offsets as a function of row, ignoring bad fits
                  gd = np.where(np.abs(a[npoly+igroup*3]) > 1)[0]
                  p=np.polyfit(rows[gd],a[npoly+igroup*3,gd],1)
                  chipafit[:,noffset]=p[0]*rows+p[1]
                  chipc[:,noffset]=a[npoly+2+igroup*3]
                  p=np.polyfit(rows[gd],a[npoly+2+igroup*3,gd],1)
                  chipcfit[:,noffset]=p[0]*rows+p[1]
                  noffset+=1
              nfile +=1
              gdfiles.append(file.split('.')[0].split('-')[2])
            except:
              pass
        # exclude bad/missing columns
        chipa=chipa[:,0:noffset]
        chipc=chipc[:,0:noffset]
        chipafit=chipafit[:,0:noffset]
        chipcfit=chipcfit[:,0:noffset]
        wave=wave[:,0:nfile]
        pix=pix[:,0:nfile]
        pixraw=pixraw[:,0:nfile]

        wmed = np.median(wave,axis=1)
        pmed = np.nanmedian(pix,axis=1)
        chipamed = np.median(chipafit,axis=1)
        chipcmed = np.median(chipcfit,axis=1)
        for ifile in range(nfile) :
            pix[:,ifile]-=pmed
            wave[:,ifile]-=wmed
            wave[:,ifile]-=np.median(wave[:,ifile])
            pixraw[:,ifile]=pix[:,ifile]
            pix[:,ifile]-=np.median(pix[:,ifile])
            plots.plotl(ax[0],x,wave[:,ifile],yr=[-0.5,0.5],xt='global pixel',yt=r'$\lambda-\lambda_{med}$')
        for ifile in range(noffset) :
            # to account for dither shifts, subtract median pixel for this solution
            # for chip gaps, get shift relative to median across all solutions
            chipa[:,ifile]-=chipamed
            chipc[:,ifile]-=chipcmed
            chipafit[:,ifile]-=chipamed
            chipcfit[:,ifile]-=chipcmed
            plots.plotl(ax[1],np.arange(300),chipa[:,ifile],yr=[-0.5,0.5],xt='Fiber',yt='chipa-chipamed')
            plots.plotl(ax[2],np.arange(300),chipc[:,ifile],yr=[-0.5,0.5],xt='Fiber',yt='chipc-chipcmed')
        if year < 0 : tit='All years'
        else : tit='Year: {:d}'.format(year)
        ytit.append(tit)
        name=root+'year{:1d}'.format(year)
        fig.savefig(name+'.png'.format(year))
        plt.close()
        t=tv.TV(aspect='auto')
        t.cmap='viridis'
        row=[name+'.png']
        names=['pixraw','pix','chipa','chipc','chipafit','chipcfit']
        for i,im in enumerate([pixraw,pix,chipa,chipc,chipafit,chipcfit]) :
            t.ax.cla()
            t.tv(im,min=-0.05,max=0.05)
            if i < 2 :
                for ifile,f in enumerate(gdfiles) :
                    t.ax.text(ifile+0.5,-1.,str(f),ha='right',rotation=90,fontsize=8)
            t.fig.suptitle(tit)
            t.fig.savefig(root+name+names[i]+'.png'.format(year))
            row.append(root+name+names[i]+'.png')
        plt.close()
        grid.append(row)
    html.htmltab(grid,file=out+'.html',ytitle=ytit)

def allplots() :
    """ Routine to put together master web page for summary plots from all years
    """ 

    # summary plots should be made in wavecal!
    fig,ax=plots.multi(1,4,hspace=0.001)
    cb_ax=fig.add_axes((0.9,0.72,0.03,0.15))
    cb_ax2=fig.add_axes((0.9,0.15,0.03,0.4))
    grid=[]
    ytit=[]
    for ical,cal in enumerate([2380000,5680000,9500000,13140000,16680000,20380000,24040000,22670000,24040000]) :
        if ical<7 :
            root='apWave-{:08d}'.format(cal)
        else :
            root='asWave-{:08d}'.format(cal)

        # read the fit parameters
        a = fits.open(root.replace('-','-b-')+'.fits'.format(cal))[3].data
        # get chip positions relative to median postion across all groups
        chipa=a[4:200:3,:]-np.median(a[4:200:3,:],axis=0)
        chipc=a[6:200:3,:]-np.median(a[6:200:3,:],axis=0)
        chipb=a[5:200:3,:]-np.median(a[5:200:3,:],axis=0)
        # image of chip b shifts
        aximage=ax[0].imshow(chipb,vmin=-2,vmax=2,cmap='viridis',interpolation='nearest',aspect='auto')
        ax[0].set_ylabel('chip loc')
        fig.colorbar(aximage,cax=cb_ax,orientation='vertical')

        # get chip b shift relative to median across all rows
        chipb=(chipb.T-np.median(chipb,axis=1)).T
        ax[1].imshow(chipb,vmin=-0.03,vmax=0.03,cmap='viridis',interpolation='nearest',aspect='auto')
        ax[1].set_ylabel('rel chip loc')
        # chip gaps 
        ax[2].imshow(chipa,vmin=-0.03,vmax=0.03,cmap='viridis',interpolation='nearest',aspect='auto')
        ax[2].set_ylabel('g-r gap')
        aximage=ax[3].imshow(chipc,vmin=-0.03,vmax=0.03,cmap='viridis',interpolation='nearest',aspect='auto')
        ax[3].set_xlabel('Row')
        ax[3].set_ylabel('b-g gap')
        fig.suptitle('{:08d}'.format(cal))
        fig.colorbar(aximage,cax=cb_ax2,orientation='vertical')
        fig.savefig('plots/'+root+'_sum.png'.format(cal))
        grid.append([root+'.png',root+'_chiploc.png',root+'_sum.png'])
        ytit.append(root)
        pdb.set_trace()
    html.htmltab(grid,file='plots/all.html',ytitle=ytit)

def gauss(x,a,x0,sig) :
    """ Evaluate Gaussian function 
    """
    return a/np.sqrt(2*np.pi)/sig*np.exp(-(x-x0)**2/2./sig**2)

def myerf(t) :
    """ Evaluate function that integrates Gaussian from -inf to t
    """
    neg = np.where(t<0.)[0]
    pos = np.where(t>=0.)[0]
    out = t*0.
    out[neg] = erfc(abs(t[neg]))/2.
    out[pos] = 0.5+erf(abs(t[pos]))/2.
    return out


def ditherplots(planfile,vers=None,inst=None) :
    """ Make HTML file for dither plots for a visit
    """
    p=yanny.yanny(planfile)
    dirname=os.path.dirname(planfile)

    if dirname == '' : dirname = '.'
    if inst is None : inst = p['instrument'].strip("'") if p.get('instrument') else 'apogee-n'
    if vers is None : vers = p['apred_vers'].strip("'") if p.get('apred_vers') else 'current'

    # set up file reader
    load=apload.ApLoad(apred=vers,instrument=inst,verbose=False)

    grid = []
    for frame in p['APEXP']['name'] :
        grid.append(['../plots/skypixshift-'+frame+'-airglow.png',
                     '../plots/skypixshift-'+frame+'-airglow_wave.png',
                     '../plots/skydithershift-'+frame+'.png',
                     '../plots/dithershift-'+frame+'.gif'])
    html.htmltab(grid,file=dirname+'/html/shift.html')

def shape(w,title=None,out=None,figax=None) :
   
    if figax is None : fig,ax=plots.multi(1,3,hspace=0.001)
    else : fig,ax = figax
    x=np.arange(300)
    for ichip,chip in enumerate(chips) :
        for col in [50,55,1995,2000] :
            y=w[chip][2].data[:,col]-w[chip][2].data[:,1000]
            plots.plotl(ax[ichip],x,y-np.median(y),yr=[-0.2,0.2])
            plots.plotl(ax[ichip],x,y-np.median(y),yr=[-0.2,0.2])
    if title is not None: fig.suptitle(title)
    if out is not None :
        fig.savefig(out)
        plt.close()
    return fig,ax

def allshape() :
    grid=[]
    r11 = apload.ApLoad(apred='r11')
    for waveid in [2380000,5680000,9500000,13140000,16680000,20380000,24040000] :
        w=r11.apWave(waveid)
        out='plots/shape_{:08d}.png'.format(waveid)
        shape(w,title='{:08d}'.format(waveid),out=out)
        grid.append(['../'+out])
    r11.settelescope('lco25m')
    for waveid in [22670000,24040000] :
        w=r11.apWave(waveid)
        out='plots/shape_{:08d}.png'.format(waveid)
        shape(w,title='{:08d}'.format(waveid),out=out)
        grid.append(['../'+out])
    r8 = apload.ApLoad(apred='r8')
    for waveid in [2420038] :
        w=r8.apWave(waveid)
        out='plots/shape_{:08d}.png'.format(waveid)
        shape(w,title='{:08d}'.format(waveid),out=out)
        grid.append(['../'+out])
    html.htmltab(grid,file='html/shape.html')

def refine(oldpars,npoly=4) :
    ''' Refine wavelength solution by averaging over groups,and smoothing over rows
    '''
    # copy parameters so as not to replace
    allpars=copy.copy(oldpars)

    # average the offsets for all of the groups. To do so, we first need to remove
    # the trends from the dither shifts, i.e. with a four-parameter fit of slope and chip offsets
    # do this relative to first group
    nframes=(allpars.shape[0]-npoly)//3
    for iframe in range(1,nframes) :
        design=np.zeros([900,4])
        y=np.zeros(900)
        # offset of each chip relative to first frame
        for ichip in range(3) :
            # global slope with rows
            design[ichip*300+np.arange(300),0] = np.arange(300)
            design[ichip*300+np.arange(300),ichip+1] = 1.
            gd = np.where(allpars[npoly+iframe*3+ichip,:] != 0.)[0]
            y[ichip*300+np.arange(300)[gd]] = allpars[npoly+iframe*3+ichip,gd]-allpars[npoly+ichip,gd]
        # reject bad fibers
        gd=np.where((abs(y) > 1.e-5) & (abs(y) < 200.))[0]
        design=design[gd,:]
        y=y[gd]
        # solve and replace offsets with offsets adjusted to first group dither position
        try : 
            w = np.linalg.solve(np.dot(design.T,design), np.dot(design.T, y))
            for ichip in range(3) : 
                gd = np.where(allpars[npoly+iframe*3+ichip,:] != 0.)[0]
                allpars[npoly+iframe*3+ichip,gd] -= (w[0]*np.arange(300)[gd] + w[ichip+1])
        except : 
            print('fit failed ....frame:',iframe)
            pdb.set_trace()
    # replace chip offsets of first group with average chip offsets
    # then calculate wavelength array for all chips and rows with this fit
    # also calculate smoothed wavelength array (across rows at each column). We will use this for a new fit
    waves={}
    swaves={}
    newwaves={}
    x = np.zeros([3,2048])
    for ichip,chip in enumerate(chips) : 
        # loop over rows so we can do outlier-rejected mean
        x[0,:] = np.arange(2048)
        x[1,:] = ichip+1
        x[2,:] = 0
        waves[chip]=np.zeros([300,2048])
        swaves[chip]=np.zeros([300,2048])
        newwaves[chip]=np.zeros([300,2048])
        for row in range(300) :
            y=allpars[4+ichip::3,row]
            # skip missing groups
            y=y[np.where(y != 0.)[0]]
            # reject outliers
            gd=np.where(abs(y-np.median(y)) < 0.1)[0]
            if len(gd) > 0 : allpars[4+ichip,row] = np.mean(y[gd])
            waves[chip][row,:] = func_multi_poly(x,*allpars[:,row],npoly=npoly)
        # to reduce noise further, fit relative wavelengths across rows and use the fit values for smoothed array
        rows=np.arange(300)
        for col in range(2048) :
            try :
                # don't include bad rows! i.e. from missing fibers
                gd = np.where( np.isfinite(waves[chip][:,col]-waves[chip][:,1024]) & (waves[chip][:,col] > 0.) )[0]
                pfit = np.polyfit(rows[gd],waves[chip][gd,col]-waves[chip][gd,1024],3)
                swaves[chip][gd,col] = np.polyval(pfit,rows[gd]) + waves[chip][gd,1024]
            except :
                print('fit across rows failed, col: ',col)
                pdb.set_trace()

    # now refit the full wavelength solutions using swaves as input
    newpars=[]
    for row in range(300) :
        if allpars[4,row] == 0. : 
            # if this was a bad row before, keep it bad
            popt = pars*0.
            newpars.append(popt)
            continue

        x = np.zeros([3,2048*3])
        y = np.zeros([2048*3])
        for ichip,chip in enumerate(chips) : 
            x[0,ichip*2048:(ichip+1)*2048] = np.arange(2048)
            x[1,ichip*2048:(ichip+1)*2048] = ichip+1
            x[2,ichip*2048:(ichip+1)*2048] = 0
            y[ichip*2048:(ichip+1)*2048] = swaves[chip][row,:]
        # we will fix the central chip position at 0, and allow wavelength to float
        pars = allpars[0:npoly+3,row]
        pars[npoly+1] = 0.
        bounds = ( np.zeros(len(pars))-np.inf, np.zeros(len(pars))+np.inf)
        bounds[0][npoly+1] = -1.e-7
        bounds[1][npoly+1] = 1.e-7
        try :
            popt,pcov = curve_fit(func_multi_poly,x,y,p0=pars,bounds=bounds)
        except :
            print('Solution failed for row: ', row)
            pdb.set_trace()
            popt = pars*0.
        newpars.append(popt)
        # calculate wavelength arrays from refined solution
        x = np.zeros([3,2048])
        for ichip,chip in enumerate(chips) : 
            x[0,:] = np.arange(2048)
            x[1,:] = ichip+1
            x[2,:] = 0
            newwaves[chip][row,:] = func_multi_poly(x,*popt,npoly=npoly)
    # return new parameters and wavelength array
    return np.array(newpars).T,newwaves
    
