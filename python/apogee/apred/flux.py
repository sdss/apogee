# routines for flux calibration a plate based on tellurics

import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
import pdb
from apogee.utils import apload

chips=['a','b','c']

def getresp(plate,mjd,apred='r12',telescope='apo25m',plot=False) :
    """ Solve for response function and apply it
    """

    # get apPlate file
    load=apload.ApLoad(apred=apred,telescope=telescope) 
    apPlate=load.apPlate(plate,mjd)

    # get rows of tellurics
    fibers=apPlate['b'][11].data['FIBERID']
    tel=np.where(apPlate['b'][11].data['OBJTYPE'] == 'HOT_STD')[0]
    rows=300-fibers[tel]

    # do polynomial fit to log(flux), with 4th order plus offset fo each star,
    # using every 10th pixel in each chip, so we have 190 pixels * 3 chips * ntelluric data points
    # and 4 + ntellurics parameters
    npix=190
    nstars=len(rows)
    design=np.zeros([3*npix*nstars,4+nstars])
    y=np.zeros([3*npix*nstars])
    for ichip,chip in enumerate(chips) :
      for irow,row in enumerate(rows):
        x=apPlate[chip][4].data - 16000.
        design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix,0] = x[row,100:2000:10]**4
        design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix,1] = x[row,100:2000:10]**3
        design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix,2] = x[row,100:2000:10]**2
        design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix,3] = x[row,100:2000:10]
        design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix,4+irow] = 1.
        y[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix] = np.log10(apPlate[chip][1].data[row,100:2000:10])
    gd=np.where(np.isfinite(y))[0]
    design=design[gd,:]
    y=y[gd]
    # do the fit
    coef = np.linalg.solve(np.dot(design.T,design), np.dot(design.T, y))

    # apply the fit. Note that norm adds a term so that response gives 1/lambda**-2 shape
    for chip in chips :
        for row in np.arange(300) :
            print(chip,row)
            w=apPlate[chip][4].data[row,:]
            spec=apPlate[chip][1].data[row,:]
            resp=norm(w,coef)
            if plot :
                plt.plot(w,spec)
                plt.plot(w,resp*1000)
                plt.plot(w,spec/resp)
                plt.show()
            apPlate[chip][1].data[row,:] /= resp
        if plot :
            pdb.set_trace()
            plt.clf()
        file=load.filename('Plate',plate=plate,mjd=mjd,apred=apred,chips=True).replace('Plate-','Plate-'+chip+'-')
        apPlate[chip].writeto(file,overwrite=True)
    # now do the apVisit files
    for row in np.arange(300) :
        try:
            print('Row: ', row)
            apVisit=load.apVisit(plate,mjd,300-row)
            for ichip in range(3) :
                print(ichip)
                w=apVisit[4].data[ichip,:]
                resp=norm(w,coef)
                apVisit[1].data[ichip,:] /= resp
            file=load.filename('Visit',plate=plate,mjd=mjd,apred=apred,fiber=300-row)
            print(file)
            apVisit.writeto(file,overwrite=True)
            print('done')
        except : pass

def norm(w,coef) :
    """ Evaluate the polynomial fit
    """
    x = w-16000.
    logflux = coef[0]*x**4 + coef[1]*x**3 + coef[2]*x**2 + coef[3]*x
    logflux += 2*np.log10(w/16000.)
    return 10.**logflux

def main(args) :

    parser = argparse.ArgumentParser(
        description='Makes a SLURM batch file')

    parser = argparse.ArgumentParser(description="flux calibration apPlate",
                                     prog=os.path.basename(args[0]),
                                     usage="apflux plate mjd --apred apred --instrument inst-s")
    parser.add_argument("plate", type=int, help="plate")
    parser.add_argument("mjd", type=int, help="MJD")
    parser.add_argument("--apred", type=str, help="apred version", default='r12')
    parser.add_argument("-t", "--telescope", type=str, required=True, help="telescope")
    args = parser.parse_args()
    getresp(args.plate,args.mjd,apred=args.apred,telescope=args.telescope)


