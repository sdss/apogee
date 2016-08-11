from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
from holtz.tools import plots
import numpy as np
import pdb

def fit1d(xdata,zdata,degree=1,reject=0,ydata=None,plot=None,plot2d=False,xr=None,yr=None,zr=None,xt=None,yt=None,zt=None) :
    """ 
    Do a 1D polynomial fit to data set and plot if requested

    Args:
        xdata : independent variable
	zdata : dependent variable to be fit 

    Keyword args:
        degree: degree of polynomial to fit (default=1 for linear fit)
        reject : single iteration rejection of points that deviate from initial by more than specified value (default=0, no rejection)
	ydata : auxiliary variable for plots (default=None)
        plot  : axes to plot into (default=None)
        plot2d (bool)  : set to make a 2D plot with auxiliary variable, rather than 1D color-coded by auxiliary variable
        xr[2] : xrange for plot
        yr[2] : yrange for plot
        zr[2] : zrange for plot
        xt    : xtitle for plot
        yt    : ytitle for plot
        zt    : ztitle for plot

    Returns :
        pfit  : 1D polynomial fit
    """

    # set up fitter and do fit
    fit_p = fitting.LinearLSQFitter()
    p_init = models.Polynomial1D(degree=degree)
    pfit = fit_p(p_init, xdata, zdata)
    # rejection of points?
    if reject > 0 :
        gd=np.where(abs(zdata-pfit(xdata)) < reject)[0]
        bd=np.where(abs(zdata-pfit(xdata)) >= reject)[0]
        #for i in bd :
        #    print x[i],10.**data[i]
        print 'rejected ',len(xdata)-len(gd),' of ',len(xdata),' points'
        pfit = fit_p(p_init, xdata[gd], zdata[gd])

    # plot if requested
    if plot is not None :
        if xr is None : xr = [xdata.min(),xdata.max()]
        if yr is None and ydata is not None : yr = [ydata.min(),ydata.max()]
        if zr is None : zr = [zdata.min(),zdata.max()]
        if ydata is None :
            # straight 1D plot
            plots.plotp(plot,xdata,zdata,xr=xr,yr=yr,zr=zr,
                   xt=xt,yt=yt,size=15)
            x = np.linspace(xr[0],xr[1],200)
            plots.plotl(plot,x,pfit(x))
        elif plot2d :
            # 2D image plot with auxiliary variable
            y, x = np.mgrid[yr[1]:yr[0]:200j, xr[1]:xr[0]:200j]
            plots.plotc(plot,xdata,ydata,zdata,xr=xr,yr=yr,zr=zr,
                   xt=xt,yt=xt,zt=yt,colorbar=True,size=15)
            plot.imshow(pfit(x),extent=[xr[1],xr[0],yr[1],yr[0]],
                aspect='auto',vmin=zr[0],vmax=zr[1], origin='lower')
        else :
            # 1D plot color-coded by auxiliary variable
            plots.plotc(plot,xdata,zdata,ydata,xr=xr,yr=zr,zr=yr,
                   xt=xt,yt=yt,size=15)
            x = np.linspace(xr[0],xr[1],200)
            plots.plotl(plot,x,pfit(x),color='k')
        pdb.set_trace()
    return pfit

def fit2d(xdata,ydata,zdata,degree=1,plot=None,xr=None,yr=None,zr=None,xt=None,yt=None,zt=None,gdrange=None) :
    """ 
    Do a 2D polynomial fit to data set and plot if requested

    Args:
        xdata : independent variable
        ydata : independent variable
	zdata : dependent variable to be fit 

    Keyword args:
        degree: degree of polynomial to fit (default=1 for linear fit)
	ydata : auxiliary variable for plots (default=None)
        plot  : axes to plot into (default=None)
        xr[2] : xrange for plot
        yr[2] : yrange for plot
        zr[2] : zrange for plot
        xt    : xtitle for plot
        yt    : ytitle for plot
        zt    : ztitle for plot

    Returns :
        pfit  : 2D polynomial fit
    """

    if gdrange is not None :
        gd = np.where((zdata > gdrange[0]) & (zdata < gdrange[1]))[0]
        xfit = xdata[gd]
        yfit = ydata[gd]
        zfit = zdata[gd]
    else :
        xfit = xdata
        yfit = ydata
        zfit = zdata

    # set up fitter and do fit
    fit_p = fitting.LinearLSQFitter()
    p_init = models.Polynomial2D(degree=degree)
    pfit = fit_p(p_init, xfit, yfit, zfit)
    
    if plot is not None :
        if xr is None : xr = [xfit.min(),xfit.max()]
        if yr is None : yr = [yfit.min(),yfit.max()]
        if zr is None : zr = [zfit.min(),zfit.max()]
        # plot data
        plots.plotc(plot,xfit,yfit,zfit,xr=xr,yr=yr,zr=zr,
                xt=xt,yt=xt,zt=zt,colorbar=True,size=15,linewidth=1)
        # create independent variable grid for model and display
        y, x = np.mgrid[yr[1]:yr[0]:200j, xr[1]:xr[0]:200j]
        plot.imshow(pfit(x,y),extent=[xr[1],xr[0],yr[1],yr[0]],
                aspect='auto',vmin=zr[0],vmax=zr[1], origin='lower')
        plt.show()

    return pfit

def linear(data,design,err=None) :
    '''
    Given data array (npts) and design matrix (npar, npts), return linear solution
    '''
    npar=design.shape[0]
    npts=design.shape[1]
    if err is None:
        err = np.ones([npts])

    ATA = np.dot(design, design.T / err[:, None]**2)
    soln = np.linalg.solve(ATA, np.dot(design, data / err**2))
    inv = np.linalg.inv(ATA)

    #alpha=np.zeros([npar,npar])
    #beta=np.zeros([npar])
    #for ipar in range(npar) :
    #    beta[ipar]+=(data*design[ipar,:]/err**2).sum()
    #    for jpar in range(npar) :
    #        alpha[ipar,jpar]+=(design[ipar,:]*design[jpar,:]/err**2).sum()
    #c=np.linalg.inv(alpha)
    #soln=np.dot(c,beta)
    #pdb.set_trace()

    return soln, inv
 
