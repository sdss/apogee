import matplotlib
import matplotlib.pyplot as plt
import os
import pdb
import pickle
import copy
import scipy.signal
import scipy.interpolate
import numpy as np
from astropy.modeling import models, fitting
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.io import ascii, fits
from astropy.convolution import convolve, Box1DKernel, Box2DKernel
import pyvista
from pyvista import image
from pyvista import tv
from tools import plots

ROOT = os.path.dirname(os.path.abspath(__file__)) + '/../../'


class SpecData(CCDData) :
    """ Class to include a wavelength array on top of CCDData, with simple read/write/plot methods
    """
    def __init__(self,data,wave=None) :
        if type(data) is str :
            hdulist=fits.open(data)
            self.meta = hdulist[0].header
            self.unit = hdulist[0].header['BUNIT']
            self.data = hdulist[1].data
            self.uncertainty = StdDevUncertainty(hdulist[2].data)
            self.mask = hdulist[3].data
            self.wave = hdulist[4].data
        elif type(data) is CCDData :
            self.unit = data.unit
            self.meta = data.meta
            self.data = data.data
            self.uncertainty = data.uncertainty
            self.mask = data.mask
            self.wave = wave
        else :
            print('Input must be a filename or CCDData object')

    def write(self,file,overwrite=True) :
        hdulist=fits.HDUList()
        hdulist.append(fits.PrimaryHDU(header=self.meta))
        hdulist.append(fits.ImageHDU(self.data))
        hdulist.append(fits.ImageHDU(self.uncertainty.array))
        hdulist.append(fits.ImageHDU(self.mask.astype(np.int16)))
        hdulist.append(fits.ImageHDU(self.wave))
        hdulist.writeto(file,overwrite=overwrite)

    def plot(self,ax,**kwargs) :
        for row in range(self.wave.shape[0]) :
            gd = np.where(self.mask[row,:] == False)[0]
            plots.plotl(ax,self.wave[row,gd],self.data[row,gd],**kwargs)
        


def get_wavecal(file) :
    """ load a wavecal object from disk file 
    """
    with open(file,'rb') as wavecal :
        return pickle.load(wavecal) 

class WaveCal() :
    """ Class for wavelength solutions
    """
    def __init__ (self,type='chebyshev',degree=2,ydegree=2,pix0=0,orders=[1]) :
        """ Initialize the wavecal object

            type : type of solution ('poly' or 'chebyshev')
            degree : polynomial degree for wavelength
            ydegree : polynomial degree for  y dimension
            pix0 : reference pixel
            orders : spectral order for each row
            spectrum : spectrum from which fit is derived
        """
        self.type = type
        self.degree = degree
        self.ydegree = ydegree
        self.pix0 = pix0
        self.orders = orders
        self.waves = None
        self.x = None
        self.y = None
        self.weights = None
        self.model = None
        self.ax = None

    def wave(self,pixels=None,image=None) :
        """ Wavelength from pixel using wavelength solution model

            pix : input pixel positions [x] or [y,x]
            image : for input image size [nrows,ncols], return wavelengths at all pixels
            returns wavelength
        """
        if pixels is not None :
            out=np.zeros(len(pixels[0]))
            for i,pixel in enumerate(pixels[0]) :
                if self.type.find('2D') > 0 :
                    order=self.orders[pixels[1][i]]
                    out[i]=self.model(pixel-self.pix0,pixels[1][i])/order
                else :
                    out[i]=self.model(pixel-self.pix0)/self.orders[0]
            return out
        else :
            out=np.zeros(image)
            cols=np.arange(out.shape[-1])
            if out.ndim == 2 :
                for row in range(out.shape[0]) : 
                    rows=np.zeros(len(cols))+row
                    try : order = self.orders[row]
                    except : order=self.orders[0]
                    out[row,:] = self.model(cols-self.pix0,rows)/order
            else :
                out= self.model(cols-self.pix0)/self.orders[0]
            return out

    def getmod(self) :
        """ Return model for current attributes
        """

        if self.type == 'poly' :
            mod=models.Polynomial1D(degree=self.degree)
        elif self.type == 'chebyshev' :
            mod=models.Chebyshev1D(degree=self.degree)
        elif self.type == 'chebyshev2D' :
            sz=self.spectrum.data.shape
            mod=models.Chebyshev2D(x_degree=self.degree,y_degree=self.ydegree,
                                   x_domain=[0,sz[1]],y_domain=[0,sz[0]])
        else :
            raise ValueError('unknown fitting type: '+self.type)
            return
        return mod

    def fit(self,plot=True) :
        """ do a wavelength fit 
        """
        print("doing wavelength fit")
        # set up fitter and model
        twod='2D' in self.type
        fitter=fitting.LinearLSQFitter()
        mod = self.getmod()

        if not hasattr(self,'ax') : self.ax = None
        if twod :
            nold=-1
            nbd=0
            while nbd != nold :
                nold=nbd
                self.model=fitter(mod,self.pix-self.pix0,self.y,self.waves*self.waves_order,weights=self.weights)
                diff=self.waves-self.wave(pixels=[self.pix,self.y])
                gd = np.where(self.weights > 0)[0]
                print('  rms: {:8.3f}'.format(diff[gd].std()))
                bd = np.where(abs(diff) > 3*diff.std())[0]
                nbd = len(bd)
                print('rejecting {:d} points from {:d} total: '.format(nbd,len(self.waves)))
                self.weights[bd] = 0.

            if self.ax is not None : 
                self.ax[1].cla()
                scat=self.ax[1].scatter(self.waves,diff,marker='o',c=self.y,s=2)
                scat=self.ax[1].scatter(self.waves[bd],diff[bd],marker='o',c='r',s=2)
                xlim=self.ax[1].get_xlim()
                self.ax[1].set_ylim(diff.min()-0.5,diff.max()+0.5)
                self.ax[1].plot(xlim,[0,0],linestyle=':')
                self.ax[1].text(0.1,0.9,'rms: {:8.3f}'.format(diff[gd].std()),transform=self.ax[1].transAxes)
                cb_ax = self.fig.add_axes([0.94,0.05,0.02,0.4])
                cb = self.fig.colorbar(scat,cax=cb_ax)
                cb.ax.set_ylabel('Row')
                plt.draw()
                self.fig.canvas.draw_idle()
                input('  See 2D wavecal fit. Hit any key to continue....')

        else :
            self.model=fitter(mod,self.pix-self.pix0,self.waves*self.waves_order,weights=self.weights)
            diff=self.waves-self.wave(pixels=[self.pix])
            print('  rms: {:8.3f} Angstroms'.format(diff.std()))
            if self.ax is not None :
                # iterate allowing for interactive removal of points
                done = False
                ymax = self.ax[0].get_ylim()[1]
                while not done :

                    # do fit
                    gd=np.where(self.weights>0.)[0]
                    bd=np.where(self.weights<=0.)[0]
                    self.model=fitter(mod,self.pix[gd]-self.pix0,self.waves[gd]*self.waves_order[gd],weights=self.weights[gd])
                    diff=self.waves-self.wave(pixels=[self.pix])
                    print('  rms: {:8.3f} Anstroms'.format(diff[gd].std()))

                    # replot spectrum with new fit wavelength scale
                    self.ax[0].cla()
                    self.ax[0].plot(self.wave(image=self.spectrum.data.shape)[0,:],self.spectrum.data[0,:])
                    # plot residuals
                    self.ax[1].cla()
                    self.ax[1].plot(self.waves[gd],diff[gd],'go')
                    self.ax[1].text(0.1,0.9,'rms: {:8.3f} Angstroms'.format(diff[gd].std()),transform=self.ax[1].transAxes)
                    self.ax[1].set_xlabel('Wavelength')
                    self.ax[1].set_ylabel('obs wave - fit wave')
                    if len(bd) > 0 : self.ax[1].plot(self.waves[bd],diff[bd],'ro')
                    self.ax[1].set_ylim(diff[gd].min()-0.5,diff[gd].max()+0.5)
                    for i in range(len(self.pix)) :
                        self.ax[1].text(self.waves[i],diff[i],'{:2d}'.format(i),va='top',ha='center')
                        if self.weights[i] > 0 :
                            self.ax[0].plot([self.waves[i],self.waves[i]],[0,ymax],'g')
                        else :
                            self.ax[0].plot([self.waves[i],self.waves[i]],[0,ymax],'r')
                    plt.draw()

                    # get input from user on lines to remove
                    for i in range(len(self.pix)) :
                        print('{:3d}{:8.2f}{:8.2f}{:8.2f}{:8.2f}'.format(
                               i, self.pix[i], self.waves[i], diff[i], self.weights[i]))
                    i = input('  enter ID of line to remove (-n for all lines<n, +n for all lines>n, O for new degree, return to continue): ')
                    if i == '' :
                        done = True
                    elif i == 'O' :
                        print('  current degree of fit: {:d}'.format(self.degree))
                        self.degree = int(input('  enter new degree of fit: '))
                        mod = self.getmod()
                    elif '+' in i :
                        self.weights[int(i)+1:] = 0.
                    elif '-' in i :
                        self.weights[0:abs(int(i))] = 0.
                    elif int(i) >= 0 :
                        self.weights[int(i)] = 0.
                    else :
                        print('invalid input')

    def set_spectrum(self,spectrum) :
        """ Set spectrum used to derive fit
        """
        self.spectrum = np.atleast_2d(spectrum)

    def get_spectrum(self) :
        """ Set spectrum used to derive fit
        """
        return self.spectrum 

    def identify(self,spectrum,file=None,wav=None,wref=None,disp=None,display=None,plot=None,rad=5,thresh=10,
                 xmin=None, xmax=None, lags=range(-300,300), nskip=1) :
        """ Given some estimate of wavelength solution and file with lines,
            identify peaks and centroid
        """

        sz=spectrum.shape
        if len(sz) == 1 : 
            spectrum.data = np.atleast_2d(spectrum.data)
            spectrum.uncertainty.array = np.atleast_2d(spectrum.uncertainty.array)
            sz=spectrum.shape
        if xmin is None : xmin=0
        if xmax is None : xmax=sz[-1]
        nrow=sz[0]

        # get initial reference wavelengths if not given
        if wav is None :
            pix=np.arange(sz[-1])
            if self.spectrum is not None :
                # cross correlate with reference image to get pixel shift
                print('  cross correlating with reference spectrum using lags: ', lags)
                fitpeak,shift = image.xcorr(self.spectrum.data,spectrum.data,lags)
                if shift.ndim == 1 :
                    pixshift=(fitpeak+lags[0])[0]
                    print('  Derived pixel shift from input wcal: ',fitpeak+lags[0])
                    if display is not None :
                        display.plotax1.cla()
                        display.plotax1.text(0.05,0.95,'spectrum and reference',transform=display.plotax1.transAxes)
                        for row in range(spectrum.data.shape[0]) :
                            display.plotax1.plot(spectrum.data[row,:],color='m')
                            display.plotax1.plot(self.spectrum.data[row,:],color='g')
                        display.plotax1.set_xlabel('Pixel')
                        display.plotax2.cla()
                        display.plotax2.text(0.05,0.95,'cross correlation: {:8.3f}'.format(pixshift),
                                             transform=display.plotax2.transAxes)
                        display.plotax2.plot(lags,shift)
                        display.plotax1.set_xlabel('Lag')
                        plt.draw()
                        input("  See spectrum and template spectrum (top), cross corrleation(bottom). hit any key to continue")
                    # single shift for all pixels
                    self.pix0 = self.pix0+fitpeak+lags[0]
                    wav=np.atleast_2d(self.wave(image=np.array(sz)))
                else :
                    # different shift for each row
                    wav=np.zeros(sz)
                    cols = np.arange(sz[-1])
                    orders=[]
                    for row in range(wav.shape[0]) : 
                        print('  Derived pixel shift from input wcal for row: {:d} {:d}'.format
                               (row,shift[row,:].argmax()+lags[0]),end='\r')
                        rows=np.zeros(len(cols))+row
                        try : order = self.orders[row]
                        except : order=self.orders[0]
                        orders.append(order)
                        pix0 = self.pix0+fitpeak[row]+lags[0]
                        wav[row,:] = self.model(cols-pix0)/order
                    # ensure we have 2D fit
                    self.type = 'chebyshev2D'
                    self.orders = orders
                    print("")
            else :
                # get dispersion guess from header cards if not given in disp
                if disp is None: disp=hd.header['DISPDW']
                if wref is not None :
                    w0=wref[0]
                    pix0=wref[1]
                else:
                    w0=hd.header['DISPWC']
                    pix0=sz[1]/2 
                wav=np.atleast_2d(w0+(pix-pix0)*disp)

        # open file with wavelengths and read
        if file is not None :
            f=open(ROOT+'/data/lamps/'+file,'r')
            lines=[]
            for line in f :
                if line[0] != '#' :
                    w=float(line.split()[0])
                    # if we have microns, convert to Angstroms
                    if w<10 : w*=10000
                    if w > wav.min() and w < wav.max() : lines.append(w)
            lines=np.array(lines)
            f.close()
        else :
            lines = self.waves
            weights = self.weights
            gd = np.where(weights >0)[0]
            lines = lines[gd]

        # get centroid around expected lines
        x=[]
        y=[]
        waves=[]
        waves_order=[]
        weight=[]
        diff=[]
        if display is not None and  isinstance(display,pyvista.tv.TV) :
            display.ax.cla()
            display.ax.axis('off')
            display.tv(spectrum.data)
        if plot is not None : 
            if type(plot) is matplotlib.figure.Figure :
                plot.clf()
                plt.draw()
                ax1=plot.add_subplot(2,1,1) 
                ax2=plot.add_subplot(2,1,2,sharex=ax1) 
                plot.subplots_adjust(left=0.05,right=0.92, hspace=1.05)
                ax=[ax1,ax2]
                self.fig = plot
                self.ax = ax
            else :
                fig,ax = plt.subplots(2,1,sharex=True,figsize=(14,7))
                fig.subplots_adjust(hspace=1.05)
                self.fig = fig
                self.ax = ax

        if plot is not None : ax[0].cla()
        for row in range(0,nrow,nskip) :
            print('  identifying lines in row: ', row,end='\r')
            if plot is not None :
                ax[0].plot(wav[row,:],spectrum.data[row,:])
                #ax[0].set_yscale('log')
                ax[0].set_ylim(1.,ax[0].get_ylim()[1])
                ax[0].text(0.1,0.9,'row: {:d}'.format(row),transform=ax[0].transAxes)
                ax[0].set_xlabel('Rough wavelength')
                ax[0].set_ylabel('Intensity')
            for line in lines :
                peak=abs(line-wav[row,:]).argmin()
                if isinstance(display,pyvista.tv.TV) :
                    if (peak > xmin+rad) and (peak < xmax-rad) : display.ax.scatter(peak,row,marker='o',color='r',s=2)
                if ( (peak > xmin+rad) and (peak < xmax-rad) and 
                     ((spectrum.data[row,peak-rad:peak+rad]/spectrum.uncertainty.array[row,peak-rad:peak+rad]).max() > thresh) ) :
                    cent = (spectrum.data[row,peak-rad:peak+rad]*np.arange(peak-rad,peak+rad)).sum()/spectrum.data[row,peak-rad:peak+rad].sum()
                    peak = int(cent)
                    cent = (spectrum.data[row,peak-rad:peak+rad]*np.arange(peak-rad,peak+rad)).sum()/spectrum.data[row,peak-rad:peak+rad].sum()
                    if display is not None and  isinstance(display,pyvista.tv.TV) :
                        display.ax.scatter(cent,row,marker='o',color='g',s=2)
                    if plot is not None :
                        ax[0].text(line,1.,'{:7.1f}'.format(line),rotation='vertical',va='top',ha='center')
                    x.append(cent)
                    y.append(row)
                    # we will fit for wavelength*order
                    waves.append(line)
                    try: order = self.orders[row]
                    except: order=self.orders[0]
                    waves_order.append(order)
                    weight.append(1.)
        if plot is not None : 
            if self.model is not None :
                # if we have a solution already, see how good it is (after shift)
                diff=self.wave(pixels=[x,y])-np.array(waves)
                ax[1].cla()
                ax[1].scatter(np.array(waves),diff,s=2,c=y)
                ax[1].text(0.1,0.9,'from previous fit, rms: {:8.3f}'.format(diff.std()),transform=ax[1].transAxes)
                xlim=ax[1].get_xlim()
                ax[1].plot(xlim,[0,0],linestyle=':')
                ax[1].set_ylim(diff.min()-0.5,diff.max()+0.5)
                print("  rms from old fit (with shift): {:8.3f}".format(diff.std()))
            plt.figure(plot.number)
            plt.draw()
            input('  See identified lines. hit any key to continue....')
        self.pix=np.array(x)
        self.y=np.array(y)
        self.waves=np.array(waves)
        self.waves_order=np.array(waves_order)
        self.weights=np.array(weight)
        self.spectrum = spectrum
        print('')

    def scomb(self,hd,wav,average=True,usemask=True) :
        """ Resample onto input wavelength grid
        """
        #output grid
        out=np.zeros(len(wav))
        sig=np.zeros(len(wav))
        mask=np.zeros(len(wav),dtype=bool)
        # raw wavelengths
        w=self.wave(image=np.array(np.atleast_2d(hd.data).shape))
        for i in range(np.atleast_2d(hd).shape[0]) :
            sort=np.argsort(w[i,:])
            if usemask : 
                gd = np.where(~hd.mask[i,sort])
                sort= sort[gd]
            wmin=w[i,sort].min()
            wmax=w[i,sort].max()
            w2=np.abs(wav-wmin).argmin()
            w1=np.abs(wav-wmax).argmin()
            if average :
                out[w2:w1] += ( np.interp(wav[w2:w1],w[i,sort],np.atleast_2d(hd.data)[i,sort]) /
                                np.interp(wav[w2:w1],w[i,sort],np.atleast_2d(hd.uncertainty.array)[i,sort])**2 )
                sig[w2:w1] += 1./np.interp(wav[w2:w1],w[i,sort],np.atleast_2d(hd.uncertainty.array)[i,sort])**2 
            else :
                out[w2:w1] += np.interp(wav[w2:w1],w[i,sort],np.atleast_2d(hd.data)[i,sort])
                sig[w2:w1] += np.interp(wav[w2:w1],w[i,sort],np.atleast_2d(hd.uncertainty.array**2)[i,sort])
        if average :
            out = out / sig
        else :
            sig = np.sqrt(sig)
        return CCDData(out,uncertainty=StdDevUncertainty(sig),mask=mask,header=hd.header,unit='adu')

    def save(self,file) :
        """ Save object to file
        """
        try : delattr(self,'fig')
        except: pass
        try : delattr(self,'ax')
        except: pass
        f=open(file,'wb')
        pickle.dump(self,f)
        f.close()

class Trace() :
    """ Class for spectral traces
    """

    def __init__ (self,inst=None, type='poly',order=2,pix0=0,rad=5,spectrum=None,model=None,sc0=None,rows=None,lags=None,channel=None) :
        self.type = type
        self.order = order
        self.pix0 = pix0
        self.spectrum = spectrum
        self.rad = rad
        if inst == 'TSPEC' :
            self.order = 3
            self.rows = [[135,235],[295,395],[435,535],[560,660],[735,830]]
            self.lags = range(-75,75) 
        elif inst == 'DIS' :
            if channel == 0 : self.rows=[[215,915]]
            elif channel == 1 : self.rows=[[100,800]]
            else : raise ValueError('need to specify channel')
            self.lags = range(-300,300) 
        elif inst == 'ARCES' :
            self.lags = range(-10,10) 
        if rows is not None : self.rows=rows
        if lags is not None : self.lags=lags
        if model is not None : self.model=model
        if sc0 is not None : self.sc0=sc0

    def trace(self,hd,srows,sc0=None,plot=None,thresh=20) :
        """ Trace a spectrum from starting position
        """

        fitter=fitting.LinearLSQFitter()
        if self.type == 'poly' :
            mod=models.Polynomial1D(degree=self.order)
        else :
            raise ValueError('unknown fitting type: '+self.type)
            return

        nrow = hd.data.shape[0]
        ncol = hd.data.shape[1]
        if sc0 is None : self.sc0 = int(ncol/2)
        else : self.sc0 = sc0
        self.spectrum = hd[:,self.sc0]
        self.spectrum.data[self.spectrum.data<0] = 0.
        rows = np.arange(nrow)
        ypos = np.zeros(ncol)
        ysum = np.zeros(ncol)
        yvar = np.zeros(ncol)
        ymask = np.zeros(ncol,dtype=bool)

        # we want to handle multiple traces, so make sure srows is iterable
        if type(srows ) is int or type(srows) is float : srows=[srows]
        oldmodel=copy.copy(self.model)
        self.model=[]
        if plot is not None : 
            plot.clear()
            plot.tv(hd)

        rad = self.rad-1
        for irow,srow in enumerate(srows) :
            print('  Tracing row: {:d}'.format(int(srow)),end='\r')
            sr=copy.copy(srow)
            sr=int(round(sr))
            sr=hd.data[sr-rad:sr+rad+1,self.sc0].argmax()+sr-rad
            # march left from center
            for col in range(self.sc0,0,-1) :
                # centroid
                cr=sr-rad+hd.data[sr-rad:sr+rad+1,col].argmax()
                ysum[col] = np.sum(hd.data[cr-rad:cr+rad+1,col]) 
                ypos[col] = np.sum(rows[cr-rad:cr+rad+1]*hd.data[cr-rad:cr+rad+1,col]) / ysum[col]
                yvar[col] = np.sum(hd.uncertainty.array[cr-rad:cr+rad+1,col]**2) 
                ymask[col] = np.any(hd.mask[cr-rad:cr+rad+1,col]) 
                # if centroid is too far from starting guess, mask as bad
                if np.abs(ypos[col]-sr) > rad/2. : ymask[col] = True
                # use this position as starting center for next if above threshold S/N
                if (not ymask[col]) & np.isfinite(ysum[col]) & (ysum[col]/np.sqrt(yvar[col]) > thresh)  : sr=int(round(ypos[col]))
            sr=copy.copy(srow)
            sr=int(round(sr))
            sr=hd.data[sr-rad:sr+rad+1,self.sc0].argmax()+sr-rad
            # march right from center
            for col in range(self.sc0+1,ncol,1) :
                # centroid
                cr=sr-rad+hd.data[sr-rad:sr+rad+1,col].argmax()
                ysum[col] = np.sum(hd.data[cr-rad:cr+rad+1,col]) 
                ypos[col] = np.sum(rows[cr-rad:cr+rad+1]*hd.data[cr-rad:cr+rad+1,col]) / ysum[col]
                yvar[col] = np.sum(hd.uncertainty.array[cr-rad:cr+rad+1,col]**2) 
                ymask[col] = np.any(hd.mask[cr-rad:cr+rad+1,col]) 
                if np.abs(ypos[col]-sr) > rad/2. : ymask[col] = True
                # use this position as starting center for next if above threshold S/N
                if (not ymask[col]) & np.isfinite(ysum[col]) & (ysum[col]/np.sqrt(yvar[col]) > thresh)  : sr=int(round(ypos[col]))

            cols=np.arange(ncol)
            gd = np.where((~ymask) & (ysum/np.sqrt(yvar)>thresh) )[0]
            model=(fitter(mod,cols[gd],ypos[gd]))

            # reject outlier points (>1 pixel) and refit
            res = model(cols)-ypos
            gd = np.where((~ymask) & (ysum/np.sqrt(yvar)>thresh) & (np.abs(res)<1))[0]
            model=(fitter(mod,cols[gd],ypos[gd]))
            if len(gd) < 10 : 
                print('  failed trace for row: {:d}, using old model'.format(irow))
                model=copy.copy(oldmodel[irow])
            self.model.append(model)

            if plot : 
                plot.ax.scatter(cols,ypos,marker='o',color='r',s=4) 
                plot.ax.scatter(cols[gd],ypos[gd],marker='o',color='g',s=4) 
                plot.ax.plot(cols,model(cols),color='m')
                #plt.pause(0.05)

        self.pix0=0
        print("")
        if plot : input('  See trace. Hit any key to continue....')

    def retrace(self,hd,plot=None,thresh=20) :
        """ Retrace starting with existing model
        """
        self.find(hd)
        srows = []
        for row in range(len(self.model)) :
            srows.append(self.model[row](self.sc0))
        self.trace(hd,srows,plot=plot,thresh=thresh)
     
    def find(self,hd,lags=None,plot=None) :
        """ Determine shift from existing trace to input frame
        """
        if lags is None : lags = self.lags
       
        im=copy.deepcopy(hd.data)
        # if we have a window, zero array outside of window
        spec=im[:,self.sc0]
        try:
            spec[:self.rows[0]] = 0.  
            spec[self.rows[1]:] = 0.  
        except: pass
        fitpeak,shift = image.xcorr(self.spectrum,spec,lags)
        pixshift=(fitpeak+lags[0])[0]
        print('  traces shift: ', fitpeak+lags[0])
        if plot is not None :
            plot.clear()
            plot.tv(im)
            plot.plotax1.cla()
            plot.plotax1.text(0.05,0.95,'obj and ref cross-section',transform=plot.plotax1.transAxes)
            plot.plotax1.plot(self.spectrum.data/self.spectrum.data.max())
            plot.plotax1.plot(im[:,self.sc0]/im[:,self.sc0].max())
            plot.plotax1.set_xlabel('row')
            plot.plotax2.cla()
            plot.plotax2.text(0.05,0.95,'cross correlation {:8.3f}'.format(pixshift),
                              transform=plot.plotax2.transAxes)
            plot.plotax2.plot(lags,shift)
            plot.plotax2.set_xlabel('lag')
            plt.draw()
            input('  See spectra and cross-correlation. Hit any key to continue....')
        self.pix0=fitpeak+lags[0]
        return fitpeak+lags[0]
 
    def extract(self,hd,rad=None,scat=False,plot=None,medfilt=None) :
        """ Extract spectrum given trace(s)
        """
        if rad is None : rad=self.rad
        nrows=hd.data.shape[0]
        ncols=hd.data.shape[-1]
        spec = np.zeros([len(self.model),hd.data.shape[1]])
        sig = np.zeros([len(self.model),hd.data.shape[1]])
        mask = np.zeros([len(self.model),hd.data.shape[1]],dtype=bool)

        if plot is not None:
            plot.clear()
            plot.tv(hd)

        for i,model in enumerate(self.model) :
            print('  extracting aperture {:d}'.format(i),end='\r')
            cr=model(np.arange(ncols))+self.pix0
            icr=np.round(cr).astype(int)
            rfrac=cr-icr+0.5   # add 0.5 because we rounded
            rlo=[]
            rhi=[]
            for col in range(ncols) :
                r1=icr[col]-rad
                r2=icr[col]+rad
                # sum inner pixels directly, outer pixels depending on fractional pixel location of trace
                if r1>=0 and r2<nrows :
                    spec[i,col]=np.sum(hd.data[r1+1:r2,col])
                    sig[i,col]=np.sum(hd.uncertainty.array[r1+1:r2,col]**2)
                    spec[i,col]+=hd.data[r1,col]*(1-rfrac[col])
                    sig[i,col]+=hd.uncertainty.array[r1,col]**2*(1-rfrac[col])
                    spec[i,col]+=hd.data[r2,col]*rfrac[col]
                    sig[i,col]+=hd.uncertainty.array[r2,col]**2*rfrac[col]
                    sig[i,col]=np.sqrt(sig[i,col])
                    mask[i,col] = np.any(hd.mask[r1:r2+1,col]) 
                if plot is not None :
                    rlo.append(r1)
                    rhi.append(r2-1)
            if medfilt is not None :
                boxcar = Box1DKernel(medfilt)
                median = convolve(spec[i,:],boxcar,boundary='extend')
                spec[i,:]/=median
                sig[i,:]/=median

            if plot is not None :
                if i%2 == 0 : color='b'
                else : color='m'
                plot.ax.plot(range(ncols),cr,color='g',linewidth=3)
                plot.ax.plot(range(ncols),rlo,color=color,linewidth=1)
                plot.ax.plot(range(ncols),rhi,color=color,linewidth=1)
                plt.draw()
        if plot is not None : input('  See extraction window(s). Hit any key to continue....')
        print("")
        return CCDData(spec,uncertainty=StdDevUncertainty(sig),mask=mask,header=hd.header,unit='adu')
  
    def extract2d(self,hd,rows=None,plot=None) :
        """  Extract 2D spectrum given trace(s)
             Assumes all requests row uses same trace, just offset, not a 2D model for traces
        """
        nrows=hd.data.shape[0]
        ncols=hd.data.shape[-1]
        out=[]
        if plot is not None:
            plot.clear()
            plot.tv(hd)
        for model in self.model :
            if plot is not None :
                plot.ax.plot([0,ncols],[self.rows[0],self.rows[0]],color='g')
                plot.ax.plot([0,ncols],[self.rows[1],self.rows[1]],color='g')
                plt.draw()
            outrows=np.arange(self.rows[0],self.rows[1])
            noutrows=len(range(self.rows[0],self.rows[1]))
            spec=np.zeros([noutrows,ncols])
            sig=np.zeros([noutrows,ncols])
            cr=model(np.arange(ncols))
            cr-=cr[self.sc0]
            for col in range(ncols) :
                spec[:,col] = np.interp(outrows+cr[col],np.arange(nrows),hd.data[:,col])
                sig[:,col] = np.sqrt(np.interp(outrows+cr[col],np.arange(nrows),hd.uncertainty.array[:,col]**2))
            out.append(CCDData(spec,StdDevUncertainty(sig),unit='adu'))
        if plot is not None: input('  enter something to continue....')

        if len(out) == 1 : return out[0]
        else : return out

    def save(self,file) :
        """ Save object to file
        """
        try : delattr(self,'ax')
        except: pass
        f=open(file,'wb')
        pickle.dump(self,f)
        f.close()

def mash(hd,sp=None,bks=None) :
    """
    Mash image into spectra using requested window
    """
    if sp is None :
        sp=[0,hd.data.shape[0]]
    obj = hd.data[sp[0]:sp[1]].sum(axis=0)
    obj = hd.data[sp[0]:sp[1]].sum(axis=0)

    if bks is not None :
        back=[]
        for bk in bks :
           tmp=np.median(data[bk[0]:bk[1]],axis=0)
           back.append(tmp)
        obj-= np.mean(back,axis=0)

    return obj

def wavecal(hd,file=None,wref=None,disp=None,wid=[3],rad=5,snr=3,degree=2,wcal0=None,thresh=100,type='poly'):
    """
    Get wavelength solution for single 1D spectrum
    """

    # choose middle row +/ 5 rows
    sz=hd.data.shape
    spec=hd.data[int(sz[0]/2)-5:int(sz[0]/2)+5,:].sum(axis=0)
    spec=spec-scipy.signal.medfilt(spec,kernel_size=101)
    pix = np.arange(len(spec))

    fig,ax = plt.subplots(2,1,sharex=True,figsize=(14,6))
    ax[0].plot(spec)

    # get wavelength guess from input WaveCal if given, else use wref and dispersion, else header
    if wcal0 is not None :
        lags=range(-300,300)
        fitpeak,shift = image.xcorr(wcal0.spectrum,spec,lags)
        wnew=copy.deepcopy(wcal0)
        wnew.pix0 = wcal0.pix0+shift.argmax()+lags[0]
        print('  Derived pixel shift from input wcal0: ',shift.argmax()+lags[0])
        wav=wnew.wave(pix)
    else :
        # get dispersion guess from header cards if not given in disp
        if disp is None: disp=hd.header['DISPDW']
        if wref is not None :
            w0=wref[0]
            pix0=wref[1]
            wav=w0+(pix-pix0)*disp
        else:
            w0=hd.header['DISPWC']
            pix0=sz[1]/2 
            wav=w0+(pix-pix0)*disp
    ax[1].plot(wav,spec)

    # open file with wavelengths and read
    f=open(file,'r')
    lines=[]
    for line in f :
        if line[0] != '#' :
            w=float(line.split()[0])
            name=line[10:].strip()
            lpix=abs(w-wav).argmin()
            if lpix > 1 and lpix < sz[1]-1 :
                ax[0].text(lpix,0.,'{:7.1f}'.format(w),rotation='vertical',va='top',ha='center')
                lines.append(w)
    lines=np.array(lines)
    f.close()

    # get centroid around expected lines
    cents=[]
    for line in lines :
        peak=abs(line-wav).argmin()
        if (peak > rad) and (peak < sz[1]-rad) and (spec[peak-rad:peak+rad].max() > thresh) :
            print(peak,spec[peak-rad:peak+rad].max())
            cents.append((spec[peak-rad:peak+rad]*np.arange(peak-rad,peak+rad)).sum()/spec[peak-rad:peak+rad].sum())
    cents=np.array(cents)
    print('  cents:', cents)

    waves=[]
    weight=[]
    print('  Centroid  W0  Wave')
    for cent in cents :
        w=wav[int(cent)]
        ax[0].plot([cent,cent],[0,10000],'k')
        print('  {:8.2f}{:8.2f}{:8.2f}'.format(cent, w, lines[np.abs(w-lines).argmin()]))
        waves.append(lines[np.abs(w-lines).argmin()])
        weight.append(1.)
    waves=np.array(waves)
    weight=np.array(weight)

    # set up new WaveCal object
    pix0 = int(sz[1]/2)
    wcal = WaveCal(order=degree,type=type,spectrum=spec,pix0=pix0)

    # iterate allowing for interactive removal of points
    done = False
    ymax = ax[0].get_ylim()[1]
    while not done :
        gd=np.where(weight>0.)[0]
        bd=np.where(weight<=0.)[0]
        wcal.fit(cents[gd],waves[gd],weights=weight[gd])

        # plot
        ax[1].cla()
        ax[1].plot(cents[gd],wcal.wave(cents[gd])-waves[gd],'go')
        if len(bd) > 0 : ax[1].plot(cents[bd],wcal.wave(cents[bd])-waves[bd],'ro')
        diff=wcal.wave(cents[gd])-waves[gd]
        ax[1].set_ylim(diff.min()-1,diff.max()+1)
        for i in range(len(cents)) :
            ax[1].text(cents[i],wcal.wave(cents[i])-waves[i],'{:2d}'.format(i),va='top',ha='center')
            if weight[i] > 0 :
              ax[0].plot([cents[i],cents[i]],[0,ymax],'g')
            else :
              ax[0].plot([cents[i],cents[i]],[0,ymax],'r')
        plt.draw()

        # get input from user on lines to remove
        for i in range(len(cents)) :
            print('  {:3d}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}'.format(
                   i, cents[i], wcal.wave(cents[i]), waves[i], waves[i]-wcal.wave(cents[i]),weight[i]))
        print('  rms: {:8.2f} Anstroms'.format(diff.std()))
        i = input('enter ID of line to remove (-n for all lines<n, +n for all lines>n, return to continue): ')
        if i is '' :
            done = True
        elif '+' in i :
            weight[int(i)+1:] = 0.
        elif '-' in i :
            weight[0:abs(int(i))] = 0.
        elif int(i) >= 0 :
            weight[int(i)] = 0.
        else :
            print('invalid input')

    plt.close()

    return wcal.wave(pix),wcal

def fluxcal(obs,wobs,file=None) :
    """
    flux calibration
    """

    fluxdata=ascii.read(file)
    stan=np.interp(wobs,fluxdata['col1'],fluxdata['col2'])
    return stan/obs
  

def trace(hd,apertures=None,pix0=1024) : 
    """ Get all traces
        apertures is a list of row numbers at pixel 1024
    """
    alltr=[]
    for i in range(len(apertures)) :
        tr=Trace()
        print('tracing aperture {:d}'.format(i),end='\r')
        sr=apertures[i]
        tr.trace(hd,pix0,sr)
        alltr.append(tr)

    return alltr

def extract(hd,apertures) :
    """ Do all extractions
    """
    spec = np.zeros([len(apertures),hd.data.shape[1]])
    for i,order in enumerate(apertures) :
        print('extracting aperture {:d}'.format(i),end='\r')
        spec[i] = order.extract(hd)

    return spec


