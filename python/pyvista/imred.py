import numpy as np
import astropy
import code
import copy
from astropy import units
from astropy.nddata import CCDData, NDData, StdDevUncertainty
from astropy.nddata import NDData
from astropy.io import fits
from astropy.io import ascii
from astropy.modeling import models, fitting
from astropy.convolution import convolve, Box1DKernel, Box2DKernel, Box2DKernel
import ccdproc
import scipy.signal
import yaml
import sys

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

import matplotlib.pyplot as plt
import glob
import bz2
import os
import pdb
from pyvista import image
from pyvista import tv
try: 
    import pyds9
except:
    print('pyds9 is not available, proceeding')

ROOT = os.path.dirname(os.path.abspath(__file__)) + '/../../'



class Reducer() :
    """ Class for reducing images of a given instrument
    """
    def __init__(self,inst=None,dir='./',root='*',formstr='{04d}',gain=1,rn=0.,verbose=True,nfowler=1) :
        """  Initialize reducer with information about how to reduce
        """
        self.dir=dir
        self.root=root
        self.verbose=verbose
        self.inst=inst
        self.badpix=None
        self.scat=None
        self.mask=None

        # we will allow for instruments to have multiple channels, so everything goes in lists
        self.channels=['']
        if type(gain) is list : self.gain=gain
        else : self.gain = [gain]
        if type(rn) is list : self.rn=rn
        else : self.rn = [rn]
        if type(formstr) is list : self.formstr=formstr
        else : self.formstr=[formstr]
       
        # Read instrument configuation from YAML configuration file 
        if inst is not None :
            if inst.find('/') < 0 :
                config = yaml.load(open(ROOT+'/data/'+inst+'/'+inst+'_config.yml','r'), Loader=yaml.FullLoader)
            else :
                config = yaml.load(open(inst+'_config.yml','r'), Loader=yaml.FullLoader)
            self.channels=config['channels']
            self.formstr=config['formstr']
            self.gain=config['gain']
            self.rn=config['rn']/np.sqrt(nfowler)
            self.crbox=config['crbox']
            self.biastype=config['biastype']
            self.biasbox=[]
            for box in config['biasbox'] :
                self.biasbox.append(image.BOX(xr=box[0],yr=box[1]) )
            self.trimbox=[]
            for box in config['trimbox'] :
                self.trimbox.append(image.BOX(xr=box[0],yr=box[1]) )
            self.normbox=[]
            for box in config['normbox'] :
                self.normbox.append(image.BOX(xr=box[0],yr=box[1]) )
            try: self.scat=config['scat']
            except : pass
           
            # Add bad pixel mask if it exists
            try: self.mask=fits.open(ROOT+'/data/'+inst+'/'+inst+'_mask.fits')[0].data.astype(bool)
            except: pass

        # save number of chips for convenience
        self.nchip = len(self.formstr)

        # output setup if verbose
        if self.verbose :
            if inst is not None : print('INSTRUMENT: {:s}'.format(inst))
            for form in self.formstr :
                print('  will use format:  {:s}/{:s}{:s}.fits*'.format(self.dir,self.root,form))
            print('         gain:  {}    rn: {}'.format(self.gain,self.rn))
            print('  Biastype : {:d}'.format(self.biastype))
            print('  Bias box: ')
            for box in self.biasbox :
                box.show()
            print('  Trim box: ')
            for box in self.trimbox :
                box.show()
            print('  Norm box: ')
            for box in self.normbox :
                box.show()

    def rd(self,num, ext=0) :
        """ Read an image

        Args :
            num (str or int) : name or number of image to read
        Returns :
            image (CCDData ) : CCDData object
        """
        out=[]
        # loop over different channels (if any)
        idet=0 
        for form,gain,rn in zip(self.formstr,self.gain,self.rn) :
            # find the files that match the directory/format
            if type(num) is int :
                search=self.dir+'/'+self.root+form.format(num)+'.fits*'
            elif type(num) is str or type(num) is np.str_ :
                if num.find('/') >= 0 :
                    search=num+'*'
                else :
                    search=self.dir+'/*'+num+'*'
            else :
                print('stopping in rd... num:',num)
                pdb.set_trace()
            file=glob.glob(search)
            if len(file) == 0 : 
                print('cannot find file matching: '+search)
                return
            elif len(file) > 1 : 
                if self.verbose : print('more than one match found, using first!',file)
            file=file[0]

            # read the file into a CCDData object
            if self.verbose : print('  Reading file: {:s}'.format(file)) 
            try : im=CCDData.read(file,hdu=ext,unit='adu')
            except : raise RuntimeError('Error reading file: {:s}'.format(file))
            im.header['FILE'] = os.path.basename(file)
            if 'OBJECT' not in im.header :
                try: im.header['OBJECT'] = im.header['OBJNAME']
                except KeyError : im.header['OBJECT'] = im.header['FILE']

            # Add uncertainty (will be in error if there is an overscan, but redo with overscan subraction later)
            data=copy.copy(im.data)
            data[data<0] = 0.
            im.uncertainty = StdDevUncertainty(np.sqrt( data/gain + (rn/gain)**2 ))

            # Add mask
            if self.mask is not None : im.mask = self.mask
            else : im.mask = np.zeros(im.data.shape,dtype=bool)
            if self.badpix is not None :
                for badpix in self.badpix[idet] :
                    badpix.setval(im.mask,True)

            out.append(im)
            idet+=1

        # return the data
        if len(out) == 1 : return out[0]
        else : return out
            
    def overscan(self,im,display=None) :
        """ Overscan subtraction
        """
        if self.biastype < 0 : return

        if type(im) is not list : ims=[im]
        else : ims = im
       
        for ichan,(im,gain,rn,biasbox) in enumerate(zip(ims,self.gain,self.rn,self.biasbox)) :
            if display is not None : 
                display.tv(im)
            if self.biastype == 0 :
                b=biasbox.mean(im.data)
                if self.verbose: print('  subtracting overscan: ', b)
                if display is not None : 
                    display.tvbox(0,0,box=biasbox)
                    if ichan %2 == 0 : ax=display.plotax1
                    else : ax=display.plotax2
                    ax.cla()
                    ax.plot(np.mean(im.data[:,biasbox.xmin:biasbox.xmax],axis=1))
                    ax.text(0.05,0.95,'Overscan mean',transform=ax.transAxes)
                    ax.set_xlabel('Row')
                    display.fig.canvas.draw_idle()
                    plt.draw()
                    get=input("  See bias box and cross section. Hit any key to continue")
                    if get == 'i' : code.interact(local=locals())
                    elif get == 'q' : quit()
                    elif get == 'p' : pdb.set_trace()
                im.data = im.data.astype(float)-b
                im.header.add_comment('subtracted overscan: {:f}'.format(b))
            #elif det.biastype == 1 :
            #    over=np.median(hdu[ext].data[:,det.biasbox.xmin:det.biasbox.xmax],axis=1)
            #    boxcar = Box1DKernel(10)
            #    over=convolve(over,boxcar,boundary='extend')
            #    over=image.stretch(over,ncol=hdu[ext].data.shape[1])
            #    hdu[ext].data -= over

            # Add uncertainty (redo from scratch after overscan)
            data=copy.copy(im.data)
            data[data<0] = 0.
            im.uncertainty = StdDevUncertainty(np.sqrt( data/gain + (rn/gain)**2 ))

    def trim(self,im) :
        """ Trim image by masking non-trimmed pixels
            Need to preserve image size to match reference/calibration frames, etc.
        """
        if type(im) is not list : ims=[im]
        else : ims = im

        for  im,trimbox in zip(ims,self.trimbox) :
            tmp = np.ones(im.mask.shape,dtype=bool)
            trimbox.setval(tmp,False)
            im.mask = np.logical_or(im.mask,tmp)


    def bias(self,im,superbias=None) :
         """ Superbias subtraction
         """
         # only subtract if we are given a superbias!
         if superbias is None : return im

         # work with lists so that we can handle multi-channel instruments
         if type(im) is not list : ims=[im]
         else : ims = im
         if type(superbias) is not list : superbiases=[superbias]
         else : superbiases = superbias
         out=[]
         for im,bias in zip(ims,superbiases) :
             if self.verbose : print('  subtracting superbias...')
             out.append(ccdproc.subtract_bias(im,bias))
         if len(out) == 1 : return out[0]
         else : return out

    def dark(self,im,superdark=None) :
         """ Superdark subtraction
         """
         # only subtract if we are given a superdark!
         if superdark is None : return im

         # work with lists so that we can handle multi-channel instruments
         if type(im) is not list : ims=[im]
         else : ims = im
         if type(superdark) is not list : superdarks=[superdark]
         else : superdarks = superdark
         out=[]
         for im,dark in zip(ims,superdarks) :
             if self.verbose : print('  subtracting superdark...')
             out.append(ccdproc.subtract_dark(im,dark,exposure_time='EXPTIME',exposure_unit=units.s))
         if len(out) == 1 : return out[0]
         else : return out

    def flat(self,im,superflat=None,display=None) :
         """ Flat fielding
         """
         # only flatfield if we are given a superflat!
         if superflat is None : return im

         if type(im) is not list : ims=[im]
         else : ims = im
         if type(superflat) is not list : superflats=[superflat]
         else : superflats = superflat
         out=[]
         for im,flat in zip(ims,superflats) :
             if self.verbose : print('  flat fielding...')
             if display is not None : 
                 display.tv(im)
             corr = ccdproc.flat_correct(im,flat)
             out.append(corr)
             if display is not None : 
                 display.tv(corr)
                 #plot central crossections
                 display.plotax1.cla()
                 dim=corr.data.shape
                 col = int(dim[1]/2)
                 row = corr.data[:,col]
                 display.plotax1.plot(row)
                 min,max=tv.minmax(row,low=5,high=5)
                 display.plotax1.set_ylim(min,max)
                 display.plotax1.set_xlabel('row')
                 display.plotax1.text(0.05,0.95,'Column {:d}'.format(col),transform=display.plotax1.transAxes)
                 display.plotax2.cla()
                 row = int(dim[0]/2)
                 col = corr.data[row,:]
                 min,max=tv.minmax(col,low=10,high=10)
                 display.plotax2.plot(col)
                 display.plotax2.set_xlabel('col')
                 display.plotax2.text(0.05,0.95,'Row {:d}'.format(row),transform=display.plotax2.transAxes)
                 display.plotax2.set_ylim(min,max)
                 input("  See flat-fielded image and original with - key. Hit any key to continue")
         if len(out) == 1 : return out[0]
         else : return out


    def scatter(self,im,scat=None,display=None,smooth=3,smooth2d=31) :
        """ Removal of scattered light (for multi-order/object spectrograph)
        """
        if scat is None : return

        print('  estimating scattered light ...')
        boxcar = Box1DKernel(smooth)
        points=[]
        values=[]
        nrows = im.data.shape[0]
        ncols = im.data.shape[-1]

        # find minima in each column, and save location and value
        for col in range(0,ncols,scat) :
            print('    column: {:d}'.format(col),end='\r')
            yscat = scipy.signal.find_peaks(-convolve(im.data[:,col],boxcar))[0]
            for y in yscat :
                if im.mask is None or not im.mask[y,col] :
                    points.append([y,col])
                    values.append(im.data[y,col])

        # fit surface to the minimum values
        print('    fitting surface ...')
        grid_x, grid_y = np.mgrid[0:nrows,0:ncols]

        # smooth and reject outlying points
        boxcar = Box2DKernel(smooth2d)
        grid_z=convolve(scipy.interpolate.griddata(points,values,(grid_x,grid_y),
                        method='cubic',fill_value=0.),boxcar)
        # go back and try to reject outliers
        print('    rejecting points ...')
        points_gd=[]
        values_gd=[]
        for point,value in zip(points,values) :
            if value < 1.1*grid_z[point[0],point[1]] :
                points_gd.append(point)
                values_gd.append(value)

        # refit surface
        print('    refitting surface ...')
        grid_z=convolve(scipy.interpolate.griddata(points_gd,values_gd,(grid_x,grid_y),
                        method='cubic',fill_value=0.),boxcar)

        if display is not None :
            display.clear()
            display.tv(im)
            points=np.array(points)
            display.ax.scatter(points[:,1],points[:,0],color='r',s=3)
            points_gd=np.array(points_gd)
            display.ax.scatter(points_gd[:,1],points_gd[:,0],color='g',s=3)
            input("  See image with scattered light points. Hit any key to continue".format(im))
            display.clear()
            display.tv(im)
            display.tv(grid_z)
            col=int(im.shape[-1]/2)
            display.plotax1.cla()
            display.plotax1.plot(im.data[:,col])
            display.plotax1.plot(grid_z[:,col])
            plt.draw()
            input("  See scattered light image. Hit any key to continue".format(im))

        im.data -= grid_z

    def crrej(self,im,crbox=None,nsig=5,display=None) :
        """ CR rejection
        """
        if crbox is None: return im
        if type(im) is not list : ims=[im]
        else : ims = im
        out=[]
        for i,im in enumerate(ims) :
            if display is not None : 
                display.tv(im)
            if crbox == 'lacosmic':
                if self.verbose : print('  zapping CRs with ccdproc.cosmicray_lacosmic')
                im= ccdproc.cosmicray_lacosmic(im)
            else :
                if self.verbose : print('  zapping CRs with filter [{:d},{:d}]...'.format(*crbox))
                image.zap(im,crbox,nsig=nsig)
            if display is not None : 
                display.tv(im)
                input("  See CR-zapped image and original with - key. Hit any key to continue")
            out.append(im)
        if len(out) == 1 : return out[0]
        else : return out

    def badpix_fix(self,im,val=0.) :
        """ Replace bad pixels
        """
        if val is None : return
        if type(im) is not list : ims=[im]
        else : ims = im
        for i, im in enumerate(ims) :
            if im.mask is not None :
                 bd=np.where(im.mask)
                 ims[i].data[bd[0],bd[1]] = val
                 ims[i].uncertainty.array[bd[0],bd[1]] = np.inf

    def display(self,display,id) :

        im = self.reduce(id)
        if type(im) is not list : ims=[im]
        else : ims = im
        for i, im in enumerate(ims) :
            display.tv(im)

    def reduce(self,num,crbox=None,superbias=None,superdark=None,superflat=None,scat=None,badpix=None,return_list=False,display=None) :
        """ Full reduction
        """
        im=self.rd(num)
        self.overscan(im,display=display)
        im=self.crrej(im,crbox=crbox,display=display)
        im=self.bias(im,superbias=superbias)
        im=self.dark(im,superdark=superdark)
        self.scatter(im,scat=scat,display=display)
        im=self.flat(im,superflat=superflat,display=display)
        self.badpix_fix(im,val=badpix)
        self.trim(im)
        if return_list and type(im) is not list : im=[im]
        return im

    def write(self,im,name,overwrite=True,trim=False) :
        """ write out image, deal with multiple channels 
        """

        if type(im) is not list : ims=[im]
        else : ims = im
        for image,trimbox in zip(ims,self.trimbox) :
            if trim :
                image.data = image.data[trimbox.ymin:trimbox.ymin+trimbox.nrow(),
                                        trimbox.xmin:trimbox.xmin+trimbox.ncol()]
                image.uncertainty.array = image.uncertainty.array[trimbox.ymin:trimbox.ymin+trimbox.nrow(),
                                                                  trimbox.xmin:trimbox.xmin+trimbox.ncol()]
                if image.mask is not None :
                    image.mask = image.mask[trimbox.ymin:trimbox.ymin+trimbox.nrow(),
                                            trimbox.xmin:trimbox.xmin+trimbox.ncol()]
        if self.nchip > 1 :
            for i,frame in enumerate(im) : frame.write(name+'_'+self.channels[i]+'.fits',overwrite=overwrite)
        else :
            im.write(name+'.fits',overwrite=overwrite)

    def getcube(self,ims,**kwargs) :
        """ Read images into data cube
        """
        # create list of images, reading and overscan subtracting
        allcube = []
        for im in ims :
            if type(im) is not astropy.nddata.CCDData :
                data = self.reduce(im, **kwargs)
            allcube.append(data)

        # if just one frame, put in 2D list anyway so we can use same code, allcube[nframe][nchip]
        if self.nchip == 1 :
            allcube=[list(i) for i in zip(*[allcube])]

        return allcube

    def sum(self,ims, return_list=False, **kwargs) :
        """ Coadd input images
        """
        allcube = self.getcube(ims, **kwargs)
        nframe = len(allcube)
        
        out=[]
        for chip in range(self.nchip) :
            datacube = []
            varcube = []
            maskcube = []
            for im in range(nframe) :
                datacube.append(allcube[im][chip].data)
                varcube.append(allcube[im][chip].uncertainty.array**2)
                maskcube.append(allcube[im][chip].mask)
            sum = np.sum(np.array(datacube),axis=0)
            sig = np.sqrt(np.sum(np.array(varcube),axis=0))
            mask = np.any(maskcube,axis=0)
            out.append(CCDData(sum,uncertainty=StdDevUncertainty(sig),mask=mask,unit='adu'))
        
        # return the frame
        if len(out) == 1 : 
           if return_list : return [out[0]]
           else : return out[0]
        else : return out

    def median(self,ims, normalize=False,display=None,div=True,return_list=False, **kwargs) :
        """ Combine images from list of images 
        """
        # create list of images, reading and overscan subtracting
        allcube = self.getcube(ims,**kwargs)
        nframe = len(allcube)

        # do the combination
        out=[] 
        for chip in range(self.nchip) :
            datacube = []
            varcube = []
            maskcube = []
            allnorm = []
            for im in range(nframe) :
                if normalize :
                    norm=self.normbox[chip].mean(allcube[im][chip].data)
                    allnorm.append(norm)
                    allcube[im][chip].data /= norm
                    allcube[im][chip].uncertainty.array /= norm
                datacube.append(allcube[im][chip].data)
                varcube.append(allcube[im][chip].uncertainty.array**2)
                maskcube.append(allcube[im][chip].mask)
            if self.verbose: print('  median combining data....')
            med = np.median(np.array(datacube),axis=0)
            if self.verbose: print('  calculating uncertainty....')
            sig = 1.253 * np.sqrt(np.mean(np.array(varcube),axis=0)/nframe)
            mask = np.any(maskcube,axis=0)
            comb=CCDData(med,header=allcube[im][chip].header,uncertainty=StdDevUncertainty(sig),mask=mask,unit='adu')
            if normalize: comb.meta['MEANNORM'] = np.array(allnorm).mean()
            out.append(comb)

            # display final combined frame and individual frames relative to combined
            if display :
                display.clear()
                display.tv(comb,sn=True)
                display.tv(comb)
                gd=np.where(comb.mask == False)
                min,max=tv.minmax(med[gd[0],gd[1]],low=10,high=10)
                display.plotax1.hist(med[gd[0],gd[1]],bins=np.linspace(min,max,100),histtype='step')
                display.fig.canvas.draw_idle()
                get = input("  See final image, use - key for S/N image. Hit any key to continue")
                if get == 'i' : code.interact(local=locals())
                elif get == 'q' : sys.exit()
                elif get == 'p' : pdb.set_trace()
                for i,im in enumerate(ims) :
                    min,max=tv.minmax(med[gd[0],gd[1]],low=5,high=5)
                    display.fig.canvas.draw_idle()
                    if div :
                        display.plotax2.hist((allcube[i][chip].data/med)[gd[0],gd[1]],bins=np.linspace(0.5,1.5,100),histtype='step')
                        display.tv(allcube[i][chip].data/med,min=0.5,max=1.5)
                        input("    see image: {} divided by master, hit any key to continue".format(im))
                    else :
                        delta=5*self.rn[chip]
                        display.plotax2.hist((allcube[i][chip].data-med)[gd[0],gd[1]],bins=np.linspace(-delta,delta,100),histtype='step')
                        display.tv(allcube[i][chip].data-med,min=-delta,max=delta)
                        input("    see image: {} minus master, hit any key to continue".format(im))

        # return the frame
        if len(out) == 1 :
           if return_list : return [out[0]]
           else : return out[0]
        else : return out

    def mksuperbias(self,ims,display=None,scat=None) :
        """ Driver for superbias combination (no superbias subraction no normalization)
        """
        return self.median(ims,display=display,div=False,scat=scat)

    def mksuperdark(self,ims,superbias=None,display=None,scat=None) :
        """ Driver for superdark combination (no normalization)
        """
        return self.median(ims,superbias=superbias,display=display,div=False,scat=scat)

    def mksuperflat(self,ims,superbias=None,superdark=None,scat=None,display=None) :
        """ Driver for superflat combination (with superbias if specified, normalize to normbox
        """
        return self.median(ims,superbias=superbias,superdark=superdark,normalize=True,scat=scat,display=display)

    def mkspecflat(self,flats,wid=101) :
        """ Spectral flat takes out variation along wavelength direction
        """
        boxcar = Box1DKernel(wid)
        for iflat,flat in enumerate(flats) : 
            nrows=flats[iflat].data.shape[0]
            med = convolve(np.median(flat,axis=0),boxcar,boundary='extend')
            for row in range(flats[iflat].data.shape[0]) :
                flats[iflat].data[row,:] /= med
                flats[iflat].uncertainty.array[row,:] /= med

        return flats


# old combine routine
def combine(ims,norm=False,bias=None,flat=None,trim=False,verbose=False,
            disp=None,min=None,max=None,div=False) :
    """ 
    Combines input list of images (names or numbers) by median, 
    optionally normalizes before combination

    Args:
      ims (list of int or str): list of images to combine

    Keyword args:
      norm (bool) : normalize images before combining? (default=False)
      bias (numpy array) : bias frame to subtract? (default=None)
      flat (numpy array) : flat frame to divide? (default=None)

    Returns:
      median of input data arrays, after reduction as requested
    """
    cube=[]
    for im in ims :
        print('Reading image: ', im)
        h=reduce(im,bias=bias,flat=flat,trim=trim,verbose=verbose) 
        if norm :
            b=det.normbox
            norm=np.median(h.data[b.ymin:b.ymax,b.xmin:b.xmax])
            print('Normalizing image by : ', norm)
            cube.append(h.data/norm)
        else :
            cube.append(h.data)
    print('Combining: ', ims)
    comb = np.median(cube,axis=0)
    if disp is not None :
        for im in ims :
            print(im)
            h=reduce(im,bias=bias,flat=flat,trim=trim,verbose=verbose) 
            if norm :
                b=det.normbox
                norm=np.median(h.data[b.ymin:b.ymax,b.xmin:b.xmax])
                h.data /= norm
            print('Normalizing image by : ', norm)
            if div :
                disp.tv(h.data/comb,min=min,max=max)
            else :
                disp.tv(h.data-comb,min=min,max=max)
        disp.tv(comb,min=min,max=max)
        pdb.set_trace()
    return comb
    
def specflat(flat,rows=True,indiv=False,wid=100) :
    """
    Removes spectral signature from a flat by dividing by smoothed version

    Args:
      flat: input flat fields
   
    Keyword args:
      rows (bool) : specifies if smoothing is along rows (default), otherwise columns
      indiv (bool) : specifies if smoothing is done row-by-row, or single for whole image (default)
      wid (int) : boxcar kernel width (default=100)

    Returns:
      flat with spectral signature removed
    """     
    boxcar = Box1DKernel(wid)
    smooth=flat
    if rows :
        if indiv :
            for row in range(flat.shape[0]) :
                smooth[row,:] /= convolve(flat[row,:],boxcar,boundary='extend')
        else : 
            c=convolve(flat.sum(axis=0),boxcar,boundary='extend')
            for row in range(flat.shape[0]) :
                smooth[row,:] /= c

    else :
        print('smoothing by columns not yet implemented!')
        pdb.set_trace()
      
    return smooth 

class DET() :
    """ 
    Defines detector class 
    """
    def __init__(self,inst=None,gain=0.,rn=0.,biastype=0,biasbox=None,normbox=None,trimbox=None,formstr='{:04d}') :
        self.gain = gain
        self.rn = rn
        self.biastype = biastype
        if biasbox is None : self.biasbox = image.BOX()
        else : self.biasbox = biasbox
        if normbox is None : self.normbox = image.BOX()
        else : self.normbox = normbox
        if trimbox is None : self.trimbox = image.BOX()
        else : self.trimbox = trimbox
        self.formstr = formstr
        if inst == 'ARCES' :
            # APO ARCES
            self.gain=3.8
            self.rn=7
            self.biasbox.set(2052,2057,20,2028)
            self.trimbox.set(200,1850,0,2047)
            self.biasbox = [ image.BOX(xr=[1030,1050],yr=[0,2047]) ,
                             image.BOX(xr=[1030,1050],yr=[0,2047]) ]
            self.trimbox = [ image.BOX(xr=[0,2047],yr=[0,1023]) ,
                             image.BOX(xr=[1030,1050],yr=[0,2047]) ]
        elif inst == 'DIS' :
            # DIS blue
            self.gain=[1.71,1.71]
            self.rn=[3.9,3.9]
            self.biasbox = [ image.BOX(xr=[1030,1050],yr=[0,2047]) ,
                             image.BOX(xr=[1030,1050],yr=[0,2047]) ]
            self.trimbox = [ image.BOX(xr=[0,2047],yr=[0,1023]) ,
                             image.BOX(xr=[1030,1050],yr=[0,2047]) ]
            self.formstr=['{:04d}b','{:04d}r']

def look(tv,pause=True,files=None,list=None,min=None, max=None) :
    """ 
    Displays series of files 
    """
    if files is None and list is None :
        files=glob.glob(indir+'/*.fits*')
    if list is not None :
        f=open(indir+list,'r')
        files=[]
        for line in f :
            files.extend(np.arange(int(line.split()[0]),int(line.split()[1])))
        f.close()
        pdb.set_trace()
 
    for file in files :
       hd=read(file,verbose=True,bias=False)
       disp(tv,hd,min=min,max=max)
       if pause :
           pdb.set_trace()

def getfiles(type,listfile=None,filter=None,verbose=False) :
    """ 
    Get all files of desired type from specified directory. If file is specified, read numbers from that file, else use IMAGETYP card
    """
    if listfile is None :
        list=[]
        if verbose: print('directory: ', indir)
        for file in glob.glob(indir+'/'+root+'*.fits*') :
            head=fits.open(file)[0].header
            #if verbose :
            #   print('file: ', file)
            #   print('IMAGETYPE: ', head['IMAGETYP'])

            try :
                if head['IMAGETYP'] == type :
                    if filter is None or head['FILTER'] == filter :
                        list.append(file)
            except :
                pass
    else :
        list=ascii.read(indir+listfile,Reader=ascii.NoHeader)['col1']
    return list

def disp(tv,hd,min=None,max=None,sky=False) :
    """ 
    Displays HDU or data array on specified ds9/tv device
    """

    if type(tv) is pyds9.DS9 :
        if isinstance(hd, (np.ndarray)) :
            tv.set_np2arr(hd)
            data=hd
        elif isinstance(hd, (astropy.io.fits.hdu.hdulist.HDUList)) :
            tv.set_pyfits(hd)
            data=hd[0].data
        elif isinstance(hd, (astropy.io.fits.hdu.image.PrimaryHDU)) :
            tv.set_np2arr(hd.data)
            data=hd.data
        else :
            print('Unrecognized data type for display: ',type(hd)) 
        if sky :
           skyval = mmm.mmm(data)
           min = skyval[0]-5*skyval[1]
           max = skyval[0]+20*skyval[1]
        if min is not None and max is not None :
            tv.set("scale limits {:5d} {:5d}".format(min,max))
        else :
            tv.set("scale histequ")

    else :
        if isinstance(hd, (astropy.io.fits.hdu.hdulist.HDUList)) :
            tv.tv(hd[0],min=min,max=max)
        else :
            tv.tv(hd,min=min,max=max)

def mkmask(inst=None) :

    if inst == 'ARCES' :
        nrow=2068
        ncol=2128
        badpix = [ image.BOX(yr=[0,2067],xr=[0,230]),      # left side
                   image.BOX(yr=[0,2067],xr=[1900,2127]), # right side
                   image.BOX(yr=[802,2000],xr=[787,787]),
                   image.BOX(yr=[663,2000],xr=[1682,1682]),
                   image.BOX(yr=[219,2067],xr=[101,101]),
                   image.BOX(yr=[1792,1835],xr=[1284,1284]),
                   image.BOX(yr=[1474,2067],xr=[1355,1355]),
                   image.BOX(yr=[1418,1782],xr=[1602,1602]),
                   image.BOX(yr=[1905,1943],xr=[1382,1382]),
                   image.BOX(yr=[1926,1974],xr=[1416,1416]),
                   image.BOX(yr=[1610,1890],xr=[981,981]),
                   image.BOX(yr=[1575,2067],xr=[490,490]),
                   image.BOX(yr=[1710,1722],xr=[568,568]),
                   image.BOX(yr=[1905,1981],xr=[653,654]),
                   image.BOX(yr=[1870,1925],xr=[853,853]) ] 
            

    elif inst == 'DIS' :
        nrow=1078
        ncol=2098
        badpix = [ image.BOX(yr=[474,1077],xr=[803,803]),
                   image.BOX(yr=[0,1077],xr=[1196,1196]),
                   image.BOX(yr=[0,1077],xr=[0,0]) ]

    mask = np.zeros([nrow,ncol],dtype=np.int16)
    for box in badpix :
        box.setval(mask,True)

    hdulist=fits.HDUList()
    hdulist.append(fits.PrimaryHDU(mask))
    hdulist.writeto(inst+'_mask.fits',overwrite=True)

    return mask

def getinput(str) :
    get = input(str)
    if get == 'i' : code.interact(local=globals())
    elif get == 'p' :
        pdb.set_trace()
    return get

class Data(object) :
    """ Experimental data class to cinclude wavelength array
    """
    def __init__(self,data,wave=None) :
        if type(data) is str :
            hdulist=fits.open(data)
            self.meta = hdulist[0].header
            self.attr_list = []
            for i in range(1,len(hdulist) ) :
                try : 
                    attr=hdulist[i].header['ATTRIBUT']
                except KeyError :
                    if i == 1 : attr='data'
                    elif i == 2 : attr='uncertainty'
                    elif i == 3 : attr='mask'
                    elif i == 4 : attr='wave'
                print('attr: {:s}'.format(attr))
                self.attr_list.append(attr)
                setattr(self,attr,hdulist[i].data) 
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
        for attr in self.attr_list :
            hdulist.append(fits.ImageHDU(getattr(self,attr)))
        hdulist.writeto(file,overwrite=overwrite)
