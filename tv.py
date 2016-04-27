import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
from astropy.wcs import wcs
import cmap
import mmm
import autopy
import pdb
 
class TV:
 
    """
    A "TV" figure
    
    Usage: import tv
           tv=TV()  to set up a new TV object (display window)
    """
 
    def __init__(self, img=np.zeros([5,5]), fig=None):
    
        """
        Initialize TV object
        """
    
        # create new figure,set title, margins, and facecolor
        tv = plt.figure(figsize=(8,8.5))
        self.fig = tv
        tv.canvas.set_window_title('Image display window')
        tv.set_facecolor('darkred')
        #ax = plt.axes()
        #ax = tv.add_subplot(111)
        rect = 0., 0.05, 1., 0.95
        ax = tv.add_axes(rect)
        self.ax = ax
        ax.axis('off')
        self.axis = False

        # set up initial img and header lists
        self.current = -1
        self.images = 0
 
        # initialize rolling buffers
        self.img = None
        self.imglist = [None, None, None, None]
        self.hdr = None
        self.hdrlist = [None, None, None, None]
        self.scale = [0.,1.]
        self.scalelist = [[0.,1.],[0.,1.],[0.,1.],[0.,1.]]
        self.cmap = 'Greys_r'
        self.axlist = [None, None, None, None]

        # set up colorbar
        self.cb = None
        self.cblist = [None, None, None, None]
        rect = 0.00, 0.03, 1., 0.06
        self.cb_ax = tv.add_axes(rect)
        #tv.subplots_adjust(left=-0.15,right=1.15,bottom=-0.10,top=1.00)
        self.bottom = 0.
        self.top = 1.

        # function to show image values, etc.
        def format_coord(x, y):
            x = int(x + 0.5)
            y = int(y + 0.5)

            try:
                self.img
                try:
                    hdr=self.hdr
                    mywcs=wcs.WCS(hdr)
                    pixcrd = np.array([[x,y]])
                    world=mywcs.wcs_pix2world(pixcrd,1)
                    try:
                       object=self.hdr['object']
                    except:
                       object=None
                    return "[x,y]=[%4d, %4d] val=%8.1f   [%s %s]=[%10.6f,%10.6f]   OBJECT: %s" % (x,y, self.img[y, x], mywcs.wcs.ctype[0],mywcs.wcs.ctype[1],world[0,0], world[0,1], object)
		except:
                    mywcs=None
                try:
                    return "[x,y]\n [%4i, %4i] val=%8.1f OBJECT: %s" % (x,y, self.img[y, x], object)
                except IndexError:
                    return ""
            except:
                return " [%4i, %4i]" % (x, y)

        # set this routine up for format
        ax.format_coord = format_coord

        #event handling 
        self.event = None
        # turn off default key bindings
        tv.canvas.mpl_disconnect(tv.canvas.manager.key_press_handler_id)
        # set up our event handler
        self.cid = tv.canvas.mpl_connect('key_press_event', self.__onEvent)
        self.cid2 = tv.canvas.mpl_connect('button_press_event', self.__onEvent)
        self.cid3 = tv.canvas.mpl_connect('button_release_event', self.__onEvent)
        self.cid4 = tv.canvas.mpl_connect('motion_notify_event', self.__onEvent)
        self.button = False
        self.blocking = 0

    def __onEvent(self, event):
        """
        Handler for all trapped events. 
        
        Args:
          event -- a KeyEvent
        """
        self.event = event
        subPlotNr = self.__getSubPlotNr(event)        

        if event.name == 'key_press_event' :
            # keypress events: '-', '+/=', 'r'
            self.key = event.key

            # function for computing "scale"
            def scale() :
                # compute screen pixels per image pixel
                p1 = self.ax.transData.transform((0.,0.))
                p2 = self.ax.transData.transform((100.,100.))
                return (p2[1]-p1[1])/100., (p2[0]-p1[0])/100.

            if event.key == '-' or event.key == '+' or event.key == '=':
                if event.key == '-' :
                    self.current = (self.current-1) % self.images
                elif event.key == '+' or event.key == '=':
                    self.current = (self.current+1) % self.images
                self.img = self.imglist[self.current]
                self.hdr = self.hdrlist[self.current]
                self.scale = self.scalelist[self.current]
                for i in range(self.images) :
                    if i == self.current :
                        self.axlist[i].set_visible(True)
                    else :
                        self.axlist[i].set_visible(False)
                self.aximage=self.axlist[self.current]
                self.cb=self.cblist[self.current]
                #self.cb.ax.clear()
                #self.cb = self.fig.colorbar(self.aximage,cax=self.cb.ax,orientation='horizontal')
                #self.cb = self.ax.get_figure().colorbar(self.aximage,cax=self.cb_ax,orientation='horizontal')
                #cm=cmap.remap(self.cmap,self.bottom,self.top)
                #self.aximage.set_cmap(cm)
                plt.draw()
                x,y= autopy.mouse.get_pos()
                autopy.mouse.move(x,y)

            elif (event.key == 'p' or event.key == 'v') and subPlotNr == 0 :
                n=7
                xdata=int(round(event.xdata))
                ydata=int(round(event.ydata))
                if event.key == 'p' :
                    py, px = np.unravel_index(np.argmax(self.img[ydata-n:ydata+n,xdata-n:xdata+n]),
                                              self.img[ydata-n:ydata+n,xdata-n:xdata+n].shape)
                else :
                    py, px = np.unravel_index(np.argmin(self.img[ydata-n:ydata+n,xdata-n:xdata+n]),
                                              self.img[ydata-n:ydata+n,xdata-n:xdata+n].shape)
                px-=n
                py-=n
                x,y= autopy.mouse.get_pos()
                xs,ys = scale()
                autopy.mouse.move(int(x+px*xs),int(y-py*ys))

            elif event.key == 'r' and subPlotNr == 0 :
                dim=np.shape(self.img)
                size=np.max([dim[0],dim[1]])
                self.ax.set_xlim(dim[1]/2.-size/2.,dim[1]/2.+size/2.)
                self.ax.set_ylim(dim[0]/2.-size/2.,dim[0]/2.+size/2.)
                #self.ax.set_xlim(-0.5,dim[1]-0.5)
                #self.ax.set_ylim(-0.5,dim[0]-0.5)
                plt.draw()

            elif event.key == 'r' and subPlotNr == 1 :
                self.bottom=0.
                self.top=1.
                cm=cmap.remap(self.cmap,self.bottom,self.top)
                self.aximage.set_cmap(cm)
                plt.draw()

            elif event.key == 'left' and subPlotNr == 0 :
                xs,ys = scale()
                x,y= autopy.mouse.get_pos()
                if xs < 1. :
                    autopy.mouse.move(x-1,y)
                else :
                    autopy.mouse.move(int(x-xs),y)

            elif event.key == 'right' and subPlotNr == 0 :
                xs,ys = scale()
                x,y= autopy.mouse.get_pos()
                if xs < 1. :
                    autopy.mouse.move(x+1,y)
                else :
                    autopy.mouse.move(int(x+xs),y)

            elif event.key == 'up' and subPlotNr == 0 :
                xs,ys = scale()
                x,y= autopy.mouse.get_pos()
                if ys < 1. :
                    autopy.mouse.move(x,y-1)
                else :
                    autopy.mouse.move(x,int(y-ys))

            elif event.key == 'down' and subPlotNr == 0 :
                xs,ys = scale()
                x,y = autopy.mouse.get_pos()
                if ys < 1. :
                    autopy.mouse.move(x,y+1)
                else :
                    autopy.mouse.move(x,int(y+ys))

            elif event.key == 'a' and subPlotNr == 0 :
                if self.axis :
                    rect = 0., 0.05, 1., 0.95
                    self.ax.axis('off')
                else :
                    rect = 0.05, 0.15, 0.95, 0.85
                    self.ax.axis('on')
                self.ax.set_position(rect)
                self.axis = not self.axis
                plt.draw()

            elif event.key == 'h' or event.key == '?' :
                print 'Asynchronous commands: '
                print 'Image window: '
                print '  mouse:'
                print '    left mouse  : zoom in, centered on cursoe'
                print '    center mouse: zoom out, centered on cursoe'
                print '    right mouse : pan, center to cursor'
                print '  keys:'
                print '    r           : redraw at default zoom'
                print '    +/=         : toggle to next image in stack'
                print '    -           : toggle to previous image in stack'
                print '    arrow keys  : move single image pixels'
                print '    a           : toggle axes on/off'
                print '    h/?         : print this help'

            if self.blocking == 1 : self.__stopBlock()

        elif event.name == 'button_press_event' :
            if subPlotNr == 0 :
                # button press in image window to zoom/pan
                xlim = self.ax.get_xlim()
                ylim = self.ax.get_ylim()
                if event.button == 1 :
                    # zoom in
                    xsize = ( xlim[1]-xlim[0] )/ 2.
                    ysize = ( ylim[1]-ylim[0] )/ 2.
                elif event.button == 2 :
                    # zoom out
                    xsize = ( xlim[1]-xlim[0] )* 2.
                    ysize = ( ylim[1]-ylim[0] )* 2.
                else :
                    # pan
                    xsize = xlim[1]-xlim[0]
                    ysize = ylim[1]-ylim[0]
                size=max([xsize,ysize])
                self.ax.set_xlim(event.xdata-size/2.,event.xdata+size/2.)
                self.ax.set_ylim(event.ydata-size/2.,event.ydata+size/2.)
                plt.draw()
            elif subPlotNr == 1 :
                # flag button press in colorbar
                self.button = True
                self.xstart = event.xdata

        elif event.name == 'button_release_event' :
            self.button = False

        elif event.name == 'motion_notify_event'  and self.button :
            # if motion in colorbar with key pressed in colorbar, adjust colorbar
            if subPlotNr == 1 :
                xend = event.xdata
                if event.button == 2 :
                    diff = (xend - self.xstart)
                    self.top = self.top + diff
                    self.bottom = self.bottom + diff
                    self.xstart = xend
                else :
                    if self.xstart > 0.5 :
                      if xend > self.bottom :
                          self.top = xend
                      else :
                          self.top = self.bottom
                    else :
                      if xend < self.top :
                          self.bottom = xend
                      else :
                          self.bottom = self.top
                cm=cmap.remap(self.cmap,self.bottom,self.top)
                self.aximage.set_cmap(cm)
                plt.draw()

        
    def __getSubPlotNr(self, event):
    
        """
        Get the nr of the subplot that has been clicked
        
        Args::
          event -- an event
        
        Returns:
          A number or None if no subplot has been clicked
        """
    
        i = 0
        axisNr = None
        for axis in self.fig.axes:
            if axis == event.inaxes:
                axisNr = i        
                break
            i += 1
        return axisNr
        
    def __stopBlock(self) :
        """
        stops blocking for keypress event
        """
        if self.blocking == 1 :
            self.fig.canvas.stop_event_loop()

    def __startBlock(self) :
        """
        starts blocking for keypress event
        """
        self.blocking = 1
        self.fig.canvas.start_event_loop(-1)

    def tv(self,img,min=None,max=None,cmap=None) :
        """
        main display routine: displays image with optional scaling

        Args:
          img: a numpy array OR a fits HDU

        Keyword args:
          min=, max= : optional scaling arguments
        """

        # load data array depending on input type
        if isinstance(img, (np.ndarray)) :
            data = img
        elif isinstance(img.data, (np.ndarray)) :
            data = img.data
        else :
            print 'input must be numpy array or have data attribute that is'
            return

        # set figure and axes
        plt.figure(self.fig.number)
        plt.axes(self.ax)
        #self.fig.clf()

        # load new image data onto rolling stack
        current= (self.current+1) % 4
        self.images += 1
        if self.images > 4 : self.images = 4
        self.current = current
        self.imglist.pop(current)
        self.imglist.insert(current,data)
        self.img = data

        # save the header if we have one
        if hasattr(img,'header') :
            self.hdrlist.pop(current)
            self.hdrlist.insert(current,img.header)
            self.hdr=img.header
        else :
            self.hdr=None
      
        # get autodisplay parameters if needed, and save display params
        if min is None : 
           min = 0.
        if max is None : 
           sky = mmm.mmm(data)
           min = sky[0]-5*sky[1]
           max = sky[0]+20*sky[1]
        self.scale = [min,max]
        self.scalelist.pop(current)
        self.scalelist.insert(current,self.scale)

        if cmap != None :
           self.cmap = cmap
 
        # display image and new colorbar 
        dim=np.shape(self.img)
        size=np.max([dim[0],dim[1]])
        self.ax.set_xlim(dim[1]/2.-size/2.,dim[1]/2.+size/2.)
        self.ax.set_ylim(dim[0]/2.-size/2.,dim[0]/2.+size/2.)
        #self.ax.set_xlim(-0.5,dim[1]-0.5)
        #self.ax.set_ylim(-0.5,dim[0]-0.5)
        self.aximage = self.ax.imshow(data,vmin=min,vmax=max,cmap=self.cmap,interpolation='nearest')
        old=self.axlist.pop(current)
        # if we had a previous image, reload the data with a single value
        # so we don't continually accumulate memory (matplotlib doesn't
        # describe how memory can be released
        z=np.zeros([1,1])
        if old is not None : old.set_data(z)
        self.axlist.insert(current,self.aximage)
        if self.cb is None :
            #self.cb = self.fig.colorbar(self.aximage,orientation='horizontal',shrink=0.7,pad=0)
            self.cb = self.fig.colorbar(self.aximage,cax=self.cb_ax,orientation='horizontal')
            #plt.subplots_adjust(left=-0.15,right=1.15,bottom=-0.10,top=1.00)
        else :
            self.cb.ax.clear()
            self.cb = self.fig.colorbar(self.aximage,cax=self.cb.ax,orientation='horizontal')
        self.cblist.pop(current)
        self.cblist.insert(current,self.cb)

        plt.draw()
        # instead of redraw color, could replace data, but not if sizes change?
        # img.set_data()
        # img.changed()
        # plt.draw()

            
    def tvbox(self,x,y,size=None,color=None) :
        """
        displays a patch (box by default) on an image

        Args:
          x,y : center position of patch

        Keyword args :
          size= :  patch size
          color= :  patch color
        """
        plt.figure(self.fig.number)
        if size == None : size = 3
        if color == None : color = 'm'
        self.ax.add_patch(patches.Rectangle((x-size/2-0.5,y-size/2-0.5),size,size,fill=False,color=color))
        plt.draw()

    def tvclear(self) :
        """
        clears patches from image
        """
        plt.figure(self.fig.number)
        for i in range(len(self.ax.patches)) : self.ax.patches[0].remove()
        plt.draw()

    def tvmark(self) :
        """
        Blocking input: waits for key press in display and returns key 
        that was pressed and data pixel location of the keypress

        Args:
          none

        Returns:
          key pressed, x data position, y data position
        """
        self.__startBlock()
        return self.event.key, self.event.xdata, self.event.ydata

    def fill(self) :
        y,x=np.mgrid[0:100,0:100]
        self.tv(x)
        self.tv(y)
        self.tv(x+y)
        self.tv(x-y)
