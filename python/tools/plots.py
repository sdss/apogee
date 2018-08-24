#import matplotlib.colors as colors
import matplotlib.pyplot as plt
from scipy import spatial
from scipy.stats import gaussian_kde
import struct
import numpy as np
import sys
import pdb

from matplotlib.widgets import Lasso
from matplotlib.collections import RegularPolyCollection
from matplotlib import colors as mcolors, path

# default values for event handling
_index = 0
_data = None
_button = None
_id_cols = ['APOGEE_ID']
_data_x = None
_data_y = None
_block = 0

def event(fig) :
    '''
    Define event handler, on a key press, will set _button (key pressed) and _index (index of nearest point), and
    if _data is not None, will list _id_cols from the structure _data[_index]
    '''
    global _block
    def onpress(event) :
        global _index, _block, _x, _y, _button
        _button = event.key
        inv = event.inaxes.transData.inverted()
        _x,_y = inv.transform((event.x,event.y))
        print(_x, _y)
        #A[spatial.KDTree(A).query([event.x,event.y])[1]]
        #distance,index = spatial.KDTree(A).query([event.x,event.y])
        if _data_x is not None and _data_y is not None :
            print('Transform', len(_data_x))
            A = event.inaxes.transData.transform(zip(_data_x,_data_y))
            print('KDTree')
            tree=spatial.KDTree(A)
            print('query')
            distance,index = tree.query([event.x,event.y])
            _index = [index]
            print('_index: ',_index,_data_x[_index],_data_y[_index])
            if _data is not None :
                struct.list(_data,ind=_index,cols=_id_cols)
        if _block == 1 :
            _block = 0
            fig.canvas.stop_event_loop()
    cid = fig.canvas.mpl_connect('key_press_event',onpress)

def mark(fig) :
    global _block, _x, _y, _button
    _block = 1
    fig.canvas.start_event_loop(-1)
    return _x, _y, _button


def plotc(ax,x,y,z,yerr=None,xr=None,yr=None,zr=None,size=5,cmap='rainbow',colorbar=False,xt=None,yt=None,zt=None,label=None,linewidth='0',marker='o',draw=True,orientation='vertical',labelcolor='k',tit=None,nxtick=None,nytick=None,rasterized=None) :
    """
    Plots a scatter plot with point color-coded by z data

    Args:
      ax (axis)  : existing axes
      x (float)  : x values
      y (float)  : y values
      z (float)  : z values, or a single color

    Keyword args:
      xr : x range  (default=None)
      yr : y range (default=None)
      zr : z range (default=None)
      xt : x axis title (default=None) (default=None)
      yt : y axis title (default=None)
      zt : z axis title (default=None)
      marker : marker type (default='o')
      cmap : colormap (default='rainbow')
      size : point size(s) in points (default=5)
      linewidth : linewidth (default='0')
      colorbar (bool) : draw a colorbar? (default=False)
      orientation : colorbar orientation (default='vertical')
      label [x,y,text] : label plot with text at position x,y in relative (0:1,0:1) coordinates
      labelcolor : color for label

    Returns:
      aximage

    """
    set_limits_ticks(ax,xr,yr,nxtick,nytick)
    if xt is not None : ax.set_xlabel(xt) 
    if yt is not None : ax.set_ylabel(yt)
    if tit is not None : ax.set_title(tit)
    if zr is None :
        scat=ax.scatter(x,y,c=z,s=size,cmap=cmap,linewidth=linewidth,marker=marker,rasterized=rasterized)
    else :
        scat=ax.scatter(x,y,c=z,vmin=zr[0],vmax=zr[1],s=size,cmap=cmap,linewidth=linewidth,marker=marker,rasterized=rasterized)

    if yerr is not None :
        ax.errorbar(x,y,yerr=yerr,fmt='none',capsize=0,ecolor='k')
    if label is not None :
        ax.text(label[0],label[1],label[2],transform=ax.transAxes,color=labelcolor) 
    if colorbar :
        cb=plt.colorbar(scat,ax=ax,orientation=orientation)
        cb.ax.set_ylabel(zt)
    if draw : plt.draw()
    _data_x = x[np.isfinite(x)]
    _data_y = y[np.isfinite(y)]
    return scat

def set_limits_ticks(ax,xr,yr,nxtick=None,nytick=None) :
    if nxtick is not None:
        if xr is not None : ax.set_xlim(xr[0],xr[1])
        ax.xaxis.set_ticks(np.linspace(ax.get_xlim()[0],ax.get_xlim()[1],nxtick)[1:-1])
    else :
        if xr is not None : ax.set_xlim(xr[0]+0.01*(xr[1]-xr[0]),xr[1]-0.01*(xr[1]-xr[0]))
    if nytick is not None:
        if yr is not None : ax.set_ylim(yr[0],yr[1])
        ax.yaxis.set_ticks(np.linspace(ax.get_ylim()[0],ax.get_ylim()[1],nytick)[1:-1])
    else :
        if yr is not None : ax.set_ylim(yr[0]+0.01*(yr[1]-yr[0]),yr[1]-0.01*(yr[1]-yr[0]))

def plotc_append(ax,x,y,z,size=25,linewidth=1,marker='o',facecolor='none',draw=True) :
    '''
    Adds points to a plot
    '''
    scat=ax.scatter(x,y,c=z,s=size,linewidth=linewidth,marker=marker,facecolor=facecolor)
    if draw : plt.draw()

def plotrow(ax,img,r,norm=True,draw=True) :
    '''
    Plots a row of an input image

    Args:
      ax (axis)    : existing axes
      img (float)  : image to plot
      r (int)      : row or rows ([rmin,rmax]) to plot

    Keyword args:
      norm         : for multiple row plots, average instead of sum (default=True)
    '''
    ax.set_xlim(xr[0],xr[1])
    ax.set_ylim(yr[0],yr[1])
    if len(r) == 1 :
        ax.plotl(img[r,:])
    elif len(r) == 2 :
        if norm :
            ax.plotl(np.average(img[r[0]:r[1],:],axis=1))
        else :
            ax.plotl(np.sum(img[r[0]:r[1],:],axis=1))
    if draw : plt.draw()

def plotp(ax,x,y,z=None,typeref=None,types=None,xr=None,yr=None,zr=None,marker='o',size=5,linewidth='0',color='r',facecolors=None,xt=None,yt=None,draw=True,xerr=None,yerr=None,label=None,text=None,labelcolor='k',linewidths=None,nxtick=None,nytick=None,tit=None,contour=None,levels=None,alpha=None) :
    '''
    Plot points, optionally with a series of different markers/sizes keyed to z data

    Args:
        ax : axes to plot in
        x : x data
        y : y data

    Keyword args:
        z=  : specifies data to be used to color points if using types to specify sizes, markers (default=None)
        typeref= : specfies array to be used to determine groupings
        types= : array of different types (from typeref) to plot with specified sizes, markers, colors
        size= : array of different sizes to plot for different types, or single size
        marker= : array of different markers to plot for different types, or single marker
        color= : array of different colors markers to plot for different types, or single color
        xr : x limits (default=None)
        yr : y limits (default=None)
        xt : x title (default=None)
        yt : y title (default=None)
        label : label for legend
        text=[x,y,text] : put text at (x,y) relative coords
        labelcolor=  : color for label
       
    '''
    set_limits_ticks(ax,xr,yr,nxtick,nytick)
    if xt is not None : ax.set_xlabel(xt) 
    if yt is not None : ax.set_ylabel(yt)
    if tit is not None : ax.set_title(tit)
    if facecolors is None: facecolors=color

    if typeref is not None and types is not None :
        # Make sure types, sizes, markers are all lists
        try :
            test = len(types)
        except :
            types = [types]
        try :
            test = len(size)
        except :
            size = [size]
        try :
            test = len(marker)
        except :
            marker = [marker]
        try :
            test = len(color)
        except :
            color = [color]

        # loop through the types
        for i in range(len(types)) :
            gd = np.where(typeref == types[i])[0]
            sz= size[i] if (len(size) > 1)  else size[0]
            mark=marker[i] if (len(marker) > 1) else marker[0]
            col=color[i] if (len(color) > 1) else color[0]
            if facecolors is 'none' : facecol = 'none'
            else : facecol = col
            if z is not None :
                if zr is not None :
                    ax.scatter(x[gd],y[gd],c=z[gd],s=sz,marker=mark,vmin=zr[0],vmax=zr[1],linewidths=linewidths)
                else :
                    ax.scatter(x[gd],y[gd],c=z[gd],s=sz,marker=mark,linewidths=linewidths)
            else :
                ax.scatter(x[gd],y[gd],s=sz,marker=mark,facecolors=facecol,edgecolors=col,linewidths=linewidths)
            if yerr is not None :
                ax.errorbar(x[gd],y[gd],marker=mark,yerr=yerr[gd],fmt='none',capsize=0,ecolor=col)
    elif contour is not None:
        if contour <= 0 :
            gd = np.where((x > xr[0]) & (x < xr[1]) & (y>yr[0]) & (y<yr[1]) )[0]
            data = np.vstack([x[gd],y[gd]])
            kde = gaussian_kde(data,bw_method=abs(contour))
            xgrid = np.linspace(xr[0], xr[1], 40)
            ygrid = np.linspace(yr[0], yr[1], 40)
            Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
            Z = kde.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))
            # Plot the result as an image
            ax.imshow(np.log(Z).reshape(Xgrid.shape), origin='lower', aspect='auto',vmin=0,
                       extent=[xr[0], xr[1], yr[0], yr[1]], cmap='Blues')
        else :
            im = np.histogram2d(y,x,range=[yr,xr],bins=20)
            if levels is None:
                levels=np.linspace(1.,im[0].max(),contour)
            ax.contour((im[2][0:-1]+im[2][1:])/2.,(im[1][0:-1]+im[1][1:])/2.,im[0],
               colors=color,levels=levels,alpha=alpha)
    else :
        ax.scatter(x,y,marker=marker,s=size,linewidth=linewidth,facecolors=facecolors,edgecolors=color,linewidths=linewidths,alpha=alpha,label=label)
        _data_x = x[np.isfinite(x)]
        _data_y = y[np.isfinite(y)]
        if xerr is not None or yerr is not None :
            ax.errorbar(x,y,marker=marker,xerr=xerr,yerr=yerr,fmt='none',capsize=0,ecolor=color)

    if text is not None :
        if labelcolor is 'k' and color is not None : labelcolor=color
        ax.text(text[0],text[1],text[2],transform=ax.transAxes,color=labelcolor)

    if draw : plt.draw()


def plotl(ax,x,y,xr=None,yr=None,color=None,xt=None,yt=None,draw=True,label=None,ls=None,semilogy=False,linewidth=1.,tit=None,nxtick=None,nytick=None) :
    '''
    Plot connected points
    '''
    set_limits_ticks(ax,xr,yr,nxtick,nytick)
    if xt is not None : ax.set_xlabel(xt) 
    if yt is not None : ax.set_ylabel(yt)
    if tit is not None : ax.set_title(tit)
    if ls is None : ls='-'
    if semilogy :
        line = ax.semilogy(x,y,color=color,label=label,ls=ls,linewidth=linewidth)
    else :
        line = ax.plot(x,y,color=color,label=label,ls=ls,linewidth=linewidth)
    if draw : plt.draw()
    return line
    
def ax(subplot=111) :
    '''
    Return axes object for a new figure and desired subplots

    Keyword args :

    subplot  : matplotlib subplot specification
    '''
    fig=plt.figure()
    return fig.add_subplot(subplot)

def multi(nx,ny,figsize=None,hspace=1,wspace=1,sharex=False,sharey=False,squeeze=True,xtickrot=None) :
    '''
    Returns figure and axes array for grid of nx by ny plots, suppressing appropriate axes if requested by hspace and wspace

    Args:
       nx : number of plots in horizontal direction
       ny : number of plots in vertical direction

    Keyword args:
       figsize  : specifies figure size
       hspace   : space (0.-1.) between vertical plots (height)
       wspace   : space (0.-1.) between horizont plots (width)
       sharex   : force subplots to have same x-axes
       sharey   : force subplots to have same y-axes
    '''
    fig,ax = plt.subplots(ny,nx,figsize=figsize,sharex=sharex,sharey=sharey,squeeze=squeeze)
    fig.subplots_adjust(hspace=hspace,wspace=wspace)
    if (hspace < 0.01) & (ny>1):
        # if we are vertical stacking, turn off xticks for all except bottom
        if nx == 1 :
            ticklabels = ax[0].get_xticklabels()
            for i in range(1,ny-1) : 
                ticklabels = ticklabels + ax[i].get_xticklabels()
        else :
            ticklabels = ax[0,0].get_xticklabels()
            for i in range(nx) :
                for j in range(0,ny-1) : 
                    ticklabels = ticklabels + ax[j,i].get_xticklabels()
        plt.setp(ticklabels, visible=False)
    if (wspace < 0.01) & (nx> 1):
        # if we are horizontal stacking, turn off yticks for all except left
        if ny == 1 :
            ticklabels = ax[1].get_yticklabels()
            for i in range(2,nx) : 
                ticklabels = ticklabels + ax[i].get_yticklabels()
        else :
            ticklabels = ax[0,1].get_yticklabels()
            for i in range(1,nx) :
                for j in range(ny) : 
                    ticklabels = ticklabels + ax[j,i].get_yticklabels()
        plt.setp(ticklabels, visible=False)

    if xtickrot is not None:
      for i in range(nx) :
        for j in range(0,ny) :
          if nx == 1 and ny == 1:
              print('setting rotation')
              for tick in ax.get_xticklabels(): tick.set_rotation( xtickrot )
          elif nx>1 and ny == 1:
              for tick in ax[i].get_xticklabels(): tick.set_rotation( xtickrot )
          elif ny>1 and nx == 1:
              for tick in ax[j].get_xticklabels(): tick.set_rotation( xtickrot )
          else :
              for tick in ax[j,i].get_xticklabels(): tick.set_rotation( xtickrot )
      fig.subplots_adjust(bottom=0.2)
    return fig,ax


def close() :
    '''
    Close open plots windows
    '''
    plt.close('all')

class LassoManager(object):
    ''' Simple Lasso manager to allow user to lasso points and get indices

        Adapted from version from Google search on matplotlib lasso
    '''
    def __init__(self, ax, x, y):
        '''  Initialize manager with axes, x, and y data arrays
        '''
        self.axes = ax
        self.canvas = ax.figure.canvas
        self.Nxy = len(x)

        self.xys=[]
        for i in range(len(x)) :
            self.xys.append((x[i],y[i]))
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)

    def callback(self, verts):
        ''' Once a Lasso is marked with mouse, get points within the path
        '''
        p = path.Path(verts)
        self.ind = p.contains_points(self.xys)
        self.canvas.widgetlock.release(self.lasso)
        del self.lasso

    def onpress(self, event):
        '''called on button_press_event, runs Lasso with callback
        '''
        if self.canvas.widgetlock.locked():
            return
        if event.inaxes is None:
            return
        self.lasso = Lasso(event.inaxes,
                           (event.xdata, event.ydata),
                           self.callback)
        # acquire a lock on the widget drawing
        self.canvas.widgetlock(self.lasso)
