import matplotlib.colors as colors
import matplotlib.pyplot as plt
from scipy import spatial
import struct
import numpy as np
import sys
import pdb

_index = 0
_data = None
_button = None
_id_cols = ['APOGEE_ID']
_data_x = None
_data_y = None

def event(fig) :
    def onpress(event) :
        global _index
        _button = event.key
        A = event.inaxes.transData.transform(zip(_data_x,_data_y))
        #A[spatial.KDTree(A).query([event.x,event.y])[1]]
        #distance,index = spatial.KDTree(A).query([event.x,event.y])
        tree=spatial.KDTree(A)
        distance,index = tree.query([event.x,event.y])
        #print distance, index
        _index = [index]
        if _data is not None :
            struct.list(_data,ind=_index,cols=_id_cols)
    cid = fig.canvas.mpl_connect('key_press_event',onpress)

def plotc(ax,x,y,z,xr=None,yr=None,zr=None,size=5,cmap='rainbow',colorbar=False,xt=None,yt=None,zt=None,label=None,linewidth='0',marker='o',draw=True) :
    """
    Plots a scatter plot with point colors

    Args:
      ax (axis)  : existing axes
      x (float)  : x values
      y (float)  : y values
      z (float)  : z values, or a single color

    Keyword args:
      xr : x range (optional)
      yr : y range (optional)
      zr : z range (optional)
      xt : x axis title (optional)
      yt : y axis title (optional)
      zt : z axis title (optional)
      marker : marker type
      cmap : colormap
      size : point size(s) in points (optional)

    Returns:

    """
    if xr is not None : ax.set_xlim(xr[0],xr[1])
    if yr is not None : ax.set_ylim(yr[0],yr[1])
    if xt is not None : ax.set_xlabel(xt) 
    if yt is not None : ax.set_ylabel(yt)
    print 'marker: ', marker, type(marker),type('o')
    if zr is None :
        scat=ax.scatter(x,y,c=z,s=size,cmap=cmap,linewidth=linewidth,marker=marker)
    else :
        scat=ax.scatter(x,y,c=z,vmin=zr[0],vmax=zr[1],s=size,cmap=cmap,linewidth=linewidth,marker=marker)
    if label is not None :
        ax.text(label[0],label[1],label[2]) 
    if colorbar :
        cb=plt.colorbar(scat)
        cb.ax.set_ylabel(zt)
    if draw : plt.draw()

def plotc_append(ax,x,y,z,size=25,linewidth=1,marker='o',facecolor='none',draw=True) :
    scat=ax.scatter(x,y,c=z,s=size,linewidth=linewidth,marker=marker,facecolor=facecolor)
    if draw : plt.draw()


def plotrow(ax,img,r,norm=True,draw=True) :
    ax.set_xlim(xr[0],xr[1])
    ax.set_ylim(yr[0],yr[1])
    if len(r) == 1 :
        ax.plotl(img[r,:])
    elif len(r) == 2 :
        ax.plotl(np.average(img[r[0]:r[1],:],axis=1))
    if draw : plt.draw()

def plotp(ax,x,y,xr=None,yr=None,marker='o',size=5,linewidth='0',color='r',facecolors=None,xt=None,yt=None,draw=True,xerr=None,yerr=None) :
    try: ax.set_xlim(xr[0],xr[1])
    except : pass
    try : ax.set_ylim(yr[0],yr[1])
    except : pass
    if xt is not None : ax.set_xlabel(xt) 
    if yt is not None : ax.set_ylabel(yt)
    #ax.plot(x,y,sym)
    if facecolors is None: facecolors=color
    print 'color: ', color
    print 'marker: ', marker
    print 'facecolors: ', facecolors
    print 'edgecolors: ', color
    ax.scatter(x,y,marker=marker,s=size,linewidth=linewidth,facecolors=facecolors,edgecolors=color)
    if xerr is not None or yerr is not None :
      ax.errorbar(x,y,marker=marker,xerr=xerr,yerr=yerr,fmt='none',capsize=0,ecolor=color)
    if draw : plt.draw()

def plotl(ax,x,y,xr=None,yr=None,color=None,xt=None,yt=None,draw=True) :
    try: ax.set_xlim(xr[0],xr[1])
    except : pass
    try : ax.set_ylim(yr[0],yr[1])
    except : pass
    if xt is not None : ax.set_xlabel(xt) 
    if yt is not None : ax.set_ylabel(yt)
    line = ax.plot(x,y,color=color)
    if draw : plt.draw()
    return line
    

