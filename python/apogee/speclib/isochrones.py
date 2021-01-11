import os
import pdb
import numpy as np
from astropy.io import ascii
import astropy.table as table
from tools import plots
import matplotlib.pyplot as plt

def read(infile,columns=None,agerange=[0,20.]) :
    """
    Routine to read a Padova isochrone file and return
    a numpy structured array with the contents

    Args: 
        file name : input data file name
        columns=[list] : list of columns to extract
        age = age : single age to extract

    Returns: 
        structured array with isochrone data
    """
    if os.getenv('ISOCHRONE_DIR') != "" and infile.find('/') < 0 :
        data=ascii.read(os.getenv('ISOCHRONE_DIR')+'/'+infile,
             names=['z','age','mini','mact','logl','logte','logg',
                    'mbol','u','b','v','r','i','j','h','k','intimf','stage'])
    else :
        data=ascii.read(infile,
             names=['z','age','mini','mact','logl','logte','logg',
                    'mbol','u','b','v','r','i','j','h','k','intimf','stage'])

    # add some derived columns
    data.add_column(table.column.Column(name='feh',data=np.log10(data['z']/0.0152)))
    data.add_column(table.column.Column(name='teff',data=10.**(data['logte'])))
    data.add_column(table.column.Column(name='jk',data=data['j']-data['k']))

    # select ages within specified age range
    gd = np.where((data['age'] >=agerange[0]) & (data['age'] <= agerange[1]))
    data=data[gd]

    # default columns
    if columns is None:
        # can set default columns here, or keep all quantitites
        #data.keep_columns(['Z','age','logte','logl','intimf','stage'])
        pass
    # option to extract specified columns
    else:
        data.keep_columns(columns)

    return data

def plot(ax,iso,x,y,xr=None,yr=None,color=None,dx=0.,dy=0.,isoadj=False,alpha=None) :
    ''' plotting routine that handles tip of RGB '''
    mdiff = iso['mini'][0:-1]-iso['mini'][1:]
    j=np.where(abs(mdiff) < 1.e-8)[0]
    if len(j) > 0 :
        j=j[0]
        if isoadj :
            a=116.
            b=155.
            c=-22.
            dx = a + b*iso['feh'][0:j] + c*iso['logg'][0:j]
        if x == 'te' : line = plots.plotl(ax,10.**(iso['logte'][0:j]+dx),iso[y][0:j]+dy,xr=xr,yr=yr,color=color,alpha=alpha)
        else : line = plots.plotl(ax,iso[x][0:j]+dx,iso[y][0:j]+dy,xr=xr,yr=yr,color=color,alpha=alpha)
        if x == 'te' : plots.plotl(ax,10.**iso['logte'][j+1:],iso[y][j+1:],color=line[0].get_color(),alpha=alpha)
        else : plots.plotl(ax,iso[x][j+1:],iso[y][j+1:],color=line[0].get_color(),alpha=alpha)
    else :
        if x == 'te' : line = plots.plotl(ax,10.**(iso['logte']+dx),iso[y]+dy,xr=xr,yr=yr,color=color,alpha=alpha)
        else : line = plots.plotl(ax,iso[x]+dx,iso[y]+dy,xr=xr,yr=yr,color=color,alpha=alpha)
    plt.draw()

