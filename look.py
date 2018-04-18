import matplotlib
matplotlib.use('Agg')
import glob
from astropy.io import ascii
from astropy.io import fits
from holtz.tools import plots
from holtz.tools import html
import matplotlib.pyplot as plt
import subprocess
import argparse
import multiprocessing
import logging
import pdb
import numpy as np

def apogee1m(root) :

    dir='/home/1m/'+root+'/'
    files=glob.glob(dir+'apogee_*.dat')
    figs=[]
    ylab=[]
    for file in files :
        print(file)
        fp = open(file,"r")
        for line in fp :
            print(line)
            obj=line.split()[0]
            #makemovie(root,obj)
            p=multiprocessing.Process(target=makemovie,args=(root,obj,))
            p.start()
            ylab.append(obj)
            col=[]
            col.append(obj+'.jpg')
            col.append(obj+'.gif')
            figs.append(col)
    html.htmltab(figs,file=root+'sum.html',ytitle=ylab)

def makemovie(root,obj,min=1500,max=2000) :

    dir='/home/1m/'+root+'/'
    logging.debug('makemovie '+dir+obj)
    try: 
        seqs=ascii.read(dir+obj+'.dat')
        gifs=['convert','-loop','0']
    except:
        logging.debug("can't open file "+dir+obj+".dat")
    else :
        for i in range(len(seqs)) :
            i1=seqs['col4'][i]-3
            i2=seqs['col5'][i]
            alt=[]
            az=[]
            dalt=[]
            daz=[]
            air=[]
            for j in range(i1,i2,3) :
                print(j)
                if j < 1000 : num='{:03d}'.format(j)
                else : num='{:04d}'.format(j)
                logging.debug(dir+root+'.'+num+'.fits')
                try :
                    try:
                        hdulist=fits.open(dir+root+'.'+num+'.fits',ignore_missing_end=True)
                    except:
                        hdulist=fits.open(dir+root+'.'+num+'.fits.bz2',ignore_missing_end=True)
                except :
                    logging.debug("can't open"+dir+root+'.'+num+'.fits')
                else :
                    data=hdulist[0].data-32768
                    alt.append(hdulist[0].header['ALT'])
                    dalt.append(hdulist[0].header['OBS_ALT'] - hdulist[0].header['ALT'])
                    az.append(hdulist[0].header['AZ'])
                    daz.append(hdulist[0].header['OBS_AZ'] - hdulist[0].header['AZ'])
                    air.append(hdulist[0].header['AIRMASS'])
                    fig=plt.figure()
                    ax=fig.add_subplot(111)
                    ax.imshow(data,vmin=min,vmax=max,interpolation='nearest',cmap='Greys_r')
                    if j == i1+3 : 
                        tot=data.astype(float)
                    elif j > i1+3 :
                        tot+=data
                    hdulist.close()
                    del data
                    fig.savefig(num+'.png')
                    plt.close(fig)
                    gifs.append(num+'.png')
                
            gifs.append(obj+'.gif')
            subprocess.call(gifs)
            subprocess.call(['rm']+gifs[3:-1])
        hdu=fits.PrimaryHDU(tot)
        hdu.writeto(dir+obj+'.fits',overwrite=True)
        fig=plt.figure()
        ax=fig.add_subplot(111)
        az=np.array(az)
        alt=np.array(alt)
        daz=np.array(daz)
        dalt=np.array(dalt)
        air=np.array(air)
       
        vmin=tot.min()-3*tot.std() 
        vmax=vmin+20*tot.std()
        ax.imshow(tot,vmin=vmin,vmax=vmax,interpolation='nearest',cmap='Greys_r')
        txt=('{:8.2f}'*4+'{:10.2f}-{:7.2f}').format(az.max()-az.min(),alt.max()-alt.min(),daz.std(),dalt.std(),air.min(),air.max())
        ax.text(0.05,0.9,txt,transform=ax.transAxes,color='w')
        fig.savefig(dir+obj+'.jpg')

if __name__ == '__main__' :
    parser = argparse.ArgumentParser(description="Process 1m images into movies")
    parser.add_argument("root", type=str,help="YYMMDD to process")
    args = parser.parse_args()
    apogee1m(args.root)
 
