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
import multiprocessing as mp
import logging
import os
import pdb
import numpy as np

def apogee1m(root) :

    dir='/home/1m/'+root+'/'
    files=glob.glob(dir+'apogee_*.dat')
    figs=[]
    ylab=[]
    pars=[]
    # loop over each program
    for file in files :
        print(file)
        fp = open(file,"r")
        # get all of the objects observed for this program
        objs=[]
        for line in fp :
            print(line)
            objs.append(line.split()[0])
        fp.close()

        for obj in set(objs) :
            #makemovie(root,obj)
            #p=multiprocessing.Process(target=makemovie,args=(root,obj,))
            #p.start()
            pars.append((root,obj))
            ylab.append(obj)
            col=[]
            col.append(obj+'.png')
            col.append(obj+'.jpg')
            col.append(obj+'.gif')
            figs.append(col)
    #outputs = []
    #for par in pars :
    #    outputs.append(makemovie(par))

    pool = mp.Pool(32)
    outputs = pool.map_async(makemovie, pars).get()
    pool.close()
    pool.join()
    html.htmltab(figs,file=root+'sum.html',ytitle=ylab)

#def makemovie(root,obj,min=1500,max=2000) :
def makemovie(pars) :
    root=pars[0]
    obj=pars[1]   
    min=1500
    max=2000

    dir='/home/1m/'+root+'/'
    print('dir: ',dir)
    logging.debug('makemovie '+dir+obj)
    try: 
        seqs=ascii.read(dir+obj+'.dat')
    except:
        logging.debug("can't open file "+dir+obj+".dat")
    else :
        for i in range(len(seqs)) :
            i1=seqs['col4'][i]-2
            i2=seqs['col5'][i]
            if i2-i1 < 3 : return
            alt=[]
            az=[]
            dalt=[]
            daz=[]
            air=[]
            gifs=['convert','-loop','0']
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
            print('obj: ',obj)    
            gifs.append(obj+'.gif')
            subprocess.call(gifs)
            os.rename(gifs[3],obj+'.png')
            subprocess.call(['rm']+gifs[4:-1])
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
 
