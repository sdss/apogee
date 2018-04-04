import glob
from astropy.io import ascii
from astropy.io import fits
from holtz.tools import plots
import matplotlib.pyplot as plt
import subprocess
import argparse
import multiprocessing
import logging
import pdb

def apogee1m(root) :

    dir='/home/1m/'+root+'/'
    files=glob.glob(dir+'apogee_*.dat')
    for file in files :
        print file
        fp = open(file,"r")
        for line in fp :
            print line
            obj=line.split()[0]
            #makemovie(root,obj)
            p=multiprocessing.Process(target=makemovie,args=(root,obj,))
            p.start()

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
            for j in range(i1,i2) :
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
                    fig=plt.figure()
                    ax=fig.add_subplot(111)
                    ax.imshow(data,vmin=min,vmax=max,interpolation='nearest',cmap='Greys_r')
                    hdulist.close()
                    del data
                    fig.savefig(num+'.png')
                    plt.close(fig)
                    gifs.append(num+'.png')
                
            gifs.append(obj+'.gif')
            subprocess.call(gifs)
            subprocess.call(['rm']+gifs[3:-1])

if __name__ == '__main__' :
    parser = argparse.ArgumentParser(description="Process 1m images into movies")
    parser.add_argument("root", type=str,help="YYMMDD to process")
    args = parser.parse_args()
    apogee1m(args.root)
 
