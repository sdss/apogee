#!/usr/bin/env python
"""
doc
./noise.py -m 56462
./noise.py -m 56462
./noise.py -m 56756 -m2 56750

"""

import glob
import pyfits, numpy, scipy
import sys, os.path
import argparse
import datetime as dt
import time
import scipy.optimize
import re

import warnings
warnings.filterwarnings('ignore')

# -- 
def curSjd():
# current mjd
  TAI_UTC =34; sjd1=(time.time() + TAI_UTC) / 86400.0 + 40587.3;  
  sjd= int (sjd1)
  return sjd

# -- 
def getMap(files): 
    map=''
    files=sorted(files)
    for i,f in enumerate(files):
        imtype=pyfits.getval(f,'IMAGETYP') 
        if imtype=="Dark": map=map+"d"
        elif imtype=="QuartzFlat": map=map+"q"
        elif imtype=="InternalFlat": map=map+"i"
        elif imtype=="DomeFlat": map=map+"m"
        elif imtype=="ArcLamp":
             map=map+"a"
#            nreads=pyfits.getval(f,'NFRAMES') 
#            if  nreads ==12:   map=map+"t" # Thar
#            elif nreads ==40:  map=map+"u" # Une
#            else: 
#                print f, " ArcLamp with wrong number of reads" 
#                map=map+"_"
        elif imtype=="Object" : map=map+"o" 
        else:
            print f, "wrong image type" 
            map=map+"_" 
    return map
    
# ---
def getMorning(map, files): 
#    regex="dddqqqtutudiiid"
    regex="dddqqqaaaadiiid"
    m = re.search(regex, map)
    if m==None:
     #   print "  no morning sequence found "
        return None
    return files[m.start():m.end()]

# ---
def getStd(f1,f2):
    # check 'NFRAMES', should be 60
    dat1 = pyfits.getdata(f1,0)/pyfits.getval(f1,'NFRAMES')
    dat2 = pyfits.getdata(f2,0)/pyfits.getval(f2,'NFRAMES')
    diff=dat2-dat1
    vlist=list();  y=1024
    for i in range(3):
        for j in range(4):
            x= i*2048+j*512+256
            sub_diff=diff[(y-50):(y+50),(x-50):(x+50)] 
            std=numpy.std(sub_diff)
            vlist.append(std)
    return vlist

def darkOneMjd(mjd): 
    fNames="/data/apogee/utr_cdr/%s/apRaw-%s.fits"%(mjd,"*")
    files = glob.glob(fNames)
    if len(files)==0: 
        print "%s  - no files for this mjd " % (mjd)
      #  sys.exit(" - no files found -- ")
        return 

    files = sorted(files)
    map=getMap(files)
    morning_files=getMorning(map,files)
    if morning_files == None:
         print  "%s - no morning sequence " % (mjd)
         return 
    mm1=morning_files[1]
    mm2=morning_files[2]
    vlist=getStd(mm1,mm2) 

    sa="  "+"".join("%5.2f"%(v) for v in vlist[0:4])   
    sb="  "+"".join("%5.2f"%(v) for v in vlist[4:8])   
    sc="  "+"".join("%5.2f"%(v) for v in vlist[8:13])   
    ss=sa+sb+sc
    print "%s %s  %s%s" % (mjd, mm1[33:41], mm2[33:41], ss)


# ---
def main(argv=None):
    if argv is None: argv = sys.argv[1:]
  
    sjd=curSjd()
    desc = 'APOGEE Dark noises evaluation'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-m', '--mjd', 
           help='mjd, default is current mjd',    
           default=int(sjd), type=int)
    parser.add_argument('-m2', '--mjd2', 
           help='end mjd if check some data range,  default is mjd',    
      #     default=args.mjd, 
           type=int)

    args = parser.parse_args()    
    mjd=args.mjd
    mjd2= args.mjd2;  
    if mjd2==None:  mjd2=mjd 
        
    if mjd > mjd2:
        dd=mjd; mjd=mjd2;  mjd2=dd;
            
    print "APOGEE Dark noises"     
    separator="#%s" % ("-"*79)
    print(separator)
    print "#mjd  Exp1      Exp2  %s%s%s%s%s%s" % (" "*5,"Chip A"," "*16,"Chip B"," "*16,"Chip C")

    for m in  range(mjd, mjd2+1):
        darkOneMjd(m)
    print(separator)
    
        
if __name__ == "__main__":
    sys.exit(main())        
        
        
#y0=1024
# 2048/512
# 0-512


 #   print dat.shape, dat.min(), dat.max() 
#    ss=ss+" "+" ".join("%5.2f"%(v) for v,r in zip(vlist, getDarkRef()))            
#    ss=ss+" "+" ".join("%5.2f"%((v-r)/r*100.) for v,r in zip(vlist, getDarkRef()))            

#APOGEE Dark noises
#mjd  Exp1      Exp2       Chip A                Chip B                Chip C
#56462 09000017  09000018   0.69 0.94 0.69 0.74   0.66 1.64 0.79 0.88   0.75 0.72 1.56 1.05


 
