#!/usr/bin/env python

""" Print list of apogee files in <mjd> directory (or current). 

list_ap1.py  (for current mjd)
list_ap1.py -m  <mjd>  (for other mjd)

Testing:
list_ap1.py -m  56494  # before  summer shakedown
list_ap1.py -m  56531  # reference  
list_ap1.py -m  56553  # after summer shakedown 

History 2013: 
06/22  created by EM
06/25  fixed mjd calculation
06/26  added type of arc lamps 
09/09  added relative position of spectral line #1
        average middle section of raw data to get spe, set  ref position 
09/18  changed spectral relative position to line #2  [941.1758, 941.6995] 
09/24  contacted with JH, he recommended to use quick-reduction data 
      /data/apogee/quickred/56531/ap1D-a-09690003.fits.fz, select one fiber
      (I took 150),  and I use average A and B central position as reference, 939.646. 
      Ref: int=44941.0,  x=939.402,  wg=1.287
      Fit-A: int=44942.1,  x=939.402,  wg=1.287
      Fit-B: int=62787.1,  x=939.890,  wg=1.262
      file-A= /data/apogee/quickred/56531/ap1D-a-09690003.fits.fz 
      file-B= /data/apogee/quickred/56531/ap1D-a-09690005.fits.fz

10/03  added _try_ statement while read fits file  if it is not available
10/04  I stopped errro messages, but I still got pyfits warnings, 
        I stopped them by  setting warnings.filterwarnings('ignore')
11/12  fixed bug with format message if function cannot read file   fits file     
12/10 fixed bug  fitting gaussian function if offset large  (copied function from apogeeThar);
    reformat - repalce ArcLams t Arc,  QuaFlat, and IntFlat; added dither set, format offset 
    to 4.1f (it was  5.2f). 

01/09/2014 Added column with  normed flux for all three flats.  
       test:  ./list_ap1.py -m 56721 -morning
June/19  added special list for dither analysis  
   possible tests: 
   list_ap1.py -m 56803 -dither
   list_ap1.py -m 56803
   list_ap1.py -m 56803 -evening
   list_ap1.py -m 56803 -morning

Sept, 08, 2014  changes A.B dithers settings: 
   A 13.0 --> 10.0
   B 13.5 --> 10.5
Oct,15, 2014  removed changes on Sept, 08, 2014

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

import apogeeThar

#p0 = scipy.c_[44941, 939.646, 1.287] #   A and B average for Line 2

zone=20
p0=apogeeThar.p0a.copy()

morning="dddqqqtutudiiid"
evening="dqtutud"
dither='t'*17

ditherA=0
regex=""


#------------------------
def curSjd():
# current mjd
  TAI_UTC =34; sjd1=(time.time() + TAI_UTC) / 86400.0 + 40587.3;  
  sjd= int (sjd1)
  return sjd

#-----
def getFiles(mask):
    files = glob.glob(mask)
    return sorted(files)

#----------------------- 
def ifLamp(f, lamp ):   #  'LAMPUNE'; 'LAMPTHAR'
    fexp=f[33:41]
    fReduct="/data/apogee/quickred/?????/ap2D-a-%s.fits.fz" % (fexp)  
    fileNames = glob.glob(fReduct)
    if len(fileNames) > 1:  sys.exit(" -- error: several reduced arc files found -- ")
    if len(fileNames) > 1:  sys.exit(" -- error: no reduced arc file found -- ")
    lamp=pyfits.getval(fileNames[0],lamp,0)
    if lamp == 1:  return True
    else: return False

#-------------  
def getMap(files):
    if len(files)==0:  return None
    files=sorted(files)
    map=''    
    for i,f in enumerate(files):
        if not os.path.isfile(f):
#            print "%3i    no file: %s " % (i,f)
            map=map+"-"
            continue
        hdr = pyfits.getheader(f)
        imtype=hdr.get('IMAGETYP')
        nreads=hdr.get('NFRAMES')
        if imtype=="Dark": map=map+"d"
        elif imtype=="QuartzFlat": map=map+"q"
        elif imtype=="InternalFlat": map=map+"i"
        elif imtype=="DomeFlat": map=map+"m"
        elif imtype=="ArcLamp" and ifLamp(f, 'LAMPTHAR'):  
               map=map+"t" # Thar
        elif imtype=="ArcLamp" and ifLamp(f, 'LAMPUNE'):  
               map=map+"u" # Une
        elif imtype=="Object" : 
               map=map+"o" 
        else:
            print "wrong image type, f=%s, imtype = %s" % (f,imtype) 
            map=map+"?"
      #  print i, f, map
    return map             

# ----
def getSequence(map, regex): 
    m = re.search(regex, map)
    if m==None:
        return None
    else:
        return  [m.start(), m.end()]

# ------
def getOffset(qrfile1):
    try :
        data1 = pyfits.getdata(qrfile1,0)
    except  :
        return None 
    success, p1, x, spe, ref, fit =apogeeThar.OneFileFitting(data1, 150, p0)
    if success==5:
        return " err"
    else:
        return "%5.2f" % (p1[1] - p0[0][1])

#----
def getStd(f): 
    hdr = pyfits.getheader(f)
    nreads=hdr.get('NFRAMES')
    dat = pyfits.getdata(f,0)/nreads
    dat=numpy.array(dat)    
    def ff(i, dat): 
        y=1024;   x=i*2048+1024 
        dat=dat[(y-200):(y+200),(x-100):(x+100)] 
        return "%4.2f" % numpy.std(dat)
#    return "(%s, %s, %s)" % (ff(0, dat), ff(1, dat), ff(2, dat))
    return "%s" % (ff(2, dat))
    
    	
# ......
def getFlux(f):
    hdr = pyfits.getheader(f)
    nreads=hdr.get('NFRAMES')
    imtype=hdr.get('IMAGETYP')
    if imtype not in ["QuartzFlat","InternalFlat","DomeFlat"]:
        return " - "
    dat = pyfits.getdata(f,0)/nreads
    dat=numpy.array(dat)
    y=1024;  x= 1024
    dat=dat[(y-200):(y+200),(x-100):(x+100)] 
    med=numpy.mean(dat)
    if imtype=="QuartzFlat":  nrm=167.05
    elif imtype=="InternalFlat":  nrm=94.05
    elif imtype=="DomeFlat":  nrm=81.525
    else: nrm=None
    if nrm != None:  
        return "%3i" % (med/nrm*100)
    else:  return " ? "

                    
def  list_one_file(i,f,mjd, dither1=False):
    path="/data/apogee/quickred/%s" % mjd
    fexp=f[33:41]
    qrfile1="%s/ap1D-a-%s.fits.fz" % (path,fexp)
    qrfile2="%s/ap2D-a-%s.fits.fz" % (path,fexp)
    
    ff=f
    if os.path.exists(qrfile2):  
        ff=qrfile2

    q=False
    for j in range(3):
      try :
        hdr = pyfits.getheader(ff)
        q=True
        break
      except :    # it was except IOError:
        continue 
    if not q:
       print "    cannot read file : %s" % ff 
       return
    
    ct=hdr.get('CARTID'); plate=hdr.get('PLATEID'); 
    if ct == None: ct="--"
    if plate==None: plate="----"

    dth= float(hdr['DITHPIX'])
    
    if dth==12.994: sdth="A  %4.1f" % dth
    elif dth==13.499: sdth="B  %4.1f" % dth
    else: sdth="?  %4.1f"% dth 

    #dth=round(dth, 1)
  #  if dth==10.0: sdth="A  %4.1f" % dth
  #  elif dth==10.5: sdth="B  %4.1f" % dth
  #  else: sdth="?  %4.1f"% dth 
        
    imtype= hdr.get('IMAGETYP')
    offset=" - "
    if imtype=="ArcLamp":
      imtype="Arc"
      if hdr.get('LAMPUNE')==1:  imtype=imtype+"-Une"
      elif hdr.get('LAMPTHAR')==1:
          imtype=imtype+"-Thar"
          offs=getOffset(qrfile1)
          if offs != None: 
              offset=offs              
      else: imtype=imtype+"----"
    if imtype=="QuartzFlat":  imtype="QuarFlat"
    if imtype=="InternalFlat": imtype="InteFlat"    
    imtype=imtype.center(10)

    arc=list("x-x-x")    
    for k,l in enumerate(["a","b","c"]):
        pp="/data/apogee/archive/%s/apR-%s-%s.apz"%(mjd,l,f[33:41])
        if os.path.exists(pp): arc[2*k]=l 
    
    flux=" - "    
    if hdr.get('IMAGETYP') in ["QuartzFlat","InternalFlat","DomeFlat"]:
            flux=getFlux(f)

#    std="  - "
##    print imtype, hdr.get('NFRAMES'), len(imtype)
#    if (imtype.strip(' \t\n\r') == "Dark") and (hdr.get('NFRAMES') == 60):
#            std=getStd(f)
    
# print information
    ss1="%s %3i "% (mjd, i+1)  #i
    ss1=ss1+"%s  " % (hdr['DATE-OBS'][11:16]) # UT time
    ss1=ss1+"%s " % (f[33:41])  # exp number
    
    if dither1:
        ss1=ss1+" %s " % (sdth)  # dither            
        ss1=ss1+" %5s " % (offset)  # offset
        dDth=[0,-2,0,+2,0,-4, 0,+4,0,-3,0,2,0,-1,0,0.5,0]
        ss1=ss1+"%5s " % (dDth[i])
        global ditherA
        if i==0:  
            ditherA=dth
        dd=float(dth - ditherA)
      #  print dd
        ss1=ss1+"    %5.2f" % (dd)
        ss1=ss1+"    %5.2f" % (dDth[i]-dd)        
    else: 
        ss1=ss1+"%s" % (imtype)  # image type
        ss1=ss1+"%2i " %  hdr.get('NFRAMES')  # nframes
        ss1=ss1+" %s " % (sdth)  # dither            
        ss1=ss1+" %2s-%4s  " % (ct, plate)
        ss1=ss1+"%s "%"".join(arc)    # archive file existence
        ss1=ss1+"%5s " % (offset)  # offset
        ss1=ss1+"%3s " % (flux)  # flux
        #   ss1=ss1+" %s " % (std)  # std        
        comm=hdr["OBSCMNT"]
        if comm  !="None": 
            ss1=ss1+ " %s" % comm[0:8] # comment
#    if (regex==dither):
#        ss1=ss1+" %5.2f" % diffDither
    
    print ss1 #, len(ss1)


def  list_one_mjd(mjd, args):
#    bs=os.path.basename(sys.argv[0])
#    bs=os.path.abspath(sys.argv[0])
#    print  "# %s -m1 %s" % (bs, mjd)
        
#    print "APOGEE data list,   mjd=%s" % mjd    
    pp="/data/apogee/utr_cdr/"
    fNames="%s%s/apRaw-%s.fits"%(pp,mjd,"*")   # raw data
#    print "   raw_data: ", fNames
    fReduct="/data/apogee/quickred/%s/ap2D-a-*.fits.fz" % (mjd)  # quick reduction
#    print "   quick_red:", fReduct
    fArch="/data/apogee/archive/%s/apR-[a,b,c]-*.apz" % (mjd)  # archive
#    print "   archive:  ", fArch
       
    files=getFiles(fNames)
    if len(files)==0: 
        print " - mjd=%s, no files found  -- " % (mjd)        
        #sys.exit(" - mjd=%s, no files found  -- " % (mjd))
        return 

    files=sorted(files)
    lastNumber=files[-1][33:41]
    dd=int(round(int(lastNumber)/10000)*10000)
    s2=int(lastNumber)-dd

    files1=[]
    for i in range(1,s2+1):
        fItem="%s%s/apRaw-%8i.fits"%(pp,mjd,dd+i)   # raw data 
        files1.append(fItem)

    files=files1
    dither1=False

    if args.morning or args.evening or args.dither:
        map=getMap(files)
        if args.morning: regex=morning
        elif args.evening:  regex=evening
        elif args.dither:  
            regex=dither
            dither1=True 
        else:   sys.exit("  error: 1")   
        mapRange=getSequence(map, regex)
        if mapRange==None:
            print  "  no sequence found"
            return              
    #    print mapRange[0], mapRange[1]
        files=files[mapRange[0]:mapRange[1]]
    #    print len(files)    
                            
    for i,f in enumerate(sorted(files)):
        if not os.path.isfile(f):
            ss1="%s %3i"% (mjd, i+1)  #i
            print "%s      no file found:  %s " % (ss1, f)
            continue
        list_one_file(i,f, mjd, dither1)   

def main(): 
 #   print ("my_test =====")
    sjd=curSjd()
    desc = 'list of files for one night of apogee observations'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-m', '--mjd',  help='enter mjd, default is current mjd',    
           default=int(sjd), type=int)
    parser.add_argument('-m2', '--mjd2',  help='and mjd if use range of dates, default is mjd', type=int)           
    parser.add_argument('-morning', '--morning', action="store_true") 
    parser.add_argument('-evening', '--evening', action="store_true") 
    parser.add_argument('-dither', '--dither', action="store_true")            
    args = parser.parse_args()    
    mjd=args.mjd
    if args.mjd2==None: mjd2=mjd
    else: mjd2=args.mjd2
    if mjd2 < mjd:
        dd=mjd; mjd=mjd2; mjd2=dd
    
    if args.morning: request="morning"
    elif args.evening: request="evening"
    elif args.dither:  request="dither"
    else: request="all"
    
    prc="%"
    if request=="dither":
        line="-"*67
        header="mjd     i   UT   File/Exp  A/B Dth  Offset Command   Moved   Error" 
    else:
        line="-"*80 
        # header=" i   UT   File/Exp   Imtype  Nread  Dth  Ct-Plate Arch Offset %sFlux Std Comm" % prc
        header="mjd     i   UT   File/Exp  Imtype Nread A/B Dth  Ct-Plate Arch Offset %sFlux Comm" % prc

    for m in range(mjd, mjd2+1):
        pp="/data/apogee/utr_cdr/"
        fNames="%s%s/apRaw-%s.fits"%(pp,m,"*")   # raw data
        fArch="/data/apogee/archive/%s/apR-[a,b,c]-*.apz" % (mjd)  # archive
        fReduct="/data/apogee/quickred/%s/ap2D-a-*.fits.fz" % (mjd)  # quickred
        print "APOGEE data list,   mjd=%s,  %s" % (m, request)  
        print "raw_data: ", fNames  
        print "archive: ", fArch  
        print "quickred: ", fArch  
        
        
        print line, "\n",  header, "\n", line        
        files=getFiles(fNames)
        if len(files)==0: 
             print " - no files found  for this mjd  %s-- " % m 
             print line, "\n" 
             continue
        list_one_mjd(m, args)
        print line, "\n" 
        
if __name__ == "__main__":
    sys.exit(main())
        