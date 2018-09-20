#!/usr/bin/env python

'''apogeeThar.py  check apogee Thar arcs (EM)

Program to check the stability of APOGEE instrument:
 - take ArcLamp THAR quick-red data 
 - takes 1-3 spectra  which is optional
 - do fitting of gaussian function to one spectral line on each chip.
 - compare results with reference, average for mjd=56531, Aug 26,2013, 
             ff=09690003, 09690005, 09690014, 09690016.
  
EM 09/01/2013
Usage:  ./apogeeThar.py <-m1 mjd1> <-m2 mjd2> 

examples:
 ./apogeeThar.py -m1 56531    # reference
 ./apogeeThar.py -m1 56560 -m2  56567
 ./apogeeThar.py -m1 56445  -m2 56537    # 06/01/2013 (56445) -- 09/01/2013 (56537)

History: 
09/11/2013: EM added to apogeeql svn bin repository
09/16/2013: EM reference is the first calibration set after 2013 summer
   mjd=56531, Aug 26,2013, ff=09690003, 09690005, 09690014, 09690016     
09/30/2013: EM:   use quick-red instead of instead of raw data; 
   reorganized output for differences only, 3-lines and 3 rows output, 
   designed for night log. 
09/02/2014  switchehd to sdss-4,  needed more libraries, commented 
    matplotlib for a while until installed.
       
 '''
import argparse
import pyfits, numpy, scipy
#from pylab import *
import scipy.optimize
import sys, os.path, glob
import time

#  Constants
p0a = scipy.c_[53864, 939.646, 1.2745]
p0b = scipy.c_[46184.2, 924.366, 1.071]
p0c = scipy.c_[31715, 1776.62, 0.803]
zone=20
path="/data/apogee/quickred/"

#---------  
def sdth(dth): 
    if dth==12.994:  return "A"
    elif dth==13.499:  return "B"
    else: return "?"    
#-------
def getFileName(path, mjd, chip, fExp="*"):
   return "%s%s/ap1D-%s-%s.fits.fz" % (path, mjd, chip,fExp)       
#------------------------
def curSjd():
# current mjd
  TAI_UTC =34; sjd1=(time.time() + TAI_UTC) / 86400.0 + 40587.3;  
  sjd= int (sjd1)
  return sjd
#----------------------- 

def OneFileFitting(data1, fiber, pRef):
  p0=pRef.copy()
# select zone around line  center
  x1=p0[0][1]-zone;   x2=p0[0][1]+zone
  x=numpy.arange(data1.shape[1])[x1:x2]  #  x-axis array in pix
  spe=data1[fiber,x1:x2]  #read spectrum in the line range 

  ll=numpy.where(spe == max(spe) )
  p0[0][1]= ll[0][0]+x1   

  fitfunc = lambda p0, x: p0[0]*scipy.exp(-(x-p0[1])**2/(2.0*p0[2]**2))
  errfunc = lambda p, x, y: fitfunc(p,x)-y
  p1, success= scipy.optimize.leastsq(errfunc, p0.copy()[0],args=(x,spe))
  
  ref= fitfunc(pRef[0], x)
  fit= fitfunc(p1, x)
  return success, p1, x, spe, ref, fit
  
# -----------------
def checkOneFile(mjd,fibers,fExp,sHeader, outfile):
    for fiber in fibers:      
         ss="%s %3i" %  (sHeader, fiber)  # row
         for chip in ["a","b","c"]: 

            ff= getFileName(path, mjd, chip, fExp)
            hdulist=pyfits.open(ff,'readonly')
            data1=hdulist[1].data
            hdulist.close()

            pp="p0%s" % chip
            pRef=eval(pp)                
            (succ,p1,x,spe,ref,fit)=OneFileFitting(data1, fiber, pRef)                 
            sProf="%5.2f %5.2f %5.2f" % \
                  (p1[0]/pRef[0][0],  # intensity
                   p1[1] - pRef[0][1],  # profile center
                   p1[2] - pRef[0][2] ) # width
            ss="%s  %s" % (ss, sProf)      
                   
   #      print "%s " % ss
         myPrint(ss, outfile)
    return True
  
#-----------------------------------------  
def checkOneMjd(mjd, dither, fibers, outfile):
  mask=getFileName(path, mjd, "a")
  files = sorted(glob.glob(mask))
  if len(files)==0:  
       return False
  for ff in files:  
     hdulist=pyfits.open(ff,'readonly')
     hdr = hdulist[0].header
     hdulist.close()
     dth=float(hdr['DITHPIX'])
#     q0= sdth(dth) in dither
     q1= hdr.get('IMAGETYP')=="ArcLamp"
     q2= hdr.get('LAMPTHAR')==1
#     if q0 and q1 and q2:
     if q1 and q2:
         fExp=ff[35:43] 
         ss="%5s %8s %s" %  (mjd,fExp,sdth(dth))                  
         checkOneFile(mjd,fibers,fExp,ss, outfile)
  return True
     
def myPrint(ss, outfile): 
  print ss
  outfile.write(ss+'\n')
  outfile.flush()

#--------------------------  
def main(argv=None):
  if argv is None: argv = sys.argv[1:]

# current mjd
#  TAI_UTC =34; sjd1=(time.time() + TAI_UTC) / 86400.0 + 40587.3;  sjd= int (sjd1)

  desc = 'apogee arc Thar check for mjd range'
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('-m1', '--mjd1', help='start mjd, default current', default=curSjd()-1,  type=int)
  parser.add_argument('-m2', '--mjd2', default=curSjd(), help='end mjd, default is the same as mjd1', type=int)
  parser.add_argument('-r', '--ref',  action="store_true",  help='references for  lines')
  parser.add_argument('-c', '--cfull',  action="store_true",  help='do three fibers [30,150,270], default one [150]')
  parser.add_argument('-f', '--fiber',  default=150,  help='select one fiber, default [150]', type=int)
  parser.add_argument('-A', '--ditherA', action="store_true",  help='dither A, default both A and B')
  parser.add_argument('-B', '--ditherB', action="store_true", help='dither B, default both A and B')
  parser.add_argument('-o', '--outfile',  help='output file name')
  parser.add_argument('-s', '--sortList',  action="store_true", help='to sort the listing')
  
  args = parser.parse_args()    
  
  mjd1=args.mjd1
  mjd2=args.mjd2
#  if mjd2==None:  mjd2=curSjd()
  mjds=sort(range(mjd1, mjd2+1))

  fibers=[args.fiber]
  if args.cfull: 
      fibers=[30, 150, 270]

  dither=[]
  if args.ditherA: dither.append("A")
  if args.ditherB: dither.append("B")
  if len(dither)==0: 
       dither=["A","B"]

  outfileName=args.outfile
  if outfileName==None:
       outfileName="apThar_%5i.outfile" % (mjd2)
  outfile = open(outfileName,'w')
    
  title="# ./apogeeThar.py -m1 %s  -m2 %s" % (mjd1, mjd2)
  separator="#%s" % ("-"*79)  
#  header1 = "# mjd  nExp  Dth fiber"+" "*8+"a"+" "*18+"b"+" "*18+"c"
#  spHed="I/Io  X-Xo  W-Wo"
#  header2 = "#%s" % (" "*22+spHed+" "*3+spHed+" "*3+spHed)
#  ss="%s\n%s\n%s\n%s" %  (separator, header1, header2, separator)
  def sBl(a):  return "(I%s    X%s    W%s)" % (a,a,a) 
  header="# mjd  nExp  Dth fiber %s   %s   %s" % (sBl("a"),sBl("b"),sBl("c"),) 
  myPrint("%s\n%s\n%s\n%s" %  (title, separator, header, separator), outfile)
  
  if args.sortList:
    for dd in dither:
       for fbr in fibers: 
          ff1=[fbr]
          for mjd in mjds: 
              checkOneMjd(mjd, dd, ff1, outfile)
       myPrint(separator, outfile)        
  else:   
     for mjd in  mjds:
         checkOneMjd(mjd, dither, fibers, outfile)
         myPrint(separator, outfile)
      
  if  args.ref:
     for line in ['a','b','c']:
        pp="p0%s" % line;   p0=eval(pp)     
        myPrint("# ref-%s:  [%6i,  %7.2f,  %4.2f]"% (line, p0[0][0], p0[0][1],p0[0][2]), outfile)

  print "" 


if __name__ == "__main__":
    sys.exit(main())
