#!/usr/bin/env python

'''apogeeTharTest.py  test #2 line fitting, plot optional (EM)

EM,  09/23/2013 -  program for the testing of fitting compare reference. 
The reference  are "/data/apogee/quickred/56531/ap1D-a-09690003.fits.fz"  (A)
                                               "ap1D-a-09690005.fits.fz"  (B)

Usage: 
./apogeeTharProfile.py   -    to run for reference data 

./apogeeTharProfile.py -f 09960005
      print fitting results 

./apogeeTharProfile.py -f 09960005 -p
    print fitting results and  them 


Use average A and B params as reference and initial values for fitting
p0_A - A profile 


 '''
import argparse
import pyfits, numpy, scipy
import scipy.optimize
import sys, os.path, glob
import time

import matplotlib.pyplot as plt
from scipy.interpolate import spline

import apogeeThar as th

#p0a = scipy.c_[53864, 939.646, 1.2745]
#p0b = scipy.c_[46184.2, 924.366, 1.071]
#p0c = scipy.c_[31715, 1776.62, 0.803]


if __name__ == "__main__":
     desc = 'apogeeThar fitting for one image, line #2, plot optional'
     nexpDef='09690003'
     parser = argparse.ArgumentParser(description=desc)
     parser.add_argument('-n', '--nexp', help='exposure number, default=%s, on mjd=56531' % nexpDef,\
             default=nexpDef)
     parser.add_argument('-l', '--line', help='chip and line number (a,b,c), default a',\
             default='a') 
     parser.add_argument('-f', '--fiber', help='fiber, default 150', type=int, default=150) 
     parser.add_argument('-p', '--plot', help='select to plot',\
             default=False, action='store_true')
             
     args = parser.parse_args()    
     pp=args.plot
     if  args.line not in ['a','b','c']:
         sys.exit("not right line")
     fiber=args.fiber 
     if fiber > 300: 
         sys.exit("fiber should be 1:300, you entered %s" % fiber)

     mask="/data/apogee/quickred/*/ap1D-%s-%s.fits.fz" % (args.line, args.nexp)      
     file = glob.glob(mask)
     if len(file)==0: 
          sys.exit("Erros: no file found %s" % file)
     print "file=", file[0]

     # read fits file   
     hdulist=pyfits.open(file[0],'readonly')
     hdr = hdulist[0].header
     data1=hdulist[1].data
     hdulist.close()

     # check is file  SrcLamp and Thar? 
     q1=hdr.get('IMAGETYP')=="ArcLamp"
     q2=hdr.get('LAMPTHAR')==1
     if not(q1 and q2):
          sys.exit("Error: the file is not Thar arc")

     pp="th.p0%s" % args.line
     pRef=eval(pp);    
     (succ,p1,x,spe,ref,fit)=th.OneFileFitting(data1, fiber, pRef)           
     print "success (0-4 ok) =",succ
     print "Relative = [%9.2f,  %7.2f,  %4.2f] " % (p1[0]/pRef[0][0],p1[1] - pRef[0][1], p1[2] - pRef[0][2]   ) 
     print "Fitting  = [%9.2f,  %7.2f,  %4.2f]" % (p1[0], p1[1], p1[2])
     print "Refferenc= [%9.2f,  %7.2f,  %4.2f]" % (pRef[0][0], pRef[0][1],pRef[0][2])
     if  args.plot: print " .. plotting"
     print ""
          
     if not args.plot:
         sys.exit()  # no plotting requested         
         
# plotting
     xnew = numpy.linspace(x[0],x[-1],100)
     spe_smooth = spline(x,spe,xnew)
     fitting_smooth = spline(x,fit,xnew)
     refFit_smooth = spline(x,ref,xnew)

     fig = plt.figure(figsize=(8,5))
     plt.subplot(1, 1, 1)   # horiz, vertical
     
     plt.title('%s' %(file[0])) 
     plt.xlabel('pixels')
     plt.ylabel('data')
     
     rr=70000;  r1=-0.1*rr;  r2=1.1*rr;  
     plt.ylim((r1,r2))   
     plt.xlim((x[0],x[len(x)-1]))   
     
     plt.plot(xnew, refFit_smooth, color='green', )   # reference 
     plt.plot(xnew, spe_smooth, color='black')   # data smooth ilne
     plt.plot(x, spe, 'o', color='black', markersize=3.5)  # data with symbols
     plt.plot(xnew, fitting_smooth, color='red', )  # fitting
     plt.plot([pRef[0][1],pRef[0][1]], [r1,r2], color='black', )  # center for reference

     plt.grid(True, which='both')
     plt.show()
