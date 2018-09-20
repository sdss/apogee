#!/usr/bin/env python

''' EM: program to  check if broken fibers in dome flat apogee data. 

Use VM algorithm to  separate fibers on Dome flat data, calulate total intensity 
on each fibers and compare with master. The selection of broken / faint fibers 
is making based on empirical threshhold.    
The output is copied/pasted to observer's night log after each apogee observations
This program started to be in use since ____ 

This program use apRaw-03720068.fits data as master to  track the 
difference. 

Usage: aptest <mjd> <exposure number>

Example: 

[observer@sos3 ~]$ aptest 56544 09820019
flat= /data/apogee/utr_cdr/56544/apRaw-09820019.fits
master= /home/observer/bin/apRaw-03720068.fits
nfibers= 300
average ratio= 0.801929857011
missing fibers (flat/master < 0.2):
  no missing fibers
faint fibers (0.2 <= flat/master < 0.7) :
  faint: fiber=  15, intFlat= 2872, intMaster= 5262, ratio=0.681
  faint: fiber=  34, intFlat= 2156, intMaster= 5896, ratio=0.456
  faint: fiber= 297, intFlat= 2847, intMaster= 5455, ratio=0.651

It use packages: 
numpy
pyfits
matplotlib
sys

History: 
09/24/2013:  re-structured to have  "__main__" call

'''

from numpy import *
import matplotlib.pyplot as plt
import sys ## ,os, string
import pyfits

if __name__ == "__main__":


    if not sys.argv[1:3]:
          print "usage:  test1.py mjd  exp"
          sys.exit("Error: program did not get mjd and filename, exit")
    #pp="/Users/elenam/idl/fiberplot/data/"
    pp="/data/apogee/utr_cdr/"
    mjd=sys.argv[1]
    exp=sys.argv[2]
    dd="--"
    if  sys.argv[3:4] :
        dd=sys.argv[3]
    #print "dd=",dd   

    flatName="%s%s/apRaw-%s.fits"%(pp,mjd,exp)
    masterName="/home/observer/bin/apRaw-03720068.fits"
    print "flat=", flatName
    print "master=", masterName

    def rdfits(fname):
       hdulist=pyfits.open(fname)
       #hdulist.info()
       dat=hdulist[0].data
       hdulist.close
       return dat
    flatData=rdfits(flatName)
    masterData=rdfits(masterName)

    nd=2048
    col=3000
    dat=zeros((nd), dtype=int64)
    flat=zeros((nd), dtype=int64)
    for i in range(nd):
       dat[i]=masterData[nd-i-1,col]   
       flat[i]=flatData[nd-i-1,col] 

    #print "-"*20
    cutoff=200.;  imax1=500;  imax2=100  # vm cutoff=120
    j=0 # fiber number
    k=0 # pix number in fiber
    ndat=zeros((imax1,imax2), dtype=int64)
    qj=False
    for i in range(nd):
      if(j >= imax1) or (k >= imax2):
             print "break", j,k
             break
      if dat[i] >= cutoff:
        ndat[j,k]=i
    #    print i, dat[i], j, k
        qj=True
        k=k+1
      else:
    #    print i, dat[i], j, "n/a"    
        if qj:
          qj=False
          j=j+1
          k=0
    #      print "   -----"
    nfibers=j
    print "nfibers=",nfibers

    def flux(dat,ndat, nfibers):
      datFlux=zeros((nfibers))
      for j in range(nfibers):
        datFlux[j]=0
        for k in range(10):
          if ndat[j,k] != 0:
            datFlux[j]=datFlux[j]+dat[ndat[j,k]]
      return datFlux

    masFlux=flux(dat,ndat,nfibers)
    #masAver=sum(masFlux)/float(len(masFlux))
    #masFlux=masFlux/masAver
    #masFlux=masFlux/(sum(masFlux)/len(masFlux))

    flatFlux=flux(flat,ndat,nfibers)
    #flatAver=sum(flatFlux)/float(len(flatFlux))
    #flatFlux=flatFlux/flatAver
    #flatFlux=float(flatFlux/(sum(flatFlux)/len(flatFlux)))

    #print "average master / flat intensity=",masAver, flatAver

    fbr=zeros((nfibers))
    for j in range(nfibers):
      fbr[j]=j+1

    ratio=flatFlux/masFlux   #  * (sum(masFlux)/sum(flatFlux) )
    print "average ratio=", sum(ratio)/len(ratio)
    ratio=ratio/(sum(ratio)/len(ratio))

    #for j in range(nfibers):
    #   print j, flatFlux[j], masFlux[j], ratio[j]

    print "missing fibers (flat/master < 0.2):"
    n=0
    for j in range(nfibers):
      if ratio[j] < 0.2: 
          print "   missing: fiber=%4i, intFlat=%5i, intMaster=%5i, ratio=%3.2f" % (fbr[j], flatFlux[j], masFlux[j], ratio[j])
          n=n+1
    if n == 0:
      print "   no missing fibers"
    #print "-"*20

    print "faint fibers (0.2 <= flat/master < 0.7) :"
    n=0
    for j in range(nfibers):
      if (ratio[j] >= 0.2) and (ratio[j] < 0.7):
          print "   faint: fiber=%4i, intFlat=%5i, intMaster=%5i, ratio=%4.3f" % (fbr[j], flatFlux[j], masFlux[j], ratio[j])
          n=n+1
    if n == 0:
      print "   no fainted fibers"
    print "-----------------"

    if dd=="plot":
        xmin=1; xmax=300; ymin=-0.02; ymax=1.35  
        plt.figure(num=None,figsize=(9,4),)
        plt.subplot(111)
        #rr=plt.plot(ratio,"o",color="black")
        rr=plt.plot(fbr,ratio,color="black")
        for j in range(nfibers):
            if ratio[j] < 0.2: 
                rr1=plt.plot(fbr[j],ratio[j],"ro")
        for j in range(nfibers):
            if (ratio[j] >= 0.2) and (ratio[j] < 0.7) : 
                rr1=plt.plot(fbr[j],ratio[j],"bo")
        plt.axis([xmin, xmax, ymin, ymax]); plt.grid(True)
        line1 = plt.plot([xmin,xmax],[0.7,0.7],color='black', linewidth=0.8)
        line2 = plt.plot([xmin,xmax],[1.25,1.25],color='black', linewidth=0.8)
        line3 = plt.plot([xmin,xmax],[0.2,0.2],color='black', linewidth=0.8)
        for j in range(0,nfibers,30):
            line4 = plt.plot([fbr[j],fbr[j]],[ymin,ymax],color='green', linewidth=0.8)
            plt.annotate("%2i"%(fix(j/30)+1), [fbr[j]+10,1.27], color="r", size=13)
        plt.xlabel('fiber number');  plt.ylabel('intensity flat / intensity master');
        plt.title('Apogee fibers relative intensity',size=15)
        plt.show()



