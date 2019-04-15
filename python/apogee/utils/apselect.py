# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: March 2018
# @Filename: synth.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

from astropy.io import fits
from astropy.io import ascii
import os
from sdss_access.path import path
from sdss_access.sync.http import HttpAccess
import pdb
import sys
from astroquery.gaia import Gaia

import numpy as np
import copy
from tools import html
from tools import plots
from tools import match
from apogee.utils import bitmask
import matplotlib.pyplot as plt

def select(data,badval=None,badstar=None,logg=[-1,10],teff=[0,10000],mh=[-100.,100.],alpha=[-100.,100.],sn=[0,1000], raw=False, 
           glon=[-1,360],glat=[-90,90],vscatter=[-1,100],grid=None,field=None,giants=None, dwarfs=None,rgb=None, rc=None,inter=None, 
           id=None, redid=None, badtarg=None, gdtarg=None) :
    '''  
    Return indices of requested subsamples from input allStar structure 

    Args:
        data : allStar structure

    Keyword args:
        badval (char) : ASPCAPFLAG value or list of values to reject
        logg  (float[2]) : min, max log g to accept (default [-1,10])
        teff  (float[2]) : min, max Teff to accept (default [0,10000])
        mh    (float[2]) : min, max [M/H] to accept (default [-100,100])
        alpha (float[2]) : min, max [alpha/H] to accept (default [-100,100])
        sn    (float[2]) : min, max S/N to accept (default [0,1000])
        glon  (float[2]) : min, max GLON to accept (default [0,360])
        glat  (float[2]) : min, max GLAT to accept (default [-90,90])
        grid  (char)     : ASPCAP_CLASS to accept, if specified (default None)
        field (char)     : FIELD to accept, if specified (default None)
        locid (char)     : FIELD to accept, if specified (default None)
        giants (bool)    : ASPCAP giants if true, per line in HR diagram, else all (default = None)
        dwarfs (bool)    : ASPCAP dwarfs if true, per line in HR diagram, else all (default = None)
        rgb    (bool)    : ASPCAP RGB per Teff, [M/H], [C/N] criterion (default = None)
        rc     (bool)    : ASPCAP RC per Teff, [M/H], [C/N] criterion (default = None)
        inter  (bool)    : ASPCAP "intermediate" RGB/RC per Teff, [M/H], [C/N] criterion (default = None)
        raw   (bool)     : specifies raw (FPARAM) values rather than calibrated (PARAM) (default=False)
        id    (char)     : APOGEE_ID to return (ignores other constraints)
        redid (char)     : REDUCTION_ID to return (ignores other constraints)
 
    '''

    if id is not None :
        gd = np.where(np.core.defchararray.find(data.APOGEE_ID,id) >=0 )
        return gd[0]
 
    if redid is not None :
        gd = np.where(np.core.defchararray.find(data.REDUCTION_ID,redid) >=0 )
        return gd[0]
    
    if alpha is None :
        alpha=[-100,100]    
    if mh is None :
        mh=[-100,100]    
    if glon is None :
        glon=[0,360]
    if glat is None :
        glat=[-90,90]

    # filter for bad ASPCAPFLAG 
    aspcapflag=bitmask.AspcapBitMask()
    if type(badval) is str :
        badval = [badval]
    badbits = 0
    if badval is not None : 
        for name in badval :
            badbits = badbits | aspcapflag.getval(name)
    try :
       bad = data['ASPCAPFLAG'] & badbits
    except :
       bad = np.zeros(len(data),dtype=np.int8)

    # filter for bad STARFLAG
    starflag=bitmask.StarBitMask()
    if type(badstar) is str :
        badstar = [badstar]
    badbits = 0
    if badstar is not None : 
        for name in badstar :
            badbits = badbits | starflag.getval(name)
    try :
       bad = bad | (data['STARFLAG'] & badbits)
    except :
       pass

    # filter for bad TARGFLAG
    targ = np.zeros(len(data),dtype=np.int8)
    if badtarg is not None :
        if type(badtarg) is str :
            badtarg=[badtarg]
        for val in badtarg :
            j=np.where(np.core.defchararray.find(data['TARGFLAGS'],val) >= 0)[0]
            targ[j]=1
            print(val, len(j))

    # filter for good TARGFLAG if not None, supercedes (default is all bad unless specified)
    if gdtarg is not None :
        if type(gdtarg) is str :
            gdtarg=[gdtarg]
        targ = np.ones(len(data),dtype=np.int8)
        for val in gdtarg :
            j=np.where(np.core.defchararray.find(data['TARGFLAGS'],val) >= 0)[0]
            targ[j]=0
   
    if raw :
        param='FPARAM'
    else :
        param='PARAM'
    t = data[param][:,0] 
    g = data[param][:,1] 
    m = data[param][:,3] 
    a = data[param][:,6] 
    vscat = data['VSCATTER']

    if giants is not None or dwarfs is not None :
        startype = (g < 2./1300.*(t-3500.)+2) & (g < 4) & (t < 7000)
        if dwarfs is not None :
           startype = np.logical_not(startype)
    else  :
        # take all stars
        startype = g < 100

    # define dt from Bovy ridgeline as f([M/H])
    dt=data['FPARAM'][:,0]-(4468+(data['FPARAM'][:,1]-2.5)/0.0018 - 382.5*data['FPARAM'][:,3])
    cn=data['FPARAM'][:,4]-data['FPARAM'][:,5]
    if rgb :
       startype = np.logical_and(startype,
                  (dt<0) | (data['FPARAM'][:,1] < 2.3) |
                  ((dt>0) & (dt<100) & (cn<-0.1-0.005*dt)))
    if rc :
       startype = np.logical_and(startype,
                  (data['FPARAM'][:,1] > 2.3) &
                  ((dt>100) | ((dt>0) & (cn>-0.075-0.0025*dt)) ))
    if inter :
       startype = np.logical_and(startype,
                  (data['FPARAM'][:,1] > 2.3) & (dt>0) & (dt<100) &
                  (cn>-0.1-0.005*dt) & (cn<-0.075-0.0025*dt))

    try :
       snr = data['SNR']
    except :
       snr = data['SNR_2']

    gd = np.where((bad == 0)  & (targ == 0) &
         (startype) &
         (t >= teff[0]) & (t <= teff[1])  &
         (g >= logg[0]) & (g <= logg[1])  &
         (m >= mh[0]) & (m <= mh[1])  &
         (a >= alpha[0]) & (a <= alpha[1])  &
         (vscat >= vscatter[0]) & (vscat <= vscatter[1])  &
         (data['GLON'] >= glon[0]) & (data['GLON'] <= glon[1])  &
         (data['GLAT'] >= glat[0]) & (data['GLAT'] <= glat[1])  &
         (snr >= sn[0]) & (snr <= sn[1])  
         )[0]

    if grid is not None :
        try: gdclass = np.where(np.core.defchararray.find(data['ASPCAP_CLASS'][gd],grid) >=0 )[0]
        except: gdclass = np.where(np.core.defchararray.find(data['CLASS'][gd],grid) >=0 )[0]
        gd=gd[gdclass]

    if field is not None :
        gdfield = np.where(np.core.defchararray.strip(data['FIELD'][gd]) ==field)[0]
        gd=gd[gdfield]

    return gd

def clustdata() :
    """
    Returns structure containing cluster data
    """

    clust=['M92','M15','M53','N5466','N4147',
        'M2','M13','M3','M5','M12','M107',
        'M71','N2243','Be29', 'N2158','M35','N2420',
        'N188','M67','N7789','Pleiades','N6819',
        'N6791',
        'N5053','M68','N6397','M55','N5634','M22','M79','N3201','M10',
        'N6752','Omegacen','M54','N6229','Pal5','N6544','N6522','N288','N362','N1851',
        'M4','N2808','Pal6','47TUC','Pal1','N6539','N6388','Terzan12','N6441','N6316',
        'N6760','Terzan5','N6553','N6528']

    out = np.recarray(len(clust),dtype=[
                       ('name','U24'),
                       ('field','U24'),
                       ('rv','f4'),
                       ('drv','f4'),
                       ('mh','f4'),
                       ('dist','f4'),
                       ('age','f4'),
                       ('giant_mass','f4'),
                       ('ra','f4'),
                       ('dec','f4'),
                       ('rad','f4'),
                        ('ebv','f4')
                        ])
    out['ebv']=[0.02,0.1,0.02,0.,0.02,
                0.06,0.02,0.01,0.03,0.19,0.33,
                0.25,0.051,0.157,0.36,0.262,0.05,
                0.08,0.04,0.217,0.03,0.16,
                0.122,
                0.01,0.05,0.18,0.08,0.05,0.34,0.01,0.24,0.28,
                0.04,0.12,0.15,0.01,0.03,0.76,0.48,0.03,0.05,0.02,
                0.35,0.22,1.46,0.04,0.15,1.02,0.37,2.06,0.47,0.54,
                0.77,2.28,0.63,0.54]

    out['name']=clust
    out['field']=['M92','M15','M53','N5466','N4147',
                'M2','M13','M3','M5PAL5','M12','M107',
                'M71','N2243','198+08', 'M35N2158','M35N2158','N2420',
                'N188','M67','N7789','Pleiades','N6819',
                'N6791',
                'N5053','M68','N6397','M55','N5634','M22','M79','N3201','M10',
                'N6752','Omegacen','M54','N6229','Pal5','N6544','N6522','N288','N362','N1851',
                'M4','N2808','Pal6','47TUC','Pal1','N6539','N6388','Terzan12','N6441','N6316',
                'N6760','Terzan5','N6553','N6528']

    out['rv']=[ -118.517, -107.508, -61.5988, 106.883, 183., 
              -3.74874, -246.589, -145.525, 54.8727, -41.4000, -35.2638,
              -23.5504, 60., 25., 26.6529, -7.02990, 74.3025, 
              -42.1051, 34.0525, -55.0546, 5.59228, 2.55553, 
              -47.4558,
               44.00,-94.70,18.80,174.70,-45.10,-146.30,205.80,494.00,75.20,
              -26.70,232.10,141.30,-154.20,-58.70,-27.30,-21.10,-45.40,223.50,320.50,
               70.70,101.60,181.00,-18.00,-82.80,31.00,80.10,94.10,16.50,71.40,
               -27.50,-93.00,-3.20,206.60]
    out['drv']=[ 12., 12., 10., 10., 10., 
               10., 12., 12., 10., 10., 8., 
               8., 10., 10., 10., 10., 10., 
               6., 6., 6., 10., 6., 
               6.,
               10.00,10.00,10.0,8.00,10.00,12.00,10.00,10.00,12.00,
               12.00,25.00,10.50,12.00,1.10,10.00,10.00,10.00,10.00,10.40,
               15.00,15.00,10.00,20.00,10.00,10.00,18.90,10.00,10.00,10.00,
               10.00,10.00, 10.00,10.00]
    out['mh']=[-2.35,-2.33,-2.06,-2.01,-1.78,
             -1.66,-1.58,-1.50,-1.33,-1.37,-1.03,
             -0.82,-0.35,-0.44,-0.21,-0.14,-0.13,
              0.04, 0.06,0.09,0.03,0.16,
              0.37,
             -2.27,-2.23,-2.02,-1.94,-1.88,-1.70,-1.60,-1.59,-1.56,
             -1.54,-1.53,-1.49,-1.47,-1.41,-1.40,-1.34,-1.32,-1.26,-1.18,
             -1.16,-1.14,-0.91,-0.72,-0.65,-0.63,-0.55,-0.50,-0.46,-0.45,
             -0.40,-0.23,-0.18,-0.11]
    out['dist']=[8.3, 10.4, 17.9, 16.0, 19.3, 
               12.59, 7.1, 10.2, 7.5, 6.33, 6.4, 
               4.0, 4.45, 14.87, 5.06, 0.816, 2.44, 
               2.04, 0.907, 2.33, 0.444, 2.36, 
               4.09,
               17.40,10.30,2.30,5.40,25.20,3.20,12.90,4.90,4.40,
               4.00,5.20,26.50,30.50,23.20,3.00,7.70,12.00,9.40,12.10,
               2.20,9.60,5.80,4.50,11.10, 7.80,9.90,4.80,11.60,10.40,
               7.40,6.90,6.00,7.90]
    out['age']=[12., 12., 12., 12., 12., 
               12., 12., 12., 12., 12., 12., 
               10.0, 4.5, 1.1, 1.0, 0.5, 1.1, 
               4.26, 6., 1.7, 0.5, 1.5, 
               4.4,
               12., 12., 12., 12., 12., 12., 12., 12., 12., 
               12., 12., 12., 12., 12., 12., 12., 12., 12., 12., 
               12., 12., 12., 12., 12., 12., 12., 12., 12., 12., 
               12., 12., 12., 12.]
    out['giant_mass'] = [ 0.85, 0.85, 0.85, 0.85, 0.85,
                          0.85, 0.85, 0.85, 0.85, 0.85, 0.85,
                          0.85, -9999.,  -9999.,  -9999., -9999., 1.6,
                          -9999., 1.36, -9999., -9999., 1.63,
                          1.15,
                          0.85, 0.85,0.85,0.85, 0.85,0.85, 0.85,0.85, 0.85,
                          0.85, 0.85,0.85, 0.85,0.85, 0.85,0.85, 0.85,0.85, 0.85,
                          0.85, 0.85,0.85, 0.85,0.85, 0.85,0.85, 0.85,0.85, 0.85,
                          0.85, 0.85,0.85, 0.85]
    out['ra']=[259.27917,322.4929,198.2292,211.3625,182.5262,
             323.3625,250.42083,205.55,229.6375,251.80908,248.1333,
             298.4438,97.3917,103.325,91.8542,92.25,114.5958,
             12.1083,132.825,359.35,56.75,295.325,
             290.22083,
             199.11288 , 189.86658 , 265.17537 , 294.99879 , 217.40541 ,
             279.09975 ,  81.046215, 154.403415, 254.28771 , 287.71713 ,
             201.69684 , 283.76388 , 251.744955, 229.021875, 271.835745,
             270.89166 ,  13.188495,  15.809415,  78.528165, 245.896755,
             138.012915, 265.925835,   6.023625,  53.333505, 271.207005,
             264.07179 , 273.065835, 267.554415, 259.15542 , 287.800035,
             267.019995, 272.323335, 271.20684 ]

    out['dec']=[43.1358,12.1669,18.1806,28.5344,18.5425,
              -0.82325,36.4611,28.3772,2.0811,-1.94853,-13.05361,
              18.7792,-31.2833,16.9167,24.0967,24.35,21.5733,
              85.255,11.8,56.7083,24.11667,40.1867,
              37.77167,
              17.700250,-26.744056, -53.674333,-30.964750, -5.976389,-23.904750,-24.524722,-46.412472, -4.100306,
             -59.984556,-47.479583,-30.479861, 47.527750, -0.111611,-24.997333,-30.033972,-26.582611,-70.848778,-40.046556,
             -26.525750,-64.863500,-26.222500,-72.081278, 79.581056, -7.585861,-44.735500,-22.741944,-37.051444,-28.140111,
              1.030472,-24.779167,-25.908694,-30.056278]
    out['rad']=[20.,18.,13.,11.,4.,
             16.,20.,24.,23.,16.,13.,
             7.,5.,5.,5.,25.,6.,
             16.,45.,16.,110.,12.,
             16.,
             11.80,13.70,15.80,16.30,8.40,29.00,8.30,28.50,21.50,
             55.30,57.00,7.50,5.40,16.30,2.05,16.40,12.90,16.10,11.70,
             32.50,15.60,8.40,42.90,9.00,21.50,6.20,5.00,12.00,5.90,
             13.00,13.30,8.20,16.60]

    return out.view(np.recarray)


def clustmember(data,cluster,logg=[-1,3.8],te=[3800,5500],rv=True,pm=True,dist=True,raw=False,firstgen=False,firstpos=True,plot=False,hard=None) :

    clust=clustdata()
    ic = np.where( np.core.defchararray.strip(clust.name) == cluster)[0]
    if len(ic) == 0 :
        print('no cluster found: ',cluster)
        return []
    print('cluster: ', cluster)
    ic=ic[0]

    # adjust ra for wraparound if needed
    ra=copy.copy(data['RA'])
    if clust[ic].ra > 300 :
        j=np.where(data['RA'] < 180)[0]
        ra[j]+=360
    if clust[ic].ra < 60 :
        j=np.where(data['RA'] > 180)[0]
        ra[j]-=360

    # select by location relative to cluster
    jc=np.where((np.abs(ra-clust[ic].ra)*np.cos(clust[ic].dec*np.pi/180.) < clust[ic].rad/60.) & 
                (np.abs(data['DEC']-clust[ic].dec) < clust[ic].rad/60.))[0]
    if len(jc) > 0 :
        j=np.where( ((ra[jc]-clust[ic].ra)*np.cos(clust[ic].dec*np.pi/180.))**2+ 
                     (data[jc]['DEC']-clust[ic].dec)**2 < (clust[ic].rad/60.)**2)[0]
        jc=jc[j]
    else :
        jc=[]
    print('{:d} stars after location criterion'.format(len(jc)))
    if len(jc) == 0 : return jc
    if plot :
        jf=np.where((np.abs(ra-clust[ic].ra)*np.cos(clust[ic].dec*np.pi/180.) < 1.5) & 
                (np.abs(data['DEC']-clust[ic].dec) < 1.5))[0]
        fig,ax=plots.multi(1,1)
        fig.suptitle('{:s} Radius: {:4.2f} arcmin'.format(cluster,clust[ic].rad))
        plots.plotp(ax,ra[jf],data['DEC'][jf],color='k',size=20,draw=False,xt='RA',yt='DEC')
        plots.plotp(ax,ra[jc],data['DEC'][jc],color='g',size=20,draw=False)
        circle = plt.Circle((clust[ic].ra,clust[ic].dec), clust[ic].rad/60., color='g', fill=False)
        ax.add_artist(circle)

        if hard is not None :
            print(hard+'/'+clust[ic].name+'_pos.jpg')
            fig.savefig(hard+'/'+clust[ic].name+'_pos.jpg')
            plt.close()
        else :
            pdb.set_trace()

    # RV criterion
    try :
        vhelio = data['VHELIO']
    except :
        vhelio = data['VHELIO_AVG']
    j=np.where(np.abs(vhelio[jc]-clust[ic].rv) < clust[ic].drv)[0]
    if plot :
        ax.cla() 
        try :
            if len(jf) > 0 : ax.hist(vhelio[jf],color='k',bins=np.arange(clust[ic].rv-100,clust[ic].rv+100,1.),histtype='step')
            if len(jc) > 0 : ax.hist(vhelio[jc],color='r',bins=np.arange(clust[ic].rv-100,clust[ic].rv+100,1.),histtype='step')
            if len(j) > 0 : ax.hist(vhelio[jc[j]],color='g',bins=np.arange(clust[ic].rv-100,clust[ic].rv+100,1.),histtype='step',linewidth=3)
        except : pass
        ymax=ax.get_ylim()[1]
        ax.plot([clust[ic].rv,clust[ic].rv],[0,ymax],color='g')
        ax.plot([clust[ic].rv+clust[ic].drv,clust[ic].rv+clust[ic].drv],[0,ymax],color='r',ls=':')
        ax.plot([clust[ic].rv-clust[ic].drv,clust[ic].rv-clust[ic].drv],[0,ymax],color='r',ls=':')
        fig.suptitle('{:s} RV: {:4.2f} +/- {:4.2f}'.format(cluster,clust[ic].rv,clust[ic].drv))
        ax.set_xlabel('RV')
        if hard is not None :
            fig.savefig(hard+'/'+clust[ic].name+'_rv.jpg')
            plt.close()
        else :
            plt.draw()
            pdb.set_trace()
    if len(j) > 0 :
        if rv: jc=jc[j]
    else :
        jc=[]
    print('{:d} stars after RV criterion'.format(len(jc)))
    if len(jc) <= 1 : return jc

    # proper motion criterion
    if dist or pm : gaia = True
    else : gaia = False
    if gaia :
      job=Gaia.launch_job_async("SELECT xm.original_ext_source_id, gaia.pmra, gaia.pmra_error, gaia.pmdec, gaia.pmdec_error, gaia.parallax, gaia.parallax_error "+
                               "FROM gaiadr2.gaia_source AS gaia, gaiadr2.tmass_best_neighbour AS xm "+
                               "WHERE gaia.source_id = xm.source_id AND "+
                                  "CONTAINS(POINT('ICRS',gaia.ra,gaia.dec),CIRCLE('ICRS',{:12.6f},{:12.6f},{:12.6f}))=1;".format(
                                                                                  clust[ic].ra,clust[ic].dec,clust[ic].rad/60.))
      # convert to velocities (note mas and kpc cancel out factors of 1000) and get median
      gaia=job.get_results()
      i1, i2 = match.match(np.core.defchararray.replace(data['APOGEE_ID'][jc],'2M',''),gaia['original_ext_source_id'])
      vra=4.74*gaia['pmra']*clust[ic].dist
      vdec=4.74*gaia['pmdec']*clust[ic].dist
      med_vra=np.median(vra[i2])
      med_vdec=np.median(vdec[i2])
      j=np.where((vra[i2]-med_vra)**2+(vdec[i2]-med_vdec)**2 < clust[ic].drv**2)[0]
  
      if plot :
        ax.cla() 
        plots.plotp(ax,vra,vdec,color='k',
                    xr=[med_vra-100,med_vra+200],xt='PMRA (km/sec at cluster dist)',
                    yr=[med_vdec-100,med_vdec+100],yt='PMDEC (km/sec at cluster dist)')
        plots.plotp(ax,vra[i2],vdec[i2],color='r',size=30)
        plots.plotp(ax,vra[i2[j]],vdec[i2[j]],color='g',size=30)
        fig.suptitle('{:s} PM (km/s): {:4.2f} +/- {:4.2f}  {:4.2f} +/ {:4.2f}'.format(cluster,med_vra,clust[ic].drv, med_vdec,clust[ic].drv))
        if hard is not None :
            fig.savefig(hard+'/'+clust[ic].name+'_pm.jpg')
            plt.close()
        else :
            pdb.set_trace()
      if len(j) > 0 :
        if pm: jc=jc[i1[j]]
      else :
        jc=[]
      print('{:d} stars after PM criterion'.format(len(jc)))
      if len(jc) <= 1 : return jc
   
      # parallaxes
      gaia=job.get_results()
      i1, i2 = match.match(np.core.defchararray.replace(data['APOGEE_ID'][jc],'2M',''),gaia['original_ext_source_id'])
      par=gaia['parallax']
      par_error=gaia['parallax_error']
      gd=np.where(np.isfinite(par[i2]))[0]
      med_par=np.median(par[i2[gd]])
      med_par_error=np.median(par_error[i2[gd]])
      j=np.where(np.isfinite(par[i2]) & (np.abs(par[i2]-med_par) < 3*med_par_error))[0]
      if plot :
        ax.cla() 
        ax.hist(par,color='k',bins=np.arange(par.min(),par.max(),0.01),histtype='step',range=(0,2))
        ax.hist(par[i2],color='r',bins=np.arange(par.min(),par.max(),0.01),histtype='step',range=(0,2))
        ax.hist(par[i2[j]],color='g',bins=np.arange(par.min(),par.max(),0.01),histtype='step',range=(0,2))
        ax.set_xlabel('Parallax')
        fig.suptitle('{:s} Parallax : {:4.2f} +/- {:4.2f} '.format(cluster,med_par, 3*med_par_error))
        if hard is not None :
            fig.savefig(hard+'/'+clust[ic].name+'_parallax.jpg')
            plt.close()
        else :
            plt.draw()
            pdb.set_trace()
      if len(j) > 0 :
        if dist: jc=jc[i1[j]]
      else :
        jc=[]
      print('{:d} stars after parallax criterion'.format(len(jc)))
      if len(jc) <= 1 : return jc

    # parameters criteria
    if raw :
        param='FPARAM'
    else :
        param='PARAM'
    try :
        j = np.where((data[param][jc,1] >= logg[0]) & (data[param][jc,1] <= logg[1]) &
                     (data[param][jc,0] >= te[0]) & (data[param][jc,0] <= te[1]) )[0]
        if len(j) > 0 :
            jc=jc[j]
        else :
            jc=[]
    except: pass
    print('{:d} stars after parameters criterion'.format(len(jc)))
    if len(jc) <= 1 : return jc

    # Remove badstars
    if plot :
        ax.cla()
        plots.plotp(ax,data['J'][jf]-data['K'][jf],data['K'][jf],color='k',size=20,xr=[-0.5,1.5],yr=[15,6],facecolors='none',linewidth=1,draw=False,xt='J-K',yt='K')
        plots.plotp(ax,data['J'][jc]-data['K'][jc],data['K'][jc],color='g',size=30,xr=[-0.5,1.5],yr=[15,6],draw=False)
    badstars = open(os.environ['APOGEE_DIR']+'/data/calib/badcal.dat')
    bad = []
    for line in badstars :
       bad.append(line.split()[0])
    jc = [x for x in jc if data[x]['APOGEE_ID'] not in bad]
    print('{:d} stars after badstars rejection'.format(len(jc)))
    if len(jc) <= 1 : return jc

    # remove non firstgen GC stars if requested
    if firstgen :
        gcstars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/gc_szabolcs.dat')
        if firstpos :
            gd=np.where(gcstars['pop'] == 1)[0]
            jc = [x for x in jc if data[x]['APOGEE_ID'] in gcstars['id'][gd]]
        else :
            bd=np.where(gcstars['pop'] != 1)[0]
            jc = [x for x in jc if data[x]['APOGEE_ID'] not in gcstars['id'][bd]]
    print('{:d} stars after firstgen rejection'.format(len(jc)))

    if plot :
        plots.plotp(ax,data['J'][jc]-data['K'][jc],data['K'][jc],color='b',size=30,draw=False)
        if hard is not None :
            fig.savefig(hard+'/'+clust[ic].name+'_cmd.jpg')
            plt.close()
        else :
            plt.draw()
            pdb.set_trace()

    return jc

def clusters(data,dir='clusters/') :
    """ Determine cluster members from input structure, make web pages with selection
    """
    try: os.mkdir(dir)
    except: pass
    # get the cluster data
    clusts=clustdata()
    # master output text file
    fstars=open(dir+'/allclust.txt','w')
    # HTML output
    f=html.head(file=dir+'/clust.html')
    f.write('<A HREF=allclust.txt> cluster stars list </a>')
    f.write('<TABLE BORDER=2>\n')
    f.write('<TR><TD>NAME<TD>RA<TD>DEC<TD>Radius<TD>RV<TD>Delta RV<TD>Position criterion<TD>RV criterion<TD>PM criterion<TD>Parallax criterion<TD> CMD')
    clust=clustdata()
    for ic in range(len(clust.name)) :
        print(clust[ic].name)
        j=clustmember(data,clust[ic].name,plot=True,hard=dir)
        print(clust[ic].name,len(j))
        # clusters to exclude here
        f.write('<TR><TD><A HREF='+clust[ic].name+'.txt>'+clust[ic].name+'</A><TD>{:12.6f}<TD>{:12.6f}<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                clust[ic].ra,clust[ic].dec,clust[ic].rad,clust[ic].rv,clust[ic].drv))
        f.write('<TD><A HREF='+clust[ic].name+'_pos.jpg><IMG SRC='+clust[ic].name+'_pos.jpg width=300></A>\n')
        f.write('<TD><A HREF='+clust[ic].name+'_rv.jpg><IMG SRC='+clust[ic].name+'_rv.jpg width=300></A>\n')
        f.write('<TD><A HREF='+clust[ic].name+'_pm.jpg><IMG SRC='+clust[ic].name+'_pm.jpg width=300></A>\n')
        f.write('<TD><A HREF='+clust[ic].name+'_parallax.jpg><IMG SRC='+clust[ic].name+'_parallax.jpg width=300></A>\n')
        f.write('<TD><A HREF='+clust[ic].name+'_cmd.jpg><IMG SRC='+clust[ic].name+'_cmd.jpg width=300></A>\n')
        np.savetxt(dir+'/'+clust[ic].name+'.txt',data[j]['APOGEE_ID'],fmt='%s')
        for star in data[j]['APOGEE_ID'] : fstars.write('{:s} {:s}\n'.format(star,clust[ic].name))
    html.tail(f)
    fstars.close()

