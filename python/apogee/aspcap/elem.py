# routines related to individual element calibration for APOGEE/ASPCAP

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
from apogee.utils import apload, apselect, bitmask
from apogee.aspcap import err, aspcap
from tools import plots, html, fit, match
import pdb
from astropy.io import fits, ascii
from astropy.table import Table, Column
try:
   import esutil
except:
   pass
import copy
import os

def read(file='allStar-testcal.fits') :
    '''
    Read allStar file, get main structure, elem_symbol, and elemtoh
    '''
    dr13load=apload.ApLoad(dr='dr13')
    #a=apload.allStar()[1].data
    c=dr13load.allStar()[3].data
    #a=fits.open('../dist/allStar+.fits')[1].data
    #x,y,z,r = galmodel.lbd2xyz(a['GLON'],a['GLAT'],a['DISO'][:,2]/1000.)
    #zone=np.where((r>9) & (r<11) & (dt<40))[0]

    a=fits.open(file)[1].data
    #c=fits.open(file)[3].data
    elem=c['ELEM_SYMBOL'][0]
    elemtoh=c['ELEMTOH'][0]

    return a, elem, elemtoh

def arctabun(el) :
    '''
    Define Arcturus abundances, and return requested abundance
    '''
    abun = { "C" : 0.090000, "CI" : 0.09, "N" : 0.400000, "O" : 0.480000, "Na" : 0.210000, "Mg" : 0.370000, "Al" : 0.400000, "Si" : 0.330000, "P" : 0.070000, "S" : 0.350000, "K" : 0.200000, "Ca" : 0.090000, "Sc" : 0.070000, "Ti" : 0.250000, "TiII" : 0.25, "V" : 0.160000, "Cr" : -0.050000, "Mn" : -0.120000, "Fe" : -0.000000, "Co" : 0.040000, "Ni" : 0.030000, "Cu" : -0.050000, "Ge" : 0.000000, "Rb" : 0.000000, "Y" : 0.000000, "Ce" : -0.190000, "Nd" : 0.130000, "Yb" : 0., "M" : 0., "alpha" : 0.3}
    return(abun[el]) 

def optabun(el) :
    '''
    ??? define abundance offsets from some optical analysis ???
    '''
    abun = {"Na" : -0.15, "Mg" :  0.06, "Al" :  0.04, "Si" : -0.21, "Ca" :  0.11, "Ti" : -0.14, "TiII" : 0.08, "V" : -0.15, "Cr" : -0.04, "Mn" : -0.36, "Fe" : 0.06, "Co" : -0.26}
    try :
        return(abun[el]) 
    except :
        return(-9999.)


def refabun(el,dwarf=False) :
   '''
   Return reference abundance: 0 if giant,  Arcturus if not?
   '''
   if dwarf :
       return 0.
   else :
       return arctabun(el)

def plot(a,elem,etoh,dwarf=False,suffix='',gcal=None,dcal=None,glon=None,glat=None,res=None,usemh=False,sn=[200,1000]) :
    '''
    Make a bunch of plots for elemental abundances
    '''

    try: os.mkdir('elem')
    except: pass

    # selection
    #dt=a['FPARAM'][:,0]-(4468+(a['FPARAM'][:,1]-2.5)/0.0018 - 382.5*a['FPARAM'][:,3])
    #gd=apselect.select(a[zone],badval='STAR_BAD',logg=[-1,3.5],sn=[200,1000],teff=[4000,4800])
    #gd=zone[gd]
    if dwarf :
        tit = 'Dwarfs, S/N>200'
        prefix = 'd'+suffix
        tmax=6500
        gd=apselect.select(a,badval='STAR_BAD',sn=sn,raw=True,glon=glon,glat=glat,dwarfs=True)
        etoh[0]=1
        etoh[1]=1
        etoh[2]=1
        ref=apselect.select(a,id='VESTA')
    else :
        tit = 'Giants, S/N>200'
        prefix = 'g'+suffix
        tmax=6500
        gd=apselect.select(a,badval='STAR_BAD',sn=sn,raw=True,glon=glon,glat=glat,giants=True)
        ref=apselect.select(a,id='alpha_Boo')
    out = open('elem/'+prefix+'.dat','w')


    # get the indices for different grids, and for stars near solar metallicity
    fgrid=apselect.select(a[gd],grid='F',raw=True)
    gkgrid=apselect.select(a[gd],grid='GK',raw=True)
    mgrid=apselect.select(a[gd],grid='M',raw=True)
    solar=apselect.select(a[gd],mh=[-0.1,0.1],raw=True)

    ifeh=17
    if len(a['FELEM'].shape) == 2: felem_feh = a['FELEM'][:,ifeh]
    else : felem_feh = a['FELEM'][:,0,ifeh]
  
    ytit=[]
    files=[]
    # loop over elements
    nelem=len(elem)
    for ielem in range(nelem+2) :
        file=[]
        if ielem < nelem :
            el = elem[ielem].strip()
            #eelem = a['ELEM'][gd,ielem]
            eelem = a['X_M'][gd,ielem]
            if len(a['FELEM'].shape) == 2: felem = a['FELEM'][gd,ielem]
            else : felem = a['FELEM'][gd,0,ielem]
            eelem_err = a['X_M_ERR'][gd,ielem]
            if len(a['FELEM'].shape) == 2: felem_err = a['FELEM_ERR'][gd,ielem]
            else: felem_err = a['FELEM_ERR'][gd,0,ielem]
            tmp=etoh[ielem]
            if ielem > 2 :
                if usemh and etoh[ielem] :
                    #eelem -= a['FPARAM'][gd,3]
                    felem -= a['FPARAM'][gd,3]
                elif not usemh and not etoh[ielem] :
                    #eelem += a['FPARAM'][gd,3]
                    felem += a['FPARAM'][gd,3]
            else :
                giants = apselect.select(a[gd],grid='g_',raw=True)
                if not usemh :
                    #eelem[giants] += a['FPARAM'][gd[giants],3] 
                    felem[giants] += a['FPARAM'][gd[giants],3] 
                dwarfs = apselect.select(a[gd],grid='d_',raw=True)
                if usemh :
                    #eelem[dwarfs] -= a['FPARAM'][gd[dwarfs],3] 
                    felem[dwarfs] -= a['FPARAM'][gd[dwarfs],3] 
        elif ielem == nelem :
            el = 'M'
            eelem = a['PARAM'][gd,0]
            felem = a['FPARAM'][gd,3]
            eelem_err = np.sqrt(a['PARAM_COV'][gd,3,3])
            felem_err = np.sqrt(a['FPARAM_COV'][gd,3,3])
            tmp = 1
        else :
            el = 'alpha'
            eelem = a['PARAM'][gd,6]
            felem = a['FPARAM'][gd,6]
            eelem_err = np.sqrt(a['PARAM_COV'][gd,6,6])
            felem_err = np.sqrt(a['FPARAM_COV'][gd,6,6])
            tmp = 0
            if not usemh :
                eelem += a['FPARAM'][gd,3]
                felem += a['FPARAM'][gd,3]

        if (tmp == 1 and not usemh) or (tmp == 0 and usemh ):
            refoffset=0
        else :
            refoffset=a['FPARAM'][ref,3]
            if usemh :
                refoffset *= -1

        name=prefix+el
        print(name)
        fname = 'elem/'+name
        # loop over plots
        xtit = []
        for iplot in range(0,8) :
        #for iplot in range(2,3) :
            if iplot == 0 :
              #x = a['ELEM'][gd,ifeh]
              x = a['X_H'][gd,ifeh]
              xr = [-1.5,1.]
              xt= '[Fe/H] (cal)'
              y = eelem
              if not usemh: y-=a['PARAM'][gd,3]
              yr=[-0.25,0.5]
              yt = '['+name+'/M](cal)'
              z = a['FPARAM'][gd,0]
              zr = [3000,tmax]
              zt='Teff'
              xtit.append('calibrated vs [Fe/H]')
            elif iplot == 1 :
              x = felem_feh[gd]
              xr = [-1.5,1.]
              xt= '[Fe/H] (raw)'
              y = felem
              if not usemh: y-=a['FPARAM'][gd,3]
              yr=[-0.25,0.5]
              yt = '['+name+'/M](raw)'
              z = a['FPARAM'][gd,0]
              zr = [3000,tmax]
              zt='Teff'
              xtit.append('raw vs [Fe/H]')
            elif iplot == 2 :
              x = a['FPARAM'][gd,0]
              xr = [2500,tmax]
              xt= 'Teff'
              y = eelem
              if not usemh: y-=a['PARAM'][gd,3]
              yr=[-0.25,0.5]
              yt = '['+name+'/M](cal)'
              #z = a['ELEM'][gd,ifeh]
              z = a['X_H'][gd,ifeh]
              zr = [-1.5,1.]
              zt='[Fe/H]'
              xtit.append('calibrated vs Teff')
            elif iplot == 3 :
              x = a['FPARAM'][gd,0]
              xr = [2500,tmax]
              xt= 'Teff'
              y = felem
              if not usemh: y-=a['FPARAM'][gd,3]
              yr=[-0.25,0.5]
              yt = '['+name+'/M](raw)'
              z = felem_feh[gd]
              zr = [-1.5,1.]
              zt='[Fe/H]'
              xtit.append('raw vs Teff')
            elif iplot == 4 :
              x = a['FPARAM'][gd,0]
              xr = [3000,tmax]
              xt = 'Teff'
              y = eelem-felem
              yr = [-0.3,0.3]
              yt = 'cal - raw'
              z = felem_feh[gd]
              zr = [-1.5,1.]
              zt='[Fe/H]'
              xtit.append('calibration')
            elif iplot == 5 :
              x = a['FPARAM'][gd,0]
              xr = [2500,tmax]
              xt = 'Teff'
              y = eelem_err
              yr= [0,0.3]
              yt = 'Empirical uncertainty'
              z = felem_feh[gd]
              zr = [-1.5,1.]
              zt='[Fe/H]'
              xtit.append('empirical uncertainty')
            elif iplot == 6 :
              x = a['FPARAM'][gd,0]
              xr = [2500,tmax]
              xt = 'Teff'
              y = felem_err
              yr= [0,0.3]
              yt = 'FERRE uncertainty'
              z = felem_feh[gd]
              zr = [-1.5,1.]
              zt='[Fe/H]'
              xtit.append('FERRE uncertainty')
            elif iplot == 7 :
              x = a['FPARAM'][gd,0]
              xr = [2500,tmax]
              xt = 'Teff'
              if ielem < nelem :
                y = a['ELEM_CHI2'][gd,ielem]
              else :
                y = x*0.
              yr= [0,50]
              yt = 'ELEM_CHI2'
              z = felem_feh[gd]
              zr = [-1.5,1.]
              zt='[Fe/H]'
              xtit.append('CHI2 from element fit')
    
            fig=plt.figure(figsize=(10,8))
            ax=fig.add_subplot(111)
            if len(x) > 0 :
                if len(fgrid) > 0 :
                    plots.plotc(ax,x[fgrid],y[fgrid],z[fgrid],xr=xr,yr=yr,zr=zr,colorbar=False,size=10,marker='s',yt=yt,xt=xt,zt=zt)
                if len(gkgrid) > 0 :
                    plots.plotc(ax,x[gkgrid],y[gkgrid],z[gkgrid],xr=xr,yr=yr,zr=zr,size=10,marker='o',yt=yt,xt=xt,zt=zt)
                if len(mgrid) > 0 :
                    plots.plotc(ax,x[mgrid],y[mgrid],z[mgrid],xr=xr,yr=yr,zr=zr,size=7,marker='^',yt=yt,xt=xt,zt=zt)
            if (iplot == 0 or iplot ==  2) : 
                if res is not None :
                  clust, = np.where(res['col3'] == ielem)
                  plots.plotp(ax,res['col4'][clust],res['col9'][clust],xr=xr,yr=yr,size=50,marker='o',facecolors='none',linewidth=1)

                # plot the reference star abundance (Arcturus or Vesta)
                if ielem < nelem-2 :
                  #refval = a['ELEM'][ref,ielem]+refoffset
                  #referr = a['ELEM_ERR'][ref,ielem]
                  refval = a['X_M'][ref,ielem]
                  referr = a['X_M_ERR'][ref,ielem]
                elif ielem == nelem-2 :
                  refval = a['PARAM'][ref,3]+refoffset
                  referr = np.sqrt(a['PARAM_COV'][ref,3,3])
                else :
                  refval = a['PARAM'][ref,6]+refoffset
                  referr = np.sqrt(a['PARAM_COV'][ref,6,6])
                if not usemh: refval -= a['PARAM'][ref,3]
                reflit = (refabun(el,dwarf=dwarf)-refabun('Fe',dwarf=dwarf))
                plots.plotl(ax,xr,[refval-reflit,refval-reflit],color='r')

                # Plot the median of solar abundance stars
                cal=np.where(y[solar] > -9000)[0]
                med = np.median(y[solar[cal]])
                plots.plotl(ax,xr,[med,med],color='y')
                plt.text(xr[0]+0.05*(xr[1]-xr[0]),yr[0]+0.85*(yr[1]-yr[0]),'solar metallicity stars: {:4.2f}'.format(med),color='y')
                if iplot == 0 :
                    out.write(el+'{:8.3f}  {:8d}\n'.format(med,len(cal)))

                # Plot the offset from the optical analysis 
                opt=optabun(el)
                plots.plotl(ax,xr,[opt,opt],color='m')
                plt.text(xr[0]+0.05*(xr[1]-xr[0]),yr[0]+0.75*(yr[1]-yr[0]),'optical offset: {:4.2f}'.format(opt),color='m')

                # Plot M67 points
                #clust=fits.open('../../cal/clust.fits')[1].data
                #m67=np.where(np.core.defchararray.find(clust['FIELD'],'M67') >= 0)
                #m1, m2 = esutil.numpy_util.match(a['APOGEE_ID'][gd],clust['APOGEE_ID'][m67])
                #plots.plotp(ax,x[m1],y[m1],size=6,color='k',xr=xr,yr=yr)
                #if iplot == 2 :
                #    print(m1, x[m1])
                if dwarf : 
                    refstar = 'VESTA'
                    if dcal is not None : 
                        refclust=dcal[ielem]  #-dcal[ifeh]
                        plots.plotl(ax,xr,[refclust,refclust],color='g')
                        plt.text(xr[0]+0.05*(xr[1]-xr[0]),yr[0]+0.95*(yr[1]-yr[0]),'M67 dwarfs: {:4.2f}'.format(refclust),color='g')
                else :
                    refstar = 'Arcturus'
                    if gcal is not None : 
                        refclust=gcal[ielem]  #-gcal[ifeh]
                        plots.plotl(ax,xr,[refclust,refclust],color='g')
                        plt.text(xr[0]+0.05*(xr[1]-xr[0]),yr[0]+0.95*(yr[1]-yr[0]),'M67 giants: {:4.2f}'.format(refclust),color='g')
                        if dcal is not None : 
                            drefclust=dcal[ielem] #-dcal[ifeh]
                            plots.plotl(ax,xr,[refclust-drefclust,refclust-drefclust],color='b')
                            plt.text(xr[0]+0.05*(xr[1]-xr[0]),yr[0]+0.90*(yr[1]-yr[0]),'M67 dwarfs: {:4.2f}'.format(drefclust),color='b')
                plt.text(xr[0]+0.05*(xr[1]-xr[0]),yr[0]+0.05*(yr[1]-yr[0]),
                    refstar+'   ASPCAP: {:4.2f}+/-{:4.2f}'.format(refval[0],referr[0])+'   lit: '+'{:4.2f}'.format(reflit),color='r')

            #plt.show()
            plt.savefig(fname+'{:1d}.png'.format(iplot))
            plt.close()
            file.append(name+'{:1d}.png'.format(iplot))

        ytit.append(name)

#        plt.show()
        files.append(file)
    out.close()
    html.htmltab(files,ytitle=ytit,file='elem/'+prefix+'elem.html',xtitle=xtit,header=tit)

def elemindex() :
    ''' Make the HTML pages (assumes plots have already been made) for individual elements '''

    a,elem,elemtoh = read()

    # loop over elements
    nelem=len(elem)
    for ielem in range(len(elem)+2) :
      if ielem < len(elem) :
          el = elem[ielem].strip()
      elif ielem == nelem :
            el = 'M'
      elif ielem == nelem+1 :
            el = 'alpha'

      ytit=[]
      files=[]
      for prefix in [ 'g','d' ] :
       for suffix in [ '', 'gal' ] :
        name=prefix+el
        file = [prefix+'1mh'+suffix+el+'2.png',prefix+'1emh'+suffix+el+'2.png',prefix+'1eemh'+suffix+el+'2.png',
                prefix+'2mh'+suffix+el+'2.png',prefix+'2emh'+suffix+el+'2.png',prefix+'2eemh'+suffix+el+'2.png',
                prefix+'3mh'+suffix+el+'2.png',prefix+'3emh'+suffix+el+'2.png',prefix+'3eemh'+suffix+el+'2.png']
        xtit = ['Linear 4000-5250','Linear 3750-5250','Linear 3500-5250',
               'Quadratic 4000-5250','Quadratic 3750-5250','Quadratic 3500-5250',
               'Cubic 4000-5250','Cubic 3750-5250','Cubic 3500-5250']
        files.append(file)

      ytit = ['Giants (full)','Giants (70&lt;l&lt;110)','Dwarfs (full)','Dwarfs (70&lt;l&lt;110)']
      html.htmltab(files,file='elem/'+el+'.html',xtitle=xtit,ytitle=ytit)

def main() :
    ''' Make series of plots and web pages for each calibration "type" '''    

    #files = ['testcal','testcal1mh','testcal2mh','testcal3mh',
    #         'testcal1emh','testcal2emh','testcal3emh', 
    #         'testcal1eemh','testcal2eemh','testcal3eemh'] 
    #dirs = ['../testcal','testcal1mh','testcal2mh','testcal3mh',
    #         'testcal1emh','testcal2emh','testcal3emh', 
    #         'testcal1eemh','testcal2eemh','testcal3eemh'] 
    #suffixes = ['','1mh','2mh','3mh','1emh','2emh','3emh','1eemh','2eemh','3eemh']
    files = ['l30e.2']
    dirs = ['../cal']
    suffixes = ['']
    for i in range(len(files)) :
      a,e,etoh = read(file='allStar-'+files[i]+'.fits')
      gcal = fits.open(dirs[i]+'/giantcal.fits')[2].data['ABUN'][0,:,17]
      dcal = fits.open(dirs[i]+'/dwarfcal.fits')[2].data['ABUN'][0,:,17]
      for d in [ False, True ] :
        if d :
          res = ascii.read(dirs[i]+'/dwarfcal.res')
        else :
          res = ascii.read(dirs[i]+'/giantcal.res')
        tmp=etoh   # since etoh gets changed for dwarfs
        plot(a,e,tmp,suffix=suffixes[i],dwarf=d,gcal=gcal,dcal=dcal,res=None,usemh=True,sn=[200,1000])
        plot(a,e,tmp,suffix=suffixes[i]+'gal',dwarf=d,gcal=gcal,dcal=dcal,glon=[70,110],glat=[-5,5],res=None,usemh=True)
    #a,e,etoh = read()
    #plot(a,e,etoh)
    #plot(a,e,etoh,dwarf=True)

def globalscatter(allstar,elems,vscatter=[0,0.2],pm=True,dist=True) :
    ''' 
    Compute scatter in clusters
    '''
    clust=apselect.clustdata()
    gd=apselect.select(allstar,badval='STAR_BAD',vscatter=vscatter)
    members=[]
    print('selecting')
    clusts = ['N2420', 'M67', 'N188', 'N7789', 'N6819', 'N6791']
    fp=open('global.dat','w')
    for cluster in clusts :
        j=np.array(apselect.clustmember(allstar[gd],str(cluster),param=None,pm=pm,dist=dist))
        print(cluster,len(j))
        members.append(j)
        for jj in j :
            fp.write('{:s} {:s} {:8.3f} {:8.1f} {:8.1f} {:8.1f} {:8.2f} {:s}\n'.format(
              cluster,allstar['APOGEE_ID'][gd[jj]],allstar['FE_H'][gd[jj]],allstar['TEFF'][gd[jj]],
              allstar['SNR'][gd[jj]],allstar['ASPCAP_CHI2'][gd[jj]],
              allstar['VSCATTER'][gd[jj]],allstar['STARFLAGS'][gd[jj]]))
    fp.close()

    iel=0
    nels=len(elems[0])+2
    fig,ax=plots.multi(2,int(round(nels/2.)),hspace=0.001,wspace=0.001,figsize=(8,10))
    plots.event(fig)
    plots._data=allstar
    plots._id_cols=['APOGEE_ID']

    color=['r','g','b','c','m','y']
    for iel,el in enumerate(np.append(elems,['M','alpha'])) :
        iclust=0
        all=np.array([])
        ix=iel%2
        iy=iel/2
        for cluster in clusts :
            i=np.where(clust.name == cluster)
            mh=clust[i].mh
            name=clust[i].name
            # get cluster members
            j=members[iclust]
            if len(j) > 0 :
                if el.strip() == 'Fe' :
                  abun=allstar['X_H'][gd[j],iel]
                  ok=np.where(((allstar['ELEMFLAG'][gd[j],iel] & 255) == 0) & (allstar['X_H_ERR'][gd[j],iel] < 0.2))[0]
                elif el.strip() == 'M' :
                  abun=allstar['M_H'][gd[j]]
                  ok=np.where(((allstar['PARAMFLAG'][gd[j],3] & 255) == 0) & (allstar['M_H_ERR'][gd[j]] < 0.2))[0]
                elif el.strip() == 'alpha' :
                  abun=allstar['ALPHA_M'][gd[j]]
                  ok=np.where(((allstar['PARAMFLAG'][gd[j],6] & 255) == 0) & (allstar['ALPHA_M_ERR'][gd[j]] < 0.2))[0]
                else :
                  abun=allstar['X_M'][gd[j],iel]
                  ok=np.where(((allstar['ELEMFLAG'][gd[j],iel] & 255) == 0) & (allstar['X_M_ERR'][gd[j],iel] < 0.2) & (allstar['X_M'][gd[j],iel] > -999) )[0]
                if len(ok) > 3 :
                    all=np.append(all,abun[ok]-abun[ok].mean())
                    plots.plotp(ax[iy,ix],allstar['TEFF'][gd[j[ok]]],abun[ok]-abun[ok].mean(),color=color[iclust],size=10,yr=[-0.5,0.5])

            iclust+=1
        print('{:s} {:10.3f} {:10.3f} {:d}\n'.format(el, all.mean(), all.std(), len(all)))
        ax[iy,ix].text(0.1,0.9,el.strip(),ha='left',va='top',transform=ax[iy,ix].transAxes)
        ax[iy,ix].text(0.9,0.9,'{:8.3f}'.format(all.std()),ha='right',va='top',transform=ax[iy,ix].transAxes)
        iel+=1

def getabun(data,elems,elemtoh,el,xh=False,terange=[-1,10000],calib=False,line=0) :
    '''
    Return the abundance of the requested element, given data array, elem array, element
    '''
    if calib :
        param = 'PARAM'
    else :
        param = 'FPARAM'
    parammask = bitmask.ParamBitMask()
    if el.strip() == 'M' :
        ok=np.where(((data['PARAMFLAG'][:,3] & parammask.badval()) == 0) & (data['FPARAM_COV'][:,3,3] < 0.2) &
                    ((data['ASPCAPFLAG']&aspcapmask.badval())== 0) &
                    (data['FPARAM'][:,0] >= terange[0]) & (data['FPARAM'][:,0] <= terange[1]) & (data[param][:,3] > -9990.) )[0]
        abun = data[param][:,3]
    elif el.strip() == 'alpha' :
        ok=np.where(((data['PARAMFLAG'][:,6] & parammask.badval()) == 0) & (data['FPARAM_COV'][:,6,6] < 0.2) &
                    ((data['ASPCAPFLAG']&aspcapmask.badval())== 0) &
                    (data['FPARAM'][:,0] >= terange[0]) & (data['FPARAM'][:,0] <= terange[1]) & (data[param][:,6] > -9990.) )[0]
        abun = data[param][:,6]
        if xh : abun+=data['FPARAM'][:,3]
    else :
        iel=np.where(np.core.defchararray.strip(elems) == el.strip())[0][0]
        if calib :
          if xh : 
              abun = data['X_H'][:,iel]
              abunerr = data['X_H_ERR'][:,iel]
          else :
              abun = data['X_M'][:,iel]
              abunerr = data['X_M_ERR'][:,iel]
        else :
          if len(data['FELEM'].shape) == 2: 
              abun = data['FELEM'][:,iel]
              abunerr = data['FELEM_ERR'][:,iel]
          else : 
              abun = data['FELEM'][:,line,iel]
              abunerr = data['FELEM_ERR'][:,line,iel]

          if xh and not elemtoh[iel] : abun+=data['FPARAM'][:,3]
          if not xh and elemtoh[iel] : abun-=data['FPARAM'][:,3]
          #if el.strip() == 'C' or el.strip() == 'CI' or el.strip() == 'N' :
          #  # special case for C and N for dwarfs, since those use [M/H] dimension
          #  try :
          #      dw = np.where((np.core.defchararray.find(data['ASPCAP_CLASS'],'GKd')>=0) | (np.core.defchararray.find(data['ASPCAP_CLASS'],'Fd')>=0)  |
          #                   (np.core.defchararray.find(data['ASPCAP_CLASS'],'Md')>=0))[0]
          #  except :
          #      dw = np.where((np.core.defchararray.find(data['CLASS'],'GKd')>=0) | (np.core.defchararray.find(data['CLASS'],'Fd')>=0)  |
          #                   (np.core.defchararray.find(data['CLASS'],'Md')>=0))[0]
          #  if xh : abun[dw]-=data['FPARAM'][dw,3]
          #  else : abun[dw]-=data['FPARAM'][dw,3]
        if calib : badflag = 255
        else : badflag = 0
       
        elemflag =  data['ELEMFLAG'][:,iel].astype(np.uint64)

        #try: ok=np.where(( (data['ELEMFLAG'][:,iel] & badflag) == 0) &
        try: ok=np.where(( (elemflag & badflag) == 0) &
                      (abunerr < 0.2) &
                      ((data['ASPCAPFLAG']&aspcapmask.badval())== 0) &
                      (data['FPARAM'][:,0] >= terange[0]) & 
                      (data['FPARAM'][:,0] <= terange[1]) & 
                      (abun > -9990.) )[0]
        except: pdb.set_trace()
    return abun, ok

def docal(allstar,elems,elemtoh,doels,xh=False,plot=True,sepplot=False,hard=None, maxvisit=100,calvers='default',dwarfs=False,inter=False,
        errpar=False,calib=False,nx=4,ny=2,maxvscatter=0.2,pm=True,dist=True, lines=False) :
    ''' 
    Determine internal calibration relations for elements
   
    Args:
        allstar : allStar-like HDUList
        elems : list of elems
        elemtoh : coresponding list of elemtoh code

    Keyword args:
        xh  : fit in [X/H]? (default=False, i.e. fit in [X/M])
        plot : show individual element plots
    '''

    # select cluster members from array that don't have STAR_BAD into data structure
    clusters=apselect.clustdata()
    calclusters=['M92','M15','M13','M3','M5','M12','M35','N2420','N188','M67','N7789','Pleiades','N6819','N6791',
                 'N6397','M55','N3201','N6752','N362','M4','N2808','47TUC']

    #calclusters=['N2420','N188','M67','N7789','Pleiades','N6819','N6791']
    errpar = False

    clusts = clusters.name
    types = np.arange(len(clusts))
    markers = np.zeros(len(clusts),dtype=str)
    colors = np.zeros(len(clusts),dtype=str)
    markers[np.where(clusters.mh > -999)[0]] = 's'
    markers[np.where(clusters.mh < -1)[0]] = 'o'
    markers[np.where(clusters.mh > 0)[0]] = '^'
    allcol=['r','g','b','c','m','y']
    for i in range(len(colors)) : colors[i] = allcol[i%6]
    if dwarfs : 
        logg=[3.8,5.5]
        reject=0.25
        glon=[0,360]
    else : 
        logg=[-1,3.8]
        reject=0.15
        glon=[70,110]
    #solar=apselect.select(allstar[1].data,badval='STAR_BAD',badtarg=['YOUNG','EMBEDDED','EMISSION','EXTENDED'],
    #                      raw=True,logg=logg,glon=glon,glat=[-5,5],sn=[200,10000])
    #solar=apselect.select(allstar[1].data,badval='STAR_BAD',badtarg=['YOUNG','EMBEDDED','EMISSION','EXTENDED'],
    #                      raw=True,logg=logg,sn=[200,10000],maxdist=500.)
    solar=apselect.select(allstar[1].data,badval='STAR_BAD',badtarg=['YOUNG','EMBEDDED','EMISSION','EXTENDED'],
                          raw=True,logg=logg,sn=[200,10000])
    try :
        gd=np.where((allstar[1].data['gaiaedr3_parallax_error'][solar]/abs(allstar[1].data['gaiaedr3_parallax'][solar]) < 0.1) )[0]   
        solar=solar[gd]
        distance = 1000./allstar[1].data['gaiaedr3_parallax'][solar]
        x,y,z,r=lbd2xyz(allstar[1].data['GLON'][solar],allstar[1].data['GLAT'][solar],distance/1000.)
        gd = np.where((abs(z) < 0.5) & (r>8) & (r<9))[0]
        solar=solar[gd]
    except:
        print('no distance information available for solar sample, using glon/glat')
        solar=apselect.select(allstar[1].data,badval='STAR_BAD',badtarg=['YOUNG','EMBEDDED','EMISSION','EXTENDED'],
                              raw=True,logg=logg,glon=glon,glat=[-5,5],sn=[200,10000])
    
    gd=apselect.select(allstar[1].data,badval='STAR_BAD',raw=True,logg=logg)
    print('ngd: ',len(gd))
    print('nsolar: ',len(solar))
    try :
        v=np.where(allstar[1].data['VISIT'][gd]<= maxvisit)[0]
        gd=gd[v]
    except :
        print('VISIT keyword does not exist')
    # preselect with fast HTM method at largest cluster radius
    #try:
    #    print('pre-selecting cluster members using HTM')
    #    h=esutil.htm.HTM()
    #    maxrad=clusters['rad'].max()
    #    m1,m2,rad=h.match(clusters['ra'],clusters['dec'],allstar[1].data['RA'][gd],allstar[1].data['DEC'][gd],maxrad,maxmatch=500)
    #    gd=gd[m2]
    #except :
    #    pass
    # now select per cluster
    print('selecting cluster members')
    all=[]
    for cluster in clusts :
        if cluster in calclusters :
            #clustdir=os.environ['APOGEE_REDUX']+'/r12/stars/junk/'
            #if clustdir :
            #    stars=ascii.read(clustdir+'/'+cluster+'.txt',names=['APOGEE_ID'],format='no_header')
            #    jsaved,j2 = match.match(allstar[1].data[gd]['APOGEE_ID'],stars['APOGEE_ID'])
            j=apselect.clustmember(allstar[1].data[gd],str(cluster),param=None,firstgen=True,firstpos=False,logg=logg,
                                   pm=pm,dist=dist)
            #print(cluster,len(j),len(jsaved))
            if len(j) < 1 :
                j=apselect.clustmember(allstar[1].data[gd],str(cluster),param=None,logg=logg,pm=pm,dist=dist)
            all=set(all).union(gd[j].tolist())
    data=allstar[1].data[list(all)]

    # in the abbreviated array, get the lists of cluster members
    members=[]
    fig,ax=plots.multi(1,1,figsize=(16,8))
    for label in ax.axes.get_xticklabels():
        label.set_visible(False)
    for label in ax.axes.get_yticklabels():
        label.set_visible(False)

    iplot=0
    for iclust,cluster in enumerate(clusts) :
        if cluster in calclusters :
            ax.scatter((iplot//12)*0.1+0.25,12-iplot%12,marker=markers[iclust],color=colors[iclust])
            ax.text((iplot//12)*0.1+0.26,12-iplot%12,clusts[iclust]+' ( '+str(clusters[iclust].mh)+')',color=colors[iclust],va='center')
            ax.set_xlim(0.23,0.8)
            j=apselect.clustmember(data,str(clusts[iclust]),param=None,firstgen=True,firstpos=False,logg=logg,pm=pm,dist=dist)
            if len(j) < 1 :
                j=apselect.clustmember(data,str(clusts[iclust]),param=None,logg=logg, pm=pm, dist=dist)
            iplot+=1
        else :
            j=[]
        # members is a list of lists of cluster members
        members.append(j)
    if hard is not None : 
        fig.savefig(hard+'clust_key.png')
        fig.savefig(hard+'clust_key.pdf')
        plt.close(fig)

    # setup output structured array
    rec = np.zeros(len(doels),dtype=[
                       ('elem','S5'),
                       ('elemfit','i4'),
                       ('mhmin','f4'),
                       ('te0','f4'),
                       ('temin','f4'),
                       ('temax','f4'),
                       ('femin','f4'),
                       ('femax','f4'),
                       ('caltemin','f4'),
                       ('caltemax','f4'),
                       ('extfit','i4'),
                       ('extpar','3f4'),
                       ('clust','{:1d}S16'.format(len(clusts))),
                       ('par','3f4'),
                       ('abun','{:1d}f4'.format(len(clusts))),
                       ('nstars','{:1d}i4'.format(len(clusts))),
                       ('mean','{:1d}f4'.format(len(clusts))),
                       ('rms','{:1d}f4'.format(len(clusts))),
                       ('rmsgd','{:1d}f4'.format(len(clusts))),
                       ('rawmean','{:1d}f4'.format(len(clusts))),
                       ('errpar','4f4'),
                       ])
    # empirical scatter bin setup: these are bin left edges
    if dwarfs :
        dmhbin=3.
        mhbins=np.arange(-2.25,0.75,dmhbin)
        nerrfit=2
        xr=[3000,7500]
    else :
        dmhbin=0.5
        mhbins=np.arange(-2.25,0.75,dmhbin)
        nerrfit=3
        xr=[3000,5500]
    dteffbin=250
    teffbins=np.arange(3500,6000,dteffbin)
    dsnbin=50
    snbins=np.arange(50,250,dsnbin)

    # plot setup
    if plot and not sepplot :
        fig,ax = plots.multi(nx,ny,hspace=0.001,wspace=0.5,figsize=(18,6))
    # plot setup for summary all-element plots
    if plot and len(doels) > 2 :
        nels=0
        for el in doels :
          # parameters for the fit for this element
          if calvers == 'dr13' :
              pars = dr13cal(el,dwarfs=dwarfs)
          elif calvers == 'dr14' :
              pars = dr14cal(el,dwarfs=dwarfs)
          elif calvers == 'dr16' :
              pars = dr16cal(el,dwarfs=dwarfs)
          else :
              pars = defaultcal(el,dwarfs=dwarfs)
          if pars['elemfit'] >=0 : nels+=1
        # to get all elements in plots 
        nels = len(doels)
        allfig,allax=plots.multi(2,(nels-1)//2+1,hspace=0.001,wspace=0.3,figsize=(12,18))
        if len(solar) > 0 : allsolarfig,allsolarax=plots.multi(2,(nels-1)//2+1,hspace=0.001,wspace=0.3,figsize=(12,18))
    if errpar :
        errfig,errax=plots.multi(len(snbins),len(doels),hspace=0.001,wspace=0.001,figsize=(3*len(snbins),2*len(doels)))

    # loop over all the elements!
    iel=0
    #iplot=0
    grid=[]
    yt=[]
    for iplot,el in enumerate(doels) :
        if lines :
            jelem = np.where(allstar[3].data['ELEM_SYMBOL'][0] == el)[0]
            nlines = len(np.where(allstar[3].data['FELEM_WIND'][0][0,:,jelem] > 0)[0])
            if nlines > 0 :
                linefig,lineax=plots.multi(2,nlines+1,hspace=0.001,wspace=0.4,figsize=(10,18))
        else : nlines = 0

        # parameters for the fit for this element
        if calvers == 'dr13' :
            pars = dr13cal(el,dwarfs=dwarfs)
        elif calvers == 'dr14' :
            pars = dr14cal(el,dwarfs=dwarfs)
        elif calvers == 'dr16' :
            pars = dr16cal(el,dwarfs=dwarfs)
        else :
            pars = defaultcal(el,dwarfs=dwarfs)
        pars['clust'] = np.array(clusts,dtype='S16')
        pars['abun'] = np.zeros(len(clusts))
        pars['par'] = np.zeros(3)
        pars['elem'] = el
        pars['errpar'] = np.zeros(4)
        elemfit = pars['elemfit']
        while elemfit >= 0 :
          # get the good abundance data for this element, load variables for fit (teff, abun, clust)
          abundata, ok = getabun(data,elems,elemtoh,el,xh=xh,calib=calib)
          snr=np.clip(data['SNR'],0.,snbins[-1]+dsnbin-0.1)
          print(el,pars['elemfit'],pars['mhmin'],len(ok))

          # get cluster members
          ind=np.array([],dtype=int)
          clust=np.array([],dtype='S16')
          apogee_id=np.array([],dtype='S16')
          jclust=[]
          for iclust,cluster in enumerate(clusts) :
              #if cluster in calclusters :
                  i=np.where(clusters.name == clusts[iclust])
                  # get cluster members: intersection of all cluster members and good ones for this element
                  j=list(set(ok).intersection(members[iclust]))
                  jclust.append(j)
                  if clusters[i].mh > pars['mhmin']  and len(j) > 3 :
                      # ind has the indices of all stars above the [M/H] threshold and good abundances
                      ind=np.append(ind,j)
                      clust=np.append(clust,[clusts[iclust]]*len(j))

          # loop if we have individual lines analysis
          for iline in range(1+nlines) :

            abundata, ok = getabun(data,elems,elemtoh,el,xh=xh,calib=calib,line=iline)

            teff=data['FPARAM'][ind,0]
            mh=data['FPARAM'][ind,3]
            vscatter=data['VSCATTER'][ind]
            abun=abundata[ind]
            try :
                visit=data['VISIT'][ind]
            except :
                visit = np.zeros(len(ind))
            # only use visits=0 and vscatter<maxvscatter[gd] for fit, but we'll plot all
            gd=np.where((visit == 0) & (vscatter<maxvscatter) & (teff>=pars['temin']) & (teff<=pars['temax']))[0]
            bd=np.where((visit > 0) | (vscatter>=maxvscatter) | (teff<pars['temin']) | (teff>pars['temax']))[0]
            if len(gd) > 2 :
                print(el,len(ind))
                for iter in range(2) :
                    print(el,iter,len(gd),pars['temin'],pars['temax'])
                    deriv=calderiv(teff[gd]-pars['te0'],abun[gd],clust[gd],order=pars['elemfit'])
                    soln,inv = fit.linear(abun[gd],deriv)
                    nclust = len(np.unique(clust[gd]))
                    pars['clust'] = np.sort(np.unique(clust[gd]))
                    pars['par'][0:pars['elemfit']] = soln[nclust:len(soln)]
                    pars['abun'] = soln[0:nclust]
                    func=calfunc(pars,teff,mh,abun,clust,order=pars['elemfit'],extcal=False)
                    if iter == 0 :
                        res=abun-func
                        gd=np.where((visit == 0) & (vscatter<maxvscatter) & (teff>=pars['temin']) & (teff<=pars['temax']) & (abs(res) <= reject))[0]
                        tmpreject=reject
                        while len(gd) < 10 and tmpreject<reject*8 :
                          tmpreject*=2.
                          gd=np.where((visit == 0) & (vscatter<maxvscatter) & (teff>=pars['temin']) & (teff<=pars['temax']) & (abs(res) <= tmpreject))[0]
        
                        bd=np.where((visit > 0) | (vscatter>=maxvscatter) | (teff<pars['temin']) | (teff>pars['temax']) | (abs(res) > tmpreject))[0]

                print('\nGlobal {:<8s} {:8.3f} (summed) {:8.3f} (with 3 visits)'.format(el, (abun[gd]-func[gd]).std(), (abun[bd]-func[bd]).std()))
                # loop through all clusters and determine mean and scatter for each cluster, and accumulate 
                #   data for scatter as f([M/H],Teff,S/N)
                print(' Clusters:  mean std (cal)  mean std (raw)')
                rmsdata=[]
                rmsderiv=[]
                if errpar and iline == 0 and hard is not None: 
                    f=open(hard+el.strip()+'_err_obj.dat','w')
                    fc=open(hard+el.strip()+'_err_clust.dat','w')
                tedata=[]
                sndata=[]
                mhdata=[]
                val=[]
                for iclust,cluster in enumerate(clusts) :
                  if cluster in calclusters and len(jclust[iclust])>3 :
                    j=np.array(jclust[iclust])
                    try:
                        cgd=np.where((data['VISIT'][j] == 0) & (data['VSCATTER'][j]<maxvscatter) & 
                                     (data['FPARAM'][j,0]>=pars['temin']) & (data['FPARAM'][j,0]<=pars['temax']))[0]
                    except:
                        cgd=np.where((data['VSCATTER'][j]<maxvscatter) & (data['FPARAM'][j,0]>=pars['temin']) & (data['FPARAM'][j,0]<=pars['temax']))[0]
                    if len(gd) > 1 :
                      rmsgd = (abundata[j[cgd]]-calfunc(pars,data['FPARAM'][j[cgd],0],data['FPARAM'][j[cgd],3],abundata[j[cgd]],''*len(j),order=pars['elemfit'])).std()
                    else :
                      rmsgd=-1.
                    rec['rms'][iel,iclust] = (abundata[j]-calfunc(pars,data['FPARAM'][j,0],data['FPARAM'][j,3],abundata[j],''*len(j),order=pars['elemfit'])).std()
                    rec['rmsgd'][iel,iclust] = rmsgd
                    rec['mean'][iel,iclust] = (abundata[j]-calfunc(pars,data['FPARAM'][j,0],data['FPARAM'][j,3],abundata[j],''*len(j),order=pars['elemfit'])).mean()
                    rec['rawmean'][iel,iclust] = abundata[j].mean()
                    rec['nstars'][iel,iclust] = len(j)
                    print('  {:<10s}{:8.3f}{:8.3f}{:8.3f}{:6d}{:6d}{:8.3f}{:8.3f}'.format(
                      clusts[iclust],rec['mean'][iel,iclust],rec['rms'][iel,iclust],rmsgd,rec['nstars'][iel,iclust],
                      len(cgd),abundata[j].mean(),abundata[j].std()))

                    # empirical uncertainties
                    if errpar and iline==0 :
                        tedata.extend(data['FPARAM'][j,0])
                        sndata.extend(snr[j])
                        mhdata.extend(data['FPARAM'][j,3])
                        val.extend(abundata[j]-rec['mean'][iel,iclust])
                        if hard is not None:
                          for jj in j : 
                              f.write('{:8.1f}{:8.2f}{:8.2f}{:8.3f}{:8.1f} {:s} {:s}\n'.format(
                                      data['FPARAM'][jj,0],snr[jj],data['FPARAM'][jj,3],abundata[jj]-rec['mean'][iel,iclust],
                                      rec['mean'][iel,iclust],clusts[iclust],data['APOGEE_ID'][jj]))
                        i=np.where(clusters.name == clusts[iclust])
                        for mhbin in mhbins :
                            if (clusters[i].mh > mhbin) and (clusters[i].mh <= mhbin+dmhbin) :
                              for teffbin in teffbins :
                                for snbin in snbins :
                                  ibin = np.where(( data['FPARAM'][j,0] > teffbin) & (data['FPARAM'][j,0] <= teffbin+dteffbin) &
                                                  ( snr[j] > snbin) & (snr[j] <= snbin+dsnbin) & (abs(abundata[j]-rec['mean'][iel,iclust]) < 0.3) )[0]
                                  if len(ibin) > 3 :
                                      if not np.isfinite(np.log(abundata[np.array(j)[ibin]].std())) : 
                                          pdb.set_trace()
                                      rmsdata.append(np.log(abundata[np.array(j)[ibin]].std()))
                                      if dwarfs :
                                        rmsderiv.append([1.,teffbin+dteffbin/2.-4500.,snbin+dsnbin/2.-100.])
                                      else :
                                        rmsderiv.append([1.,teffbin+dteffbin/2.-4500.,snbin+dsnbin/2.-100.,mhbin+dmhbin/2.])
                                      if hard is not None:
                                        fc.write('{:8.1f}{:8.2f}{:8.2f}{:8.2f}{:5d}{:8.3f} {:s}\n'.format(
                                              teffbin+dteffbin/2.,snbin+dsnbin/2.,mhbin+dmhbin/2.,clusters[i].mh[0],len(ibin),abundata[np.array(j)[ibin]].std(),clusts[iclust]))
                                      iplt = np.where(snbins == snbin)[0][0]
                                      try: plots.plotc(errax[iel,iplt],clusters[i].mh.tolist(),[teffbin+dteffbin/2.],[abundata[np.array(j)[ibin]].std()],
                                                    size=30,zr=[0,0.1],xr=[-2.5,0.5],yr=[3500,5500],linewidth=1)
                                      except: pdb.set_trace()
                   
                if errpar and iline==0 :
                    if hard is not None: 
                        f.close()
                        fc.close()
                    #empirical uncertainties
                    rmsdata=np.array(rmsdata)
                    rmsderiv=np.array(rmsderiv)
                    if len(rmsdata) > 5 :
                        soln,inv = fit.linear(rmsdata,rmsderiv.transpose())
                        y, x = np.mgrid[3500:5500:200j,-2.5:0.5:200j]
                        for iplt in range(len(snbins)) :
                              sn = snbins[iplt]+dsnbin/2.
                              errax[iel,iplt].imshow(elemerr(soln,y-4500.,sn-100.,x),extent=[-2.5,0.5,3500,5500], aspect='auto',vmin=0,vmax=0.1, origin='lower',cmap='rainbow')
                              errax[iel,iplt].text(0.98,0.98,el+' S/N={:4.0f}'.format(sn),va='top',ha='right',transform=errax[iel,iplt].transAxes)
 
                        pars['errpar'] = soln
                        # send all points to generic errfit function (not rms within each bin) for alternative approach and to get plots
                        try:
                            soln2 = err.errfit(np.array(tedata),np.array(sndata),np.array(mhdata),np.array(val),out=hard+el.strip(),mkhtml=False)
                            grid.append([os.path.basename(hard+el.strip()+'_err.png'),os.path.basename(hard+el.strip()+'_err_sn.png')])
                            yt.append(el.strip())
                        except: 
                            print('errfit failed: ',el)
 
                # get calibrated values before external calibration 
                print('getting calibrated values....')
                func_cal=calfunc(pars,teff,abun,mh,clust,order=pars['elemfit'],extcal=False)
                func_uncal=calfunc(pars,teff,abun,mh,clust,order=0,extcal=False)

                # get the abundances of the "solar circle" stars
                if len(solar) > 0 and len(doels) > 2 :
                    print('getting abundances of solar neighborhood stars....')
                    solar_teff=allstar[1].data['FPARAM'][solar,0]
                    solar_mh=allstar[1].data['FPARAM'][solar,3]
                    solar_abun,solar_ok= getabun(allstar[1].data[solar],elems,elemtoh,el,xh=xh,calib=calib)
                    solar_func=calfunc(pars,solar_teff,solar_mh,solar_abun,np.array(['']*len(solar_teff)),order=pars['elemfit'],calib=calib)
                    # get mean and scatter of solar metallicity stars, rejecting points more than 0.2 from mean
                    ss=np.where((solar_mh[solar_ok] > -0.05) & (solar_mh[solar_ok] < 0.05) & 
                                (solar_teff[solar_ok] > pars['temin']) & (solar_teff[solar_ok] < pars['temax']))[0]
                    median=np.median(solar_abun[solar_ok[ss]]-solar_func[solar_ok[ss]])
                    ss=np.where((solar_mh[solar_ok] > -0.05) & (solar_mh[solar_ok] < 0.05) & 
                                (solar_teff[solar_ok] > pars['temin']) & (solar_teff[solar_ok] < pars['temax']) &
                                (np.abs(solar_abun[solar_ok]-solar_func[solar_ok])<0.2))[0]
                    std=(solar_abun[solar_ok[ss]]-solar_func[solar_ok[ss]]).std()
                    if pars['extfit'] == 4 :
                        pars['extpar'] = np.array([median,0.,0.])
                    median_uncal=np.median(solar_abun[solar_ok[ss]])
                    std_uncal=solar_abun[solar_ok[ss]].std()
                if pars['extfit'] == 10 :
                    j=np.where(rec['nstars'][iel]>0)[0]
                    pars['extpar'][0] = np.median(rec['mean'][iel][j]-clusters[j].mh)
                elif pars['extfit'] == 11 :
                    j=np.where((clusters.mh < -1) & (rec['nstars'][iel]>0))[0]
                    pars['extpar'][0] = np.median(rec['mean'][iel][j]-clusters[j].mh)
                    j=np.where((clusters.mh > -0.5) & (rec['nstars'][iel]>0))[0]
                    pars['extpar'][1] = np.median(rec['mean'][iel][j]-clusters[j].mh)
 
                # make plots!
                print('making plots ...')
                if plot :
                    if sepplot :
                        fig,ax = plots.multi(nx,ny,hspace=0.001,wspace=0.5,figsize=[12,6])
                        fig1,ax1 = plots.multi(1,1,figsize=[12,4])
                        fig2,ax2 = plots.multi(1,1,figsize=[12,4])
                    else :
                        for iy in range(ny) :
                            for ix in range(nx) :
                                ax[iy,ix].cla()
                    if iline == 0 :
                        #after calibration
                        print('plots 1 ....')
                        plots.plotp(ax[0,0],teff[gd],abun[gd]-func_cal[gd], typeref=clust[gd],yr=[-0.29,0.29],xr=xr,
                                    types=clusts,color=colors,marker=markers,size=16,yt=el)
                        plots.plotp(ax[0,0],teff[bd],abun[bd]-func_cal[bd],typeref=clust[bd],yr=[-0.29,0.29],xr=xr,
                                types=clusts,color=colors,marker=markers,size=16,facecolors='none',linewidths=0.2)
                        ax[0,0].text(0.98,0.98,'{:5.3f}'.format((abun[gd]-func_cal[gd]).std()),transform=ax[0,0].transAxes,va='top',ha='right')
                        #before calibration
                        plots.plotp(ax[1,0],teff[gd],abun[gd]-func_uncal[gd],typeref=clust[gd],yr=[-0.29,0.29],xr=xr,
                                    types=clusts,color=colors,marker=markers,size=16,xt='Teff',yt=el)
                        plots.plotp(ax[1,0],teff[bd],abun[bd]-func_uncal[bd],typeref=clust[bd],yr=[-0.29,0.29],xr=xr,
                                    types=clusts,color=colors,marker=markers,size=16,facecolors='none',linewidths=0.2)
                        if sepplot:
                            plots.plotp(ax1,teff[gd],abun[gd]-func_uncal[gd],typeref=clust[gd],yr=[-0.29,0.29],xr=xr,
                                    types=clusts,color=colors,marker=markers,size=16,xt='Teff',yt=el)
                            plots.plotp(ax1,teff[bd],abun[bd]-func_uncal[bd],typeref=clust[bd],yr=[-0.29,0.29],xr=xr,
                                    types=clusts,color=colors,marker=markers,size=16,facecolors='none',linewidths=0.2)
                        ax[1,0].text(0.98,0.98,'{:5.3f}'.format((abun[gd]-func_uncal[gd]).std()),transform=ax[1,0].transAxes,va='top',ha='right')
                    # figure with all elements on same plot
                    if len(doels) > 2 :
                        print('plots 2 ....')
                        if iline == 0 :
                            plots.plotp(allax[iplot//2,iplot%2],teff[gd],abun[gd]-func_uncal[gd],typeref=clust[gd],yr=[-0.29,0.29],xr=[3500,5500],
                                        types=clusts,color=colors,marker=markers,size=8,xt='Teff',yt=el)
                            plots.plotp(allax[iplot//2,iplot%2],teff[bd],abun[bd]-func_uncal[bd],typeref=clust[bd],yr=[-0.29,0.29],xr=[3500,5500],
                                        types=clusts,color=colors,marker=markers,size=8,facecolors='none',linewidths=0.2)
                            allax[iplot//2,iplot%2].text(0.98,0.98,'{:5.3f}'.format(
                                        (abun[gd]-func_uncal[gd]).std()),transform=allax[iplot//2,iplot%2].transAxes,va='top',ha='right')
                            m67 = np.where(clusts == 'M67')[0][0]
                            allax[iplot//2,iplot%2].text(0.98,0.75,'{:5.3f}'.format(
                                        rec['rms'][iel,m67]),transform=allax[iplot//2,iplot%2].transAxes,va='top',ha='right',color='r')
                            allax[iplot//2,iplot%2].yaxis.set_major_locator(MultipleLocator(0.2))
                            allax[iplot//2,iplot%2].yaxis.set_minor_locator(MultipleLocator(0.05))
                            label = allax[iplot//2,iplot%2].yaxis.get_label()
                            if len(label.get_text()) < 5 : label.set_rotation(0)
                        if nlines > 0 :
                            plots.plotp(lineax[iline,0],teff[gd],abun[gd]-func_uncal[gd],typeref=clust[gd],yr=[-0.29,0.29],xr=xr,
                                        types=clusts,color=colors,marker=markers,size=16,xt='Teff',yt=el)
                            plots.plotp(lineax[iline,0],teff[bd],abun[bd]-func_uncal[bd],typeref=clust[bd],yr=[-0.29,0.29],xr=xr,
                                        types=clusts,color=colors,marker=markers,size=16,facecolors='none',linewidths=0.2)
                            plots.plotp(lineax[iline,1],teff[gd],abun[gd]-func_uncal[gd],typeref=clust[gd],yr=[-2.,2.],xr=xr,
                                        types=clusts,color=colors,marker=markers,size=16,xt='Teff',yt=el)
                            plots.plotp(lineax[iline,1],teff[bd],abun[bd]-func_uncal[bd],typeref=clust[bd],yr=[-2,2],xr=xr,
                                        types=clusts,color=colors,marker=markers,size=16,facecolors='none',linewidths=0.2)
                            if iline > 0 :
                                w=np.squeeze(allstar[3].data['FELEM_WIND'][0][:,iline-1,jelem])
                                lineax[iline,0].text(0.05,0.8,'{:8.2f}-{:8.2f}   {:8.2f}'.format(w[0],w[1],w[2]),transform=lineax[iline,0].transAxes,fontsize=10)
                                lineax[iline,1].text(0.05,0.8,'{:8.2f}-{:8.2f}   {:8.2f}'.format(w[0],w[1],w[2]),transform=lineax[iline,1].transAxes,fontsize=10)

                    # stuff for interactive plots
                    print('plots 3 ....')
                    plots._id_cols=['APOGEE_ID','VISIT']
                    plots._id_cols=['APOGEE_ID']
                    plots._data=data[ind]
                    plots._data_x=teff
                    plots._data_y=abun-func
                    # plot fits
                    x=np.linspace(pars['caltemin'],pars['caltemax'],200)
                    func=calfunc(pars,x,x*0.,x*0,np.array(['']*len(x)),order=pars['elemfit'],extcal=False)
                    plots.plotl(ax[1,0],x,func)
                    if sepplot: plots.plotl(ax1,x,func)
                    if len(doels) > 2 :
                        print('plots 4 ....')
                        # figure with all elements on same plot
                        if iline==0 : plots.plotl(allax[iplot//2,iplot%2],x,func)
                        # solar circle stars
                        if iline==0 and len(solar) > 0 :
                            plots.plotc(allsolarax[iplot//2,iplot%2],solar_teff[solar_ok],solar_abun[solar_ok]-solar_func[solar_ok],solar_mh[solar_ok],
                                        xr=xr,yr=[-0.5,0.5],zr=[-1,0.5],xt='Teff',yt=el)
                            plots.plotl(allsolarax[iplot//2,iplot%2],[pars['temin'],pars['temax']],[median,median],color='k')
                            plots.plotl(allsolarax[iplot//2,iplot%2],xr,[median,median],color='k',ls=':')
                            allsolarax[iplot//2,iplot%2].text(0.98,0.98,'{:5.3f}'.format(std),ha='right',va='top',transform=allsolarax[iplot//2,iplot%2].transAxes)
                            allsolarax[iplot//2,iplot%2].text(0.98,0.02,'{:5.3f}'.format(median),ha='right',va='bottom',transform=allsolarax[iplot//2,iplot%2].transAxes)
                            label = allsolarax[iplot//2,iplot%2].yaxis.get_label()
                            if len(label.get_text()) < 5 : label.set_rotation(0)
                            plots.plotc(ax[0,2],solar_teff[solar_ok],solar_abun[solar_ok]-solar_func[solar_ok],solar_mh[solar_ok],xr=xr,yr=[-0.5,0.5],zr=[-1,0.5])
                            plots.plotl(ax[0,2],xr,[median,median],color='orange')
                            ax[0,2].text(0.98,0.98,'{:5.3f}'.format(std),ha='right',va='top',transform=ax[0,2].transAxes)
                            ax[0,2].text(0.98,0.02,'{:5.3f}'.format(median),ha='right',va='bottom',transform=ax[0,2].transAxes)
                            plots.plotc(ax[0,3],solar_mh[solar_ok],solar_abun[solar_ok]-solar_func[solar_ok],solar_teff[solar_ok],yr=[-0.5,0.5],zr=xr)
                            #uncalibrated
                            plots.plotc(ax[1,2],solar_teff[solar_ok],solar_abun[solar_ok],solar_mh[solar_ok],xr=xr,yr=[-0.5,0.5],zr=[-1,0.5])
                            plots.plotl(ax[1,2],xr,[median_uncal,median_uncal],color='orange')
                            ax[1,2].text(0.98,0.98,'{:5.3f}'.format(std_uncal),ha='right',va='top',transform=ax[1,2].transAxes)
                            ax[1,2].text(0.98,0.02,'{:5.3f}'.format(median_uncal),ha='right',va='bottom',transform=ax[1,2].transAxes)
                            plots.plotc(ax[1,3],solar_mh[solar_ok],solar_abun[solar_ok],solar_teff[solar_ok],yr=[-0.5,0.5],zr=xr)
                    if xh or el == 'M' :
                      gdplt=np.where(rec['nstars'][iel]>0)[0]
                      plots.plotp(ax[0,1],clusters[gdplt].mh,rec['rawmean'][iel][gdplt]-clusters[gdplt].mh,
                                  typeref=clusters[gdplt].name,types=clusts,color=colors,marker=markers,size=16,
                                  xr=[-2.5,0.5],yr=[-0.6,0.6],xt='Lit [M/H]',yt='ASPCAP-lit [M/H]',yerr=rec['rms'][iel])
                      plots.plotp(ax[1,1],clusters[gdplt].mh,rec['mean'][iel][gdplt]-clusters[gdplt].mh,
                                  typeref=clusters[gdplt].name,types=clusts,color=colors,marker=markers,size=16,
                                  xr=[-2.5,0.5],yr=[-0.6,0.6],xt='Lit [M/H]',yt='ASPCAP-lit [M/H]',yerr=rec['rms'][iel])
                      if sepplot :
                          plots.plotp(ax2,clusters[gdplt].mh,rec['mean'][iel][gdplt]-clusters[gdplt].mh,
                                      typeref=clusters[gdplt].name,types=clusts,color=colors,marker=markers,size=16,
                                      xr=[-2.5,0.5],yr=[-0.6,0.6],xt='Lit [M/H]',yt='ASPCAP-lit [M/H]',yerr=rec['rms'][iel])
                          ax2.plot([-2.5,-1.0],[0.108797,0.108797],color='k')
                          ax2.plot([-1.0,-0.5],[0.108797,-0.0272657],color='k')
                          ax2.plot([-0.5,0.5],[-0.0272657,-0.0272657],color='k')


                    else :
                      gdplt=np.where(rec['nstars'][iel]>0)[0]
                      plots.plotp(ax[0,1],clusters[gdplt].mh,rec['rawmean'][iel][gdplt],
                                  typeref=clusters[gdplt].name,types=clusts,color=colors,marker=markers,size=16,
                                  xr=[-2.5,0.5],yr=[-0.6,0.6],xt='Lit [M/H]',yt='ASPCAP-lit [M/H]',yerr=rec['rms'][iel])
                      plots.plotp(ax[1,1],clusters[gdplt].mh,rec['mean'][iel][gdplt],
                                  typeref=clusters[gdplt].name,types=clusts,color=colors,marker=markers,size=16,
                                  xr=[-2.5,0.5],yr=[-0.6,0.6],xt='Lit [M/H]',yt='ASPCAP',yerr=rec['rms'][iel])
                    plots.event(fig)
                    #if iline == nlines : iplot+=1
                    #if not sepplot and cal != 'inter' : pdb.set_trace()
                    if iline == nlines and hard is not None : 
                        fig.savefig(hard+el.strip()+'.png')
                        plt.close(fig)
                        if sepplot: 
                            fig1.savefig(hard+el+'.pdf')
                            fig2.savefig(hard+el+'_lit.pdf')
                            plt.close(fig1)
                            plt.close(fig2)
                        if nlines > 0 : 
                            linefig.savefig(hard+el+'_lines.png')
                            linefig.savefig(hard+el+'_lines.pdf')
                            plt.close(linefig)     
           
            if inter :
                # with interactive options, can adjust fit order and limits and redo
                plt.draw()
                plt.pause(1)
                print('elemfit: ',elemfit)
                s = raw_input('enter new elemfit (-1 to continue to next element, l for new fit limits): ')
                try:
                    elemfit = int(s)
                except:
                    s = raw_input('enter new lower and upper fit limits in Teff: ')
                    pars['temin'] = int(s.split()[0])
                    pars['temax'] = int(s.split()[1])
                if elemfit >=0 : pars['elemfit'] = elemfit
            else :
                elemfit = -1
        # transfer results for this element to output summary array 
        for key in ['elem','elemfit','mhmin','te0','temin','temax','caltemin','caltemax','extfit','extpar','clust','abun','par','errpar'] :
            print(key, pars[key], pars['elem'], pars['elemfit'])
            if key == 'clust' or key == 'abun' or key == 'errpar':
                n=len(pars[key])
                rec[iel][key][0:n]=pars[key]
            elif key == 'par' :
                # reverse for aspcap_correct
                rec[iel][key][:]=pars[key][::-1]
            elif key == 'extpar' :
                print(pars[key])
                rec[iel][key][:]=pars[key][:]
            else :
                rec[iel][key]=pars[key]
        rec[iel]['femin'] = -99.999
        rec[iel]['femax'] = 99.999
        iel+=1
    #if plot and iplot%2 == 1 : 
    if plot and len(doels)%2 == 1 : 
        allax[iplot//2,iplot%2].set_visible(False)
        ticklabels = allax[iplot//2-1,iplot%2].get_xticklabels()
        plt.setp(ticklabels, visible=True)
 
    if plot and hard is not None and len(doels) > 2: 
        allfig.savefig(hard+'all.png')
        if len(solar) > 0 : allsolarfig.savefig(hard+'allsolar.png')
    if errpar and hard is not None :
        try: html.htmltab(grid,ytitle=yt,file=hard+'err_all.html')
        except: pass
        errfig.savefig(hard+'err_all.png')
        plt.close(errfig)

    return rec

def calfunc(pars,teff,mh,abun,clust,order=1,calib=False,extcal=True) :
    '''
    Apply calibration function. If clust is not '', then include the mean abundance for the cluster as determined from the fit,
    otherwise only apply the temperature correction

    '''
    npts=len(teff)
    func=np.zeros([npts])
    # if we are given clusters that are not part of the calibration, set them to -999
    j=np.where(clust != '')[0]
    func[j]=-999.
    # start with the cluster mean abundances if requested
    for iclust in range(len(pars['clust'])) :
        j=np.where(clust == pars['clust'][iclust].strip())[0]
        func[j] = pars['abun'][iclust]
    # add the temperature terms, truncating at temin and temax
    if calib == False :
      if order >= 1:
        temp=copy.copy(teff)
        bd=np.where(temp < pars['temin'])[0]
        temp[bd]=pars['temin']
        bd=np.where(temp > pars['temax'])[0]
        temp[bd]=pars['temax']
        for iorder in range(0,order) :
            func += pars['par'][iorder]*(temp-pars['te0'])**(iorder+1)
      if extcal :
        if pars['extfit'] == 4 :
          func += pars['extpar'][0]
        elif pars['extfit'] == 10 :
          func += pars['extpar'][0]+pars['extpar'][1]*mh+pars['extpar'][2]*mh**2
        elif pars['extfit'] == 11 :
          mhclip=np.clip(mh,-1.,-0.5)
          func += pars['extpar'][0] + (mhclip-(-1.))*(pars['extpar'][1]-pars['extpar'][0])/0.5
    return func

def calderiv(teff,abun,clust,order=1) :
    '''
    Function/derivatives for abundance calibration
    '''
    uclust=np.sort(np.unique(clust))
    npar=order+len(uclust)
    npts=len(teff)
    deriv=np.zeros([npar,npts])
    for iclust in range(len(uclust)) :
        j=np.where(clust == uclust[iclust])[0]
        deriv[iclust,j] = 1.
    if order >= 1:
        for iorder in range(0,order) :
            deriv[len(uclust)+iorder,:] = teff**(iorder+1)
    return deriv
        

def defaultcal(el,dwarfs=False) :
    '''
    Return default parameters for abundance calibrtion
    '''
    te0=4500
    temin=4000
    if dwarfs : temax=6000
    else : temax=5000
    elemfit=1
    extfit=0
    caltemin=3532.5
    caltemax=6500
    extpar=[0.,0.,0.]
    mhmin=-1
    return {'elemfit': elemfit, 'mhmin' : mhmin, 'te0': te0, 'temin': temin, 'temax': temax, 
            'caltemin': caltemin, 'caltemax' : caltemax, 'extfit' : extfit, 'extpar' : np.array(extpar)}

def dr16cal(el,dwarfs=False) :
    '''
    Return default parameters for abundance calibrtion
    '''
    te0=4500
    # values for WARN and to use for fits, if any
    temin=0
    if dwarfs : temax=100000
    else : temax=10000
    # default method/order for fit with Teff (0=none)
    elemfit=0
    # default method for zeropoint (4=solar neighborhood)
    extfit=4
    # values for BAD, i.e. no calibration
    caltemin=3032.5
    caltemax=7500
    extpar=[0.,0.,0.]
    # minimum metallicity to use in clusters
    mhmin=-1

    if el.strip() == 'Ge' : elemfit=-1
    if el.strip() == 'Rb' : elemfit=-1
    if el.strip() == 'Nd' : elemfit=-1
    if el.strip() == 'Yb' : elemfit=-1

    if not dwarfs :
        if el.strip() == 'C' : 
            extfit=0
        if el.strip() == 'CI' : 
            extfit=0
        if el.strip() == 'N' : 
            extfit=0
        if el.strip() == 'O' : 
            temax=5000
        if el.strip() == 'Na' : 
            temin=3750
        elif el.strip() == 'Al' : 
            temin=3400
        elif el.strip() == 'K' : 
            temin=3900
        elif el.strip() == 'P' : 
            temax=6000
        elif el.strip() == 'Ti' : 
            temin=4200
        elif el.strip() == 'TiII' : 
            temin=4000
        elif el.strip() == 'V' : 
            temax=4800
        elif el.strip() == 'Mn' : 
            temin=4000
        elif el.strip() == 'Fe' : 
            extfit=0
        elif el.strip() == 'Co' : 
            temin=3300
            temax=6500
        elif el.strip() == 'Cu' : 
            temin=4000
        elif el.strip() == 'Ce' : 
            temin=4000
            temax=5000

    else :

        if el.strip() == 'O' : 
            temax=5000
        elif el.strip() == 'Na' : 
            temin=5500
            temax=5500
        elif el.strip() == 'P' : 
            temin=4300
            temin=5500
            temax=5500
        elif el.strip() == 'S' : 
            temin=4260
        elif el.strip() == 'K' : 
            temin=4000
            temax=6500
        elif el.strip() == 'Ti' : 
            temin=4000
            temax=6000
        elif el.strip() == 'TiII' : 
            temin=5500
            temax=6000
        elif el.strip() == 'V' : 
            temin=4800
            temax=5500
        elif el.strip() == 'Cr' : 
            temin=3800
            temax=6200
        elif el.strip() == 'Mn' : 
            temin=3800
        elif el.strip() == 'Fe' : 
            extfit=0
        elif el.strip() == 'Co' : 
            temax=6500
        elif el.strip() == 'Cu' : 
            temax=6200
        elif el.strip() == 'Ce' : 
            temin=4200
            temin=5500
            temax=5500

    return {'elemfit': elemfit, 'mhmin' : mhmin, 'te0': te0, 'temin': temin, 'temax': temax, 
            'caltemin': caltemin, 'caltemax' : caltemax, 'extfit' : extfit, 'extpar' : np.array(extpar)}

def dr14cal(el,dwarfs=False) :
    '''
    Return calibration parameters for requested element for DR14 choices

    elemfit gives order/type of polynomial in cluster fit: 1 (linear), 2 (quadratic), 3 (cubic)
    temin/temax gives range over which fit is performed
    caltemin/caltemax gives range over which calibration can be applied (bad outside range)
    extfit gives source of external calibration:  1 (Arcturus), 2 (Vesta), 3 (M67), 4 (solar sequence), 10 (quadratic fit to clusters), 11(piecewise fit to clusters)
    extpar gives the values of the external calibration
    '''

    # defaults
    te0=4500
    temin=4000
    temax=5000
    elemfit=1
    extfit=0
    caltemin=3532.5
    caltemax=6500
    extpar=[0.,0.,0.]
    mhmin=-1
    if el.strip() == 'Ca' : mhmin = -2.
    if el.strip() == 'C' : mhmin = -0.6
    if el.strip() == 'Fe' : mhmin = -3.
    if el.strip() == 'K' : mhmin = -0.6
    if el.strip() == 'Mn' : mhmin = -2.0
    if el.strip() == 'Na' : mhmin = -0.6
    if el.strip() == 'Ni' : mhmin = -3.0
    if el.strip() == 'N' : mhmin = -0.6
    if el.strip() == 'O' : mhmin = -0.6
    if el.strip() == 'Si' : mhmin = -3.0
    if el.strip() == 'V' : mhmin = -0.6

    # nothing below -1
    if mhmin < -1 : mhmin=-1.

    if not dwarfs :
        # calibration parameters for giants
        if el.strip() == 'C' :
            elemfit= 0
        elif el.strip() == 'CI' :
            elemfit= 0
        elif el.strip() == 'N' :
            elemfit= 0
        elif el.strip() == 'O' :
            elemfit= 2
            temin= 3750
            extfit= 4
        elif el.strip() == 'Na' :
            elemfit= 1
            extfit= 4
        elif el.strip() == 'Mg' :
            elemfit= 1
            temax= 5250
            extfit= 4
        elif el.strip() == 'Al' :
            elemfit= 1
            extfit= 4
        elif el.strip() == 'Si' :
            elemfit= 2
            temin= 3750
            temax= 5250
            extfit= 4
        elif el.strip() == 'P' :
            elemfit= 1
            extfit= 4
        elif el.strip() == 'S' :
            elemfit= 1
            extfit= 4
        elif el.strip() == 'K' :
            elemfit= 1
            extfit= 4
        elif el.strip() == 'Ca' :
            elemfit= 1
            temin= 3750
            extfit= 4
        elif el.strip() == 'Ti' :
            elemfit= 1
            temin= 3750
            extfit= 4
        elif el.strip() == 'TiII' :
            elemfit= 1
            extfit= 4
        elif el.strip() == 'V' :
            elemfit= 1
            extfit= 4
        elif el.strip() == 'Cr' :
            elemfit= 1
            temin= 3750
            extfit= 4
        elif el.strip() == 'Mn' :
            elemfit= 1
            extfit= 4
        elif el.strip() == 'Fe' :
            elemfit= 1
            extfit= 4
        elif el.strip() == 'Co' :
            elemfit= 1
            extfit= 4
        elif el.strip() == 'Ni' :
            elemfit= 1
            extfit= 4
        elif el.strip() == 'Cu' :
            elemfit= -1
        elif el.strip() == 'Ge' :
            elemfit= -1
        elif el.strip() == 'Ce' :
            elemfit= -1
        elif el.strip() == 'Rb' :
            elemfit= -1
            #elemfit= 1
            #extfit= 4
        elif el.strip() == 'Y' :
            elemfit= -1
        elif el.strip() == 'Nd' :
            elemfit= -1
        elif el.strip() == 'M' :
            elemfit= 1
            extfit=11
        elif el.strip() == 'alpha' :
            elemfit= 1
            temax=5250
            extfit= 4
    else :

        # default values for dwarfs
        temin=3200
        temax=6250
        elemfit=3
        caltemin=-1
        caltemax=6500
        extfit=0
        extpar=[0.,0.,0.]


        # manual overrides for each element, dwarfs
        if el.strip() == 'C' :
            elemfit=1
            extfit=4
        elif el.strip() == 'CI' :
            elemfit=1
            caltemin=3500
            caltemax=5000
            extfit=4
        elif el.strip() == 'N' :
            elemfit=0
            caltemin=3500
            extfit=4
        elif el.strip() == 'O' :
            elemfit=2
            temin=3500
            temax=4500
            extfit=4
        elif el.strip() == 'Na' :
            elemfit=-1 #0
            temin=3750
            temax=5500
            caltemin=3750
            extfit=4
        elif el.strip() == 'Mg' :
            elemfit=1
            temin=3750
            extfit=4
        elif el.strip() == 'Al' :
            elemfit=2
            temin=3750
            caltemin=3500
            extfit=4
        elif el.strip() == 'Si' :
            elemfit=1
            temin=3500
            extfit=4
        elif el.strip() == 'P' :
            elemfit=0
            caltemin=3750
            caltemax=5000
            extfit=0
        elif el.strip() == 'S' :
            elemfit=1
            temin=3750
            caltemin=3532
            extfit=4
        elif el.strip() == 'K' :
            elemfit=2
            temin=3750
            caltemin=3750
            extfit=4
        elif el.strip() == 'Ca' :
            elemfit=1
            temin=3750
            caltemin=3750
            extfit=4
        elif el.strip() == 'Ti' :
            elemfit=3
            temin=3750
            temax=5250
            caltemin=3750
            extfit=4
        elif el.strip() == 'TiII' :
            elemfit=-1
            caltemax=-1
            extfit=0
        elif el.strip() == 'V' :
            elemfit=2
            temax=5250
            caltemin=3750
            extfit=4
        elif el.strip() == 'Cr' :
            elemfit=0
            temax=5250
            caltemin=3750
            extfit=4
        elif el.strip() == 'Mn' :
            elemfit=3
            temin=3500
            caltemin=3500
            extfit=4
        elif el.strip() == 'Fe' :
            elemfit=2
            temin=3500
            extfit=4
        elif el.strip() == 'Co' :
            elemfit=-1
        elif el.strip() == 'Ni' :
            elemfit=1
            temin=3500
            caltemin=3500
            extfit=4
        elif el.strip() == 'Cu' :
            elemfit=-1  #2
            temin=3750
            caltemin=3750
            extfit=4
        elif el.strip() == 'Ge' :
            elemfit=-1
        elif el.strip() == 'Ce' :
            elemfit=-1
        elif el.strip() == 'Rb' :
            elemfit=-1 #1
            caltemin=3500
            temin=3200
            temax=5250
            extfit=4
        elif el.strip() == 'Y' :
            elemfit=-1
        elif el.strip() == 'Nd' :
            elemfit=-1
        elif el.strip() == 'M' :
            elemfit=1
            temin=3200
            extfit=10
        elif el.strip() == 'alpha' :
            elemfit=2
            temin=3500
            caltemin=3500
            extfit=4
        
    return {'elemfit': elemfit, 'mhmin' : mhmin, 'te0': te0, 'temin': temin, 'temax': temax, 
            'caltemin': caltemin, 'caltemax' : caltemax, 'extfit' : extfit, 'extpar' : np.array(extpar)}

def dr13cal(el,dwarfs=False) :
    '''
    Return calibration parameters for requested element for DR13 choices

    elemfit gives order/type of polynomial in cluster fit: 1 (linear), 2 (quadratic), 3 (cubic)
    temin/temax gives range over which fit is performed
    caltemin/caltemax gives range over which calibration can be applied (bad outside range)
    extfit gives source of external calibration:  1 (Arcturus), 2 (Vesta), 3 (M67), 4 (solar sequence), 10 (quadratic fit to clusters)
    extpar gives the values of the external calibration
    '''

    # defaults
    te0=4500
    temin=4000
    temax=5000
    elemfit=1
    extfit=0
    caltemin=3532.5
    caltemax=6500
    extpar=[0.,0.,0.]
    mhmin=-1
    if el.strip() == 'Ca' : mhmin = -2.
    if el.strip() == 'C' : mhmin = -0.6
    if el.strip() == 'Fe' : mhmin = -3.
    if el.strip() == 'K' : mhmin = -0.6
    if el.strip() == 'Mn' : mhmin = -2.0
    if el.strip() == 'Na' : mhmin = -0.6
    if el.strip() == 'Ni' : mhmin = -3.0
    if el.strip() == 'N' : mhmin = -0.6
    if el.strip() == 'O' : mhmin = -0.6
    if el.strip() == 'Si' : mhmin = -3.0
    if el.strip() == 'V' : mhmin = -0.6

    # nothing below -1
    if mhmin < -1 : mhmin=-1.

    if not dwarfs :
        # calibration parameters for giants
        if el.strip() == 'C' :
            elemfit= 0
        elif el.strip() == 'CI' :
            elemfit= 0
        elif el.strip() == 'N' :
            elemfit= 0
        elif el.strip() == 'O' :
            elemfit= 2
            temin= 3750
            extfit= 4
            extpar= [0.060,0.,0.]
        elif el.strip() == 'Na' :
            elemfit= 2
            extfit= 4
            extpar= [0.186,0.,0.]
        elif el.strip() == 'Mg' :
            elemfit= 3
            temin= 3500
            extfit= 4
            extpar= [0.045,0.,0.]
        elif el.strip() == 'Al' :
            elemfit= 3
            extfit= 4
            extpar= [0.108,0.,0.]
        elif el.strip() == 'Si' :
            elemfit= 3
            temin= 3500
            extfit= 4
            extpar= [0.107,0.,0.]
        elif el.strip() == 'P' :
            elemfit= 2
            extfit= 4
            extpar= [-0.008,0.,0.]
        elif el.strip() == 'S' :
            elemfit= 2
            extfit= 4
            extpar= [-0.092,0.,0.]
        elif el.strip() == 'K' :
            elemfit= 1
            extfit= 4
            extpar= [-0.026,0.,0.]
        elif el.strip() == 'Ca' :
            elemfit= 3
            temin= 3750
            extfit= 4
            extpar= [-0.021,0.,0.]
        elif el.strip() == 'Ti' :
            elemfit= 3
            temin= 3500
            extfit= 4
            extpar= [-0.014,0.,0.]
        elif el.strip() == 'TiII' :
            elemfit= 2
            extfit= 4
            extpar= [0.166,0.,0.]
        elif el.strip() == 'V' :
            elemfit= 3
            temin= 3750
            extfit= 4
            extpar= [0.110,0.,0.]
        elif el.strip() == 'Cr' :
            elemfit= 2
            temin= 3500
            extfit= 4
            extpar= [-0.057,0.,0.]
        elif el.strip() == 'Mn' :
            elemfit= 1
            extfit= 4
            extpar= [0.041,0.,0.]
        elif el.strip() == 'Fe' :
            elemfit= 2
            temin= 3500
            extfit= 4
            extpar= [-0.005,0.,0.]
        elif el.strip() == 'Co' :
            elemfit= 3
            extfit= 4
            extpar= [0.003,0.,0.]
        elif el.strip() == 'Ni' :
            elemfit= 2
            temin= 3750
            extfit= 4
            extpar= [-0.001,0.,0.]
        elif el.strip() == 'Cu' :
            elemfit= 3
            temin= 3
            extfit= 4
            extpar= [0.452,0.,0.]
        elif el.strip() == 'Ge' :
            elemfit= 2
            extfit= 4
            extpar= [0.354,0.,0.]
        elif el.strip() == 'Ce' :
            elemfit= -1
        elif el.strip() == 'Rb' :
            elemfit= 2
            temin= 3750
            extfit= 4
            extpar= [-0.105,0.,0.]
        elif el.strip() == 'Y' :
            elemfit= -1
        elif el.strip() == 'Nd' :
            elemfit= -1
        elif el.strip() == 'M' :
            elemfit= 1
        elif el.strip() == 'alpha' :
            elemfit= 2
            extfit= 4
            extpar = [0.056,0.,0.]
    else :

        # default values for dwarfs
        temin=3200
        temax=6250
        elemfit=3
        caltemin=-1
        caltemax=6500
        extfit=0
        extpar=[0.,0.,0.]


        # manual overrides for each element, dwarfs
        if el.strip() == 'C' :
            elemfit=1
            extfit=4
            extpar=[-0.019,0.,0.]
        elif el.strip() == 'CI' :
            extfit=4
            extpar=[-0.026,0.,0.]
        elif el.strip() == 'N' :
            extfit=4
            extpar=[-0.01,0.,0.]
        elif el.strip() == 'O' :
            elemfit=3
            temin=3500
            temax=4500
            extfit=4
            extpar=[0.068,0.,0.]
        elif el.strip() == 'Na' :
            elemfit=1
            temin=3750
            temax=5500
            caltemin=3750
            extfit=4
            extpar=[0.096,0.,0.]
        elif el.strip() == 'Mg' :
            elemfit=3
            temin=3750
            extfit=4
            extpar=[-0.003,0.,0.]
        elif el.strip() == 'Al' :
            elemfit=2
            temin=3750
            caltemin=3500
            extfit=4
            extpar=[0.043,0.,0.]
        elif el.strip() == 'Si' :
            elemfit=1
            temin=3500
            extfit=4
            extpar=[-0.023,0.,0.]
        elif el.strip() == 'P' :
            caltemax=-1
            extfit=0
            extpar=[0.,0.,0.]
        elif el.strip() == 'S' :
            elemfit=1
            temin=3750
            caltemin=5500
            extfit=4
            extpar=[-0.017,0.,0.]
        elif el.strip() == 'K' :
            elemfit=2
            temin=3750
            caltemin=3750
            extfit=4
            extpar=[-0.029,0.,0.]
        elif el.strip() == 'Ca' :
            elemfit=1
            temin=3750
            caltemin=3750
            extfit=4
            extpar=[0.023,0.,0.]
        elif el.strip() == 'Ti' :
            elemfit=3
            temin=3750
            temax=5250
            caltemin=3750
            extfit=4
            extpar=[-0.002,0.,0.]
        elif el.strip() == 'TiII' :
            caltemax=-1
            extfit=0
            extpar=[0.,0.,0.]
        elif el.strip() == 'V' :
            elemfit=2
            temax=5250
            caltemin=3750
            extfit=4
            extpar=[0.002,0.,0.]
        elif el.strip() == 'Cr' :
            elemfit=1
            temax=5250
            caltemin=3750
            extfit=4
            extpar=[-0.044,0.,0.]
        elif el.strip() == 'Mn' :
            elemfit=3
            temin=3500
            caltemin=3500
            extfit=4
            extpar=[-0.077,0.,0.]
        elif el.strip() == 'Fe' :
            elemfit=2
            temin=3500
            extfit=4
            extpar=[0.016,0.,0.]
        elif el.strip() == 'Co' :
            elemfit=-1
        elif el.strip() == 'Ni' :
            elemfit=1
            temin=3500
            caltemin=3500
            extfit=4
            extpar=[0.03,0.,0.]
        elif el.strip() == 'Cu' :
            elemfit=2
            temin=3750
            caltemin=4500
            extfit=4
            extpar=[0.026,0.,0.]
        elif el.strip() == 'Ge' :
            elemfit=-1
        elif el.strip() == 'Ce' :
            elemfit=-1
        elif el.strip() == 'Rb' :
            elemfit=1
            temin=3200
            temax=5250
            extfit=4
            extpar=[-0.217,0.,0.]
        elif el.strip() == 'Y' :
            elemfit=-1
        elif el.strip() == 'Nd' :
            elemfit=-1
        elif el.strip() == 'M' :
            elemfit=3
            temin=3200
            extfit=0
            extpar=[0.0,0.,0.]
        elif el.strip() == 'alpha' :
            elemfit=1
            extfit=4
            extpar=[-0.004,0.,0.]
        
    return {'elemfit': elemfit, 'mhmin' : mhmin, 'te0': te0, 'temin': temin, 'temax': temax, 
            'caltemin': caltemin, 'caltemax' : caltemax, 'extfit' : extfit, 'extpar' : np.array(extpar)}

def elemerr(soln,te,sn,fe) :
    '''
    Function to evaluate function for empirical uncertainties
    '''
    out=soln[0]+soln[1]*te+soln[2]*sn
    if len(soln) > 3: out+= soln[3]*fe
    return np.exp(out)

if __name__ == '__main__' :
    main()



def lbd2xyz(l,b,d,R0=8.5) :
    ''' Angular coordinates + distance -> galactocentry x,y,z '''

    brad = b*np.pi/180.
    lrad = l*np.pi/180.

    x = d*np.sin(0.5*np.pi-brad)*np.cos(lrad)-R0
    y = d*np.sin(0.5*np.pi-brad)*np.sin(lrad)
    z = d*np.cos(0.5*np.pi-brad)
    r = np.sqrt(x**2+y**2)
    return x, y, z, r


def zerocal(tab,j,elems,elemtoh,doels,calvers='dr17',extfit=None,calib=False,dwarfs=False) :
    """ Get zeropoint calibration only from abundances of specified indices
    """
    rec = np.zeros(len(doels),dtype=[
                       ('elem','S5'),
                       ('elemfit','i4'),
                       ('mhmin','f4'),
                       ('te0','f4'),
                       ('temin','f4'),
                       ('temax','f4'),
                       ('femin','f4'),
                       ('femax','f4'),
                       ('caltemin','f4'),
                       ('caltemax','f4'),
                       ('gdtemin','f4'),
                       ('gdtemax','f4'),
                       ('extfit','i4'),
                       ('extpar','3f4'),
                       ('exterr','f4'),
                       ('par','3f4'),
                       ('errpar','4f4'),
                       ])
    data = tab[j]
    for iel,el in enumerate(doels) :
        if calvers == 'dr17' :
            pars=dr17cal(el,dwarfs=dwarfs)
        elif calvers == 'dr16' :
            pars=dr16cal(el,dwarfs=dwarfs)
        else :
            pdb.set_trace()           

        pars['elem'] = el
        for key in ['elem','elemfit','mhmin','te0','temin','temax','caltemin','caltemax','gdtemin','gdtemax','extfit','extpar'] :
            rec[iel][key] = pars[key]

        if pars['extfit'] > 0 :
            abun,ok = getabun(data,elems,elemtoh,el,xh=False,terange=[-1,10000],calib=calib,line=0) 
            zero = np.median(abun[ok])
            zero_err = np.median(np.abs(abun[ok]-zero))
            rec[iel]['extpar'][0] = zero
            rec[iel]['exterr'] = zero_err

        if extfit is not None : rec[iel]['extfit'] = extfit
        if rec[iel]['extfit'] == 2 : 
            if el == 'M' :
                rec[iel]['exterr'] = np.sqrt(data['PARAM_COV'][0,3,3])
            elif el == 'alpha' :
                rec[iel]['exterr'] = np.sqrt(data['PARAM_COV'][0,6,6])
            else :
                jel = np.where(elems == el)[0]
                rec[iel]['exterr'] = data['X_H_ERR'][0,jel]

    return rec

def zeroplot(tab3) :
    fig,ax=plots.multi(1,1,figsize=(12,4))
    x=np.arange(len(tab3['SOLAR_ZERO'][0])) 
    plots.plotp(ax,x-0.2,tab3['GIANT_SOLARNEIGH_ZERO'][0],yerr=tab3['GIANT_SOLARNEIGH_ZERO_ERR'][0],color='r',size=10,label='giants')
    plots.plotp(ax,x,tab3['DWARF_SOLARNEIGH_ZERO'][0],yerr=tab3['DWARF_SOLARNEIGH_ZERO_ERR'][0],color='b',size=10,label='dwarfs')
    plots.plotp(ax,x+0.05,tab3['DWARF2_SOLARNEIGH_ZERO'][0],yerr=tab3['DWARF2_SOLARNEIGH_ZERO_ERR'][0],color='c',size=10,label='dwarfs45-50')
    plots.plotp(ax,x+0.2,tab3['SOLAR_ZERO'][0],yerr=tab3['SOLAR_ZERO_ERR'][0],color='g',size=10,label='Vesta')
    ax.set_ylim(-0.4,0.4)
    ax.legend()
    ax.grid()
    ax.set_xticks(x)
    ax.set_xticklabels(tab3['ELEM_SYMBOL'][0].astype(str))

def cal(a,tab3,caldir='cal/') :
    """ Calibrate abundances 
    """

    aspcapmask=bitmask.AspcapBitMask()
    parammask=bitmask.ParamBitMask()

    # calibrate ALL stars with >=0 
    gd=np.where( ((a['ASPCAPFLAG']&aspcapmask.badval()) >= 0) )[0]
    gd=np.where( ((a['ASPCAPFLAG']&aspcapmask.getval('NO_ASPCAP_RESULT')) == 0) )[0]

    giant = np.where( (a['FPARAM'][gd,1] < 2./1300.*(a['FPARAM'][gd,0]-3600)+2.) &
                      (a['FPARAM'][gd,1] < 4.0) & ((a['FPARAM'][gd,0] < 5500) | (a['FPARAM'][gd,1]<2.8)) )[0]
    tmp = np.zeros(len(gd),dtype=bool)
    tmp[giant] = True
    dwarf = np.where(~tmp)[0]

    # initialize calibrated arrays and flag
    # turn on CALRANGE_BAD to start
    # turn off TEFF_CUT to start
    for i in [3,4,5,6] :
        a['PARAM'][:,i] = np.nan
        a['PARAMFLAG'][gd,i] |= parammask.getval('CALRANGE_BAD')
        a['PARAMFLAG'][gd,i] &= ~parammask.getval('TEFF_CUT')

    a['X_H'][:,:] = np.nan
    a['X_M'][:,:] = np.nan
    a['X_H_SPEC'] = a['X_H']
    a['X_M_SPEC'] = a['X_M']
    for iel,el in enumerate(aspcap.elems()[0]) :
        a['ELEMFLAG'][gd,iel] |= parammask.getval('CALRANGE_BAD')
        a['ELEMFLAG'][gd,iel] &= ~parammask.getval('TEFF_CUT')

    # populate [N/M] and [C/M] parameters with uncalibrated values
    for i in [4,5] : 
        a['PARAM'][gd,i] = a['FPARAM'][gd,i]
        a['PARAMFLAG'][gd,i] &= ~parammask.getval('CALRANGE_BAD')

    # calibrate [M/H], [alpha/M], and individual elemets
    els = ['M','alpha']
    els.extend(aspcap.elems()[0])
    elemtoh = aspcap.elems()[1]

    if caldir == 'none' :
        a['PARAM'][gd,3] = a['FPARAM'][gd,3]
        a['PARAM'][gd,6] = a['FPARAM'][gd,6]
        a['PARAMFLAG'][gd,3] &= ~parammask.getval('CALRANGE_BAD')
        a['PARAMFLAG'][gd,6] &= ~parammask.getval('CALRANGE_BAD')
        for iel,el in enumerate(aspcap.elems()[0]) :
            if elemtoh[iel] :
                a['X_M'][gd,iel] = a['FELEM'][gd,iel]-a['FPARAM'][gd,3]
                a['X_H'][gd,iel] = a['FELEM'][gd,iel]
            else :
                a['X_M'][gd,iel] = a['FELEM'][gd,iel]
                a['X_H'][gd,iel] = a['FELEM'][gd,iel]+a['FPARAM'][gd,3]
            a['ELEMFLAG'][gd,iel] &= ~parammask.getval('CALRANGE_BAD')
        a['X_H_SPEC'] = a['X_H']
        a['X_M_SPEC'] = a['X_M']
        return

    for group in ['dwarf','giant'] :
        cal=fits.open(caldir+'/'+group+'_abuncal.fits')[1].data
        if group == 'giant' :
            ok = gd[giant]
        else :
            ok = gd[dwarf]

        for el in els :
            print(el)
            iel = np.where(cal['elem'] == el)[0][0]

            # no calibration if elemfit<0
            if cal['elemfit'][iel] < 0 : continue

            calteffmin=cal['caltemin'][iel]
            calteffmax=cal['caltemax'][iel]
            gdteffmin=cal['gdtemin'][iel]
            gdteffmax=cal['gdtemax'][iel]
            print(el,calteffmin,calteffmax)

            # only calibrate within calteffmin-calteffmax 
            # allow GRIDEDGE here (flag>=0 allows all), but not in named tag
            if el == 'M' :
                gdel=np.where( (a['FPARAM'][ok,0] >= calteffmin)  &
                               (a['FPARAM'][ok,0] <= calteffmax)  &
                               ((a['PARAMFLAG'][ok,3]&parammask.badval()) >= 0) ) [0]
            elif el == 'alpha' :
                gdel=np.where( (a['FPARAM'][ok,0] >= calteffmin)  &
                               (a['FPARAM'][ok,0] <= calteffmax)  &
                               ((a['PARAMFLAG'][ok,6]&parammask.badval()) >= 0) ) [0]
            else :
                jel = np.where(aspcap.elems()[0] == el)[0]
                gdel=np.where( (a['FPARAM'][ok,0] >= calteffmin)  &
                               (a['FPARAM'][ok,0] <= calteffmax)  &
                               ((a['ELEMFLAG'][ok,jel]&parammask.badval()) >= 0) ) [0]

            # flag stars outside of "gd" range, but still populate here
            bdel=np.where((a['FPARAM'][ok,0]<gdteffmin) |
                          (a['FPARAM'][ok,0]>gdteffmax) )[0]

            teff=np.clip(a['FPARAM'][ok[gdel],0],cal['temin'][iel],cal['temax'][iel])
            mh=np.clip(a['FPARAM'][ok[gdel],3],cal['femin'][iel],cal['femax'][iel])
            try: snr=np.clip(a['SNREV'][ok[gdel]],0,200.)
            except:
                print('No SNREV, continue with SNR?')
                pdb.set_trace()
                snr=np.clip(a['SNR'][ok[gdel]],0,200.)

            x = teff-cal['te0'][iel]

            # "internal" calibration
            fit=0
            for iorder in range(0,cal['elemfit'][iel]) :
                fit+=cal['par'][iel,iorder] * x**(iorder+1)

            # "external" calibration
            for iorder in range(2) :
                fit+=cal['extpar'][iel,iorder] * x**(iorder)

            if el == 'M' :
                a['PARAM'][ok[gdel],3] = a['FPARAM'][ok[gdel],3]-fit
                # populate uncertainties with err.apply()
                #a['PARAM_COV'][ok[gdel],3,3] = err.elemerr(cal['errpar'][iel],
                #    a['FPARAM'][ok[gdel],0]-4500,snr-100,a['FPARAM'][ok[gdel],3],quad=True)**2
                a['PARAMFLAG'][ok[gdel],3] &= ~parammask.getval('CALRANGE_BAD')
                a['PARAMFLAG'][ok[bdel],3] |= parammask.getval('TEFF_CUT')
            elif el == 'alpha' :
                a['PARAM'][ok[gdel],6] = a['FPARAM'][ok[gdel],6]-fit
                # populate uncertainties with err.apply()
                #a['PARAM_COV'][ok[gdel],6,6] = err.elemerr(cal['errpar'][iel],
                #    a['FPARAM'][ok[gdel],0]-4500,snr-100,a['FPARAM'][ok[gdel],3],quad=True)**2
                a['PARAMFLAG'][ok[gdel],6] &= ~parammask.getval('CALRANGE_BAD')
                a['PARAMFLAG'][ok[bdel],6] |= parammask.getval('TEFF_CUT')
            else :
                jel = np.where(aspcap.elems()[0] == el)[0]
                print(iel,jel,elemtoh[jel])
                if elemtoh[jel] :
                    a['X_M'][ok[gdel],jel] = a['FELEM'][ok[gdel],jel]-fit-a['FPARAM'][ok[gdel],3]
                    a['X_M_SPEC'][ok[gdel],jel] = a['FELEM'][ok[gdel],jel]-a['FPARAM'][ok[gdel],3]
                    a['X_M_ERR'][ok[gdel],jel] = np.sqrt(a['X_M_ERR'][ok[gdel],jel]**2 + a['PARAM_COV'][ok[gdel],3,3])
                    # [X/H] calculated with all calibrated parameters
                    a['X_H'][ok[gdel],jel] = a['X_M'][ok[gdel],jel]+a['PARAM'][ok[gdel],3]
                    a['X_H_SPEC'][ok[gdel],jel] = a['X_M_SPEC'][ok[gdel],jel]+a['FPARAM'][ok[gdel],3]
                else :
                    a['X_M'][ok[gdel],jel] = a['FELEM'][ok[gdel],jel]-fit
                    a['X_M_SPEC'][ok[gdel],jel] = a['FELEM'][ok[gdel],jel]
                    # [X/H] calculated with all calibrated parameters
                    a['X_H'][ok[gdel],jel] = a['X_M'][ok[gdel],jel]+a['PARAM'][ok[gdel],3]
                    a['X_H_SPEC'][ok[gdel],jel] = a['X_M_SPEC'][ok[gdel],jel]+a['FPARAM'][ok[gdel],3]
                    a['X_H_ERR'][ok[gdel],jel] = np.sqrt(a['X_H_ERR'][ok[gdel],jel]**2 + a['PARAM_COV'][ok[gdel],3,3])
                a['ELEMFLAG'][ok[gdel],jel] &= ~parammask.getval('CALRANGE_BAD')
                a['ELEMFLAG'][ok[bdel],jel] |= parammask.getval('TEFF_CUT')

                # populate uncertainties with err.apply()
                #a['X_H_ERR'][ok[gdel],jel] = err.elemerr(cal['errpar'][jel],
                #    a['FPARAM'][ok[gdel],0]-4500,snr-100,a['FPARAM'][ok[gdel],3],quad=True)
                #a['X_M_ERR'][ok[gdel],jel] = err.elemerr(cal['errpar'][jel],
                #    a['FPARAM'][ok[gdel],0]-4500,snr-100,a['FPARAM'][ok[gdel],3],quad=True)
                # use FERRE uncertainty if larger
                tmp=ok[gdel]
                #j=np.where(a['FELEM_ERR'][tmp,jel] > a['X_H_ERR'][tmp,jel])[0]
                #a['X_H_ERR'][tmp[j],jel] = a['FELEM_ERR'][tmp[j],jel]
                #a['ELEMFLAG'][tmp[j],jel] |= parammask.getval('FERRE_ERR_USED')
                #print(el,'FERRE ERR used: ',len(j))

    return

def dr17cal(el,dwarfs=False) :

    elemfit = 0
    extfit = 4
    te0 = 4500
    temin=0
    temax=10000 
    caltemin=3000
    caltemax=7000
    gdtemin=3000
    gdtemax=7000
    extpar = [0.,0.,0.]
    mhmin=-1

    if el.strip() == 'P' : elemfit=-1
    if el.strip() == 'Ge' : elemfit=-1
    if el.strip() == 'Rb' : elemfit=-1
    if el.strip() == 'Cu' : elemfit=-1
    if el.strip() == 'Yb' : elemfit=-1
    if el.strip() == 'Nd' : elemfit=-1
    if el.strip() == 'C13' : elemfit=-1

    if dwarfs :
        if el.strip() == 'C' : 
            extfit = 0
        elif el.strip() == 'CI' : 
            extfit = 0
        elif el.strip() == 'N' : 
            extfit = 0
            gdtemin=4500
        elif el.strip() == 'Na' : 
            gdtemin=4400
        elif el.strip() == 'Al' : 
            gdtemin=4400
        elif el.strip() == 'S' : 
            gdtemin=4500
        elif el.strip() == 'K' : 
            gdtemin=3300
        elif el.strip() == 'Ca' : 
            gdtemin=3700
        elif el.strip() == 'Ti' : 
            gdtemin=3700
        elif el.strip() == 'TiII' : 
            gdtemin=100000
        elif el.strip() == 'V' : 
            gdtemin=4800
        elif el.strip() == 'Cr' : 
            gdtemin=4400
        elif el.strip() == 'Mn' : 
            gdtemin=4000
        elif el.strip() == 'Fe' : 
            extfit = 0
        elif el.strip() == 'Co' : 
            gdtemin=100000
        elif el.strip() == 'Ni' : 
            gdtemin=0
        elif el.strip() == 'Ce' : 
            gdtemin=100000
            gdtemax=6800
        elif el.strip() == 'Nd' : 
            gdtemin=4500
        elif el.strip() == 'M' : 
            extfit = 0
        gdtemin=np.max([3500,gdtemin])
    else :
        if el.strip() == 'C' : 
            extfit = 0
        elif el.strip() == 'CI' : 
            extfit = 0
            gdtemin=0
        elif el.strip() == 'N' : 
            extfit = 0
        elif el.strip() == 'Na' : 
            gdtemin=3600
        elif el.strip() == 'Al' : 
            gdtemin=3700
        elif el.strip() == 'S' : 
            gdtemin=3700
        elif el.strip() == 'K' : 
            gdtemin=3600
        elif el.strip() == 'Ca' : 
            gdtemin=3300
        elif el.strip() == 'Ti' : 
            gdtemin=3800
        elif el.strip() == 'TiII' : 
            gdtemin=3800
        elif el.strip() == 'V' : 
            gdtemin=3500
        elif el.strip() == 'Cr' : 
            gdtemin=3700
        elif el.strip() == 'Mn' : 
            gdtemin=3800
        elif el.strip() == 'Fe' : 
            extfit = 0
        elif el.strip() == 'Co' : 
            gdtemin=3600
        elif el.strip() == 'Ni' : 
            gdtemin=3500
        elif el.strip() == 'Ce' : 
            gdtemin=3900
        elif el.strip() == 'Nd' : 
            gdtemin=4400
        elif el.strip() == 'M' : 
            extfit = 0
  
    return {'elemfit': elemfit, 'mhmin' : mhmin, 'te0': te0, 'temin': temin, 'temax': temax, 
            'caltemin': caltemin, 'caltemax' : caltemax, 'gdtemin' : gdtemin, 'gdtemax': gdtemax,
            'extfit' : extfit, 'extpar' : np.array(extpar)}
