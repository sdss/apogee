from apogee.utils import apload
from apogee.utils import apselect
from apogee.aspcap import elem
from apogee.aspcap import teffcomp
from apogee.aspcap import loggcomp
from tools import html
from tools import match
from tools import struct
from tools import plots
from tools import fit
try: from tools import vfit
except: pass
import os
import shutil
import pdb
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from astropy.io import fits
from astropy.io import ascii

def allField(files=['apo*/*/apField-*.fits','apo*/*/apFieldC-*.fits'],out='allField.fits',verbose=False) :
    '''
    Concatenate set of apField files
    '''
    # concatenate the structures
    all=struct.concat(files,verbose=verbose)

    # write out the file
    if out is not None:
        print('writing',out)
        struct.wrfits(all,out)

    return all

def allCal(files=['clust???/aspcapField-*.fits','cal???/aspcapField-*.fits'],nelem=15,out='allCal.fits',allfield=None) :
    '''
    Concatenate aspcapField files, adding ELEM tags if not there
    '''
    # concatenate the structures
    all=struct.concat(files,verbose=True)

    # add elements tags if we don't have them
    try :
        test=all['FELEM'][0]
    except :
        n=len(all)
        form='{:<d}f4'.format(nelem)
        iform='{:<d}f4'.format(nelem)
        all=struct.add_cols(all,np.zeros(n,dtype=[('FELEM',form),('FELEM_ERR',form),
                                                  ('ELEM',form),('ELEM_ERR',form),('ELEMFLAG',iform)]))

    # add in NINST information from allField file
    if allfield is not None:
        a=fits.open(allfield)[1].data
        i1,i2=match.match(np.core.defchararray.add(all['APOGEE_ID'],all['LOCATION_ID'].astype(str)),
                    np.core.defchararray.add(a['APOGEE_ID'],a['LOCATION_ID'].astype(str)))
        n=len(all)
        all=struct.add_cols(all,np.zeros(n,dtype=[('NINST','3i4')]))
        all['NINST'][i1]=a['NINST'][i2]

    # write out the file
    if out is not None:
        print('writing',out)
        struct.wrfits(all,out)

    return all

def concat(files,hdu=1) :
    '''
    Create concatenation of all apField files
    '''
    if type(files) == str:
        files=[files]
    allfiles=[]
    for file in files :   
        allfiles.extend(glob.glob(file))
    if len(allfiles) == 0 :
        print('no files found!',file)
        return

    for file in files :
        print(file)
        a=fits.open(file)[hdu].data
        try:
            all=struct.append(all,a)
        except :
            all=a
        print(len(all), len(a))
    return all

def hrsample(indata,hrdata,maxbin=50,raw=True) :
    ''' 
    selects stars covering HR diagram as best as possible from input sample
    '''
    i1,i2 = match.match(indata['APOGEE_ID'],hrdata['APOGEE_ID'])
    gd=[]
    for teff in np.arange(3500,6000,500) :
        gdteff=apselect.select(hrdata[i2],badval=['STAR_BAD'],badtarg=['EMBEDDED','EXTENDED'],teff=[teff,teff+500],sn=[100,1000],raw=raw)
        for logg in np.arange(0,5,1) :
            j=apselect.select(hrdata[i2[gdteff]],logg=[logg,logg+1],raw=True)
            gdlogg=gdteff[j]
            for mh in np.arange(-2.5,0.5,0.5) :
                j=apselect.select(hrdata[i2[gdlogg]],mh=[mh,mh+0.5],raw=True)
                j=gdlogg[j]
                js=np.argsort(hrdata[i2[j]]['SNR'])[::-1]
                x = j[js] if len(j) < maxbin else j[js[0:maxbin]]
                gd.extend(x)
    return i1[gd],i2[gd]

def calsample(indata=None,file='clust.html',plot=True,clusters=True,apokasc='APOKASC_cat_v3.6.0',cal1m=True,galcen=True,dir='cal',hrdata=None,optical=True) :
    '''
    selects a calibration subsample from an input apField structure, including several calibration sub-classes: 
        cluster, APOKASC stars, 1m calibration stars. Creates cluster web pages/plots if requested
    '''

    if indata is None :
        indata=allField(files=['apo25m/*/apField-*.fits','apo1m/calibration/apField-*.fits'],out=None,verbose=True)

    j=np.where(indata['COMMISS'] == 0)[0]
    data=indata[j]

    try: os.mkdir(dir)
    except: pass
    if clusters :
        jc=[]
        clusts=apselect.clustdata()
        f=html.head(file=dir+'/'+file)
        f.write('<TABLE BORDER=2>\n')
        clust=apselect.clustdata()
        for ic in range(len(clust.name)) :
            j=apselect.clustmember(data,clust[ic].name,plot=plot,hard=dir)
            print(clust[ic].name,len(j))
            jc.extend(j)
            f.write('<TR><TD>'+clust[ic].name+'<TD>{:12.6f}<TD>{:12.6f}<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                    clust[ic].ra,clust[ic].dec,clust[ic].rad,clust[ic].rv,clust[ic].drv))
            f.write('<TD><A HREF='+clust[ic].name+'_pos.jpg><IMG SRC='+clust[ic].name+'_pos.jpg width=300></A>\n')
            f.write('<TD><A HREF='+clust[ic].name+'_pos.jpg><IMG SRC='+clust[ic].name+'_rv.jpg width=300></A>\n')
            f.write('<TD><A HREF='+clust[ic].name+'_pos.jpg><IMG SRC='+clust[ic].name+'_cmd.jpg width=300></A>\n')
        html.tail(f)
        print('Number of cluster stars: ',len(jc))
        mklinks(data,jc,dir+'/clust')
    
    if apokasc is not None :
        jc=[]
        apokasc = fits.open(os.environ['APOGEE_DIR']+'/data/apokasc/'+apokasc+'.fits')[1].data
        i1,i2=match.match(data['APOGEE_ID'],apokasc['2MASS_ID'])
        rgb=np.where(apokasc['CONS_EVSTATES'][i2] == 'RGB')[0]
        print('Number of APOKASC RGB stars (every 3rd): ',len(rgb[0:-1:3]))
        jc.extend(i1[rgb][0:-1:3])
        rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC')[0]
        print('Number of APOKASC RC stars (every 2nd): ',len(rc[0:-2:2]))
        jc.extend(i1[rc][0:-1:2])
        lowg=np.where((apokasc['LOGG_SYD_SCALING'][i2] < 2) & (apokasc['LOGG_SYD_SCALING'][i2] > 0.1))[0]
        print('Number of APOKASC low log g  stars: ',len(lowg))
        jc.extend(i1[lowg])
        highg=np.where((apokasc['LOGG_SYD_SCALING'][i2] > 3.8) & (apokasc['LOGG_SYD_SCALING'][i2] < 5.5))[0]
        print('Number of APOKASC high log g  stars: ',len(highg))
        jc.extend(i1[highg])
        lowz=np.where((apokasc['FE_H_ADOP_COR'][i2] < -1.) & (apokasc['FE_H_ADOP_COR'][i2] > -90.))[0]
        print('Number of APOKASC low [Fe/H] stars: ',len(lowz))
        jc.extend(i1[lowz])
        mklinks(data,jc,dir+'/apokasc')
    
    if galcen :
        j=np.where(data['FIELD'] == 'GALCEN')[0]
        print('Number of GALCEN stars: ',len(j))
        mklinks(data,j,dir+'/galcen')

    if cal1m :
        j=np.where(data['FIELD'] == 'calibration')[0]
        print('Number of 1m calibration stars: ',len(j))
        mklinks(data,j,dir+'/cal')

    if optical :
        stars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/validation_stars_DR16.txt',names=['id'],format='fixed_width_no_header')
        i1,i2=match.match(data['APOGEE_ID'],stars['id'])
        print('Number of optical validation stars: ',len(i1))
        mklinks(data,i1,dir+'/optical')

    if hrdata is not None:
        i1, i2 = hrsample(data,hrdata)
        print('Number of HR sample stars: ',len(i1))
        mklinks(data,i1,dir+'/hr')

    return indata

def mklinks(data,j,out,n=48) :
    """ Create links in n different output directories for requested indices
    """

    # remove existing output directories, create new ones
    cleandir(out,n)
    # create symbolic links in output directories
    nsplit=len(j)//n+1
    for i in range(len(j)) :
        symlink(data[j[i]],out,i//nsplit)


def cleandir(out,n) :
    '''
    auxiliary routine to clean and remake calibration directories
    '''
    for i in range(n) : 
        try:
            shutil.rmtree('{:s}{:03d}'.format(out,i))
        except : pass
        try:
            os.mkdir('{:s}{:03d}'.format(out,i))
        except : pass

def symlink(data,out,idir) :
    '''
    auxiliary routine to create symlinks to appropriate files from calibration directories
    '''
    outfile='{:s}{:03d}/{:s}.{:s}.fits'.format(
            out,idir,os.path.splitext(os.path.basename(data['FILE']))[0],data['FIELD'])
    if data['TELESCOPE'] == 'apo25m' or data['TELESCOPE'] == 'lco25m' :
        infile='../../{:s}/{:s}/{:s}'.format(data['TELESCOPE'],data['FIELD'],data['FILE'])
    else :
        infile='../../{:s}/calibration/{:s}'.format(data['TELESCOPE'],data['FILE'])
    os.symlink(infile,outfile)

def docal(vers,clobber=False,allstar=True,hr=True,teff=True,logg=True,vmicro=True,vmacro=True,elemcal=True,out=None,stp=False,cal='dr14',calib=False) :
    '''
    Derives all calibration relations and creates plots of them, as requested
    '''

    # output subdirectory
    try:
        os.mkdir('cal')
    except:
        pass
    print(os.getcwd())

    # combine aspcapField files into allCal
    if clobber :
        allc=allCal(['hr???/aspcapField*.fits','cal???/aspcapField*.fits','clust???/aspcapField*.fits'],out='allCal-'+vers+'.fits')
    else :
        try:
            allc=fits.open('allCal-'+vers+'.fits')[1].data
        except:
            allc=allCal(['hr???/aspcapField*.fits','cal???/aspcapField*.fits','clust???/aspcapField*.fits'],out='allCal-'+vers+'.fits')
    if allstar :
        c=fits.open('allStar-'+vers+'.fits')[1].data
        cc=fits.open('allStar-'+vers+'.fits')[3].data
    else :
        c=allc
        cc=apl.allStar()[3].data
    print('Total stars:',len(c))

    figs=[]
    ytitle=[]
    # HR diagram
    if hr :
        reload(apselect)
        fig,ax=plots.multi(1,1)
        if calib : param='PARAM'
        else : param='FPARAM'
        plots.plotc(ax,c[param][:,0],c[param][:,1],c[param][:,3],xr=[6000,3000],yr=[5,-1],zr=[-2,0.5])
        plt.savefig(out+'hr.jpg')                                                                                                 
    figs.append(['hr.jpg','hr.jpg'])
    ytitle.append('HR')

    allcal={}
    # Teff vs photometric
    if teff :
        allcal['giant_teffcal'] = teffcomp.ghb(c,ebvmax=0.02,glatmin=10,out=out+'giant_teffcomp',yr=[-750,750],dwarf=False,calib=calib)
        if out is not None :
            struct.wrfits(struct.dict2struct(allcal['giant_teffcal']),out+'giant_tecal.fits')
        allcal['dwarf_teffcal'] = teffcomp.ghb(c,ebvmax=0.02,glatmin=10,trange=[4000,7500],out=out+'dwarf_teffcomp',yr=[-750,750],dwarf=True,calib=calib)
        if out is not None :
            struct.wrfits(struct.dict2struct(allcal['dwarf_teffcal']),out+'dwarf_tecal.fits')
        if stp : pdb.set_trace()
    figs.append(['giant_teffcomp.jpg','dwarf_teffcomp.jpg'])
    ytitle.append('Teff')
    figs.append(['giant_teffcomp_b.jpg','dwarf_teffcomp_b.jpg'])
    ytitle.append('Teff')

    # log g vs asteroseismic
    if logg :
        allcal['rgbrcsep' ] = loggcomp.rcrgb(c,out=out+'rcrgbsep')
        allcal['loggcal'] = loggcomp.apokasc(c,plotcal=False,out=out+'loggcomp',calib=calib)
        if out is not None :
            struct.wrfits(struct.dict2struct(dict(allcal['rgbrcsep'].items()+allcal['loggcal'].items())),out+'loggcal.fits')
        if stp : pdb.set_trace()
    figs.append(['rcrgbsep.jpg','none.jpg'])
    ytitle.append('log g')
    figs.append(['loggcomp_b.jpg','loggcomp.jpg'])
    ytitle.append('log g')

    # vmicro calibration
    if vmicro :
        print("vmicro fit, cubic in log g, linear in [M/H]")
        print("sample limited to FERRE vmicro error <0.01 ")
        vfit.fit_vmicro(c,degree=3,reject=0.15,mhrange=[-2,1],loggrange=[-0.3,4.9],vmrange=[0,7],teffrange=[3550,6500],vrange=[0.55,4],maxerr=0.01,func=vfit.vm3_1,out=out+'vmicro3_1')

        print("full sample (maxerr=0.1)")
        vfit.fit_vmicro(c,degree=3,reject=0.15,mhrange=[-2,1],loggrange=[-0.3,4.9],vmrange=[0,7],teffrange=[3500,6500],vrange=[0.55,4],maxerr=0.1,func=vfit.vm3_1,out=out+'vmicro3_1all')

        #dwarfs only
        print("dwarfs only, fit as f(Teff)")
        dw=np.where(c['FPARAM'][:,1] > 4)[0]
        vfit.fit_vmicro(c,reject=0.15,mhrange=[-2,1],loggrange=[4,5],vmrange=[0,7],teffrange=[3500,8000],vrange=[0.55,4],maxerr=0.1,func=vfit.vm1t,out=out+'vmicro1t')
        fig,ax=plots.multi(1,1)
        plots.plotc(ax,c['FPARAM'][dw,0],10**c['FPARAM'][dw,2],c['FPARAM'][dw,3],xr=[3500,8000],xt='Teff',
           yr=[0,4],yt='vmicro',zr=[-2,0.5],zt='[M/H]',colorbar=True)

    # vmacro
    if vmacro :
        vfit.fit_vmacro(c,mhrange=[-2.5,1],reject=0.3,maxerr=0.1,out=out+'vmacro_2d')

    # elemental abundances
    if elemcal :
        elems=np.append(cc['ELEM_SYMBOL'][0],['M','alpha'])
        # use allCal file for uncertainty calibration, so we have multiple visits
        # use allstar for calibration, so we have full solar circle sample
        errcal=elem.cal(allc,cc['ELEM_SYMBOL'][0],cc['ELEMTOH'][0],elems,hard=out+'giants_',cal=cal,errpar=True,plot=False,calib=calib)
        allcal['giantcal']=elem.cal(c,cc['ELEM_SYMBOL'][0],cc['ELEMTOH'][0],elems,hard=out+'giants_',cal=cal,errpar=True,calib=calib)
        allcal['giantcal']['errpar']=errcal['errpar']
        if out is not None :
            struct.wrfits(allcal['giantcal'],out+'giantcal.fits')
        errcal=elem.cal(allc,cc['ELEM_SYMBOL'][0],cc['ELEMTOH'][0],elems,hard=out+'dwarfs_',dwarfs=True,errpar=True,plot=False,cal=cal,calib=calib)
        allcal['dwarfcal']=elem.cal(c,cc['ELEM_SYMBOL'][0],cc['ELEMTOH'][0],elems,hard=out+'dwarfs_',dwarfs=True,cal=cal,calib=calib)
        allcal['dwarfcal']['errpar']=errcal['errpar']
        if out is not None :
            struct.wrfits(allcal['dwarfcal'],out+'dwarfcal.fits')
        if stp : pdb.set_trace()
    figs.append(['giants_all.jpg','dwarfs_all.jpg'])
    ytitle.append('clusters')
    figs.append(['giants_allsolar.jpg','dwarfs_allsolar.jpg'])
    ytitle.append('solar circule')
    figs.append(['giants_M.jpg','dwarfs_M.jpg'])
    ytitle.append('cluster [M/H]')

    html.htmltab(figs,xtitle=['giants','dwarfs'],ytitle=ytitle,file=out+vers+'.html')
    return allcal

def comp(plots=['hr','giant_teffcomp','dwarf_teffcomp','rcrgbsep','loggcomp_b','loggcomp','giants_all','clust_key','dwarfs_all','giants_allsolar','dwarfs_allsolar','giants_M','dwarfs_M','giants_err_all','dwarfs_err_all'],runs=['l31a','l31b','l30b_vm4','l31b_vm4','l31a_asset'],out=None) :
    '''
    Generate web page with (existing) calibration plots for multiple runs
   
    Keyword args:
         plots=[list of plot names]
         runs=[list of runs]    
         out=name of output HTML file
    '''
    grid = []
    for plot in plots :
        y=[]
        for run in runs :
            y.append(run+'/'+run+out+plot+'.jpg')
        grid.append(y)
            
    html.htmltab(grid,file=out,ytitle=plots,xtitle=runs)


def compstars(d1,d2,out=None) :
    '''
    Creates plots to compare 2 different version
    '''
    v1=fits.open(d1+'/'+d1+'/allCal-'+d1+'.fits')[1].data
    v2=fits.open(d2+'/'+d2+'/allCal-'+d2+'.fits')[1].data
    i1,i2=match.match(v1['APOGEE_ID'],v2['APOGEE_ID'])

    fig,ax=plots.multi(1,7,hspace=0.001,figsize=(8,20))
    plots.plotc(ax[0],v1['FPARAM'][i1,0],v2['FPARAM'][i2,0]-v1['FPARAM'][i1,0],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta Teff$',yr=[-1000,1000],xt='Teff')
    plots.plotc(ax[1],v1['FPARAM'][i1,0],v2['FPARAM'][i2,1]-v1['FPARAM'][i1,1],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta log g$',yr=[-1,1],xt='Teff')
    plots.plotc(ax[2],v1['FPARAM'][i1,0],10.**v2['FPARAM'][i2,2]-10.**v1['FPARAM'][i1,2],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta vmicro$',yr=[-1,1],xt='Teff')
    plots.plotc(ax[3],v1['FPARAM'][i1,0],v2['FPARAM'][i2,3]-v1['FPARAM'][i1,3],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta [M/H]$',yr=[-0.75,0.75],xt='Teff')
    plots.plotc(ax[4],v1['FPARAM'][i1,0],v2['FPARAM'][i2,4]-v1['FPARAM'][i1,4],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta [C/M]$',yr=[-0.75,0.75],xt='Teff')
    plots.plotc(ax[5],v1['FPARAM'][i1,0],v2['FPARAM'][i2,5]-v1['FPARAM'][i1,5],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta [N/M]$',yr=[-0.75,0.75],xt='Teff')
    plots.plotc(ax[6],v1['FPARAM'][i1,0],v2['FPARAM'][i2,6]-v1['FPARAM'][i1,6],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta \alpha/M]$',yr=[-0.75,0.75],xt='Teff')
    if out is not None:
        plt.savefig(out+'.jpg')

    # plots as a function of delta logvmicro
    fig,ax=plots.multi(1,7,hspace=0.001,figsize=(8,20))
    plots.plotc(ax[0],v1['FPARAM'][i1,2]-v2['FPARAM'][i2,2],v2['FPARAM'][i2,0]-v1['FPARAM'][i1,0],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta Teff$',yr=[-1000,1000],xt=r'$\Delta log vmicro$' )
    plots.plotc(ax[1],v1['FPARAM'][i1,2]-v2['FPARAM'][i2,2],v2['FPARAM'][i2,1]-v1['FPARAM'][i1,1],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta log g$',yr=[-1,1],xt=r'$\Delta log vmicro$' )
    plots.plotc(ax[2],v1['FPARAM'][i1,2]-v2['FPARAM'][i2,2],10.**v2['FPARAM'][i2,2]-10.**v1['FPARAM'][i1,2],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta vmicro$',yr=[-1,1],xt=r'$\Delta log vmicro$' )
    plots.plotc(ax[3],v1['FPARAM'][i1,2]-v2['FPARAM'][i2,2],v2['FPARAM'][i2,3]-v1['FPARAM'][i1,3],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta [M/H]$',yr=[-0.75,0.75],xt=r'$\Delta log vmicro$' )
    plots.plotc(ax[4],v1['FPARAM'][i1,2]-v2['FPARAM'][i2,2],v2['FPARAM'][i2,4]-v1['FPARAM'][i1,4],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta [C/M]$',yr=[-0.75,0.75],xt=r'$\Delta log vmicro$' )
    plots.plotc(ax[5],v1['FPARAM'][i1,2]-v2['FPARAM'][i2,2],v2['FPARAM'][i2,5]-v1['FPARAM'][i1,5],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta [N/M]$',yr=[-0.75,0.75],xt=r'$\Delta log vmicro$' )
    plots.plotc(ax[6],v1['FPARAM'][i1,2]-v2['FPARAM'][i2,2],v2['FPARAM'][i2,6]-v1['FPARAM'][i1,6],v1['FPARAM'][i1,3],zr=[-2,0.5],colorbar=True,zt='[M/H]',yt=r'$\Delta \alpha/M]$',yr=[-0.75,0.75],xt=r'$\Delta log vmicro$' )
    if out is not None:
        plt.savefig(out+'_dvmicro.jpg')

def starcomp(ref='l31a',comps=['l31b','l31b_vm4','l30b_vm4','l31a_asset'],out=None) :
    '''
    Creates web page for comparison of mulitple version
    '''
    grid=[]
    for comp in comps :
        compstars(ref,comp,out='comp/'+ref+'_'+comp)
        y=[ref+'_'+comp+'.jpg',ref+'_'+comp+'_dvmicro.jpg']
        grid.append(y)
    html.htmltab(np.asarray(grid).T.tolist(),file=out,xtitle=comps)

def errcomp(vers=['dr14','dr13','dr12'],els=['alpha','O','Mg','Ni','M'],out='comp.html') :
    '''
    Creates web page for comparison of mulitple version
    '''
    grid=[]
    for ver in vers :
        y=[]
        ytit=[]
        for el in els :
          for plot in ['err','err_sn','clusterr_all'] :
             y.append('../'+ver+'/cal/'+el+'_'+plot+'.jpg')
             ytit.append(el+'_'+plot)
        grid.append(y)
    html.htmltab(np.asarray(grid).T.tolist(),file=out,xtitle=vers,ytitle=ytit)

def errplots(tags=['ALPHA_M','O_FE','MG_FE','NI_FE','M_H'],cannon=None) :
    a=apl.allStar()[1].data
    a3=apl.allStar()[3].data
    if cannon is not None :
        cannon=fits.open(cannon)[1].data
        abun=cannon
    else :
        abun = a
    solar=apselect.select(a,badval='STAR_BAD',logg=[-1,3.8],glon=[70,110],glat=[-5,5],alpha=[-0.15,0.15])
    for tag in tags :
        if tag == 'ALPHA_M' : el = 'alpha'
        elif tag == 'PARAM_ALPHA_M' : el = 'alpha'
        elif tag == 'PARAM_M_H' : el = 'M'
        else : el = tag.split('_')[0].capitalize()
        print(tag,el)
        # scatter in solar circle stars
        errfit(a['TEFF'][solar],a['SNR'][solar],a['PARAM'][solar,3],abun[tag][solar],out='cal/'+el,title=el,snbins=np.arange(50,250,25),tebins=np.arange(3600,5200,400),mhbins=np.arange(-1,0.75,0.25))
        if cannon is not None:
            a[tag] = cannon[tag]
            if tag == 'M_H' :
                a['PARAM'][:,3] = cannon[tag]
            elif tag == 'ALPHA_M' :
                a['PARAM'][:,6] = cannon[tag]
            else :
                j=np.where(a3['ELEM_SYMBOL'][0] == el)[0]
                a['X_M'][:,j[0]] = cannon[tag]
        # scatter in clusters
        elem.cal(a,a3['ELEM_SYMBOL'][0],a3['ELEMTOH'][0],[el,el],hard='cal/'+el+'_clust',errpar=True,calib=True)

def allplots() :
    global apl

    os.chdir( '../dr14')
    apl=apload.apLoad(dr='dr14')
    errplots()
    os.chdir('../dr13')
    apl=apload.apLoad(dr='dr13')
    errplots()
    os.chdir('../dr12')
    apl=apload.apLoad(dr='dr12')
    errplots(tags=['PARAM_ALPHA_M','O_H','MG_H','NI_H','PARAM_M_H'])

