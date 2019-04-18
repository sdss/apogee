from apogee.utils import apload
from apogee.utils import apselect
from apogee.aspcap import elem
from apogee.aspcap import teffcomp
from apogee.aspcap import loggcomp
from apogee.aspcap import aspcap
from apogee.aspcap import qa
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

def allField(files=['apo*/*/a?Field-*.fits','apo*/*/a?FieldC-*.fits','lco*/*/a?Field-*.fits'],out='allField.fits',verbose=False) :
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

def allCal(files=['clust???/aspcapField-*.fits','cal???/aspcapField-*.fits'],nelem=15,out='allCal.fits',allfield=None,elemcal=True) :
    '''
    Concatenate aspcapField files, adding ELEM tags if not there
    '''
    # concatenate the structures
    all=struct.concat(files,verbose=True,fixfield=True)

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
        #struct.wrfits(all,out)
        hdulist=fits.HDUList()
        hdu=fits.BinTableHDU.from_columns(all)
        hdulist.append(hdu)
        filelist=glob.glob(files[0])
        hdu=fits.open(filelist[0])[3]
        hdulist.append(hdu)
        hdu=fits.open(filelist[0])[3]
        hdulist.append(hdu)
        hdulist.writeto(out,overwrite=True)

def summary(hdulist,elemcal=False,out='allCal.fits') :
    """ Create QA summary page and plots
    """
    all=hdulist[1].data
    try: os.mkdir('plots/')
    except: pass

    # HR diagrams
    aspcap.hr(all,hard='plots/hr.png',xr=[8000,3000],grid=True)
    aspcap.hr(all,hard='plots/hrhot.png',xr=[20000,3000],iso=True)
    aspcap.multihr(all,hard='plots/multihr.png')
    #teffcomp.ghb(all,ebvmax=0.02,glatmin=10,out='plots/giant_teffcomp',yr=[-750,750],dwarf=False,calib=False)
    #loggcomp.apokasc(all,plotcal=False,out='plots/loggcomp',calib=False)
    grid=[['plots/hr.png','plots/multihr.png','plots/hrhot.png']] 

    # Master summary HTML file
    f=html.head(file=out.replace('.fits','.html'))
    f.write(html.table(grid))
    f.write('<br>Uncalibrated parameters:<br>')
    ids = ['VESTA','alpha_Boo']
    j=[]
    for id in ids: j.extend( np.where( (np.core.defchararray.strip(all['APOGEE_ID']) == id) & (all['VISIT'] == 0)) [0] )
    f.write(html.table(all['FPARAM'][j],plots=False,ytitle=ids,xtitle=aspcap.params()[1]))
    # table of abundances (relative to M)
    f.write('<br>Uncalibrated abundances:<br>')
    abun=all['FELEM'][j,0,:]
    xtit=[]
    for i in range(len(hdulist[3].data['ELEM_SYMBOL'][0])) : 
        if hdulist[3].data['ELEMTOH'][0][i] == 1 : abun[:,i]-=all['FPARAM'][j,3]
        xtit.append('['+hdulist[3].data['ELEM_SYMBOL'][0][i]+'/M]')
    f.write(html.table(abun,plots=False,ytitle=ids,xtitle=xtit))
    f.write('<p> <a href=calib/'+out.replace('.fits','.html')+'> Calibration plots</a>')
    f.write('<p> <a href=qa/elem_chem.html> Chemistry plots</a>')
    f.write('<p> <a href=qa/repeat.html> Duplicate observations plots, including APO/LCO</a>')
    f.write('<p> <a href=qa/dr14_diffs.html> DR14 comparison plots</a>')
    html.tail(f)

    # do the calibration and calibration and QA plots
    try: os.mkdir('calib/')
    except: pass
    docal(out,clobber=False,allstar=False,hr=False,teff=True,logg=True,vmicro=False,vmacro=False,elemcal=elemcal,out='calib/',stp=False,cal='default',calib=False) 
    pdb.set_trace()
    try: os.mkdir('calibrated/')
    except: pass
    try: docal(out,clobber=False,allstar=False,hr=False,teff=True,logg=True,vmicro=False,vmacro=False,elemcal=elemcal,out='calibrated/',stp=False,cal='default',calib=True) 
    except: pass
    try: os.mkdir('qa/')
    except: pass
    hdulist=fits.open(out)
    qa.dr14comp(hdulist,out='qa/',elem=elemcal)
    if elemcal : 
        qa.plotelems(hdulist,out='qa/')
        qa.plotcn(hdulist,out='qa/')
    qa.repeat(hdulist,out='repeat/giant_',elem=elemcal,logg=[-1,3.8])
    qa.repeat(hdulist,out='repeat/dwarf_',elem=elemcal,logg=[3.8,5.5])

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
    for teff in np.arange(3000,6000,500) :
        gdteff=apselect.select(hrdata[i2],badval=['STAR_BAD'],badtarg=['EMBEDDED','EXTENDED'],teff=[teff,teff+500],sn=[100,1000],raw=raw)
        for logg in np.arange(-0.5,5.5,1) :
            j=apselect.select(hrdata[i2[gdteff]],logg=[logg,logg+1],raw=True)
            gdlogg=gdteff[j]
            for mh in np.arange(-2.5,0.5,0.5) :
                j=apselect.select(hrdata[i2[gdlogg]],mh=[mh,mh+0.5],raw=True)
                j=gdlogg[j]
                js=np.argsort(hrdata[i2[j]]['SNR'])[::-1]
                x = j[js] if len(j) < maxbin else j[js[0:maxbin]]
                gd.extend(x)
    return i1[gd],i2[gd]

def calsample(indata=None,file='clust.html',plot=True,clusters=True,apokasc='APOKASC_cat_v4.4.2',
              cal1m=True,coolstars=True,dir='cal',hrdata=None,optical='cal_stars_20190329.txt',ns=True,special=None,Ce=True,ebvmax=None,snmin=75,mkall=True) :
    '''
    selects a calibration subsample from an input apField structure, including several calibration sub-classes: 
        cluster, APOKASC stars, 1m calibration stars. Creates cluster web pages/plots if requested
    '''

    if indata is None :
        indata=allField(files=['apo25m/*/a?Field-*.fits','lco25m/*/a?Field-*.fits','apo1m/calibration/a?Field-*.fits'],out=None,verbose=True)

    j=np.where((indata['COMMISS'] == 0) & (indata['SNR'] > 75) )[0]
    data=indata[j]

    try: os.mkdir(dir)
    except: pass
    all=[]
    if clusters :
        jc=[]
        clusts=apselect.clustdata()
        fstars=open(dir+'/allclust.txt','w')
        f=html.head(file=dir+'/'+file)
        f.write('<A HREF=allclust.txt> cluster stars list </a>')
        f.write('<TABLE BORDER=2>\n')
        f.write('<TR><TD>NAME<TD>RA<TD>DEC<TD>Radius<TD>RV<TD>Delta RV<TD>Position criterion<TD>RV criterion<TD>PM criterion<TD>Parallax criterion<TD> CMD')
        clust=apselect.clustdata()
        for ic in range(len(clust.name)) :
            print(clust[ic].name)
            j=apselect.clustmember(data,clust[ic].name,plot=plot,hard=dir)
            print(clust[ic].name,len(j))
            # clusters to exclude here
            if (clust[ic].name not in ['OmegaCen','Pal1','Pal6','Pal5','Terzan12'])  and (len(j) >= 5): jc.extend(j)
            f.write('<TR><TD><A HREF='+clust[ic].name+'.txt>'+clust[ic].name+'</A><TD>{:12.6f}<TD>{:12.6f}<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                    clust[ic].ra,clust[ic].dec,clust[ic].rad,clust[ic].rv,clust[ic].drv))
            f.write('<TD><A HREF='+clust[ic].name+'_pos.jpg><IMG SRC='+clust[ic].name+'_pos.jpg width=300></A>\n')
            f.write('<TD><A HREF='+clust[ic].name+'_rv.jpg><IMG SRC='+clust[ic].name+'_rv.jpg width=300></A>\n')
            f.write('<TD><A HREF='+clust[ic].name+'_pm.jpg><IMG SRC='+clust[ic].name+'_pm.jpg width=300></A>\n')
            f.write('<TD><A HREF='+clust[ic].name+'_parallax.jpg><IMG SRC='+clust[ic].name+'_parallax.jpg width=300></A>\n')
            f.write('<TD><A HREF='+clust[ic].name+'_cmd.jpg><IMG SRC='+clust[ic].name+'_cmd.jpg width=300></A>\n')
            np.savetxt(dir+'/'+clust[ic].name+'.txt',data[jc]['APOGEE_ID'],fmt='%s')
            for star in data[jc]['APOGEE_ID'] : fstars.write('{:s} {:s}\n'.format(star,clust[ic].name))
        html.tail(f)
        fstars.close()
        print('Number of cluster stars: ',len(jc))
        mklinks(data,jc,dir+'_clust')
        all.extend(jc)
    
    if apokasc is not None :
        jc=[]
        apokasc = fits.open(os.environ['APOGEE_DIR']+'/data/apokasc/'+apokasc+'.fits')[1].data
        i1,i2=match.match(data['APOGEE_ID'],apokasc['2MASS_ID'])
        rgb=np.where(apokasc['CONS_EVSTATES'][i2] == 'RGB')[0]
        print('Number of APOKASC RGB stars (every 4th): ',len(rgb[0:-1:3]))
        jc.extend(i1[rgb][0:-1:4])
        rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC')[0]
        print('Number of APOKASC RC stars (every 2nd): ',len(rc[0:-2:2]))
        jc.extend(i1[rc][0:-1:2])
        rc=np.where(apokasc['CONS_EVSTATES'][i2] == '2CL')[0]
        print('Number of APOKASC 2CL stars: ',len(rc))
        jc.extend(i1[rc])
        rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC/2CL')[0]
        print('Number of APOKASC RC/2CL stars: ',len(rc))
        jc.extend(i1[rc])
        logg = 'APOKASC2_LOGG'  # older version used LOGG_SYSD_SCALING
        lowg=np.where((apokasc[logg][i2] < 2) & (apokasc[logg][i2] > 0.1))[0]
        print('Number of APOKASC low log g  stars: ',len(lowg))
        jc.extend(i1[lowg])
        highg=np.where((apokasc[logg][i2] > 3.8) & (apokasc[logg][i2] < 5.5))[0]
        print('Number of APOKASC high log g  stars: ',len(highg))
        jc.extend(i1[highg])
        highg=np.where((apokasc['LOGG_DW'][i2] > 3.8) & (apokasc['LOGG_DW'][i2] < 5.5))[0]
        print('Number of APOKASC high LOGG_DW stars: ',len(highg))
        jc.extend(i1[highg])
        lowz=np.where((apokasc['FE_H_ADOP_COR'][i2] < -1.) & (apokasc['FE_H_ADOP_COR'][i2] > -90.))[0]
        print('Number of APOKASC low [Fe/H] stars: ',len(lowz))
        jc.extend(i1[lowz])
        mklinks(data,jc,dir+'_apokasc')
        all.extend(jc)
    
    if coolstars :
        jc=[]
        j=np.where(data['FIELD'] == 'GALCEN')[0]
        print('Number of GALCEN stars: ',len(j))
        jc.extend(j)
        stars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/coolstars.txt',names=['id'],format='fixed_width_no_header')
        i1,i2=match.match(data['APOGEE_ID'],stars['id'])
        print('Number of cool stars: ',len(i1))
        jc.extend(i1)
        mklinks(data,jc,dir+'_galcen')
        all.extend(jc)

    if cal1m :
        j=np.where(data['FIELD'] == 'calibration')[0]
        print('Number of 1m calibration stars: ',len(j))
        all.extend(j)
        j=np.where(data['FIELD'] == 'RCB')[0]
        print('Number of 1m RCB stars: ',len(j))
        all.extend(j)
        mklinks(data,j,dir+'_cal1m')

    if optical is not None:
        #stars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/validation_stars_DR16.txt',names=['id'],format='fixed_width_no_header')
        stars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/'+optical,names=['id'],format='fixed_width_no_header')
        i1,i2=match.match(data['APOGEE_ID'],stars['id'])
        print('Number of optical validation stars: ',len(i1))
        mklinks(data,i1,dir+'_optical')
        all.extend(i1)

    if Ce :
        stars = np.loadtxt(os.environ['APOGEE_DIR']+'/data/calib/Ce_test_stars.txt',dtype=str)[:,0]
        i1,i2=match.match(data['APOGEE_ID'],stars)
        print('Number of Ce test stars: ',len(i1))
        mklinks(data,i1,dir+'_Ce')
        all.extend(i1)

    if ns :
        #north-south overlap
        jn=np.where((data['FIELD'] == 'N2243') | (data['FIELD'] == '000+08') |
                    (data['FIELD'] == '300+75') | (data['FIELD'] == 'M12-N') )[0]
        js=np.where((data['FIELD'] == 'N2243-S') | (data['FIELD'] == '000+08-S') |
                    (data['FIELD'] == '300+75-S') | (data['FIELD'] == 'M12-S') )[0]
        i1,i2=match.match(data['APOGEE_ID'][jn], data['APOGEE_ID'][js])
        jc=list(jn[i1])
        jc.extend(js[i2])
        stars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/apogee_overlap.txt',names=['id'],format='fixed_width_no_header')
        jn=np.where(data['TELESCOPE'] == 'apo25m')[0]
        js=np.where(data['TELESCOPE'] == 'lco25m')[0]
        i1,i2=match.match(data['APOGEE_ID'][jn],stars['id'])
        jc.extend(jn[i1])
        i1,i2=match.match(data['APOGEE_ID'][js],stars['id'])
        jc.extend(js[i2])
        mklinks(data,jc,dir+'_ns')
        print('Number of N/S overlap stars: ',len(jc))
        all.extend(jc)

    if special is not None:
        stars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/'+special,names=['id'],format='fixed_width_no_header')
        jn=np.where(data['TELESCOPE'] == 'apo25m')[0]
        js=np.where(data['TELESCOPE'] == 'lco25m')[0]
        i1,i2=match.match(data['APOGEE_ID'][jn],stars['id'])
        jc=list(jn[i1])
        i1,i2=match.match(data['APOGEE_ID'][js],stars['id'])
        jc.extend(js[i2])
        print('Number of '+special+' stars: ',len(i1))
        mklinks(data,jc,dir+'_special')
        all.extend(jc)

    if ebvmax is not None:
        jc=np.where( (data['SFD_EBV']>0) & (data['SFD_EBV'] < ebvmax) )[0]
        print('Number of low E(B-V) stars: ',len(jc))
        mklinks(data,jc,dir+'_ebv')
        all.extend(jc)

    if hrdata is not None:
        i1, i2 = hrsample(data,hrdata)
        print('Number of HR sample stars: ',len(i1))
        mklinks(data,i1,dir+'_hr')
        all.extend(i1)

    # create "all" directories with all stars, removing duplicates
    print('Total number of stars: ',len(list(set(all))))
    if mkall: mklinks(data,list(set(all)),dir+'_all')

    return indata

def mklinks(data,j,out,n=48) :
    """ Create links in n different output directories for requested indices
    """

    for tel in ['apo1m','apo25m','lco25m'] :
        # create symbolic links in output directories, separate for each instrument
        outdir=tel+'/'+out+'_'+tel
        # remove existing output directories, create new ones
        cleandir(outdir,n)
        gd=np.where(data['TELESCOPE'][j] == tel )[0]
        nsplit=len(gd)//n+1
        for i in range(len(gd)) :
            symlink(data[j[gd[i]]],outdir,i//nsplit)


def cleandir(out,n) :
    '''
    auxiliary routine to clean and remake calibration directories
    '''
    for i in range(n) : 
        try:
            shutil.rmtree('{:s}{:03d}'.format(out,i))
        except : pass
        try:
            os.makedirs('{:s}{:03d}'.format(out,i))
        except : pass

def symlink(data,out,idir) :
    '''
    auxiliary routine to create symlinks to appropriate files from calibration directories
    '''
    outfile='{:s}{:03d}/{:s}.{:s}.fits'.format(
            out,idir,os.path.splitext(os.path.basename(data['FILE']))[0],data['FIELD'])
    if data['TELESCOPE'] == 'apo25m' or data['TELESCOPE'] == 'lco25m' :
        infile='{:s}/{:s}/{:s}'.format(data['TELESCOPE'],data['FIELD'],data['FILE'])
        if not os.path.exists(infile) :
            infile='{:s}/{:d}/{:s}'.format(data['TELESCOPE'],data['LOCATION_ID'],data['FILE'])
    else :
        infile='{:s}/calibration/{:s}'.format(data['TELESCOPE'],data['FILE'])
    os.symlink('../../'+infile,outfile)

def docal(calfile,clobber=False,allstar=True,hr=True,teff=True,logg=True,vmicro=True,vmacro=True,elemcal=True,out=None,stp=False,cal='dr14',calib=False) :
    '''
    Derives all calibration relations and creates plots of them, as requested
    '''

    # output subdirectory
    try: os.mkdir('cal')
    except: pass
    print(os.getcwd())

    # combine aspcapField files into allCal
    if clobber or not os.path.isfile(calfile) :
        allc=allCal(['hr???/aspcapField*.fits','cal???/aspcapField*.fits','clust???/aspcapField*.fits'],out=calfile)
     
    allc=fits.open(calfile)
    if allstar :
        c=fits.open(calfile)
    else :
        c=allc
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
        allcal['giant_teffcal'] = teffcomp.ghb(c[1].data,ebvmax=0.02,glatmin=10,out=out+'giant_teffcomp',yr=[-750,750],dwarf=False,calib=calib)
        if out is not None :
            struct.wrfits(struct.dict2struct(allcal['giant_teffcal']),out+'giant_tecal.fits')
        allcal['dwarf_teffcal'] = teffcomp.ghb(c[1].data,ebvmax=0.02,glatmin=10,trange=[4000,7500],out=out+'dwarf_teffcomp',yr=[-750,750],dwarf=True,calib=calib)
        if out is not None :
            struct.wrfits(struct.dict2struct(allcal['dwarf_teffcal']),out+'dwarf_tecal.fits')
        if stp : pdb.set_trace()
    figs.append(['giant_teffcomp.jpg','dwarf_teffcomp.jpg'])
    ytitle.append('Teff')
    figs.append(['giant_teffcomp_b.jpg','dwarf_teffcomp_b.jpg'])
    ytitle.append('Teff')

    # log g vs asteroseismic
    if logg :
        allcal['rgbrcsep' ] = loggcomp.rcrgb(c[1].data,out=out+'rcrgbsep')
        allcal['loggcal'] = loggcomp.apokasc(c[1].data,plotcal=False,out=out+'loggcomp',calib=calib)
        if out is not None :
            struct.wrfits(struct.dict2struct(dict(allcal['rgbrcsep'].items()+allcal['loggcal'].items())),out+'loggcal.fits')
        if stp : pdb.set_trace()
    figs.append(['rcrgbsep.jpg','none.jpg'])
    ytitle.append('log g')
    figs.append(['loggcomp_b.jpg','loggcomp.jpg'])
    ytitle.append('log g')
    loggcomp.dwarf(c[1].data,out=out+'logg')
    figs.append(['logg_dwarfs.png','logg_all.png'])
    ytitle.append('log g ')

    # vmicro calibration
    if vmicro :
        print("vmicro fit, cubic in log g, linear in [M/H]")
        print("sample limited to FERRE vmicro error <0.01 ")
        vfit.fit_vmicro(c[1].data,degree=3,reject=0.15,mhrange=[-2,1],loggrange=[-0.3,4.9],vmrange=[0,7],
                        teffrange=[3550,6500],vrange=[0.55,4],maxerr=0.01,func=vfit.vm3_1,out=out+'vmicro3_1')

        print("full sample (maxerr=0.1)")
        vfit.fit_vmicro(c[1].data,degree=3,reject=0.15,mhrange=[-2,1],loggrange=[-0.3,4.9],vmrange=[0,7],
                        teffrange=[3500,6500],vrange=[0.55,4],maxerr=0.1,func=vfit.vm3_1,out=out+'vmicro3_1all')

        #dwarfs only
        print("dwarfs only, fit as f(Teff)")
        dw=np.where(c[1].data['FPARAM'][:,1] > 4)[0]
        vfit.fit_vmicro(c[1].data,reject=0.15,mhrange=[-2,1],loggrange=[4,5],vmrange=[0,7],teffrange=[3500,8000],
                        vrange=[0.55,4],maxerr=0.1,func=vfit.vm1t,out=out+'vmicro1t')
        fig,ax=plots.multi(1,1)
        plots.plotc(ax,c[1].data['FPARAM'][dw,0],10**c[1].data['FPARAM'][dw,2],c[1].data['FPARAM'][dw,3],xr=[3500,8000],xt='Teff',
           yr=[0,4],yt='vmicro',zr=[-2,0.5],zt='[M/H]',colorbar=True)

    # vmacro
    if vmacro :
        vfit.fit_vmacro(c[1].data,mhrange=[-2.5,1],reject=0.3,maxerr=0.1,out=out+'vmacro_2d')

    # elemental abundances
    if elemcal :
        elems=np.append(c[3].data['ELEM_SYMBOL'][0],['M','alpha'])
        # use allCal file for uncertainty calibration, so we have multiple visits
        # use allstar for calibration, so we have full solar circle sample
        errcal=elem.cal(allc,c[3].data['ELEM_SYMBOL'][0],c[3].data['ELEMTOH'][0],elems,hard=out+'giants_',cal=cal,errpar=True,plot=False,calib=calib)
        allcal['giantcal']=elem.cal(c,c[3].data['ELEM_SYMBOL'][0],c[3].data['ELEMTOH'][0],elems,hard=out+'giants_',cal=cal,errpar=True,calib=calib)
        allcal['giantcal']['errpar']=errcal['errpar']
        if out is not None :
            struct.wrfits(allcal['giantcal'],out+'giantcal.fits')
        errcal=elem.cal(allc,c[3].data['ELEM_SYMBOL'][0],c[3].data['ELEMTOH'][0],elems,hard=out+'dwarfs_',dwarfs=True,errpar=True,plot=False,cal=cal,calib=calib)
        allcal['dwarfcal']=elem.cal(c,c[3].data['ELEM_SYMBOL'][0],c[3].data['ELEMTOH'][0],elems,hard=out+'dwarfs_',dwarfs=True,cal=cal,calib=calib)
        allcal['dwarfcal']['errpar']=errcal['errpar']
        if out is not None :
            struct.wrfits(allcal['dwarfcal'],out+'dwarfcal.fits')
        if stp : pdb.set_trace()
    figs.append(['giants_all.jpg','dwarfs_all.jpg'])
    ytitle.append('clusters')
    figs.append(['giants_allsolar.jpg','dwarfs_allsolar.jpg'])
    ytitle.append('solar circle')
    figs.append(['giants_M.jpg','dwarfs_M.jpg'])
    ytitle.append('cluster [M/H]')

    html.htmltab(figs,xtitle=['giants','dwarfs'],ytitle=ytitle,file=out+calfile.replace('.fits','.html'))
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

