from apogee.utils import apload
from apogee.utils import apselect
from apogee.utils import bitmask
from apogee.aspcap import elem
from apogee.aspcap import teffcomp
from apogee.aspcap import loggcomp
from apogee.aspcap import aspcap
from apogee.aspcap import err
from apogee.aspcap import qa
from tools import html
from tools import match
from tools import struct
from tools import plots
from tools import fit
try: from tools import vfit
except: pass
import copy
import os
import shutil
import pdb
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy.table import Column
from astropy.coordinates import SkyCoord
import astropy.units as units

os.environ['ISOCHRONE_DIR']='/uufs/chpc.utah.edu/common/home/apogee/isochrones/'

def allField(search=['apo*/*/a?Field-*.fits','apo*/*/a?FieldC-*.fits','lco*/*/a?Field-*.fits'],out='allField.fits',verbose=False) :
    '''
    Concatenate set of apField files
    '''

    if type(search) == str:
        search=[search]
    allfiles=[]
    for path in search :
        allfiles.extend(glob.glob(path))

    a=[]
    for file in allfiles :
        if 'Field-cal_' not in file :
            dat=fits.open(file)[1].data
            if type(dat['STARFLAG'][0]) is not np.uint64 :
                print(file, type(dat['STARFLAG'][0]))
            a.append(dat)
    all = np.hstack(a)

    bd=np.where((all['RA'] < 0.0001) & (all['DEC'] < 0.0001) )[0]
    starmask=bitmask.StarBitMask()
    rvfail=starmask.getval('RV_FAIL')
    for i in bd:
        all['STARFLAG'][i] |= np.int64(rvfail)

    # concatenate the structures
    #all=struct.concat(search,verbose=verbose)

    # write out the file
    if out is not None:
        print('writing',out)
        struct.wrfits(all,out)

    return all

def m67(all) :

    m67=[]
    for a in all :
        m67.append(np.array(apselect.clustmember(a,'M67',raw=True)))

    els = aspcap.elems()[0]
    fig,ax=plots.multi(1,len(all),hspace=0.001)
    for iel,el in enumerate(els) :
        for i,(a,j) in enumerate(zip(all,m67)) :
            ax[i].cla()
            plots.plotc(ax[i],a['FPARAM'][j,1],a['FELEM'][j,iel],a['FPARAM'][j,0],yt=el,
                    xr=[0,5],yr=[-0.3,0.3])
        pdb.set_trace()    
    return m67
    
def allCal(search=['clust???/aspcapField-*.fits','cal???/aspcapField-*.fits'],nelem=15,out='allCal.fits',allfield=None) :
    '''
    Concatenate aspcapField files, adding ELEM tags if not there
    '''
    # concatenate the structures
    if type(search) == str:
        search=[search]
    allfiles=[]
    for path in search :
        allfiles.extend(glob.glob(path))

    a=[]
    for file in allfiles :
        print(file)
        a.append(fits.open(file)[1].data)
    all = np.hstack(a)

    return all

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

def dr16(allstar='allStar-r12-l33cal.fits',allcal='allCal-r12-l33.fits') :
    """ Run the summary routines for DR16 allStar and allCal
    """
    if allstar is not None : summary(out=allstar,prefix='allStar/',cal='dr16',repeat=False)
    if allcal is not None : summary(out=allcal,prefix='allCal/',cal='dr16',calib=False)

def writecal(allstar='allStar',allcal='allCal') :
    """ Write the calibration FITS files into output calibration directory
        Replace errpar parameters with those from repeats
    """
    # replace the scatter coefficients in calibration relations with these
    giant_errfit=fits.open(allcal+'/repeat/giant_errfit.fits')
    dwarf_errfit=fits.open(allcal+'/repeat/dwarf_errfit.fits')
    giant_abuncal=Table.read(allstar+'/calib/giant_abuncal.fits')
    dwarf_abuncal=Table.read(allstar+'/calib/dwarf_abuncal.fits')
    giant_abuncal.remove_column('errpar')
    dwarf_abuncal.remove_column('errpar')
    errfit=np.vstack([giant_errfit[2].data['ERRFIT'],giant_errfit[1].data['ERRFIT'][3:7:3,:]])
    giant_abuncal.add_column(Column(errfit, name='errpar'))
    errfit=np.vstack([dwarf_errfit[2].data['ERRFIT'],dwarf_errfit[1].data['ERRFIT'][3:7:3,:]])
    dwarf_abuncal.add_column(Column(errfit, name='errpar'))
    shutil.copy(allstar+'/calib/all_tecal.fits','cal/')
    shutil.copy(allstar+'/calib/giant_loggcal.fits','cal/')
    shutil.copy(allstar+'/calib/dwarf_loggcal.fits','cal/')
    giant_abuncal.write('cal/giant_abuncal.fits',overwrite=True)
    dwarf_abuncal.write('cal/dwarf_abuncal.fits',overwrite=True)


def summary(out='allCal.fits',prefix='allcal/',cal='dr16',hr=True,repeat=True,dr14comp=True,teff=True,logg=True,elemcal=True, calib=True, doqa=True, calibrate=True) :
    """ Create QA summary page and plots
    """
    hdulist=fits.open(out)
    all=hdulist[1].data

    try: os.mkdir(prefix)
    except: pass
    try: os.mkdir(prefix+'hr/')
    except: pass
    try: os.mkdir(prefix+'calib/')
    except: pass
    try: os.mkdir(prefix+'calibrated/')
    except: pass
    try: os.mkdir(prefix+'optical/')
    except: pass
    try: os.mkdir(prefix+'qa/')
    except: pass
    try: os.mkdir(prefix+'qa_calibrated/')
    except: pass
    try: os.mkdir(prefix+'repeat/')
    except: pass

    # HR diagrams
    if hr :
        qa.hr(all,hard=prefix+'hr/hr.png',xr=[8000,3000],grid=True,iso=[9.0,10.0],alpha=1.0,snrbd=5,target=prefix+'hr/hr',size=1)
        qa.hr(all,hard=prefix+'hr/hrhot.png',xr=[20000,3000],iso=[8.0,10.0],snrbd=30,size=1)
        qa.multihr(all,hard=prefix+'hr/multihr.png',size=1)
        qa.hr(all,hard=prefix+'hr/hr_cal.png',xr=[8000,3000],grid=True,iso=[9.0,10.0],alpha=1.0,snrbd=5,param='PARAM',target=prefix+'hr/hr_cal',size=1)
    grid=[[prefix+'hr/hr.png',prefix+'hr/multihr.png',prefix+'hr/hrhot.png'],
          [prefix+'hr/hr_main.png',prefix+'hr/hr_targ.png',''],
          [prefix+'hr/hr_cal_main.png',prefix+'hr/hr_cal.png','']] 

    # Master summary HTML file
    f=html.head(file=out.replace('.fits','.html'))
    f.write(html.table(grid,ytitle=['uncalibrated','uncalibrated','calibrated']))
    f.write('<br>Uncalibrated parameters:<br>')
    ids = ['VESTA','2M14153968+1910558']
    j=[]
    try: 
        for id in ids: j.extend( np.where( (np.core.defchararray.strip(all['APOGEE_ID']) == id) & (all['VISIT'] == 0)) [0] )
    except: 
        for id in ids: j.extend( np.where( (np.core.defchararray.strip(all['APOGEE_ID']) == id) ) [0] )
    ids = ['VESTA','Arcturus']
    f.write(html.table(all['FPARAM'][j],plots=False,ytitle=ids,xtitle=aspcap.params()[1]))
    f.write('<br>calibrated parameters:<br>')
    f.write(html.table(all['PARAM'][j],plots=False,ytitle=ids,xtitle=aspcap.params()[1]))
    # table of abundances (relative to M)
    f.write('<br>Uncalibrated abundances:<br>')
    try: abun=all['FELEM'][j,0,:]
    except: abun=all['FELEM'][j,:]
    xtit=[]
    for i in range(len(hdulist[3].data['ELEM_SYMBOL'][0])) : 
        if hdulist[3].data['ELEMTOH'][0][i] == 1 : abun[:,i]-=all['FPARAM'][j,3]
        xtit.append('['+hdulist[3].data['ELEM_SYMBOL'][0][i]+'/M]')
    f.write(html.table(abun,plots=False,ytitle=ids,xtitle=xtit))

    f.write('<br>calibrated abundances:<br>')
    xtit=[]
    for i in range(len(hdulist[3].data['ELEM_SYMBOL'][0])) : 
        xtit.append('['+hdulist[3].data['ELEM_SYMBOL'][0][i]+'/M]')
    f.write(html.table(all['X_M'][j],plots=False,ytitle=ids,xtitle=xtit))

    f.write('<p> Calibration relations<ul>\n')
    f.write('<li> <a href='+prefix+'calib/'+out.replace('.fits','.html')+'> Calibration plots</a>\n')
    f.write('<li> <a href='+prefix+'calibrated/'+out.replace('.fits','.html')+'> Calibration check (calibration plots from calibrated values) </a>\n')
    f.write('<li> <a href='+prefix+'qa/calib.html> Calibrated-uncalibrated plots</a>\n')
    f.write('</ul>\n')
    f.write('<p> Comparisons<ul>\n')
    f.write('<li> <a href='+prefix+'optical/optical.html> Comparison with optical abundances</a>\n')
    f.write('<li> <a href='+prefix+'qa/apolco.html> APO-LCO comparison</a>\n')
    f.write('<li> <a href='+prefix+'qa/m67.html> M67 abundances</a>\n')
    f.write('</ul>\n')
    f.write('<p> Chemistry plots<ul>\n')
    f.write('<li> <a href='+prefix+'qa/elem_chem.html> Chemistry plots with uncalibrated abundances</a>\n')
    f.write('<li> <a href='+prefix+'qa_calibrated/elem_chem.html> Chemistry plots with calibrated abundances, main sample</a>\n')
    f.write('<li> <a href='+prefix+'qa_calibrated/all_elem_chem.html> Chemistry plots with calibrated abundances, all</a>\n')
    f.write('<li> <a href='+prefix+'qa_calibrated/named_elem_chem.html> Chemistry plots with calibrated abundances, named tags</a>\n')
    f.write('</ul>\n')
    f.write('<p> QA checks<ul>\n')
    f.write('<li> <a href='+prefix+'qa/cn.html> C,N parameters vs abundances</a>\n')
    f.write('<li> <a href='+prefix+'qa/flags.html> Bitmasks</a>\n')
    f.write('<li> <a href='+prefix+'qa/elem_errs.html> Abundance uncertainties</a>\n')
    f.write('<li> <a href='+prefix+'qa/dr14_diffs.html> DR14 comparison plots</a>')
    f.write('</ul>\n')
    f.write('<br> Duplicates/repeats, for empirical uncertainties: <ul>\n')
    f.write('<li> <a href='+prefix+'repeat/giant_repeat_elem.html> Elemental abundances, giants</a>\n')
    f.write('<li> <a href='+prefix+'repeat/giant_repeat_param.html> Parameters,  giants</a>\n')
    f.write('<li> <a href='+prefix+'repeat/dwarf_repeat_elem.html> Elemental abundances, dwarfs</a>\n')
    f.write('<li> <a href='+prefix+'repeat/dwarf_repeat_param.html> Parameters,  dwarfs</a>\n')
    f.write('</ul>\n')
    html.tail(f)

    # optical comparison index
    grid=[]
    grid.append(['r12_uncal_paramcomp.png','r12_cal_paramcomp.png'])
    grid.append(['dr14_uncal_paramcomp.png','dr14_cal_paramcomp.png'])
    yt=['DR16 Parameters','DR14 parameters']
    for el in hdulist[3].data['ELEM_SYMBOL'][0] :
        grid.append(['r12_uncal_abundcomp_{:s}.png'.format(el),'r12_cal_abundcomp_{:s}.png'.format(el)])
        yt.append(el)
    html.htmltab(grid,file=prefix+'optical/optical.html',ytitle=yt,xtitle=['uncalibrated','calibrated'])
   
    # do the calibration and calibration and QA plots
    if calibrate :
        allcal=docal(out,clobber=False,hr=False,teff=teff,logg=logg,vmicro=False,vmacro=False,elemcal=elemcal,
              out=prefix+'calib/',stp=False,cal=cal,calib=False) 
    if repeat :
        # get scatter from repeat observations
        giant_param_errfit,giant_elem_errfit=err.repeat(hdulist,out=prefix+'repeat/giant_',elem=elemcal,logg=[-1,3.8])
        dwarf_param_errfit,dwarf_elem_errfit=err.repeat(hdulist,out=prefix+'repeat/dwarf_',elem=elemcal,logg=[3.8,5.5])

    # now get plots for calibrated data
    if calib : 
        docal(out,clobber=False,hr=False,teff=teff,logg=logg,vmicro=False,vmacro=False,elemcal=elemcal,
              out=prefix+'calibrated/',stp=False,cal=cal,calib=True) 
    hdulist=fits.open(out)
    if dr14comp :
        qa.dr14comp(hdulist,out=prefix+'qa/',elem=elemcal)
    if doqa : 
        qa.plotelems(hdulist,out=prefix+'qa/')
        qa.plotelem_errs(hdulist,out=prefix+'qa/')
        qa.plotelems(hdulist,calib=True,out=prefix+'qa_calibrated/')
        qa.plotelems(hdulist,out=prefix+'qa_calibrated/all_',main=False)
        qa.plotelems(hdulist,out=prefix+'qa_calibrated/named_',named=True)
        qa.plotcn(hdulist,out=prefix+'qa/')
        qa.calib(hdulist,out=prefix+'qa/')
        qa.m67(hdulist,out=prefix+'qa/')
        qa.apolco(hdulist,out=prefix+'qa/')
        qa.flags(hdulist,out=prefix+'qa/')

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

def calsample(indata=None,file='clust.html',plot=True,clusters=True,apokasc='APOKASC_cat_v4.4.2',apred=None,
              cal1m=True,coolstars=True,dir='cal',hrdata=None,optical='cal_stars_20190329.txt',ns=True,special=None,Ce=True,ebvmax=None,snmin=75,mkindiv=False,mkall=True) :
    '''
    selects a calibration subsample from an input apField structure, including several calibration sub-classes: 
        cluster, APOKASC stars, 1m calibration stars. Creates cluster web pages/plots if requested
    '''

    if apred is None : 
        print('need to provide apred')
        pdb.set_trace()

    if indata is None :
        indata=allField(['apo25m/*/a?Field-*.fits','lco25m/*/a?Field-*.fits','apo1m/calibration/a?Field-*.fits'],out=None,verbose=True)

    j=np.where((indata['COMMISS'] == 0) & (indata['SNR'] > 75) )[0]
    data=indata[j]

    try: os.mkdir(dir)
    except: pass
    all=[]
    if clusters :
        jc=[]
        fstars=open(dir+'/allclust.txt','w')
        f=html.head(file=dir+'/'+file)
        f.write('<A HREF=allclust.txt> cluster stars list </a>')
        f.write('<TABLE BORDER=2>\n')
        f.write('<TR><TD>NAME<TD>RA<TD>DEC<TD>Radius<TD>RV<TD>Delta RV<TD>Position criterion<TD>RV criterion<TD>PM criterion<TD>Parallax criterion<TD> CMD')
        clust=apselect.clustdata(gals=False)
        for ic in range(len(clust.name)) :
            print(clust[ic].name)
            j=apselect.clustmember(data,clust[ic].name,plot=False,hard=None,gals=False)
            print(clust[ic].name,len(j))
            # clusters to exclude here
            if (clust[ic].name not in ['OmegaCen','Pal1','Pal6','Pal5','Terzan12'])  and (len(j) >= 5): jc.extend(j)
            f.write('<TR><TD><A HREF='+clust[ic].name+'.txt>'+clust[ic].name+'</A><TD>{:12.6f}<TD>{:12.6f}<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                    clust[ic].ra,clust[ic].dec,clust[ic].rad,clust[ic].rv,clust[ic].drv))
            f.write('<TD><A HREF='+clust[ic].name+'_pos.png><IMG SRC='+clust[ic].name+'_pos.png width=300></A>\n')
            f.write('<TD><A HREF='+clust[ic].name+'_rv.png><IMG SRC='+clust[ic].name+'_rv.png width=300></A>\n')
            f.write('<TD><A HREF='+clust[ic].name+'_pm.png><IMG SRC='+clust[ic].name+'_pm.png width=300></A>\n')
            f.write('<TD><A HREF='+clust[ic].name+'_parallax.png><IMG SRC='+clust[ic].name+'_parallax.png width=300></A>\n')
            f.write('<TD><A HREF='+clust[ic].name+'_cmd.png><IMG SRC='+clust[ic].name+'_cmd.png width=300></A>\n')
            np.savetxt(dir+'/'+clust[ic].name+'.txt',data[jc]['APOGEE_ID'],fmt='%s')
            for star in data[jc]['APOGEE_ID'] : fstars.write('{:s} {:s}\n'.format(star,clust[ic].name))
        html.tail(f)
        fstars.close()
        print('Number of cluster stars: ',len(jc))
        if mkindiv: mklinks(data,jc,dir+'_clust',apred=apred)
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
        if mkindiv: mklinks(data,jc,dir+'_apokasc',apred=apred)
        all.extend(jc)
    
    if coolstars :
        jc=[]
        j=np.where((data['FIELD'] == 'GALCEN') or (data['FIELD'] == b'GALCEN'))[0]
        print('Number of GALCEN stars: ',len(j))
        jc.extend(j)
        stars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/coolstars.txt',names=['id'],format='fixed_width_no_header')
        i1,i2=match.match(data['APOGEE_ID'],stars['id'])
        print('Number of cool stars: ',len(i1))
        jc.extend(i1)
        if mkindiv: mklinks(data,jc,dir+'_galcen',apred=apred)
        all.extend(jc)

    if cal1m :
        j=np.where((data['FIELD'] == 'calibration') or (data['FIELD'] == b'calibration'))[0]
        print('Number of 1m calibration stars: ',len(j))
        all.extend(j)
        j=np.where(data['FIELD'] == 'RCB')[0]
        print('Number of 1m RCB stars: ',len(j))
        all.extend(j)
        if mkindiv: mklinks(data,j,dir+'_cal1m',apred=apred)

    if optical is not None:
        #stars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/validation_stars_DR16.txt',names=['id'],format='fixed_width_no_header')
        stars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/'+optical,names=['id'],format='fixed_width_no_header')
        i1,i2=match.match(data['APOGEE_ID'],stars['id'])
        print('Number of optical validation stars: ',len(i1))
        if mkindiv: mklinks(data,i1,dir+'_optical',apred=apred)
        all.extend(i1)

    if Ce :
        stars = np.loadtxt(os.environ['APOGEE_DIR']+'/data/calib/Ce_test_stars.txt',dtype=str)[:,0]
        i1,i2=match.match(data['APOGEE_ID'],stars)
        print('Number of Ce test stars: ',len(i1))
        if mkindiv: mklinks(data,i1,dir+'_Ce',apred=apred)
        all.extend(i1)

    if ns :
        #north-south overlap
        jn=np.where((data['FIELD'] == 'N2243') | (data['FIELD'] == b'N2243') | (data['FIELD'] == '000+08')  | (data['FIELD'] == b'000+08') |
                    (data['FIELD'] == '300+75') | (data['FIELD'] == b'300+75') | (data['FIELD'] == 'M12-N') | (data['FIELD'] == b'300+75') )[0]
        js=np.where((data['FIELD'] == 'N2243-S') | (data['FIELD'] == b'N2243-S') | (data['FIELD'] == '000+08-S') | (data['FIELD'] == b'N2243-S') | 
                    (data['FIELD'] == '300+75-S') | (data['FIELD'] == b'300+75-S') |  (data['FIELD'] == 'M12-S') | (data['FIELD'] == b'M12-S') )[0]
        i1,i2=match.match(data['APOGEE_ID'][jn], data['APOGEE_ID'][js])
        jc=list(jn[i1])
        jc.extend(js[i2])
        stars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/apogee_overlap.txt',names=['id'],format='fixed_width_no_header')
        if type(data['TELESCOPE'][0]) is str :
            jn=np.where(data['TELESCOPE'] == 'apo25m')[0]
            js=np.where(data['TELESCOPE'] == 'lco25m')[0]
        else :
            jn=np.where(data['TELESCOPE'] == b'apo25m')[0]
            js=np.where(data['TELESCOPE'] == b'lco25m')[0]
        i1,i2=match.match(data['APOGEE_ID'][jn],stars['id'])
        jc.extend(jn[i1])
        i1,i2=match.match(data['APOGEE_ID'][js],stars['id'])
        jc.extend(js[i2])
        if mkindiv: mklinks(data,jc,dir+'_ns',apred=apred)
        print('Number of N/S overlap stars: ',len(jc))
        all.extend(jc)

    if special is not None:
        stars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/'+special,names=['id'],format='fixed_width_no_header')
        if type(data['TELESCOPE'][0]) is str :
            jn=np.where(data['TELESCOPE'] == 'apo25m')[0]
            js=np.where(data['TELESCOPE'] == 'lco25m')[0]
        else :
            jn=np.where(data['TELESCOPE'] == b'apo25m')[0]
            js=np.where(data['TELESCOPE'] == b'lco25m')[0]
        i1,i2=match.match(data['APOGEE_ID'][jn],stars['id'])
        jc=list(jn[i1])
        i1,i2=match.match(data['APOGEE_ID'][js],stars['id'])
        jc.extend(js[i2])
        print('Number of '+special+' stars: ',len(i1))
        if mkindiv: mklinks(data,jc,dir+'_special',apred=apred)
        all.extend(jc)

    if ebvmax is not None:
        jc=np.where( (data['SFD_EBV']>0) & (data['SFD_EBV'] < ebvmax) )[0]
        print('Number of low E(B-V) stars: ',len(jc))
        if mkindiv: mklinks(data,jc,dir+'_ebv',apred=apred)
        all.extend(jc)

    if hrdata is not None:
        i1, i2 = hrsample(data,hrdata)
        print('Number of HR sample stars: ',len(i1))
        if mkindiv: mklinks(data,i1,dir+'_hr',apred=apred)
        all.extend(i1)

    # create "all" directories with all stars, removing duplicates
    print('Total number of stars: ',len(list(set(all))))
    if mkall: mklinks(data,list(set(all)),dir+'_all',apred=apred)

    return indata,all

def mklinks(data,j,out,n=24,apred=None) :
    """ Create links in n different output directories for requested indices
    """

    for tel in ['apo1m','apo25m','lco25m'] :
        # create symbolic links in output directories, separate for each instrument
        outdir=tel+'/'+out+'_'
        # remove existing output directories, create new ones
        cleandir(outdir,n)
        if type(data['TELESCOPE'][0]) is str :
            gd=np.where(data['TELESCOPE'][j] == tel )[0]
        else: 
            gd=np.where(data['TELESCOPE'][j] == tel.encode() )[0]
        nsplit=len(gd)//n+1
        load=apload.ApLoad(apred=apred,telescope=tel)
        if len(gd) > 0 :
            # create apField files
            ii=0
            for i in range(n) :
                field=os.path.basename(outdir)+'{:03d}'.format(i)
                apfield=load.filename('Field',field=field)
                tmp=Table(data[np.array(j)[gd[ii:ii+nsplit]]])
                # we will add original FIELD to APOGEE_ID, so increase column width
                tmp['APOGEE_ID']=tmp['APOGEE_ID'].astype('U40')
                for star in tmp :
                    star['APOGEE_ID']=star['APOGEE_ID']+'.'+star['FIELD']
                    star['FIELD'] = field
                tmp.write(apfield)
                ii+=nsplit
            # create apStar links
            for i in range(len(gd)) :
                symlink(data[j[gd[i]]],outdir,i//nsplit,load=load)


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

def tostr(dat) :
    if type(dat) is str : return dat
    else : return dat.decode()

def symlink(data,out,idir,load=None) :
    '''
    auxiliary routine to create symlinks to appropriate files from calibration directories
    '''
    if data['FILE'] == '' :
        data['FILE'] = os.path.basename(load.filename('Star',field=data['FIELD'],obj=data['APOGEE_ID']))
    try:
        outfile='{:s}{:03d}/{:s}.{:s}.fits'.format(
                out,idir,os.path.splitext(os.path.basename(tostr(data['FILE'])))[0],tostr(data['FIELD']))
    except :
        outfile='{:s}{:03d}/{:s}.{:s}.fits'.format(
                out,idir,os.path.splitext(os.path.basename(tostr(data['FILE'])))[0],tostr(data['FIELD']))
    if tostr(data['TELESCOPE']) == 'apo25m' or tostr(data['TELESCOPE']) == 'lco25m' :
        infile='{:s}/{:s}/{:s}'.format(tostr(data['TELESCOPE']),tostr(data['FIELD']),tostr(data['FILE']))
        if not os.path.exists(infile) :
            infile='{:s}/{:d}/{:s}'.format(tostr(data['TELESCOPE']),data['LOCATION_ID'],tostr(data['FILE']))
    else :
        infile='{:s}/calibration/{:s}'.format(tostr(data['TELESCOPE']),tostr(data['FILE']))
    os.symlink('../../'+infile,outfile)

def docal(infile,clobber=False,hr=True,teff=True,logg=True,vmicro=True,vmacro=True,elemcal=True,out=None,stp=False,cal='dr14',calib=False,ebvmax=0.02) :
    '''
    Derives all calibration relations and creates plots of them, as requested
    '''

    # output subdirectory
    try: os.mkdir('cal')
    except: pass
    print(os.getcwd())

    c=fits.open(infile)

    # if we don't have GLON/GLAT, add it
    ebvmax=0.05
    if (c[1].data['GLON'].max()<1) :
        coords=SkyCoord(ra=c[1].data['RA']*units.degree,dec=c[1].data['DEC']*units.degree) 
        c[1].data['GLON'] = coords.galactic.l
        c[1].data['GLAT'] = coords.galactic.b


    figs=[]
    ytitle=[]
    # HR diagram
    if hr :
        print('HR diagram...')
        reload(apselect)
        fig,ax=plots.multi(1,1)
        if calib : param='PARAM'
        else : param='FPARAM'
        plots.plotc(ax,c[param][:,0],c[param][:,1],c[param][:,3],xr=[6000,3000],yr=[5,-1],zr=[-2,0.5])
        plt.savefig(out+'hr.png')                                                                                                 
        figs.append(['hr.png','hr.png'])
        ytitle.append('HR')

    allcal={}
    # Teff vs photometric
    if teff :
        print('Teff calibration...')
        allcal['teffcal'] = teffcomp.ghb(c[1].data,ebvmax=ebvmax,glatmin=10,out=out+'tecal',yr=[-750,750],trange=[4500,7000],loggrange=[-1,6],calib=calib)
        allcal['giant_teffcal'] = teffcomp.ghb(c[1].data,ebvmax=ebvmax,glatmin=10,out=out+'giant_tecal',yr=[-750,750],loggrange=[-1,3.8],calib=calib)
        allcal['dwarf_teffcal'] = teffcomp.ghb(c[1].data,ebvmax=ebvmax,glatmin=10,trange=[4500,7000],out=out+'dwarf_tecal',yr=[-750,750],loggrange=[3.8,6],calib=calib)
        if out is not None :
            struct.wrfits(struct.dict2struct(allcal['teffcal']),out+'all_tecal.fits')
            struct.wrfits(struct.dict2struct(allcal['giant_teffcal']),out+'giant_tecal.fits')
            struct.wrfits(struct.dict2struct(allcal['dwarf_teffcal']),out+'dwarf_tecal.fits')
        if stp : pdb.set_trace()
    figs.append(['tecal.png','tecal_b.png'])
    ytitle.append('Teff all together')
    figs.append(['giant_tecal.png','dwarf_tecal.png'])
    ytitle.append('Teff, giants and dwarfs')
    figs.append(['giant_tecal_b.png','dwarf_tecal_b.png'])
    ytitle.append('Teff, giants and dwarfs')

    # log g vs asteroseismic
    if logg :
        print('log g calibration...')
        allcal['rgbrcsep' ] = loggcomp.rcrgb(c[1].data,out=out+'rcrgbsep')
        allcal['giant_loggcal'] = loggcomp.apokasc(c[1].data,plotcal=False,out=out+'rcrgb_loggcal',calib=calib)
        allcal['dwarf_loggcal'] = loggcomp.dwarf(c[1].data,out=out+'logg',calib=calib)
        if out is not None :
            # following is Python 2
            #struct.wrfits(struct.dict2struct(dict(allcal['rgbrcsep'].items()+allcal['giant_loggcal'].items())),
            struct.wrfits(struct.dict2struct({**allcal['rgbrcsep'],**allcal['giant_loggcal']}),
                            out+'giant_loggcal.fits')
            struct.wrfits(struct.dict2struct(allcal['dwarf_loggcal']),out+'dwarf_loggcal.fits')
        if stp : pdb.set_trace()
    figs.append(['rcrgbsep.png','rcrgb_loggcal.png'])
    ytitle.append('log g, RGB/RC')
    figs.append(['rcrgb_loggcal_b.png','logg_dwarfs.png'])
    ytitle.append('log g, RGB/RC and dwarfs')
    figs.append(['logg_all.png','logg_all.png'])
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
        print('vmacro relation...')
        vfit.fit_vmacro(c[1].data,mhrange=[-2.5,1],reject=0.3,maxerr=0.1,out=out+'vmacro_2d')

    # elemental abundances
    if elemcal :
        print('abundances ...')
        elems=np.append(c[3].data['ELEM_SYMBOL'][0],['M','alpha'])
        allcal['giant_abuncal']=elem.cal(c,c[3].data['ELEM_SYMBOL'][0],c[3].data['ELEMTOH'][0],elems,hard=out+'giants_',cal=cal,errpar=True,calib=calib)
        allcal['dwarf_abuncal']=elem.cal(c,c[3].data['ELEM_SYMBOL'][0],c[3].data['ELEMTOH'][0],elems,hard=out+'dwarfs_',dwarfs=True,cal=cal,calib=calib)
        if out is not None :
            struct.wrfits(allcal['giant_abuncal'],out+'giant_abuncal.fits')
            struct.wrfits(allcal['dwarf_abuncal'],out+'dwarf_abuncal.fits')
        if stp : pdb.set_trace()
    figs.append(['giants_all.png','dwarfs_all.png'])
    ytitle.append('clusters')
    figs.append(['giants_allsolar.png','dwarfs_allsolar.png'])
    ytitle.append('solar circle')
    figs.append(['giants_M.png','dwarfs_M.png'])
    ytitle.append('cluster [M/H]')
    figs.append(['giants_clust_key.png','dwarfs_clust_key.png'])
    ytitle.append('cluster ID')

    html.htmltab(figs,xtitle=['giants','dwarfs'],ytitle=ytitle,file=out+infile.replace('.fits','.html'))
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
            y.append(run+'/'+run+out+plot+'.png')
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
        plt.savefig(out+'.png')

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
        plt.savefig(out+'_dvmicro.png')

def starcomp(ref='l31a',comps=['l31b','l31b_vm4','l30b_vm4','l31a_asset'],out=None) :
    '''
    Creates web page for comparison of mulitple version
    '''
    grid=[]
    for comp in comps :
        compstars(ref,comp,out='comp/'+ref+'_'+comp)
        y=[ref+'_'+comp+'.png',ref+'_'+comp+'_dvmicro.png']
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
             y.append('../'+ver+'/cal/'+el+'_'+plot+'.png')
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

def logg(a,caldir='cal/') :
    """ Check IDL corrections for log g
    """
    #a=fits.open('allStar-r12-l33.fits')[1].data

    aspcapmask=bitmask.AspcapBitMask()
    parammask=bitmask.ParamBitMask()
    starmask=bitmask.StarBitMask()
    gd=np.where( ((a['ASPCAPFLAG']&aspcapmask.badval()) == 0) )[0]

    cal=fits.open(caldir+'/giant_loggcal.fits')[1].data
    rgbsep=cal['rgbsep'][0]
    cnsep=cal['cnsep'][0]
    rclim=cal['rclim'][0]
    rcfit2=cal['rcfit2'][0]
    rgbfit2=cal['rgbfit2'][0]
    calloggmin=cal['calloggmin']
    calloggmax=cal['calloggmax']
    calteffmin=cal['calteffmin']
    calteffmax=cal['calteffmax']

    # for stars that aren't bad, get cn and dt
    cn=a['FPARAM'][gd,4]-a['FPARAM'][gd,5]
    dt=a['FPARAM'][gd,0] - (rgbsep[0] + rgbsep[1]*(a['FPARAM'][gd,1]-2.5) +rgbsep[2]*a['FPARAM'][gd,3])

    new=np.zeros(len(a))-9999.99

    # select RC
    rc=np.where((a['FPARAM'][gd,1]<rclim[1])&(a['FPARAM'][gd,1]>rclim[0])&
                (cn>cnsep[0]+cnsep[1]*a['FPARAM'][gd,3] + cnsep[2]*dt)&
                (a['FPARAM'][gd,1]<calloggmax)&(a['FPARAM'][gd,1]>calloggmin) &
                (a['FPARAM'][gd,0]<calteffmax)&(a['FPARAM'][gd,0]>calteffmin))[0]
    rccorr=rcfit2[0] + rcfit2[1]*a['FPARAM'][gd,1] + rcfit2[2]*a['FPARAM'][gd,1]**2
    new[gd[rc]]=a['FPARAM'][gd[rc],1]-rccorr[rc]
    #rcidl=np.where( (a['PARAMFLAG'][gd,1]&parammask.getval('LOGG_CAL_RC')) >0)[0]

    # select RGB
    rgb=np.where(((a['FPARAM'][gd,1]>rclim[1])|(a['FPARAM'][gd,1]<rclim[0])|
                (cn<cnsep[0]+cnsep[1]*a['FPARAM'][gd,3] + cnsep[2]*dt)) &
                (a['FPARAM'][gd,1]<calloggmax)&(a['FPARAM'][gd,1]>calloggmin) &
                (a['FPARAM'][gd,0]<calteffmax)&(a['FPARAM'][gd,0]>calteffmin))[0]
    #clip logg at loggmin and loggmax
    logg=clip(a['FPARAM'][gd,1],cal['loggmin'],cal['loggmax'])
    mh=clip(a['FPARAM'][gd,3],cal['mhmin'],cal['mhmax'])
    # get correction
    rgbcorr=(rgbfit2[0] + rgbfit2[1]*logg + rgbfit2[2]*logg**2 +
                       rgbfit2[3]*logg**3 + rgbfit2[4]*mh )
    new[gd[rgb]]=a['FPARAM'][gd[rgb],1]-rgbcorr[rgb]
    #rgbidl=np.where( (a['PARAMFLAG'][gd,1]&parammask.getval('LOGG_CAL_RGB')) >0)[0]

    cal=fits.open(caldir+'/dwarf_loggcal.fits')[1].data
    teff=clip(a['FPARAM'][gd,0],cal['temin'],cal['temax'])
    logg=clip(a['FPARAM'][gd,1],cal['loggmin'],cal['loggmax'])
    mh=clip(a['FPARAM'][gd,3],cal['mhmin'],cal['mhmax'])
    msfit=cal['msfit'][0]
    mscorr=msfit[0]+msfit[1]*teff+msfit[2]*mh
    ms=np.where(a['FPARAM'][gd,1] > cal['calloggmin'])[0]
    new[gd[ms]]=a['FPARAM'][gd[ms],1]-mscorr[ms]
    #msidl=np.where( (a['PARAMFLAG'][gd,1]&parammask.getval('LOGG_CAL_MS')) >0)[0]

    trans=np.where((a['FPARAM'][gd,1] < 4) & (a['FPARAM'][gd,1] > 3.5) &
                (a['FPARAM'][gd,0] < calteffmax) )[0]
    ms_weight=(a['FPARAM'][gd[trans],1]-3.5)/0.5
    new[gd[trans]] = a['FPARAM'][gd[trans],1]-(mscorr[trans]*ms_weight+rgbcorr[trans]*(1-ms_weight))

    diff =a['PARAM'][:,1]-new
    bd = np.where (np.isclose(diff,0.,1.e-6,0.01) == False)[0]
    return new

def clip(x,xmin,xmax) :
    """ clip input array so that x<xmin becomes xmin, x>xmax becomes xmax, return clipped array
    """
    new=copy.copy(x)
    bd=np.where(x<xmin)[0]
    new[bd]=xmin
    bd=np.where(x>xmax)[0]
    new[bd]=xmax
    return new



