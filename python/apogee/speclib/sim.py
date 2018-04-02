import pdb
import os
import math
import numpy as np
from astropy.io import fits
from numpy import isclose

def cval(x) :
    ''' routine to convert value to "Kurucz-style" string'''
    if x < 0 :
      prefix = 'm'
    else :
      prefix = 'p'
    return prefix+'{:02d}'.format(int(round(abs(x)*10.)))

def filename(a,c,n,v) :
    return 'a'+cval(a)+'c'+cval(c)+'n'+cval(n)+'v'+cval(v)

def vector(header,axis) :
    caxis='{:1d}'.format(axis)
    return header['CRVAL'+caxis]+header['CDELT'+caxis]*np.arange(header['NAXIS'+caxis])

rootdir=os.environ['APOGEE_SPECLIB']+'/synth/turbospec/kurucz/giantisotopes/tgGK_150714_lsfcombo5_l31c/'

write=True
random=0.01

allcol_apogee_id=[]
allcol_telescope=[]
allcol_field=[]
allcol_afe=[]
allcol_cfe=[]
allcol_nfe=[]
allcol_vmicro=[]
allcol_teff=[]
allcol_logg=[]
allcol_mh=[]
param=[]
w0=4.179
dw=6.e-6
nw=8575
w=w0+np.arange(nw)*dw
ntot=0
for afe in np.arange(-0.5,1.00,0.25)  :
  for cfe in np.arange(-0.5,0.75,0.25) :
    for nfe in np.arange(-0.5,1.1,0.5) :
      for logvmicro in np.arange(math.log10(0.5),math.log10(8),math.log10(2)) :
        vmicro=10.**logvmicro
        file = filename(afe,cfe,nfe,vmicro) 
        try: os.mkdir(file)    
        except: pass
        synth=fits.open(rootdir+file+'.fits')[0]
        wsynth=vector(synth.header,1)/math.log(10)
        col_apogee_id=[]
        col_telescope=[]
        col_field=[]
        col_afe=[]
        col_cfe=[]
        col_nfe=[]
        col_vmicro=[]
        col_teff=[]
        col_logg=[]
        col_mh=[]
        n=0
        for it,teff in enumerate(vector(synth.header,2)) :
            for ig,logg in enumerate(vector(synth.header,3)) :
              if abs((teff-3500)*4/2000. - logg)  < 2 or (isclose(cfe,0.) and isclose(nfe,0.) and isclose(afe,0.) and isclose(vmicro,2.)) :
                for im,mh in enumerate(vector(synth.header,4)) :
                  if np.random.random() < random or (isclose(cfe,0.) and isclose(nfe,0.) and isclose(afe,0.) and isclose(vmicro,2.)) :
                    #print(teff,logg,mh)
                    out=np.interp(w,wsynth,synth.data[im,ig,it,:])
                    err=np.ones(nw)*0.005
                    mask=np.zeros(nw).astype(int)
                    bd=np.where(np.isnan(out))[0]
                    err[bd]=float('nan')
                    apogee_id='{:s}t{:4d}g{:s}m{:s}'.format(file,int(teff),cval(logg),cval(mh))
                    outfile='apStar-sim-'+apogee_id
                    #print(outfile)
                    hdulist=fits.HDUList()
                    hdu=fits.PrimaryHDU()
                    hdu.header['APOGEE_ID']=apogee_id
                    hdu.header['SNR']=1000.
                    hdu.header['SNRVIS1']=1000.
                    hdu.header['FIELD']=file
                    hdu.header['LOCID']=0
                    hdu.header['TEFF']=teff
                    hdu.header['LOGG']=logg
                    hdu.header['MH']=mh
                    hdulist.append(hdu)
                    for data in [out,err,mask] :
                        hdu=fits.ImageHDU(data)
                        hdu.header['CRVAL1']=w0
                        hdu.header['CRPIX1']=1
                        hdu.header['CDELT1']=dw
                        hdulist.append(hdu)
                    if write : hdulist.writeto(file+'/'+outfile+'.fits',overwrite=True)
                    n+=1
                    ntot+=1
                    col_apogee_id.append(apogee_id)
                    col_telescope.append('')
                    col_field.append(file)
                    col_afe.append(afe)
                    col_cfe.append(cfe)
                    col_nfe.append(nfe)
                    col_vmicro.append(vmicro)
                    col_teff.append(teff)
                    col_logg.append(logg)
                    col_mh.append(mh)
                    allcol_apogee_id.append(apogee_id)
                    allcol_telescope.append('')
                    allcol_field.append(file)
                    allcol_afe.append(afe)
                    allcol_cfe.append(cfe)
                    allcol_nfe.append(nfe)
                    allcol_vmicro.append(vmicro)
                    allcol_teff.append(teff)
                    allcol_logg.append(logg)
                    allcol_mh.append(mh)
                    param.append([teff,logg,vmicro,mh,cfe,nfe,afe])
        print(file,' n: ',n)
        hdu = fits.BinTableHDU.from_columns([
                fits.Column(name='APOGEE_ID',format='A64',array=col_apogee_id),
                fits.Column(name='TELESCOPE',format='A10',array=col_telescope),
                fits.Column(name='FIELD',format='A20',array=col_field),
                fits.Column(name='ALPHA_FE',format='E',array=col_afe),
                fits.Column(name='C_FE',format='E',array=col_cfe),
                fits.Column(name='N_FE',format='E',array=col_nfe),
                fits.Column(name='VMICRO',format='E',array=col_vmicro),
                fits.Column(name='TEFF',format='E',array=col_teff),
                fits.Column(name='LOGG',format='E',array=col_logg),
                fits.Column(name='M_H',format='E',array=col_mh) 
              ] )
        if write: hdu.writeto(file+'/apField-'+file+'.fits',overwrite=True)
hdu = fits.BinTableHDU.from_columns([
        fits.Column(name='APOGEE_ID',format='A64',array=allcol_apogee_id),
        fits.Column(name='TELESCOPE',format='A10',array=allcol_telescope),
        fits.Column(name='SNR',format='E',array=np.ones(len(allcol_teff))*1000),
        fits.Column(name='GLON',format='E',array=np.zeros(len(allcol_teff))),
        fits.Column(name='GLAT',format='E',array=np.zeros(len(allcol_teff))),
        fits.Column(name='FPARAM',format='7E',array=param),
        fits.Column(name='PARAM',format='7E',array=param),
        fits.Column(name='FIELD',format='A20',array=allcol_field),
        fits.Column(name='ALPHA_FE',format='E',array=allcol_afe),
        fits.Column(name='C_FE',format='E',array=allcol_cfe),
        fits.Column(name='N_FE',format='E',array=allcol_nfe),
        fits.Column(name='VMICRO',format='E',array=allcol_vmicro),
        fits.Column(name='TEFF',format='E',array=allcol_teff),
        fits.Column(name='LOGG',format='E',array=allcol_logg),
        fits.Column(name='M_H',format='E',array=allcol_mh) 
      ] )
if write : hdu.writeto('allStar.fits',overwrite=True)
print('Ntot: ', ntot)
