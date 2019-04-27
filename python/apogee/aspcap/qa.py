import numpy as np
import os
import pdb
from astropy.io import fits
from tools import plots
from tools import html
from tools import match
import matplotlib.pyplot as plt
from apogee.aspcap import aspcap
from apogee.aspcap import err
from apogee.utils import apload
from apogee.utils import apselect

def plotparams(a,title=None,hard=None) :
    """ Plot parameters vs Teff
    """
    fig,ax=plots.multi(1,8,hspace=0.001)

    paramnames,tagnames,flagnames = params()

    for i in range(8) :
        plots.plotc(ax[i],a['FPARAM'][:,0],a['FPARAM'][:,i],a['FPARAM'][:,3],yt=tagnames[i],xt='Teff')
    if title is not None : fig.suptitle(title)
    if hard is not None : fig.savefig(hard)

def plotcn(hdulist,title=None,out=None) :
    """ Compare parameter level C/M and N/M with element level
    """
    a=hdulist[1].data
    grid=[]
    yt=[]
    for iel,el in enumerate(['C','N']) :
        row=[]
        xt=[]
        yt.append(el)
        for icol in range(5) :
            if icol == 0 :
                te=[3000,4000]
                logg=[0,3.8]
            elif icol == 1 :
                te=[4000,4500]
                logg=[0,3.8]
            elif icol == 2 :
                te=[4500,8000]
                logg=[0,3.8]
            elif icol == 3 :
                te=[3000,4000]
                logg=[3.8,5.5]
            elif icol == 4 :
                te=[4000,8000]
                logg=[3.8,5.5]
            gd = np.where((a['FPARAM'][:,0] >= te[0]) & (a['FPARAM'][:,0] <= te[1]) &
                          (a['FPARAM'][:,1] >= logg[0]) & (a['FPARAM'][:,1] <= logg[1]) )[0]
            if el == 'C' : 
               try: abun=a['FELEM'][gd,0,0]
               except: abun=a['FELEM'][gd,0]
               pabun=a['FPARAM'][gd,4]
            elif el == 'N' :
               try: abun=a['FELEM'][gd,0,2]
               except: abun=a['FELEM'][gd,2]
               pabun=a['FPARAM'][gd,5]
   
            fig,ax=plots.multi(1,2,hspace=0.001)
            plots.plotc(ax[0],a['FPARAM'][gd,3],abun-pabun,a['FPARAM'][gd,0],yr=[-0.5,0.5],zr=te,xt='[M/H]',colorbar=True,zt='Teff',yt=r'$\Delta$['+el+'/M] (el-param)')
            ax[0].text(0.1,0.9,'uncalibrated params',transform=ax[0].transAxes)
            plots.plotc(ax[1],a['FPARAM'][gd,3],abun-pabun,a['SNR'][gd],yr=[-0.5,0.5],zr=[50,200],xt='[M/H]',colorbar=True,zt='S/N',yt=r'[$\Delta$'+el+'/M] (el-param)')
            if out is not None :
                outfile=out+'param_elem_'+el+'_{:1d}.png'.format(icol)
                fig.savefig(outfile)
                row.append(os.path.basename(outfile))
            else: pdb.set_trace()
            plt.close(fig)
            xt.append('{:6.0f}&lt;Teff&lt;{:6.0f} {:6.1f}&lt;logg&lt;{:6.1f}'.format(te[0],te[1],logg[0],logg[1]))
        grid.append(row)
    html.htmltab(grid,file=out+'cn.html',ytitle=yt,xtitle=xt)

def plotelems(hdulist,title=None,out=None) :
    """ Make [X/M] vs [M/H] plots for all elements as f(Teff, logg)
    """
    a=hdulist[1].data
    if 'EXTRATARG' in a.columns.names :
        gd=np.where(a['EXTRATARG'] == 0)[0]
        a=a[gd]
        comment=', main sample only'
    else :
        comment=', full sample'

    els=hdulist[3].data['ELEM_SYMBOL'][0]
    etoh=hdulist[3].data['ELEMTOH'][0]
    grid=[]
    yt=[]
    for iel,el in enumerate(els) :
        row=[]
        xt=[]
        yt.append(el)
        for icol in range(5) :
            if icol == 0 :
                te=[3000,4000]
                logg=[0,3.8]
            elif icol == 1 :
                te=[4000,4500]
                logg=[0,3.8]
            elif icol == 2 :
                te=[4500,8000]
                logg=[0,3.8]
            elif icol == 3 :
                te=[3000,4000]
                logg=[3.8,5.5]
            elif icol == 4 :
                te=[4000,8000]
                logg=[3.8,5.5]
            gd = np.where((a['FPARAM'][:,0] >= te[0]) & (a['FPARAM'][:,0] <= te[1]) &
                          (a['FPARAM'][:,1] >= logg[0]) & (a['FPARAM'][:,1] <= logg[1]) )[0]
            print(el,te,logg,len(gd))
            try:
                if etoh[iel] == 1 : abun=a['FELEM'][gd,0,iel]-a['FPARAM'][gd,3]
                else : abun = a['FELEM'][gd,0,iel]
            except:
                if etoh[iel] == 1 : abun=a['FELEM'][gd,iel]-a['FPARAM'][gd,3]
                else : abun = a['FELEM'][gd,iel]
            fig,ax=plots.multi(1,2,hspace=0.001)
            plots.plotc(ax[0],a['FPARAM'][gd,3],abun,a['FPARAM'][gd,0],yr=[-0.5,1],zr=te,xt='[M/H]',colorbar=True,zt='Teff',yt='['+el+'/M]')
            ax[0].text(0.1,0.9,'uncalibrated params'+comment,transform=ax[0].transAxes)
            plots.plotc(ax[1],a['FPARAM'][gd,3],abun,a['SNR'][gd],yr=[-0.5,1],zr=[50,200],xt='[M/H]',colorbar=True,zt='S/N',yt='['+el+'/M]')
            if out is not None :
                outfile=out+el+'_{:1d}.png'.format(icol)
                fig.savefig(outfile)
                row.append(os.path.basename(outfile))
            else: pdb.set_trace()
            plt.close(fig)
            xt.append('{:6.0f}&lt;Teff&lt;{:6.0f} {:6.1f}&lt;logg&lt;{:6.1f}'.format(te[0],te[1],logg[0],logg[1]))
        grid.append(row)
    html.htmltab(grid,file=out+'elem_chem.html',ytitle=yt,xtitle=xt)

def plotparamdiffs(data,bdata,title=None,cal=False,out=None,elem=True) :
    """ Plot parameter differences between two different runs
    """

    a=data[1].data
    b=bdata[1].data
    paramnames,tagnames,flagnames = aspcap.params()
    i1,i2=match.match(a['APOGEE_ID'],b['APOGEE_ID'])
    print('number of matching stars: ',len(i1))

    # parameters first
    if cal : param = 'PARAM'
    else : param = 'FPARAM'
    grid=[]
    yt=[]
    for i in range(8) :
        if i == 0 : 
            yr=[-200,200]
        elif i==1  : 
            xr=[-0.5,5]
            yr=[-0.5,0.5]
        else :
            xr=[-2.5,1.0]
            yr=[-0.5,0.5]
        row=[]
        yt.append(paramnames[i])
        for j in range(3) :
            fig,ax=plots.multi(1,1)
            if j == 0 :
                plots.plotc(ax,a[param][i1,0],b[param][i2,i]-a[param][i1,i],a[param][i1,3],yt=r'$\Delta$'+tagnames[i],xt='Teff',yr=yr,xr=[3000,8000],zr=[-2,0.5])
            elif j == 1 :
                plots.plotc(ax,a[param][i1,1],b[param][i2,i]-a[param][i1,i],a[param][i1,3],yt=r'$\Delta$'+tagnames[i],xt='log g',yr=yr,xr=[-1,6],zr=[-2,0.5])
            elif j == 2 :
                plots.plotc(ax,a[param][i1,3],b[param][i2,i]-a[param][i1,i],a[param][i1,3],yt=r'$\Delta$'+tagnames[i],xt='[M/H]',yr=yr,xr=[-2.5,1.0],zr=[-2,0.5])
            if out is not None:
                outfile = out+'paramdiffs_{:1d}_{:1d}.png'.format(i,j)
                fig.savefig(outfile)
                row.append(os.path.basename(outfile))
            else: 
                pdb.set_trace()
            plt.close(fig)
        grid.append(row)
    ptab=html.table(grid,ytitle=yt)
    html.htmltab(grid,file=out+'paramdiffs.html',ytitle=yt)

    # now elements
    if elem :
        grid=[]
        yt=[]
        elemnames = data[3].data['ELEM_SYMBOL'][0]
        belemnames = bdata[3].data['ELEM_SYMBOL'][0]
        if cal : elem = 'ELEM'
        else : elem = 'FELEM'
        for i,el in enumerate(elemnames) :
          ii = np.where(belemnames == el)[0]
          if len(ii) > 0 : 
            yr=[-0.5,0.5]
            row=[]
            yt.append(el)
            if len(a[elem].shape) == 3 : abun=a[elem][i1,0,i]
            else : abun=a[elem][i1,i]
            if len(b[elem].shape) == 3 : abun_b=b[elem][i2,0,ii]
            else : abun_b=b[elem][i2,i]
            for j in range(3) :
                fig,ax=plots.multi(1,1)
                if j == 0 :
                    plots.plotc(ax,a[param][i1,0],abun_b-abun,a[param][i1,3],yt=r'$\Delta$'+el,xt='Teff',yr=yr,xr=[3000,8000],zr=[-2,0.5])
                elif j == 1 :
                    plots.plotc(ax,a[param][i1,1],abun_b-abun,a[param][i1,3],yt=r'$\Delta$'+el,xt='log g',yr=yr,xr=[-1,6],zr=[-2,0.5])
                elif j == 2 :
                    plots.plotc(ax,a[param][i1,3],abun_b-abun,a[param][i1,3],yt=r'$\Delta$'+el,xt='[M/H]',yr=yr,xr=[-2.5,1.0],zr=[-2,0.5])
                if out is not None:
                    outfile = out+el+'_diff_{:1d}.png'.format(j)
                    fig.savefig(outfile)
                    row.append(os.path.basename(outfile))
                else: 
                    pdb.set_trace()
                plt.close(fig)
            grid.append(row)
        etab=html.table(grid,ytitle=yt)
        html.htmltab(grid,file=out+'elemdiffs.html',ytitle=yt)

    # HR diagrams
    grid=[]
    row=[]
    aspcap.hr(a[i1],hard=out+'hr_match1.png',grid=True,size=1)
    row.append(os.path.basename(out+'hr_match1.png'))
    aspcap.hr(b[i2],hard=out+'hr_match2.png',grid=True,size=1)
    row.append(os.path.basename(out+'hr_match2.png'))
    grid.append(row)
    row=[]
    aspcap.hr(a[i1],hard=out+'hrcal_match1.png',param='PARAM',grid=True,size=1)
    row.append(os.path.basename(out+'hrcal_match1.png'))
    aspcap.hr(b[i2],hard=out+'hrcal_match2.png',param='PARAM',grid=True,size=1)
    row.append(os.path.basename(out+'hrcal_match2.png'))
    grid.append(row)
    hrtab=html.table(grid)

    fp=html.head(out+'diffs.html')
    fp.write('<h2> HR diagrams, raw(top), calibrated(bottom)')
    fp.write(hrtab)
    fp.write('<h2> Parameter differences as f(Teff, logg, [M/H]')
    fp.write(ptab)
    fp.write('<h2> Abundance differences as f(Teff, logg, [M/H]')
    if elem: fp.write(etab)
    html.tail(fp)
    return 
    
def dr14comp(a,out=None,elem=True) :
    """ Comparisons to DR14
    """
    apl=apload.ApLoad(dr='dr14')
    dr14=apl.allStar()
    plotparamdiffs(a,dr14,out=out+'dr14_',elem=elem)

def repeat(data,out='./',elem=True,logg=[-1,6]) :
    """ Comparison of repeat observations of objects
    """

    gd=apselect.select(data[1].data,badval='STAR_BAD',raw=True,logg=logg)

    a=data[1].data[gd]
    stars = set(a['APOGEE_ID'])
    colors=['r','g','b']
    tels=['apo1m','apo25m','lco25m'] 
    diff=[]
    ediff=[]
    teff=[]
    mh=[]
    sn=[]

    snbins=np.arange(50,300,50)
    # looop over stars looking for duplications
    for star in stars :
        j = np.where(a['APOGEE_ID'] == star)[0]
        n = len(j)
        if n > 1 :
            print(star,n,a['APOGEE_ID'][j])
            # if we have a duplicate, we need duplicate observations at comparable S/N
            for isn in range(len(snbins)-1) :
                jj = np.where((a['SNR'][j] > snbins[isn]) & (a['SNR'][j] < snbins[isn+1]) )[0]
                if len(jj) > 1 :
                    for jjj in jj :
                        teff.append(a['FPARAM'][j,0].mean())
                        mh.append(a['FPARAM'][j,3].mean())
                        sn.append(a['SNR'][j[jj]].mean())
                        diff.append(a['FPARAM'][j[jjj],:]-a['FPARAM'][j[jj],:].mean(axis=0))
                        try:
                            if elem:  ediff.append(a['FELEM'][j[jjj],0,:]-a['FELEM'][j[jj],0,:].mean(axis=0))
                        except:
                            if elem:  ediff.append(a['FELEM'][j[jjj],:]-a['FELEM'][j[jj],:].mean(axis=0))
    teff=np.array(teff)
    mh=np.array(mh)
    sn=np.array(sn)
    j=np.where(sn>250)[0]
    sn[j]=249.999
    diff=np.array(diff)
    grid=[]
    params=data[3].data['PARAM_SYMBOL'][0]
    for i,param in enumerate(params) :
        if param == 'TEFF' : 
            zr=[0,100] 
            gd=np.where(abs(diff[:,i]) < 500.)[0]
        else : 
            zr=[0,0.2]
            gd=np.where(abs(diff[:,i]) < 1.)[0]
        err.errfit(teff[gd],sn[gd],mh[gd],diff[gd,i],out=out+param,zr=zr,verbose=False,mkhtml=False)
        grid.append([os.path.basename(out+param+'_err.jpg'),os.path.basename(out+param+'_err_sn.jpg')])
    html.htmltab(grid,file=out+'repeat_params.html',ytitle=params)
    if elem : 
        els=data[3].data['ELEM_SYMBOL'][0]
        ediff=np.array(ediff)
        grid=[]
        for i,el in enumerate(els) :
            gd=np.where(abs(ediff[:,i]) < 1.)[0]
            err.errfit(teff[gd],sn[gd],mh[gd],ediff[gd,i],out=out+el,quad=True)
            grid.append([os.path.basename(out+el+'_err.jpg'),os.path.basename(out+el+'_err_sn.jpg')])
        html.htmltab(grid,file=out+'repeat_elem.html',ytitle=els,mkhtml=False)

#l33=fits.open('allCal-r12-l33p.fits')
#notie=fits.open('../l33notie/allCal-r12-l33notie.fits')
#i1,i2=match.match(l33[1].data['APOGEE_ID'],notie[1].data['APOGEE_ID'])
#bd=np.where((l33[1].data['FPARAM'][i1,0] < 4000) & (np.abs(l33[1].data['FELEM'][i1,0,17]-notie[1].data['FELEM'][i2,0,17]) > 0.2) )[0]
#pdb.set_trace()
#
#els=l33[3].data['ELEM_SYMBOL'][0]
#etoh=l33[3].data['ELEMTOH'][0]
#
#grid=[]
#ytit=[]
#for iel,el in enumerate(els) :
#  for col in [0,1] :
#    if col == 0 : gd=np.where((l33[1].data['FPARAM'][:,1] < 3.8) & (l33[1].data['FPARAM'][:,0] > 3050) )[0]
#    else :  gd=np.where((l33[1].data['FPARAM'][:,1] > 3.8) & (l33[1].data['FPARAM'][:,0] > 3050) )[0]
#    fig,ax=plots.multi(1,3,hspace=0.001)
#    if etoh[iel] == 1 : abun=l33[1].data['FELEM'][gd,0,iel]-l33[1].data['FPARAM'][gd,3]
#    else : abun = l33[1].data['FELEM'][gd,0,iel]
#    plots.plotc(ax[0],l33[1].data['FPARAM'][gd,3],abun,l33[1].data['FPARAM'][gd,0],yr=[-0.5,1],zr=[3500,5500],xt='[M/H]')
#    ax[0].text(0.1,0.9,'uncalibrated params',transform=ax[0].transAxes)
#
#    if etoh[iel] == 1 : abun=l33[1].data['FELEM_CAL'][gd,0,iel]-l33[1].data['FPARAM'][gd,3]
#    else : abun = l33[1].data['FELEM_CAL'][gd,0,iel]
#    plots.plotc(ax[1],l33[1].data['FPARAM'][gd,3],abun,l33[1].data['FPARAM'][gd,0],yr=[-0.5,1],zr=[3500,5500],xt='[M/H]')
#    ax[1].text(0.1,0.9,'calibrated params',transform=ax[1].transAxes)
#
#    if col == 0 : gd=np.where((notie[1].data['FPARAM'][:,1] < 3.8) & (l33[1].data['FPARAM'][:,0] > 3050) )[0]
#    else : gd=np.where((notie[1].data['FPARAM'][:,1] > 3.8) & (l33[1].data['FPARAM'][:,0] > 3050) )[0]
#    if etoh[iel] == 1 : abun=notie[1].data['FELEM'][gd,0,iel]-notie[1].data['FPARAM'][gd,3]
#    else : abun = notie[1].data['FELEM'][gd,0,iel]
#    plots.plotc(ax[2],notie[1].data['FPARAM'][gd,3],abun,notie[1].data['FPARAM'][gd,0],yr=[-0.5,1],zr=[3500,5500],xt='[M/H]')
#    ax[2].text(0.1,0.9,'old notie',transform=ax[2].transAxes)
#
#    fig.savefig('elem/comp_'+el+'_{:d}.png'.format(col))
#    plt.close(fig)
#  grid.append(['comp_'+el+'_0.png','comp_'+el+'_1.png'])
#  ytit.append(el)
#
#xtit=['giants','dwarfs']
#html.htmltab(grid,file='elem/comp.html',ytitle=ytit,xtitle=xtit)
