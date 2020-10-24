import numpy as np
import os
import pdb
from astropy.io import fits
from tools import plots
from tools import html
from tools import match
from tools import fit
import matplotlib
import matplotlib.pyplot as plt
from apogee.aspcap import aspcap
from apogee.aspcap import err
from apogee.speclib import isochrones
from apogee.utils import apload
from apogee.utils import apselect
from apogee.utils import bitmask

def plotparams(a,title=None,hard=None) :
    """ Plot parameters vs Teff
    """
    fig,ax=plots.multi(1,8,hspace=0.001)

    paramnames,tagnames,flagnames = params()

    for i in range(8) :
        plots.plotc(ax[i],a['FPARAM'][:,0],a['FPARAM'][:,i],a['FPARAM'][:,3],yt=tagnames[i],xt='Teff')
    if title is not None : fig.suptitle(title)
    if hard is not None : fig.savefig(hard)

te_ranges=[[3000,4000],[4000,4500],[4500,8000],[3000,4000],[4000,8000]]
logg_ranges=[[0,3.8],[0,3.8],[0,3.8],[3.8,5.5],[3.8,5.5]]

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
        icol=0
        for te,logg in zip(te_ranges,logg_ranges) :
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
            icol+=1
            xt.append('{:6.0f}&lt;Teff&lt;{:6.0f} {:6.1f}&lt;logg&lt;{:6.1f}'.format(te[0],te[1],logg[0],logg[1]))
        grid.append(row)
    html.htmltab(grid,file=out+'cn.html',ytitle=yt,xtitle=xt)

def elemvslogg(hdulist,title=None,out=None,calib=False,main=True,named=False) :
    """ Abundances vs log g for solar sample
    """
    a=hdulist[1].data
    if main and ('EXTRATARG' in a.columns.names) :
        gd=np.where(a['EXTRATARG'] == 0)[0]
        a=a[gd]
        comment=', main sample only'
    else :
        comment=', full sample'

    if calib : param = 'PARAM'
    else : param = 'FPARAM'

    els=hdulist[3].data['ELEM_SYMBOL'][0]
    etoh=hdulist[3].data['ELEMTOH'][0]

    # for plots vs logg of solar neighborhood solar stars
    solar=np.where((a['gaia_parallax_error']/np.abs(a['gaia_parallax']) < 0.1) )[0]
    distance = 1000./a['gaia_parallax'][solar]
    x,y,z,r=lbd2xyz(a['GLON'][solar],a['GLAT'][solar],distance/1000.)
    gd = np.where((abs(z) < 0.5) & (r>8) & (r<9))[0]
    solar=solar[gd]
    gd = np.where((a[param][solar,3] >= -0.1) & (a[param][solar,3] <= 0.1) )[0]
    solar=solar[gd]
    gels=list(els)
    for el in ['Fe','Ge','Rb','Nd','Yb','Ce'] : gels.remove(el)
    igel=0
    if len(gels)%2 == 0 : gfig,gax=plots.multi(2,len(gels)//2,wspace=0.4,hspace=0.001,figsize=(9,12))
    else : gfig,gax=plots.multi(2,len(gels)//2+1,hspace=0.001,wspace=0.001,figsize=(9,12))
    
    for iel,el in enumerate(els) :
        if named:
            if el == 'Fe' : tag = 'FE_H'
            else : tag=(el+'_FE').upper()
            abun=a[tag]
            ytit='['+el+'/Fe]'
            rabun=a['X_M'][:,iel]+a['M_H']-a['FE_H']
        elif calib :
            abun=a['X_M'][:,iel]
            rabun=a['X_M'][:,iel]
            ytit='['+el+'/M] (cal)'
        else :
            try:
                if etoh[iel] == 1 : abun=a['FELEM'][:,0,iel]-a['FPARAM'][:,3]
                else : abun = a['FELEM'][:,0,iel]
            except:
                if etoh[iel] == 1 : abun=a['FELEM'][:,iel]-a['FPARAM'][:,3]
                else : abun = a['FELEM'][:,iel]
            rabun=a['X_M'][:,iel]
            ytit='['+el+'/M] (uncal)'
        if el in gels :
            print(el,igel)
            plots.plotp(gax[igel//2,igel%2],a[param][solar,1],rabun[solar],color='r',
                        xr=[5.9,-0.9],yr=[-0.39,0.39], size=1,xt='log g',yt=ytit,alpha=0.2,rasterized=False)
            plots.plotp(gax[igel//2,igel%2],a[param][solar,1],abun[solar],color='k',
                        xr=[5.9,-0.9],yr=[-0.39,0.39], size=1,xt='log g',yt=ytit,rasterized=False)
            igel+=1

    gfig.savefig(out+'elem_solar_logg.png',dpi=150)
    # C/N
    fig,ax=plots.multi(1,1,figsize=(6,4))
    gd=np.where(a['C_FE'][solar]>-9)[0]
    plots.plotp(ax,a[param][solar[gd],1],a['C_FE'][solar[gd]]-a['N_FE'][solar[gd]],color='k',
                xr=[5.9,-0.9],yr=[-0.59,0.39],
                size=3,xt='log g',yt='[C/N]')
    plt.tight_layout()
    fig.savefig(out+'cn_logg.png',dpi=150)


def plotelems(hdulist,title=None,out=None,calib=False,main=True,named=False) :
    """ Make [X/M] vs [M/H] plots for all elements as f(Teff, logg)
    """
    a=hdulist[1].data
    if main and ('EXTRATARG' in a.columns.names) :
        gd=np.where(a['EXTRATARG'] == 0)[0]
        a=a[gd]
        comment=', main sample only'
    else :
        comment=', full sample'

    if calib : param = 'PARAM'
    else : param = 'FPARAM'
    try: vhelio = a['VHELIO_AVG']
    except: vhelio = a['VHELIO']

    els=hdulist[3].data['ELEM_SYMBOL'][0]
    etoh=hdulist[3].data['ELEMTOH'][0]
    grid=[]
    yt=[]

    # for plots vs logg of solar neighborhood solar stars
    solar=np.where((a['gaia_parallax_error']/np.abs(a['gaia_parallax']) < 0.1) )[0]
    distance = 1000./a['gaia_parallax'][solar]
    x,y,z,r=lbd2xyz(a['GLON'][solar],a['GLAT'][solar],distance/1000.)
    gd = np.where((abs(z) < 0.5) & (r>8) & (r<9))[0]
    solar=solar[gd]
    gd = np.where((a[param][solar,3] >= -0.1) & (a[param][solar,3] <= 0.1) )[0]
    solar=solar[gd]
    gels=list(els)
    for el in ['Fe','Ge','Rb','Nd','Yb','Ce'] : gels.remove(el)
    igel=0
    if len(gels)%2 == 0 : gfig,gax=plots.multi(2,len(gels)//2,wspace=0.5,hspace=0.001,figsize=(8,16))
    else : gfig,gax=plots.multi(2,len(gels)//2+1,hspace=0.001,wspace=0.001,figsize=(8,16))
    
    for iel,el in enumerate(els) :
        if named:
            if el == 'Fe' : tag = 'FE_H'
            else : tag=(el+'_FE').upper()
            abun=a[tag]
            ytit='['+el+'/Fe]'
        elif calib :
            abun=a['X_M'][:,iel]
            ytit='['+el+'/M] (cal)'
        else :
            try:
                if etoh[iel] == 1 : abun=a['FELEM'][:,0,iel]-a['FPARAM'][:,3]
                else : abun = a['FELEM'][:,0,iel]
            except:
                if etoh[iel] == 1 : abun=a['FELEM'][:,iel]-a['FPARAM'][:,3]
                else : abun = a['FELEM'][:,iel]
            ytit='['+el+'/M] (uncal)'
        row=[]
        xt=[]
        yt.append(el)
        icol=0
        for te,logg in zip(te_ranges,logg_ranges) :
            gd = np.where((a[param][:,0] >= te[0]) & (a[param][:,0] <= te[1]) &
                          (a[param][:,1] >= logg[0]) & (a[param][:,1] <= logg[1]) )[0]
            print(el,te,logg,len(gd))

            fig,ax=plots.multi(1,3,hspace=0.001,figsize=(8,8))
            plots.plotc(ax[0],a[param][gd,3],abun[gd],a[param][gd,0],xr=[-2.5,1.0],yr=[-0.5,1],zr=te,
                        xt='[M/H]',colorbar=True,zt='Teff',yt=ytit)
            ax[0].text(0.1,0.9,'uncalibrated params'+comment,transform=ax[0].transAxes)
            plots.plotc(ax[1],a[param][gd,3],abun[gd],a['SNR'][gd],xr=[-2.5,1.0],yr=[-0.5,1],zr=[50,200],
                        xt='[M/H]',colorbar=True,zt='S/N',yt=ytit)
            plots.plotc(ax[2],a[param][gd,3],abun[gd],vhelio[gd],xr=[-2.5,1.0],yr=[-0.5,1],zr=[-200,200],
                        xt='[M/H]',colorbar=True,zt='vhelio',yt=ytit)
            if out is not None :
                outfile=out+el+'_{:1d}.png'.format(icol)
                fig.savefig(outfile)
                row.append(os.path.basename(outfile))
            else: pdb.set_trace()
            plt.close(fig)
            icol+=1
            xt.append('{:6.0f}&lt;Teff&lt;{:6.0f} {:6.1f}&lt;logg&lt;{:6.1f}'.format(te[0],te[1],logg[0],logg[1]))
        #plot as f(Teff) for giants
        gd = np.where((a[param][:,1] >= -1) & (a[param][:,1] <= 3.8) )[0]
        fig,ax=plots.multi(1,3,hspace=0.001,figsize=(8,8))
        plots.plotc(ax[0],a[param][gd,0],abun[gd],a[param][gd,3],xr=[3000,5000],yr=[-0.5,1],zr=[-2,0.5],
                    xt='Teff',colorbar=True,zt='[M/H]',yt=ytit)
        plots.plotc(ax[1],a[param][gd,0],abun[gd],a['SNR'][gd],xr=[3000,5000],yr=[-0.5,1],zr=[50,200],
                    xt='Teff',colorbar=True,zt='S/N',yt=ytit)
        plots.plotc(ax[2],a[param][gd,0],abun[gd],vhelio[gd],xr=[3000,5000],yr=[-0.5,1],zr=[-200,200],
                    xt='Teff',colorbar=True,zt='vhelio',yt=ytit)
        xt.append('{:6.1f}&lt;logg&lt;{:6.1f}'.format(-0.5,3.8))
        outfile=out+el+'_teff.png'
        fig.savefig(outfile)
        plt.close(fig)
        row.append(os.path.basename(outfile))
        grid.append(row)

        if el in gels :
            print(el,igel)
            plots.plotc(gax[igel//2,igel%2],a[param][solar,1],abun[solar],a[param][solar,3],xr=[5.9,-0.9],yr=[-0.39,0.39],zr=[-0.1,0.1],
                        size=1,xt='log g',yt=ytit,zt='[M/H]')
            igel+=1

    gfig.savefig(out+'elem_solar_logg.png')
    if out is not None : html.htmltab(grid,file=out+'elem_chem.html',ytitle=yt,xtitle=xt)

def plotparam_errs(hdulist,title=None,out=None) :
    """ Plot uncertainties
    """
    a=hdulist[1].data
    els=hdulist[3].data['PARAM_SYMBOL'][0]
    grid=[]
    yt=[]
    dwarfs=np.where(a['FPARAM'][:,1] > 3.8)[0]
    giants=np.where(a['FPARAM'][:,1] < 3.8)[0]

    for iel,el in enumerate(els) :
        gfig,gax=plots.multi(1,2,hspace=0.001)
        plots.plotc(gax[0],a['PARAM'][giants,0],a['X_H_ERR'][giants,iel],a[param][giants,3],
                    xr=[3000,7000],zr=[-2.5,1.0],yr=[0.0,0.5],xt='Teff',colorbar=True,zt='[M/H]',yt=el,size=1)
        try :
            plots.plotc(gax[1],a[param][giants,0],a['FELEM_ERR'][giants,0,iel],a[param][giants,3],
                        xr=[3000,7000],zr=[-2.5,1.0],yr=[0.0,0.5],xt='Teff',colorbar=True,zt='[M/H]',yt=el,size=1)
        except :
            plots.plotc(gax[1],a[param][giants,0],a['FELEM_ERR'][giants,iel],a[param][giants,3],
                        xr=[3000,7000],zr=[-2.5,1.0],yr=[0.0,0.5],xt='Teff',colorbar=True,zt='[M/H]',yt=el,size=1)
        dfig,dax=plots.multi(1,2,hspace=0.001)
        plots.plotc(dax[0],a[param][dwarfs,0],a['X_H_ERR'][dwarfs,iel],a[param][dwarfs,3],
                    xr=[3000,7000],zr=[-2.5,1.0],yr=[0.0,0.5],xt='Teff',colorbar=True,zt='[M/H]',yt=el,size=1)
        try :
            plots.plotc(dax[1],a[param][dwarfs,0],a['FELEM_ERR'][dwarfs,0,iel],a[param][dwarfs,3],
                        xr=[3000,7000],zr=[-2.5,1.0],yr=[0.0,0.5],xt='Teff',colorbar=True,zt='[M/H]',yt=el,size=1)
        except :
            plots.plotc(dax[1],a[param][dwarfs,0],a['FELEM_ERR'][dwarfs,iel],a[param][dwarfs,3],
                        xr=[3000,7000],zr=[-2.5,1.0],yr=[0.0,0.5],xt='Teff',colorbar=True,zt='[M/H]',yt=el,size=1)
        if out is not None :
            goutfile=out+'giant_'+el+'_err.png'
            gfig.savefig(goutfile)
            plt.close(gfig)
            doutfile=out+'dwarf_'+el+'_err.png'
            dfig.savefig(doutfile)
            plt.close(dfig)
            grid.append([os.path.basename(goutfile),os.path.basename(doutfile)])
            yt.append(el) 
    if out is not None : html.htmltab(grid,file=out+'elem_errs.html',ytitle=yt,xtitle=['giants','dwarfs'])

def plotelem_errs(hdulist,title=None,out=None,calib=False) :
    """ Plot uncertainties
    """
    a=hdulist[1].data
    els=hdulist[3].data['ELEM_SYMBOL'][0]
    etoh=hdulist[3].data['ELEMTOH'][0]
    grid=[]
    yt=[]
    dwarfs=np.where(a['FPARAM'][:,1] > 3.8)[0]
    giants=np.where(a['FPARAM'][:,1] < 3.8)[0]

    if calib: param='PARAM'
    else :param='FPARAM'
    for iel,el in enumerate(els) :
        gfig,gax=plots.multi(1,2,hspace=0.001)
        plots.plotc(gax[0],a[param][giants,0],a['X_H_ERR'][giants,iel],a[param][giants,3],
                    xr=[3000,7000],zr=[-2.5,1.0],yr=[0.0,0.5],xt='Teff',colorbar=True,zt='[M/H]',yt=el,size=1)
        try :
            plots.plotc(gax[1],a[param][giants,0],a['FELEM_ERR'][giants,0,iel],a[param][giants,3],
                        xr=[3000,7000],zr=[-2.5,1.0],yr=[0.0,0.5],xt='Teff',colorbar=True,zt='[M/H]',yt=el,size=1)
        except :
            plots.plotc(gax[1],a[param][giants,0],a['FELEM_ERR'][giants,iel],a[param][giants,3],
                        xr=[3000,7000],zr=[-2.5,1.0],yr=[0.0,0.5],xt='Teff',colorbar=True,zt='[M/H]',yt=el,size=1)
        dfig,dax=plots.multi(1,2,hspace=0.001)
        plots.plotc(dax[0],a[param][dwarfs,0],a['X_H_ERR'][dwarfs,iel],a[param][dwarfs,3],
                    xr=[3000,7000],zr=[-2.5,1.0],yr=[0.0,0.5],xt='Teff',colorbar=True,zt='[M/H]',yt=el,size=1)
        try :
            plots.plotc(dax[1],a[param][dwarfs,0],a['FELEM_ERR'][dwarfs,0,iel],a[param][dwarfs,3],
                        xr=[3000,7000],zr=[-2.5,1.0],yr=[0.0,0.5],xt='Teff',colorbar=True,zt='[M/H]',yt=el,size=1)
        except :
            plots.plotc(dax[1],a[param][dwarfs,0],a['FELEM_ERR'][dwarfs,iel],a[param][dwarfs,3],
                        xr=[3000,7000],zr=[-2.5,1.0],yr=[0.0,0.5],xt='Teff',colorbar=True,zt='[M/H]',yt=el,size=1)
        if out is not None :
            goutfile=out+'giant_'+el+'_err.png'
            gfig.savefig(goutfile)
            plt.close(gfig)
            doutfile=out+'dwarf_'+el+'_err.png'
            dfig.savefig(doutfile)
            plt.close(dfig)
            grid.append([os.path.basename(goutfile),os.path.basename(doutfile)])
            yt.append(el) 
    if out is not None : html.htmltab(grid,file=out+'elem_errs.html',ytitle=yt,xtitle=['giants','dwarfs'])

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
            fig,ax=plots.multi(1,2,hspace=0.001)
            if j == 0 :
                plots.plotc(ax[0],a['FPARAM'][i1,0],b['FPARAM'][i2,i]-a['FPARAM'][i1,i],a['FPARAM'][i1,3],
                            yt=r'$\Delta$'+tagnames[i],xt='Teff',yr=yr,xr=[3000,8000],zr=[-2,0.5])
                plots.plotc(ax[1],a['PARAM'][i1,0],b['PARAM'][i2,i]-a['PARAM'][i1,i],a['PARAM'][i1,3]
                            ,yt=r'$\Delta$'+tagnames[i],xt='Teff',yr=yr,xr=[3000,8000],zr=[-2,0.5])
            elif j == 1 :
                plots.plotc(ax[0],a['FPARAM'][i1,1],b['FPARAM'][i2,i]-a['FPARAM'][i1,i],a['FPARAM'][i1,3],
                            yt=r'$\Delta$'+tagnames[i],xt='log g',yr=yr,xr=[-1,6],zr=[-2,0.5])
                plots.plotc(ax[1],a['PARAM'][i1,1],b['PARAM'][i2,i]-a['PARAM'][i1,i],a['PARAM'][i1,3],
                            yt=r'$\Delta$'+tagnames[i],xt='log g',yr=yr,xr=[-1,6],zr=[-2,0.5])
            elif j == 2 :
                plots.plotc(ax[0],a['FPARAM'][i1,3],b['FPARAM'][i2,i]-a['FPARAM'][i1,i],a['FPARAM'][i1,3],
                            yt=r'$\Delta$'+tagnames[i],xt='[M/H]',yr=yr,xr=[-2.5,1.0],zr=[-2,0.5])
                plots.plotc(ax[1],a['PARAM'][i1,3],b['PARAM'][i2,i]-a['PARAM'][i1,i],a['PARAM'][i1,3],
                            yt=r'$\Delta$'+tagnames[i],xt='[M/H]',yr=yr,xr=[-2.5,1.0],zr=[-2,0.5])
            ax[0].text(0.1,0.9,'Uncalibrated',transform=ax[0].transAxes)
            ax[1].text(0.1,0.9,'Calibrated',transform=ax[1].transAxes)
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
                fig,ax=plots.multi(1,2)
                if j == 0 :
                    plots.plotc(ax[0],a['FPARAM'][i1,0],abun_b-abun,a['FPARAM'][i1,3],yt=r'$\Delta$'+el,xt='Teff',yr=yr,xr=[3000,8000],zr=[-2,0.5])
                    plots.plotc(ax[1],a['PARAM'][i1,0],abun_b-abun,a['PARAM'][i1,3],yt=r'$\Delta$'+el,xt='Teff',yr=yr,xr=[3000,8000],zr=[-2,0.5])
                elif j == 1 :
                    plots.plotc(ax[0],a['FPARAM'][i1,1],abun_b-abun,a['FPARAM'][i1,3],yt=r'$\Delta$'+el,xt='log g',yr=yr,xr=[-1,6],zr=[-2,0.5])
                    plots.plotc(ax[1],a['PARAM'][i1,1],abun_b-abun,a['PARAM'][i1,3],yt=r'$\Delta$'+el,xt='log g',yr=yr,xr=[-1,6],zr=[-2,0.5])
                elif j == 2 :
                    plots.plotc(ax[0],a['FPARAM'][i1,3],abun_b-abun,a['FPARAM'][i1,3],yt=r'$\Delta$'+el,xt='[M/H]',yr=yr,xr=[-2.5,1.0],zr=[-2,0.5])
                    plots.plotc(ax[1],a['PARAM'][i1,3],abun_b-abun,a['PARAM'][i1,3],yt=r'$\Delta$'+el,xt='[M/H]',yr=yr,xr=[-2.5,1.0],zr=[-2,0.5])
                ax[0].text(0.1,0.9,'Uncalibrated',transform=ax[0].transAxes)
                ax[1].text(0.1,0.9,'Calibrated',transform=ax[1].transAxes)
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
    hr(a[i1],hard=out+'hr_match1.png',grid=True,size=1)
    row.append(os.path.basename(out+'hr_match1.png'))
    hr(b[i2],hard=out+'hr_match2.png',grid=True,size=1)
    row.append(os.path.basename(out+'hr_match2.png'))
    grid.append(row)
    row=[]
    hr(a[i1],hard=out+'hrcal_match1.png',param='PARAM',grid=True,size=1)
    row.append(os.path.basename(out+'hrcal_match1.png'))
    hr(b[i2],hard=out+'hrcal_match2.png',param='PARAM',grid=True,size=1)
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
    
def dr14comp(a,out=None,elem=True,domiss=False) :
    """ Comparisons to DR14
    """
    apl=apload.ApLoad(dr='dr14')
    dr14=apl.allStar()
    plotparamdiffs(a,dr14,out=out+'dr14_',elem=elem)

    if domiss :
        miss=set(dr14[1].data['APOGEE_ID'])-set(a[1].data['APOGEE_ID'])
        print('{:d} stars in DR14 missing from current data'.format(len(miss)))

        bad=[]
        bad1m=[]
        for m in miss :
            j=np.where(dr14[1].data['APOGEE_ID'] == m)[0]
            for v,jj in zip(dr14[1].data['VISITS'][j],j) :
                for vv in v.split(',') :
                    mjd=vv.split('-')[2]
                    try: 
                        if int(mjd) > 55800 : 
                            bad.append(m)
                            #print(m,vv,mjd,dr14[1].data['FIELD'][jj],dr14[1].data['LOCATION_ID'][jj])   
                    except: 
                        bad1m.append(m)
                        #print(m,vv)
        print('not 1m',len(bad),len(set(bad)))
        print('1m',len(bad1m),len(set(bad1m)))


def m67(allstar,out='./') :
    """ M67 abundances
    """
    gd=apselect.select(allstar[1].data,badval='STAR_BAD')
    m67=np.array(apselect.clustmember(allstar[1].data[gd],'M67',raw=True,pm=True,dist=True))
    m67=gd[m67]
    els=allstar[3].data['ELEM_SYMBOL'][0]
    grid=[]
    for iel,el in enumerate(els) :
        fig,ax=plots.multi(1,2,figsize=(6,4),hspace=0.001)
        plots.plotc(ax[0],allstar[1].data['LOGG'][m67],allstar[1].data['X_M'][m67,iel],allstar[1].data['TEFF'][m67],
                    yerr=allstar[1].data['X_M_ERR'][m67,iel],
                    xr=[6,0],yr=[-0.75,0.75],xt='log g',yt='['+el+'/M]')
        plots.plotc(ax[1],allstar[1].data['LOGG'][m67],allstar[1].data['X_H'][m67,iel],allstar[1].data['TEFF'][m67],
                    yerr=allstar[1].data['X_H_ERR'][m67,iel],
                    xr=[6,0],yr=[-0.75,0.75],xt='log g',yt='['+el+'/H]')
        ax[0].text(0.1,0.9,'rms: {:8.3f}'.format(allstar[1].data['X_M'][m67,iel].std()),transform=ax[0].transAxes)
        ax[1].text(0.1,0.9,'rms: {:8.3f}'.format(allstar[1].data['X_H'][m67,iel].std()),transform=ax[1].transAxes)
        outfile=out+'m67_{:s}.png'.format(el.strip())
        fig.savefig(outfile)
        plt.close(fig)
        grid.append([os.path.basename(outfile)])
    html.htmltab(grid,file=out+'m67.html',ytitle=els)

def calib(allstar,out='./') :
    """ Plot calibration relations
    """

    grid=[]
    yt=[]
    # Teff f([M/H], Teff)
    fig,ax=plots.multi(1,1)
    plots.plotc(ax,allstar[1].data['FPARAM'][:,3],allstar[1].data['PARAM'][:,0]-allstar[1].data['FPARAM'][:,0],allstar[1].data['FPARAM'][:,0],
                xt='[M/H]',yt=r'$\Delta$ Teff',zr=[3000,8000],colorbar=True,zt='Teff')
    outfile1=out+'calib_Teff_all.png'
    fig.savefig(outfile1)
    ax.set_xlim(-3,1)
    ax.set_ylim(-500,500)
    outfile2=out+'calib_Teff.png'
    fig.savefig(outfile2)
    plt.close(fig)
    grid.append([os.path.basename(outfile2),os.path.basename(outfile1)])
    yt.append('Teff')

    # logg f(logg,[M/H])
    fig,ax=plots.multi(1,1)
    plots.plotc(ax,allstar[1].data['FPARAM'][:,1],allstar[1].data['PARAM'][:,1]-allstar[1].data['FPARAM'][:,1],allstar[1].data['FPARAM'][:,3],
                xt='log g',yt=r'$\Delta$logg',zr=[-2,0.5],colorbar=True,zt='[M/H]')
    outfile1=out+'calib_logg_all.png'
    fig.savefig(outfile1)
    ax.set_xlim(-1,6)
    ax.set_ylim(-1,1)
    outfile2=out+'calib_logg.png'
    fig.savefig(outfile2)
    plt.close(fig)
    grid.append([os.path.basename(outfile2),os.path.basename(outfile1)])
    yt.append('log g')

    els=allstar[3].data['ELEM_SYMBOL'][0]
    for iel,el in enumerate(els) :
        print(el)
        fig,ax=plots.multi(1,1)
        if allstar[3].data['ELEMTOH'][0][iel] == 0 : abun = allstar[1].data['FELEM'][:,iel] 
        else : abun=allstar[1].data['FELEM'][:,iel]-allstar[1].data['FPARAM'][:,3]
        plots.plotc(ax,allstar[1].data['FPARAM'][:,0],allstar[1].data['X_M'][:,iel]-abun,allstar[1].data['FPARAM'][:,1],
                    xt='Teff',yt=r'$\Delta$['+el+'/M]',zr=[0,5.],colorbar=True,zt='log g')
        outfile1=out+'calib_{:s}_all.png'.format(el.strip())
        fig.savefig(outfile1)
        ax.set_xlim(3000,8000)
        ax.set_ylim(-0.5,0.5)
        outfile2=out+'calib_{:s}.png'.format(el.strip())
        fig.savefig(outfile2)
        plt.close(fig)
        grid.append([os.path.basename(outfile2),os.path.basename(outfile1)])
        yt.append(el.strip())

    html.htmltab(grid,file=out+'calib.html',ytitle=yt)

def flags(hdulist,out='./',alpha=0.005) :
    """ Tabulate number of objects with different flag bits set, and make HR diagrams showing these
    """
    f=html.head(out+'flags.html')
    mask=bitmask.AspcapBitMask()
    xt=[]
    for i in range(32) : 
        if mask.name[i] != '' : xt.append(mask.name[i])
    data=[]
    row=[]
    yt=['ASPCAPFLAG']
    for i in range(32) :
        j=np.where(hdulist[1].data['ASPCAPFLAG'] & 2**i)[0]
        if mask.name[i] == '' and len(j) > 0 :
            print('Unnamed bits are set!',i,len(j))
            pdb.set_trace()
        elif mask.name[i] != '' :
            print(mask.name[i])
            if len(j) > 0 :
                if np.core.defchararray.find(mask.name[i],'COLORTE') >= 0:
                    jk0=hdulist[1].data['J']-hdulist[1].data['K']-1.5*hdulist[1].data['AK_TARG']
                    fig,ax=plots.multi(1,2,hspace=0.001)
                    plots.plotc(ax[0],hdulist[1].data['FPARAM'][:,0],hdulist[1].data['FPARAM'][:,1],hdulist[1].data['FPARAM'][:,3],alpha=alpha,zr=[-2,0.5])
                    plots.plotc(ax[1],hdulist[1].data['FPARAM'][:,0],jk0,hdulist[1].data['FPARAM'][:,3],zr=[-2,0.5],alpha=alpha)
                    plots.plotc(ax[0],hdulist[1].data['FPARAM'][j,0],hdulist[1].data['FPARAM'][j,1],hdulist[1].data['FPARAM'][j,3],
                                xr=[10000,3000],yr=[6,-1],zr=[-2,0.5],xt='Teff (raw)',yt='logg (raw)')
                    plots.plotc(ax[1],hdulist[1].data['FPARAM'][j,0],jk0[j],hdulist[1].data['FPARAM'][j,3],
                                xr=[10000,3000],yr=[-2,4],zr=[-2,0.5],xt='Teff (raw)',yt='J-K')
                else :
                    fig,ax=plots.multi(1,1)
                    plots.plotc(ax,hdulist[1].data['FPARAM'][:,0],hdulist[1].data['FPARAM'][:,1],hdulist[1].data['FPARAM'][:,3],alpha=alpha,zr=[-2,0.5])
                    plots.plotc(ax,hdulist[1].data['FPARAM'][j,0],hdulist[1].data['FPARAM'][j,1],hdulist[1].data['FPARAM'][j,3],
                               xr=[10000,3000],yr=[6,-1],zr=[-2,0.5],xt='Teff (raw)',yt='logg (raw)')
                outfile=out+'flag_aspcapflag_{:d}.png'.format(i)
                fig.savefig(outfile)
                plt.close()
                row.append('<a href={:s}> {:d} </a>'.format(os.path.basename(outfile),len(j)))
            else : 
                row.append('{:d}'.format(len(j)))
    data.append(row)
    f.write(html.table(data,xtitle=xt,ytitle=yt,plots=False,formstr=':s'))

    mask=bitmask.ParamBitMask()
    xt=[]
    for i in range(32) : 
        if mask.name[i] != '' : xt.append(mask.name[i])
    data=[]
    yt=[]
    for iparam,param in enumerate(hdulist[3].data['PARAM_SYMBOL'][0]) :
        yt.append(param)
        row=[]
        for i in range(32) :
            print(param,mask.name[i])
            j=np.where(hdulist[1].data['PARAMFLAG'][:,iparam] & 2**i)[0]
            if mask.name[i] == '' and len(j) > 0 :
                print('Unnamed bits are set!',i,len(j))
                pdb.set_trace()
            elif mask.name[i] != '' :
                if len(j) > 0 :
                    fig,ax=plots.multi(1,1)
                    plots.plotc(ax,hdulist[1].data['FPARAM'][:,0],hdulist[1].data['FPARAM'][:,1],hdulist[1].data['FPARAM'][:,3],alpha=alpha,zr=[-2,0.5])
                    plots.plotc(ax,hdulist[1].data['FPARAM'][j,0],hdulist[1].data['FPARAM'][j,1],hdulist[1].data['FPARAM'][j,3],
                                xr=[10000,3000],yr=[6,-1],zr=[-2,0.5],xt='Teff (raw)',yt='logg (raw)')
                    outfile=out+'flag_param_{:d}_{:d}.png'.format(iparam,i)
                    fig.savefig(outfile)
                    plt.close()
                    row.append('<a href={:s}> {:d} </a>'.format(os.path.basename(outfile),len(j)))
                else : 
                    row.append('{:d}'.format(len(j)))
        data.append(row)
    f.write(html.table(data,xtitle=xt,ytitle=yt,plots=False,formstr=':s'))

    xt=[]
    for i in range(32) : 
        if mask.name[i] != '' : xt.append(mask.name[i])
    data=[]
    yt=[]
    for ielem,el in enumerate(hdulist[3].data['ELEM_SYMBOL'][0]) :
        if hdulist[3].data['ELEMTOH'][0][ielem] == 1 : fzr=[-2,0.5]
        else : fzr = [-1,1]
        yt.append(el)
        row=[]
        for i in range(32) :
            j=np.where(hdulist[1].data['ELEMFLAG'][:,ielem] & 2**i)[0]
            if mask.name[i] == '' and len(j) > 0 :
                print('Unnamed bits are set!',i,len(j))
                pdb.set_trace()
            elif mask.name[i] != '' :
                print(el,mask.name[i])
                if len(j) > 0 :
                    fig,ax=plots.multi(1,2,hspace=0.001)
                    plots.plotc(ax[0],hdulist[1].data['FPARAM'][:,0],hdulist[1].data['FPARAM'][:,1],hdulist[1].data['FPARAM'][:,3],alpha=alpha,zr=[-2,0.5])
                    plots.plotc(ax[1],hdulist[1].data['FPARAM'][:,0],hdulist[1].data['FPARAM'][:,1],hdulist[1].data['FELEM'][:,ielem],alpha=alpha,zr=fzr)
                    plots.plotc(ax[0],hdulist[1].data['FPARAM'][j,0],hdulist[1].data['FPARAM'][j,1],hdulist[1].data['FPARAM'][j,3],
                                xr=[10000,3000],yr=[6,-1],zr=[-2,0.5],xt='Teff (raw)',yt='logg (raw)',colorbar=True,zt='[M/H]')
                    plots.plotc(ax[1],hdulist[1].data['FPARAM'][j,0],hdulist[1].data['FPARAM'][j,1],hdulist[1].data['FELEM'][j,ielem],
                                xr=[10000,3000],yr=[6,-1],zr=fzr,xt='Teff (raw)',yt='logg (raw)',colorbar=True,zt='FELEM')
                    outfile=out+'flag_{:s}_{:d}.png'.format(el,i)
                    fig.savefig(outfile)
                    plt.close()
                    row.append('<a href={:s}> {:d} </a>'.format(os.path.basename(outfile),len(j)))
                else : 
                    row.append('{:d}'.format(len(j)))
        data.append(row)
    f.write(html.table(data,xtitle=xt,ytitle=yt,plots=False,formstr=':s'))

    html.tail(f)

def apolco(hdulist,out='./',snmin=150) :
    """ histograms of LCO-APO  parameters and abundances
    """

    gd=apselect.select(hdulist[1].data,badval='STAR_BAD',sn=[snmin,10000])
    a=hdulist[1].data[gd]

    apo=np.where(a['TELESCOPE'] == 'apo25m')[0]
    lco=np.where(a['TELESCOPE'] == 'lco25m')[0]
    i1,i2=match.match(a['APOGEE_ID'][apo],a['APOGEE_ID'][lco])
    grid=[]
    yt=[]
    for iparam,param in enumerate(hdulist[3].data['PARAM_SYMBOL'][0]) :
        fig,ax=plots.multi(1,1,figsize=(6,4.5))
        diff=a['FPARAM'][lco[i2],iparam]-a['FPARAM'][apo[i1],iparam]
        if iparam == 0 : ax.hist(diff,bins=np.arange(-100.,100.,1.))
        else : ax.hist(diff,bins=np.arange(-0.5,0.5,0.01))
        ax.set_xlabel('{:s} (LCO-APO)'.format(param))
        ax.text(0.1,0.9,'S/N> {:d}'.format(snmin),transform=ax.transAxes)
        ax.text(0.1,0.8,'mean: {:8.3f}'.format(diff.mean()),transform=ax.transAxes)
        ax.text(0.1,0.7,'std: {:8.3f}'.format(diff.std()),transform=ax.transAxes)
        outfile=out+'apolco_param_{:d}.png'.format(iparam)
        fig.savefig(outfile)
        plt.close()
        grid.append([os.path.basename(outfile)])
        yt.append(param)
    for ielem,el in enumerate(hdulist[3].data['ELEM_SYMBOL'][0]) :
        fig,ax=plots.multi(1,1,figsize=(6,4.5))
        diff=a['FELEM'][lco[i2],ielem]-a['FELEM'][apo[i1],ielem]
        ax.hist(diff,bins=np.arange(-0.5,0.5,0.01))
        ax.set_xlabel('{:s} (LCO-APO)'.format(el))
        ax.text(0.1,0.9,'S/N> {:d}'.format(snmin),transform=ax.transAxes)
        ax.text(0.1,0.8,'mean: {:8.3f}'.format(diff.mean()),transform=ax.transAxes)
        ax.text(0.1,0.7,'std: {:8.3f}'.format(diff.std()),transform=ax.transAxes)
        outfile=out+'apolco_{:s}.png'.format(el)
        fig.savefig(outfile)
        plt.close()
        grid.append([os.path.basename(outfile)])
        yt.append(el)

    html.htmltab(grid,file=out+'apolco.html',ytitle=yt)


def lbd2xyz(l,b,d,R0=8.5) :
    ''' Angular coordinates + distance -> galactocentry x,y,z '''

    brad = b*np.pi/180.
    lrad = l*np.pi/180.

    x = d*np.sin(0.5*np.pi-brad)*np.cos(lrad)-R0
    y = d*np.sin(0.5*np.pi-brad)*np.sin(lrad)
    z = d*np.cos(0.5*np.pi-brad)
    r = np.sqrt(x**2+y**2)
    return x, y, z, r

def elemsens(els=None,plot=None,ylim=[0.1,-0.3],teff=4750,logg=2.,feh=-1.,smooth=None) :
    '''
    Returns and optionally plots wavelength sensitivity to changes in elemental abundances for specified elements from MOOG mini-elem grid
    '''
    elem=fits.open(os.environ['APOGEE_REDUX']+'/speclib/moog/elemsens.fits')
    if els is None :
        els = elems()[0]
    elif type(els) == str :
        els = [els]
    wave=[]
    out=[]
    for el in els :
        for i in range(1,25) :
            card='HDU{:02d}'.format(i)
            try :
              if elem[0].header[card].strip().upper() == el.strip().upper() :
                it=int(round((teff-elem[i].header['CRVAL2'])/elem[i].header['CDELT2']))
                ig=int(round((logg-elem[i].header['CRVAL3'])/elem[i].header['CDELT3']))
                ife=int(round((feh-elem[i].header['CRVAL4'])/elem[i].header['CDELT4']))
                diff=elem[i].data[ife,ig,it,:]
                if smooth is not None:
                    diff=scipy.ndimage.filters.gaussian_filter(diff,smooth)
                wave=elem[i].header['CRVAL1']+np.arange(elem[i].header['NAXIS1'])*elem[i].header['CDELT1']
                if plot is not None:
                    #plot.plot(wave,diff,color='g')
                    plot.plot(wave,diff)
                    plot.set_ylim(ylim[0],ylim[1])
                out.append(diff)
            except: pass
    if len(out) == 1 :
        return wave, out[0]
    else :
        return wave, out

def sensplot(ax=None,offset=0) :
    if ax is None :
        fig,ax=plots.multi(1,2,hspace=0.001,sharex=True)
    els=['O','Mg','Si','S','Ca','Ti','Na','Al','K','P']
    cols=['r','g','b','c','y','m','r','g','b','c']
    ls=['-','-','-','-','-','-',':',':',':',':']
    for i in range(len(els)) :
        w,s=elemsens(els=els[i])
        plots.plotl(ax[0],w,s+offset,label=els[i],color=cols[i],ls=ls[i])
    ax[0].legend(fontsize='small')

    elems=['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Ce','Rb','Y','Nd']
    els=['V','Cr','Mn','Co','Ni','Cu','Ge','Ce','Rb','Nd']
    cols=['r','g','b','c','y','m','r','g','b','c']
    ls=['-','-','-','-','-','-',':',':',':',':']
    for i in range(len(els)) :
        w,s=elemsens(els=els[i])
        plots.plotl(ax[1],w,s+offset,label=els[i],color=cols[i],ls=ls[i])
    ax[1].legend(fontsize='small')
    #elems=['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Ce','Rb','Y','Nd']
    

def intplot(a=None,param='FPARAM',indir='cal',apred='r10',aspcap='t33b',verbose=False) :
    """ Given input structure, plot HR diagram, and enter event loop to mark stars to plot spectra
    """

    load=apload.ApLoad(apred=apred,aspcap=aspcap,verbose=verbose)
    if a is None : a=load.allCal()[1].data

    fig,ax=hr(a)
    plots.event(fig)
    sf,sa=plots.multi(1,1)
    nplot = 11
    hf,ha=plots.multi(1,nplot,figsize=(8.5,11),hspace=0.2)
    ha2=[]
    for i in range(nplot) : ha2.append(ha[i].twinx())
    print('hit any key near object to plot in HR diagram, q to quit')
    while (1) :
        ret=plots.mark(fig)
        if ret[2] == 'q' : break
        ind=plots._index[0]
        load.settelescope('apo25m')
        try : f=load.aspcapField(a['ALTFIELD'][ind])
        except : f=load.aspcapField(a['FIELD'][ind])
        if f is None :
            load.settelescope('lco25m')
            try : f=load.aspcapField(a['ALTFIELD'][ind])
            except : f=load.aspcapField(a['FIELD'][ind])
        if f is None :
            f=glob.glob(indir+'/*'+a['APOGEE_ID'][plots._index[0]]+'*')
            dir=os.path.dirname(f[0])
            f=glob.glob(dir+'/*aspcapField*.fits')
            f=fits.open(f[0])
        print('f: ',f)
        data=f[1].data
        j=np.where(data['APOGEE_ID'] == a['APOGEE_ID'][plots._index[0]])[0][0]
        sa.cla()
        for i in range(11) : 
            ha[i].cla()
            ha2[i].cla()
            ha[i].set_ylabel('Flux')
            ha2[i].set_ylabel(r'$\chi^2$')
            ha[i].set_ylim(0.5,1.3)
            ha2[i].set_ylim(0.,20.)
        plot(10.**f[3].data['WAVE'][0],f[2].data['SPEC'][j,:],ax=ha,sum=True,color='k')
        plot(10.**f[3].data['WAVE'][0],f[2].data['SPEC_BESTFIT'][j,:],ax=ha,sum=True,color='b')
        chi2 = (f[2].data['SPEC'][j,:]-f[2].data['SPEC_BESTFIT'][j,:])**2/f[2].data['ERR'][j,:]**2
        plot(10.**f[3].data['WAVE'][0],chi2,ax=ha2,sum=True,alpha=0.4)
        plots.plotl(sa,10.**f[3].data['WAVE'][0],f[2].data['SPEC'][j,:])
        plots.plotl(sa,10.**f[3].data['WAVE'][0],f[2].data['SPEC_BESTFIT'][j,:])
        text1=r'ID: {:s} FIELD: {:s} SNR: {:6.1f} $\chi^2$: {:6.1f}'.format(
             data['APOGEE_ID'][j],data['FIELD'][j],data['SNR'][j],data['PARAM_CHI2'][j])
        text2=r'Teff: {:5.0f} logg: {:5.1f} [M/H]: {:5.2f} [$\alpha$/M]: {:5.2f} [C/M]: {:5.2f} [N/M]: {:5.2f}'.format(
             data[param][j,0],data[param][j,1],data[param][j,3],data[param][j,6],data[param][j,4],data[param][j,5])
        sf.suptitle(text1+'\n'+text2)
        hf.suptitle(text1+'\n'+text2)
        plt.draw()
        plt.show()
    plt.close(hf)
    plt.close(sf)
    plt.close(fig)

def hr(a,param='FPARAM',colorbar=False,zt='[M/H]',zr=None,iso=None, alpha=0.3,hard=None, gridclass=None,xr=[8000,3000],yr=[6,-1],grid=False,contour=False,snrbd=0,target=None,size=5) :
    """ Plot an HR diagram from input structure

        Args:
            all  : structure that includes stellar parameter array with (Teff, logg, ...)
            param : tag to use (default='FPARAM')
            colorbar : show colorbar? (default= False)
    """
    fig,ax = plots.multi(1,2,figsize=(8,12),hspace=0.001)
    if gridclass is None :
        teff=a[param][:,0]
        logg=a[param][:,1]
    else :
        teff=a['FPARAM_CLASS'][:,gridclass,0]
        logg=a['FPARAM_CLASS'][:,gridclass,1]
    if zt == '[M/H]' : 
        z=a[param][:,3]
        if zr is None : zr=[-2,0.5]
    elif zt == 'chi2' : 
        z=a['PARAM_CHI2']
        if zr is None : zr=[0,10]
    aspcapmask=bitmask.AspcapBitMask()
    starmask=bitmask.StarBitMask()
    bd=np.where( ((a['ASPCAPFLAG']&aspcapmask.badval()) > 0) |
                 ((a['STARFLAG']&starmask.badval()) > 0) |
                  (a['SNR']<snrbd) ) [0]
    gd=np.where( ((a['ASPCAPFLAG']&aspcapmask.badval()) == 0) &
                 ((a['STARFLAG']&starmask.badval()) == 0) &
                  (a['SNR']>=snrbd) ) [0]
    if contour :
        plots.plotp(ax[0],teff,logg,xr=xr,yr=yr,
                    xt='Teff',yt='log g',contour=-1)
        plots.plotp(ax[1],teff,logg,xr=xr,yr=yr,
                    xt='Teff',yt='log g',contour=-1)
    else :
        plots.plotc(ax[0],teff[gd],logg[gd],z[gd],xr=xr,yr=yr,zr=zr,
                    xt='Teff',yt='log g',zt=zt,colorbar=colorbar,size=size)
        plots.plotc(ax[1],teff[gd],logg[gd],z[gd],xr=xr,yr=yr,zr=zr,
                    xt='Teff',yt='log g',zt=zt,colorbar=colorbar,size=size)
        plots.plotp(ax[1],teff[bd],logg[bd],color='k',size=2)
    if grid: 
        ax[0].grid()
        ax[1].grid()
    ax[0].text(0.05,0.9,'{:d} stars, S/N>{:5.0f}, no GRIDEDGE_BAD'.format(len(gd),snrbd),transform=ax[0].transAxes)
    ax[1].text(0.05,0.9,'{:d} stars'.format(len(gd)+len(bd)),transform=ax[1].transAxes)
    plots._data = a
    if iso is not None:
        cmap = matplotlib.cm.get_cmap('rainbow')
        for mh in [-2.0,-1.0,0.0,0.5] :
            if mh < -0.01 : name = 'zm{:02d}'.format(int(abs(mh)*10.))
            else : name = 'zp{:02d}'.format(int(abs(mh)*10.))
            rgba=cmap((mh-zr[0])/(zr[1]-zr[0]))
            for age in iso :
                isodata=isochrones.read(os.environ['ISOCHRONE_DIR']+'/'+name+'.dat',agerange=[age-0.01,age+0.01])
                isochrones.plot(ax[0],isodata,'teff','logg',color=rgba,alpha=alpha)
    if hard is not None: 
        fig.savefig(hard)
        plt.close()

    if (target is not None) and ('EXTRATARG' in a.columns.names) :
        tfig,tax=plots.multi(1,1)
        main =np.where(a['EXTRATARG'][gd] == 0)[0]
        plots.plotc(tax,teff[gd[main]],logg[gd[main]],z[gd[main]],xr=xr,yr=yr,zr=zr,
                    xt='Teff',yt='log g',zt=zt,colorbar=colorbar,size=size)
        tax.text(0.05,0.9,'{:d} stars, main sample, S/N>{:5.0f}, no GRIDEDGE_BAD'.format(len(main),snrbd),transform=tax.transAxes)
        tax.grid()

        if iso is not None:
            cmap = matplotlib.cm.get_cmap('rainbow')
            for mh in [-2.0,-1.0,0.0,0.5] :
                if mh < -0.01 : name = 'zm{:02d}'.format(int(abs(mh)*10.))
                else : name = 'zp{:02d}'.format(int(abs(mh)*10.))
                rgba=cmap((mh-zr[0])/(zr[1]-zr[0]))
                for age in iso :
                    isodata=isochrones.read(os.environ['ISOCHRONE_DIR']+'/'+name+'.dat',agerange=[age-0.01,age+0.01])
                    isochrones.plot(tax,isodata,'teff','logg',color=rgba,alpha=alpha)
        tfig.savefig(target+'_main.png')
        plt.close(tfig)

        t1=bitmask.Apogee2Target1()
        t2=bitmask.Apogee2Target2()
        t3=bitmask.Apogee2Target3()
        t1_1=bitmask.ApogeeTarget1()
        t2_1=bitmask.ApogeeTarget2()
        tfig,tax=plots.multi(1,1)
        plots.plotp(tax,teff[gd[main]],logg[gd[main]],color='b',xr=xr,yr=yr,
                    xt='Teff',yt='log g',label='MAIN')

        mc =np.where(((a['APOGEE2_TARGET1'][gd] & t1.getval(['APOGEE2_MAGCLOUD_MEMBER'])) > 0) |
                     ((a['APOGEE2_TARGET1'][gd] & t1.getval(['APOGEE2_MAGCLOUD_CANDIDATE'])) > 0) )[0]
        plots.plotp(tax,teff[gd[mc]],logg[gd[mc]],color='g',xr=xr,yr=yr,
                    xt='Teff',yt='log g',label='MC')
        dsph =np.where(((a['APOGEE2_TARGET1'][gd] & t1.getval(['APOGEE2_DSPH_MEMBER'])) > 0) |
                     ((a['APOGEE2_TARGET1'][gd] & t1.getval(['APOGEE2_DSPH_CANDIDATE'])) > 0) |
                     ((a['APOGEE2_TARGET1'][gd] & t1.getval(['APOGEE2_SGR_DSPH'])) > 0)  )[0]
        plots.plotp(tax,teff[gd[dsph]],logg[gd[dsph]],color='c',xr=xr,yr=yr,
                    xt='Teff',yt='log g',label='DSPH')
        rr =np.where((a['APOGEE2_TARGET1'][gd] & t1.getval(['APOGEE2_RRLYR'])) > 0 )[0]
        plots.plotp(tax,teff[gd[rr]],logg[gd[rr]],color='r',xr=xr,yr=yr,
                    xt='Teff',yt='log g',label='RRLYR')
        young =np.where(((a['APOGEE2_TARGET3'][gd] & t3.getval(['APOGEE2_YOUNG_CLUSTER'])) > 0) |
                        ((a['APOGEE_TARGET2'][gd] & t2_1.getval(['APOGEE_EMBEDDEDCLUSTER_STAR'])) > 0) )[0]
        plots.plotp(tax,teff[gd[young]],logg[gd[young]],color='m',xr=xr,yr=yr,
                    xt='Teff',yt='log g',label='YOUNG')
        emission =np.where( ((a['APOGEE_TARGET2'][gd] & t2_1.getval(['APOGEE_EMISSION_STAR'])) > 0) )[0]
        plots.plotp(tax,teff[gd[emission]],logg[gd[emission]],color='y',xr=xr,yr=yr,
                    xt='Teff',yt='log g',label='EMISSION')
        extended =np.where( ((a['APOGEE_TARGET1'][gd] & t1_1.getval(['APOGEE_EXTENDED'])) > 0) |
                     ((a['APOGEE_TARGET1'][gd] & t1_1.getval(['APOGEE_M31_CLUSTER'])) > 0) |
                     ((a['APOGEE2_TARGET3'][gd] & t3.getval(['APOGEE2_M31'])) > 0) |
                     ((a['APOGEE2_TARGET3'][gd] & t3.getval(['APOGEE2_M33'])) > 0) )[0]
        plots.plotp(tax,teff[gd[extended]],logg[gd[extended]],color='orange',xr=xr,yr=yr,
                    xt='Teff',yt='log g',label='EXTENDED')
        tax.grid()
        tax.legend(loc='upper left')
        tfig.savefig(target+'_targ.png')
        plt.close(tfig)

    return fig,ax

def multihr(a,param='FPARAM',colorbar=False,hard=None,xr=[8000,3000],yr=[6,-1],size=5) :
    """ Series of HR diagram plots, color-coded by different quantities
    """
    fig,ax = plots.multi(3,4,hspace=0.001,wspace=0.001,figsize=(12,10))

    aspcapmask=bitmask.AspcapBitMask()
    starmask=bitmask.StarBitMask()
    bd=np.where( ((a['ASPCAPFLAG']&aspcapmask.badval()) > 0) |
                 ((a['STARFLAG']&starmask.badval()) > 0) ) [0]
    gd=np.where( ((a['ASPCAPFLAG']&aspcapmask.badval()) == 0) &
                 ((a['STARFLAG']&starmask.badval()) == 0) ) [0]

    z=a[param][gd,3]
    zr=[-2,0.5]
    zt='[M/H] (-2:0.5)'
    plots.plotc(ax[0,0],a[param][gd,0],a[param][gd,1],z,xr=xr,yr=yr,zr=zr,
                xt='Teff',yt='log g',zt=zt,colorbar=colorbar,size=size)
    ax[0,0].text(0.05,0.9,zt,transform=ax[0,0].transAxes)

    z=10.**a[param][gd,2]
    zr=[0.3,4]
    zt='vmicro (0.3:4)'
    plots.plotc(ax[0,1],a[param][gd,0],a[param][gd,1],z,xr=xr,yr=yr,zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar,size=size)
    ax[0,1].text(0.05,0.9,zt,transform=ax[0,1].transAxes)

    z=10.**a[param][gd,7]
    zr=[0,10]
    zt='vrot (0:10)'
    plots.plotc(ax[0,2],a[param][gd,0],a[param][gd,1],z,xr=xr,yr=yr,zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar,size=size)
    ax[0,2].text(0.05,0.9,zt,transform=ax[0,2].transAxes)

    z=a[param][gd,4]
    zr=[-0.5,0.5]
    zt='[C/M] (-0.5:0.5)'
    plots.plotc(ax[1,0],a[param][gd,0],a[param][gd,1],z,xr=xr,yr=yr,zr=zr,
                xt='Teff',yt='log g',zt=zt,colorbar=colorbar,size=size)
    ax[1,0].text(0.05,0.9,zt,transform=ax[1,0].transAxes)

    z=a[param][gd,5]
    zr=[-0.5,0.5]
    zt='[N/M] (-0.5:0.5)'
    plots.plotc(ax[1,1],a[param][gd,0],a[param][gd,1],z,xr=xr,yr=yr,zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar,size=size)
    ax[1,1].text(0.05,0.9,zt,transform=ax[1,1].transAxes)

    z=a[param][gd,4]-a[param][gd,5]
    zr=[-0.5,0.5]
    zt='[C/N] (-0.5:0.5)'
    plots.plotc(ax[1,2],a[param][gd,0],a[param][gd,1],z,xr=xr,yr=yr,zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar,size=size)
    ax[1,2].text(0.05,0.9,zt,transform=ax[1,2].transAxes)

    z=a[param][gd,6]
    zr=[-1,1.0]
    zt=r'[$\alpha$/M] (-1:1)'
    plots.plotc(ax[2,0],a[param][gd,0],a[param][gd,1],z,xr=xr,yr=yr,zr=zr,
                xt='Teff',yt='log g',zt=zt,colorbar=colorbar,size=size)
    ax[2,0].text(0.05,0.9,zt,transform=ax[2,0].transAxes)

    z=a['VSCATTER']
    zr=[0,5]
    zt='VSCATTER'
    plots.plotc(ax[2,1],a[param][gd,0],a[param][gd,1],z[gd],xr=xr,yr=yr,zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar,size=size)
    ax[2,1].text(0.05,0.9,zt,transform=ax[2,1].transAxes)

    z=a['MEANFIB']
    zr=[0,300]
    zt='MEANFIB'
    plots.plotc(ax[2,2],a[param][gd,0],a[param][gd,1],z[gd],xr=xr,yr=yr,zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar,size=size)
    ax[2,2].text(0.05,0.9,zt,transform=ax[2,2].transAxes)

    try: z=np.log10(a['PARAM_CHI2'])
    except: z=np.log10(a['ASPCAP_CHI2'])
    zr=[0,2]
    zt='log(CHI2) (0:2)'
    plots.plotc(ax[3,0],a[param][gd,0],a[param][gd,1],z[gd],xr=xr,yr=yr,zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar,size=size)
    ax[3,0].text(0.05,0.9,zt,transform=ax[3,0].transAxes)

    z=a['SNR']
    zr=[20,200]
    zt='SNR (20:200)'
    plots.plotc(ax[3,1],a[param][gd,0],a[param][gd,1],z[gd],xr=xr,yr=yr,zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar,size=size)
    ax[3,1].text(0.05,0.9,zt,transform=ax[3,1].transAxes)

    z=a['SNR']
    zr=[20,100]
    zt='SNR (20:100)'
    plots.plotc(ax[3,2],a[param][gd,0],a[param][gd,1],z[gd],xr=xr,yr=yr,zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar,size=size)
    ax[3,2].text(0.05,0.9,zt,transform=ax[3,2].transAxes)

    if hard is not None: 
        fig.savefig(hard)
        plt.close()

def plot(wave,spec,color=None,figax=None,ax=None,hard=None,sum=False,title=None,alpha=None,yr=None,lineids=None,multipage=False, refline=None, figsize=(8,11), textsize=8) :
    """  Multipanel plots of APOGEE spectra
    """
    # set up plots
    if sum : ny=11
    else : ny=10
    if figax is not None : fig,ax=figax
    if ax is None : fig,ax=plots.multi(1,ny,figsize=figsize,hspace=0.2)

    # get line labels if requested
    if lineids is not None :
        file = os.environ['APOGEE_DIR']+'/data/lines/atlas_line_ids_apogee.txt'
        lines = np.loadtxt(file,delimiter=';',dtype={'names' : ('wave','label'), 'formats': ('f4','S24') } )
        lines['wave'] = spectra.airtovac(lines['wave'])

    # plot chunks of 200 A
    for i in range(10) :
        plots.plotl(ax[i],wave,spec,xr=[15000+i*200,15200+i*200],color=color,linewidth=0.3,alpha=alpha,yr=yr)
        ax[i].xaxis.label.set_size(6)
        ax[i].yaxis.label.set_size(6)
        ax[i].tick_params(axis = 'both', which = 'major', labelsize = 6)
        ax[i].xaxis.set_minor_locator(plt.MultipleLocator(10.))
        if lineids is not None :
            gd = np.where( (lines['wave'] > 15000+i*200) & ( lines['wave'] < 15200+i*200) ) [0]
            for line in lines[gd] :
                ax[i].text(line['wave'],lineids,line['label'],rotation=90,size=textsize,ha='center',va='bottom')
        if refline is not None :
            # plot a reference horizontal line
            plots.plotl(ax[i],wave,spec*0.+refline,xr=[15000+i*200,15200+i*200],color=color,linewidth=0.3,alpha=alpha,yr=yr,ls=':')

        if multipage :
            mfig,multiax = plots.multi(1,1)
            plots.plotl(multiax,wave,spec,xr=[15000+i*200,15200+i*200],color=color,linewidth=0.3,alpha=alpha,yr=yr)
            multiax.xaxis.label.set_size(6)
            multiax.yaxis.label.set_size(6)
            multiax.tick_params(axis = 'both', which = 'major', labelsize = 6)
            multiax.xaxis.set_minor_locator(plt.MultipleLocator(10.))
            if lineids is not None :
                gd = np.where( (lines['wave'] > 15000+i*200) & ( lines['wave'] < 15200+i*200) ) [0]
                for line in lines[gd] :
                    multiax.text(line['wave'],1.,line['label'],rotation=90,size=4,ha='left',va='bottom')

    # final panel with full wavelength range if requested
    if sum :
        plots.plotl(ax[10],wave,spec,xr=[15100,17000],color=color,linewidth=0.3,alpha=alpha,yr=yr)
        ax[10].xaxis.label.set_size(6)
        ax[10].yaxis.label.set_size(6)
        ax[10].tick_params(axis = 'both', which = 'major', labelsize = 6)
        ax[10].xaxis.set_minor_locator(plt.MultipleLocator(100.))

    try: 
        if title is not None : fig.suptitle(title)
    except: pass
    if hard is not None : fig.savefig(hard)

    try: return fig,ax
    except: return

def multiwind(data,apred='r10',aspcap='t33w',out='plots/') :
    """ Plot results from different windows for each element
    """
    #load=apload.ApLoad(apred=apred,aspcap=aspcap)
    #data=load.allCal()
    els=data[3].data['ELEM_SYMBOL'][0]
    elemtoh=data[3].data['ELEMTOH'][0]
    grid=[]
    ytit=[]
    for iel,el in enumerate(els) :
        if os.path.exists(os.environ['APOGEE_ASPCAP']+'/'+apred+'/'+aspcap+'/config/apogee-n/'+el+'.wind') :
            w=np.loadtxt(os.environ['APOGEE_ASPCAP']+'/'+apred+'/'+aspcap+'/config/apogee-n/'+el+'.wind') 
            nwind=w.shape[0]
            fig,ax=plots.multi(1,nwind,hspace=0.001,figsize=(6,nwind))
            fig.suptitle(el)
            for i in range(1,nwind+1) : 
                plots.plotc(ax[i-1],data[1].data['FPARAM'][:,0],data[1].data['FELEM'][:,i,iel]-data[1].data['FELEM'][:,0,iel],
                            data[1].data['FPARAM'][:,3],xr=[3000,6000],yr=[-1,1],zr=[-2,0.5],xt='Teff',size=5,yt=r'$\Delta$(line-global)')
                ax[i-1].text(0.05,0.8,'{:8.2f}-{:8.2f}   {:8.2f}'.format(w[i-1,0],w[i-1,1],w[i-1,2]),transform=ax[i-1].transAxes,fontsize=10)
                ax[i-1].yaxis.label.set_size(6)
            fig.savefig(out+el+'.png')
            fig,ax=plots.multi(1,nwind,hspace=0.001,figsize=(6,nwind))
            fig.suptitle(el)
            for i in range(1,nwind+1) : 
                if elemtoh[iel] == 1: y = data[1].data['FELEM'][:,i,iel]-data[1].data['FPARAM'][:,3]
                else : y = data[1].data['FELEM'][:,i,iel]
                plots.plotc(ax[i-1],data[1].data['FPARAM'][:,3],y,
                            data[1].data['FPARAM'][:,0],xr=[-2.5,1.],yr=[-0.5,0.5],zr=[3500,5500],xt='[M/H]',size=5,yt='[X/M]')
                ax[i-1].text(0.05,0.8,'{:8.2f}-{:8.2f}   {:8.2f}'.format(w[i-1,0],w[i-1,1],w[i-1,2]),transform=ax[i-1].transAxes,fontsize=10)
                ax[i-1].yaxis.label.set_size(6)
            fig.savefig(out+el+'_2.png')
            grid.append([el+'.png',el+'_2.png'])
            ytit.append(el)
            plt.close()
            plt.close()
    html.htmltab(grid,ytitle=ytit,file=out+'wind.html')

def compspec(a,b,j=0,hard=None) :
    """ compare spectra and best fits from two different input aspcapField files for specified object
    """
    nplot=11
    hf,ha=plots.multi(1,nplot,figsize=(33,44),hspace=0.2)
    ha2=[]
    for i in range(nplot) : ha2.append(ha[i].twinx())
    print(hf.get_size_inches())
    chi2_a = (a[2].data['SPEC'][j,:]-a[2].data['SPEC_BESTFIT'][j,:])**2/a[2].data['ERR'][j,:]**2
    chi2_b = (b[2].data['SPEC'][j,:]-b[2].data['SPEC_BESTFIT'][j,:])**2/b[2].data['ERR'][j,:]**2
    plot(10.**a[3].data['WAVE'][0],a[2].data['SPEC'][j,:] ,sum=True,color='k',figax=(hf,ha),yr=[0.8,1.3],lineids=1.1)
    plot(10.**a[3].data['WAVE'][0],a[2].data['SPEC_BESTFIT'][j,:] ,sum=True,color='r',figax=(hf,ha),yr=[0.8,1.3])
    plot(10.**b[3].data['WAVE'][0],b[2].data['SPEC_BESTFIT'][j,:] ,sum=True,color='b',figax=(hf,ha),yr=[0.8,1.3])
    plot(10.**a[3].data['WAVE'][0],a[2].data['SPEC'][j,:]-b[2].data['SPEC_BESTFIT'][j,:]+0.95 ,sum=True,color='c',
         figax=(hf,ha),yr=[0.8,1.3],refline=0.9)
    plot(10.**a[3].data['WAVE'][0],a[2].data['SPEC_BESTFIT'][j,:]-b[2].data['SPEC_BESTFIT'][j,:]+0.9 ,sum=True,color='m',
         figax=(hf,ha),yr=[0.8,1.3],refline=0.9)
    plot(10.**a[3].data['WAVE'][0],chi2_a-chi2_b ,sum=True,color='g',multipage=False,figax=(hf,ha2),yr=[-10,10],refline=0.)
    if hard is not None : hf.savefig(hard+'.pdf')

    fig,ax=plots.multi(1,2,hspace=0.001)
    plots.plotl(ax[0],10.**a[3].data['WAVE'][0],np.cumsum(chi2_a),xt='Wavelength',yt='Cumulative chi^2')
    plots.plotl(ax[0],10.**a[3].data['WAVE'][0],np.cumsum(chi2_b))
    plots.plotl(ax[1],10.**a[3].data['WAVE'][0],np.cumsum(chi2_a)-np.cumsum(chi2_b),yt='Difference in cumulative chi^2')
    if hard is not None : fig.savefig(hard+'_chi2.pdf')
    plt.close()
    plt.close()

def vesta() :
    """ series of plots for Vesta for DR16 tests
    """
    l33_fit = fits.open(os.environ['APOGEE_ASPCAP']+'/r10/test/apo1m/fit/aspcapField-standards.fits')
    l33_fixed = fits.open(os.environ['APOGEE_ASPCAP']+'/r10/test/apo1m/fixed/aspcapField-standards.fits')
    l31c_fit = fits.open(os.environ['APOGEE_ASPCAP']+'/r10/testl31c/apo1m/fit/aspcapField-standards.fits')
    l31c_fixed = fits.open(os.environ['APOGEE_ASPCAP']+'/r10/testl31c/apo1m/fixed/aspcapField-standards.fits')
    l31crenorm_fit = fits.open(os.environ['APOGEE_ASPCAP']+'/r10/testl31crenorm/apo1m/fit/aspcapField-standards.fits')
    l31crenorm_fixed = fits.open(os.environ['APOGEE_ASPCAP']+'/r10/testl31crenorm/apo1m/fixed/aspcapField-standards.fits')
    compspec(l33_fit,l33_fixed,hard='Vesta_l33_fitvsfixed')
    compspec(l31c_fit,l31c_fixed,hard='Vesta_l31c_fitvsfixed')
    #compspec(l31c_fixed,l33_fixed,hard='Vesta_fixed_l31cvl33')
    compspec(l31crenorm_fit,l31crenorm_fixed,hard='Vesta_l31crenorm_fitvsfixed')

def repeat(data,out=None) :
    """ Comparison of repeat observations of objects
    """

    a=data[1].data
    els=data[3].data['ELEM_SYMBOL'][0]
    stars = set(a['APOGEE_ID'])
    fig,ax=plots.multi(2,7,figsize=(8,18),wspace=0.4)
    efig,eax=plots.multi(2,len(els),hspace=0.001,figsize=(8,36),wspace=0.4)
    telescope=np.zeros(len(a),dtype='S6')
    colors=['r','g','b']
    tels=['apo1m','apo25m','lco25m'] 
    for tel in tels :
        j=np.where(np.core.defchararray.find(a['ALTFIELD'],tel) >= 0)[0]
        telescope[j] = tel
    diff=[]
    ediff=[]
    tdiff=[]
    teldiff=[]
    for star in stars :
        j = np.where(a['APOGEE_ID'] == star)[0]
        n = len(j)
        if n > 1 :
            print(star,n,telescope[j])
            for i in j : 
                ediff.append(a['FELEM'][i,0,:]-a['FELEM'][j,0,:].mean(axis=0))
                diff.append(a['FPARAM'][i,:]-a['FPARAM'][j,:].mean(axis=0))
                tdiff.append(a['FPARAM'][j,0].mean())
                teldiff.append(telescope[i])
            for i in range(7)  :
                plots.plotp(ax[i,0],np.repeat(a['FPARAM'][j,0].mean(),n),a['FPARAM'][j,i]-a['FPARAM'][j,i].mean(),typeref=telescope[j],
                            types=tels,color=colors)
                #for itel,tel in enumerate(tels) :
                #    gd=np.where(np.core.defchararray.find(a['ALTFIELD'][j],tel) >= 0)[0]
                #    gd=j[gd]
            for i in range(len(els))  :
                plots.plotp(eax[i,0],np.repeat(a['FPARAM'][j,0].mean(),n),a['FELEM'][j,0,i]-a['FELEM'][j,0,i].mean(),typeref=telescope[j],
                            types=tels,color=colors,yr=[-0.2,0.2])
                #for itel,tel in enumerate(tels) :
                #    gd=np.where(np.core.defchararray.find(a['ALTFIELD'][j],tel) >= 0)[0]
                #    gd=j[gd]
                #    ax[i,1].hist(a['FELEM'][gd,i]-a['FELEM'][gd,i].mean(),color=colors[itel])
    diff=np.array(diff) 
    ediff=np.array(ediff) 
    tdiff=np.array(tdiff) 
    teldiff=np.array(teldiff) 
    for itel,tel in enumerate(tels) :
        gd=np.where(teldiff == tel)[0]
        if len(gd) > 0 :
            for i in range(7)  :
                if i == 0 : bins=np.arange(-200,200,10)
                elif i ==1 : bins=np.arange(-0.5,0.5,0.025)
                else : bins=np.arange(-0.2,0.2,0.01)
                ax[i,1].hist(diff[gd,i],color=colors[itel],histtype='step',bins=bins)
            for i in range(len(els))  :
                bins=np.arange(-0.2,0.2,0.01)
                try: eax[i,1].hist(ediff[gd,i],color=colors[itel],histtype='step',bins=bins)
                except: pass
    for i,el in enumerate(els) : 
        eax[i,0].set_ylabel(el)
        eax[i,1].text(0.1,0.9,el,transform=eax[i,1].transAxes)
    for i,param in enumerate(data[3].data['PARAM_SYMBOL'][0]) : 
        if i < 7 :
            ax[i,0].set_ylabel(param)
            ax[i,1].text(0.1,0.9,param,transform=ax[i,1].transAxes)
    if out is not None :
        fig.savefig(out+'param_diff.png')
        efig.savefig(out+'elem_diff.png')
    else :
        pdb.set_trace()
    plt.close(fig)
    plt.close(efig)


def average(a,ind,apred='r12',aspcap='l33', median=False) :
    ''' Average together multiple aspcapStar spectra
    '''

    load=apload.ApLoad(apred=apred,aspcap=aspcap)

    spec=np.zeros([8575])
    err=np.zeros([8575])
    ratio=np.zeros([8575])
    spec=[]
    err=[]
    ratio=[]
    for j in ind :
        load.settelescope(a['TELESCOPE'][j])
        z=load.aspcapStar(a['FIELD'][j],a['APOGEE_ID'][j])
        spec.append(z[1].data)
        err.append(z[2].data)
        ratio.append(z[1].data/z[3].data)
    spec=np.array(spec)
    err=np.array(err)
    ratio=np.array(ratio)

    if median: return np.median(spec,axis=0), np.median(err,axis=0), np.median(ratio,axis=0)
    else : return spec.mean(axis=0) ,err.mean(axis=0), ratio.mean(axis=0)


def resid(data,spec,teff=[3000,3500,4000,4500,5000,5500],logg=[0,2,4,6],out='resid') :
    """ Make average residual plots for different ranges of Teff, logg
    """
    wave = np.hstack(aspcap.gridWave()).flatten()

    grid=[]
    for ite in range(len(teff[0:-1])) :
        xg=[]
        for ilogg in range(len(logg[0:-1])) :
            fig,ax=plots.multi(1,2,hspace=0.001,wspace=0.001,figsize=(24,8))
            gd = np.where((data['FPARAM'][:,0]>=teff[ite]) & (data['FPARAM'][:,0]<teff[ite+1])  &
                          (data['FPARAM'][:,1]>=logg[ilogg]) & (data['FPARAM'][:,1]<logg[ilogg+1])  ) [0]
            res=np.median( spec['SPEC'][gd,:]-spec['SPEC_BESTFIT'][gd,:],axis=0 )
            plots.plotl(ax[0],wave,res)
            med=np.median( spec['SPEC'][gd,:],axis=0 )
            plots.plotl(ax[1],wave,med)
            name='{:s}_{:d}_{:d}.png'.format(out,ite,ilogg)
            fig.savefig(name)
            plt.close()
            xg.append(name)
        grid.append(xg)

    html.htmltab(grid,file='resid.html')
