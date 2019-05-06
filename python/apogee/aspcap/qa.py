import numpy as np
import os
import pdb
from astropy.io import fits
from tools import plots
from tools import html
from tools import match
from tools import fit
import matplotlib.pyplot as plt
from apogee.aspcap import aspcap
from apogee.aspcap import err
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

def plotelems(hdulist,title=None,out=None,calib=False) :
    """ Make [X/M] vs [M/H] plots for all elements as f(Teff, logg)
    """
    a=hdulist[1].data
    if 'EXTRATARG' in a.columns.names :
        gd=np.where(a['EXTRATARG'] == 0)[0]
        a=a[gd]
        comment=', main sample only'
    else :
        comment=', full sample'

    if calib : param = 'PARAM'
    else : param = 'FPARAM'

    els=hdulist[3].data['ELEM_SYMBOL'][0]
    etoh=hdulist[3].data['ELEMTOH'][0]
    grid=[]
    yt=[]
    for iel,el in enumerate(els) :
        row=[]
        xt=[]
        yt.append(el)
        icol=0
        for te,logg in zip(te_ranges,logg_ranges) :
            gd = np.where((a[param][:,0] >= te[0]) & (a[param][:,0] <= te[1]) &
                          (a[param][:,1] >= logg[0]) & (a[param][:,1] <= logg[1]) )[0]
            print(el,te,logg,len(gd))
            if calib :
                abun=a['X_M'][gd,iel]
            else :
                try:
                    if etoh[iel] == 1 : abun=a['FELEM'][gd,0,iel]-a['FPARAM'][gd,3]
                    else : abun = a['FELEM'][gd,0,iel]
                except:
                    if etoh[iel] == 1 : abun=a['FELEM'][gd,iel]-a['FPARAM'][gd,3]
                    else : abun = a['FELEM'][gd,iel]

            fig,ax=plots.multi(1,2,hspace=0.001)
            plots.plotc(ax[0],a[param][gd,3],abun,a[param][gd,0],xr=[-2.5,1.0],yr=[-0.5,1],zr=te,xt='[M/H]',colorbar=True,zt='Teff',yt='['+el+'/M]')
            ax[0].text(0.1,0.9,'uncalibrated params'+comment,transform=ax[0].transAxes)
            plots.plotc(ax[1],a[param][gd,3],abun,a['SNR'][gd],xr=[-2.5,1.0],yr=[-0.5,1],zr=[50,200],xt='[M/H]',colorbar=True,zt='S/N',yt='['+el+'/M]')
            if out is not None :
                outfile=out+el+'_{:1d}.png'.format(icol)
                fig.savefig(outfile)
                row.append(os.path.basename(outfile))
            else: pdb.set_trace()
            plt.close(fig)
            icol+=1
            xt.append('{:6.0f}&lt;Teff&lt;{:6.0f} {:6.1f}&lt;logg&lt;{:6.1f}'.format(te[0],te[1],logg[0],logg[1]))
        grid.append(row)
    if out is not None : html.htmltab(grid,file=out+'elem_chem.html',ytitle=yt,xtitle=xt)

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
    if out is not None : html.htmltab(grid,file=out+'elem_err.html',ytitle=yt,xtitle=['giants','dwarfs'])

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
                    outfile=out+'flag_param_{:d}.png'.format(i)
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
