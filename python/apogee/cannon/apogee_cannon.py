# Routines for use of The Cannon with APOGEE data
# 
# Basic sequence:
#    apogee_cannon.norm   : creates normalized spectra
#    apogee_cannon.train  : trains Cannon model
#    apogee_cannon.fit    : applies Cannon model
#    apogee_cannon.merge  : merges Cannon result with allStar file to create allStarCannon file
#  all routines run with input planfile that identifies version/field (so run in parallel for multiple fields with multiple planfiles)
#  and desired input parameters
#
# Auxiliary
#    apogee_cannon.send   : plots of Cannon sensitivies to individual label changes
#    apogee_cannon.compare : makes plots that show ASPCAP and Cannon results



import cPickle as pickle
import logging
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
import copy
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, TableColumns, Column
from glob import glob
from collections import OrderedDict
import subprocess


#import Cannon.AnniesLasso as tc
import AnniesLasso as tc
from apogee.cannon import continuum

from apogee.utils import apload
from apogee.utils import apselect
from apogee.aspcap import aspcap
from apogee.aspcap import cal
from apogee.aspcap import elem
from tools import match
from tools import struct
from tools import plots
from tools import html
from sdss_access.path import path

import argparse
import pdb
from apogee.utils import yanny

PICKLE_PROTOCOL = -1

def norm(planfile,threads=8,inter=False,sim=False) :
    '''
    Do Cannon normalization of files
    '''
    print('planfile: ', planfile)
    print('threads: ', threads)

    # construct list of apStar files to normalize and corresponding output files
    p=yanny.yanny(planfile,np=True)
    apred = p['apred_vers'].strip("'")
    apstar = p['apstar_vers'].strip("'")
    aspcap = p['aspcap_vers'].strip("'")
    results = p['results_vers'].strip("'")
    cannon = getval(p,'cannon_vers','cannon_aspcap')
    threads=int(getval(p,'ncpus','16'))
    

    # loop over fields in planfile
    for field in p['ASPCAP']['field'] :
        try:
            paths = getfiles(apred,apstar,aspcap,results,cannon,field)
            aspcappaths = getfiles(apred,apstar,aspcap,results,cannon,field,aspcapStar=True)
        except:
            return
        if sim :
            aspcappaths=paths
 
        N_individual_visits = len(paths)
        # create output directories
        for path in paths :
            outdir=os.path.dirname(path[2])
            if not os.path.exists(outdir):
                print('makedirs output_path',outdir)
                os.makedirs(outdir)

        # Enable logging.
        logger = logging.getLogger("apogee.dr14.tc")
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)-8s] %(message)s"))
        logger.addHandler(handler)

        # inspection
        if inter :
            norm_kwds = normalization_kwds()
            for path,aspcappath in zip(paths,aspcappaths): 
                stacked, visits, metadata = continuum.normalize_individual_visits(
                    path[1], full_output=True, **norm_kwds)
                plt.clf()
                plt.plot(np.arange(8575),stacked[0])
                plt.plot(np.arange(8575),1./np.sqrt(stacked[1]))
                stacked, visits, metadata = continuum.normalize_individual_visits(
                    aspcappath[1], full_output=True, apStar_sum=True, aspcapStar=True,ignore_bitmask_values=None,**norm_kwds)
                plt.plot(np.arange(8575),stacked[0])
                plt.plot(np.arange(8575),1./np.sqrt(stacked[1]))
                plt.ylim([0,1.3])
                plt.draw()
                pdb.set_trace()
            return       

        # Process the normalization in parallel.
        logger.info(
            "Starting {} threads to do pseudo-continuum-normalization of {} stars"\
            .format(threads, N_individual_visits))

        pool = mp.Pool(threads)
        normalized_result = \
            pool.map_async(_process_normalization, aspcappaths).get()
        pool.close()
        pool.join()

        # If there were input spectra that failed, then show a summary.
        _failed = []
        for result, input_path, output_path in normalized_result:
            if result not in (True, None):
                _failed.append(input_path)

        if _failed:
            logger.info("Summary of failures ({}):".format(len(_failed)))
            for input_path in _failed:
                logger.info("\t{}".format(input_path))

        # Save the dispersion from one spectrum.
        for result, input_path, output_path in normalized_result:
            if result in (True, None):
            
                with fits.open(input_path) as image:
                    dispersion = 10**(image[1].header["CRVAL1"] \
                        + np.arange(image[1].data.size) * image[1].header["CDELT1"])
            
                with open(os.path.join(outdir, "dispersion.pkl"), "wb") as fp:
                    pickle.dump(dispersion, fp, PICKLE_PROTOCOL)
            
                # We only need to save the dispersion once: break the loop.
                break

def normalization_kwds() :
    continuum_file = os.environ['IDLWRAP_DIR']+'/python/Cannon/continuum.list'
    return {
        "regions": [
            (15140, 15812),
            (15857, 16437),
            (16472, 16960),
        ],
        "conservatism": (2.0, 0.1),
        "normalized_ivar_floor": 1e-4,
        "continuum_pixels": np.loadtxt(continuum_file, dtype=int)
    }

# Create a function so that we can process the normalization in parallel.
def _process_normalization(path):
    """
    Produce pseudo-continuum-normalized data products for APOGEE DR14 apStar
    spectra.

    :param input_path:
        The local path of an apStar spectrum.

    :returns:
        A three-length tuple indicating: (1) whether the normalization was
        successful, (2) the `input_path`, and (3) the `output_path` if the
        normalization was successful. If `None` is provided in (1), it is 
        because the output file already exists and we were not instructed 
        to clobber it.
    """

    apogee_id, input_path,output_path = path

    print input_path, output_path
    #normalization set up
    clobber_normalization = False
    norm_kwds = normalization_kwds()

    print apogee_id

    # Check if this output file already exists.
    #if os.path.exists(output_path) and not clobber_normalization:
    #    logger.info("Skipping normalization of {}..".format(input_path))
    #    return (None, input_path, output_path)

    try:
        stacked, visits, metadata = continuum.normalize_individual_visits(
            input_path, full_output=True, aspcapStar=True, apStar_sum = True, **norm_kwds)

    except:
        #logger.exception("Normalization failed on {}".format(input_path))
        print("Normalization failed on {}".format(input_path))
        return (False, input_path, None)

    metadata.update(APOGEE_ID=apogee_id)

    stacked = np.vstack(stacked)
    visits = np.vstack(visits)
     
    with open(output_path, "wb") as fp:
        pickle.dump((metadata,stacked), fp, PICKLE_PROTOCOL)

    with open("{}.visits".format(output_path), "wb") as fp:
        pickle.dump(visits, fp, PICKLE_PROTOCOL)

    with open("{}.meta".format(output_path), "wb") as fp:
        pickle.dump(metadata, fp, PICKLE_PROTOCOL)

    #logger.info("Normalized spectra in {} successfully".format(input_path))
    print("Normalized spectra in {} successfully".format(input_path))

    return (True, input_path, output_path)

def getrange(val):
    '''
    from character string, return list of limits
    '''
    return [float(val.split()[0]),float(val.split()[1])]

def train(planfile,skip=1,threads=8,xh=None,model_name=None,censor=None,sim=False,gb=None,mh=None) :
    '''
    Define training set and train Cannon
    '''

    p=yanny.yanny(planfile,np=True)
    apred = p['apred_vers'].strip("'")
    apstar = getval(p,'apstar_vers','stars').strip("'")
    aspcap_vers = getval(p,'aspcap_vers','aspcap').strip("'")
    results = getval(p,'results_vers','results').strip("'")
    cannon = getval(p,'cannon_vers','cannon_aspcap')
    if model_name is None : model_name = getval(p,'model_name','apogee-dr14-giants')
    model_order = int(getval(p,'model_order','2'))
    model_scale_factor = float(getval(p,'model_scale_factor','1.0'))
    model_regularization = float(getval(p,'model_regularization','0.0'))
    threads=int(getval(p,'ncpus',threads))
    if xh is None: xh=int(getval(p,'xh',False))
    if censor is None: censor=int(getval(p,'censor',False))
    logg=getrange(getval(p,'logg','-1 3.9'))
    teff=getrange(getval(p,'teff','3500 5500'))
    if gb is None: gb=getval(p,'gb',0)
    if mh is None: mh=getrange(getval(p,'mh','-3. 1.'))
    alpha=getrange(getval(p,'alpha','-0.5 1.'))

    # label names
    elems = aspcap.elems()[0]
    #model_labels = ['TEFF','LOGG','M_H']
    model_labels = ['TEFF','LOGG','M_H','ALPHA_M','FE_H']
    input_labels = ['TEFF','LOGG','M_H','ALPHA_M','FE_H']
    for el in elems :
        d=elem.dr14cal(el)
        if el is not 'Fe' and d['elemfit'] >= 0:
            if xh :
                model_labels.append(el.upper()+'_H')
            else :
                model_labels.append(el.upper()+'_FE')
            input_labels.append(el.upper()+'_FE')

    apload.apred = apred
    apload.apstar = apstar
    apload.aspcap = aspcap_vers
    apload.results = results

    if sim :
        allstar=fits.open('allStar.fits')[1].data
        gd=apselect.select(allstar,logg=logg,teff=teff,mh=mh,alpha=alpha,sn=[100,10000])
        model_labels = ['TEFF','LOGG','M_H']
        input_labels = ['TEFF','LOGG','M_H']
        model_labels = sim
        input_labels = sim
        if gb :
          gd2 = np.where( np.abs( (allstar['TEFF'][gd]-3500)*4/2000. - allstar['LOGG'][gd])  < float(gb) )[0]
          gd = gd[gd2]
    else :
        allstar=apload.allStar()[1].data
        gd=apselect.select(allstar,badval=['STAR_BAD'],sn=[100,10000],logg=logg,teff=teff,mh=mh,alpha=alpha,badstar=['PERSIST_HIGH','PERSIST_MED','PERSIST_LOW'],gb=gb)

        gcstars = ascii.read(os.environ['IDLWRAP_DIR']+'/data/gc_szabolcs.dat')
        bd=np.where(gcstars['pop'] != 1)[0]
        jc = [x for x in gd if allstar[x]['APOGEE_ID'] not in gcstars['id'][bd]]
        gd=jc

        # down select stars using HR+[M/H] sampling
        i1,i2=cal.hrsample(allstar,allstar[gd],raw=False)
        # make sure all labels are good
        gd=[]
        for i in i1 :
           good = True
           for label in input_labels :
                # special handling for NA in DR14
                if label == 'NA_FE' and allstar[label][i] < -5 and allstar['FE_H'][i] < -1 : allstar[label][i] = 0.
                if allstar[label][i] < -5 : 
                    good = False
                    print('reject',allstar['APOGEE_ID'][i],label,allstar[label][i])
                    break
           if good : gd.append(i)
        
    print('selected ',len(gd),' training set stars')
    root = os.environ['APOGEE_ASPCAP']+'/'+apred+'/'+cannon+'/'
    training_set = os.path.join(root,"{}-training-set.fits".format(model_name))
    if not os.path.exists(os.path.dirname(training_set)):
           os.makedirs(os.path.dirname(training_set))
    struct.wrfits(np.array(allstar[gd]),training_set)

    # The label names to use in the model.

    model_filename = os.path.join(root, "{}.model".format(model_name))
    initial_filename = os.path.join(root, "{}.initial".format(model_name))

    clobber_model = True
    labelled_set = Table.read(training_set)[0:-1:skip]
    N_labelled = len(labelled_set)
    if xh :
        for el in elems :
            d=elem.dr14cal(el)
            if el is not 'Fe' and d['elemfit'] >= 0:
                labelled_set[el.upper()+'_H'] = labelled_set[el.upper()+'_FE']+labelled_set['FE_H']

    # TODO: something's wrong with our dispersion that we extracted.
    #with open(os.path.join(CANNON_DATA_DIR, "dispersion.pkl"), "rb") as fp:
    #    dispersion = pickle.load(fp)
    #P = dispersion.size
    dispersion = None
    P = 8575 # MAGIC

    # These defaults (flux = 1, ivar = 0) will mean that even if we don't find a
    # spectrum for a single star in the training set, then that star will just have
    # no influence on the training (since ivar = 0 implies infinite error on flux).

    normalized_flux = np.ones((N_labelled, P), dtype=float)
    normalized_ivar = np.zeros((N_labelled, P), dtype=float)

    # Enable logging.
    logger = logging.getLogger("apogee.dr14.tc")
    logger.setLevel(logging.INFO)

    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)-8s] %(message)s"))
    logger.addHandler(handler)

    sdss_path=path.Path()
    ngd = 0
    for i, row in enumerate(labelled_set):

        logger.info(
            "Reading labelled set spectra ({}/{})".format(i + 1, N_labelled))
        if row['TELESCOPE'] == 'apo1m' :
            filename = sdss_path.full('cannonStar-1m',apred=apred,apstar=apstar,aspcap=aspcap_vers,results=results,cannon=cannon,
                   field=row['FIELD'],reduction=row['REDUCTION_ID'],telescope=row['TELESCOPE'])
        else :
            filename = sdss_path.full('cannonStar',apred=apred,apstar=apstar,aspcap=aspcap_vers,results=results,cannon=cannon,
                   field=row['FIELD'],obj=row['APOGEE_ID'],telescope=row['TELESCOPE'])

        if not os.path.exists(filename):
            logger.warn("Could not find filename for labelled set star {}: {}"\
                .format(row["APOGEE_ID"], filename))
            continue

        with open(filename, "rb") as fp:
            #flux, ivar = pickle.load(fp)
            metadata, data = pickle.load(fp)
            flux,ivar = data

        if (np.isfinite(flux).all()) & (np.isfinite(ivar).all()) :
            normalized_flux[i, :] = flux
            normalized_ivar[i, :] = ivar
        else :
            print('non-finite values in',row['APOGEE_ID'])
            normalized_flux[i, :] = 0.
            normalized_ivar[i, :] = 0.
            #pdb.set_trace()

    # TODO: Cache the normalized_flux and normalized_ivar into a single file so that
    #       it is faster to read in next time?
    assert  np.isfinite(normalized_flux).all(), \
            "Non-finite values in normalized_flux!"
    assert  np.isfinite(normalized_ivar).all(), \
            "Non-finite values in normalized_ivar!"

    # Exclude labelled set stars where there is no spectrum, only because it
    # will get annoying later on when we are doing 1-to-1 and cross-validation
    keep = np.any(normalized_ivar > 0, axis=1)
    if not np.all(keep):
        logger.info(
            "Excluding {} labelled set stars where there was no information in "
            "the spectrum".format(np.sum(~keep)))
        labelled_set = labelled_set[keep]
        normalized_flux = normalized_flux[keep]
        normalized_ivar = normalized_ivar[keep]

    # Construct and train a model. #
    model = tc.L1RegularizedCannonModel(
        labelled_set, normalized_flux, normalized_ivar, dispersion, 
        threads=threads)

    model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(
        labelled_set, 
        tc.vectorizer.polynomial.terminator(model_labels, model_order),
        scale_factor=model_scale_factor)

    if censor :
        for label in model_labels :
            for el in elems :
                d=elem.dr14cal(el)
                if xh :
                    lab = el.upper()+'_H'
                else :
                    lab = el.upper()+'_FE'
                if lab == label :
                    model.censors[label] = getcensor(el,maskdir=os.environ['SPECLIB_DIR']+'/lib/filters_26042016/',length=P)
                    print(label,getcensor(el,maskdir=os.environ['SPECLIB_DIR']+'/lib/filters_26042016/'))

    model.s2 = 0
    model.regularization = model_regularization

    model.train()
    model._set_s2_by_hogg_heuristic()

    model.save(
        model_filename, include_training_data=False, overwrite=clobber_model)
    model.save(
        model_filename+'.full', include_training_data=True, overwrite=clobber_model)


    # Make some 1-to-1 plots just to show sensible behaviour.
    #X = model.labels_array()
    X = model.labels_array
    Y = model.fit(model.normalized_flux, model.normalized_ivar)
    out=Table(np.transpose([np.mean(model.labels_array, axis=0), np.nanmean(Y-X,axis=0), np.nanstd(Y-X,axis=0), np.array(model_labels)]) )
    out.write(initial_filename,overwrite=True,format='ascii') 
    #np.savetxt(initial_filename, [np.mean(model.labels_array, axis=0).reshape(-1, 1), np.nanmean(Y-X,axis=0), np.nanstd(Y-X,axis=0), model_labels] )

    try:
        os.makedirs(os.path.join(root,'plots'))
    except:
        pass
    it = model_labels.index('TEFF')
    ig = model_labels.index('LOGG')
    iz = model_labels.index('M_H')

    def plotit(ax,x,y,z,label) :
        plots.plotc(ax,x,y-x,z,xt=label,yt='inferred-labelled')
        lims = ax.get_xlim()
        ax.plot(lims, [0.,0.], c="#666666", zorder=-1, linestyle=":")
        mean, rms = np.nanmean(y - x), np.nanstd(y - x)
        title = "{}: ({:.2f}, {:.2f})".format(label, mean, rms)
        ax.set_title(title)

    fig,ax=plots.multi(2,3)
    plotit(ax[0,0],X[:,it],Y[:,it],X[:,iz],'TEFF')
    plotit(ax[0,1],X[:,ig],Y[:,ig],X[:,it],'LOGG')
    plotit(ax[1,0],X[:,iz],Y[:,iz],X[:,it],'M_H')
    plots.plotc(ax[1,1],X[:,it],X[:,ig],X[:,iz],xr=[6000,3500],yr=[5,-0.5],xt='TEFF',yt='LOGG')
    gd=np.where(model.normalized_ivar.flatten() >1)[0]
    diff=np.abs(model.normalized_flux-model.predict(Y))
    ax[2,0].hist(diff.flatten()[gd],cumulative=True,normed=True,bins=10.**np.arange(-8,0,0.05),histtype='step')
    ax[2,0].set_xlabel('|Model-true|')
    ax[2,0].set_ylim(0.,1.)
    ax[2,0].set_xscale('log')

    figure_path = os.path.join(root, "plots/{}-1to1.png".format(model_name))
    fig.tight_layout()
    fig.savefig(figure_path, dpi=300)
    plt.close()

    for i, label_name in enumerate(model_labels):

        x = X[:, i]
        y = Y[:, i]

        fig, ax = plt.subplots()
        if label_name == 'TEFF' :
            ax.scatter(x, y, c=X[:,iz], alpha=0.5)
        else :
            ax.scatter(x, y, c=X[:,it],alpha=0.5)

        lims = np.array([ax.get_xlim(), ax.get_ylim()])
        lims = (lims.min(), lims.max())
        ax.plot(lims, lims, c="#666666", zorder=-1, linestyle=":")
        ax.set_xlim(lims)
        ax.set_ylim(lims)

        ax.set_xlabel("Labelled")
        ax.set_ylabel("Inferred")

        mean, rms = np.nanmean(y - x), np.nanstd(y - x)
        title = "{}: ({:.2f}, {:.2f})".format(label_name, mean, rms)
        ax.set_title(title)
        logger.info("Mean and RMS for {}".format(title))

        figure_path = os.path.join(root, "plots/{}-{}-1to1.png".format(
            model_name, label_name))
        fig.tight_layout()
        fig.savefig(figure_path, dpi=300)
        plt.close()

        logger.info( "Created 1-to-1 figure for {} at {}".format(label_name, figure_path))

def getval(p,key,default) :
    try :
        return p[key].strip("'")
    except:
        return default

def fit(planfile, model_name=None, spectrum_filenames=None, threads=8, clobber=True, from_filename=False ,fit_velocity=False, chunk_size=1000,output_suffix=None,
    **kwargs):
    """
    Fit a series of spectra.
    """
    p=yanny.yanny(planfile,np=True)
    apred=p['apred_vers'].strip("'")
    apstar=p['apstar_vers'].strip("'")
    aspcap=p['aspcap_vers'].strip("'")
    results=p['results_vers'].strip("'")
    threads=int(getval(p,'ncpus','16'))
    cannon = getval(p,'cannon_vers','cannon_aspcap')
    if model_name is None: model_name = getval(p,'model_name','apogee-dr14-giants')
    if output_suffix is None: output_suffix = getval(p,'output_suffix','result')
    logg=getrange(getval(p,'logg','-1 3.9'))
    teff=getrange(getval(p,'teff','3500 5500'))
    mh=getrange(getval(p,'mh','-3. 1.'))
    alpha=getrange(getval(p,'alpha','-0.5 1.'))

    root = os.environ['APOGEE_REDUX']+'/'+apred+'/'+apstar+'/'+aspcap+'/'+results+'/'+cannon+'/'
    model = tc.load_model(os.path.join(root, "{}.model".format(model_name)), threads=threads)
    assert model.is_trained
    label_names = model.vectorizer.label_names
    mean_labels = Table.read(os.path.join(root, "{}.initial".format(model_name)),format='ascii')['col0']
    sig_labels = Table.read(os.path.join(root, "{}.initial".format(model_name)),format='ascii')['col2']
    #mean_labels = np.loadtxt(os.path.join(root, "{}.initial".format(model_name)))

    logger = logging.getLogger("AnniesLasso")

    # get allStar file for initial labels
    apload.apred = apred
    apload.apstar = apstar
    apload.aspcap = aspcap
    apload.results = results
    allstar=apload.allStar()[1].data

    # loop over fields in planfile
    for field in p['ASPCAP']['field'] :
        metadatas = []
        fluxes = []
        ivars = []
        output_filenames = []
        apogee_names = []
        failures = 0

        # get file names to fit 
        try:
            paths = getfiles(apred,apstar,aspcap,results,cannon,field)
        except:
            return
        
        spectrum_filenames = []
        initial_labels = []
        apogee_ids = []
        for apogee_id,inpath,outpath in paths :
            # only take stars within certain parameter ranges
            print(apogee_id)
            #j=apselect.select(allstar,redid=apogee_id)[0]
            j=np.where(((allstar['REDUCTION_ID'] == apogee_id) | (allstar['APOGEE_ID'] == apogee_id) ) & (allstar['COMMISS'] == 0) )[0]
            if (len(j) == 0) : 
                print('missing target',apogee_id)
            else:
                if len(j)>1 : j=j[0]
                if ((allstar['FPARAM'][j,1] >= logg[0]) & (allstar['FPARAM'][j,1] <= logg[1]) &
                    (allstar['FPARAM'][j,0] >= teff[0]) & (allstar['FPARAM'][j,0] <= teff[1]) &
                    (allstar['FPARAM'][j,3] >= mh[0]) & (allstar['FPARAM'][j,3] <= mh[1]) &
                    (allstar['FPARAM'][j,6] >= alpha[0]) & (allstar['FPARAM'][j,6] <= alpha[1]) ) :
                    spectrum_filenames.append(outpath)
                    apogee_names.append(apogee_id)
                    #labels=[]
                    #for i,label in enumerate(label_names) :
                    #    if allstar[label][j][0] > -9 :
                    #        labels.append(allstar[label][j][0])
                    #    else :
                    #        labels.append(mean_labels[i])
                    #initial_labels.append(labels)

        if len(apogee_names) == 0 : return

        #initial_labels=np.array(initial_labels)
        initial_labels = mean_labels
        # MAGIC HACK
        delete_meta_keys = ("fjac", ) # To save space...

        #output_suffix = kwargs.get("output_suffix", None)
        #output_suffix = "result" if output_suffix is None else str(output_suffix)
        summary_file =  root+field+'/cannonField-'+os.path.basename(field)+'-'+output_suffix+'.fits'
        N = len(spectrum_filenames)
        for i, names in enumerate(zip(apogee_names,spectrum_filenames)):
            apogee_id = names[0]
            filename = names[1]
            logger.info("At spectrum {0}/{1}: {2}".format(i + 1, N, filename))

            basename, _ = os.path.splitext(filename)
            output_filename = "-".join([basename, output_suffix]) + ".pkl"
        
            if os.path.exists(output_filename) and not clobber:
                logger.info("Output filename {} already exists and not clobbering."\
                    .format(output_filename))
                continue

            try:
                with open(filename, "rb") as fp:
                    metadata, data = pickle.load(fp)
                    metadatas.append(metadata)
                    flux,ivar = data
                    fluxes.append(flux)
                    ivars.append(ivar)

                output_filenames.append(output_filename)
                apogee_ids.append(apogee_id)

            except:
                logger.exception("Error occurred loading {}".format(filename))
                failures += 1

            else:
                if len(output_filenames) >= chunk_size:
                
                    results, covs, metas = model.fit(fluxes, ivars,
                        initial_labels=initial_labels, model_redshift=fit_velocity,
                        full_output=True)

                    for result, cov, meta, output_filename \
                    in zip(results, covs, metas, output_filenames):

                        for key in delete_meta_keys:
                            if key in meta:
                                del meta[key]

                        with open(output_filename, "wb") as fp:
                            pickle.dump((result, cov, meta), fp, 2) # For legacy.
                        logger.info("Saved output to {}".format(output_filename))
                
                    del output_filenames[0:], fluxes[0:], ivars[0:]


        if len(output_filenames) > 0:

            results, covs, metas = model.fit(fluxes, ivars, 
                initial_labels=initial_labels, model_redshift=fit_velocity,
                full_output=True)

            # Create an ordered dictionary of lists for all the data.
            data_dict = OrderedDict([("FILENAME", [])])
            data_dict['APOGEE_ID'] = []
            data_dict['LOCATION_ID'] = []
            data_dict['FIELD'] = []
            for label_name in label_names:
                data_dict[label_name] = []
            for label_name in label_names:
                data_dict["{}_RAWERR".format(label_name)] = []
            for label_name in label_names:
                data_dict["{}_ERR".format(label_name)] = []
            #data_dict["COV"] = []
            #meta_keys=metas[0].keys()
            meta_keys=['chi_sq','r_chi_sq','model_flux']
            for key in meta_keys:
                data_dict[key] = []
            data_dict['flux'] = []
            data_dict['ivar'] = []

            # loop over spectra, output individual files, and accumulate for summary file
            for result, cov, meta, output_filename,apogee_id,metadata,flux,ivar \
            in zip(results, covs, metas, output_filenames, apogee_ids,metadatas,fluxes,ivars):

              if np.isfinite(result).all() :
                outlist=[os.path.basename(output_filename),apogee_id,metadata['LOCATION_ID'],metadata['FIELD']]+result.tolist()
                try:
                    rawerr=np.diag(cov)**0.5
                    outlist.extend(rawerr)
                except:
                    pdb.set_trace()
                outlist.extend(np.max([rawerr,sig_labels],axis=0))
                #outlist.append(cov.tolist())
                for key in delete_meta_keys:
                    if key in meta:
                        del meta[key]
                #outlist += [meta.get(k, v) for k, v in meta.items()]
                outlist += [meta.get(k) for k in meta_keys]
                outlist.append(flux)
                outlist.append(ivar)
                for key, value in zip(data_dict.keys(), outlist):
                    data_dict[key].append(value)

                # save to pkl file?
                #with open(output_filename, "wb") as fp:
                #    pickle.dump((result, cov, meta), fp, 2) # For legacy.
                #logger.info("Saved output to {}".format(output_filename))

                # save to FITS cannonStar file
                hdr = fits.Header()
                hdr['HISTORY'] = 'IDLWRAP_VERSION: '+subprocess.check_output('idlwrap_version').strip('\n')
                hdr['OBJ'] = apogee_id
                hdr['LOCID'] = metadata['LOCATION_ID']
                hdr['FIELD'] = metadata['FIELD']
                hdr['CHI2'] = meta.get('r_chi_sq')
                for i,label_name in enumerate(label_names) :
                    hdr[label_name] = result[i]
                hdulist = fits.HDUList(fits.PrimaryHDU(header=hdr))
                hdr = fits.Header()
                hdr['OBSERVER'] = 'Edwin Hubble'
                hdr['CRVAL1'] = 4.179e0
                hdr['CDELT1'] = 6.e-6
                hdr['CRPIX1'] = 1
                hdr['CTYPE1'] = 'LOG-LINEAR'
                hdr['DC-FLAG'] = 1
                hdulist.append(fits.ImageHDU(flux,header=hdr))
                hdulist.append(fits.ImageHDU(1./np.sqrt(ivar),header=hdr))
                hdulist.append(fits.ImageHDU(meta.get('model_flux'),header=hdr))
                hdulist.writeto(output_filename.replace('-result','').replace('.pkl','.fits'),overwrite=True)
        
            del output_filenames[0:], fluxes[0:], ivars[0:]


        logger.info("Number of failures: {}".format(failures))
        logger.info("Number of successes: {}".format(N - failures))
        table = Table(TableColumns(data_dict))
        table.write(summary_file.replace('-result',''), overwrite=clobber)
        logger.info("Written to {}".format(summary_file))


    return None

def getfiles(apred,apstar,aspcap,results,cannon,field,aspcapStar=False) :

    # construct list of apStar files to normalize and corresponding output files
    root = os.environ['APOGEE_REDUX']+'/'+apred+'/'+apstar
    sdss_path=path.Path()
    paths = []

    apfieldfile = root+'/'+field+'/apField-'+os.path.basename(field)+'.fits'
    try :
        apfield=fits.open(apfieldfile)[1].data
    except:
        print('ERROR reading file',apfieldfile)
        raise
        return
    if aspcapStar :
        root = 'aspcapStar'
    else :
        root = 'apStar'

    for star in apfield['APOGEE_ID'] :
        if apfield['TELESCOPE'][0] == 'apo1m' :
           infile = sdss_path.full(root+'-1m',apred=apred,apstar=apstar,aspcap=aspcap,results=results,prefix='ap',
               field=apfield['FIELD'][0],reduction=star,telescope=apfield['TELESCOPE'][0])
           outfile = sdss_path.full('cannonStar-1m',apred=apred,apstar=apstar,aspcap=aspcap,results=results,cannon=cannon,
               field=apfield['FIELD'][0],reduction=star,telescope=apfield['TELESCOPE'][0])
        else :
           infile = sdss_path.full(root,apred=apred,apstar=apstar,aspcap=aspcap,results=results,prefix='ap',
               field=apfield['FIELD'][0],obj=star,telescope=apfield['TELESCOPE'][0])
           outfile = sdss_path.full('cannonStar',apred=apred,apstar=apstar,aspcap=aspcap,results=results,cannon=cannon,
               field=apfield['FIELD'][0],obj=star,telescope=apfield['TELESCOPE'][0])
        paths.append((star,infile,outfile))
    return paths

def merge(planfile,fields=None,outfile=None,clobber=True) :
    '''
    Match Cannon results to existing allStar file to make new table
    '''
    p=yanny.yanny(planfile,np=True)
    apred=p['apred_vers'].strip("'")
    apstar=p['apstar_vers'].strip("'")
    aspcap_vers=p['aspcap_vers'].strip("'")
    results=p['results_vers'].strip("'")
    apload.apred = apred
    apload.apstar = apstar
    apload.aspcap = aspcap_vers
    apload.results = results
    a=apload.allStar()[1].data
    t=Table(a)

    out = Table()
    out['APOGEE_ID'] = t['APOGEE_ID']
    length=len(out)

    if fields is None :
        fields=glob('*/cannonField*.fits')
    else :
        fields=glob(fields)

    c=fits.open(fields[0])[1].data
    for i,name in enumerate(c.names) :
        print name
        if name != 'APOGEE_ID'  and name != 'model_flux' and name != 'fvec' and name != 'flux' and name != 'ivar' :
            out.add_column(Column(name=name,dtype=c.dtype[i],length=length))
            print(name,type(out[name][0]))
            if type(out[name][0]) is np.string_: 
                print 'str!'
                out[name] = ''
            else :
                out[name] = -9999.
    # add X_M tag
    #out.add_column(Column(name='X_M',dtype='{:d}f4'.format(len(a['X_M'])),length=length))

    for field in fields :
        print('field',field)
        c=fits.open(field)[1].data
        j1=np.where(a['FIELD'] == c['FIELD'][0])[0]
        i1,i2=match.match(a['APOGEE_ID'][j1],c['APOGEE_ID'])
        bd = np.where(c['chi_sq'][i2] <= 0.)[0]
        print(len(bd))
        for name in out.columns :
            out[name][j1[i1]] = c[name][i2]
            if  type(c[name][0]) is np.string_: 
                bad = ''
            else :
                bad = -9999.
            if name is not 'APOGEE_ID' : out[name][j1[i1[bd]]] = bad

    # for DR14, "fix" NA_H
    bd = np.where(out['NA_H'] < -1)[0]
    out['NA_H'][bd] = -9999.

    out['CANNON_ID'] = t['ASPCAP_ID']
    if outfile is None :
        outfile='allStarCannon-'+results+'.fits'

    prihdr=fits.Header()
    prihdr['HISTORY'] = 'IDLWRAP_VERSION: '+subprocess.check_output('idlwrap_version').strip('\n')
    prihdu=fits.PrimaryHDU(header=prihdr)
    hdu=fits.BinTableHDU.from_columns(np.array(out))
    hdulist=fits.HDUList([prihdu,hdu])
    hdulist.writeto(outfile,overwrite=clobber)
    ##out.write(outfile,overwrite=clobber)
    return out

def compare(planfile,model_name=None,outfile=None,xh=False,output_suffix='') :
    '''
    Make some plots with results
    ''' 
    p=yanny.yanny(planfile,np=True)
    apred=p['apred_vers'].strip("'")
    apstar=p['apstar_vers'].strip("'")
    aspcap_vers=p['aspcap_vers'].strip("'")
    results=p['results_vers'].strip("'")
    if model_name is None: model_name = getval(p,'model_name','apogee-dr14-giants')
    apload.apred = apred
    apload.apstar = apstar
    apload.aspcap = aspcap_vers
    apload.results = results
    a=apload.allStar()[1].data

    if outfile is None :
        outfile='allStarCannon-'+results+'.fits'
    c = fits.open(outfile)[1].data

    elems = aspcap.elems()[0]
    elems = ['Ca','Ni']
    model_labels = ['TEFF','LOGG','M_H','ALPHA_M','FE_H']
    figs=[]
    ytit=[]
    for el in elems :
        d=elem.dr14cal(el)
        if el is not 'Fe' and d['elemfit'] >= 0:
          print(el)
          tag = el.upper()+'_FE'
          if xh :
              ctag = el.upper()+'_H'
          else :
              ctag = el.upper()+'_FE'
          model_labels.append(ctag)
          f=[]
          xtit=[]
          for tmin in range(3500,5500,500) :
            fig,ax=plots.multi(4,2,hspace=0.001,wspace=0.001,figsize=(10,4.5))
            for i,snmin in enumerate(range(50,250,50)) :
                #j = apselect.select(a,badval='STAR_BAD',sn=[snmin,100000],teff=[tmin,tmin+500],logg=[-1,3.9],mh=[-3.,1],alpha=[-0.5,1.],raw=True)
                j = np.where((np.chararray.find(a['ASPCAPFLAGS'],'STAR_BAD')<0) & (c['TEFF']>tmin) & (c['TEFF']<tmin+500) & (a['SNR']>snmin) & (a['SNR']<snmin+50))[0]
                if i == 0 : 
                  yt='['+el+'/Fe]'
                else : 
                  yt=None
                ax1=plots.plotc(ax[0,i],a['FE_H'][j],a[tag][j],a['ASPCAP_CHI2'][j],xr=[-2.5,1],yr=[-0.5,0.75],zr=[0,10],xt='[Fe/H]',yt=yt,nxtick=6,rasterized=True)
                ax[0,i].text(0.1,0.9,'ASPCAP:',transform=ax[0,i].transAxes)
                ax[0,i].text(0.2,0.8,'S/N > {:d}'.format(snmin),transform=ax[0,i].transAxes)
                if xh :
                    ax2=plots.plotc(ax[1,i],c['FE_H'][j],c[ctag][j]-c['FE_H'][j],c['r_chi_sq'][j],xr=[-2.5,1],yr=[-0.5,0.75],zr=[0,5],xt='[Fe/H]',yt=yt,nxtick=6,rasterized=True)
                else :
                    ax2=plots.plotc(ax[1,i],c['FE_H'][j],c[ctag][j],c['r_chi_sq'][j],xr=[-2.5,1],yr=[-0.5,0.75],zr=[0,5],xt='[Fe/H]',yt=yt,nxtick=6,rasterized=True)
                ax[1,i].text(0.1,0.9,'Cannon:',transform=ax[1,i].transAxes)
                ax[1,i].text(0.2,0.8,'{:d} < S/N < {:d}'.format(snmin,snmin+50),transform=ax[1,i].transAxes)
            file = '{:s}_{:04d}{:s}'.format(el,tmin+250,output_suffix)
            cbaxes = fig.add_axes([0.91, 0.55, 0.01, 0.3])
            cb = plt.colorbar(ax1, cax = cbaxes)
            cb.set_label('CHI2')
            cbaxes.tick_params(axis='both',labelsize=8)
            cbaxes = fig.add_axes([0.91, 0.15, 0.01, 0.3])
            cb = plt.colorbar(ax2, cax = cbaxes)
            cb.set_label('CHI2')
            cbaxes.tick_params(axis='both',labelsize=8)

            #fig.savefig('newplots/'+file+'.jpg',dpi=300)
            fig.savefig('newplots/'+file+'.pdf')
            plt.close()
            f.append(file)
            xtit.append('{:04d} < Teff < {:04d}'.format(tmin,tmin+500))
          figs.append(f)
          ytit.append(el)
    html.htmltab(figs,file='plots/cannon'+output_suffix+'.html',ytitle=ytit,xtitle=xtit)        
    figs=[]
    ytit=[]
    for label in model_labels :
        print(label)
        fig,ax=plots.multi(1,1,figsize=(12,8))
        try :
            gd=np.where((a[label]>-999) & (a[label]<999) & (c[label] > -999))[0]
            plots.plotc(ax,a[label][gd],c[label][gd],a['M_H'][gd],zr=[-2,0.5],colorbar=False,xt=label)
        except :
            alabel=label.replace('_H','_FE')
            print(alabel)
            gd=np.where((a[alabel]>-999) & (a[alabel]<999) & (c[label] > -999))[0]
            plots.plotc(ax,a[alabel][gd]+a['FE_H'][gd],c[label][gd],a['M_H'][gd],zr=[-2,0.5],colorbar=False,xt=label)
        file = 'fullcomp_'+label+output_suffix+'.jpg'
        fig.tight_layout()
        fig.savefig('newplots/'+file,dpi=300)
        plt.close()
        f=[file,model_name+'-'+label+'-1to1.png']
        figs.append(f)
        ytit.append(label)
    html.htmltab(figs,file='newplots/compare'+output_suffix+'.html',ytitle=ytit)

def sens(model='apogee-dr14-giants.model',teff=4500,logg=3,mh=0.,xh=False):
    '''
    Plot sensitivities of spectra to changes in Cannon labels
    '''
    model=tc.load_model(model)
    # get label names
    label_names = model.vectorizer.label_names
    
    # set default parameters
    if xh :
      d0 = mh
    else :
      d0 = 0.
    pars = [teff,logg,mh,0.,mh,d0,d0,d0,d0,d0,d0,d0,d0,d0,d0,d0,d0,d0,d0,d0,d0,d0,d0,d0] #,d0]
    synth0=model.predict(pars)[0]
    fig,ax=plots.multi(1,1)
    axalt=ax.twinx()
    wave=10.**(4.179+6.e-6*np.arange(8575))

    figs=[] 
    for i,name in enumerate(label_names) :
        print(name)
        delta = np.zeros(len(pars)) 
        if i == 0 : 
            delta[i]+=100.
        else  :
            delta[i]+=0.1
        synth=model.predict(pars+delta)[0]
        ax.cla()
        axalt.cla()
        plots.plotl(ax,wave,synth-synth0,yr=[-0.10,0.05])
        ax.text(0.1,0.9,name,transform=ax.transAxes)
        if name not in ['TEFF','LOGG','M_H','ALPHA_M'] :
            elem=name.split('_')[0]
            print(elem,' ',name)
            aspcap.elemsens([elem],plot=axalt,teff=teff,logg=logg,feh=mh,smooth=2,ylim=[0.05,-0.10])
            file=name+'_sens.jpg'
            fig.savefig('plots/'+file)
            figs.append([file])
            pdb.set_trace()
        #pdb.set_trace()
    html.htmltab(figs,file='plots/sens.html')

def getcensor(el,maskdir='./',length=8575) :
    '''
    Create censor list from ASPCAP mask files and wave.dat file that gives pixel number association
    '''
    pix=ascii.read(maskdir+'/wave.dat',Reader=ascii.NoHeader)['col2']
    filt=ascii.read(maskdir+'/'+el+'.filt',Reader=ascii.NoHeader)['col1']

    start=0
    cens=[]
    for i in range(len(filt) ):
      if filt[i] > 0 and start >= 0:
        skip=[start,pix[i]]
        cens.append(skip)
        start=-1
      elif filt[i] == 0. and start <0 :
        start = pix[i]

    cens.append([start,length-1])
    return cens

def allsim() :
    train('grid.par',sim=True,xh=False,model_name='grid')
    train('grid.par',sim=True,xh=False,model_name='grid_gb1',gb=1)
    train('grid.par',sim=True,xh=False,model_name='grid_gb2',gb=2)
    train('grid.par',sim=True,xh=False,model_name='grid_no_low_z',mh=[-1.01,1.])
    train('grid.par',sim=True,xh=False,model_name='grid_gb1_no_low_z',mh=[-1.01,1.],gb=1)
    train('grid.par',sim=True,xh=False,model_name='grid_gb2_no_low_z',mh=[-1.01,1.],gb=2)
