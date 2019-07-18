import glob
import os
import pdb

def link(oldver='r12',dates='58*',newver='r13',fields='*', tels=['lco25m','apo1m','apo25m'] ) :
    """ Create links from existing reduction version into new version under current working directory
    """

    # exposure/TELESCOPE/MJD directories
    dirs=glob.glob(os.environ['APOGEE_REDUX']+'/'+oldver+'/exposures/*/'+dates+'/')
    mklinks(dirs,-4,-2,oldver=oldver)

    # cal/TELESCOPE/MJD directories
    dirs=glob.glob(os.environ['APOGEE_REDUX']+'/'+oldver+'/cal/*/'+dates+'/')
    mklinks(dirs,-4,-2,oldver=oldver)

    # visit/TELESCOPE/FIELD/PLATE/MJD directories and visit/TELESCOPE/FIELD/*VisitSum files
    for tel in tels :
        if tel == 'apo1m' :
            dirs=glob.glob(os.environ['APOGEE_REDUX']+'/'+oldver+'/visit/'+tel+'/*/'+dates+'/*')
            mklinks(dirs,-5,-1,oldver=oldver,newver=newver)
        else :
            dirs=glob.glob(os.environ['APOGEE_REDUX']+'/'+oldver+'/visit/'+tel+'/'+fields+'/*/'+dates+'/*')
            mklinks(dirs,-6,-1,oldver=oldver,newver=newver)
        files=glob.glob(os.environ['APOGEE_REDUX']+'/'+oldver+'/visit/'+tel+'/'+fields+'/*VisitSum*'+dates+'*')
        mklinks(files,-4,-1,oldver=oldver)

    # stars/TELESCOPE/FIELD/apStar and apField
    for tel in tels :
        files=glob.glob(os.environ['APOGEE_REDUX']+'/'+oldver+'/stars/'+tel+'/'+fields+'/a?Star*')
        mklinks(files,-4,-1,oldver=oldver,newver=newver)
        files=glob.glob(os.environ['APOGEE_REDUX']+'/'+oldver+'/stars/'+tel+'/'+fields+'/a?Field*')
        mklinks(files,-4,-1,oldver=oldver,newver=newver)
        files=glob.glob(os.environ['APOGEE_REDUX']+'/'+oldver+'/stars/'+tel+'/'+fields+'/plots/*.gif')
        mklinks(files,-5,-1,oldver=oldver,newver=newver)
        files=glob.glob(os.environ['APOGEE_REDUX']+'/'+oldver+'/stars/'+tel+'/'+fields+'/plots/*.jpg')
        mklinks(files,-5,-1,oldver=oldver,newver=newver)

    # calibration files
    for caldir in ['bpm', 'darkcorr','detector','flatcorr','flux','littrow','lsf','persist','psf','telluric','trace','wave'] :
        try : os.makedirs('cal/'+caldir)
        except : pass
        files =glob.glob(os.environ['APOGEE_REDUX']+'/'+oldver+'/cal/'+caldir+'/*')
        mklinks(files,-3,-1,oldver=oldver)

def aspcap(oldver='r12',newver='r13',oldaspcap='l33',newaspcap='l33',fields='*',tels=['lco25m','apo1m','apo25m']) :

    # visit/TELESCOPE/FIELD/PLATE/MJD directories and visit/TELESCOPE/FIELD/*VisitSum files
    for tel in tels :
        files=glob.glob(os.environ['APOGEE_ASPCAP']+'/'+oldver+'/'+oldaspcap+'/'+tel+'/'+fields+'/*')
        mklinks(files,-4,-1,oldver=oldver,newver=newver)

def mklinks(dirs,start,last,oldver='r12',absolute=False,newver=None,test=False) :
    """ Routine that actually does the linking
    """
    
    for dir in dirs :
        # get output name relative to top level
        out='/'.join(dir.split('/')[start:last])
        try: os.makedirs(out)
        except: pass

        if absolute :
            ref=dir
        else :
            # create relative link
            nlevels=len(out.split('/'))
            refdir='../'+oldver+'/'
            for i in range(nlevels) : refdir='../'+refdir
            ref=refdir+out+'/'+dir.split('/')[last]
        new =  out+'/'+dir.split('/')[last]
        new=new.replace('-'+oldver+'-','-'+newver+'-')
        if newver is not None and not test:
            try: os.remove(new)
            except: pass
            print(new)
        print('linking: ',ref, new,oldver,newver)
        if not test :
            try: os.remove(new)
            except: pass
            os.symlink(ref,new)
