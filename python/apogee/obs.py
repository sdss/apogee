import glob
import os
import pdb
import string
from astropy.io import fits

def done1m(program,vers='t9') :

    print('program: ', program)
    stars=fits.open(os.environ['APOGEEREDUCEPLAN_DIR']+'/data/1m/'+program+'.fits')[1].data
    print('Number of stars: ',len(stars))
    fp = open(program+'.html','w')
    fbad = open(program+'.bad','w')
    fp.write('<HTML><BODY><h2>'+program+'</h2><TABLE BORDER=2>')
    for star in stars :
        name = star['NAME']
        try :
            comment = star['COMMENT']
        except :
            comment = ''
        allobs = glob.glob(os.environ['APOGEE_REDUX']+'/'+vers+'/visit/apo1m/'+program+'/*/apVisit-*'+name+'*')
        snmax=-1
        for iobs,obs in enumerate(allobs) :
            v=fits.open(obs)[0].header
            mjd=obs.split('/')[-2]
            link=program+'/'+mjd+'/html/'+mjd+'-'+name+'sum.html'
            print(name,mjd,v['SNR'],link)
            if iobs == 0 :
                fp.write('<TR><TD>{:s}<TD>{:s}<TD><A HREF={:s}>{:s}</A><TD>{:8.2f}</A>\n'.format(name,comment,link,mjd,v['SNR']))
            else :
                fp.write('<TR><TD><TD><TD><A HREF={:s}>{:s}<TD>{:8.2f}</A>\n'.format(link,mjd,v['SNR']))
            if v['SNR'] > snmax : snmax = v['SNR']
        if snmax>0 and  snmax < 75 : fbad.write(name+'\n')
        if len(allobs) == 0 : 
                fp.write('<TR><TD BGCOLOR=red>{:s}<TD bgcolor=red>{:s}\n'.format(name,comment))
    fp.write('</TABLE></BODY></HTML>')
    fp.close()
    fbad.close()

def all() :
    progs = glob.glob('*')
    for prog in progs: 
       
        if os.path.isdir(prog) and prog != 'hip' and prog != 'DQTAU' : done1m(prog)
