import glob
import os
from sdss import yanny
import numpy as np

def check(apred='r12',mjdmax=58360,mjdmin=55800) :

  for telescope in ['apo25m','lco25m'] :

    print(telescope)
    files= glob.glob(os.environ['APOGEEREDUCEPLAN_DIR']+'/pro/'+telescope+'/'+telescope+'_?????.pro')

    dir=os.environ['APOGEE_REDUX']+'/r12/visit/'+telescope+'/'
    plateplans = yanny.yanny(os.environ['PLATELIST_DIR']+'/platePlans.par')['PLATEPLANS']
    n=0
    for file in files :
        #print(file)
        fp=open(file,'r')
        mjd='0'
        for line in fp :
            if line.find('mjd=') >= 0 : mjd = line.split('=')[1].strip('\n')
            if int(mjd) <= mjdmax and int(mjd) >= mjdmin :
                if line.find('plate=') >= 0 and line[0] != ';' : 
                    plate = line.split('=')[1].strip('\n')
                    if plate != '0' : 
                        n+=1
                        try :
                            j = np.where(np.array(plateplans['plateid']) == int(plate))[0][0]
                            if plateplans['survey'][j] == 'manga-apogee2' :
                                field = plateplans['comments'][j]
                            else :
                                field = plateplans['name'][j].replace('APG_','')
                                field = field.replace('APGS_','')
                                field = field.replace('MC-','MC')

                            if not os.path.isfile(dir+field+'/apVisitSum-'+plate+'-'+mjd+'.fits') :
                                print('\n'+dir+field+'/apVisitSum-'+plate+'-'+mjd+'.fits') 
                                print(os.listdir(dir+field+'/'+plate+'/'+mjd))
                        except :
                            print('problem: ', line)
    print('ntotal: ', n)

