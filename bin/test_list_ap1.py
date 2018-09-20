#!/usr/bin/env python
'''
test for list_ap1.py



'''

import unittest
from list_ap1 import *

def ff(mjd, fexp):
    path="/data/apogee/quickred/%s" % mjd
    qrfile1="%s/ap1D-a-%s.fits.fz" % (path,fexp)
    return qrfile1

def ffraw(mjd, fexp):
    path="/data/apogee/utr_cdr/"
    fName="%s%s/apRaw-%s.fits"%(path,mjd,fexp)
    return fName 

class TestList(unittest.TestCase): 

    def test_offset(self):
        self.assertEqual(getOffset(ff("56660", "10980001")), " err") 
        self.assertEqual(getOffset(ff("56660", "10980002")), " 2.82")  # ok
        self.assertEqual(getOffset(ff("56660", "10980003")), " 0.77")  # ok
        self.assertEqual(getOffset(ff("56660", "10980004")), " 4.88")  # ok 
        self.assertEqual(getOffset(ff("56660", "10980005")), "-1.22")  # ok
        self.assertEqual(getOffset(ff("56660", "10980006")), " 6.93")  # ok
        self.assertEqual(getOffset(ff("56660", "10980007")), " 0.80")  # ok
        self.assertEqual(getOffset(ff("56660", "10980008")), " 5.88")  # ok
        self.assertEqual(getOffset(ff("56660", "10980009")), " 1.82")  # ok
# ... 
        self.assertEqual(getOffset(ff("56636", "10740003")), " 3.71")   # ok

        self.assertEqual(getOffset(ff("56531", "09690003")), "-0.24") 
        self.assertEqual(getOffset(ff("56531", "09690005")), " 0.24")
        

    def test_flux(self):
        self.assertEqual(getFlux(ffraw("56531", "09690019")), "100")  # InteFlat
        self.assertEqual(getFlux(ffraw("56531", "09690011")), " 99")  # QuarFlat
        self.assertEqual(getFlux(ffraw("56534", "09720044")), " 95")  # DomeFlat
        self.assertEqual(getFlux(ffraw("56534", "09720044")), " 95")  # DomeFlat
        self.assertEqual(getFlux(ffraw("56730", "11680022")), "479")  # DomeFlat
        self.assertEqual(getFlux(ffraw("56730", "11680025")), "  1")  # DomeFlat
        self.assertEqual(getFlux(ffraw("56730", "11680025")), "  1")  # DomeFlat
        self.assertEqual(getFlux(ffraw("56754", "11920031")), "149")  # DomeFlat


    def test_map(self):
    
        fNames="/data/apogee/utr_cdr/56721/apRaw-*.fits"   # raw data
        files=getFiles(fNames)
        self.assertEqual(len(files), 41)  # number of files for this mjd  
    
        map56721="ddqtutuddtttttttttttttttttdddqqqtutudiiid"         
        map=getMap(files)
        self.assertEqual(map, map56721)  # total mjd map
        
        mapRange=getSequence(map, morning)
        self.assertEqual(mapRange, [26,41])  # morning map
        files=files[mapRange[0]:mapRange[1]]
        self.assertEqual(len(files), 15)  # get 15 files for morning 15 calibrations
        self.assertEqual(files[0], '/data/apogee/utr_cdr/56721/apRaw-11590027.fits')  # 1st file
        self.assertEqual(files[14], '/data/apogee/utr_cdr/56721/apRaw-11590041.fits') # last file

        fNames="/data/apogee/utr_cdr/99999/apRaw-*.fits"   # raw data
        files=getFiles(fNames)
        self.assertEqual(getMap(files), None)  # no files found for map
        
        # test dither 
        fNames="/data/apogee/utr_cdr/56721/apRaw-*.fits"   # raw data
        files=getFiles(fNames)
        map=getMap(files)
        mapRange=getSequence(map, dither)
        self.assertEqual(mapRange, [9,26])  
            # dither map - start count from 0, files 10-26
        files=files[mapRange[0]:mapRange[1]]
        self.assertEqual(len(files), 17)
        
if __name__ == "__main__": 
    unittest.main() 


# I checked visual profiles and they are ok