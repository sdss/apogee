"""
Routines for handling ASPCAP bitmasks
"""

import numpy as np
import pdb

class BitMask():
    '''
    Base class for bitmasks, define common methods
    
    At a minimum, a BitMask will have a set of name, level, descrip

    BitMask provides 3 methods:
        getname(val,level=level) : returns name(s) of all set bits
                                  (optionally, of requested level)
        getval(name) : returns value of bit with input name
        badval()     : returns value of all bits that are marked bad (level=1)
        warnval()    : returns value of all bits that are marked warn (level=2)
    '''
    def getname(self,val,level=0,strip=True):
        '''
        Given input value, returns names of all set bits, optionally of a given level
        '''
        strflag=''
        for ibit,name in enumerate(self.name) :
            if (name != 'RESERVED' ) :
                try:
                    if ( val & 2**ibit ) > 0 and ( level == 0 or self.level == level ) :
                      strflag = strflag + name +','
                except: 
                    print('bit problem: ', ibit)
                    pdb.set_trace()
        if strip : return strflag.strip(',')
        else : return strflag

    def getval(self,name) :
        """
        Get the numerical bit value of a given character name(s)
        """
        if type(name) is str :
            name = [name]
        bitval = np.int64(0)
        for n in name :
            try:
                j=self.name.index(n.strip())
                bitval|=np.int64(2**j)
            except :
                print('WARNING: undefined name: ',n)
        return bitval


    def badval(self) :
        """
        Return bitmask value of all bits that indicate BAD in input bitmask
        """
        val=np.int64(0)
        for i,level in enumerate(self.level) :
            if level == 1 :
                try: val=val | np.int64(2**i)
                except: pdb.set_trace()
        return val

    def warnval(self) :
        """
        Return bitmask value of all bits that indicate BAD in input bitmask
        """
        val=np.int64(0)
        for i,level in enumerate(self.level) :
            if level == 2 :
                val=val | np.int64(2**i)
        return val

    def print(self,fmt='txt') :
        """ Formatted output of bit definitions
        """
        if fmt == 'txt' :
            print('{:25s}{:>6s}  {:s}'.format('Name','Bit','Description'))
        elif fmt == 'wiki' :
            print('||{:25s}||{:>6s}||{:s}||'.format('Name','Bit','Description'))
        elif fmt == 'latex' :
            print('{:25s}&{:>6s}&{:s}\\\\'.format('Name','Bit','Description'))
        for ibit,name in enumerate(self.name) :
            if (name != 'RESERVED' and name != '' ) :
                if fmt == 'txt' :
                    print('{:25s}{:6d}  {:s}'.format(name,ibit,self.descrip[ibit]))
                elif fmt == 'wiki' :
                    print('||{:25s}||{:6d}||{:s}||'.format(name,ibit,self.descrip[ibit]))
                elif fmt == 'latex' :
                    print('{:25s}&{:6d}&{:s}\\\\'.format(name,ibit,self.descrip[ibit]))

class StarBitMask(BitMask):
    '''
    BitMask class for APOGEE star bitmask (APOGEE_STARFLAG)
    '''

    name=(['BAD_PIXELS','COMMISSIONING','BRIGHT_NEIGHBOR','VERY_BRIGHT_NEIGHBOR','LOW_SNR','','','',
          '','PERSIST_HIGH','PERSIST_MED','PERSIST_LOW','PERSIST_JUMP_POS','PERSIST_JUMP_NEG','','',
          'SUSPECT_RV_COMBINATION','SUSPECT_BROAD_LINES','BAD_RV_COMBINATION','RV_REJECT','RV_SUSPECT','MULTIPLE_SUSPECT','RV_FAIL','',
          'MTPFLUX_LT_75','MTPFLUX_LT_50','','','','','','RESERVED'])
    level=([1,0,0,1,0,0,0,0,
             0,0,0,0,0,0,0,0,
             0,0,1,0,0,0,1,0,
             0,0,0,0,0,0,0,0])

    descrip=([
     'Spectrum has many bad pixels (>20%):  BAD',                                                         
     'Commissioning data (MJD<55761), non-standard configuration, poor LSF: WARN',                       
     'Star has neighbor more than 10 times brighter: WARN',
     'Star has neighbor more than 100 times brighter: BAD',
     'Spectrum has low S/N (S/N<5)',                                                                    
     '',
     '',
     '',
     '',
     'Spectrum has significant number (>20%) of pixels in high persistence region: WARN',               
     'Spectrum has significant number (>20%) of pixels in medium persistence region: WARN',
     'Spectrum has significant number (>20%) of pixels in low persistence region: WARN',
     'Spectrum show obvious positive jump in blue chip: WARN',
     'Spectrum show obvious negative jump in blue chip: WARN',                                         
     '',
     '',
     'RVs from synthetic template differ significantly (~2 km/s) from those from combined template: WARN', 
     'Cross-correlation peak with template significantly broader than autocorrelation of template: WARN',
     'RVs from synthetic template differ very significatly (~10 km/s) from those from combined template: BAD',
     'Rejected visit because cross-correlation RV differs significantly from least squares RV',
     'Suspect visit (but used!) because cross-correlation RV differs slightly from least squares RV',
     'Suspect multiple components from Gaussian decomposition of cross-correlation',
     'RV failure',
     '',
     'Spectrum falls on fiber in MTP block with relative flux < 0.75',
     'Spectrum falls on fiber in MTP block with relative flux < 0.5',
     '',
     '',
     '',
     '',
     '',
     ''
     ])
    def persist(self) :
        '''
        Returns bitwise OR of all persistence bits
        '''
        return self.getval(['PERSIST_HIGH','PERSIST_MED','PERSIST_LOW',
                             'PERSIST_JUMP_POS','PERSIST_JUMP_NEG'] )

class AspcapBitMask(BitMask):
    '''
    BitMask class for APOGEE ASPCAP bitmask (APOGEE_ASPCAPFLAG)
    '''
    name=(['TEFF_WARN','LOGG_WARN','VMICRO_WARN','M_H_WARN','ALPHA_M_WARN','C_M_WARN','N_M_WARN','STAR_WARN',
          'CHI2_WARN','COLORTE_WARN','ROTATION_WARN','SN_WARN','SPEC_HOLE_WARN','ATMOS_HOLE_WARN','VSINI_WARN','',
          'TEFF_BAD','LOGG_BAD','VMICRO_BAD','M_H_BAD','ALPHA_M_BAD','C_M_BAD','N_M_BAD','STAR_BAD',
          'CHI2_BAD','COLORTE_BAD','ROTATION_BAD','SN_BAD','SPEC_HOLE_BAD','ATMOS_HOLE_BAD','VSINI_BAD','NO_ASPCAP_RESULT',
          'MISSING_APSTAR','NO_GRID','BAD_FRAC_LOWSNR','BAD_FRAC_BADPIX','FERRE_FAIL','','','',
          'PROBLEM_TARGET','','','','','','','',
          '','','','','','','','',
          '','','','','','','','RESERVED'])
    level=([2,2,0,0,0,0,0,2,
            2,2,2,2,2,2,0,0,
            1,1,0,0,0,0,0,1,
            1,1,1,1,1,2,0,1,
            1,1,1,1,1,0,0,0,
            1,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0])

    descrip=([
     'WARNING on effective temperature (see PARAMFLAG[0] for details) ',
     'WARNING on log g (see PARAMFLAG[1] for details) ',
     'WARNING on vmicro (see PARAMFLAG[2] for details) ',
     'WARNING on metals (see PARAMFLAG[3] for details) ',
     'WARNING on [alpha/M] (see PARAMFLAG[4] for details) ',
     'WARNING on [C/M] (see PARAMFLAG[5] for details) ',
     'WARNING on [N/M] (see PARAMFLAG[6] for details) ',
     'WARNING overall for star: set if any of TEFF, LOGG, CHI2, COLORTE, ROTATION, SN warn are set ',
     'high chi^2 (> 2*median at ASPCAP temperature (WARN)',
     'effective temperature more than 500K from photometric temperature for dereddened color (WARN)',
     'Spectrum has broad lines, with possible bad effects: FWHM of cross-correlation of spectrum with best RV template relative to auto-correltion of template > 1.5 (WARN)',
     'S/N<70 (WARN)',
     'Grid point within 2 grid steps of hole-filled synthesis ',
     'Grid point within 2 grid steps of hole-filled atmosphere ',
     ' ',
     ' ',
     'BAD effective temperature (see PARAMFLAG[0] for details) ',
     'BAD log g (see PARAMFLAG[1] for details) ',
     'BAD vmicro (see PARAMFLAG[2] for details) ',
     'BAD metals (see PARAMFLAG[3] for details) ',
     'BAD [alpha/M] (see PARAMFLAG[4] for details) ',
     'BAD [C/M] (see PARAMFLAG[5] for details) ',
     'BAD [N/M] (see PARAMFLAG[6] for details) ',
     'BAD overall for star: set if any of TEFF, LOGG, CHI2, COLORTE, ROTATION, SN error are set, or any GRIDEDGE_BAD ',
     'high chi^2 (> 5*median at ASPCAP temperature (BAD)',
     'effective temperature more than 1000K from photometric temperature for dereddened color (BAD)',
     'Spectrum has broad lines, with possible bad effects: FWHM of cross-correlation of spectrum with best RV template relative to auto-correltion of template > 2 (BAD)',
     'S/N<50 (BAD)',
     'Grid point within 1 grid steps of hole-filled synthesis ',
     'Grid point within 1 grid steps of hole-filled atmosphere ',
     ' ',' ',
     'Missing apStar file','Not processed by any ASPCAP grid','Fraction low SNR pixels > 0.5',
     'Fraction bad pixels > 0.5 or 0.33 in any chip','FERRE failure (bad input?)','','','',
     'Target extended or not suitable for standard star analysis','','','','','','','',
     '','','','','','','','',
     '','','','','','','',''
     ])

class ParamBitMask(BitMask):
    '''
    BitMask class for APOGEE ASPCAP bitmask (APOGEE_ASPCAPFLAG)
    '''
    name =['GRIDEDGE_BAD','CALRANGE_BAD','OTHER_BAD','FERRE_FAIL','PARAM_MISMATCH_BAD','FERRE_ERR_USED','','',
           'GRIDEDGE_WARN','CALRANGE_WARN','OTHER_WARN','FERRE_WARN','PARAM_MISMATCH_WARN','OPTICAL_WARN','ERR_WARN','FAINT_WARN',
           'PARAM_FIXED','RV_WARN','','','','','','',
           'LOGG_CAL_RC','LOGG_CAL_RGB','LOGG_CAL_MS','LOGG_CAL_RGB_MS','','','','RESERVED']


    level=[1,1,1,1,1,0,0,0,
             0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0]

    descrip=[
     'Parameter within 1/8 grid spacing of grid edge ',
     'Parameter outside valid range of calibration determination ',
     'Other error condition ',
     'Failed solution in FERRE ',
     'Elemental abundance from window differs significantly from parameter abundance ',
     'FERRE uncertainty used (larger than parametric uncertainty) ',
     ' ',
     ' ',
     'Parameter within 1 grid spacing of grid edge ',
     'Parameter in possibly unreliable range of calibration determination ',
     'Other warning condition ',
     'FERRE warning (not implemented?) ',
     'Elemental abundance from window differs from parameter abundance ',
     'Comparison with optical abundances suggests problem ',
     'Large expected uncertainty or upper limit based on location in parameter space (Teff, [M/H], S/N) ',
     'Warning based on faint star/RV combination ',
     'Parameter set at fixed value, not fit',
     'RV puts important line off of chip ',
     ' ',
     ' ',
     ' ',
     ' ',
     ' ',
     ' ',
     'Use RC gravity calibration ',
     'Use RGB gravity calibration ',
     'Use MS gravity calibration ',
     'Use RBG/MS transition gravity calibration ',
     ' ',
     ' ',
     ' ',
     ' '
     ]

class PixelBitMask(BitMask) :
    '''
    BitMask class for APOGEE pixel bitmask (APOGEE_PIXMASK)
    '''
    name=(['BADPIX','CRPIX','SATPIX','UNFIXABLE','BADDARK','BADFLAT','BADERR','NOSKY',
          'LITTROW_GHOST','PERSIST_HIGH','PERSIST_MED','PERSIST_LOW','SIG_SKYLINE','SIG_TELLURIC','NOT_ENOUGH_PSF','',
          'FERRE_MASK','','','','','','','',
          '','','','','','','','RESERVED'])

    level=([1,1,1,1,1,1,1,1,
            0,0,0,0,0,0,1,0,
            0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0])

    maskcontrib=([0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
                 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.0,
                 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.0,
                 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.0])

    descrip=([
     'Pixel marked as BAD in bad pixel mask or from strong persistence jump',
     'Pixel marked as cosmic ray in ap3d',
     'Pixel marked as saturated in ap3d',
     'Pixel marked as unfixable in ap3d',
     'Pixel marked as bad as determined from dark frame',
     'Pixel marked as bad as determined from flat frame',
     'Pixel set to have very high error (not used)',
     'No sky available for this pixel from sky fibers',
     'Pixel falls in Littrow ghost, may be affected',
     'Pixel falls in high persistence region, may be affected',
     'Pixel falls in medium persistence region, may be affected',
     'Pixel falls in low persistence region, may be affected',
     'Pixel falls near sky line that has significant flux compared with object',
     'Pixel falls near telluric line that has significant absorption',
     'Less than 50 percent PSF in good pixels',
     '',
     'Pixel masked by FERRE mask < 0.001',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     ''
    ])

class Apogee2Target1(BitMask) :
    '''
    BitMask class for APOGEE2_TARGET1
    '''

    name = ([ 'APOGEE2_ONEBIN_GT_0_5','APOGEE2_TWOBIN_0_5_TO_0_8','APOGEE2_TWOBIN_GT_0_8','APOGEE2_IRAC_DERED',
              'APOGEE2_WISE_DERED','APOGEE2_SFD_DERED','APOGEE2_NO_DERED','APOGEE2_WASH_GIANT',
              'APOGEE2_WASH_DWARF','APOGEE2_SCI_CLUSTER','APOGEE2_CLUSTER_CANDIDATE','APOGEE2_SHORT',
              'APOGEE2_MEDIUM','APOGEE2_LONG','APOGEE2_NORMAL_SAMPLE','APOGEE2_MANGA_LED',
              'APOGEE2_ONEBIN_GT_0_3','APOGEE2_WASH_NOCLASS','APOGEE2_STREAM_MEMBER','APOGEE2_STREAM_CANDIDATE',
              'APOGEE2_DSPH_MEMBER','APOGEE2_DSPH_CANDIDATE','APOGEE2_MAGCLOUD_MEMBER','APOGEE2_MAGCLOUD_CANDIDATE',
              'APOGEE2_RRLYR','APOGEE2_BULGE_RC','APOGEE2_SGR_DSPH','APOGEE2_APOKASC_GIANT',
              'APOGEE2_APOKASC_DWARF','APOGEE2_FAINT_EXTRA','APOGEE2_APOKASC',''])

class Apogee2Target2(BitMask) :
    '''
    BitMask class for APOGEE2_TARGET2
    '''
    name = ([ 'APOGEE2_K2_GAP','APOGEE2_CCLOUD_AS4','APOGEE2_STANDARD_STAR','APOGEE2_RV_STANDARD',
              'APOGEE2_SKY','APOGEE2_EXTERNAL_CALIB','APOGEE2_INTERNAL_CALIB','APOGEE2_DISK_SUBSTRUCTURE_MEMBER',
              'APOGEE2_DISK_SUBSTRUCTURE_CANDIDATE','APOGEE2_TELLURIC','APOGEE2_CALIB_CLUSTER','APOGEE2_K2_PLANETHOST',
              'APOGEE2_TIDAL_BINARY','APOGEE2_LITERATURE_CALIB','APOGEE2_GES_OVERLAP','APOGEE2_ARGOS_OVERLAP',
              'APOGEE2_GAIA_OVERLAP','APOGEE2_GALAH_OVERLAP','APOGEE2_RAVE_OVERLAP','APOGEE2_COMMIS_SOUTH_SPEC',
              'APOGEE2_HALO_MEMBER','APOGEE2_HALO_CANDIDATE','APOGEE2_1M_TARGET','APOGEE2_MOD_BRIGHT_LIMIT',
              'APOGEE2_CIS','APOGEE2_CNTAC','APOGEE2_EXTERNAL','APOGEE2_CVZ_AS4_OBAF',
              'APOGEE2_CVZ_AS4_GI','APOGEE2_CVZ_AS4_CTL','APOGEE2_CVZ_AS4_GIANT',''])

class Apogee2Target3(BitMask) :
    '''
    BitMask class for APOGEE2_TARGET3
    '''
    name = ([ 'APOGEE2_KOI','APOGEE2_EB','APOGEE2_KOI_CONTROL','APOGEE2_MDWARF',
              'APOGEE2_SUBSTELLAR_COMPANIONS','APOGEE2_YOUNG_CLUSTER','APOGEE2_K2','APOGEE2_OBJECT',
              'APOGEE2_ANCILLARY','APOGEE2_MASSIVE_STAR','APOGEE2_QSO','APOGEE2_CEPHEID',
              'APOGEE2_LOW_AV_WINDOWS','APOGEE2_BE_STAR','APOGEE2_YOUNG_MOVING_GROUP','APOGEE2_NGC6791',
              'APOGEE2_LABEL_STAR','APOGEE2_FAINT_KEPLER_GIANTS','APOGEE2_W345','APOGEE2_MASSIVE_EVOLVED',
              'APOGEE2_REDDENING_TARGETS','APOGEE2_KEPLER_MDWARF_KOI','APOGEE2_AGB','APOGEE2_M33',
              'APOGEE2_ULTRACOOL','APOGEE2_DISTANT_SEGUE_GIANTS','APOGEE2_CEPHEID_MAPPING','APOGEE2_SA57',
              'APOGEE2_K2_MDWARF','APOGEE2_RVVAR','APOGEE2_M31','APOGEE2_'])

class Apogee2Target4(BitMask) :
    '''
    BitMask class for APOGEE2_TARGET4
    '''
    name = ([ '','','','',
              '','','','',
              '','','','',
              '','','','',
              '','','','',
              '','','','',
              '','','','',
              '','','',''])

class ApogeeTarget1(BitMask) :
    '''
    BitMask class for APOGEE_TARGET1
    '''
    name = ([ 'APOGEE_FAINT','APOGEE_MEDIUM','APOGEE_BRIGHT','APOGEE_IRAC_DERED',
              'APOGEE_WISE_DERED','APOGEE_SFD_DERED','APOGEE_NO_DERED','APOGEE_WASH_GIANT',
              'APOGEE_WASH_DWARF','APOGEE_SCI_CLUSTER','APOGEE_EXTENDED','APOGEE_SHORT',
              'APOGEE_INTERMEDIATE','APOGEE_LONG','APOGEE_DO_NOT_OBSERVE','APOGEE_SERENDIPITOUS',
              'APOGEE_FIRST_LIGHT','APOGEE_ANCILLARY','APOGEE_M31_CLUSTER','APOGEE_MDWARF',
              'APOGEE_HIRES','APOGEE_OLD_STAR','APOGEE_DISK_RED_GIANT','APOGEE_KEPLER_EB',
              'APOGEE_GC_PAL1','APOGEE_MASSIVE_STAR','APOGEE_SGR_DSPH','APOGEE_KEPLER_SEISMO',
              'APOGEE_KEPLER_HOST','APOGEE_FAINT_EXTRA','APOGEE_SEGUE_OVERLAP',''])

class ApogeeTarget2(BitMask) :
    '''
    BitMask class for APOGEE_TARGET2
    '''
    name = ([ 'LIGHT_TRAP','APOGEE_FLUX_STANDARD','APOGEE_STANDARD_STAR','APOGEE_RV_STANDARD',
              'APOGEE_SKY','APOGEE_SKY_BAD','APOGEE_GUIDE_STAR','APOGEE_BUNDLE_HOLE',
              'APOGEE_TELLURIC_BAD','APOGEE_TELLURIC','APOGEE_CALIB_CLUSTER','APOGEE_BULGE_GIANT',
              'APOGEE_BULGE_SUPER_GIANT','APOGEE_EMBEDDEDCLUSTER_STAR','APOGEE_LONGBAR','APOGEE_EMISSION_STAR',
              'APOGEE_KEPLER_COOLDWARF','APOGEE_MIRCLUSTER_STAR','APOGEE_RV_MONITOR_IC348','APOGEE_RV_MONITOR_KEPLER',
              'APOGEE_GES_CALIBRATE','APOGEE_BULGE_RV_VERIFY','APOGEE_1MTARGET','APOGEE_',
              'APOGEE_','APOGEE_','APOGEE_','APOGEE_',
              'APOGEE_','APOGEE_','APOGEE_',''])

class ApogeeTarget3(BitMask) :
    '''
    BitMask class for APOGEE_TARGET3
    '''
    name = ([ '','','','','','','','',
              '','','','','','','','',
              '','','','','','','','',
              '','','','','','','',''])

class ApogeeTarget4(BitMask) :
    '''
    BitMask class for APOGEE_TARGET4
    '''
    name = ([ '','','','','','','','',
              '','','','','','','','',
              '','','','','','','','',
              '','','','','','','',''])

class RVBitMask(BitMask) :
    '''
    BitMask class for RVBitMask
    '''
    name = ([ 'RV_BCFIT','RV_BCFIT_FAIL','RV_FAINT_FIT','RV_WINDOW_MASK','RV_VALUE_ERROR','RV_RUNTIME_ERROR','RV_ERROR','',
              'NO_GOOD_VISITS','ALL_VISITS_REJECTED','RV_REJECT','RV_SUSPECT','','','','',
              '','','','','','','','',
              '','','','','','','','RESERVED'])
    descrip=([
     'Initial fit on BC combined spectra, then small RV range',
     'Failed fit on BC combined spectra',
     'Faint star, RV fit on reduced range',
     'Regions masked in RV fit'
     'Jointfit failed with ValueError',
     'Jointfit failed with RuntimeError',
     'Jointfit failed with exception',
     '',
     'No good visits for RV',
     'All visits rejected',
     'RV rejected based on fit vs xcorr RV',
     'RV suspect based on fit vs xcorr RV',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     '',
     ''
    ])

class MembersBitMask(BitMask) :
    ''' Membership '''
    name=['M92','M15','M53','N5466','N4147',
        'M2','M13','M3','M5','M12','M107',
        'M71','N2243','Be29', 'N2158','M35','N2420',
        'N188','M67','N7789','Pleiades','N6819',
        'ComaBer','N6791',
        'N5053','M68','N6397','M55','N5634','M22','M79','N3201','M10',
        'N6752','Omegacen','M54','N6229','Pal5','N6544','N6522','N288','N362','N1851',
        'M4','N2808','Pal6','47TUC','Pal1','N6539','N6388','N6441','N6316',
        'N6760','N6553','N6528',
        'DRACO','URMINOR','BOOTES1','SEXTANS','FORNAX','SCULPTOR','CARINA','','RESERVED']
    descrip = []
    for n in name: 
        descrip.append('Likely member of '+n)


def targflags(targ1,targ2,targ3,targ4,survey='apogee2') :

    if 'apogee2' in survey :
        mask1=Apogee2Target1()
        mask2=Apogee2Target2()
        mask3=Apogee2Target3()
        mask4=Apogee2Target4()
        return ','.join([mask1.getname(targ1),mask2.getname(targ2),mask3.getname(targ3),mask4.getname(targ4)]).strip(',')
    else :
        mask1=ApogeeTarget1()
        mask2=ApogeeTarget2()
        return ','.join([mask1.getname(targ1),mask2.getname(targ2)]).strip(',')
