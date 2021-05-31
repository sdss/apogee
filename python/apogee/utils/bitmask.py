"""
Routines for handling ASPCAP bitmasks
"""

import numpy as np
import pdb
import sys

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

    def print(self,fmt='txt',fp=sys.stdout) :
        """ Formatted output of bit definitions
        """
        if fmt == 'txt' :
            fp.write('{:25s}{:>6s}  {:s}'.format('Name','Bit','Description'))
        elif fmt == 'wiki' :
            fp.write('||{:25s}||{:>6s}||{:s}||'.format('Name','Bit','Description'))
        elif fmt == 'latex' :
            fp.write('{:25s}&{:>6s}&{:s}\\\\'.format('Name','Bit','Description'))
        elif fmt == 'par' :
            fp.write('masktype {:s} {:d}\n'.format(self.flagname,len(self.name)))
        elif fmt == 'html' or fmt == 'sdsshtml' :
            if fmt == 'html' :
                fp.write("[SDSS_GROUP TITLE='<h2 id={:s}>{:s}']\n".format(self.shorttitle,self.title))
                fp.write("{:s}\n".format(self.blurb))
            else :
                fp.write('<div id="{:s}" class="panel panel-default">\n'.format(self.flagname))
                fp.write('<div class="panel-heading">\n')
                fp.write('<h3 class="panel-title"><a class="accordion-toggle" href="#collapse{:s}" data-toggle="collapse" data-parent="#accordion-bitmask">{:s}&nbsp; </a></h3>\n'.format(self.shorttitle,self.title))
                fp.write('</div>\n')
                fp.write('<div id="collapse{:s}" class="panel-collapse collapse">\n'.format(self.flagname))
                fp.write('<div class="panel-body">\n')

            fp.write('<table class="table table-bordered table-condensed"\n')
            fp.write('<thead>\n')
            fp.write('<tr><th style="white-space:nowrap;">Bit&nbsp;Name</th><th style="white-space:nowrap;">Binary&nbsp;Digit</th><th>Description</th></tr>\n')
            fp.write('</thead>\n')
            fp.write('<tbody>\n')
        for ibit,name in enumerate(self.name) :
            if (name != 'RESERVED' and name != '' ) :
                if fmt == 'txt' :
                    fp.write('{:25s}{:6d}  {:s}'.format(name,ibit,self.descrip[ibit]))
                elif fmt == 'wiki' :
                    fp.write('||{:25s}||{:6d}||{:s}||'.format(name,ibit,self.descrip[ibit]))
                elif fmt == 'latex' :
                    fp.write('{:25s}&{:6d}&{:s}\\\\'.format(name,ibit,self.descrip[ibit]))
                elif fmt == 'par' :
                    fp.write('maskbits {:s} {:d} {:s} "{:s}"\n'.format(self.flagname,ibit,name,self.descrip[ibit]))
                elif fmt == 'html' :
                    fp.write('<tr><td style="white-space:nowrap;">{:s}<td>{:d}<td>{:s}\n'.format(name,ibit,self.descrip[ibit]))
        if fmt == 'html' or fmt == 'sdsshtml' :
            fp.write('</tbody>\n')
            fp.write('</table>\n')
            if fmt == 'html' :
                fp.write("[/SDSS_GROUP]\n")
            else :
                fp.write("</div>\n")
                fp.write("</div>\n")
                fp.write("</div>\n")
class StarBitMask(BitMask):
    '''
    BitMask class for APOGEE star bitmask (APOGEE_STARFLAG)
    '''

    flagname='APOGEE_STARFLAG'
    shorttitle='StarBitMask'
    title='APOGEE_STARFLAG, APOGEE_ANDFLAG : APOGEE star level bitmask '
    blurb='This bitmask is used to provide information and identify issues associated with <a href=”/dr17/irspec/apred/”>spectral processing</a>, radial velocity measurement, and <a href=”/dr17/irspec/spectral_combination/”>spectral combination</a>. At the visit level, it conveys information relevant to a specific spectrum. It is also used for combined spectra in both STARFLAG and ANDFLAG, where the former is a bitwise OR of the visit STARFLAG for all of the visits, and the latter is a bitwise AND of the visit STARFLAG for all of the visits. The bitmask <code>RV_FLAG</code> provides more details relevant to the RV determination.'
    name=(['BAD_PIXELS','COMMISSIONING','BRIGHT_NEIGHBOR','VERY_BRIGHT_NEIGHBOR','LOW_SNR','','','',
          '','PERSIST_HIGH','PERSIST_MED','PERSIST_LOW','PERSIST_JUMP_POS','PERSIST_JUMP_NEG','','',
          'SUSPECT_RV_COMBINATION','SUSPECT_BROAD_LINES','BAD_RV_COMBINATION','RV_REJECT','RV_SUSPECT','MULTIPLE_SUSPECT','RV_FAIL','SUSPECT_ROTATION',
          'MTPFLUX_LT_75','MTPFLUX_LT_50','','','','','','RESERVED'])
    level=([1,0,0,1,0,0,0,0,
             0,0,0,0,0,0,0,0,
             0,0,1,0,0,0,1,0,
             0,0,0,0,0,0,0,0])

    descrip=([
         'Spectrum has many bad pixels (>20%):  BAD',                                                         
         'Commissioning data (MJD<55761), non-standard configuration, poor LSF: WARN',                       
         'Star has neighbor more than 10 times brighter: WARN',
         'Star has neighbor more than 100 times brighter: possibly BAD',
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
         'Autocorrelation peak width large',
         'RVs from synthetic template differ very significatly (~10 km/s) from those from combined template: potentially BAD',
         'Rejected visit because cross-correlation RV differs significantly from least squares RV',
         'Suspect visit (but used!) because cross-correlation RV differs slightly from least squares RV',
         'Suspect multiple components from Gaussian decomposition of cross-correlation',
         'RV failure',
         'Suspect rotation: cross-correlation peak with template significantly broader than autocorretion of template',
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

    flagname='APOGEE_ASPCAPFLAG'
    shorttitle='AspcapBitMask'
    title='APOGEE_ASPCAPFLAG : ASPCAP star level bitmask '
    blurb='This bitmask is used to provide information and indicate possible issues associated with the <a href=”/dr17/irspec/aspcap/”>ASPCAP fits</a>. The flags refer to issues in the spectra impacting the fit (inherited from <code>STARFLAG</code>),  issues incurred during the fitting process, and issues incurred when evaluating the parameters. For flags associated with a stellar parameter or chemical abundance measurement, more detailed information is available in the <code>PARAMFLAG</code> or <code>ELEMFLAG</code> bitmask for that measurement. '

    name=(['TEFF_WARN','LOGG_WARN','VMICRO_WARN','M_H_WARN','ALPHA_M_WARN','C_M_WARN','N_M_WARN','STAR_WARN',
          'CHI2_WARN','COLORTE_WARN','ROTATION_WARN','SN_WARN','SPEC_HOLE_WARN','ATMOS_HOLE_WARN','VSINI_WARN','',
          'TEFF_BAD','LOGG_BAD','VMICRO_BAD','M_H_BAD','ALPHA_M_BAD','C_M_BAD','N_M_BAD','STAR_BAD',
          'CHI2_BAD','COLORTE_BAD','ROTATION_BAD','SN_BAD','SPEC_HOLE_BAD','ATMOS_HOLE_BAD','VSINI_BAD','NO_ASPCAP_RESULT',
          'MISSING_APSTAR','NO_GRID','BAD_FRAC_LOWSNR','BAD_FRAC_BADPIX','FERRE_FAIL','','','',
          'PROBLEM_TARGET','MULTIPLE_SUSPECT','','','','','','',
          '','','','','','','','',
          '','','','','','','','RESERVED'])
    level=([2,2,0,0,0,0,0,2,
            2,2,2,2,2,2,0,0,
            1,1,0,0,0,0,0,1,
            0,1,1,0,1,2,0,1,
            1,1,1,1,1,0,0,0,
            1,1,0,0,0,0,0,0,
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
         'higher than typical chi^2 (> 30*(SNR/100)**2)',
         'effective temperature more than 500K from photometric temperature for dereddened color (WARN)',
         'Spectrum has broad lines, with possible bad effects: FWHM of cross-correlation of spectrum with best RV template relative to auto-correltion of template > 1.5 (WARN)',
         'S/N<70 (WARN)',
         'Grid point within 2 grid steps of hole-filled synthesis ',
         'Grid point within 2 grid steps of hole-filled atmosphere ',
         ' ',
         ' ',
         'potentially BAD effective temperature (see PARAMFLAG[0] for details) ',
         'potentially BAD log g (see PARAMFLAG[1] for details) ',
         'potentially BAD vmicro (see PARAMFLAG[2] for details) ',
         'potentially BAD metals (see PARAMFLAG[3] for details) ',
         'potentially BAD [alpha/M] (see PARAMFLAG[4] for details) ',
         'potentially BAD [C/M] (see PARAMFLAG[5] for details) ',
         'potentially BAD [N/M] (see PARAMFLAG[6] for details) ',
         'BAD overall for star: set if any of TEFF, LOGG, CHI2, COLORTE, ROTATION, SN error are set, or any GRIDEDGE_BAD ',
         'significantly higher than typical chi^2 (> 50*(SNR/100)**2)',
         'effective temperature more than 1000K from photometric temperature for dereddened color',
         'Spectrum has broad lines, with possible bad effects: FWHM of cross-correlation of spectrum with best RV template relative to auto-correltion of template > 2 (BAD)',
         'S/N<50',
         'Grid point within 1 grid steps of hole-filled synthesis ',
         'Grid point within 1 grid steps of hole-filled atmosphere ',
         ' ',' ',
         'Missing apStar file','Not processed by any ASPCAP grid','Fraction low SNR pixels > 0.5',
         'Fraction bad pixels > 0.5 or 0.33 in any chip','FERRE failure (bad input?)','','','',
         'Target probably not suitable for standard ASPCAP analysis: extended, embedded, galaxy, etc',
         'Likely composite spectrum, possibly not suitable for standard ASPCAP analysis','','','','','','',
         '','','','','','','','',
         '','','','','','','',''
         ])

class ParamBitMask(BitMask):
    '''
    BitMask class for APOGEE ASPCAP bitmask (APOGEE_ASPCAPFLAG)
    '''

    flagname='APOGEE_PARAMFLAG'
    shorttitle='ParamBitMask'
    title='APOGEE_PARAMFLAG, APOGEE_ELEMFLAG : ASPCAP bitmask for individual parameters/abundances'
    blurb='These bitmasks are used to provide information indicate possible issues associated with individual measurements from <a href=”/dr17/irspec/aspcap/”>ASPCAP fits</a>. This bit provides more context for some of the bits in <code>ASPCAPFLAG</code>. A <code>PARAMFLAG</code> or <code>ELEMFLAG</code> is produced for each of the stellar parameters and chemical abundances measured in DR17, and each of the bitmasks has the same format. '
    name =['GRIDEDGE_BAD','CALRANGE_BAD','OTHER_BAD','FERRE_FAIL','PARAM_MISMATCH_BAD','FERRE_ERR_USED','TEFF_CUT','',
           'GRIDEDGE_WARN','CALRANGE_WARN','OTHER_WARN','FERRE_WARN','PARAM_MISMATCH_WARN','OPTICAL_WARN','ERR_WARN','FAINT_WARN',
           'PARAM_FIXED','RV_WARN','','','','','','',
           'SPEC_RC','SPEC_RGB','LOGG_CAL_MS','LOGG_CAL_RGB_MS','','','','RESERVED']

    level=[1,0,1,1,0,0,1,0,
           0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,0]

    descrip=([
         'Parameter within 1/8 grid spacing of grid edge : true value may be off grid ',
         'Parameter outside valid range of calibration determination ',
         'Other error condition ',
         'Failed solution in FERRE ',
         'Elemental abundance from window differs significantly from parameter abundance ',
         'FERRE uncertainty used (larger than parametric uncertainty) ',
         'Star in region of parameter space where abundances do not appear valid for this element ',
         ' ',
         'Parameter within 1 grid spacing of grid edge (not necessarily bad) ',
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
         'Spectroscopically identified RC star ',
         'Spectroscopically identified RGB star ',
         'DR16: use MS gravity calibration ',
         'DR16: use RBG/MS transition gravity calibration ',
         ' ',
         ' ',
         ' ',
         ' '
         ])

class PixelBitMask(BitMask) :
    '''
    BitMask class for APOGEE pixel bitmask (APOGEE_PIXMASK)
    '''
    flagname='APOGEE_PIXMASK'
    shorttitle='PixelBitMask'
    title='APOGEE_PIXMASK : APOGEE bitmask for individual pixels in a spectrum'
    blurb='This bitmask is used to provide information associated with individual pixels in a one-dimensional spectrum. At the visit level, <code>PIXMASK</code> refers to the individual spectrum. At the combined level, PIXMASK attempts to appropriately combine the PIXMASKs of the visit spectrum level. '
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

class ExtratargBitMask(BitMask) :
    '''
    BitMask class for EXTRATARG
    '''
    flagname='APOGEE2_EXTRATARG'
    shorttitle='ExtratargBitMask'
    title='EXTRATARG : basic targeting information'
    blurb='This bitmask is used to indicate some basic useful targeting information. Some bits are used to tie together similar samples that are indicated through distinct APOGEE-1 and APOGEE-2 targeting bitmasks. Another bit is used to indicate targets that have “duplicate” observations (the same <code>APOGEE_ID</code> observed in distinct <code>FIELD</code>).  '
    name=(['NOT_MAIN','COMMISSIONING','TELLURIC','APO1M','DUPLICATE','','',''])
    descrip=([
     'Not a main sample target',
     'Commissioning observation',
     'Targeted as telluric',
     'APO/NMSU 1M observation',
     'Non-primary (not highest S/N) duplicate',
     '','',''])

class Apogee2Target1(BitMask) :
    '''
    BitMask class for APOGEE2_TARGET1
    '''

    flagname='APOGEE2_TARGET1'
    shorttitle='Apogee2Target1'
    title='APOGEE2_TARGET1 : APOGEE2 targeting bitmask (1 of 3)'
    blurb='APOGEE-2 targeting information is conveyed in three bitmasks: <code>APOGEE2_TARGET1</code>, <code>APOGEE2_TARGET2</code>, and <code>APOGEE2_TARGET3</code> (<code>APOGEE2_TARGET4</code> was added to the data model, but is not used). The bitmasks provide context for the reason that a target was selected for observations on a plate. Descriptions can be found on <a href="/dr17/irspec/targets/">Targeting Information</a>.  '
    name = ([ 'APOGEE2_ONEBIN_GT_0_5','APOGEE2_TWOBIN_0_5_TO_0_8','APOGEE2_TWOBIN_GT_0_8','APOGEE2_IRAC_DERED',
              'APOGEE2_WISE_DERED','APOGEE2_SFD_DERED','APOGEE2_NO_DERED','APOGEE2_WASH_GIANT',
              'APOGEE2_WASH_DWARF','APOGEE2_SCI_CLUSTER','APOGEE2_CLUSTER_CANDIDATE','APOGEE2_SHORT',
              'APOGEE2_MEDIUM','APOGEE2_LONG','APOGEE2_NORMAL_SAMPLE','APOGEE2_MANGA_LED',
              'APOGEE2_ONEBIN_GT_0_3','APOGEE2_WASH_NOCLASS','APOGEE2_STREAM_MEMBER','APOGEE2_STREAM_CANDIDATE',
              'APOGEE2_DSPH_MEMBER','APOGEE2_DSPH_CANDIDATE','APOGEE2_MAGCLOUD_MEMBER','APOGEE2_MAGCLOUD_CANDIDATE',
              'APOGEE2_RRLYR','APOGEE2_BULGE_RC','APOGEE2_SGR_DSPH','APOGEE2_APOKASC_GIANT',
              'APOGEE2_APOKASC_DWARF','APOGEE2_FAINT_EXTRA','APOGEE2_APOKASC',''])

    descrip=([
     'Selected in single (J-Ks)o > 0.5 color bin',
     'Selected in "blue" 0.5 < (J-Ks)o < 0.8 color bin',
     'Selected in "red" (J-Ks)o > 0.8 color bin',
     'Selected with RJCE-IRAC dereddening',
     'Selected with RJCE-WISE dereddening',
     'Selected with SFD_EBV dereddening',
     'Selected with no dereddening',
     'Selected as Wash+DDO51 photometric giant',
     'Selected as Wash+DDO51 photometric dwarf',
     'Science cluster candidate member',
     'Selected as Globular Cluster candidate',
     'Selected as part of a short cohort',
     'Selected as part of a medium cohort',
     'Selected as part of a long cohort',
     'Selected as part of the random sample',
     'Star on a shared MaNGA-led design',
     'Selected in single (J-Ks)o > 0.3 color bin',
     'Selected because it has no W+D classification',
     'Selected as confirmed halo tidal stream member',
     'Selected as potential halo tidal stream member (based on photometry)',
     'Selected as confirmed dSph member (non Sgr)',
     'Selected as potential dSph member (non Sgr) (based on photometry)',
     'Selected as confirmed Mag Cloud member',
     'Selected as potential Mag Cloud member (based on photometry)',
     'Selected as an RR Lyrae star',
     'Selected as a bulge candidate RC star',
     'Selected as confirmed Sgr core/stream member',
     'Selected as part of APOKASC "giant" sample',
     'Selected as part of APOKASC "dwarf" sample',
     '"Faint" star (fainter than cohort limit; not required to reach survey S/N requirement)',
     'Selected as part of the APOKASC program (incl. seismic/gyro targets and others, both the Cygnus field and K2)'
    ])


class Apogee2Target2(BitMask) :
    '''
    BitMask class for APOGEE2_TARGET2
    '''
    flagname='APOGEE2_TARGET2'
    shorttitle='Apogee2Target2'
    title='APOGEE2_TARGET2 : APOGEE2 targeting bitmask (2 of 3)'
    blurb=''
    name = ([ 'APOGEE2_K2_GAP','APOGEE2_CCLOUD_AS4','APOGEE2_STANDARD_STAR','APOGEE2_RV_STANDARD',
              'APOGEE2_SKY','APOGEE2_EXTERNAL_CALIB','APOGEE2_INTERNAL_CALIB','APOGEE2_DISK_SUBSTRUCTURE_MEMBER',
              'APOGEE2_DISK_SUBSTRUCTURE_CANDIDATE','APOGEE2_TELLURIC','APOGEE2_CALIB_CLUSTER','APOGEE2_K2_PLANETHOST',
              'APOGEE2_TIDAL_BINARY','APOGEE2_LITERATURE_CALIB','APOGEE2_GES_OVERLAP','APOGEE2_ARGOS_OVERLAP',
              'APOGEE2_GAIA_OVERLAP','APOGEE2_GALAH_OVERLAP','APOGEE2_RAVE_OVERLAP','APOGEE2_COMMIS_SOUTH_SPEC',
              'APOGEE2_HALO_MEMBER','APOGEE2_HALO_CANDIDATE','APOGEE2_1M_TARGET','APOGEE2_MOD_BRIGHT_LIMIT',
              'APOGEE2_CIS','APOGEE2_CNTAC','APOGEE2_EXTERNAL','APOGEE2_CVZ_AS4_OBAF',
              'APOGEE2_CVZ_AS4_GI','APOGEE2_CVZ_AS4_CTL','APOGEE2_CVZ_AS4_GIANT',''])

    descrip = ([
     'K2 Galactic Archeology Program Star',
     'California Cloud target',
     'Stellar parameters/abundance standard',
     'Stellar RV standard',
     'Sky fiber',
     'External survey calibration target (generic flag; others below dedicated to specific surveys)',
     'Internal survey calibration target (observed in at least 2 of: APOGEE-1, -2N, -2S)',
     'Bright time extension: outer disk substructure (Triand, GASS, and A13) members',
     'Bright time extension: outer disk substructure (Triand, GASS, and A13) candidates',
     'Telluric calibrator target',
     'Selected as calibration cluster member',
     'Planet host in the K2 field',
     'Ancillary KOI Program (Simonian)',
     'Overlap with high-resolution literature studies',
     'Overlap with Gaia-ESO',
     'Overlap with ARGOS',
     'Overlap with Gaia',
     'Overlap with GALAH',
     'Overlap with RAVE',
     'Commissioning special targets for APOGEE2S',
     'Halo Member',
     'Halo Candidate',
     'Selected as a 1-m target',
     'Selected in a cohort with H>10 rather than H>7',
     'Carnegie program target',
     'Chilean community target',
     'Proprietary external target',
     'OBAF stars selected for multi-epoc observations Andrew T.',
     'Submitted program to be on CVZ plate (Known Planets, ATL, Tayar-Subgiant, Canas-Cool-dwarf)',
     'Filler CTL star selected from the TESS Input Catalog',
     'Filler Giant selected with RPMJ'
    ])

class Apogee2Target3(BitMask) :
    '''
    BitMask class for APOGEE2_TARGET3
    '''
    flagname='APOGEE2_TARGET3'
    shorttitle='Apogee2Target2'
    title='APOGEE2_TARGET3 : APOGEE2 targeting bitmask (3 of 3)'
    blurb=''
    name = ([ 'APOGEE2_KOI','APOGEE2_EB','APOGEE2_KOI_CONTROL','APOGEE2_MDWARF',
              'APOGEE2_SUBSTELLAR_COMPANIONS','APOGEE2_YOUNG_CLUSTER','APOGEE2_K2','APOGEE2_OBJECT',
              'APOGEE2_ANCILLARY','APOGEE2_MASSIVE_STAR','APOGEE2_QSO','APOGEE2_CEPHEID',
              'APOGEE2_LOW_AV_WINDOWS','APOGEE2_BE_STAR','APOGEE2_YOUNG_MOVING_GROUP','APOGEE2_NGC6791',
              'APOGEE2_LABEL_STAR','APOGEE2_FAINT_KEPLER_GIANTS','APOGEE2_W345','APOGEE2_MASSIVE_EVOLVED',
              'APOGEE2_REDDENING_TARGETS','APOGEE2_KEPLER_MDWARF_KOI','APOGEE2_AGB','APOGEE2_M33',
              'APOGEE2_ULTRACOOL','APOGEE2_DISTANT_SEGUE_GIANTS','APOGEE2_CEPHEID_MAPPING','APOGEE2_SA57',
              'APOGEE2_K2_MDWARF','APOGEE2_RVVAR','APOGEE2_M31',''])

    descrip = ([
     'Selected as part of the long cadence KOI study',
     'Selected as part of the EB program',
     'Selected as part of the long cadence KOI "control sample"',
     'Selected as part of the M dwarf study',
     'Selected as part of the substellar companion search',
     'Selected as part of the young cluster study (IN-SYNC)',
     'Selected as part of the K2 program (BTX and Main Survey)',
     'This object is an APOGEE-2 target',
     'Selected as an ancillary target',
     'Selected as part of the Massive Star program',
     'ancillary QSO pilot program (Albareti)',
     'ancillary Cepheid sparse targets (Beaton)',
     'ancillary Deep Disk sample (Bovy)',
     'ancillary ASHELS sample (Chojnowski)',
     'ancillary young moving group members (Downes)',
     'ancillary NGC 6791 star (Geisler)',
     'ancillary Cannon calibrator Sample (Ness)',
     'ancillary APOKASC faint giants (Pinsonneault)',
     'ancillary W3/4/5 star forming complex (Roman-Lopes)',
     'ancillary massive/evolved star targets (Stringfellow)',
     'ancillary extinction targets (Schlafly)',
     'ancillary M dwarf targets (Smith)',
     'ancillary AGB sample (Zamora)',
     'Ancillary M33 Program (Anguiano)',
     'Ancillary Ultracool Dwarfs Program (Burgasser)',
     'Ancillary Distant SEGUE Giants Program (Harding)',
     'Ancillary Cepheid Mapping Program (Inno)',
     'Ancillary SA57 Kapteyn Field Program (Majewski)',
     'Ancillary K2 M dwarf Program (Smith)',
     'Ancillary RV Variables Program (Troup)',
     'Ancillary M31 Program (Zasowski)'
    ])

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
    flagname='APOGEE_TARGET1'
    shorttitle='ApogeeTarget1'
    title='APOGEE_TARGET1 : APOGEE1 targeting bitmask (1 of 2)'
    blurb='APOGEE-1 targeting information is conveyed in two bitmasks, <code>APOGEE_TARGET1</code> and <code>APOGEE_TARGET2</code>. The bitmasks provide context for the reason that a target was selected for observations on a plate. Detailed descriptions can be found on <a href="/dr17/irspec/targets/">Targeting Information</a>.  '
    name = ([ 'APOGEE_FAINT','APOGEE_MEDIUM','APOGEE_BRIGHT','APOGEE_IRAC_DERED',
              'APOGEE_WISE_DERED','APOGEE_SFD_DERED','APOGEE_NO_DERED','APOGEE_WASH_GIANT',
              'APOGEE_WASH_DWARF','APOGEE_SCI_CLUSTER','APOGEE_EXTENDED','APOGEE_SHORT',
              'APOGEE_INTERMEDIATE','APOGEE_LONG','APOGEE_DO_NOT_OBSERVE','APOGEE_SERENDIPITOUS',
              'APOGEE_FIRST_LIGHT','APOGEE_ANCILLARY','APOGEE_M31_CLUSTER','APOGEE_MDWARF',
              'APOGEE_HIRES','APOGEE_OLD_STAR','APOGEE_DISK_RED_GIANT','APOGEE_KEPLER_EB',
              'APOGEE_GC_PAL1','APOGEE_MASSIVE_STAR','APOGEE_SGR_DSPH','APOGEE_KEPLER_SEISMO',
              'APOGEE_KEPLER_HOST','APOGEE_FAINT_EXTRA','APOGEE_SEGUE_OVERLAP',''])

    descrip = ([
     'Star selected in faint bin1 of its cohort',
     'Star selected in medium bin1 of its cohort',
     'Star selected in bright bin1 of its cohort',
     'Selected w/ RJCE-IRAC dereddening',
     'Selected w/ RJCE-WISE dereddening',
     'Selected w/ SFD dereddening',
     'Selected w/ no dereddening',
     'Selected as giant using Washington photometry',
     'Selected as dwarf using Washington photometry',
     'Selected as probable cluster member',
     'Extended object',
     'Selected as "short" (~3-visit) cohort target (includes 1-visit samples)',
     'Selected as "intermediate" cohort (~6-visit) target',
     'Selected as "long" cohort (~12- or 24-visit) target',
     'Do not observe (again) -- undesired dwarf, galaxy, etc',
     'Serendipitous interesting target to be re-observed',
     'First Light plate target',
     'Ancillary target',
     'M31 Clusters (Allende Prieto, Schiavon, Bizyaev, OConnell, Shetrone)',
     'RVs of M Dwarfs (Blake, Mahadevan, Hearty, Deshpande, Nidever, Bender, Crepp, Carlberg, Terrien, Schneider) -- both original list and second-round extension',
     'Stars with Optical Hi-Res Spectra (Fabbian, Allende Prieto, Smith, Cunha)',
     'Oldest Stars in Galaxy (Harding, Johnson)',
     'Ages/Compositions? of Disk Red Giants (Johnson, Epstein, Pinsonneault, Lai, Bird, Schonrich, Chiappini)',
     'Kepler EBs (Mahadevan, Fleming, Bender, Deshpande, Hearty, Nidever, Terrien)',
     'Globular Cluster Pops in the MW (Simmerer, Ivans, Shetrone)',
     'Massive Stars in the MW (Herrero, Garcia-Garcia, Ramirez-Alegria)',
     'Sgr (dSph) member',
     'Kepler asteroseismology program target (Epstein)',
     'planet-host program target (Epstein)',
     'as "faint" target in low-target-density fields',
     'SEGUE overlap',
    ])


class ApogeeTarget2(BitMask) :
    '''
    BitMask class for APOGEE_TARGET2
    '''
    flagname='APOGEE_TARGET2'
    shorttitle='ApogeeTarget2'
    title='APOGEE_TARGET2 : APOGEE1 targeting bitmask (2 of 2)'
    blurb=''
    name = ([ 'LIGHT_TRAP','APOGEE_FLUX_STANDARD','APOGEE_STANDARD_STAR','APOGEE_RV_STANDARD',
              'APOGEE_SKY','APOGEE_SKY_BAD','APOGEE_GUIDE_STAR','APOGEE_BUNDLE_HOLE',
              'APOGEE_TELLURIC_BAD','APOGEE_TELLURIC','APOGEE_CALIB_CLUSTER','APOGEE_BULGE_GIANT',
              'APOGEE_BULGE_SUPER_GIANT','APOGEE_EMBEDDEDCLUSTER_STAR','APOGEE_LONGBAR','APOGEE_EMISSION_STAR',
              'APOGEE_KEPLER_COOLDWARF','APOGEE_MIRCLUSTER_STAR','APOGEE_RV_MONITOR_IC348','APOGEE_RV_MONITOR_KEPLER',
              'APOGEE_GES_CALIBRATE','APOGEE_BULGE_RV_VERIFY','APOGEE_1MTARGET','',
              '','','','',
              '','','',''])
    descrip = ([
     'Light trap',
     'Flux standard',
     'Stellar abundance/parameters standard',
     'RV standard',
     'Sky',
     'Selected as sky but IDed as bad (via visual examination or observation)',
     'Guide star',
     'Bundle hole',
     'Selected as telluric std but IDed as too red (via SIMBAD or observation)',
     'Hot (telluric) standard',
     'Known calibration cluster member',
     'Selected as probable giant in bulge',
     'Selected as probable supergiant in bulge',
     'Young Nebulous Clusters (Covey, Tan)',
     'Milky Way Long Bar (Zasowski)',
     'Be Emission Line Stars (Chojnowski, Whelan)',
     'Kepler Cool Dwarfs (van Saders)',
     'Outer Disk MIR Clusters (Beaton)',
     'RV Variability in IC348 (Nidever, Covey(?))',
     'RV Variability for Kepler Planet Hosts and Binaries (Deshpande, Fleming, Mahadevan)',
     'GAIA-ESO calibration targets',
     'RV Verification (Nidever)',
     'Selected as a 1-m target (Holtzman)',
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
    flagname='APOGEE_RV_FLAG'
    shorttitle='RVBitMask'
    title='APOGEE_RV_FLAG : Bitmask for radial velocity information'
    blurb='This bitmask is used to provide information ms associated with radial velocity measurements.  The bitmask provides more context to RV-related bits in <code>STARFLAG</code>. This bitmask appears at the visit- and combined- level; some bits are set at the visit level (e.g., bits 10, 11) and others at the combined level. '
    name = ([ 'RV_BCFIT','RV_BCFIT_FAIL','RV_FAINT_FIT','RV_WINDOW_MASK','RV_VALUE_ERROR','RV_RUNTIME_ERROR','RV_ERROR','',
              'NO_GOOD_VISITS','ALL_VISITS_REJECTED','RV_REJECT','RV_SUSPECT','','','','',
              '','','','','','','','',
              '','','','','','','','RESERVED'])
    descrip=([
     'Initial fit on BC combined spectra, then small RV range',
     'Failed fit on BC combined spectra',
     'Faint star, RV fit on reduced range',
     'Regions masked in RV fit',
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
    flagname='APOGEE_MEMBERFLAG'
    shorttitle='MembersBitMask'
    title='APOGEE_MEMBERFLAG : Bitmask to identify likely members of clusters/dwarf Spheroidals'
    blurb='This bitmask is produced to indicate the likely membership of a given target in intentionally targeted star clusters and dwarf galaxies. Each bit refers to a specific star cluster or dwarf galaxy. A bit is set if a given target meets a set of membership criteria, including its sky position, mean radial velocity, and proper motions. <code>MEMBER</code> is convenient but may not be applicable for all scientific applications. '

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
        return ','.join([mask1.getname(targ1),mask2.getname(targ2),mask3.getname(targ3),mask4.getname(targ4)]).strip(',').replace(',,',',')
    else :
        mask1=ApogeeTarget1()
        mask2=ApogeeTarget2()
        return ','.join([mask1.getname(targ1),mask2.getname(targ2)]).strip(',').replace(',,',',')

def print_bitmasks(fmt='html',out=None) :
  
    if out is not None : fp = open(out,'w') 
    else : fp  = sys.stdout

    for mask in [ AspcapBitMask(), ParamBitMask(), StarBitMask(), RVBitMask(), PixelBitMask(), ExtratargBitMask(),  
                  Apogee2Target1(), Apogee2Target2(), Apogee2Target3(), ApogeeTarget1(), ApogeeTarget2(), MembersBitMask()] :
        mask.print(fmt=fmt,fp=fp)
