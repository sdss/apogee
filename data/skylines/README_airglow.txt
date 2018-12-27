Description of the airglow linelists
------------------------------------
David Nidever
2/29/2012

-------
CRIRES
-------
There are some linelists from the CRIRES static calibration tables:
http://www.eso.org/observing/dfo/quality/CRIRES/pipeline/pipe_calib.html#other
It's not clear where these linelists originated from, but probably HITRAN
http://www.cfa.harvard.edu/hitran/

CR_GCAT_061130A_lines_oh.fits
CR_GCAT_061130A_lines_hitran.fits
These files are available in apogeereduce/lib/linelists/

The OH list was used as a starting point for the "original" APOGEE airglow
linelist (see below).

----------
Rousselot
----------
These are OH lines from Rousselot et al. (2000).  These were not used
for APOGEE linelists partly because of the low accuracy (~0.05A) of
the wavelengths.

rousselot_list_v2.0.dat
rousselot.vac.txt
These files are available in apogeereduce/lib/linelists/

---------------------------------
Original APOGEE airglow linelist
---------------------------------
This airglow linelist was made with
apogeereduce/lib/linelists/make_ohcrires_apogee_linelist.pro
and used for the original APOGEE reductions (through v0.92).
Started by David Nidever in March 2011.

-The linelist was started with the brightest lines in the
 CRIRES OH linelist.
-Close doublets are occupy a single line in the list with
 doublet=1, dbl_wsep=the wavelength separation of the two lines
 and wave=the wavelength midway betwee the two lines.  Most of
 the doublet information was used from the CRIRES linelist.
-For a handful of lines the CRIRES doublet wavelength separation
 was found to not work well with the APOGEE spectra and were refit
 using a two-Gaussian component (equal heights).
-A large number of weak lines (~100) that were missing from the
 linelist but seen in the APOGEE spectra were added by hand.
 The wavelengths of these lines are quite inaccurate.
-USELSF and USEWAVE columns were added to show which lines are good
 for determining the LSF and wavelength solution in the APOGEE
 pipeline.  This was vetted manually using trial and error.

The original linelist (without USELSF and USEWAVE) information
is called crires.oh.vac.apogee, but the final version was called
airglow.txt (the default pipeline linelist to use).  The "permanant"
name is now original_apogee_airglow_linelist.txt.  These are all
available in: apogeereduce/lib/linelists/


---------------------------------------
Dmitry Bizyaev APOGEE airglow Linelist
---------------------------------------
Dmitry Bizyaev created a new APOGEE airglow linelist using a large
stack of APOGEE sky spectra (in January 2012).

These files are in apogeereduce/lib/linelists/
dmitry_a_supersky.fits  red chip summed sky spectrum
dmitry_b_supersky.fits  green/middle chip summed sky spectrum
dmitry_c_supersky.fits  blue chip summed sky spectrum
dmitry_airglow_list_vac.txt  airglow linelist
dmitry_airglow_list_vac.txt.orig (original file from Dmitry)

The linelist was created by manually fitting Gaussians to lines in
IRAF and is fairly complete (~400 lines).  The linelist suffers from
wavelength errors because it was calibrated using the Rousselot list
which has wavelength errors (see above).


-------------------------------------
End-of-Night APOGEE Airglow linelist
-------------------------------------
It is not straightforward to create a new airglow linelist using
APOGEE data because all of the data are spectrally dithered by
half a pixel and lamp exposures are only taken in the afternoon
and end of night.  The APOGEE pipeline uses the airglow lines to
wavelength calibrate the science exposures, and therefore the
"normal" APOGEE science frames are not themselves suitable for
creating a new airglow linelist.  However, there are often no dithers
between the object exposure of the night and the first end-of-night
lamp exposures.

The 35 sky fibers from the last object frame of the night
(55968/ap1D-04060081) was used to create a new linelist (Feb 2012
by D.Nidever)  The lamps taken shortly afterwards (04060091+92) were
used for the wavelength calibration.  Lines are identified and fit
with a Gaussian. A second round of fitting does simultaneous fitting
of multiple Gaussians for lines that are close together.  The fitted
parameters are then averaged over the 35 sky fibers.

Lines with close neighbors (less than 1A) then they have doublet=1
and dbl_wsep=the wavelength separation between the lines.  However,
in contrast to the original and final APOGEE linelists, close
doublet take up two lines and wave=the wavelength of the individual
component lines.

Lines are flagged as being good for determining the LSF and wavelength
calibration based on the scatter in wavelength across the 35 fibers,
their height, width and if they have a close neighbor or are part of
a doublet.

This program makes the linelist:
make_airglow_linelist_endofnight.pro
apogeereduce/pro/lib/linelists/

The origingl linelist is called
airglow_linelist_endofnight.txt
apogeereduce/lib/linelists/

The lines were then manually checked and vetted for LSF/wavelength
determination by hand using trial and error.  This list is called
airglow_linelist_endofnight_vetted.txt

All lines are called OH even though some lines migth be from other
species.

This whole process can be redone better in the future if a sky flat is
taken at the very end of the night.


-------------------------------
Final APOGEE Airglow Linelist
-------------------------------
A final APOGEE airglow linelist was created using the best parts of the
linelists mentioned above (Feb 2012 by D.Nidever).

The endofnight lines are as a starting point.
-doublet information is taken from the original APOGEE airglow linelist
 (most originally coming from the CRIRES OH linelist)
-lines in the Dmitry linelist that are not in the endofnight list
 are added (the wavelengths of the entire Dmitry linelist were first
 calibrated to the endofnight list)
-some manual tweaking and deleting of lines
-USELSF and USEWAVE comes from the endofnight list

All lines are called OH even though some lines migth be from other
species.

The program that makes the linelist is called
combine_airglow_linelists.pro
apogeereduce/pro/lib/linelists/

and the final linelist is called
combine_airglow_linelists.txt
apogeereduce/lib/linelists/

The pipeline uses the file apogeereduce/lib/linelists/airglow.txt
which is a soft link to the most current linelist.  At the
writing of this document that is combine_airglow_linelists.txt.


------------------------
October 2018 revisit
-------------------------
Here are the results from 6 nights.

Here's the information that I'm using:
daynum = [2231, 2239, 2244, 2245, 2247, 2255]
ndaynum = n_elements(daynum)
plates = [9684, 9684, 9050, 9050, 8443, 8309]
fields = ['330+60', '330+60', '040+43_MGA', '040+43_MGA', '029+77_MGA', '102+61_MGA']

Here's the program I wrote that's calculating this:
/uufs/chpc.utah.edu/common/home/u0914350/wave/check_airglow_lines.pro


   ID  NLINES  WORIG       WMED      WDIFF   WSIG   WRMS    WERR    WDIFF2
    6  204  15187.1453  15187.1340  0.0113  0.0177  0.0241  0.0012  ------
    8  204  15240.9612  15240.9567  0.0045  0.0114  0.0222  0.0008  0.0078
   17  204  15287.7944  15287.7915  0.0029  0.0103  0.0217  0.0007 -0.0022
   24  204  15332.4069  15332.4044  0.0025  0.0156  0.0225  0.0011  0.0004
   42  204  15395.3348  15395.3346  0.0002  0.0097  0.0214  0.0007 -0.0027
   61  202  15500.8735  15500.8695  0.0040  0.0120  0.0224  0.0008  ------
   62  204  15509.7906  15509.7882  0.0024  0.0143  0.0226  0.0010  0.0103
   63  202  15517.7986  15517.8215 -0.0229  0.0168  0.0224  0.0012  ------
   71  204  15546.1479  15546.1417  0.0062  0.0114  0.0219  0.0008  0.0084
   86  204  15597.6376  15597.6318  0.0058  0.0128  0.0229  0.0009  0.0072
   94  204  15631.5221  15631.5155  0.0066  0.0136  0.0231  0.0009  0.0034
   99  204  15654.9717  15654.9614  0.0103  0.0176  0.0281  0.0012 -0.2109
  118  204  15702.5377  15702.5283  0.0094  0.0120  0.0233  0.0008  0.0086
  181  204  15972.6102  15972.6099  0.0003  0.0163  0.0235  0.0011  0.0031
  191  200  16030.8482  16030.8442  0.0040  0.0155  0.0233  0.0011  0.0035
  205  202  16079.7690  16079.7670  0.0020  0.0145  0.0211  0.0010 -0.0006
  239  202  16194.6264  16194.6272 -0.0008  0.0161  0.0215  0.0011 -0.0029
  270  204  16317.1734  16317.1635  0.0099  0.0145  0.0240  0.0010  0.0372
  275  204  16350.6348  16350.6238  0.0110  0.0167  0.0261  0.0012  0.0444
  285  204  16388.5017  16388.4885  0.0132  0.0141  0.0239  0.0010  0.0694
  291  204  16502.3666  16502.3632  0.0034  0.0121  0.0226  0.0008 -0.0068
  303  204  16553.8217  16553.8173  0.0044  0.0132  0.0235  0.0009  0.0195
  350  202  16840.4948  16840.5043 -0.0095  0.0179  0.0214  0.0013  ------
  363  204  16903.6824  16903.6819  0.0005  0.0149  0.0223  0.0010  0.1241

On Mon, Oct 8, 2018 at 4:51 PM David Nidever <dnidever@email.noao.edu> wrote:
I've added a final column that shows the difference between the new wavelengths and the airglow.new file.
I've bolded the 5 problematic ones

   ID  NLINES  WORIG       WMED      WDIFF   WSIG   WRMS    WERR  WDIFF2
    6   34  15187.1453  15187.1271  0.0182  0.0170  0.0191  0.0029  ------
    8   34  15240.9612  15240.9567  0.0045  0.0126  0.0170  0.0022  0.0078
   17   34  15287.7944  15287.7933  0.0011  0.0096  0.0158  0.0016 -0.0022
   24   34  15332.4069  15332.4058  0.0011  0.0137  0.0178  0.0023  0.0004
   42   34  15395.3348  15395.3354 -0.0006  0.0081  0.0166  0.0014 -0.0027
   61   33  15500.8735  15500.8693  0.0042  0.0102  0.0168  0.0018  ------
   62   34  15509.7906  15509.7859  0.0047  0.0112  0.0182  0.0019  0.0103
   63   34  15517.7986  15517.8124 -0.0138  0.0159  0.0180  0.0027  ------
   71   34  15546.1479  15546.1407  0.0072  0.0108  0.0176  0.0018  0.0084
   86   34  15597.6376  15597.6311  0.0065  0.0115  0.0180  0.0020  0.0072
   94   34  15631.5221  15631.5182  0.0039  0.0098  0.0181  0.0017  0.0034
   99   34  15654.9717  15654.9608  0.0109  0.0174  0.0212  0.0030 -0.2109
  118   34  15702.5377  15702.5273  0.0104  0.0100  0.0184  0.0017  0.0086
  181   34  15972.6102  15972.6099  0.0003  0.0167  0.0212  0.0029  0.0031
  191   34  16030.8482  16030.8439  0.0043  0.0170  0.0242  0.0029  0.0035
  205   34  16079.7690  16079.7653  0.0037  0.0163  0.0203  0.0028 -0.0006
  239   34  16194.6264  16194.6260  0.0004  0.0175  0.0208  0.0030 -0.0029
  270   34  16317.1734  16317.1621  0.0113  0.0138  0.0203  0.0024  0.0372
  275   34  16350.6348  16350.6220  0.0128  0.0135  0.0221  0.0023  0.0444
  285   34  16388.5017  16388.4891  0.0126  0.0148  0.0204  0.0025  0.0694
  291   34  16502.3666  16502.3617  0.0049  0.0108  0.0185  0.0018 -0.0068
  303   34  16553.8217  16553.8163  0.0054  0.0154  0.0200  0.0026  0.0195
  350   34  16840.4948  16840.5128 -0.0180  0.0191  0.0202  0.0033  ------
  363   34  16903.6824  16903.6802  0.0022  0.0179  0.0206  0.0031  0.1241


On Mon, Oct 8, 2018 at 4:19 PM David Nidever <dnidever@email.noao.edu> wrote:
Okay, here are the first results from the 1st night.  This compares the "new" median wavelengths across the 34/35 sky fibers.
WORIG are the wavelengths from airglow.txt. WDIFF is the difference between the new and old wavelengths.
This is for the 24 lines which are a combination of the 20 lines I identified and the ones in airglow.new.

   ID  NLINES   WORIG     WMED   WDIFF     WSIG   WRMS    WERR
    6   34  15187.1453  15187.1271  0.0182  0.0170  0.0191  0.0029
    8   34  15240.9612  15240.9567  0.0045  0.0126  0.0170  0.0022
   17   34  15287.7944  15287.7933  0.0011  0.0096  0.0158  0.0016
   24   34  15332.4069  15332.4058  0.0011  0.0137  0.0178  0.0023
   42   34  15395.3348  15395.3354 -0.0006  0.0081  0.0166  0.0014
   61   33  15500.8735  15500.8693  0.0042  0.0102  0.0168  0.0018
   62   34  15509.7906  15509.7859  0.0047  0.0112  0.0182  0.0019
   63   34  15517.7986  15517.8124 -0.0138  0.0159  0.0180  0.0027
   71   34  15546.1479  15546.1407  0.0072  0.0108  0.0176  0.0018
   86   34  15597.6376  15597.6311  0.0065  0.0115  0.0180  0.0020
   94   34  15631.5221  15631.5182  0.0039  0.0098  0.0181  0.0017
   99   34  15654.9717  15654.9608  0.0109  0.0174  0.0212  0.0030
  118   34  15702.5377  15702.5273  0.0104  0.0100  0.0184  0.0017
  181   34  15972.6102  15972.6099  0.0003  0.0167  0.0212  0.0029
  191   34  16030.8482  16030.8439  0.0043  0.0170  0.0242  0.0029
  205   34  16079.7690  16079.7653  0.0037  0.0163  0.0203  0.0028
  239   34  16194.6264  16194.6260  0.0004  0.0175  0.0208  0.0030
  270   34  16317.1734  16317.1621  0.0113  0.0138  0.0203  0.0024
  275   34  16350.6348  16350.6220  0.0128  0.0135  0.0221  0.0023
  285   34  16388.5017  16388.4891  0.0126  0.0148  0.0204  0.0025
  291   34  16502.3666  16502.3617  0.0049  0.0108  0.0185  0.0018
  303   34  16553.8217  16553.8163  0.0054  0.0154  0.0200  0.0026
  350   34  16840.4948  16840.5128 -0.0180  0.0191  0.0202  0.0033
  363   34  16903.6824  16903.6802  0.0022  0.0179  0.0206  0.0031

On Mon, Oct 8, 2018 at 3:48 PM David Nidever <dnidever@email.noao.edu> wrote:
The WAVE column comes directly from the end-of-night wavelength solution at the same dither position.

On Mon, Oct 8, 2018 at 3:47 PM David Nidever <dnidever@email.noao.edu> wrote:
Thanks.

I have a meeting in 15 min and am then heading home.  Here's the first linelist if you want to take a look.
I'm trying to figure out the best way to plot these:
/uufs/chpc.utah.edu/common/home/u0914350/wave/22310093_airglow.dat

-------------
oct_18a taken as subset of David's lines, but refined from literature from http://adsabs.harvard.edu/abs/2016JQSRT.168..142B

Dear David,

I am a post-doc working for Chad Bender and Chad mentioned that you are interested in finding a list of accurate OH sky emission line wavelengths to help wavelength calibrating APOGEE data.  I have been using wavelengths calculated in the theoretical paper Brooke et al. (2016) http://adsabs.harvard.edu/abs/2016JQSRT.168..142B for a project where I am trying to forward model OH sky lines in HPF, NEID, and (hopefully) APOGEE data and the wavelengths appear to be pretty good.

I have attached the list that I use, which is from the supplementary electronic files from the paper which I have also attached.  To get the wavelengths, I simply just grab the calculated transition wavenumber from column 10 of OH-XX-Line_list_201804001442519260.txt and convert it to microns  as follows:
wavelength [um] = 1e4 / wavenumber [cm^-1]

Let me know if this helps.

Regards,

-Kyle Kaplan


---------------
December 2018
--------------

wavecal routines rewritten in Python. Multiframe solutions derived for each year using 4th order polynomial
allowing for chip shifts between different groups.

skycal written to correct frames based on sky lines, using 4 parameter model: global line fiber-dependent shift
plus chip offsets. Use with two frames to refine positions of airglow lines, which needed to be adjusted for green
chip, from end of night frames 16930029 and 22430033 (because these nights were included in the annual group
solutions). This was done by writing temporary wavelength solution for the group save_Wavecal(group=) and then
runing skycal, but looking at output wavelengths.


