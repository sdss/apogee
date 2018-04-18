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


