.. _apogeereduce:

apogeereduce
===============================

apogeereduce provides for the main APOGEE data reduction pipeline.
It is written in IDL.

The pipeline processes images from raw data cubes to flux and wavelength calibrated, sky and telluric corrected 1D spectra with measured RVs.  The main steps are:

APRED: Reduces images to 1D spectra

AP3D:  Collapse the data cube from 3D→2D.  Reference pixel, linearity, dark current correction.  Cosmic ray and saturation flagging/repair.  Collapse to 2D using Fowler/Up-the-ramp.  Flat field correction and construct error array.

AP2D: Extract the 300 spectra.  Flux/throughput correction and initial wavelength calibration.

AP1DVisit: Dither measurement, sky/telluric correction, dither combination, flux calibration and initial RV determination.

APSTAR: Combined multiple spectra of the same star, measures accurate RVs and puts final spectrum on a common rest wavelength scale.

There are also programs to create the various calibration products (PSF, LSF, Wave, Dark, etc.). There are also a number of scripts to automate the reduction process.



