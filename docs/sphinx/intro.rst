
.. _intro:

Introduction to apogee
===============================

The apogee package comprises a wide range of routines related to the reduction and analysis
of APOGEE data. It comprises several different components that were previously located in
separate SDSS SVN pacakages:

- routines related to the realtime collection of data (SVN apogeeql and apgquicklook)
- routines related to the reduction of APOGEE data (SVN apogeereduce)
- routines related to the creation of synthetic spectra and the bundling of these 
spectra into libraries for use with FERRE (SVN speclib)
- routines that run the spectral parameter and abundance analysis (ASPCAP) pipeline (SVN idlwrap)
- the FERRE program that interpolates within the spectral grids and finds the best solution
given input spectra and uncertainties
