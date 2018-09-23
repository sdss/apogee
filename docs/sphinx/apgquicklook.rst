.. _apgquicklook :

apgquicklook
===============================

apgquickook comprises a set of IDL routines for "quicklook"
and "quickred". The former processes frames as each read comes
in, and the latter processes frames after all reads for an
exposure are complete.

The main purpose of the quickred software is to bundle the annotated frames into 3D cubes and do a quick reduction of the exposure that allows for quick diagnostics on the mountain.  It is started once and exposure is finished.  The main steps are:

Bundle
Collapse datacube (3D->2D)
Extract all 300 fiberes
Diagnostics including determination of S/N vs. mag and at fiducial mag
Database update with diagnostic information and spectra
Compress the cubes, 2D images and 1D spectra




