.. _speclib:

speclib
===============================

atmos.py

Routines for working with model atmospheres, including
finding and filling "holes" in a rectangular atmospheres grid.
Main routine is fill_holes(). Takes as input command line arguments with
the grid parameters (start, end, step), and
looks for existing atmospheres files. If file not found, adopt "nearest" atmosphere
with algorithm as given in find_filler(), and copy this into the missing
file location. For Kurucz models, need to then "fix" this file for appropriate
parameters; MARCS models are left as is.

Also outputs holefile that gives location of holes and how they were filled.

dw.py

Development routine to create output plots and web pages comparing syntheses
computed with different wavelength spacings

isochrones.py

Routine to read a Padova isochrone file and return data in structured array

lsf.py

Routines for constructing kernels and convolving spectrumn with a broadening
profile using fast matrix operations, written largely by Jo Bovy. Includes
LSF kernels, macroturbulence (Gaussian) kernels, and rotation kernels. Main
routine is convolve, which does the convolution.

pca.py

Routines for creating and testing PCA-compresed synthetic spectra grids. Main
routine is pca, which does the PCA compression. Also includes test, which
will compare the PCA spectra to the raw ones for a representative sample of
spectra (from existing test.ipf file), and run some FERRE fits to see how well 
parameters are recovered.

rbf.py

sample.py

Routines for creating sample of test parameters/abundances to be used for
testing of FERRE and/or for construction of training sample for an interpolator,
etc.

sim.py

synth.py : 

Routines related making synthetic spectra using Turbospectrum. The low
level routine mkturbospec does a single Turbospectrum synthesis and returns the
resulting spectra (both raw and normalized). This is used by several higher
level routines: 

 - mkgrid will make a regular grid of synthetic spectra, as specified
by an input planfile. These are used as the basis for the libraries used
by FERRE, after (optionally) hole filling, convolution with an LSF, PCA compression, and 
bundling into a FERRE format library.

 - mksynth will make a series of synthetic spectra as specified by
an input file with parameters and abundances on each line. These might
be used as an input into a neural network interpolator..

synth also includes a number of utility routines to accomplish the above
tasks.

welem.py

window.py



