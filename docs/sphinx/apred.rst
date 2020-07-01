.. _apred:

apred
===============================

apred provides Python routines related to APOGEE data reduction. To date,
these are limited because the reduction is done in an IDL environment, but
there are a few utility routines used for quality assurance. It is expected
that additional functionality will be implemented in Python in the future.

In late 2018, routines for doing wavelength calibration were added, 
see  :ref:`Wavelength Calibration`

lsfmap.py

routines to create maps of FWHM and resolution across the APOGEE detectors,
given pre-existing apLSF files

.. automodule:: apogee.apred.lsfmap
         :members:

through.py

routines to make throughput-related plots, given pre-existing reduction
summary files

.. automodule:: apogee.apred.through
         :members:

mjdcube.py

routine to make a cube of simply-reduced (CDS?) frames for an entire
night, used by the apogeereduce persistence correction routines.

.. automodule:: apogee.apred.mjdcube
         :members:

rv.py

routines related to radial velocity (RV) measurement.

These include new routines to replace the old IDL routines for
RV determination and visit combination. The new routines are
based on the use of David Nidever's doppler package, which is
included in the apogee software data product as a sub-package.

RVs are processed at the field level. For each field, information
about all of the visits are saved by the ap1dvisit processing
in apVisitSum files. The rv.doppler_rv() routine concatenates
these files into an allvisit table. Unique objects are identified
in this table, and a list of visit file names for each objecct
is accumulated. These lists are then passed to the wrapper
routine rv.dorv(), which calls doppler.rv.jointfit() to deterine
RVs for each visit.

doppler.rv.jointfit() determines RV in several iterations. First
each input spectra is cross-correlated against a pre-determined
set of model templates. Stellar parameters from the template
that produces the highest cross-correlation peak are stored,
along with the RV corresponding to this peak. After all visits
are processed, S/N-weighted average parameters and RV are determined,
and these are used to provide first guesses for a least-squares
fit to the stack of images, in which stellar parameters and
individual visit RVs are the free parameters. The process
goes through a few iterations, masking out poorly fitting
pixels before determining a final set of parameters and RVs.
The best fitting model is then used to provide a final set
of cross-correlation functions, although these are not used
for the RVs, as the least-squares fit provides the final RVs.

The output CCFs, however, are used to attempt to identify
multi-component systems. This is done using autonomous gaussian
decomposition as implemented in the (external) gausspy package,
also included as a sub-package. This is implemented in the
driver routine gauss_decomp(), which also includes a simple
algorithm for filtering out likely spurious peaks.

The rv.dorv() wrapper does both the doppler.rv.jointfit()
and the rv.gauss_decomp(), as well as a rv.dop_plot() routine
to make a series of plots. The wrapper routine is used so
that the driver routine rv.doppler_rv() can run multiple
stars in different threads. Note that, for this to be effective,
the doppler.rv and lower level routines must be single-threaded,
accomplished by setting environment variables:

setenv OMP_NUM_THREADS 1
setenv OPENBLAS_NUM_THREADS 1
setenv MKL_NUM_THREADS 1
setenv VECLIB_MAXIMUM_THREADS 1
setenv NUMEXPR_NUM_THREADS 1

before loading the packages. This is significantly more efficient
that running a single object at a time with multi-threaded
lower levels.

After all of the objects have been run, the RV results are loaded
into the allvisit structure. In addition, a new allstar structure
is created with an entry for each star that includes the average
heliocentric (actually, barycentric) RVs, and the VSCATTER measurement
of scatter.
