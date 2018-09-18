.. _apred:

apred
===============================

apred provides Python routines related to APOGEE data reduction. To date,
these are limited because the reduction is done in an IDL environment, but
there are a few utility routines used for quality assurance. It is expected
that additional functionality will be implemented in Python in the future.

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

