.. _aspcap:

aspcap
===============================

ASPCAP is the APOGEE Stellar Parameters and
Chemical Abundances pipeline. The main engine of
ASPCAP is the code FERRE. Here, we refer to
ASPCAP as the wrapper which runs FERRE on APOGEE
data, orchestrating the different spectral libraries,
and producing output FITS files.

Through DR14, ASPCAP was mostly written in IDL and was
archived in the idlwrap SVN product. Some Python routines
were included for the post-calibration of the FERRE
output.

Various ASPCAP
routines are being ported to Python with a goal
of having ASPCAP written entirely in Python evenutally.

The migration has started with basic utility routines:

**aspcap.py**

.. automodule:: apogee.aspcap.aspcap
         :members:

**err.py**

.. automodule:: apogee.aspcap.err
         :members:

**loggcomp.py**

.. automodule:: apogee.aspcap.loggcomp
         :members:

**teffcomp.py**

.. automodule:: apogee.aspcap.teffcomp
         :members:

**cal.py**

.. automodule:: apogee.aspcap.cal
         :members:

**elem.py**

.. automodule:: apogee.aspcap.elem
         :members:

**ferre.py**

.. automodule:: apogee.aspcap.ferre
         :members:

**norm.py**

.. automodule:: apogee.aspcap.norm
         :members:


**persist.py**

.. automodule:: apogee.aspcap.persist
         :members:




