.. _apogeeql:

apogeeql
===============================

The apogeeql actor orchestrates APOGEE data taking at both APO and LCO.
It interacts with the APOGEE instrument software, quicklook, quickred and
the hub to handle the images that are taken.  Once a read is finished the
apogeeql actor "annotates" the FITS header with extra information, copies the
annotated frame to a certain location on disk, and lets the quicklook know
that a read is ready to be checked.  At the beginning of an exposure
the actor also reads in necessary information (plugmap, calibration
data, etc.), creates a temporary plugmap file (plPlugMapA).  At the end of an exposure it
lets quickred know that an exposure is done and it should process it.
When the apogeeql actor is started, it starts processes for quicklook
and quickred and initializes sockets to communicate with them.




