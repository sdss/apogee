09/11/2013
I created this file to list the programs and  supporting files in this
directory and place a very short description of them. 
The detailed description and usage might be found inside programs. 
Fill free to edit this file. EM. 

readme.txt  
list of programs and files with short description

-------------------------------------------------
apogee1m  
?  no docs

apogeeql  
?  no docs

apogeeThar.py
EM: Program to check the stability of APOGEE instrument. 
It uses ArcLamp THAR data  and fits a gaussian function
to one spectral line x-center=43 and output fitting parameters
relative to reference data apRaw-09000004.fits fitting. 
 
apogeeTharPlot.py
EM: This program is supplemented for apogeeThar.py. 
It uses apogeeThar.py output saved as  apogeeThar.outfile, 
read this file  as a table and plot 
x-centers of gaussian function from this table. 

apRaw-03720068.fits 
EM: DomeFlat raw data using as a 
reference for aptest (see aptest description)

aptest
EM: program to  check if broken fibers in dome flat apogee data. 

fix_missing_bzero.py 
To fix the files in MJD 56532-56540 where the bzero field wasn't written
during the annotation, due to the pyfits 2.4->3.1 change.

list_ap1
EM: python code to list the headers of apogee data. 
Observers use it during observations to monitor of progress 
and record to night log.  

plugmapm2a.py
This script will transform a plPlugMapM file to a plPlugMapA file for
APOGEE, adding the 2-Mass JHK magnitudes and 2-MAS target name to the
table.

talk2ql.py
server to communicate with the apql_wrapper on port 10033
