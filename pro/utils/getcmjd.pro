; getcmjd converts numerical APOGEE file number into character MJD
function getcmjd,num,mjd=mjd
 imjd=num/10000
 mjd=55562+imjd
 cmjd=string(format='(i5.5)',mjd)
 return,cmjd
end
