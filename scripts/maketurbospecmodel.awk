# awk command file for making input atmosphere model for Turbospec from Kurucz atmosphere input
# will print out all atmosphere levels by default:
#    if tmin=X is specified on command line, will only print out levels with t>X

{
# loop through file counting lines with t>tmin and storing output lines in array
count = count + 1

if (count == 1) teff = $2
if (count == 1) logg = $4

if (start > 0 && start <= lines )  {
    drho = $1 - rho
    #tau = tau + ($5+abross)/2*drho
    tau = tau + $5*drho
    rho = $1
    abross = $5
    #print $0, tau, log(tau)/log(10)
    if ( $2 > tmin && tau > taumin && start > trim) {
      nlines = nlines + 1
      out[j++]=$0
    }
    micro = $7 * 1.0
    start = start + 1
    }

if ($1 == "READ") lines = $3 
if ($1 == "READ") start = 1
#if ($1 == "READ") printf("KURUCZ   %s 5000. %s 0 0.\n",nlines,logg)

}
# now go back and output all of the desired layers, with the correct header
END{
  printf("KURUCZ   %s 5000. %s 0 0.\n",nlines,logg)
  for (j=0; j<nlines; j++) print out[j]
}
