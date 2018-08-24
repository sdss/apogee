#!/bin/csh -f

set cmd = $1

if ( ! $?QUERYHOST ) setenv QUERYHOST `hostname`

if ( $cmd == apred ) then
  setenv QUERYPORT 1050
else if ( $cmd == aspcap ) then
  setenv QUERYPORT 1051
else
  setenv QUERYPORT 1052
endif

setenv IDL_CPU_TPOOL_NTHREADS 1
if ( $cmd == apred ) then
  setenv APOGEE_MAXRUN 8
  if ( $?TACC_SYSTEM) setenv APOGEE_MAXRUN 10
  setenv IDL_CPU_TPOOL_NTHREADS 4
else if ( $cmd == mkrbf || $cmd == pca) then
  setenv APOGEE_MAXRUN 1
else if ( $cmd == mkgrid || $cmd == bundle ) then
  setenv APOGEE_MAXRUN 32
  if ( $?TACC_SYSTEM) setenv APOGEE_MAXRUN 40
else if ( $cmd == mkgridlsf ) then
  setenv APOGEE_MAXRUN 12
else 
  echo Unrecognized command in mkslurm
endif

set top = `pwd`

set prog = `echo $cmd | awk '{print $1}'`
shift

# Create the batch file
mkdir slurm
foreach node ( $cmd )
  set file = slurm/"$node"
  rm -f $file
  touch $file
  echo '#!'/bin/csh >> $file

  if ( $?TACC_SYSTEM ) then
    echo "#SBATCH --partition=gpu" >>$file
    echo "#SBATCH --time=240:00:00" >>$file
    echo "#SBATCH --ntasks-per-node=20" >>$file
    echo "#SBATCH --nodes=1" >>$file
  endif
  if ( $?UUFSCELL ) then
    echo "#SBATCH --account=sdss-kp" >> $file
    echo "#SBATCH --partition=sdss-kp" >> $file
    echo "#SBATCH --time=240:00:00" >> $file
    echo "#SBATCH --ntasks=16" >> $file
    echo "#SBATCH --nodes=1" >> $file
  endif

  echo "#SBATCH -o $prog.out" >>$file
  echo "#SBATCH -e $prog.out" >>$file
  echo "setenv QUERYHOST " $QUERYHOST >>$file
  echo "setenv QUERYPORT " $QUERYPORT >>$file
  echo "setenv APOGEE_MAXRUN " $APOGEE_MAXRUN >>$file
  echo "setenv IDL_CPU_TPOOL_NTHREADS " $IDL_CPU_TPOOL_NTHREADS >>$file

  echo cd $top >> $file
  echo runplans $cmd "$*" >> $file
  echo wait >> $file
  echo echo "DONE" >> $file
end

endif
exit
