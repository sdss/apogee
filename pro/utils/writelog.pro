pro writelog,logfile,line
openw,log,/get_lun,logfile,/append
printf,log,line
free_lun,log
end

