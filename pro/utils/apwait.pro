
pro apwait,file,time,silent=silent
  if not keyword_set(silent) then print,'waiting for file lock: ', file
  wait,time
end
