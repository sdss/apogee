function getlocaldir
  local=getenv('APOGEE_LOCALDIR')
  if local eq '' then return,0 else return,addslash(local)
end

