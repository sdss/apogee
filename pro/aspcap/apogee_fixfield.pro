function apogee_fixfield,input

   if size(input,/type) eq 3 then cname='Unknown' else cname=input

   ;; remove APG_ if necessary
   if(strmatch(cname, 'APG_*')) then $
      cname= strmid(cname, 4)

   if(strmatch(cname, 'APGS_*')) then $
      cname= strmid(cname, 5)

   case cname of
      'BULGE_04+00': cname='004+00'
      'HALO_M13': cname='M13'
      'DISK_27+00': cname='027+00'
      'DISK_30+08': cname='030+08'
      else:
   endcase

   return, cname

end
