;======================================================================
pro flatplot,flat,file

 ; make a jpg of the flat
 set_plot,'ps'
 device,file='a.eps',/encap,xsize=15,ysize=15
 low=0.5
 high=1.5
 disp=flat>low<high
 disp=255*(disp-low)/(high-low)
 tv,nint(disp)
 ;tvscl,flat>low<high
 device,/close
 spawn,'convert a.eps '+file+'.jpg'
 file_delete,'a.eps'
 set_plot,'X'

end
