{

count = count + 1
if (count == 1) print "simple moog format linelist from ",FILENAME
if (substr($0,1,1) != "#") {

 lam = $1 * 10.0
 EP1 = substr($0,55,12) * 1.0
 EP2 = substr($0,83,12) * 1.0
 orggf = substr($0,11,7)
 newgf = substr($0,19,7)
 astgf = substr($0,35,7)
 isof = substr($0,146,6) * 1.0
 hyp = substr($0,137,6) * 1.0
 specid1 = substr($0,47,5) * 1.0
 specid2 = substr($0,52,3) 
 vdW1 = substr($0,123,6)
 vdW2 = substr($0,181,6)
 if (specid1 < 100) specid = specid1 + specid2 * 10.0
 else {
   if (specid2 == ".12" && specid1 == 606) { dis = 6.210; specid2 = .01212 } 
   else if (specid2 == ".00" && specid1 == 606) { dis = 6.210; specid2 = .01212 } 
   else if (specid2 == ".13" && specid1 == 606) { dis = 6.210; specid2 = .01213 } 
   else if (specid2 == ".33" && specid1 == 606) { dis = 6.210; specid2 = .01313 } 
   else if (specid2 == ".12" && specid1 == 607) { dis = 7.760; specid2 = .01214 } 
   else if (specid2 == ".00" && specid1 == 607) { dis = 7.760; specid2 = .01214 } 
   else if (specid2 == ".13" && specid1 == 607) { dis = 7.760; specid2 = .01314 } 
   else if (specid2 == ".15" && specid1 == 607) { dis = 7.760; specid2 = .01215 } 
   else if (specid2 == ".12" && specid1 == 608) { dis =11.092; specid2 = .01216 } 
   else if (specid2 == ".00" && specid1 == 608) { dis =11.092; specid2 = .01216 } 
   else if (specid2 == ".13" && specid1 == 608) { dis =11.092; specid2 = .01316 } 
   else if (specid2 == ".17" && specid1 == 608) { dis =11.092; specid2 = .01217 } 
   else if (specid2 == ".18" && specid1 == 608) { dis =11.092; specid2 = .01218 } 
   else if (specid2 == ".16" && specid1 == 108) { dis = 4.392; specid2 = .00116 } 
   else if (specid2 == ".00" && specid1 == 108) { dis = 4.392; specid2 = .00116 } 
   else if (specid2 == ".28" && specid1 == 114) { dis = 3.060; specid2 = .00128 } 
   else if (specid2 == ".00" && specid1 == 114) { dis = 3.060; specid2 = .00128 } 
   else if (specid2 == ".29" && specid1 == 114) { dis = 3.060; specid2 = .00129 } 
   else if (specid2 == ".30" && specid1 == 114) { dis = 3.060; specid2 = .00130 } 
   else if (specid2 == ".11" && specid1 == 101) { dis = 4.478; specid2 = .00101 } 
   else if (specid2 == ".12" && specid1 == 101) { dis = 4.478; specid2 = .00102 } 
   else if (specid2 == ".01" && specid1 == 101) { dis = 4.478; specid2 = .00101 } 
   else if (specid2 == ".02" && specid1 == 101) { dis = 4.478; specid2 = .00102 } 
   else if (specid2 == ".56" && specid1 == 126) { dis = 1.590; specid2 = .00156 }
   specid = specid2 + specid1
   }

 if (EP1 * 1 < 0) ep = EP1 * -1
 else ep = EP1
 if (EP2 * 1 > 0 && EP2 * 1 < ep ) ep = EP2
 if (EP2 * 1 < 0 && EP2 * -1 < ep ) ep = -EP2
 ep = ep * 1.2398e-4

 if (astgf != "       ") gf = astgf + hyp + isof 
 else if (newgf != "       ") gf = newgf + hyp + isof 
 else gf = orggf + hyp + isof

 if (vdW2 == "      ") vdW = vdW1
 else if ( vdW2 * 1.0 < -1.0) vdW = vdW2
 else vdW = vdW1

 if (specid < 100 && specid > 1.9 && ep < 40.0) {
    if (specid2 * 1.0 < .03) printf("%10.3f%10.2f%10.3f%10.3f%10.3f\n",lam,specid,ep,gf,vdW)
    }
 else if (specid > 100) printf("%10.3f%10.5f%10.3f%10.3f%10.3f%10.3f\n",lam,specid,ep,gf,vdW,dis)
 } 
} 
