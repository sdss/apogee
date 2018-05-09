{

 count = count + 1

 if (count == 1) printf("%s\n", "WEBMARCS")
 if (count == 1) printf("%s\n", $0)
 if (count == 2) Teff = $1
 if (count == 4) logg = log($1)/log(10.0)
 if (count == 5) micro = $1
 if (count == 7) zero = $1
 if (count == 13) { c = $6;n=$7;o=$8}
 if (count == 14) { na = $1;mg=$2;al=$3;si=$4;s=$6;k=$9;ca=$10}
 if (count == 15) { ti = $2;v=$3;cr=$4;mn=$5;fe=$6;co=$7;ni=$8}
 
   if (layer >= 1 && layer < number + 3 ) {
      if (layer >= 3 ) {
        lgTau5[layer-2] =  $3
        T[layer-2] = $5
        Pe[layer-2] = $6
        Pg[layer-2] = $7  
      }
      layer = layer + 1
    }
 if ($2 == "Number") {
   number = $1
   printf("             %s\n",number)
   layer = layer + 1
 }
    if ( second >= 1 ) {
       second = second + 1
       if (second >= 2 && second < number + 2 ) {
          RHOX[second-1] = $8
          }
     }
    if ( $8 == "RHOX" ) second = second + 1

  micro = (2.24 - 0.3 * logg) * 1E5   # substitute your favorite formula here
  if (length(vmicro)>0) micro=vmicro * 1E5


}
END{ printf("%s\n","5000.0")
     for (i=1;i<=number;i++)
       {
         print i,lgTau5[i],T[i],Pe[i],Pg[i],RHOX[i]
       }
    printf("%10.3f\n",micro )
    printf("NATOMS     17 %6.3f\n",zero)
    printf("%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f\n",6.0,c,7.0,n,8.0,o)
    printf("%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f\n",11.0,na,12.0,mg,13.0,al,14.0,si)
    printf("%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f\n",16.0,s,19.0,k,20.0,ca)
    printf("%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f\n",22.0,ti,23.0,v,24.0,cr,25.0,mn,26.0,fe)
    printf("%10.1f%10.3f%10.1f%10.3f\n",27.0,co,28.0,ni)
    printf("NMOL          20\n")
printf("       1.1     107.0     108.0     607.0     608.0     708.0       6.1\n")
printf("       7.1       8.1      12.1     112.0     101.0     106.0     101.0\n")
printf("      22.1     822.0      14.1     114.0      26.1     126.0\n")
    }   

