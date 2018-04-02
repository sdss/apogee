#!/bin/csh -f

#  B. Plez
# exctraction of line lists from my database, and reformatting for turbospectrum. 

# The selection of lines is made at log(gf*lambda)-chie*theta > strongest*10^CUT, 
# with theta=5040./tempselect,
# chie is the lower energy level of the line,
# and strongest is the largest log(gf*lambda)-chie*theta of the ful line list.
# Chose a very low cut (e.g. CUT=-99.) if you want the full line list.

set CUT = -4
set tempselect = 3500.
set weakratio = 1.e${CUT}

set LAMBDAMIN = '15000'
set LAMBDAMAX = '17000'

set DIRECTORY = ${LAMBDAMIN}-${LAMBDAMAX}_cut${CUT}_new

mkdir ${DIRECTORY}

#------------------------------------------------------------------
# Selection de raies LaO pour Bsyn. BPz 22/06-09
#------------------------------------------------------------------

#../PROGRAMS/translatelinelists_autocount_identif <<eof
#../LaO/A-X/linelist_LaO_AX-Re2_equal_1.dat
#LaO
#${DIRECTORY}/LaO_AX-Re2_equal_1-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
#${LAMBDAMIN}
#${LAMBDAMAX}
#${tempselect}
#${weakratio}
#eof

#------------------------------------------------------------------
# Selection de raies C2 pour Bsyn. BPz 11/03-99
#------------------------------------------------------------------

translatelinelists_autocount_identif <<eof
../C2/linelist/PLEZ_QUERCI/C2_Ballik-Ramsay_Querci.dat
C1212Q
${DIRECTORY}/C12C12-BR_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../C2/linelist/PLEZ_QUERCI/C2_Ballik-Ramsay_Querci.dat
C1213Q
${DIRECTORY}/C12C13-BR_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../C2/linelist/PLEZ_QUERCI/C2_Ballik-Ramsay_Querci.dat
C1313Q
${DIRECTORY}/C13C13-BR_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../C2/linelist/PLEZ_QUERCI/C2_Phillips_Querci.dat
C1212Q
${DIRECTORY}/C12C12-P_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../C2/linelist/PLEZ_QUERCI/C2_Phillips_Querci.dat
C1213Q
${DIRECTORY}/C12C13-P_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../C2/linelist/PLEZ_QUERCI/C2_Phillips_Querci.dat
C1313Q
${DIRECTORY}/C13C13-P_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../C2/linelist/PLEZ_QUERCI/C2_Swan_Querci.dat
C1212Q
${DIRECTORY}/C12C12-S_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../C2/linelist/PLEZ_QUERCI/C2_Swan_Querci.dat
C1213Q
${DIRECTORY}/C12C13-S_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../C2/linelist/PLEZ_QUERCI/C2_Swan_Querci.dat
C1313Q
${DIRECTORY}/C13C13-S_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

#------------------------------------------------------------------
# Selection de raies CH pour Bsyn. BPz 27/02-99
#------------------------------------------------------------------

#translatelinelists_autocount_identif <<eof
#../CH/linelist/scan_ch.dat
#C12H
#${DIRECTORY}/C12H-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
#${LAMBDAMIN}
#${LAMBDAMAX}
#${tempselect}
#${weakratio}
#eof
#
#translatelinelists_autocount_identif <<eof
#../CH/linelist/scan_ch.dat
#C13H
#${DIRECTORY}/C13H-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
#${LAMBDAMIN}
#${LAMBDAMAX}
#${tempselect}
#${weakratio}
#eof
#

#------------------------------------------------------------------
# Selection de raies CN pour Bsyn. BPz 18/02-99
#------------------------------------------------------------------

translatelinelists_autocount_identif <<eof
../CN/linelist/linelistCN_alliso_R240511.dat
C2N4sg
${DIRECTORY}/C12N14R_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../CN/linelist/linelistCN_alliso_R240511.dat
C2N5sg
${DIRECTORY}/C12N15R_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../CN/linelist/linelistCN_alliso_R240511.dat
C3N4sg
${DIRECTORY}/C13N14R_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../CN/linelist/linelistCN_alliso_R240511.dat
C3N5sg
${DIRECTORY}/C13N15R_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../CN/linelist/linelistCN1214V130710.dat
C12N14
${DIRECTORY}/C12N14V130710_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../CN/linelist/linelistCN1215V130710.dat
C12N15
${DIRECTORY}/C12N15V130710_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../CN/linelist/linelistCN1314V130710.dat
C13N14
${DIRECTORY}/C13N14V130710_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../CN/linelist/linelistCN1315V130710.dat
C13N15
${DIRECTORY}/C13N15V130710_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof


#------------------------------------------------------------------
# Selection de raies CaH pour Bsyn. BPz 27/02-99
#------------------------------------------------------------------

translatelinelists_autocount_identif <<eof
../CaH/linelist/linelist_CaH_A-X.dat
CaH
${DIRECTORY}/CaH_A-X-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../CaH/linelist/linelist_CaH_B-X.dat
CaH
${DIRECTORY}/CaH_B-X-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof


#------------------------------------------------------------------
# Selection de raies MgH Kurucz pour Bsyn. BPz 27/02-99
#------------------------------------------------------------------
#
#foreach ku (24 25 26) 
#
#translatelinelists_autocount_identif <<eof
#../MgH/linelist/MgH_kurucz.dat
#${ku}MgH
#${DIRECTORY}/${ku}MgH-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
#${LAMBDAMIN}
#${LAMBDAMAX}
${tempselect}
${weakratio}
#eof
#
#end


#------------------------------------------------------------------
# Selection de raies MgH Skory et al. pour Bsyn. BPz 30/09-09
#------------------------------------------------------------------

translatelinelists_autocount_identif <<eof
../MgH/NEWLIST/MgH_Skory_Weck_Stancil.list
MgH
${DIRECTORY}/MgH-Skory-Weck_Stancil-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

#------------------------------------------------------------------
# Selection de raies NH pour Bsyn. BPz 27/02-99
#------------------------------------------------------------------

translatelinelists_autocount_identif <<eof
../NH/linelist/NH_kurucz.dat
14NH
${DIRECTORY}/14NH-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../NH/linelist/NH_kurucz.dat
15NH
${DIRECTORY}/15NH-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof


#------------------------------------------------------------------
# Selection de raies OH A-X Kurucz pour Bsyn. BPz 27/02-99
#------------------------------------------------------------------

#translatelinelists_autocount_identif <<eof
#../OH_A-X/linelist/OH_A-X_kurucz.dat
#16OH
#${DIRECTORY}/16OH_A-X-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
#${LAMBDAMIN}
#${LAMBDAMAX}
#${tempselect}
#${weakratio}
#eof

#translatelinelists_autocount_identif <<eof
#../OH_A-X/linelist/OH_A-X_kurucz.dat
#18OH
#${DIRECTORY}/18OH_A-X-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
#${LAMBDAMIN}
#${LAMBDAMAX}
#${tempselect}
#${weakratio}
#eof


#------------------------------------------------------------------
# Selection de raies SiH pour Bsyn. BPz 27/02-99
#------------------------------------------------------------------

translatelinelists_autocount_identif <<eof
../SiH/linelist/SiH_kurucz.dat
28SiH
${DIRECTORY}/28SiH-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof


translatelinelists_autocount_identif <<eof
../SiH/linelist/SiH_kurucz.dat
29SiH
${DIRECTORY}/29SiH-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof


translatelinelists_autocount_identif <<eof
../SiH/linelist/SiH_kurucz.dat
30SiH
${DIRECTORY}/30SiH-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof


#------------------------------------------------------------------
# Selection de raies FeH pour Bsyn. BPz 23/08-99
#------------------------------------------------------------------

translatelinelists_autocount_identif <<eof
../FeH/linelist/FeH_scaledLanghoff.list
FeH
${DIRECTORY}/FeH-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof



#------------------------------------------------------------------
# Selection de raies TiO VALD pour Bsyn. BPz 14/12-11
#------------------------------------------------------------------

foreach iso (48 )

../PROGRAMS/translatelinelists_autocount_identif <<eof
../TiOVALD/linelist/linelist_reduced${iso}_all_deltacorr_lab.dat
${iso}TiOV
${DIRECTORY}/${iso}TiOVALD-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

end

foreach iso (46 47 49 50)

../PROGRAMS/translatelinelists_autocount_identif <<eof
../TiOVALD/linelist/linelist_reduced${iso}_all_deltacorr.dat
${iso}TiOV
${DIRECTORY}/${iso}TiOVALD-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

end

#------------------------------------------------------------------
# Selection de raies TiO pour Bsyn. BPz 20/01-99
#------------------------------------------------------------------

#foreach iso (46 47 48 49 50)
#
#../PROGRAMS/translatelinelists_autocount_identif <<eof
#../TiOnew/linelist/TiO${iso}_lab.dat
#${iso}TiO
#${DIRECTORY}/${iso}TiO-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
#${LAMBDAMIN}
#${LAMBDAMAX}
#${tempselect}
#${weakratio}
#eof
#
#end

#------------------------------------------------------------------
# Selection de raies VO pour Bsyn. BPz 27/02-99
#------------------------------------------------------------------

translatelinelists_autocount_identif <<eof
../VO/linelist/linelistVO_ALL.dat
VO
${DIRECTORY}/VO-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof


#------------------------------------------------------------------
# Selection de raies ZrO pour Bsyn. BPz 20/01-99
#------------------------------------------------------------------

foreach iso (90 91 92 94 96)

gunzip ../ZrO/linelist/linelist${iso}ZrO.dat.gz
translatelinelists_autocount_identif <<eof
../ZrO/linelist/linelist${iso}ZrO.dat
${iso}ZrO
${DIRECTORY}/${iso}ZrO-bsyn_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof
#gzip ../ZrO/linelist/linelist${iso}ZrO.dat

end

#------------------------------------------------------------------
# Selection de raies CO pour Bsyn. BPz 08/09-99
#------------------------------------------------------------------

foreach co (26 27 28 36 37 38 46)

translatelinelists_autocount_identif <<eof
../CO/linelist/CO_all_Goorvitch.linelist
${co}CO
${DIRECTORY}/${co}CO_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

end


#------------------------------------------------------------------
# Selection de raies OH AX Goldman pour Bsyn. BPz 
#------------------------------------------------------------------


translatelinelists_autocount_identif <<eof
../OH_A-X/linelist/GOLDMAN/4600.fix.sorted
OH-AX
${DIRECTORY}/OH_HITRAN-AX.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

#------------------------------------------------------------------
# Selection de raies OH IR Goldman pour Bsyn. BPz 
#------------------------------------------------------------------

translatelinelists_autocount_identif <<eof
../OH/linelist/HITRAN/hitran.combined.acm.eav.ME4957
OH-IR
${DIRECTORY}/OH_HITRAN-IR.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

#------------------------------------------------------------------
# Selection de raies SiS pour Bsyn. BPz 24/02-2009
#------------------------------------------------------------------

translatelinelists_autocount_identif <<eof
../SiS/SiS_Camietal_2009.dat
SiS
${DIRECTORY}/SiS_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

#------------------------------------------------------------------
# Selection de raies SiO pour Bsyn. BPz 05/05-2009
#------------------------------------------------------------------

translatelinelists_autocount_identif <<eof
../SiO/linelist/SiO_all_Langhoff.linelist
28SiO
${DIRECTORY}/28SiO_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../SiO/linelist/SiO_all_Langhoff.linelist
29SiO
${DIRECTORY}/29SiO_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
../SiO/linelist/SiO_all_Langhoff.linelist
30SiO
${DIRECTORY}/30SiO_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

#------------------------------------------------------------------

