#!/bin/csh -f
##################################################################################################
# Output turbospectrum/babsma format compatible 
# A control plot (interpol_check.ps) is displayed at the end.
# Extrapolation is not advised, even if allowed by this program.
# Requires an "cubic" set of 8 MARCS binary format models,
# in other words 
# !!!!!   MODELS MUST DIFFER 2 BY 2 BY ONLY ONE PARAMETER !!!!!!  
# !!!!!!! ORDER OF THE INPUT MODELS MATTERS !!!!!!!
# here is the order of the files
# model1: Tefflow logglow zlow
# model2: Tefflow logglow zup
# model3: Tefflow loggup zlow
# model4: Tefflow loggup zup
# model5: Teffup logglow zlow
# model6: Teffup logglow zup
# model7: Teffup loggup zlow
# model8: Teffup loggup zup
######################################################################################################

set model_path = Testwebformat
#set model_path = /u16/marcs

#MARCS binary format (.true.) or MARCS ASCII web format (.false.)?
set marcs_binary = '.false.'
#set marcs_binary = '.true.'

set model1 = p5500:g+4.0:m0.0:t01:ap:z-0.50:a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod
set model2 = p5500:g+4.0:m0.0:t01:st:z+0.00:a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod
set model3 = p5500:g+5.0:m0.0:t01:ap:z-0.50:a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod
set model4 = p5500:g+5.0:m0.0:t01:st:z+0.00:a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod
set model5 = p6000:g+4.0:m0.0:t01:ap:z-0.50:a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod
set model6 = p6000:g+4.0:m0.0:t01:st:z+0.00:a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod
set model7 = p6000:g+5.0:m0.0:t01:ap:z-0.50:a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod
set model8 = p6000:g+5.0:m0.0:t01:st:z+0.00:a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod

#enter here the values requested for the interpolated model 
foreach Tref   ( 5750 )
foreach loggref ( +4.5 )
foreach zref ( -0.25 )
set modele_out = Testout/${Tref}g${loggref}z${zref}.interpol
set modele_out2 = Testout/${Tref}g${loggref}z${zref}.alt
#set modele_out = scratch

#### the test option is set to .true. if you want to plot comparison model (model_test)
set test = '.true.'
set model_test = 'Testwebformat/p5750:g+4.5:m0.0:t01:ap:z-0.25:a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod'

# interpolation program (for further details see interpol_modeles.f)
interpol_modeles <<EOF
'${model_path}/${model1}'
'${model_path}/${model2}'
'${model_path}/${model3}'
'${model_path}/${model4}'
'${model_path}/${model5}'
'${model_path}/${model6}'
'${model_path}/${model7}'
'${model_path}/${model8}'
'${modele_out}'
'${modele_out2}'
${Tref}
${loggref}
${zref}
${test}
${marcs_binary}
'${model_test}'
EOF

echo 'control plot loading...'
#---------------------------------------------------------------------------------------------------
# control plot with sm
if ( ${test} == ".true.") then
set testsm = '1'
else
set testsm = '0'
endif

sm << eof
verbose 0
data "${modele_out}"
lines 2 0
read {tau 1 T 2 Pe 3 Pg 4 taur 6}
define ndp  (DIMEN(tau))
data "modele.sm" 

define beg (1)
define end(\$ndp*1)
lines \$beg \$end
read  {tau1 1 T1 2 Pe1 3 Pg1 4 taur1 5}
define beg (\$ndp+1)
define end(\$ndp*2)
lines \$beg \$end
read  {tau2 1 T2 2 Pe2 3 Pg2 4 taur2 5}
define beg (\$end+1)
define end(\$ndp*3)
lines \$beg \$end
read  {tau3 1 T3 2 Pe3 3 Pg3 4 taur3 5}
define beg (\$end+1)
define end(\$ndp*4)
lines \$beg \$end
read  {tau4 1 T4 2 Pe4 3 Pg4 4 taur4 5}
define beg (\$end+1)
define end(\$ndp*5)
lines \$beg \$end
read  {tau5 1 T5 2 Pe5 3 Pg5 4 taur5 5}
define beg (\$end+1)
define end(\$ndp*6)
lines \$beg \$end
read  {tau6 1 T6 2 Pe6 3 Pg6 4 taur6 5}
define beg (\$end+1)
define end(\$ndp*7)
lines \$beg \$end
read  {tau7 1 T7 2 Pe7 3 Pg7 4 taur7 5}
define beg (\$end+1)
define end(\$ndp*8)
lines \$beg \$end
read  {tau8 1 T8 2 Pe8 3 Pg8 4 taur8 5}
if (${testsm}) {
define beg (\$ndp*8+1)
lines \$beg 0
read  {tau9 1 T9 2 Pe9 3 Pg9 4 taur9 5}
}
device postportfile interpol_check.ps
ctype green
toplabel T_{eff}=${Tref} logg=${loggref} z=${zref}
ctype black
location 3000 15000 17000 32000
ylabel log(Pe)
limits tau Pe
box 0 2 0 0

lweight 4
ltype 0
ctype blue
relocate (4000 31000)
draw (4600 31000) 
putlabel 6 low T_{eff}
ctype red
relocate (4000 30500)
draw (4600 30500)
putlabel 6 up T_{eff}
ctype black
lweight 1
relocate (4000 30000)
draw (4600 30000)
putlabel 6 low logg
lweight 4
relocate (4000 29500)
draw (4600 29500)
putlabel 6 up logg
ltype 0
relocate (4000 29000)
draw (4600 29000)
putlabel 6 low [Fe/H]
ltype 1
relocate (4000 28500)
draw (4600 28500)
putlabel 6 up [Fe/H]

ctype blue
lweight 1
ltype 0
connect tau1 Pe1
ltype 1
connect tau2 Pe2
lweight 2
ltype 0
connect tau3 Pe3
ltype 1
connect tau4 Pe4
ctype red
lweight 1
ltype 0
connect tau5 Pe5
ltype 1
connect tau6 Pe6
lweight 2
ltype 0
connect tau7 Pe7
ltype 1
connect tau8 Pe8
ltype 0
lweight 1
ctype green
connect tau Pe
ctype black
if (${testsm}) {
connect tau9 Pe9   
}

location 3000 15000 3000 17000
xlabel log (\tau_{5000})
ylabel T
limits tau T
box
ctype blue
connect tau1 T1
ltype 1
connect tau2 T2
lweight 2
ltype 0
connect tau3 T3
ltype 1
connect tau4 T4
ctype red
lweight 1
ltype 0
connect tau5 T5
ltype 1
connect tau6 T6
ltype 0
lweight 2
connect tau7 T7
ltype 1
connect tau8 T8
ltype 0
ctype green
connect tau T
ctype black
if (${testsm}) {
connect tau9 T9  
}
lweight 1

location 19000 32000 17000 32000
xlabel log (\tau_{5000})
ylabel log(Pg)
limits tau Pg
box  
ctype blue
connect tau1 Pg1
ltype 1
connect tau2 Pg2
ltype 0
lweight 2
connect tau3 Pg3
ltype 1
connect tau4 Pg4
ctype red
ltype 0
lweight 1
connect tau5 Pg5
ltype 1
connect tau6 Pg6
ltype 0
lweight 2
connect tau7 Pg7
ltype 1
connect tau8 Pg8
ltype 0
ctype green
connect tau Pg
ctype black
if (${testsm}) {
connect tau9 Pg9 
}
lweight 1

if (${testsm}) {
location 19300 32000 3000 11000
ylabel actual error (%)
location 20000 32000 3000 6000
spline tau  T tau9 Tbis
limits tau9 (abs(Tbis-T9)/T9*100)
connect  tau9 (abs(Tbis-T9)/T9*100)
box
ylabel T
xlabel log (\tau_{5000})

location 20000 32000 6000 9000
spline tau Pe tau9 Pebis
limits tau9 (abs(10**Pebis-10**Pe9)/10**Pe9*100)
box 0 2 0 0
connect  tau9  (abs(10**Pebis-10**Pe9)/10**Pe9*100)
ylabel Pe

location 20000 32000 9000 12000
spline tau Pg tau9 Pgbis
limits tau9 (abs(10**Pgbis-10**Pg9)/10**Pg9*100)
connect  tau9 (abs(10**Pgbis-10**Pg9)/10**Pg9*100)
box 0 2 0 0
ylabel Pg
}

eof


gv interpol_check.ps&

rm -f modele.sm
end
end 
end
