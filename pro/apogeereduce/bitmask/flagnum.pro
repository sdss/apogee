function flagnum,flag

if flag eq 'APOGEE_FAINT' then num=0
if flag eq 'APOGEE_MEDIUM' then num=1
if flag eq 'APOGEE_BRIGHT' then num=2
if flag eq 'APOGEE_IRAC_DERED' then num=3
if flag eq 'APOGEE_WISE_DERED' then num=4
if flag eq 'APOGEE_NICE_DERED' then num=5
if flag eq 'APOGEE_NO_DERED' then num=6
if flag eq 'APOGEE_WASH_GIANT' then num=7
if flag eq 'APOGEE_WASH_DWARF' then num=8
if flag eq 'APOGEE_SCI_CLUSTER' then num=9
if flag eq 'APOGEE_EXTENDED' then num=10
if flag eq 'APOGEE_SHORT' then num=11
if flag eq 'APOGEE_INTERMEDIATE' then num=12
if flag eq 'APOGEE_LONG' then num=13
if flag eq 'APOGEE_DO_NOT_OBSERVE' then num=14
if flag eq 'APOGEE_SERENDIPITOUS' then num=15
if flag eq 'APOGEE_FIRST_LIGHT' then num=16
if flag eq 'APOGEE_ANCILLARY' then num=17
if flag eq 'APOGEE_M31_CLUSTER' then num=18
if flag eq 'APOGEE_MDWARF' then num=19
if flag eq 'APOGEE_HIRES' then num=20
if flag eq 'APOGEE_OLD_STAR' then num=21
if flag eq 'APOGEE_DISK_RED_GIANT' then num=22
if flag eq 'APOGEE_KEPLER_EB' then num=23
if flag eq 'APOGEE_GC_PAL1' then num=24
if flag eq 'APOGEE_MASSIVE_STAR' then num=25
if flag eq 'APOGEE_SGR_DSPH' then num=26
if flag eq 'APOGEE_KEPLER_SEISMO' then num=27
if flag eq 'APOGEE_KEPLER_HOST' then num=28
if flag eq 'APOGEE_FAINT_EXTRA' then num=29
if flag eq 'APOGEE_CHECKED' then num=31
if flag eq 'LIGHT_TRAP' then num=0
if flag eq 'APOGEE_FLUX_STANDARD' then num=1
if flag eq 'APOGEE_STANDARD_STAR' then num=2
if flag eq 'APOGEE_RV_STANDARD' then num=3
if flag eq 'APOGEE_SKY' then num=4
if flag eq 'APOGEE_SKY_BAD' then num=5
if flag eq 'APOGEE_GUIDE_STAR' then num=6
if flag eq 'APOGEE_BUNDLE_HOLE' then num=7
if flag eq 'APOGEE_TELLURIC_BAD' then num=8
if flag eq 'APOGEE_TELLURIC' then num=9
if flag eq 'APOGEE_CALIB_CLUSTER' then num=10
if flag eq 'APOGEE_GC_GIANT' then num=11
if flag eq 'APOGEE_GC_SUPER_GIANT' then num=12
if flag eq 'APOGEE_EMBEDDEDCLUSTER_STAR' then num=13
if flag eq 'APOGEE_LONGBAR' then num=14
if flag eq 'APOGEE_EMISSION_STAR' then num=15
if flag eq 'APOGEE_KEPLER_COOLDWARF' then num=16
if flag eq 'APOGEE_MIRCLUSTER_STAR' then num=17

if flag eq 'APOGEE_1MTARGET' then num=22
if flag eq 'APOGEE_CHECKED' then num=31

return,num

end
