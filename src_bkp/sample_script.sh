#sample_script.sh
#!/bin/bash
outputfile=outputfilename
wdir=/home/mjh92/picomax/picomax/demo
picomax_exe=/home/mjh92/project/picomax/src/picomax
nband=30
nchi=9
encut=400
nfreq=1001
dfreq=0.05
kk=1
delta=0
epsilon=0.3
qnum=20,20,10,10,20
qpath=L,G,X,W,K,G
crystal=zb
a=4.35
vg=3,4,8,11,-0.21,0,0.04,0.08,-0.21,0,0.04,0.08
refpoint=0,0,0



${picomax_exe} -switch epsij\
    -outputfile ${outputfile} -wdir ${wdir}\
    -nband ${nband} -neps ${nchi} -encut ${encut}\
    -crystal {crystal} -a ${a} -vg ${vg}\
    -kpoint 1 -kpointorder 7\
    -nfreq ${nfreq} -dfreq ${dfreq} -kk ${kk} -delta ${delta} -epsilon ${epsilon}\
    -qnum ${qnum} -qpath ${qpath}\
    >> ${wdir}/${outputfile}.log