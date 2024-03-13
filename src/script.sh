#!/bin/bash

PICOMAX_EXE="/home/mjh92/projects/picomax/src/picomax"


# -numpw : number of plane waves {124, 342, 728, 1330, 2196}
# -numband : number of bands to consider {8}
# -gridscheme : k-grid scheme (0: Cohen-Cahdi, 1: Monkhorst-Pack)
# -gridorder : Cohen-Chadi order
# -b1/b2/b3 : reciprocal lattice vectors


# -freq : frequency [eV]
# -kk : kramers-kronig


# -q : q-vector
# -nq : total number of q-vectors
# -numq : 10,10,5,5,10
# -pathq L,G,X,W,K,G
# ==> 0->L 10->G 20->X 25->W 30->K 40->G (41 points)

# -epm : parameters for the empirical pseudopotential method


# -saveband : save electronic bandstructure
# -saveeps : save permittivity
# -savepol : save optical bandstructure

# -switch : band,eps,pol


# -name : project name

# -a : lattice constant [ang]
# -a1/a2/a3 : hmmf?


# -q 0.1,0.5,0.23 \


# -h, -help : help
# -v, -version : version


# -return report?
# -display input parameters (for debugging purposes)

# -print {cmd},file


KGRID="-gridscheme 0 -gridorder 2"
# KGRID="-gridscheme 0 -gridorder 2"
BAND="-numband 8"
RECVECTOR="-b1 1,1,-1 -b2 1,-1,1 -b3 -1,1,1"
MATERIAL="-a 5.43 -epm -0.257,0.0,-0.040,0.033,0,0,0,0,0.55,0.32,0,1.06,0"
QVECTOR="-numq 10,10,5,5,10 -pathq L,G,X,W,K,G"
FREQ="-numfreq 1 -dfreq 1"



# /home/mjh92/projects/picomax_mod/src/a.out -switch eps \
#     -debug 4 \
#     -numpw 124 \
#     -numband 8 \
#     -gridscheme 0 \
#     -gridorder 2 \
#     -b1 1,1,-1 \
#     -b2 1,-1,1 \
#     -b3 -1,1,1 \
#     -a 5.43 \
#     -epm -0.257,0.0,-0.040,0.033,0,0,0,0,0.55,0.32,0,1.06,0 \
#     -pathq L,G,X,W,K,G \
#     -numq 10,10,5,5,10 \
#     -numfreq 1 \
#     -dfreq 1

# /home/mjh92/projects/picomax_mod/src/picomax -switch eps \
#     -debug 0 \
#     -numpw 124 \
#     -numband 8 \
#     -gridscheme 0 \
#     -gridorder 1 \
#     -b1 -1,1,1 \
#     -b2 1,-1,1 \
#     -b3 1,1,-1 \
#     -a 5.43 \
#     -epm -0.257,0.0,-0.040,0.033,0,0,0,0,0.55,0.32,0,1.06,0 \
#     -pathq L,G,X,W,K,G \
#     -numq 10,10,5,5,10 \
#     -numfreq 201 \
#     -dfreq 0.05



# /home/mjh92/projects/picomax_mod/src/picomax -switch eps \
#     -debug 0 \
#     -numpw 124 \
#     -numband 8 \
#     -gridscheme 0 \
#     -gridorder 1 \
#     -b1 -1,1,1 \
#     -b2 1,-1,1 \
#     -b3 1,1,-1 \
#     -a 5.43 \
#     -epm -0.257,0.0,-0.040,0.033,0,0,0,0,0.55,0.32,0,1.06,0 \
#     -q 0.01,0,0 \
#     -numfreq 1 \
#     -dfreq 0.05

# /home/mjh92/projects/picomax_mod/src/picomax -switch eps \
#     -debug 0 \
#     -numpw 124 \
#     -numband 8 \
#     -gridscheme 0 \
#     -gridorder 3 \
#     -b1 -1,1,1 \
#     -b2 1,-1,1 \
#     -b3 1,1,-1 \
#     -a 5.43 \
#     -epm -0.257,-0.040,0.033,0,0,0,0.55,0.32,0,1.06,0 \
#     -q 0.001,0,0 \
#     -numfreq 1 \
#     -dfreq 0.05


# 728

/home/mjh92/projects/picomax_try/src/picomax -switch eps \
    -debug 1 \
    -GCUT 7 \
    -nband 20 \
    -neps 9 \
    -a 4.35 \
    -epm -0.4280,0.1017,0.1108,0.0010,0.0800,0.0277 \
    -q 0.5,0,0 \
    -nfreq 1 \
    -dfreq 0.05 \
    -kk 0 \
    -gridscheme 1 \
    -kgridfile /home/mjh92/projects/picomax_try/KPOINTS_sample
    # -numq 10,10,5,5,10 \
    # -pathq L,G,X,W,K,G



    # -kgridfile /home/mjh92/projects/picomax_mod/input/FCC_FBZ/KPOINTS_10


    
# /home/mjh92/projects/picomax_try/src/picomax -switch eps \
#     -debug 1 \
#     -GCUT 7 \
#     -nband 20 \
#     -neps 9 \
#     -a 4.35 \
#     -epm -0.4280,0.1017,0.1108,0.0010,0.0800,0.0277 \
#     -q 0.5,0,0 \
#     -nfreq 1 \
#     -dfreq 0.05 \
#     -kk 0 \
#     -gridscheme 0 \
#     -gridorder 1


  


    # 100 -> 0.948151



    # FBZ_10 -> 0.987845
    # FBZ_15 -> 1.0071
    # FBZ_20 -> 0.983902
    # IBZ_20 -> 1.11323



    # FCC_FBZ/10 -> 0.987795


picomax 

-debug 3 # debug level

-wdir /tmp/picomax/ # working directory
-outputfile filename # writes outputfile "filename.dat" in the working directory
-logfile 1  # writes logfile "filename.log" in the working directory

-kpointfile ... # just use kpointfile if kpointfile is provided
-kpoint 0 -kpointorder 11 # 11x11x11 Full BZ
-kpoint 1 -kpointorder 11 # 11x11x11 Monkhorst-Pack grid
-kpoint 2 -kpointorder 3 # Cohen-Chadi grid with order 3
# option to directly provide kpoints? hmmf...


-switch eps # calculates dielectric matrix
-switch el # calculates electronic bands
-encut 400 # planewave grids energy cutoff [eV]
-neps 9 # size of dielectric matrix

-qvec 0,0,0 # q-vector [2*pi/a]
-qvec 0,0,0,1,0,0 # q-vectors {q1x,q1y,q1z,q2x,q2y,q2z,...}
-qnum 10,10,5,5,10
-qpath L,G,X,W,K,G

-nfreq 1001 # number of frequency
-dfreq 0.05 # frequency interval [eV]
-kk 1   # kramers-kronig transformation for real parts of dielectric constants
-epsilon 0.1 # delta function
-delta 1 # Lorentzian or gaussian?

-a 4.5  # lattice constant [angstrom]
-epm 1,2,3,4,5,6 # empirical pseudopotential form factors [rydberg]
-nband 20 # number of electronic bands
-originshift 0 # 


-gammazero 1e-6,0,0 # shift for the gamma-point to avoid singularity





-crystal fcc # will other crystal structures be supported?


## subroutine
* generate q-vectors
* assign symmetries to the q-vector paths
    0 -> none
    1 -> [100] C4v
    2 -> [111] C3v
    3 -> [110] C2v
    4 -> Rz


* how to print progress? - log?




#### Demo 1
# picomax -help
# picomax -h


# picomax -version
# picomax -v

####
# picomax -debug 3 -a 4.5 -epm 1,2,3,4,5,6 -kpointfile /home/mjh92/projects/picomax/include/kpoints/KPOINT_5x5x5.txt -encut 200 -nband 20 -nfreq 101 -dfreq 0.1