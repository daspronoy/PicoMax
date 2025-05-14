# picomax



## Completed tasks:
* add output path arg
    picomax ... -output 1 -filename .../output.dat
* add encut arg
    picomax ... -encut 100 ...
* hardware output
* 

## Todo:
* add logfile output
    `picomax ... -logfile true` ? (check other software)
* linearize (m,n) dimensions
* hdf5 output support (this will be important for reusing results)
* `-crystal zb` -> a1,a2,a3,b1,b2,b3,t1,t2
* `-crystal wz`
* -crystal 2d??
* `-bvector -1,1,1,1,-1,1,1,1,-1` (zb)
* `-avector 0,0.5,0.5,0.5,0,0.5,0.5,0.5,0` (zb)
* `-atomicbasis 0,0,0,0.25,0.25,0.25` (zb)
* `-Vformfac 3,4,8,11,v1_3,v1_4,v1_8,...` (zb)
*   automatically set `vcut` from the maximum `|g|^2` value
* `-vtol 0.1`


Silicon:
`-a 5.43 -vg 3,4,8,11,-0.21,0,0.04,0.08,-0.21,0,0.04,0.08`
SiC:
`-a 4.35 -vg 3,4,8,11,-0.225,0.005,0.118,0.127,-0.545,-0.005,0.118,0.121`


`-a 4.35 -vg 3,4,8,11,12,16,-0.166,-0.020,0.130,0.002,0,0,-0.450,-0.412,-0.110,0.128,0.182,0.130`


## 
sudo apt-get install libhdf5-dev






## Features
* electron bandstructure calculation
    EPM
    spin-orbit?
* dielectric matrix calculation
    longitudinal/transverse/all?
* optical (pico)polarization calculation
* macroscopic dielectric constant
* static dielectric constant
    static/high-frequency (ion-clamped) dielectric constants
* electron energy-loss spectroscopy (EELS)
* energy-loss function (ELF)
* plasmon dispersion
* optical wave dispersion
* dynamic structure factor
    inelastic X-ray scattering (IXS)
* absorption
    ~Im(eps)
* crystal structure and symmetry
* effective mass?
* (optical) phonon modes???


## Subroutines
* g-vector grids
* kpoint grids (BZ integration grid)



## Future plans
* Estimation of memory requirements
* Estimation of calculation time
* Multi-node calculation (cluster), MPI
* leave a log file (on/off)
* Wannier function
* integration with VASP/Quantum Espresso/... for general materials
* Monkhorst-Pack grids, Cohen-Chadi grids, ...
* topological properties?
* Matlab/python wrapper?
* socket programming for full integration?


## Hmmf
* analysis for finite nanostructures? (slabs, nanoparticles, nanorods, ...)
* interaction with light with OAM?
* extension to nonlinearity?




## Engineering perspectives?
* Noninvasive screening at EUV/DUV regimes?


## Relevant papers?

* Probing correlated states with plasmons





## Meaning of input arguments

`-switch`:

`-outputfile`: 

`-wdir`: working directory

`-nband`: the number of electron bands

`-nchi`: the size of dielectric matrix

`-encut`: the energy cutoff for g-vectors

`-nfreq`: the number of frequency points

`-dfreq`: the frequency interval [eV]

`-kk`: switch for Kramers-Kronig transform

`-delta`: model of delta function

`-epsilon`: the smearing parameter of delta function [eV]

`-qpath`: 

`-qnum`: 

`-qvec`: 

`-a`:

`-f`:

`-u`:

`-vg`:


















## Appendix A. Sample scripts

#### Sample script for running picomax calculations
```bash
#runscript.sh
#!/bin/bash

outputfile=filename
picomax_exe=
wdir=


```


#### Sample script for slurm submission

```bash
#slurm_cmd.sh
#!/bin/bash
#SBATCH --partition=hcpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 48 
#SBATCH --exclusive
#SBATCH --job-name=picomax
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.log

outputfile=pmx20240328_si_freq
wdir=/home/mjh92/picomax/picomax/demo
picomax_exe=/home/mjh92/project/picomax/src/picomax
kpointfile=/home/mjh92/project/picomax/include/kpoint/KPOINT_11x11x11.txt
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
a=4.35
epm=

srun --unbuffered\
    ${picomax_exe} -switch epsij\
    -outputfile ${outputfile} -wdir ${wdir}\
    -nband ${nband} -neps ${nchi} -encut ${encut}\
    -a ${a} -epm ${epm}\
    -kpoint 0 -kpointfile ${kpointfile}\
    -nfreq ${nfreq} -dfreq ${dfreq} -kk ${kk} -delta ${delta} -epsilon ${epsilon}\
    -qnum ${qnum} -qpath ${qpath}\
    >> ${wdir}/${outputfile}.log
```

#### Print out picomax version

```bash
./picomax -v
```

#### Print out picomax doc

```bash
./picomax -h
```

#### Print out picomax inputs (for debugging purposes)

```bash
./picomax -debug 1 ... # dry-run??
```






#### progress output file example

```
*********************************
***PicoMax 0.1 progress output***
*********************************
Date
PicoMax 0.1 starting
The output filename is ...

number of cores? spec of CPU?, hostname
Running on ...
Available memory: 30 GB
Parameters:
encut=** eV
nchi=9 (neps?)
nband=20
nfreq=1
dfreq=0.1

Generated G-vectors
Elapsed time: * s
Number of planewaves: *

Generated K-vectors
Simulation size in kpoint: 15 x 15 x 15
(Imported K-vectors from kpointfile: )
Number of K-vectors: 
Elapsed time: * s

Generated Q-vectors
Number of Q-vectors:
Path
Elapsed time: * s

Electronic bandstructure solver
...

Susceptibility matrix solver
Solving for only the longitudinal-longitudinal components
Solving for the full-tensor components

Progress 
0% complete. Estimated remaining time: +++ s
10% complete. Estimated remaining time: +++ s
Completed!
Elapsed time: * s

Finished exporting data to /path/to/file
Total wall time: * s
```





