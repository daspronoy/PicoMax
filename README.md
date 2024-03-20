# picomax



## Completed tasks:
* add output path arg
    picomax ... -output 1 -filename .../output.dat
* add encut arg
    picomax ... -encut 100 ...



## Todo:
* add KPOINT command input arg
    picomax ... -kpoint 100,100,0,0,0,1,0,0,1,1,...
* add input parameter pattern file?
    picomax ... -inputfile .../inputparam.txt
* add logfile output
    picomax ... -logfile true? (check other software languages)
* add 
* linearize (m,n) dimensions

#### progress output file example

-----

**************************************
***PicoMax 0.1 progress output file***
**************************************
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
-----


## Appendix A. Sample scripts

#### PicoMax script
---
#picomax_cmd.sh
#!/bin/bash
outputfile=filename

picomax_exe=
wdir=


---

#### Slurm script
---
#slurm_cmd.sh
#!/bin/bash
#SBATCH --partition=hcpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 48 
#SBATCH --exclusive
#SBATCH --job-name=filename
#SBATCH --output=%x.log
#SBATCH --error=%x.log
outputfile=%x
wdir=
picomax_exe=
kpointfile=
srun --unbuffered\
    ${picomax_exe} -switch epsij\
    -outputfile ${outputfile} -wdir ${wdir}\
    -nband 20 -neps 9 -encut 400\
    -a 4.35 -epm \
    -kpoint 0 -kpointfile ${kpointfile}\
    -nfreq 1 -dfreq 0.05 -kk 0 -delta 0 -epsilon 0.3\
    -qnum 10,10,5,5,10 -qpath L,G,X,W,K,G
---

#### dd
---
./picomax -v
---

---
./picomax -h
---

---
./picomax -debug 1 ... # dry-run??
---






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





* Progress bar
    Starting computation...
    5% (100 s)
    10% (200 s)
    ...
    100% (1000 s)
    Complete

* Estimation of memory requirements
* Estimation of calculation time
* Multi-node calculation (cluster), MPI
* leave a log file (on/off)
* specify output file
    picomax ... -output asdf.dat

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








## GIT terminal commands:
git commit -a -m "message"
git push