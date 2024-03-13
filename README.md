# picomax


## Todo:
* add OUTPUT path arg
    picomax ... -output 1 -filename .../output.dat
* add ECUT arg
    picomax ... -ecut 100 ...
* add KPOINT command input arg
    picomax ... -kpoint 100,100,0,0,0,1,0,0,1,1,...
* add input parameter pattern file?
    picomax ... -inputfile .../inputparam.txt
* add logfile output
    picomax ... -logfile true? (check other software languages)
* add 







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