# PicoMax

PicoMax is a C++ code for deep-microscopic optical and dielectric response calculations of
crystalline solids. It solves the electronic band structure with the empirical pseudopotential
method (EPM) and builds the microscopic susceptibility tensor `chi(q, omega)`, including
local-field effects, from the resulting electronic states. From this tensor a range of optical
and spectroscopic quantities can be derived, such as the macroscopic dielectric constant,
energy-loss functions, and plasmon dispersion.

The name reflects the goal of the project: resolving the optical response of a material down to
the picoscale (sub-unit-cell) length, where the response varies within the unit cell and the
local-field (crystal-local-field) effects become important.


## Features

* electron bandstructure calculation
    * empirical pseudopotential method (EPM)
    * spin-orbit coupling (symmetric and antisymmetric form factors)
* dielectric / susceptibility matrix calculation `chi(q, omega)`
    * longitudinal-longitudinal, transverse-transverse, and mixed tensor components
* optical (pico)polarization calculation
* macroscopic dielectric constant
* static and high-frequency (ion-clamped) dielectric constants
* electron energy-loss spectroscopy (EELS)
* energy-loss function (ELF)
* plasmon dispersion
* optical wave dispersion
* dynamic structure factor / inelastic X-ray scattering (IXS)
* absorption, `~Im(eps)`
* crystal structure and symmetry
* twisted bilayer graphene module (experimental)


## Methodology

The calculation proceeds in a few well-defined stages:

1. **Plane-wave basis.** G-vectors are generated and truncated by a cutoff `gcut` that is
   derived from the plane-wave energy cutoff `encut`.
2. **Electronic structure.** For each k-point the EPM Hamiltonian `H(k)` is diagonalized to get
   the eigenvalues `E_n(k)` and eigenvectors `|psi_n(k)>`. Optional spin-orbit coupling is added
   through symmetric and antisymmetric form factors.
3. **Susceptibility tensor.** The response `chi(q, omega)` is built from Fermi's golden rule over
   inter-band transitions (valence to conduction), with energy conservation
   `E_c(k) - E_v(k+q) = hbar*omega`. Overlap integrals of the form
   `<k,c| j exp(i (q+G).r) |k+q,v>` are accumulated and broadened by a delta-function model.
4. **Post-processing.** The real part of the response can be reconstructed with a Kramers-Kronig
   transform of the imaginary part.

The linear response definition used is `P(q, omega) = chi(q, omega) E(q, omega)`. The tensor is
stored per q-vector, per pair of local-field G-vectors `(m, n)`, per Cartesian pair `(i, j)`, and
per frequency point.


## Repository layout

```
PicoMax/
  src/                 main C++ sources (see below)
  include/
    eigen-3.4.0/       bundled Eigen linear algebra library
    kpoints/           precomputed k-point grid files (KPOINTS_NxNxN.txt)
    matlab_wrapper/    MATLAB class wrapper for driving picomax
    python_wrapper/    Python wrapper and demo scripts
  output_files/        example output and plotting scripts
  pmx_dos.cpp          density-of-states helper
  README.md
```

Core source files in `src/`:

| File | Purpose |
| --- | --- |
| `pmx.cpp` | main entry point, global variables, calculation dispatch |
| `pmx.hpp` | main header, physical constants, data structures, input parser |
| `pmx_io.cpp` / `pmx_io.hpp` | initialization, input parsing, output (documentation, version, exporters) |
| `pmx_basis.cpp` / `pmx_basis.hpp` | G-vector, k-vector, and q-vector generation, high-symmetry points |
| `pmx_epm.cpp` / `pmx_epm.hpp` | electronic band structure solver (EPM, spin-orbit coupling) |
| `pmx_chi.cpp` / `pmx_chi.hpp` | susceptibility tensor `chi(q, omega)` calculation |
| `pmx_tbg.cpp` / `pmx_tbg.hpp` | twisted bilayer graphene bands (experimental) |


## Dependencies

* a C++17 compiler with OpenMP support (for example `g++`)
* the Eigen linear algebra library (a copy is bundled in `include/eigen-3.4.0`)

Optional, for the wrappers and planned features:

```bash
# HDF5 output (planned)
sudo apt-get install libhdf5-dev
```


## Building

Compilation is a single `g++` command that links the core sources against Eigen and OpenMP.
Adjust the `BASE` path in `src/compile_cmd.sh` to point at your checkout, then run it:

```bash
cd src
BASE=/path/to/PicoMax
g++ -O3 -fopenmp -ffast-math \
    -I ${BASE}/include/eigen-3.4.0 \
    -o picomax pmx.cpp pmx_io.cpp pmx_chi.cpp pmx_basis.cpp pmx_epm.cpp
```

This produces the `picomax` executable in `src/`.


## Usage

```bash
picomax [-option] [value] ...
```

The kind of calculation is selected with `-switch`:

* `-switch band` : electronic band structure. Results are written to `${outputfile}.dat`.
* `-switch epsij` : susceptibility tensor `chi(q, omega)`. Results are written to
  `${outputfile}.dat`, along with `${outputfile}.gpt` (G-vectors).
* `-switch epsLL` : longitudinal-longitudinal response only.

Utility commands:

```bash
./picomax -v          # print PicoMax version
./picomax -h          # print built-in documentation
./picomax -debug 2 .. # print parsed inputs and stop (dry-run for debugging)
```


## Command-line arguments

General:

`-h`, `-help`: print help

`-v`, `-version`: print PicoMax version

`-debug`: debug level (default 0). A level above 1 prints the parsed inputs and exits.

`-outputfile`: output filename (without extension)

`-wdir`: working directory (default is the present working directory)

`-logfile`: redirect progress output to `${wdir}/${outputfile}.log`

`-switch`: calculation switch (`band`, `epsij`, `epsLL`)

Crystal and material:

`-crystal`: crystal structure (`zb`, `wz`, `rs`, `te`, `hex`)

`-dim`: crystal dimension (default 3)

`-a`: lattice constant [angstrom] (mandatory)

`-f`: lattice parameter (for example the `c/a` ratio for hexagonal cells)

`-u`: atomic basis parameter

`-vg`: empirical pseudopotential form factors

`-gsym`: enforce point-group symmetry on the G-vector grid by truncation

Spin-orbit coupling:

`-us_so`: symmetric spin-orbit form factors, as `G^2:value,G^2:value,...`

`-use_uniform_us_so`: use a single symmetric spin-orbit value for all `G^2` (flag)

`-uniform_us_so_val`: the uniform symmetric spin-orbit value [Rydberg]

`-ua_so`: antisymmetric spin-orbit form factors, as `G^2:value,G^2:value,...`

`-ua_so_val_g0`: antisymmetric spin-orbit form factor at `G^2 = 0`

Plane-wave basis and k-points:

`-encut`: plane-wave energy cutoff [eV]

`-gcut`: G-vector cutoff [2*pi/a]; normally derived from `encut`, but can be set directly

`-kpoint`: k-point grid mode (`0` import from file, `1` full BZ, `2` Monkhorst-Pack, `3` Cohen-Chadi)

`-kpointfile`: k-point file to import when `-kpoint 0`

`-kpointorder`: k-point density `NxNxN`

`-nband`: number of electronic bands

`-nchi`, `-neps`: size of the dielectric (local-field) matrix `NxN`

`-nvalence`: number of valence electron bands

q-vectors:

`-qvec`: an explicit q-vector

`-qpath`: path of high-symmetry points (default `L,G,X,W,K,G`)

`-qnum`: number of q-vectors in each path segment (default `20,20,10,10,20`)

`-refpoint`: reference q-vector used for the electronic energy reference (default `0,0,0`)

`-gammazero`: the Gamma-point (default `0,0,0`); shift it slightly away from zero to avoid the
singularity in `epsLL` calculations

Frequency and broadening:

`-nfreq`: number of frequency points (default 1)

`-dfreq`: frequency interval [eV]

`-delta`: delta-function model (`0` Lorentzian, `1` Gaussian)

`-epsilon`: delta-function smearing [eV]

`-kk`: switch for the Kramers-Kronig transform


## Crystal structures

The following structures set up their lattice, reciprocal, and atomic basis vectors automatically
from `-a` (and `-f`, `-u` where relevant):

* `zb`: zinc-blende / sphalerite (2 atoms)
* `rs`: rock-salt / halite (2 atoms)
* `wz`: wurtzite (4 atoms)
* `te`: tellurium (3 atoms)
* `hex`: 2D hexagonal (2 atoms)


## K-point grids

The k-point grid used for Brillouin-zone integration is chosen with `-kpoint`:

* `0`: import from a k-point file (`-kpointfile`). Ready-made grids are provided in
  `include/kpoints/` as `KPOINTS_NxNxN.txt`.
* `1`: generate the full Brillouin zone at density `-kpointorder`.
* `2`: generate a Monkhorst-Pack grid.
* `3`: generate a Cohen-Chadi grid.


## Output files

All output is written to the working directory:

* `${outputfile}.dat`: the main binary result. For `-switch band` it holds the band energies as
  `[NQ][NBAND]` doubles. For `-switch epsij` it holds the real part followed by the imaginary part
  of the susceptibility tensor, indexed as `[i][j][q][m][n>=n][f]`.
* `${outputfile}.gpt`: the local-field G-vectors, as `[NEPS][3]` doubles.
* `${outputfile}.kpt`: the k-vectors used in the calculation.
* `${outputfile}.log`: the progress log, when `-logfile` is set.

Binary output is read back conveniently through the MATLAB and Python wrappers.


## Wrappers

### MATLAB

`include/matlab_wrapper/` provides a `picomax` class that assembles the command line, launches the
executable (locally, or over WSL from Windows), and reads the binary output back into MATLAB.

```matlab
c = picomax;
c.setopt('encut', 300, 'nband', 20, 'nvalence', 4, 'qpath', 'L,G,X,U,K,G', 'qnum', [20 20 10 0 20]);
c.set('crystal', 'zincblende', 'a', 4.35, 'vg', [3 4 8 11 -0.225 0.005 0.118 0.127 -0.545 -0.005 0.118 0.121]);
c.setopt('picomax_exe', '/path/to/picomax');
c.run('eband');   % solve for the band structure
```

See `include/matlab_wrapper/README.md` for WSL setup and further examples.

### Python

`include/python_wrapper/` provides `picomax.py` and demo scripts in `demo/`
(`demo_eband.py`, `demo_optical_band.py`, `demo_refractive_index.py`). See
`include/python_wrapper/README.md` for environment setup.


## Example materials

Silicon:

```
-a 5.43 -vg 3,4,8,11,-0.21,0,0.04,0.08,-0.21,0,0.04,0.08
```

Silicon carbide (SiC):

```
-a 4.35 -vg 3,4,8,11,-0.225,0.005,0.118,0.127,-0.545,-0.005,0.118,0.121
```

Tellurium (Te), with spin-orbit coupling, is driven through longer form-factor lists; see
`src/code.sh` for a complete, ready-to-run example.


## Appendix A. Sample scripts

#### Sample script for running picomax calculations

```bash
#sample_script.sh
#!/bin/bash
outputfile=outputfilename
wdir=/path/to/demo
picomax_exe=/path/to/picomax
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

${picomax_exe} -switch epsij \
    -outputfile ${outputfile} -wdir ${wdir} \
    -nband ${nband} -neps ${nchi} -encut ${encut} \
    -crystal ${crystal} -a ${a} -vg ${vg} \
    -kpoint 1 -kpointorder 7 \
    -nfreq ${nfreq} -dfreq ${dfreq} -kk ${kk} -delta ${delta} -epsilon ${epsilon} \
    -qnum ${qnum} -qpath ${qpath} \
    >> ${wdir}/${outputfile}.log
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

outputfile=pmx_si_freq
wdir=/path/to/demo
picomax_exe=/path/to/picomax
kpointfile=/path/to/include/kpoints/KPOINTS_11x11x11.txt
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

srun --unbuffered \
    ${picomax_exe} -switch epsij \
    -outputfile ${outputfile} -wdir ${wdir} \
    -nband ${nband} -neps ${nchi} -encut ${encut} \
    -a ${a} \
    -kpoint 0 -kpointfile ${kpointfile} \
    -nfreq ${nfreq} -dfreq ${dfreq} -kk ${kk} -delta ${delta} -epsilon ${epsilon} \
    -qnum ${qnum} -qpath ${qpath} \
    >> ${wdir}/${outputfile}.log
```

#### Progress output example

```
*********************************
***PicoMax 0.1 progress output***
*********************************
Date
PicoMax v0.1
  Running on ...
  Using N cores
  Total memory: ... MB
  Avail memory: ... MB
  Working dir: ...
  Output file: ....dat
---------------------------------
Parameters:
  encut = ... eV
  gcut = ...
  nchi = 9
  nband = 20
  nfreq = 1
  dfreq = 0.1
  epsilon = 0.3
---------------------------------
Generated Q-vector
Generated G-vectors
Generated K-vectors
Susceptibility tensor matrix solver
  Elapsed time: ... s
```


## Roadmap

* logfile output polish
* linearize the `(m, n)` local-field dimensions
* HDF5 output support (important for reusing results)
* automatic setup of the potential cutoff `vcut` from the maximum `|g|^2` value
* effective mass and optical phonon modes
* estimation of memory and runtime before a run
* multi-node (cluster) calculations with MPI
* Wannier functions
* integration with VASP / Quantum Espresso for general materials
* analysis of finite nanostructures (slabs, nanoparticles, nanorods)


## Authors

* Jungho Mun
* Pronoy Das
