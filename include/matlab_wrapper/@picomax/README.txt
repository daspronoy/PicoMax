

----------
Matlab wrapper






#### Miscellaneous
picomax.version
picomax.help
picomax.set
picomax.setopt


#### Solvers
picomax.electronband
picomax.permittivity
picomax.polarizationband
picomax.polarizationdist



#### Subroutines
picomax.generate_qvector
picomax.highsym_fcc
picomax.highsym_hex
picomax.highsym


export?
plot?




#### Properties
picomax.cache
picomax.lattice {sym, lattice vectors, lattice constant, ...}
	sym
	lattice_vec
	rec_vec
	qvector
	qpath
	qnum
	dim

picomax.numeric
	epsilon
	nband
	nfreq
	dfreq
	nq
	neps

picomax.result
	electronband
	electrondos
	omega [eV]
	freq [THz]
	epsilonL
	epsilonT
	alphaL
	alphaT
	qvector

picomax.unit
	kB
	...




//

optical band dispersions?
-> macroscopic local
-> macroscopic nonlocal
-> Maxwell Hamiltonian

electron energy loss function
-> 1/Im(epsL)
with/without loca-field effects

plasmon dispersion
det(epsL)=0 ???


symmetry...
irreducible BZ?
truncate KPOINTS?






kpoint grids
freq

program -> path, ...



pm = picomax;
pm.set('ECUT',100
'NBAND'
'CRYSTAL'
'a'

'NFREQ'
'DFREQ'
'KK'
'EPSILON'
'QVEC'
'NEPS'
'EPM_PARAM'

'QPATH'
'QNUM'



Monkhorst-Pack
Chadi-Cohen



