/*
    epm.cpp | includes functions related to calculation of electronic bands

    Author: Jungho Mun
    Date: June 6, 2024
*/

/*
    header files
*/
#include "pmx_epm.hpp"


/*
    function definitions
*/
namespace pmx{


/*
    Returns the Hamiltonian matrix for the empirical pseudopotential method

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions 
    [2] (Cohen and Bergstresser) "Band Structures and Pseudopotnetial Form Factors 
        for Fourteen Semiconductors of the Diamond and Zinc-blende Structures", 
        Phys. Rev. 141, 789 (1966)
    [3] (Bergstresser and Cohen) "Electronic Structure and Optical Properties of 
        Hexagonal CdSe, CdS, and ZnS", Phys. Rev. 164, 1069 (1967)
*/
Eigen::MatrixXcd HamiltonianEPM (std::vector<Eigen::Vector3d> G, Eigen::Vector3d K, std::vector<Eigen::Vector3d> T, mater mat){
    Eigen::MatrixXcd H(NPW, NPW);
    double KINETIC_CONST = KCON*pow(2.0*pi/a,2);

    for (int i=0; i<NPW; i++){ for (int j=0; j<=i; j++){
        if (i==j){ // diagonal kinetic part
            H(i,j) = KINETIC_CONST * (K+G[i]).squaredNorm();
        }else{ // local pseudopotential
            Eigen::Vector3d dG = G[i]-G[j];
            if (dG.squaredNorm()<mat.g2max)
                H(i,j) = pseudopotential(dG,T,mat);
            else
                H(i,j) = 0;
        }
    }}
    // only lower triangle is needed for eigenvalue calculation (j<=i)
    return H;
}


/*
    Returns the complex pseudopotential Vg, at the given reciprocal vector |G|

    Vg = \sum_{i=1}^{NATOM} V_{ig} exp(1i*dot(t_i,g)) / NATOM

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions 
    [2] (Cohen and Bergstresser) "Band Structures and Pseudopotnetial Form Factors 
        for Fourteen Semiconductors of the Diamond and Zinc-blende Structures", 
        Phys. Rev. 141, 789 (1966)
    [3] (Bergstresser and Cohen) "Electronic Structure and Optical Properties of 
        Hexagonal CdSe, CdS, and ZnS", Phys. Rev. 164, 1069 (1967)
*/

double a1 = -0.8175 ;
double a2 = 1.117;
double a3 = 1.056 ;
double a4 = 0.3885; 
double a5 = -0.4087; 
double a6 = 0.6047; 
double a7 = 0.4071; 
double a8 = -0.4445;
double b1 = 0.5456; 
double b2 = -0.09914; 
double b3 = -0.5045; 
double b4 = 0.8948; 
double b5 = 1.032;
double b6 = 0.9786; 
double b7 = -0.2666;
double b8 = 0.675; 
double c1 = 0.9065; 
double c2 = 1.083 ;
double c3 = -0.1143 ;
double c4 = -0.0653; 
double c5 = 1.126; 
double c6 = -0.8528; 
double c7 = 1.199; 
double c8 = 0.773 ;
double d1 = 0.6527 ;
double e1 = 0.8691;

std::complex<double> pseudopotential(
        Eigen::Vector3d G, std::vector<Eigen::Vector3d> T, mater mat){
    
    double GMAG = G.squaredNorm();
    
    std::complex<double> Vm = (a1*sin(b1*sqrt(GMAG) +c1)+a2*sin(b2*sqrt(GMAG) +c2)+a3*sin(b3*sqrt(GMAG) +c3)+a4*sin(b4*sqrt(GMAG) +c4)+a5*sin(b5*sqrt(GMAG) +c5)+a6*sin(b6*sqrt(GMAG) +c6)+a7*sin(b7*sqrt(GMAG) +c7)+a8*sin(b8*sqrt(GMAG) +c8))*exp(-d1*pow(sqrt(GMAG),e1));
    
    std::complex<double> Vg = 0;
    for (int i=0; i<NATOM; i++){
        Vg += exp(-im*2.0*pi*G.dot(T[i]));
    }
    Vg *= Vm*Ry2eV/double(NATOM);
    return Vg;

}


/*
    Calculates electron energy band using the empirical pseudopotential method

    todo: change function name

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions 
    [2] (Cohen and Bergstresser) "Band Structures and Pseudopotnetial Form Factors 
        for Fourteen Semiconductors of the Diamond and Zinc-blende Structures", 
        Phys. Rev. 141, 789 (1966)
    [3] (Bergstresser and Cohen) "Electronic Structure and Optical Properties of 
        Hexagonal CdSe, CdS, and ZnS", Phys. Rev. 164, 1069 (1967)
*/


void empiricalpseudopotentialmethod(env &dat){
    dat.eband = new double*[NQ];

    // obtain energy offset
    setRefEnergy(dat);

    #pragma omp parallel for
    for (int q=0; q<NQ; q++){
        
        // Hamiltonian
        Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,dat.lat.Q[q],dat.lat.atomic,dat.mat);

        // diagonalize the Hamiltonian
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver(H);
        if (eigsolver.info() != Eigen::Success){ std::cout << "ERROR: Eigensolver failed to solve the matrix." << std::endl; abort(); }

        // truncate to NBAND and energy offset
        dat.eband[q] = new double[NBAND];
        for (int n=0; n<NBAND; n++){
            dat.eband[q][n] = double(eigsolver.eigenvalues()(n))-dat.energyoffset;
        }
    }
    #pragma omp barrier

    return;
}


/*
    Returns the nonlocal pseudopotential Vg^NL

    Gi, Gj : reciprocal lattice vectors
    T : atomic basis vectors
    mat : material parameters

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions 
    [2] 
*/



// complex nonlocal-pseudopotential for the given reciprocal vector
std::complex<double> nonlocalpseudopotential(Eigen::Vector3d Gm, Eigen::Vector3d Gn, 
                    std::vector<double> v){
    double a0 = v[6];
    double b0 = v[7];
    double rl0 = v[9];

    //define variables for symmetric and asymmetric part of the potential
	double sym = 0.0e0, asym = 0.0e0;
	std::complex<double> potential = 0.0e0;

	//Each cell has 2 atoms.  Use mid point between them as offset in units of $a$
	Eigen::Vector3d vec_tau;
    vec_tau << 0.125, 0.125, 0.125;

	//Find norm of the reciprocal vector
	double Km = Gm.norm()*(2.0*pi/a);
	double Kn = Gn.norm()*(2.0*pi/a);
	double kfermi = pow(6.0*pi*pi/(pow(a,3)), double(1.0/3.0));
	double a0E = a0*Ry2eV + (b0*Ry2eV*KCON*(sqrt(Km*Km*Kn*Kn)-pow(kfermi,2)));
	double Fl;
	double Kmr = Km*rl0;
	double Knr = Kn*rl0;

	if(abs(Kmr-Knr) < 1e-9){
		Fl = 0.5*pow(rl0,3) * (pow(sphbesj(0,Kmr),2) - sphbesj(-1,Kmr)*sphbesj(1,Kmr));
	}else{
		Fl = (pow(rl0,2)/(Km*Km-Kn*Kn)) * (Km*sphbesj(1,Kmr)*sphbesj(0,Knr) - Kn*sphbesj(1,Knr)*sphbesj(0,Kmr));
	}	
	sym = (4.0*pi/pow(a,3))*a0E*Fl;

	//Complex Pseudopotential value
	potential = std::complex<double>(sym*cos(2*pi*((Gm-Gn).dot(vec_tau))), 
                                    asym*sin(2*pi*((Gm-Gn).dot(vec_tau))));
	return potential;
}

/*
    Returns the spherical Bessel function of the first kind
*/
double sphbesj(int n, double r){
	double jn = 0.0e0;
	if(abs(r) < 1e-8){
		switch(n){
			case -1:
				jn = 1.0e8;	break;
			case 0: 
				jn = 1.0e0;	break;
			case 1: 
				jn = 0.0e0;	break;
		}
	}
	else{
		switch(n){
			case -1:
				jn = cos(r)/r;	break;
			case 0: 
				jn = sin(r)/r;	break;
			case 1: 
				jn = (sin(r)/(r*r))-(cos(r)/r);	break;
		}
	}
	return jn;
}

/*
    set reference energy [eV] at the reference kpoint

    default reference point is the gamma point

    you can change the reference point using the option `-refpoint x,y,z`
    ex) `picomax ... -refpoint 0.5,0.5,0.5`
*/
void setRefEnergy(env &dat){

    // the global Hamiltonian matrix H at gamma point
    Eigen::Vector3d K_REF = dat.refpoint;// reference kpoint, where the bandgap occurs

    Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K_REF,dat.lat.atomic,dat.mat);

    // compute the eigenvalue problem
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver(H);
    if (eigsolver.info() != Eigen::Success){ std::cout << "ERROR::setRefEnergy:: Eigensolver failed to solve the matrix at Gamma point." << std::endl; abort(); }

    // zero level of energy at the valence band maximum at reference point
    // dat.energyoffset = -1.6;
    dat.energyoffset = eigsolver.eigenvalues()(dat.nvalence-1);
    // dat.energyoffset = (eigsolver.eigenvalues()(dat.nvalence)+eigsolver.eigenvalues()(dat.nvalence-1))/2.0;

    return;
}

}/* namespace pmx */
















// // hamiltonian matrix for EPM
// Eigen::MatrixXcd HamiltonianEPM (std::vector<Eigen::Vector3d> G, Eigen::Vector3d K, std::vector<double> v){
//     Eigen::MatrixXcd H(NPW, NPW);
//     double KINETIC_CONST = KCON*pow(2.0*pi/a,2);

//     for (int i=0; i<NPW; i++){ for (int j=0; j<=i; j++){
//         if (i==j){ // diagonal kinetic part
//             H(i,j) = KINETIC_CONST * (K+G[i]).squaredNorm();
//         }else{ // local pseudopotential
//             Eigen::Vector3d dG = G[i]-G[j];
//             if (dG.squaredNorm()<11.1)
//                 H(i,j) = pseudopotential(dG,v);
//             else
//                 H(i,j) = 0;
//         }
//         if (EPM_NONLOCAL){ // nonlocal pseudopotential
//             H(i,j) += nonlocalpseudopotential(K+G[i],K+G[j],v);
//         }
//     }}
//     // only lower triangle is needed for eigenvalue calculation
//     return H;
// }



// // complex pseudopotential for the given reciprocal vector
// std::complex<double> pseudopotential(Eigen::Vector3d G, std::vector<double> v){
//     //define variables for symmetric and asymmetric part of the potential
// 	double vs = 0.0e0, va= 0.0e0;
// 	std::complex<double> potential = 0.0e0;

// 	//Each cell has 2 atoms.  Use mid point between them as offset
// 	Eigen::Vector3d vec_tau;
//     vec_tau << 0.125, 0.125, 0.125;

// 	//Find the magnitude of the reciprocal vector
//     double GMAG = G.squaredNorm();

// 	//assign symmetric and asymmetric part of the potential 
// 	if (abs(double(GMAG - 3.0)) < 1e-6){ vs = v[0]; va = v[3]; }
// 	else if (abs(double(GMAG - 4.0)) < 1e-6){ vs = 0.0; va = v[4]; }
// 	else if (abs(double(GMAG - 8.0)) < 1e-6){ vs = v[1]; va = 0.0; }
// 	else if (abs(double(GMAG - 11.0)) < 1e-6){ vs = v[2]; va = v[5]; }
// 	else { vs = 0.0e0; va = 0.0e0; }

// 	//Complex Pseudopotential value, potentials are scaled to rydberg units
//     if (ORIGINSHIFT){
//         potential = std::complex<double>(
//                         vs*cos(2.0*pi*(G.dot(vec_tau)))*Ry2eV, 
//                         va*sin(2.0*pi*(G.dot(vec_tau)))*Ry2eV );
//     }else{
//         double V1, V2;
//         V1 = (vs+va)/2;
//         V2 = (vs-va)/2;
//         potential = std::complex<double>(
//                         (V1 + V2*cos(-4.0*pi*(G.dot(vec_tau))))*Ry2eV,
//                         (V2*sin(-4.0*pi*(G.dot(vec_tau))))*Ry2eV );
//     }
// 	return potential;
// }






    // //define variables for symmetric and asymmetric part of the potential
	// double vs = 0.0e0, va = 0.0e0;
	// std::complex<double> potential = 0.0e0;

	// //Each cell has 2 atoms.  Use mid point between them as offset
	// Eigen::Vector3d vec_tau;
    // vec_tau << 0.125, 0.125, 0.125;

	// //Find the magnitude of the reciprocal vector
    // double GMAG = G.squaredNorm();

	// //assign symmetric and asymmetric part of the potential 
	// if (abs(double(GMAG - 3.0)) < 1e-6){ vs = v[0]; va = v[3]; }
	// else if (abs(double(GMAG - 4.0)) < 1e-6){ vs = 0.0; va = v[4]; }
	// else if (abs(double(GMAG - 8.0)) < 1e-6){ vs = v[1]; va = 0.0; }
	// else if (abs(double(GMAG - 11.0)) < 1e-6){ vs = v[2]; va = v[5]; }
	// else { vs = 0.0e0; va = 0.0e0; }

	// //Complex Pseudopotential value, potentials are scaled to rydberg units
    // if (ORIGINSHIFT){
    //     potential = std::complex<double>(
    //                     vs*cos(2.0*pi*(G.dot(vec_tau)))*Ry2eV, 
    //                     va*sin(2.0*pi*(G.dot(vec_tau)))*Ry2eV );
    // }else{
    //     double V1, V2;
    //     V1 = (vs+va)/2;
    //     V2 = (vs-va)/2;
    //     potential = std::complex<double>(
    //                     (V1 + V2*cos(-4.0*pi*(G.dot(vec_tau))))*Ry2eV,
    //                     (V2*sin(-4.0*pi*(G.dot(vec_tau))))*Ry2eV );
    // }
	// return potential;