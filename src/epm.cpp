/*
    epm.cpp
    contains functions related to calculation of electronic bands


    Author: 
*/
#include "pmx.hpp"

namespace pmx{

// empirical pseudopotential method
void empiricalpseudopotentialmethod(env &dat){
    dat.band = new double*[NQ];

    // obtain energy offset
    setRefEnergy(dat);

    #pragma omp parallel for
    for (int q=0; q<NQ; q++){
        
        // Hamiltonian
        Eigen::MatrixXcd H = HamiltonianEPM(dat.G,dat.Q[q],dat.epm);

        // diagonalize the Hamiltonian
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver(H);
        if (eigsolver.info() != Eigen::Success){ std::cout << "ERROR: Eigensolver failed to solve the matrix." << std::endl; abort(); }

        // truncate to NBAND and energy offset
        dat.band[q] = new double[NBAND];
        for (int n=0; n<NBAND; n++){
            dat.band[q][n] = double(eigsolver.eigenvalues()(n))-dat.energyoffset;
        }
    }
    #pragma omp barrier

    return;
}






// hamiltonian matrix for EPM
Eigen::MatrixXcd HamiltonianEPM (std::vector<Eigen::Vector3d> G, Eigen::Vector3d K, std::vector<double> v){
    Eigen::MatrixXcd H(NPW, NPW);
    double KINETIC_CONST = KCON*pow(2.0*pi/a,2);

    for (int i=0; i<NPW; i++){ for (int j=0; j<=i; j++){
        if (i==j){ // diagonal kinetic part
            H(i,j) = KINETIC_CONST * ((K+G[i]).squaredNorm());
        }else{ // local pseudopotential
            Eigen::Vector3d dG = G[i]-G[j]; // TEST Gj-Gi or Gi-Gj
            if (dG.norm()<3.32) // dG.squaredNorm()<11.1
                H(i,j) = pseudopotential(dG,v);
            else
                H(i,j) = 0;
        }
        if (EPM_NONLOCAL){ // nonlocal pseudopotential
            H(i,j) += nonlocalpseudopotential(K+G[i],K+G[j],v);
        }
    }}
    // only lower triangle is needed for eigenvalue calculation
    return H;
}



// complex pseudopotential for the given reciprocal vector
std::complex<double> pseudopotential(Eigen::Vector3d G, std::vector<double> v){
    //define variables for symmetric and asymmetric part of the potential
	double vs = 0.0e0, va = 0.0e0;
	std::complex<double> potential = 0.0e0;

	//Each cell has 2 atoms.  Use mid point between them as offset
	Eigen::Vector3d vec_tau;
    vec_tau << 0.125, 0.125, 0.125;

	//Find the magnitude of the reciprocal vector
    double GMAG = G.squaredNorm();

	//assign symmetric and asymmetric part of the potential 
	if (abs(double(GMAG - 3.0)) < 1e-6){ vs = v[0]; va = v[3]; }
	else if (abs(double(GMAG - 4.0)) < 1e-6){ vs = 0.0; va = v[4]; }
	else if (abs(double(GMAG - 8.0)) < 1e-6){ vs = v[1]; va = 0.0; }
	else if (abs(double(GMAG - 11.0)) < 1e-6){ vs = v[2]; va = v[5]; }
	else { vs = 0.0e0; va = 0.0e0; }

	//Complex Pseudopotential value, potentials are scaled to rydberg units
    if (ORIGINSHIFT){
        potential = std::complex<double>(
                        vs*cos(2.0*pi*(G.dot(vec_tau)))*Ry2eV, 
                        va*sin(2.0*pi*(G.dot(vec_tau)))*Ry2eV );
    }else{
        double V1, V2;
        V1 = (vs+va)/2;
        V2 = (vs-va)/2;
        potential = std::complex<double>(
                        (V1 + V2*cos(-4.0*pi*(G.dot(vec_tau))))*Ry2eV,
                        (V2*sin(-4.0*pi*(G.dot(vec_tau))))*Ry2eV );
    }
	return potential;
}

// complex nonlocal-pseudopotential for the given reciprocal vector
//a : lattice constant
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


// spherical Bessel function of the first kind
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

//set reference energyoffset 
void setRefEnergy(env &dat){

    // the global Hamiltonian matrix H at gamma point
    Eigen::Vector3d K_GAMMA = {0,0,0};

    Eigen::MatrixXcd H = HamiltonianEPM(dat.G,K_GAMMA,dat.epm);

    // compute the eigenvalue problem
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver(H);
    if (eigsolver.info() != Eigen::Success){ std::cout << "ERROR::setRefEnergy:: Eigensolver failed to solve the matrix at Gamma point." << std::endl; abort(); }

    // zero level of energy at the valence band maximum at gamma point
    dat.energyoffset = eigsolver.eigenvalues()(NVALENCE-1);
    // dat.eigerror = abs(eigsolver.eigenvalues()(2)-eigsolver.eigenvalues()(1))/2.0;

    return;
}

}