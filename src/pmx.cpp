#include "pmx.h"

//======================================================//
//                      main                            //
//======================================================//
int main (int argc, char **argv) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;

    // input parser
    pmx::inputParser inp(argc, argv);

    // help
    if(inp.existOption("-h") || inp.existOption("-help")){
        pmx::printHelp();
        return 0;
    }

    // version
    if(inp.existOption("-v") || inp.existOption("-version")){
        pmx::printVersion();
        return 0;
    }

    // parameters
    pmx::env dat;
    pmx::inputProcess(inp, dat);
    // check parameters before proceeding calculations


    //------------------------------------------
    // construct k-grids (not required for bandstructure calculation)
    if (dat.gridscheme==0){
        pmx::genkgrid_cohen_chadi(dat);
    }else if (dat.gridscheme==1){
        pmx::genkgrid_MP(dat);
    }
    
    pmx::generate_planewaves(dat);

    //------------------------------------------
    // print out input parameter (for debugging purposes)
    if(dat.debug>1){
        pmx::printInput(dat);
        return 0;
    }

    //------------------------------------------
    // electronic bandstructure calculation
    if(inp.existOption("-switch") 
        && (inp.valueOption("-switch")=="band")){

        start = std::chrono::system_clock::now();
        pmx::empiricalpseudopotentialmethod(dat);
        end = std::chrono::system_clock::now();
        elapsed_seconds = end-start;

        // print bandstructure results
        pmx::printBand(dat);
    }

    //------------------------------------------
    // permittivity calculation
    if(inp.existOption("-switch") 
        && (inp.valueOption("-switch")=="eps")){
        start = std::chrono::system_clock::now();

        // obtain reference energy level
        pmx::setRefEnergy(dat);

        // permittivity calculation
        pmx::permittivity(dat);

        // pmx::printEpsilon(dat);
        pmx::writeEpsilon(dat);

        end = std::chrono::system_clock::now();
        elapsed_seconds = end-start;
    }
    


    // std::cout << "permittivity calculation...\n";
    // start = std::chrono::system_clock::now();
    // for (int q=0; q<dat.nq; q++){
    //     pmx::susceptibility(dat, dat.q[q]);
    // }
    // std::cout << "completed!\n";
    // end = std::chrono::system_clock::now();
    // elapsed_seconds = end-start;
    // std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";

    return 0;
}

namespace pmx{
//======================================================//
//                      preprocessing                   //
//======================================================//

// generate plane waves G
void generate_planewaves (env &dat){
    std::vector<int> g2 {0,3,4,8,11,12,16,19,20,24,27,32,35,36,40,43,44,48,51,52,56,59,64,67,68,72,75,76,80,83,84,88,91,96,99,100};
    std::vector<int> nnn {1,8,6,12,24,8,6,24,24,24,32,12,48,30,24,24,24,8,48,24,48,72,6,24,48,36,56,24,24,72,48,24,48,24,72,30};
    std::vector<int>::iterator itr = std::lower_bound(g2.begin(), g2.end(), pow(dat.GCUT,2));
    int N = int(itr - g2.begin());
    std::vector<int> cnn(N);

    NPW = 0;
    for (int i=0; i<N; i++){
        NPW += nnn[i];
        if (i==0) cnn[i] = 1;
        else cnn[i] = cnn[i-1]+nnn[i];
    }

    dat.G.reserve(NPW);
    for (int i=0; i<NPW; i++){
        dat.G.push_back({0,0,0});
    }

    #pragma omp parallel for
    for (int i=0; i<N; i++){
        std::vector<Eigen::Vector3d> iG = find_G(g2[i],nnn[i]);
        int in;
        if (i==0) in = 0;
        else in = cnn[i-1];
        for (int ic=0; ic<nnn[i]; ic++){
            dat.G[in+ic] = iG[ic];
        }
    }
    #pragma omp barrier

    find_locGi(dat);

    return;
}

// find n reciprocal lattice vectors with squared magnitude g2
std::vector<Eigen::Vector3d> find_G (int g2, int n){
    std::vector<Eigen::Vector3d> G(n);

    int ia = int(sqrt(double(g2)));
    int gx2, gy2, gz2, gxy2, gxyz2;
    int ic = 0;

    for (int gx=-ia; gx<=ia; gx++){
        gx2 = gx*gx;
        for (int gy=-ia; gy<=ia; gy++){
            gy2 = gy*gy;
            gxy2 = gx2+gy2;
            for (int gz=-ia; gz<=ia; gz++){
                gz2 = gz*gz;
                gxyz2 = gxy2+gz2;
                if (gxyz2==g2){
                    G[ic] << double(gx), double(gy), double(gz);
                    ic +=1;
                    if (ic==n) return G;
                }
            }
        }
    }
    return G;
}

// find local-field index
void find_locGi (env &dat){

    // initialize
    dat.loci = new int *[NEPS];
    for (int m=0; m<NEPS; m++){
        dat.loci[m] = new int [NPW];
    }

    // g(0)
    for (int i=0; i<NPW; i++){
        dat.loci[0][i] = i;
    }
    #pragma omp parallel for
    for (int m=1; m<NEPS; m++){

        // local-field vector Gm
        Eigen::Vector3d Gm = dat.G[m];
        for (int i=0; i<NPW; i++){
            // bra vector Gi
            Eigen::Vector3d Gi = dat.G[i];

            // if no such vector found, -1
            dat.loci[m][i] = -1;
            for (int j=0; j<NPW; j++){
                // ket vector Gj
                Eigen::Vector3d Gj = dat.G[j];

                // find j, such that Gi = Gj-Gm
                if ((Gi-Gj+Gm).norm()<1e-10){
                    dat.loci[m][i] = j;
                    break;
                }
            }
        }
    }
    #pragma omp barrier
    return;
}




void genkgrid_cohen_chadi (env &dat) {
    //Here we will use Cohen and Chadi scheme: Phys. Rev. B 8, 5757 (1973).
    int i, j, k, l; //dummy indices
    int N = dat.gridorder;
    NKPT = int(pow(2.0,N-1)*(pow(2.0,N)+1.0)*(pow(2.0,N-1.0)+1.0)/3.0);
    dat.kgrid = new double*[NKPT];// xi [N x 3]
    dat.kweight = new double[NKPT];// wi [N]
    for (i=0; i<NKPT; i++){
        dat.kgrid[i] = new double[3];
    }
    
    //2-point Cohen and Chadi scheme: Phys. Rev. B 8, 5757 (1973).
    if (N==1){
        // 0 point
        dat.kweight[0] = 0.7500; dat.kgrid[0][0] = 0.75; dat.kgrid[0][1] = 0.25; dat.kgrid[0][2] = 0.25;
        // 1 point
        dat.kweight[1] = 0.2500; dat.kgrid[1][0] = 0.25; dat.kgrid[1][1] = 0.25; dat.kgrid[1][2] = 0.25;
    }

    //10-point Cohen and Chadi scheme: Phys. Rev. B 8, 5757 (1973).
    else if (N==2){
        // 0 point
        dat.kweight[0] = 0.18750; dat.kgrid[0][0] = 0.875; dat.kgrid[0][1] = 0.375; dat.kgrid[0][2] = 0.125;
        // 1 point
        dat.kweight[1] = 0.09375; dat.kgrid[1][0] = 0.875; dat.kgrid[1][1] = 0.125; dat.kgrid[1][2] = 0.125;
        // 2 point
        dat.kweight[2] = 0.09375; dat.kgrid[2][0] = 0.625; dat.kgrid[2][1] = 0.625; dat.kgrid[2][2] = 0.125;
        // 3 point
        dat.kweight[3] = 0.09375; dat.kgrid[3][0] = 0.625; dat.kgrid[3][1] = 0.375; dat.kgrid[3][2] = 0.375;
        // 4 point
        dat.kweight[4] = 0.18750; dat.kgrid[4][0] = 0.625; dat.kgrid[4][1] = 0.375; dat.kgrid[4][2] = 0.125;
        // 5 point
        dat.kweight[5] = 0.09375; dat.kgrid[5][0] = 0.625; dat.kgrid[5][1] = 0.125; dat.kgrid[5][2] = 0.125;
        // 6 point
        dat.kweight[6] = 0.03125; dat.kgrid[6][0] = 0.375; dat.kgrid[6][1] = 0.375; dat.kgrid[6][2] = 0.375;
        // 7 point
        dat.kweight[7] = 0.09375; dat.kgrid[7][0] = 0.375; dat.kgrid[7][1] = 0.375; dat.kgrid[7][2] = 0.125;
        // 8 point
        dat.kweight[8] = 0.09375; dat.kgrid[8][0] = 0.375; dat.kgrid[8][1] = 0.125; dat.kgrid[8][2] = 0.125;
        // 9 point
        dat.kweight[9] = 0.03125; dat.kgrid[9][0] = 0.125; dat.kgrid[9][1] = 0.125; dat.kgrid[9][2] = 0.125;
    }

    //Higher-order Cohen and Chadi scheme: PHILOSOPHICAL MAGAZINE B 81, NO. 6, 551--559 (2001).
    else if (N>2){
        k = 0;
        for (i=1; i<=pow(2,N-1); i++)	for (j =1; j<=i; j++)	for (l=1; l<=j; l++){
            //std::cout << "Loop 1\n" << std::endl;
            dat.kweight[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.kgrid[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.kgrid[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.kgrid[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }	
        for (i=pow(2,N-1)+2; i<=pow(2,N); i+=2)	for (j=1; j<=3*pow(2,N-2)-int(i*0.5); j++)	for (l=1; l<=j; l++){
            //std::cout << "Loop 2\n" << std::endl;
            dat.kweight[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.kgrid[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.kgrid[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.kgrid[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }	
        for (i=pow(2,N-1)+2; i<=3*pow(2,N-2); i+=2)	for (j=3*pow(2,N-2)-int(i*0.5)+1; j<=i; j++)	for (l=1; l<=3*pow(2,N-1)+1-i-j; l++){
            //std::cout << "Loop 3\n" << std::endl;
            dat.kweight[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.kgrid[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.kgrid[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.kgrid[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }
        for (i=3*pow(2,N-2)+2; i<=pow(2,N); i+=2)	for (j=3*pow(2,N-2)-int(i*0.5)+1; j<=3*pow(2,N-1)-i; j++)	for (l=1; l<=3*pow(2,N-1)+1-i-j; l++){
            //std::cout << "Loop 4\n" << std::endl;
            dat.kweight[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.kgrid[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.kgrid[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.kgrid[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }
        for (i=pow(2,N-1)+1; i<=pow(2,N)-1; i+=2)	for (j=1; j<=3*pow(2,N-2)-int((i-1)*0.5); j++)	for (l=1; l<=j; l++){
            //std::cout << "Loop 5\n" << std::endl;
            dat.kweight[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.kgrid[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.kgrid[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.kgrid[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }
        for (i=pow(2,N-1)+1; i<=3*pow(2,N-2)-1; i+=2)	for (j=3*pow(2,N-2)-int((i-1)*0.5)+1; j<=i; j++)	for (l=1; l<=3*pow(2,N-1)+1-i-j; l++){
            //std::cout << "Loop 6\n" << std::endl;
            dat.kweight[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.kgrid[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.kgrid[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.kgrid[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }
        for (i=3*pow(2,N-2)+1; i<=pow(2,N)-1; i+=2)	for (j=3*pow(2,N-2)-int((i-1)*0.5)+1; j<=3*pow(2,N-1)-i; j++)	for (l=1; l<=3*pow(2,N-1)+1-i-j; l++){
            //std::cout << "Loop 7\n" << std::endl;
            dat.kweight[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.kgrid[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.kgrid[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.kgrid[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }
    }
    return;
}

void genkgrid_MP (env &dat) {
    std::ifstream kdata(dat.kgridfile);
    int i = 0; //dummy index
    int NKPT0;

    std::string line;
    if (kdata.is_open())
    {
        while(getline (kdata,line)){
        if (!(line[0]=='#' || line=="")){
            //Line 0: Number of integration k-points in entire BZ
            if (i == 0){
                std::stringstream(line) >> NKPT0;
                i++;
                continue;	
            }
            //Line 1: Number of irreducible k-points in BZ
            else if (i==1){
                std::stringstream(line) >> NKPT;  
                dat.kgrid = new double*[NKPT];
                dat.kweight = new double[NKPT];
                for (int j=0; j<NKPT; j++){
                    dat.kgrid[j] = new double[3];
                }//end of for j
                i++;
                continue;
            }
            else if (i>1 && i-2<NKPT){
                //Line 1: Number of plane waves
                std::stringstream(line) >> dat.kgrid[i-2][0] >>  dat.kgrid[i-2][1] >>  dat.kgrid[i-2][2] >>  dat.kweight[i-2];  
                dat.kweight[i-2] = dat.kweight[i-2]/NKPT0;
                i++;
                continue; 
            
            }//else if i
            }//end of if line
        }//end of while 
    }//end of if kdata open
    kdata.close();
    return;
}

// Krnoecker Delta Function
inline bool kroneckerDelta(int m, int n) {
	return (m == n);
}


//======================================================//
//                      bandstructure                   //
//======================================================//

// hamiltonian
Eigen::MatrixXcd HamiltonianEPM (std::vector<Eigen::Vector3d> G, Eigen::Vector3d K, std::vector<double> v){
    Eigen::MatrixXcd H(NPW, NPW);
    double KINETIC_CONST = KCON*(4.0*pi*pi)/(a*a);
    bool NONLOCAL = false;
    if (v.size()==11) NONLOCAL = true;

    for (int i=0; i<NPW; i++){ for (int j=0; j<=i; j++){
        if (i==j){ // diagonal kinetic part
            // H(i,j) = KINETIC_CONST * ((K+G[i]).squaredNorm());
            H(i,j) = KINETIC_CONST * pow((K+G[i]).norm(),2);
        }else{ // local pseudopotential
            Eigen::Vector3d dG = G[i]-G[j]; // TEST Gj-Gi or Gi-Gj
            if (dG.norm()<3.32) // TEST dG.squaredNorm()<11.1
                H(i,j) = pseudopotential(dG,v);
            else
                H(i,j) = 0;
        }
        if (NONLOCAL){ // nonlocal pseudopotential
            H(i,j) += nonlocalpseudopotential(K+G[i],K+G[j],v);
        }
    }}
    // only lower triangle is needed for eigenvalue calculation
    return H;
}








// empirical pseudopotential method
void empiricalpseudopotentialmethod(env &dat){
    // int N = dat.numpw;
    dat.band = new double*[dat.nq];

    setRefEnergy(dat);

    #pragma omp parallel for
    for (int q=0; q<dat.nq; q++){
        
        // Hamiltonian
        Eigen::MatrixXcd H = HamiltonianEPM(dat.G,dat.q[q],dat.epm);

        // diagonalize the Hamiltonian
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver(H);
        if (eigsolver.info() != Eigen::Success){ std::cout << "ERROR: Eigensolver failed to solve the matrix." << std::endl; abort(); }

        //
        dat.band[q] = new double[NBAND];
        for (int n=0; n<NBAND; n++){
            dat.band[q][n] = double(eigsolver.eigenvalues()(n))-dat.energyoffset;
        }
    }
    #pragma omp barrier

    return;
}

// // 
// inline std::vector<double> millerindex(int n, int N){
//     int m = int(cbrt(N));
//     std::vector<double> i {
//         floor(double(n + ((N-1)/2.0))/(m*m)) - floor(m/2.0),
// 	    floor(((n + (N-1)/2)%(m*m))/m) - floor(m/2.0),
// 	    ((n + (N-1)/2)%m) - floor(m/2.0)};
//     return i;
// }

// // Generate reciprocal lattice vector
// inline Eigen::Vector3d recvec(int n, int N, std::vector<Eigen::Vector3d> b){
//     std::vector<double> i = millerindex(n,N);
// 	//Calculate the reciprocal lattice vector
// 	Eigen::Vector3d Kvec = i[0]*b[0] + i[1]*b[1] + i[2]*b[2];
// 	return Kvec;
// }

// // Generate matrix element m, n
// std::complex<double> matrixelement(int m, int n, int N, Eigen::Vector3d kvec, std::vector<Eigen::Vector3d> b, std::vector<double> v, double a){
// 	std::complex<double> Hmn;

//     if (v.size()==6){
//         if (m!=n)
//             Hmn =  pseudopotential(recvec(m-n,N,b),v);
//         else
//             Hmn = 0.0;
//     }else {
//         Hmn =  pseudopotential(recvec(m-n,N,b),v) 
//             + nonlocalpseudopotential(kvec+recvec(m,N,b),kvec+recvec(n,N,b),v,a);
//     }
// 	// if( m==n ) Hmn += cscale*(4.0*pi*pi)/(a*a) * pow((kvec+recvec(m,N,b)).norm(),2);

// 	if( m==n ) Hmn += cscale*(4.0*pi*pi)/(a*a)*abs((kvec+recvec(m,N,b)).dot(kvec+recvec(m,N,b)));
// 	return Hmn;
// }

// complex pseudopotential for the given reciprocal vector
std::complex<double> pseudopotential(Eigen::Vector3d G, std::vector<double> v){
    //define variables for symmetric and asymmetric part of the potential
	double vs = 0.0e0, va = 0.0e0;
	std::complex<double> potential = 0.0e0;

	//Each cell has 2 atoms.  Use mid point between them as offset
	Eigen::Vector3d vec_tau;
    vec_tau << 0.125, 0.125, 0.125;

	//Find the magnitude of the reciprocal vector
	// double GMAG = pow(G.norm(),2);
    double GMAG = G.squaredNorm();

	//assign symmetric and asymmetric part of the potential 
	if (abs(double(GMAG - 3.0)) < 1e-6){ vs = v[0]; va = v[3]; }
	else if (abs(double(GMAG - 4.0)) < 1e-6){ vs = 0.0; va = v[4]; }
	else if (abs(double(GMAG - 8.0)) < 1e-6){ vs = v[1]; va = 0.0; }
	else if (abs(double(GMAG - 11.0)) < 1e-6){ vs = v[2]; va = v[5]; }
	else { vs = 0.0e0; va = 0.0e0; }

	//Complex Pseudopotential value, potentials are scaled to rydberg units
	// potential = std::complex<double>(
    //                 vs*cos(2.0*pi*(G.dot(vec_tau)))*Ry2eV, 
    //                 va*sin(2.0*pi*(G.dot(vec_tau)))*Ry2eV );
    // std::complex<double> I = (0.0,1.0);
    // potential *= exp(-I*2.0*pi*(G.dot(vec_tau)));
    // potential *= std::complex<double> (cos(2.0*pi*(G.dot(vec_tau))),-sin(2.0*pi*(G.dot(vec_tau))));

    double V1, V2;
    V1 = (vs+va)/2;
    V2 = (vs-va)/2;
    potential = std::complex<double>(
                    (V1 + V2*cos(-4.0*pi*(G.dot(vec_tau))))*Ry2eV,
                    (V2*sin(-4.0*pi*(G.dot(vec_tau))))*Ry2eV );


    // potential = std::complex<double>(
    //                 ((V1+V2)*cos(2.0*pi*(G.dot(vec_tau))))*Ry2eV,
    //                 ((V1-V2)*sin(2.0*pi*(G.dot(vec_tau))))*Ry2eV );
    // potential *= Ry2eV;


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

//======================================================//
//                 density of states                    //
//======================================================//

// Calculate density of states
// void dos(env &dat){
//     double epsilon = 0.001;

//     std::ofstream dosfile;
// 	dosfile.open("dos.out");
// 	dosfile.setf(std::ios::left); //setf sets format flags
// 	dosfile.setf(std::ios::scientific);
// 	dosfile.setf(std::ios::showpoint);
// 	dosfile.setf(std::ios::showpos);
// 	dosfile.precision(8);
// 	double energy = -5.0;
// 	double dosvalue;
//     int N = dat.numpw;
// 	Eigen::MatrixXcd A(N, N); //Global Matrix

//     std::complex<double> eigenvec[dat.numk][30][N];
// 	Eigen::MatrixXd eigenvalues(dat.numk, 30);
// 	int i, j, k, n; //dummy indices

//     setRefEnergy(dat);

//     Eigen::Vector3d kvector;
//     for (int k=0; k<dat.numk; k++){
//         kvector << dat.kgrid[k][0], dat.kgrid[k][1], dat.kgrid[k][2];

//         for(i=-N/2; i<N/2; i++){for(j=-N/2; j<N/2; j++){
// 			A(i+N/2, j+N/2) = matrixelement(i,j,N,kvector,dat.b,dat.epm,dat.a); 		
// 		}}

//         Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(A);
// 		if (eigensolver.info() !=Eigen::Success){ std::cout << "ERROR: Eigensolver failed to solve the matrix." << std::endl; abort(); }
// 		for(n = 0; n < 30; n++){
// 			eigenvalues(k, n) = eigensolver.eigenvalues()(n) - dat.energyoffset;	

// 			for(j=-N/2; j<N/2; j++)	eigenvec[k][n][j+N/2] = eigensolver.eigenvectors().col(n)(j+N/2);
// 		}//end of n
//     }
// }



//======================================================//
//                      permittivity                    //
//======================================================//

// Calculate lattice longitudinal and transverse permittivity matrix elements at G1, G2
void permittivity(env &dat){

    // scale factors
    double scaleL = 8.0*pow(e,2)/(pi*eV*a*angstrom);
    double scaleT = 128.0*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));

    std::complex<double> ***eigvector_k; // eigenvector at k
    double **eigvalue_k; // eigenvalue at k
    std::complex<double> ***eigvector_kq; // eigenvector at k+q
    double **eigvalue_kq; // eigenvalue at k

    // overlap integral L, <k,v|e^{-i*G_m}|k+q,c> <k+q,c|e^{i*G_n}|k,v>
    std::complex<double> *****ointL;
    // overlap integral T, t_{Q+G_m} * <k,v|e^{-i*G_m}j_0|k+q,c> <k+q,c|e^{i*G_n}j_0|k,v> * t_{Q+G_n}
    std::complex<double> *****ointT;

    int NBAND_V [NKPT]; // number of valence bands
    int NBAND_C [NKPT]; // number of conduction bands

    dat.realepsL = new double ***[dat.nq];
    dat.imagepsL = new double ***[dat.nq];
    dat.realepsT = new double ***[dat.nq];
    dat.imagepsT = new double ***[dat.nq];



    //======================================================================
    eigvalue_k = new double *[NKPT];
    eigvector_k = new std::complex<double> **[NKPT];
    #pragma omp parallel for
    for (int k=0; k<NKPT; k++){
        // vector{k}
        Eigen::Vector3d K = {dat.kgrid[k][0], dat.kgrid[k][1], dat.kgrid[k][2]};

        // Hamiltonian at k
        Eigen::MatrixXcd H = HamiltonianEPM(dat.G,K,dat.epm);   
        
        // compute eigenvalue problem at k
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);
        if (eigsolver_k.info() != Eigen::Success){ printf("ERROR: Eigensolver failed to solve the matrix."); abort(); }

        eigvalue_k[k] = new double [NBAND];
        for (int n=0; n<NBAND; n++){
            eigvalue_k[k][n] = eigsolver_k.eigenvalues()(n)-dat.energyoffset;
            if (eigvalue_k[k][n] > 0) {
                NBAND_V[k] = n;
                break;
            }
        }
        eigvector_k[k] = new std::complex<double> *[NBAND_V[k]];
        for (int n=0; n<NBAND_V[k]; n++){
            eigvector_k[k][n] = new std::complex<double> [NPW];
            for (int j=0; j<NPW; j++){
                eigvector_k[k][n][j] = eigsolver_k.eigenvectors().col(n)(j);
            }
        }
    }
    #pragma omp barrier
    

    //======================================================================
    ointL = new std::complex<double> ****[NKPT];
    ointT = new std::complex<double> ****[NKPT];
    eigvalue_kq = new double *[NKPT];
    eigvector_kq = new std::complex<double> **[NKPT];
    for (int k=0; k<NKPT; k++){
        ointL[k] = new std::complex<double>***[NEPS];
        ointT[k] = new std::complex<double>***[NEPS];
        for (int m=0; m<NEPS; m++){
            ointL[k][m] = new std::complex<double> **[NEPS];
            ointT[k][m] = new std::complex<double> **[NEPS];
            for (int n=0; n<=m; n++){
                ointL[k][m][n] = new std::complex<double> *[NBAND_V[k]];
                ointT[k][m][n] = new std::complex<double> *[NBAND_V[k]];
            }
        }
        eigvalue_kq[k] = new double [NBAND];
    }

    for (int q=0; q<dat.nq; q++){
        Eigen::Vector3d Q = dat.q[q];
        #pragma omp parallel for
        for (int k=0; k<NKPT; k++){
            
            // vector{k}
            Eigen::Vector3d K = {dat.kgrid[k][0], dat.kgrid[k][1], dat.kgrid[k][2]};

            // Hamiltonian at k+q
            Eigen::MatrixXcd H = HamiltonianEPM(dat.G,K+Q,dat.epm);   

            // compute eigenvalue problem at k+q
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_kq(H);
            if(eigsolver_kq.info() != Eigen::Success){ printf("ERROR: Eigensolver failed to solve the matrix."); abort(); }
            
            
            int tempindex;
            for (int n=NBAND-1; n>=0; n--){
                tempindex = NBAND-1-n;
                eigvalue_kq[k][tempindex] = double(eigsolver_kq.eigenvalues()(n))-dat.energyoffset;
                if (eigvalue_kq[k][tempindex] <=0) {
                    NBAND_C[k] = tempindex;
                    break;
                }
            }
            eigvector_kq[k] = new std::complex<double> *[NBAND_C[k]];
            for (int n=0; n<NBAND_C[k]; n++){
                tempindex = NBAND-1-n;
                eigvector_kq[k][n] = new std::complex<double> [NPW];
                for (int j=0; j<NPW; j++){
                    eigvector_kq[k][n][j] = eigsolver_kq.eigenvectors().col(tempindex)(j);
                }
            }

            // overlap integrals
            for (int m=0; m<NEPS; m++){
                for (int n=0; n<=m; n++){
                    // Eigen::Vector3d vecqGm = Q+dat.G[m]; // vector q+Gm     
                    // Eigen::Vector3d vecqGn = Q+dat.G[n]; // vector q+Gn    
                    // Eigen::Vector3d vecqGm_t = perpvec(vecqGm); // normal vector perpendicular to q+Gm
                    // Eigen::Vector3d vecqGn_t = perpvec(vecqGn); // normal vector perpendicular to q+Gn


                    // vecqGm_t = vecqGm_t.cross(vecqGm);
                    // vecqGn_t = vecqGn_t.cross(vecqGn);

                    // Eigen::Vector3d vecqGm = Q; // vector q+Gm     
                    // Eigen::Vector3d vecqGn = Q; // vector q+Gn    
                    // Eigen::Vector3d vecqGm_t = perpvec(vecqGm); // normal vector perpendicular to q+Gm
                    // Eigen::Vector3d vecqGn_t = perpvec(vecqGm); // normal vector perpendicular to q+Gn
                    // vecqGm_t = vecqGm_t.cross(vecqGm);
                    // vecqGn_t = vecqGn_t.cross(vecqGn);
                    
                    //====================================================================
                    // // LL
                    // Eigen::Vector3d vecqGm_t = Q+dat.G[m];
                    // Eigen::Vector3d vecqGn_t = Q+dat.G[n];
                    // vecqGm_t = vecqGm_t/vecqGm_t.norm();
                    // vecqGn_t = vecqGn_t/vecqGn_t.norm();

                    // TT
                    Eigen::Vector3d vecqGm_t = perpvec(Q+dat.G[m]);
                    Eigen::Vector3d vecqGn_t = perpvec(Q+dat.G[n]);
                    // if (m==3 || m==4 || m==7 || m==8){
                    //     Eigen::Vector3d normQG = Q+dat.G[m];
                    //     normQG = normQG/normQG.norm();
                    //     vecqGm_t = vecqGm_t.cross(normQG);
                    // }
                    // if (n==3 || n==4 || n==7 || n==8){
                    //     Eigen::Vector3d normQG = Q+dat.G[n];
                    //     normQG = normQG/normQG.norm();
                    //     vecqGn_t = vecqGn_t.cross(normQG);
                    // }


                    // // T'T'
                    // Eigen::Vector3d vecqGm_t = perpvec(Q+dat.G[m]);
                    // Eigen::Vector3d vecqGn_t = perpvec(Q+dat.G[n]);
                    // Eigen::Vector3d normQG = Q+dat.G[m];
                    // normQG = normQG/normQG.norm();
                    // vecqGm_t = vecqGm_t.cross(normQG);
                    // normQG = Q+dat.G[n];
                    // normQG = normQG/normQG.norm();
                    // vecqGn_t = vecqGn_t.cross(normQG);

                    // // TT'
                    // Eigen::Vector3d vecqGm_t = perpvec(Q+dat.G[m]);
                    // Eigen::Vector3d vecqGn_t = perpvec(Q+dat.G[n]);
                    // Eigen::Vector3d normQG = Q+dat.G[n];
                    // normQG = normQG/normQG.norm();
                    // vecqGn_t = vecqGn_t.cross(normQG);


                    // // LT
                    // Eigen::Vector3d vecqGm_t2 = Q+dat.G[m];
                    // vecqGm_t2 = vecqGm_t2/vecqGm_t.norm();
                    // Eigen::Vector3d vecqGn_t2 = perpvec(Q+dat.G[n]);







                    // non-G ==================================================

                    // LL
                    // Eigen::Vector3d vecqGm_t = Q/Q.norm();
                    // Eigen::Vector3d vecqGn_t = Q/Q.norm();

                    // LT
                    // Eigen::Vector3d vecqGm_t = Q/Q.norm();
                    // Eigen::Vector3d vecqGn_t = perpvec(Q);

                    // // TT
                    // Eigen::Vector3d vecqGm_t = perpvec(Q);
                    // Eigen::Vector3d vecqGn_t = perpvec(Q);

                    // TT'
                    // Eigen::Vector3d vecqGm_t = perpvec(Q);
                    // Eigen::Vector3d vecqGn_t = perpvec(Q);
                    // vecqGn_t = vecqGn_t.cross(Q)/Q.norm();


                    for (int v=0; v<NBAND_V[k]; v++){ // valence bands
                        ointL[k][m][n][v] = new std::complex<double> [NBAND_C[k]];
                        ointT[k][m][n][v] = new std::complex<double> [NBAND_C[k]];
                        for (int c=0; c<NBAND_C[k]; c++){ // conduction bands

                            std::complex<double> ointL_k_kq = 0.0; // <k,v|e^{-i*G_m}|k+q,c>
                            std::complex<double> ointT_k_kq = 0.0; // <k,v|e^{-i*G_m}j_0|k+q,c> * t_{G_m+Q}
                            std::complex<double> ointL_kq_k = 0.0; // <k+q,c|e^{i*G_n}|k,v>
                            std::complex<double> ointT_kq_k = 0.0; // <k+q,c|e^{i*G_n}j_0|k,v> * t_{G_n+Q}

                            for (int i=0; i<NPW; i++){
                                if (dat.loci[m][i]!=-1){
                                    ointL_k_kq += conj(eigvector_k[k][v][i]) * eigvector_kq[k][c][dat.loci[m][i]];
                                    // ointL_k_kq += conj(eigvector_k[k][v][i]) * eigvector_kq[k][c][dat.loci[m][i]] 
                                    //                             * (dat.G[i]+K +dat.G[m]/2+Q/2).dot(vecqGm_t2);
                                    ointT_k_kq += conj(eigvector_k[k][v][i]) * eigvector_kq[k][c][dat.loci[m][i]] 
                                                                * (dat.G[i]+K +dat.G[m]/2+Q/2).dot(vecqGm_t); // 
                                }
                                if (dat.loci[n][i]!=-1){
                                    ointL_kq_k += conj(eigvector_kq[k][c][dat.loci[n][i]]) * eigvector_k[k][v][i];
                                    // ointL_kq_k += conj(eigvector_kq[k][c][dat.loci[n][i]]) * eigvector_k[k][v][i] 
                                    //                             * (dat.G[i]+K +dat.G[n]/2+Q/2).dot(vecqGn_t2);
                                    ointT_kq_k += conj(eigvector_kq[k][c][dat.loci[n][i]]) * eigvector_k[k][v][i] 
                                                                * (dat.G[i]+K +dat.G[n]/2+Q/2).dot(vecqGn_t); //
                                }
                            }

                            ointL[k][m][n][v][c] = ointL_k_kq * ointL_kq_k;
                            ointT[k][m][n][v][c] = ointT_k_kq * ointT_kq_k;
                        }
                    }
                }
            }
        }
        #pragma omp barrier


        dat.realepsL[q] = new double **[NEPS];
        dat.imagepsL[q] = new double **[NEPS];
        dat.realepsT[q] = new double **[NEPS];
        dat.imagepsT[q] = new double **[NEPS];
        
        for (int m=0; m<NEPS; m++){
            dat.realepsL[q][m] = new double *[NEPS];
            dat.imagepsL[q][m] = new double *[NEPS];
            dat.realepsT[q][m] = new double *[NEPS];
            dat.imagepsT[q][m] = new double *[NEPS];
            for (int n=0; n<=m; n++){
                dat.realepsL[q][m][n] = new double [NFREQ];
                dat.imagepsL[q][m][n] = new double [NFREQ];
                dat.realepsT[q][m][n] = new double [NFREQ];
                dat.imagepsT[q][m][n] = new double [NFREQ];
                #pragma omp parallel for //collapse(3)
                for (int f=0; f<NFREQ; f++){
                    // temporary variables
                    std::complex<double> realL1 = 0.0e0;
                    std::complex<double> imagL1 = 0.0e0;
                    std::complex<double> realT1 = 0.0e0;
                    std::complex<double> imagT1 = 0.0e0;
                    std::complex<double> depsilon = 0;
                    double realdL, imagdL, realdT, imagdT;

                    // sum over k-grids & bands
                    for (int k=0; k<NKPT; k++){
                        for (int v=0; v<NBAND_V[k]; v++){ for (int c=0; c<NBAND_C[k]; c++){
                            depsilon = eigvalue_kq[k][c]-eigvalue_k[k][v];
                            imagL1 += dat.kweight[k] * (ointL[k][m][n][v][c]) 
                                                    * diracdelta(depsilon-dat.freq[f], dat.epsilon);
                            imagT1 += dat.kweight[k] * (ointT[k][m][n][v][c]) 
                                                    * diracdelta(depsilon-dat.freq[f], dat.epsilon);
                            if (dat.kk==0){
                                realL1 += dat.kweight[k] * (ointL[k][m][n][v][c]) 
                                                        * (depsilon) / (depsilon*depsilon-dat.freq[f]*dat.freq[f]);
                                // realL1 += dat.kweight[k] * (ointL[k][m][n][v][c]) 
                                //                         / (depsilon) / (depsilon*depsilon-dat.freq[f]*dat.freq[f]);
                                realT1 += dat.kweight[k] * (ointT[k][m][n][v][c]) 
                                                        / (depsilon) / (depsilon*depsilon-dat.freq[f]*dat.freq[f]);
                            }
                        }}
                    }

                    // longitudinal permittivity
                    Eigen::Vector3d vecqGm = Q+dat.G[m];
                    Eigen::Vector3d vecqGn = Q+dat.G[n];
                    imagdL = 0.0e0;
                    if (dat.freq[f]>0) imagdL = scaleL * pi * imagL1.real() / (vecqGm.norm()*vecqGn.norm());
                    dat.imagepsL[q][m][n][f] = imagdL;
                    
                    // transverse permittivity
                    imagdT = 0.0e0;
                    if (dat.freq[f]>0) imagdT = scaleT * pi * imagT1.real() / (dat.freq[f]*dat.freq[f]);
                    dat.imagepsT[q][m][n][f] = imagdT;

                    if (dat.kk==0){
                        realdL = scaleL * 2 * realL1.real() / (vecqGm.norm()*vecqGn.norm());
                        // realdL = scaleT * 2 * realL1.real();
                        if (m==n) realdL += 1.0;
                        dat.realepsL[q][m][n][f] = realdL;
                        realdT = scaleT * 2 * realT1.real();
                        if (m==n) realdT += 1.0;
                        dat.realepsT[q][m][n][f] = realdT;
                    }
                }
            }
        }
        #pragma omp barrier

        // kramer-kronig transform
        if (dat.kk){
            for (int m=0; m<NEPS; m++){
                for (int n=0; n<=m; n++){
                    kramerskronigtransform(dat.realepsL[q][m][n],dat.imagepsL[q][m][n],dat.freq,dat.dfreq);
                    kramerskronigtransform(dat.realepsT[q][m][n],dat.imagepsT[q][m][n],dat.freq,dat.dfreq);
                    if (m==n){
                        for (int f=0; f<NFREQ; f++){
                            dat.realepsL[q][m][n][f] += 1;
                            dat.realepsT[q][m][n][f] += 1;
                        }
                    }
                }
            }
        }

        // free dynamically allocated memory
        for (int k=0; k<NKPT; k++){
            for (int m=0; m<NEPS; m++){
                for (int n=0; n<=m; n++){
                    for (int v=0; v<NBAND_V[k]; v++){
                        delete [] ointL[k][m][n][v];
                        delete [] ointT[k][m][n][v];
                    }
                }
            }
            for (int m=0; m<NBAND_C[k]; m++){
                delete [] eigvector_kq[k][m];
            }
            delete [] eigvector_kq[k];
        }
    }// loop q

    // free dynamically allocated memory
    for (int k=0; k<NKPT; k++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                delete [] ointL[k][m][n];
                delete [] ointT[k][m][n];
            }
            delete [] ointL[k][m];
            delete [] ointT[k][m];
        }
        delete [] ointL[k];
        delete [] ointT[k];

        for (int m=0; m<NBAND_V[k]; m++){
            delete [] eigvector_k[k][m];
        }
        delete [] eigvector_k[k];
        delete [] eigvalue_k[k];
        delete [] eigvalue_kq[k];
    }
    delete [] ointL;
    delete [] ointT;
    delete [] eigvalue_k;
    delete [] eigvalue_kq;
    delete [] eigvector_k;
    delete [] eigvector_kq;

    return;
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
    dat.energyoffset = eigsolver.eigenvalues()(3);
    dat.eigerror = abs(eigsolver.eigenvalues()(2)-eigsolver.eigenvalues()(1))/2.0;

    return;
}



// return a perpendicular vector
// Eigen::Vector3d perpvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (abs(v(1))<1e-10 && abs(v(2))<1e-10){
//         p << 0, 1, 1;
//         p = p/p.norm();
//     }else if (abs(v(1)-v(2))<1e-10){
//         p << 1, -v(0)/(v(1)+v(2)), -v(0)/(v(1)+v(2));
//         p = p/p.norm();
//     }else if (abs(v(1)+v(2))<1e-10){
//         p << 1, -v(0)/(v(1)-v(2)), v(0)/(v(1)-v(2));
//         p = p/p.norm();
//     }else if (v(2)!=0){
//         p << 0, -v(2), v(1);
//         p = p/p.norm();
//     }else if (v(2)==0 && v(1)!=0){
//         p << -v(1), v(0), 0;
//         p = p/p.norm();
//     }else if (v(1)==0 && v(2)==0){
//         p << 0, 0, 1;
//         p = p/p.norm();
//     }
//     return p;
// }

// Eigen::Vector3d perpvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (v(2)!=0){
//         p << 0, -v(2), v(1);
//         p = p/p.norm();
//     }else if (v(2)==0 && v(1)!=0){
//         p << -v(1), v(0), 0;
//         p = p/p.norm();
//     }else if (v(1)==0 && v(2)==0){
//         p << 0, 0, 1;
//         p = p/p.norm();
//     }
//     return p;
// }

// Eigen::Vector3d perpvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (v(2)!=0){
//         p << -v(2), 0, v(0);
//         p = p/p.norm();
//     }else if (v(2)==0 && v(1)!=0){
//         p << v(1), -v(0), 0;
//         p = p/p.norm();
//     }else if (v(1)==0 && v(2)==0){
//         p << 0, 1, 0;
//         // p = p/p.norm();
//     }
//     return p;
// }

Eigen::Vector3d perpvec(Eigen::Vector3d v){
    Eigen::Vector3d p;
    if (v(1)!=0){
        p << -v(1), v(0), 0;
        p = p/p.norm();
    }else if (v(1)==0 && v(2)!=0){
        p << v(2), 0, -v(0);
        p = p/p.norm();
    }else if (v(1)==0 && v(2)==0){
        p << 0, 0, 1;
        // p = p/p.norm();
    }
    return p;
}


// dirac delta function
// std::complex<double> diracdelta(std::complex<double> x, double eps){
// 	std::complex<double> delta_x = exp(-x*x/(eps*eps))/(eps*sqrt(pi));
// 	return delta_x;
// }

std::complex<double> diracdelta(std::complex<double> x, double eps){
	std::complex<double> delta_x = eps/pi/(pow(eps,2)+pow(x,2));
	return delta_x;
}




// // kramers-kronig
// void kramerskronig(env &dat){

//     double realepsL1;
//     double realepsT1;
//     double ds = dat.dfreq*2.0/pi;

    
//     for (int m=0; m<NEPS; m++){ for(int n=0; n<=m; n++){
//         // #pragma omp parallel for
//         for (int f=0; f<NFREQ; f++){
//             realepsL1 = 0.0;
//             realepsT1 = 0.0;
//             for (int s=0; s<NFREQ; s++){
//                 if (f!=s){
//                     realepsL1 += dat.freq[s]*dat.imagepsL[m][n][s]/(dat.freq[s]*dat.freq[s]-dat.freq[f]*dat.freq[f]);
//                     realepsT1 += dat.freq[s]*dat.imagepsT[m][n][s]/(dat.freq[s]*dat.freq[s]-dat.freq[f]*dat.freq[f]);
//                 }
//             }
//             realepsL1 *= ds;
//             realepsT1 *= ds;
//             if (m==n){
//                 realepsL1 += 1.0;
//                 realepsT1 += 1.0;
//             }
//             dat.realepsL[m][n][f] = realepsL1;
//             dat.realepsT[m][n][f] = realepsT1;
            
//         }
        
//     }}
//     // #pragma omp barrier

//     return;
// }


void kramerskronigtransform(double *ReX, double *ImX, double *w, double dw){

double ds = dw*2.0/pi;

#pragma omp parallel for
for (int f=0; f<NFREQ; f++){
    double ReXf = 0.0;
    for (int s=0; s<NFREQ; s++){
        if (s!=f)
            ReXf += w[s]*ImX[s]/(w[s]*w[s]-w[f]*w[f]);
    }
    ReX[f] = ds*ReXf;
}
return;
}





//======================================================//


// print electronic bandstructure
void printBand(env &dat){
    for(int q=0; q<dat.nq; q++){
        for(int n=0; n<NBAND; n++){
            std::cout << dat.band[q][n];
            if(n<NBAND-1) std::cout << ",";
        }
        std::cout << std::endl;
    }
    return;
}

// print permittivity
void printEpsilon(env &dat){

    for (int q=0; q<dat.nq; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    std::cout << dat.realepsL[q][m][n][f];
                    if (f<NFREQ-1) std::cout << ",";
                }
                std::cout << std::endl;
            }
        }
    }
    

    for (int q=0; q<dat.nq; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    std::cout << dat.imagepsL[q][m][n][f];
                    if (f<NFREQ-1) std::cout << ",";
                }
                std::cout << std::endl;
            }
        }
    }

    for (int q=0; q<dat.nq; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    std::cout << dat.realepsT[q][m][n][f];
                    if (f<NFREQ-1) std::cout << ",";
                }
                std::cout << std::endl;
            }
        }
    }

    for (int q=0; q<dat.nq; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    std::cout << dat.imagepsT[q][m][n][f];
                    if (f<NFREQ-1) std::cout << ",";
                }
                std::cout << std::endl;
            }
        }
    }

    return;
}

// export
void writeEpsilon(env &dat){
    std::ofstream file("/tmp/tmpepsilon.dat", std::ios::binary);

    

    for (int q=0; q<dat.nq; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    file.write(reinterpret_cast<char*>(&dat.realepsL[q][m][n][f]), sizeof dat.realepsL[q][m][n][f]);
                }
            }
        }
    }

    for (int q=0; q<dat.nq; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    file.write(reinterpret_cast<char*>(&dat.imagepsL[q][m][n][f]), sizeof dat.imagepsL[q][m][n][f]);
                }
            }
        }
    }

    for (int q=0; q<dat.nq; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    file.write(reinterpret_cast<char*>(&dat.realepsT[q][m][n][f]), sizeof dat.realepsT[q][m][n][f]);
                }
            }
        }
    }

    for (int q=0; q<dat.nq; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    file.write(reinterpret_cast<char*>(&dat.imagepsT[q][m][n][f]), sizeof dat.imagepsT[q][m][n][f]);
                }
            }
        }
    }

    file.close();
    return;

}


//======================================================//
//                      polarization                    //
//======================================================//

// void polarization(env &dat){

//     Eigen::MatrixXcd eps_L(NEPS, NEPS), eps_T(NEPS, NEPS); //Global Matrix

//     double **eigenlongreal, **eigenlongimg, **eigentransreal, **eigentransimg;
// 	eigenlongreal = new double*[21];
// 	eigenlongimg = new double*[21];
// 	eigentransreal = new double*[21];
// 	eigentransimg = new double*[21];
// 	for(int q=0; q<dat.nq; q++){ 
//         eigenlongreal[q] = new double[NEPS]; 	
//         eigenlongimg[q] = new double[NEPS]; 	
//         eigentransreal[q] = new double[NEPS]; 	
//         eigentransimg[q] = new double[NEPS]; 	


//         // calculate permittivity
//         permittivity(dat,dat.q[q]);

//         kramerskronig(dat);
//         for (int f=0; f<NFREQ; f++){
//             for (int m=0; m<NEPS; m++){for (int n=0; n<NEPS; n++){
//                 eps_L(m,n) = std::complex(dat.realepsL[f][m][n],dat.imagepsL[f][m][n]);
//                 if (m!=n) eps_L(n,m) = conj(eps_L(m,n));
//             }}
//         }

//         Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigsolverL(eps_L);

//         // for (int m=0; m<NEPS; m++){
//         //     dat.realpolL[f][m] = real(eigsolverL.eigenvalues()(m));
//         //     eigenlongreal[f][m] = real(eigsolverL.eigenvalues()(m));
//         // }
// 	}
// }



}
