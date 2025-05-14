/*
    pmx_chi.cpp | includes susceptibility calculation functinalities of picomax

    Author: Jungho Mun
    Date: June 6, 2024
*/


/*
    header files
*/
#include "pmx_chi.hpp"
#include "pmx_epm.hpp"
#include "pmx_basis.hpp"


/*
    function definitions
*/
namespace pmx{


/*
    dirac delta function, indirect function call
*/
inline double (*diracdelta) (double value);

/*
    Lorentzian dirac delta function

    eps : the smearing parameter of delta function
    delta(x) = eps/pi/(eps^2+x^2)
*/
inline double diracdelta_lorentzian(double x){
    return EPSILON/pi/(pow(EPSILON,2)+pow(x,2));
}

/*
    Gaussian dirac delta function

    eps : the smearing parameter of delta function
    delta(x) = (eps*sqrt(2*pi))^-1 * exp(-0.5*x^2/eps^2)
*/
inline double diracdelta_gaussian(double x){
    return exp(-0.5*x*x/(EPSILON*EPSILON))/(EPSILON*sqrt(2.0*pi));
}

/*
    triangular dirac delta function

    eps : the smearing parameter of delta function
    delta(x) = 4/eps * (0.5 - abs(x)/eps)       (abs(x)<0.5*eps)
             = 0                                (otherwise)
*/
inline double diracdelta_triangle(double x){
    if (abs(x)<0.5*EPSILON){
        return 4.0/EPSILON*(0.5-abs(x)/EPSILON);
    }else{
        return 0;
    }
}

/*
    rectangular dirac delta function

    eps : the smearing parameter of delta function
    delta(x) = 1/eps        (abs(x)<0.5*eps)
             = 0            (otherwise)
*/
inline double diracdelta_rect(double x){
    if (abs(x)<0.5*EPSILON){
        return 1.0/EPSILON;
    }else{
        return 0;
    }
}

/*
    Hilbert transformation of susceptibility tensor matrix elements

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions
    [2] Phys. Rev. B 74, 035101 (2006)
*/
void kramerskronigtransform(double *ReX, double *ImX, double *w, double dw){

    // copy real part
    double ReX1 [NFREQ];
    for (int f=0; f<NFREQ; f++)
        ReX1[f] = ReX[f];

    // 
    double ds = dw*2.0/pi;

    // generate real part from imaginary part
    for (int f=0; f<NFREQ; f++){
        double ReXf = 0.0;
        for (int s=0; s<NFREQ; s++){
            if (s!=f)
                ReXf += w[s]*ImX[s]/(w[s]*w[s]-w[f]*w[f]);
        }
        ReX[f] += ds*ReXf;
    }

    // generate imaginary part from real part
    for (int f=0; f<NFREQ; f++){
        double ImXf = 0.0;
        for (int s=0; s<NFREQ; s++){
            if (s!=f)
                ImXf += w[s]*ReX1[s]/(w[s]*w[s]-w[f]*w[f]);
        }
        ImX[f] += ds*ImXf;
    }

    return;
}

/*
    calculate the susceptibility tensor matrix using an alternative expression
    which treats the longitudinal terms as charge

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions
    [2] (Adler) “Quantum Theory of the Dielectric Constant in Real Solids”, Phys. Rev. 126, 413 (1962)
*/
void chi_LL(env &dat){
    std::chrono::time_point<std::chrono::system_clock> time_0, time_1;
    std::chrono::duration<double> elapsed_time;

    NTSR = 3;

    // use indirect function calls for diracdelta function
    switch (dat.delta) {
        case 0:
            diracdelta = &diracdelta_gaussian;
            break;
        case 1:
            diracdelta = &diracdelta_lorentzian;
            break;
        case 2:
            diracdelta = &diracdelta_triangle;
            break;
        case 3:
            diracdelta = &diracdelta_rect;
            break;
    }

    // scalefactor
    double SCALEFACTOR = 2.0*dat.lat.bz_volume*pow(e,2)/(pi*eV*a*angstrom);
    double SCALEFACTOR_LL = 2.0*dat.lat.bz_volume*pow(e,2)/(pi*eV*a*angstrom);
    double SCALEFACTOR_TT = 32.0*dat.lat.bz_volume*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));
    double SCALEFACTOR_LT = 8.0*dat.lat.bz_volume*pi*pow(hbar,2)*pow(e,2)/(pow(eV,2)*e_m*pow(a*angstrom,3));

    // initialize susceptibility tensor matrix elements
    dat.ReXij = new double *****[NQ];
    dat.ImXij = new double *****[NQ];
    for (int q=0; q<NQ; q++){
        dat.ReXij[q] = new double ****[NTSR];
        dat.ImXij[q] = new double ****[NTSR];
        for (int i=0; i<NTSR; i++){
            dat.ReXij[q][i] = new double ***[NTSR];
            dat.ImXij[q][i] = new double ***[NTSR];
            for (int j=0; j<NTSR; j++){
                dat.ReXij[q][i][j] = new double **[NEPS];
                dat.ImXij[q][i][j] = new double **[NEPS];
                for (int m=0; m<NEPS; m++){
                    dat.ReXij[q][i][j][m] = new double *[NEPS];
                    dat.ImXij[q][i][j][m] = new double *[NEPS];
                    for (int n=0; n<=m; n++){
                        dat.ReXij[q][i][j][m][n] = new double [NFREQ];
                        dat.ImXij[q][i][j][m][n] = new double [NFREQ];
                    }
                }
            }
        }
    }

    //======================================================================
    std::cout << "  Solving for |k,v>..." << std::endl;
    time_0 = std::chrono::system_clock::now();
    int NBAND_C [NKPT]; // number of conduction bands
    double **E_k; // eigenvalue at k [NK x NC]
    E_k = new double *[NKPT];
    std::complex<double> ***C_k; // eigenvector at k [NK x NC]
    C_k = new std::complex<double> **[NKPT];
    #pragma omp parallel for
    for (int k=0; k<NKPT; k++){
        // vector{k}
        Eigen::Vector3d K = dat.lat.K[k];

        // Hamiltonian at k
        Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K,dat.lat.atomic,dat.mat);   
        
        // compute eigenvalue problem at k
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);

        E_k[k] = new double [NBAND];
        for (int c=0; c<NBAND; c++){
            E_k[k][c] = eigsolver_k.eigenvalues()(NBAND-1-c)-dat.energyoffset;
            if (E_k[k][c]<=0) {
                NBAND_C[k] = c;
                break;
            }
        }
        C_k[k] = new std::complex<double> *[NBAND_C[k]];
        for (int c=0; c<NBAND_C[k]; c++){
            C_k[k][c] = new std::complex<double> [NPW];
            for (int p=0; p<NPW; p++){
                C_k[k][c][p] = eigsolver_k.eigenvectors().col(NBAND-1-c)(p);
            }// loop over p
        }// loop over c
    }
    #pragma omp barrier
    time_1 = std::chrono::system_clock::now();
    elapsed_time = time_1-time_0;
    std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;

    //======================================================================
    // initialize large temporary variables
    // overlap integral, t_{Q+G_m} * <k,v|e^{-i*G_m}j_0|k+q,c> <k+q,c|e^{i*G_n}j_0|k,v> * t_{Q+G_n}
    std::complex<double> *****oint;
    oint = new std::complex<double> ****[NKPT]; // overlap integral [NK x NTSR x NEPS x NC x NV]
    for (int k=0; k<NKPT; k++){
        oint[k] = new std::complex<double> ***[3];
        for (int i=0; i<NTSR; i++){
            oint[k][i] = new std::complex<double> **[NEPS];
            for (int m=0; m<NEPS; m++){
                oint[k][i][m] = new std::complex<double> *[NBAND_C[k]];
            }
        }
    }

    std::cout << "  Solving for susceptibility tensor matrix elements (alternative form)..." << std::endl;
    int NBAND_V [NKPT]; // number of valence bands
    std::complex<double> ***C_kq; // eigenvector at k+q
    double **E_kq; // eigenvalue at k
    E_kq = new double *[NKPT];
    C_kq = new std::complex<double> **[NKPT];
    for (int k=0; k<NKPT; k++){
        E_kq[k] = new double [NBAND];
    }
    // #pragma omp parallel
    for (int q=0; q<NQ; q++){
        std::cout << "  (" << (q+1) << "/" << NQ << ")..." << std::endl; // (q/NQ)
        Eigen::Vector3d Q = dat.lat.Q[q];

        std::cout << "    Solving for overlap integrals..." << std::endl;
        time_0 = std::chrono::system_clock::now();
        #pragma omp parallel for
        for (int k=0; k<NKPT; k++){
            
            // Hamiltonian at k+q
            Eigen::Vector3d K = dat.lat.K[k];
            Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K+Q,dat.lat.atomic,dat.mat);   

            // compute eigenvalue problem at k+q
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_kq(H);
            for (int v=0; v<NBAND; v++){
                E_kq[k][v] = double(eigsolver_kq.eigenvalues()(v))-dat.energyoffset;
                if (E_kq[k][v]>0) {
                    NBAND_V[k] = v;
                    break;
                }
            }
            C_kq[k] = new std::complex<double> *[NBAND_V[k]];
            for (int v=0; v<NBAND_V[k]; v++){
                C_kq[k][v] = new std::complex<double> [NPW];
                for (int p=0; p<NPW; p++){
                    C_kq[k][v][p] = eigsolver_kq.eigenvectors().col(v)(p);
                }
            }
            
            // overlap integrals
            for (int m=0; m<NEPS; m++){
                std::vector<Eigen::Vector3cd> uvec_m = uvec_LT(Q+dat.lat.G[m]);
                for (int i=0; i<NTSR; i++){
                    for (int c=0; c<NBAND_C[k]; c++){
                        oint[k][i][m][c] = new std::complex<double> [NBAND_V[k]];
                        for (int v=0; v<NBAND_V[k]; v++){
                            oint[k][i][m][c][v] = 0;
                            if (i==0){// L, <k,c|e^{-i*(q+g_m)*r}|k+q,v>
                                for (int p=0; p<NPW; p++){if (dat.lat.loci[m][p]!=-1){
                                    oint[k][i][m][c][v] += (conj(C_k[k][c][p]) * C_kq[k][v][dat.lat.loci[m][p]]);
                                }}
                            }else{// T, u^i_{q+g_m} * <k,c|e^{-i*(q+g_m)*r} \hat{j}_0 |k+q,v>
                                for (int p=0; p<NPW; p++){if (dat.lat.loci[m][p]!=-1){
                                    oint[k][i][m][c][v] += (conj(C_k[k][c][p]) * C_kq[k][v][dat.lat.loci[m][p]]) 
                                                * (uvec_m[i].dot(dat.lat.G[p]+K+Q/2+dat.lat.G[m]/2));
                                }}
                            }
                        }//loop over v
                    }//loop over c
                }//loop over i
            }//loop over m
        }//loop over k
        #pragma omp barrier
        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;

        std::cout << "    Solving for Xijmn..." << std::endl;
        time_0 = std::chrono::system_clock::now();
        #pragma omp parallel for collapse(4)
        for (int i=0; i<NTSR; i++){for (int j=0; j<NTSR; j++){for (int m=0; m<NEPS; m++){for (int n=0; n<=m; n++){
            Eigen::Vector3d vecqgm = Q+dat.lat.G[m];
            Eigen::Vector3d vecqgn = Q+dat.lat.G[n];
            for (int f=0; f<NFREQ; f++){
                double tmp_real_1 = 0;
                double tmp_imag_1 = 0;
                double tmp_real_2 = 0;
                double tmp_imag_2 = 0;
                // sum over k,c,v
                for (int k=0; k<NKPT; k++){for (int c=0; c<NBAND_C[k]; c++){for (int v=0; v<NBAND_V[k]; v++){
                    double dE = E_k[k][c]-E_kq[k][v];
                    std::complex<double> Oij = oint[k][i][m][c][v] * conj(oint[k][j][n][c][v]);
                    
                    // 
                    tmp_imag_1 += dat.lat.KW[k] * Oij.real()
                                          * (*diracdelta)(dE-dat.freq[f]);
                    tmp_real_1 -= dat.lat.KW[k] * Oij.imag()
                                          * (*diracdelta)(dE-dat.freq[f]);

                    if (dat.kk==0){
                        if (i==0 && j==0){// LL
                            tmp_real_2 += dat.lat.KW[k] * Oij.real()
                                         * (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                            tmp_imag_2 += dat.lat.KW[k] * Oij.imag()
                                         * (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                        }else if (i==0 || j==0){// LT
                            tmp_real_2 += dat.lat.KW[k] * Oij.real()
                                          / (dE*dE-dat.freq[f]*dat.freq[f]);
                            tmp_imag_2 += dat.lat.KW[k] * Oij.imag()
                                          / (dE*dE-dat.freq[f]*dat.freq[f]);
                        }else{// TT
                            tmp_real_2 += dat.lat.KW[k] * Oij.real()
                                         / (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                            tmp_imag_2 += dat.lat.KW[k] * Oij.imag()
                                         / (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                        }
                    }
                }}}
                if (f==0){
                    tmp_imag_1 = 0;
                    tmp_real_1 = 0;
                }else{
                    if (i==0 && j==0){
                        tmp_imag_1 *= SCALEFACTOR_LL / (vecqgm.norm()*vecqgn.norm());
                        tmp_real_1 *= SCALEFACTOR_LL / (vecqgm.norm()*vecqgn.norm());
                    }else if (i==0){
                        tmp_imag_1 *= SCALEFACTOR_LT / (vecqgm.norm()*dat.freq[f]);
                        tmp_real_1 *= SCALEFACTOR_LT / (vecqgm.norm()*dat.freq[f]);
                    }else if (j==0){
                        tmp_imag_1 *= SCALEFACTOR_LT / (vecqgn.norm());
                        tmp_real_1 *= SCALEFACTOR_LT / (vecqgn.norm());
                    }else{
                        tmp_imag_1 *= SCALEFACTOR_TT / (dat.freq[f]*dat.freq[f]);
                        tmp_real_1 *= SCALEFACTOR_TT / (dat.freq[f]*dat.freq[f]);
                    }
                }
                if (dat.kk==0){
                    if (i==0 && j==0){
                        tmp_imag_2 *= SCALEFACTOR_LL / (vecqgm.norm()*vecqgn.norm());
                        tmp_real_2 *= SCALEFACTOR_LL / (vecqgm.norm()*vecqgn.norm());
                    }else if (i==0){
                        tmp_imag_2 *= SCALEFACTOR_LT / (vecqgm.norm());
                        tmp_real_2 *= SCALEFACTOR_LT / (vecqgm.norm());
                    }else if (j==0){
                        tmp_imag_2 *= SCALEFACTOR_LT / (vecqgn.norm());
                        tmp_real_2 *= SCALEFACTOR_LT / (vecqgn.norm());
                    }else{
                        tmp_imag_2 *= SCALEFACTOR_TT;
                        tmp_real_2 *= SCALEFACTOR_TT;
                    }
                    
                }

                dat.ImXij[q][i][j][m][n][f] = (pi*tmp_imag_1+2.0*tmp_imag_2);
                dat.ReXij[q][i][j][m][n][f] = (pi*tmp_real_1+2.0*tmp_real_2);
            }
        }}}}
        #pragma omp barrier
        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;
        
        // kramer-kronig transform
        if(dat.kk){
            std::cout << "    Kramers-Kronig transform..." << std::endl;
            time_0 = std::chrono::system_clock::now();
            #pragma omp parallel for collapse(4)
            for (int i=0; i<NTSR; i++){ for (int j=0; j<NTSR; j++){ for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
                kramerskronigtransform(dat.ReXij[q][i][j][m][n],dat.ImXij[q][i][j][m][n],dat.freq,dat.dfreq);
            }}}}
            #pragma omp barrier
            time_1 = std::chrono::system_clock::now();
            elapsed_time = time_1-time_0;
            std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;
        }
        

        // free dynamically allocated memory
        for (int k=0; k<NKPT; k++){for (int i=0; i<NTSR; i++){for (int m=0; m<NEPS; m++){for (int c=0; c<NBAND_C[k]; c++){
            delete [] oint[k][i][m][c];
        }}}}
        for (int k=0; k<NKPT; k++){
            for (int v=0; v<NBAND_V[k]; v++){
                delete [] C_kq[k][v];
            }
            delete [] C_kq[k];
        }
    }// loop q

    // free dynamically allocated memory
    for (int k=0; k<NKPT; k++){
        for (int i=0; i<NTSR; i++){
            for (int m=0; m<NEPS; m++){
                delete [] oint[k][i][m];
            }
            delete [] oint[k][i];
        }
        delete [] oint[k];
    }
    delete [] oint;

    for (int k=0; k<NKPT; k++){
        for (int c=0; c<NBAND_C[k]; c++){
            delete [] C_k[k][c];
        }
        delete [] C_k[k];
        delete [] E_k[k];
        delete [] E_kq[k];
    }
    delete [] E_k;
    delete [] E_kq;
    delete [] C_k;
    delete [] C_kq;

    return;
}

/*
    calculate the susceptibility tensor matrix using the current-current correlation function

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions
    [2] (Adler) “Quantum Theory of the Dielectric Constant in Real Solids”, Phys. Rev. 126, 413 (1962)
    [3] 
*/
void chi_tensor(env &dat){
    std::chrono::time_point<std::chrono::system_clock> time_0, time_1;
    std::chrono::duration<double> elapsed_time;

    // use indirect function calls for diracdelta function
    switch (dat.delta) {
        case 0:
            diracdelta = &diracdelta_gaussian;
            break;
        case 1:
            diracdelta = &diracdelta_lorentzian;
            break;
        case 2:
            diracdelta = &diracdelta_triangle;
            break;
        case 3:
            diracdelta = &diracdelta_rect;
            break;
    }
    NTSR = 3;

    // scalefactor
    double SCALEFACTOR = 32.0*dat.lat.bz_volume*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));

    // initialize susceptibility tensor matrix elements
    dat.ReXij = new double *****[NQ];
    dat.ImXij = new double *****[NQ];
    for (int q=0; q<NQ; q++){
        dat.ReXij[q] = new double ****[NTSR];
        dat.ImXij[q] = new double ****[NTSR];
        for (int i=0; i<NTSR; i++){
            dat.ReXij[q][i] = new double ***[NTSR];
            dat.ImXij[q][i] = new double ***[NTSR];
            for (int j=0; j<NTSR; j++){
                dat.ReXij[q][i][j] = new double **[NEPS];
                dat.ImXij[q][i][j] = new double **[NEPS];
                for (int m=0; m<NEPS; m++){
                    dat.ReXij[q][i][j][m] = new double *[NEPS];
                    dat.ImXij[q][i][j][m] = new double *[NEPS];
                    for (int n=0; n<=m; n++){
                        dat.ReXij[q][i][j][m][n] = new double [NFREQ];
                        dat.ImXij[q][i][j][m][n] = new double [NFREQ];
                    }
                }
            }
        }
    }

    //======================================================================
    std::cout << "  Solving for |k,v>..." << std::endl;
    time_0 = std::chrono::system_clock::now();
    int NBAND_C [NKPT]; // number of conduction bands
    double **E_k; // eigenvalue at k [NK x NC]
    E_k = new double *[NKPT];
    std::complex<double> ***C_k; // eigenvector at k [NK x NC x NPW]
    C_k = new std::complex<double> **[NKPT];
    #pragma omp parallel for
    for (int k=0; k<NKPT; k++){
        // vector{k}
        Eigen::Vector3d K = dat.lat.K[k];

        // Hamiltonian at k
        Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K,dat.lat.atomic,dat.mat);   
        
        // compute eigenvalue problem at k
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);

        E_k[k] = new double [NBAND];
        for (int c=0; c<NBAND; c++){
            E_k[k][c] = eigsolver_k.eigenvalues()(NBAND-1-c)-dat.energyoffset;
            if (E_k[k][c]<=0) {
                NBAND_C[k] = c;
                break;
            }
        }
        C_k[k] = new std::complex<double> *[NBAND_C[k]];
        for (int c=0; c<NBAND_C[k]; c++){
            C_k[k][c] = new std::complex<double> [NPW];
            for (int p=0; p<NPW; p++){
                C_k[k][c][p] = eigsolver_k.eigenvectors().col(NBAND-1-c)(p);
            }// loop over p
        }// loop over c
    }
    #pragma omp barrier
    time_1 = std::chrono::system_clock::now();
    elapsed_time = time_1-time_0;
    std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;

    //======================================================================
    // initialize large temporary variables
    // overlap integral, t_{Q+G_m} * <k,v|e^{-i*G_m}j_0|k+q,c> <k+q,c|e^{i*G_n}j_0|k,v> * t_{Q+G_n}
    std::complex<double> *****oint;
    oint = new std::complex<double> ****[NKPT]; // overlap integral [NK x NTSR x NEPS x NC x NV]
    for (int k=0; k<NKPT; k++){
        oint[k] = new std::complex<double> ***[3];
        for (int i=0; i<NTSR; i++){
            oint[k][i] = new std::complex<double> **[NEPS];
            for (int m=0; m<NEPS; m++){
                oint[k][i][m] = new std::complex<double> *[NBAND_C[k]];
            }
        }
    }

    std::cout << "  Solving for susceptibility tensor matrix elements..." << std::endl;
    int NBAND_V [NKPT]; // number of valence bands
    std::complex<double> ***C_kq; // eigenvector at k+q
    double **E_kq; // eigenvalue at k
    E_kq = new double *[NKPT];
    C_kq = new std::complex<double> **[NKPT];
    for (int k=0; k<NKPT; k++){
        E_kq[k] = new double [NBAND];
    }
    // #pragma omp parallel
    for (int q=0; q<NQ; q++){
        std::cout << "  (" << (q+1) << "/" << NQ << ")..." << std::endl; // (q/NQ)
        Eigen::Vector3d Q = dat.lat.Q[q];

        std::cout << "    Solving for overlap integrals..." << std::endl;
        time_0 = std::chrono::system_clock::now();
        #pragma omp parallel for
        for (int k=0; k<NKPT; k++){
            
            // Hamiltonian at k+q
            Eigen::Vector3d K = dat.lat.K[k];
            Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K+Q,dat.lat.atomic,dat.mat);   

            // compute eigenvalue problem at k+q
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_kq(H);
            for (int v=0; v<NBAND; v++){
                E_kq[k][v] = double(eigsolver_kq.eigenvalues()(v))-dat.energyoffset;
                if (E_kq[k][v]>0) {
                    NBAND_V[k] = v;
                    break;
                }
            }
            C_kq[k] = new std::complex<double> *[NBAND_V[k]];
            for (int v=0; v<NBAND_V[k]; v++){
                C_kq[k][v] = new std::complex<double> [NPW];
                for (int p=0; p<NPW; p++){
                    C_kq[k][v][p] = eigsolver_kq.eigenvectors().col(v)(p);
                }
            }
            
            // overlap integrals
            for (int m=0; m<NEPS; m++){
                std::vector<Eigen::Vector3cd> uvec_m = uvec_LT(Q+dat.lat.G[m]);
                for (int i=0; i<NTSR; i++){
                    for (int c=0; c<NBAND_C[k]; c++){
                        oint[k][i][m][c] = new std::complex<double> [NBAND_V[k]];
                        for (int v=0; v<NBAND_V[k]; v++){
                            // u^i_{q+g_m} * <k,c|e^{-i*(q+g_m)*r}j_0|k+q,c>
                            oint[k][i][m][c][v] = 0;
                            for (int p=0; p<NPW; p++){if (dat.lat.loci[m][p]!=-1){
                                oint[k][i][m][c][v] += (conj(C_k[k][c][p]) * C_kq[k][v][dat.lat.loci[m][p]]) 
                                            * (uvec_m[i].dot(dat.lat.G[p]+K+Q/2+dat.lat.G[m]/2));
                            }}
                        }//loop over v
                    }//loop over c
                }//loop over i
            }//loop over m
        }//loop over k
        #pragma omp barrier
        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;

        std::cout << "    Solving for Xijmn..." << std::endl;
        time_0 = std::chrono::system_clock::now();
        #pragma omp parallel for collapse(4)
        for (int i=0; i<NTSR; i++){for (int j=0; j<NTSR; j++){for (int m=0; m<NEPS; m++){for (int n=0; n<=m; n++){
            for (int f=0; f<NFREQ; f++){
                double tmp_real_1 = 0;
                double tmp_imag_1 = 0;
                double tmp_real_2 = 0;
                double tmp_imag_2 = 0;
                // sum over k,c,v
                for (int k=0; k<NKPT; k++){for (int c=0; c<NBAND_C[k]; c++){for (int v=0; v<NBAND_V[k]; v++){
                    double dE = E_k[k][c]-E_kq[k][v];
                    std::complex<double> Oij = oint[k][i][m][c][v] * conj(oint[k][j][n][c][v]);
                    
                    // 
                    tmp_imag_1 += dat.lat.KW[k] * Oij.real()
                                          * (*diracdelta)(dE-dat.freq[f]);
                    tmp_real_1 -= dat.lat.KW[k] * Oij.imag()
                                          * (*diracdelta)(dE-dat.freq[f]);

                    if (dat.kk==0){
                        tmp_real_2 += dat.lat.KW[k] * Oij.real()
                                         / (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                        tmp_imag_2 += dat.lat.KW[k] * Oij.imag()
                                         / (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                    }
                }}}
                if (f==0){
                    tmp_imag_1 = 0;
                    tmp_real_1 = 0;
                }else{
                    tmp_imag_1 /= (dat.freq[f]*dat.freq[f]);
                    tmp_real_1 /= (dat.freq[f]*dat.freq[f]);
                }
                dat.ImXij[q][i][j][m][n][f] = SCALEFACTOR * (pi*tmp_imag_1+2*tmp_imag_2);
                dat.ReXij[q][i][j][m][n][f] = SCALEFACTOR * (pi*tmp_real_1+2*tmp_real_2);
            }
        }}}}
        #pragma omp barrier
        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;


        
        // kramer-kronig transform
        if(dat.kk){
            std::cout << "    Kramers-Kronig transform..." << std::endl;
            time_0 = std::chrono::system_clock::now();
            #pragma omp parallel for collapse(4)
            for (int i=0; i<NTSR; i++){ for (int j=0; j<NTSR; j++){ for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
                kramerskronigtransform(dat.ReXij[q][i][j][m][n],dat.ImXij[q][i][j][m][n],dat.freq,dat.dfreq);
            }}}}
            #pragma omp barrier
            time_1 = std::chrono::system_clock::now();
            elapsed_time = time_1-time_0;
            std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;
        }
        

        // free dynamically allocated memory
        for (int k=0; k<NKPT; k++){for (int i=0; i<NTSR; i++){for (int m=0; m<NEPS; m++){for (int c=0; c<NBAND_C[k]; c++){
            delete [] oint[k][i][m][c];
        }}}}
        for (int k=0; k<NKPT; k++){
            for (int v=0; v<NBAND_V[k]; v++){
                delete [] C_kq[k][v];
            }
            delete [] C_kq[k];
        }
    }// loop q

    // free dynamically allocated memory
    for (int k=0; k<NKPT; k++){
        for (int i=0; i<NTSR; i++){
            for (int m=0; m<NEPS; m++){
                delete [] oint[k][i][m];
            }
            delete [] oint[k][i];
        }
        delete [] oint[k];
    }
    delete [] oint;

    for (int k=0; k<NKPT; k++){
        for (int c=0; c<NBAND_C[k]; c++){
            delete [] C_k[k][c];
        }
        delete [] C_k[k];
        delete [] E_k[k];
        delete [] E_kq[k];
    }
    delete [] E_k;
    delete [] E_kq;
    delete [] C_k;
    delete [] C_kq;

    return;
}

/*
    calculate the susceptibility tensor matrix using the current-current correlation function,
    where the polarization unit vectors are given by the Cartesian unit vectors {nx,ny,nz}

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions
    [2] (Adler) “Quantum Theory of the Dielectric Constant in Real Solids”, Phys. Rev. 126, 413 (1962)
*/
void chi_xyz(env &dat){
    std::chrono::time_point<std::chrono::system_clock> time_0, time_1;
    std::chrono::duration<double> elapsed_time;


    // use indirect function calls for diracdelta function
    switch (dat.delta) {
        case 0:
            diracdelta = &diracdelta_gaussian;
            break;
        case 1:
            diracdelta = &diracdelta_lorentzian;
            break;
        case 2:
            diracdelta = &diracdelta_triangle;
            break;
        case 3:
            diracdelta = &diracdelta_rect;
            break;
    }
    NTSR = 3;

    // scalefactor
    double SCALEFACTOR = 32.0*dat.lat.bz_volume*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));

    // initialize susceptibility tensor matrix elements
    dat.ReXij = new double *****[NQ];
    dat.ImXij = new double *****[NQ];
    for (int q=0; q<NQ; q++){
        dat.ReXij[q] = new double ****[NTSR];
        dat.ImXij[q] = new double ****[NTSR];
        for (int i=0; i<NTSR; i++){
            dat.ReXij[q][i] = new double ***[NTSR];
            dat.ImXij[q][i] = new double ***[NTSR];
            for (int j=0; j<NTSR; j++){
                dat.ReXij[q][i][j] = new double **[NEPS];
                dat.ImXij[q][i][j] = new double **[NEPS];
                for (int m=0; m<NEPS; m++){
                    dat.ReXij[q][i][j][m] = new double *[NEPS];
                    dat.ImXij[q][i][j][m] = new double *[NEPS];
                    for (int n=0; n<=m; n++){
                        dat.ReXij[q][i][j][m][n] = new double [NFREQ];
                        dat.ImXij[q][i][j][m][n] = new double [NFREQ];
                        // for (int f=0; f<NFREQ; f++){
                        //     dat.ReXij[q][i][j][m][n][f] = 0;
                        //     dat.ImXij[q][i][j][m][n][f] = 0;
                        // }
                    }
                }
            }
        }
    }

    //======================================================================
    std::cout << "  Solving for |k,v>..." << std::endl;
    time_0 = std::chrono::system_clock::now();
    int NBAND_V [NKPT]; // number of conduction bands
    double **E_k; // eigenvalue at k [NK x NC]
    E_k = new double *[NKPT];
    std::complex<double> ***C_k; // eigenvector at k [NK x NC]
    C_k = new std::complex<double> **[NKPT];
    #pragma omp parallel for
    for (int k=0; k<NKPT; k++){
        // vector{k}
        Eigen::Vector3d K = dat.lat.K[k];

        // Hamiltonian at k
        Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K,dat.lat.atomic,dat.mat);   
        
        // compute eigenvalue problem at k
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);

        E_k[k] = new double [NBAND];
        for (int v=0; v<NBAND; v++){
            E_k[k][v] = eigsolver_k.eigenvalues()(v)-dat.energyoffset;
            if (E_k[k][v]>0) {
                NBAND_V[k] = v;
                break;
            }
        }
        C_k[k] = new std::complex<double> *[NBAND_V[k]];
        for (int v=0; v<NBAND_V[k]; v++){
            C_k[k][v] = new std::complex<double> [NPW];
            for (int p=0; p<NPW; p++){
                C_k[k][v][p] = eigsolver_k.eigenvectors().col(v)(p);
            }// loop over p
        }// loop over v
    }
    #pragma omp barrier
    time_1 = std::chrono::system_clock::now();
    elapsed_time = time_1-time_0;
    std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;

    //======================================================================
    // initialize large temporary variables
    // overlap integral, t_{Q+G_m} * <k,v|e^{-i*G_m}j_0|k+q,c> <k+q,c|e^{i*G_n}j_0|k,v> * t_{Q+G_n}
    std::complex<double> *****oint;
    oint = new std::complex<double> ****[NKPT]; // overlap integral [NK x NTSR x NEPS x NC x NV]
    for (int k=0; k<NKPT; k++){
        oint[k] = new std::complex<double> ***[3];
        for (int i=0; i<NTSR; i++){
            oint[k][i] = new std::complex<double> **[NEPS];
            for (int m=0; m<NEPS; m++){
                oint[k][i][m] = new std::complex<double> *[NBAND_V[k]];
            }
        }
    }

    std::cout << "  Solving for susceptibility tensor matrix elements (xyz-basis)..." << std::endl;
    int NBAND_C [NKPT]; // number of valence bands
    std::complex<double> ***C_kq; // eigenvector at k+q
    double **E_kq; // eigenvalue at k
    E_kq = new double *[NKPT];
    C_kq = new std::complex<double> **[NKPT];
    for (int k=0; k<NKPT; k++){
        E_kq[k] = new double [NBAND];
    }
    // #pragma omp parallel
    for (int q=0; q<NQ; q++){
        std::cout << "  (" << (q+1) << "/" << NQ << ")..." << std::endl; // (q/NQ)
        Eigen::Vector3d Q = dat.lat.Q[q];

        std::cout << "    Solving for overlap integrals..." << std::endl;
        time_0 = std::chrono::system_clock::now();
        #pragma omp parallel for
        for (int k=0; k<NKPT; k++){
            
            // Hamiltonian at k+q
            Eigen::Vector3d K = dat.lat.K[k];
            Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K+Q,dat.lat.atomic,dat.mat);   

            // compute eigenvalue problem at k+q
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_kq(H);
            for (int c=0; c<NBAND; c++){
                E_kq[k][c] = double(eigsolver_kq.eigenvalues()(NBAND-1-c))-dat.energyoffset;
                if (E_kq[k][c]<=0) {
                    NBAND_C[k] = c;
                    break;
                }
            }
            C_kq[k] = new std::complex<double> *[NBAND_C[k]];
            for (int c=0; c<NBAND_C[k]; c++){
                C_kq[k][c] = new std::complex<double> [NPW];
                for (int p=0; p<NPW; p++){
                    C_kq[k][c][p] = eigsolver_kq.eigenvectors().col(NBAND-1-c)(p);
                }
            }
            
            std::vector<Eigen::Vector3d> uvec;
            uvec.reserve(3);
            uvec[0] = {1,0,0};
            uvec[1] = {0,1,0};
            uvec[2] = {0,0,1};

            // overlap integrals
            for (int m=0; m<NEPS; m++){
                for (int i=0; i<NTSR; i++){
                    for (int v=0; v<NBAND_V[k]; v++){
                        oint[k][i][m][v] = new std::complex<double> [NBAND_C[k]];
                        for (int c=0; c<NBAND_C[k]; c++){
                            // u^i_{q+g_m} * <k,v|e^{-i*(q+g_m)*r}j_0|k+q,c>
                            oint[k][i][m][v][c] = 0;
                            for (int p=0; p<NPW; p++){if (dat.lat.loci[m][p]!=-1){
                                oint[k][i][m][v][c] += (conj(C_k[k][v][p]) * C_kq[k][c][dat.lat.loci[m][p]]) 
                                            * (uvec[i].dot(dat.lat.G[p]+K+Q/2+dat.lat.G[m]/2));
                            }}
                        }//loop over v
                    }//loop over c
                }//loop over i
            }//loop over m
        }//loop over k
        #pragma omp barrier
        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;

        std::cout << "    Solving for Xijmn..." << std::endl;
        time_0 = std::chrono::system_clock::now();
        #pragma omp parallel for collapse(4)
        for (int i=0; i<NTSR; i++){for (int j=0; j<NTSR; j++){for (int m=0; m<NEPS; m++){for (int n=0; n<=m; n++){
            for (int f=0; f<NFREQ; f++){
                double tmp_real_1 = 0;
                double tmp_imag_1 = 0;
                double tmp_real_2 = 0;
                double tmp_imag_2 = 0;
                // sum over k,c,v
                for (int k=0; k<NKPT; k++){for (int v=0; v<NBAND_V[k]; v++){for (int c=0; c<NBAND_C[k]; c++){
                    double dE = E_kq[k][c]-E_k[k][v];
                    std::complex<double> Oij = oint[k][i][m][v][c] * conj(oint[k][j][n][v][c]);
                    
                    // 
                    tmp_imag_1 += dat.lat.KW[k] * Oij.real()
                                          * (*diracdelta)(dE-dat.freq[f]);
                    tmp_real_1 -= dat.lat.KW[k] * Oij.imag()
                                          * (*diracdelta)(dE-dat.freq[f]);

                    if (dat.kk==0){
                        tmp_real_2 += dat.lat.KW[k] * Oij.real()
                                         / (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                        tmp_imag_2 += dat.lat.KW[k] * Oij.imag()
                                         / (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                    }
                }}}
                if (f==0){
                    tmp_imag_1 = 0;
                    tmp_real_1 = 0;
                }else{
                    tmp_imag_1 /= (dat.freq[f]*dat.freq[f]);
                    tmp_real_1 /= (dat.freq[f]*dat.freq[f]);
                }
                dat.ImXij[q][i][j][m][n][f] = SCALEFACTOR * (pi*tmp_imag_1+2*tmp_imag_2);
                dat.ReXij[q][i][j][m][n][f] = SCALEFACTOR * (pi*tmp_real_1+2*tmp_real_2);
            }
        }}}}
        #pragma omp barrier
        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;


        
        // kramer-kronig transform
        if(dat.kk){
            std::cout << "    Kramers-Kronig transform..." << std::endl;
            time_0 = std::chrono::system_clock::now();
            #pragma omp parallel for collapse(4)
            for (int i=0; i<NTSR; i++){ for (int j=0; j<NTSR; j++){ for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
                kramerskronigtransform(dat.ReXij[q][i][j][m][n],dat.ImXij[q][i][j][m][n],dat.freq,dat.dfreq);
            }}}}
            #pragma omp barrier
            time_1 = std::chrono::system_clock::now();
            elapsed_time = time_1-time_0;
            std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;
        }
        

        // free dynamically allocated memory
        for (int k=0; k<NKPT; k++){for (int i=0; i<NTSR; i++){for (int m=0; m<NEPS; m++){for (int v=0; v<NBAND_V[k]; v++){
            delete [] oint[k][i][m][v];
        }}}}
        for (int k=0; k<NKPT; k++){
            for (int c=0; c<NBAND_C[k]; c++){
                delete [] C_kq[k][c];
            }
            delete [] C_kq[k];
        }
    }// loop q

    // free dynamically allocated memory
    for (int k=0; k<NKPT; k++){
        for (int i=0; i<NTSR; i++){
            for (int m=0; m<NEPS; m++){
                delete [] oint[k][i][m];
            }
            delete [] oint[k][i];
        }
        delete [] oint[k];
    }
    delete [] oint;

    for (int k=0; k<NKPT; k++){
        for (int v=0; v<NBAND_V[k]; v++){
            delete [] C_k[k][v];
        }
        delete [] C_k[k];
        delete [] E_k[k];
        delete [] E_kq[k];
    }
    delete [] E_k;
    delete [] E_kq;
    delete [] C_k;
    delete [] C_kq;

    return;
}







// void chi_tensor_LF(env &dat){

//     // use indirect function calls for diracdelta function
//     switch (DELTA) {
//         case 0:
//             diracdelta = &diracdelta_gaussian;
//             break;
//         case 1:
//             diracdelta = &diracdelta_lorentzian;
//             break;
//         case 2:
//             diracdelta = &diracdelta_triangle;
//             break;
//         case 3:
//             diracdelta = &diracdelta_rect;
//             break;
//     }

//     // scalefactor
//     double SCALEFACTOR = 32.0*dat.lat.bz_volume*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));
    
//     // initialize susceptibility tensor matrix
//     dat.realepsij = new double *****[3];
//     dat.imagepsij = new double *****[3];
//     for (int i=0; i<3; i++){
//         dat.realepsij[i] = new double ****[3];
//         dat.imagepsij[i] = new double ****[3];
//         for (int j=0; j<3; j++){
//             dat.realepsij[i][j] = new double ***[NQ];
//             dat.imagepsij[i][j] = new double ***[NQ];
//             for (int q=0; q<NQ; q++){
//                 dat.realepsij[i][j][q] = new double **[NEPS];
//                 dat.imagepsij[i][j][q] = new double **[NEPS];
//                 for (int m=0; m<NEPS; m++){
//                     dat.realepsij[i][j][q][m] = new double *[NEPS];
//                     dat.imagepsij[i][j][q][m] = new double *[NEPS];
//                     for (int n=0; n<=m; n++){
//                         dat.realepsij[i][j][q][m][n] = new double [NFREQ];
//                         dat.imagepsij[i][j][q][m][n] = new double [NFREQ];
//                     }
//                 }
//             }
//         }
//     }

//     //======================================================================
//     int NBAND_V [NKPT]; // number of valence bands
//     double **E_k; // Energy at k
//     std::complex<double> ***C_k; // eigenvector at k
//     E_k = new double *[NKPT];
//     C_k = new std::complex<double> **[NKPT];
//     #pragma omp parallel for
//     for (int k=0; k<NKPT; k++){
//         // vector{k}
//         Eigen::Vector3d K = dat.lat.K[k];

//         // Hamiltonian at k
//         Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K,dat.lat.atomic,dat.mat);   
        
//         // compute eigenvalue problem at k
//         Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);
//         if (eigsolver_k.info() != Eigen::Success){ printf("ERROR: Eigensolver failed to solve the matrix."); abort(); }

//         E_k[k] = new double [NBAND];
//         for (int n=0; n<NBAND; n++){
//             E_k[k][n] = eigsolver_k.eigenvalues()(n)-dat.energyoffset;
//             if (E_k[k][n] > 0) {
//                 NBAND_V[k] = n;
//                 break;
//             }
//         }
//         C_k[k] = new std::complex<double> *[NBAND_V[k]];
//         for (int n=0; n<NBAND_V[k]; n++){
//             C_k[k][n] = new std::complex<double> [NPW];
//             for (int j=0; j<NPW; j++){
//                 C_k[k][n][j] = eigsolver_k.eigenvectors().col(n)(j);
//             }
//         }
//     }
//     #pragma omp barrier
    

//     //======================================================================
//     // initialize temporary variables for k+q
//     // overlap integral, O(k,q,m,v,c) = \vec{u}_{Q+G_m} * <k,v|e^{-i*g_m}j_0|k+q+g_m,c>
//     std::complex<double> *****oint;
//     oint = new std::complex<double> ****[NKPT]; // overlap integral [NK x 3 x NEPS x NV x NC]
//     for (int k=0; k<NKPT; k++){
//         oint[k] = new std::complex<double> ***[3];
//         for (int i=0; i<3; i++){
//             oint[k][i] = new std::complex<double> **[NEPS];
//             for (int m=0; m<NEPS; m++){
//                 oint[k][i][m] = new std::complex<double> *[NBAND_V[k]];
//             }
//         }
//     }

//     // initialize E_kqm [NKPT x NEPS x NPW]
//     int NBAND_CB [NKPT]; // number of conduction bands
//     double E_kq [NKPT][NBAND]; // energy at (k+q+g_m,c)
//     // double **E_kqm; // energy at (k+q+g_m,c)
//     // E_kqm = new double *[NKPT];
//     // for (int k=0; k<NKPT; k++){
//     //     E_kqm[k] = new double [NBAND];
//     // }
    

//     for (int q=0; q<NQ; q++){
//         Eigen::Vector3d Q = dat.lat.Q[q];

//         // construct overlap integrals
//         #pragma omp parallel for
//         for (int k=0; k<NKPT; k++){
            
//             // get eigenvectors
//             Eigen::Vector3d K = dat.lat.K[k]; // vector{k}
//             std::complex<double> C_kqgm [NEPS][NBAND][NPW];
//             for (int m=0; m<NEPS; m++){
//                 Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K+Q,dat.lat.atomic,dat.mat);
//                 Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_kq(H);
//                 if (m==0){
//                     for (int c=0; c<NBAND; c++){
//                         E_kq[k][c] = double(eigsolver_kq.eigenvalues()(NBAND-1-c))-dat.energyoffset; // "NBAND-1-c" from behind
//                         if (E_kq[k][c]<=0){
//                             NBAND_CB[k] = c; // number of conduction bands
//                             break;
//                         }
//                     }// loop over c
//                 }
//                 for (int c=0; c<NBAND_CB[k]; c++){
//                     for (int p=0; p<NPW; p++){
//                         C_kqgm[m][c][p] = eigsolver_kq.eigenvectors().col(NBAND-1-c)(p); // "NBAND-1-c" from behind
//                     }
//                 }// loop over c
//             }// loop over m

//             // phase correction (if possible?)



//             // get overlap integrals
//             for (int m=0; m<NEPS; m++){
//                 std::vector<Eigen::Vector3cd> uvec_m = unitvec(dat.lat.G[m]+Q);
//                 for (int i=0; i<3; i++){
//                     for (int v=0; v<NBAND_V[k]; v++){ // valence bands
//                         oint[k][i][m][v] = new std::complex<double> [NBAND_CB[k]];
//                         for (int c=0; c<NBAND_CB[k]; c++){ // conduction bands
//                             std::complex<double> oint_k_kq = 0.0; // <k,v|e^{-i*G_m}j_0|k+q,c> * t_{G_m+Q}
//                             for (int p=0; p<NPW; p++){
//                                 oint_k_kq += conj(C_k[k][v][p]) * C_kqgm[m][c][p] * uvec_m[i].dot(dat.lat.G[p]+K+Q/2+dat.lat.G[m]/2); //+dat.G[m]/2
//                             }// loop over p
//                             oint[k][i][m][v][c] = oint_k_kq;
//                         }// loop over c
//                     }// loop over v
//                 }// loop over i
//             }// loop over m
//         }// loop over k (parallel)
//         #pragma omp barrier



//         // construct the susceptibility tensor matrix
//         for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
//             #pragma omp parallel for //collapse(3)
//             for (int f=0; f<NFREQ; f++){
//                 // temporary variables
//                 double tmpval_real[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
//                 double tmpval_imag[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
//                 double depsilon = 0; // E_{k+q,c}-E_{k,v}
                
//                 for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                    
//                     // sum over k-grids & bands
//                     // \sum_{k,c,v} <k,v|e^{-i*(q+g_m)*r}|k+q+g_m,c> <k+q+g_n,c|e^{i*(q+g_n)*r)}|k,v>
//                     for (int k=0; k<NKPT; k++){
//                         for (int v=0; v<NBAND_V[k]; v++){ for (int c=0; c<NBAND_CB[k]; c++){
//                             depsilon = E_kq[k][c]-E_k[k][v];
//                             // imaginary part
//                             tmpval_imag[i][j] += dat.lat.KW[k] * (oint[k][i][m][v][c] * conj(oint[k][j][n][v][c])).real()
//                                                 * (*diracdelta)(depsilon-dat.freq[f]);
//                             // real part
//                             if (dat.kk==0){
//                                 tmpval_real[i][j] += dat.lat.KW[k] * (oint[k][i][m][v][c] * conj(oint[k][j][n][v][c])).real()
//                                                 / (depsilon) / (depsilon*depsilon-dat.freq[f]*dat.freq[f]);
//                             }
//                         }}
//                     }// loop over k
//                 }}// loop over i,j

//                 // permittivity
//                 for (int i=0; i<3; i++){for (int j=0; j<3; j++){
//                     if (dat.freq[f]>0){
//                         dat.imagepsij[i][j][q][m][n][f] = SCALEFACTOR * pi * tmpval_imag[i][j] / (dat.freq[f]*dat.freq[f]);
//                     }else{
//                         dat.imagepsij[i][j][q][m][n][f] = 0;
//                     }
//                 }}
//                 if (dat.kk==0){
//                     for (int i=0; i<3; i++){for (int j=0; j<3; j++){
//                         dat.realepsij[i][j][q][m][n][f] = SCALEFACTOR * 2 * tmpval_real[i][j];
//                         // if (m==n && i==j) dat.realepsij[i][j][q][m][n][f] += 1;
//                     }}
//                 }
//             }// loop over f
//         }}// loop over m,n
//         #pragma omp barrier

//         // kramer-kronig transform
//         if (dat.kk){
//             for (int i=0; i<3; i++){for (int j=0; j<3; j++){
//                 for (int m=0; m<NEPS; m++){for (int n=0; n<=m; n++){
//                     kramerskronigtransform(dat.realepsij[i][j][q][m][n],dat.imagepsij[i][j][q][m][n],dat.freq,dat.dfreq);
//                     // if (m==n && i==j){
//                     //     for (int f=0; f<NFREQ; f++){
//                     //         dat.realepsij[i][j][q][m][n][f] += 1;
//                     //     }
//                     // }
//                 }}// loop over m,n
//             }}// loop over i,j
//         }

//         // free dynamically allocated memory
//         for (int k=0; k<NKPT; k++){
//             for (int i=0; i<3; i++){
//                 for (int m=0; m<NEPS; m++){
//                     for (int v=0; v<NBAND_V[k]; v++){
//                         delete [] oint[k][i][m][v];
//                     }
//                 }
//             }
//         }
//     }// loop q

//     // free dynamically allocated memory
//     for (int k=0; k<NKPT; k++){
//         for (int i=0; i<3; i++){
//             for (int m=0; m<NEPS; m++){
//                 delete [] oint[k][i][m];
//             }
//             delete [] oint[k][i];
//         }
//         delete [] oint[k];
//     }
//     delete [] oint;

//     for (int k=0; k<NKPT; k++){
//         for (int m=0; m<NBAND_V[k]; m++){
//             delete [] C_k[k][m];
//         }
//         delete [] C_k[k];
//         delete [] E_k[k];
//     }
//     delete [] E_k;
//     delete [] C_k;

//     return;
// }


}/* namespace pmx */










// // Calculate lattice longitudinal and transverse permittivity matrix elements at G1, G2
// void permittivity(env &dat){

//     // use indirect function calls for diracdelta function
//     if (DELTA) {
//         diracdelta = &diracdelta_lorentzian;
//     }
//     else {
//         diracdelta = &diracdelta_gaussian;
//     }



//     // scale factors
//     double scaleL = 8.0*pow(e,2)/(pi*eV*a*angstrom);
//     // double scaleL = 2.0*dat.lat.bz_volume*pow(e,2)/(pi*eV*a*angstrom);
//     double scaleT = 32.0*dat.lat.bz_volume*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));

//     std::complex<double> ***eigvector_k; // eigenvector at k
//     double **eigvalue_k; // eigenvalue at k
//     std::complex<double> ***eigvector_kq; // eigenvector at k+q
//     double **eigvalue_kq; // eigenvalue at k

//     // overlap integral L, <k,v|e^{-i*G_m}|k+q,c> <k+q,c|e^{i*G_n}|k,v>
//     std::complex<double> *****ointL;
//     // overlap integral T, t_{Q+G_m} * <k,v|e^{-i*G_m}j_0|k+q,c> <k+q,c|e^{i*G_n}j_0|k,v> * t_{Q+G_n}
//     std::complex<double> *****ointT;

//     int NBAND_V [NKPT]; // number of valence bands
//     int NBAND_C [NKPT]; // number of conduction bands

//     dat.realepsL = new double ***[NQ];
//     dat.imagepsL = new double ***[NQ];
//     dat.realepsT = new double ***[NQ];
//     dat.imagepsT = new double ***[NQ];



//     //======================================================================
//     eigvalue_k = new double *[NKPT];
//     eigvector_k = new std::complex<double> **[NKPT];
//     #pragma omp parallel for
//     for (int k=0; k<NKPT; k++){
//         // vector{k}
//         Eigen::Vector3d K = dat.lat.K[k];

//         // Hamiltonian at k
//         Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K,dat.lat.atomic,dat.mat);   
        
//         // compute eigenvalue problem at k
//         Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);
//         if (eigsolver_k.info() != Eigen::Success){ printf("ERROR: Eigensolver failed to solve the matrix."); abort(); }

//         eigvalue_k[k] = new double [NBAND];
//         for (int n=0; n<NBAND; n++){
//             eigvalue_k[k][n] = eigsolver_k.eigenvalues()(n)-dat.energyoffset;
//             if (eigvalue_k[k][n] > 0) {
//                 NBAND_V[k] = n;
//                 break;
//             }
//         }
//         eigvector_k[k] = new std::complex<double> *[NBAND_V[k]];
//         for (int n=0; n<NBAND_V[k]; n++){
//             eigvector_k[k][n] = new std::complex<double> [NPW];
//             for (int j=0; j<NPW; j++){
//                 eigvector_k[k][n][j] = eigsolver_k.eigenvectors().col(n)(j);
//             }
//         }
//     }
//     #pragma omp barrier
    

//     //======================================================================
//     ointL = new std::complex<double> ****[NKPT];
//     ointT = new std::complex<double> ****[NKPT];
//     eigvalue_kq = new double *[NKPT];
//     eigvector_kq = new std::complex<double> **[NKPT];
//     for (int k=0; k<NKPT; k++){
//         ointL[k] = new std::complex<double>***[NEPS];
//         ointT[k] = new std::complex<double>***[NEPS];
//         for (int m=0; m<NEPS; m++){
//             ointL[k][m] = new std::complex<double> **[NEPS];
//             ointT[k][m] = new std::complex<double> **[NEPS];
//             for (int n=0; n<=m; n++){
//                 ointL[k][m][n] = new std::complex<double> *[NBAND_V[k]];
//                 ointT[k][m][n] = new std::complex<double> *[NBAND_V[k]];
//             }
//         }
//         eigvalue_kq[k] = new double [NBAND];
//     }

//     for (int q=0; q<NQ; q++){
//         Eigen::Vector3d Q = dat.lat.Q[q];
//         #pragma omp parallel for
//         for (int k=0; k<NKPT; k++){
            
//             // vector{k}
//             Eigen::Vector3d K = dat.lat.K[k];

//             // Hamiltonian at k+q
//             Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K+Q,dat.lat.atomic,dat.mat);   

//             // compute eigenvalue problem at k+q
//             Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_kq(H);
//             if(eigsolver_kq.info() != Eigen::Success){ printf("ERROR: Eigensolver failed to solve the matrix."); abort(); }
            
            
//             int tempindex;
//             for (int n=NBAND-1; n>=0; n--){
//                 tempindex = NBAND-1-n;
//                 eigvalue_kq[k][tempindex] = double(eigsolver_kq.eigenvalues()(n))-dat.energyoffset;
//                 if (eigvalue_kq[k][tempindex] <=0) {
//                     NBAND_C[k] = tempindex;
//                     break;
//                 }
//             }
//             eigvector_kq[k] = new std::complex<double> *[NBAND_C[k]];
//             for (int n=0; n<NBAND_C[k]; n++){
//                 tempindex = NBAND-1-n;
//                 eigvector_kq[k][n] = new std::complex<double> [NPW];
//                 for (int j=0; j<NPW; j++){
//                     eigvector_kq[k][n][j] = eigsolver_kq.eigenvectors().col(tempindex)(j);
//                 }
//             }

//             // overlap integrals
//             for (int m=0; m<NEPS; m++){
//                 for (int n=0; n<=m; n++){
//                     // Eigen::Vector3d vecqGm = Q+dat.G[m]; // vector q+Gm     
//                     // Eigen::Vector3d vecqGn = Q+dat.G[n]; // vector q+Gn    
//                     // Eigen::Vector3d vecqGm_t = perpvec(vecqGm); // normal vector perpendicular to q+Gm
//                     // Eigen::Vector3d vecqGn_t = perpvec(vecqGn); // normal vector perpendicular to q+Gn


//                     // vecqGm_t = vecqGm_t.cross(vecqGm);
//                     // vecqGn_t = vecqGn_t.cross(vecqGn);

//                     // Eigen::Vector3d vecqGm = Q; // vector q+Gm     
//                     // Eigen::Vector3d vecqGn = Q; // vector q+Gn    
//                     // Eigen::Vector3d vecqGm_t = perpvec(vecqGm); // normal vector perpendicular to q+Gm
//                     // Eigen::Vector3d vecqGn_t = perpvec(vecqGm); // normal vector perpendicular to q+Gn
//                     // vecqGm_t = vecqGm_t.cross(vecqGm);
//                     // vecqGn_t = vecqGn_t.cross(vecqGn);
                    
//                     //====================================================================
//                     // // LL
//                     // Eigen::Vector3d vecqGm_t = normvec(Q+dat.G[m]);
//                     // Eigen::Vector3d vecqGn_t = normvec(Q+dat.G[n]);

//                     // vecqGm_t = vecqGm_t/vecqGm_t.norm();
//                     // vecqGn_t = vecqGn_t/vecqGn_t.norm();

//                     // TT
//                     Eigen::Vector3d vecqGm_t = perpvec(Q+dat.lat.G[m]);
//                     Eigen::Vector3d vecqGn_t = perpvec(Q+dat.lat.G[n]);


//                     // Eigen::Vector3d vecqGm_t = tanvec(Q+dat.G[m],Q);
//                     // Eigen::Vector3d vecqGn_t = tanvec(Q+dat.G[n],Q);
//                     // if (m==3 || m==4 || m==7 || m==8){
//                     //     Eigen::Vector3d normQG = Q+dat.G[m];
//                     //     normQG = normQG/normQG.norm();
//                     //     vecqGm_t = vecqGm_t.cross(normQG);
//                     // }
//                     // if (n==3 || n==4 || n==7 || n==8){
//                     //     Eigen::Vector3d normQG = Q+dat.G[n];
//                     //     normQG = normQG/normQG.norm();
//                     //     vecqGn_t = vecqGn_t.cross(normQG);
//                     // }


//                     // // T'T'
//                     // Eigen::Vector3d vecqGm_t = perpvec(Q+dat.lat.G[m]);
//                     // Eigen::Vector3d vecqGn_t = perpvec(Q+dat.lat.G[n]);
//                     // Eigen::Vector3d normQG = normvec(Q+dat.lat.G[m]);
//                     // vecqGm_t = vecqGm_t.cross(normQG);
//                     // normQG = normvec(Q+dat.lat.G[n]);
//                     // vecqGn_t = vecqGn_t.cross(normQG);

//                     // // TT'
//                     // Eigen::Vector3d vecqGm_t = perpvec(Q+dat.G[m]);
//                     // Eigen::Vector3d vecqGn_t = perpvec(Q+dat.G[n]);
//                     // Eigen::Vector3d normQG = Q+dat.G[n];
//                     // normQG = normQG/normQG.norm();
//                     // vecqGn_t = vecqGn_t.cross(normQG);


//                     // // LT
//                     // Eigen::Vector3d vecqGm_t2 = Q+dat.G[m];
//                     // vecqGm_t2 = vecqGm_t2/vecqGm_t.norm();
//                     // Eigen::Vector3d vecqGn_t2 = perpvec(Q+dat.G[n]);







//                     // non-G ==================================================

//                     // LL
//                     // Eigen::Vector3d vecqGm_t = Q/Q.norm();
//                     // Eigen::Vector3d vecqGn_t = Q/Q.norm();

//                     // LT
//                     // Eigen::Vector3d vecqGm_t = Q/Q.norm();
//                     // Eigen::Vector3d vecqGn_t = perpvec(Q);

//                     // // TT
//                     // Eigen::Vector3d vecqGm_t = perpvec(Q);
//                     // Eigen::Vector3d vecqGn_t = perpvec(Q);

//                     // TT'
//                     // Eigen::Vector3d vecqGm_t = perpvec(Q);
//                     // Eigen::Vector3d vecqGn_t = perpvec(Q);
//                     // vecqGn_t = vecqGn_t.cross(Q)/Q.norm();


//                     for (int v=0; v<NBAND_V[k]; v++){ // valence bands
//                         ointL[k][m][n][v] = new std::complex<double> [NBAND_C[k]];
//                         ointT[k][m][n][v] = new std::complex<double> [NBAND_C[k]];
//                         for (int c=0; c<NBAND_C[k]; c++){ // conduction bands

//                             std::complex<double> ointL_k_kq = 0.0; // <k,v|e^{-i*G_m}|k+q,c>
//                             std::complex<double> ointT_k_kq = 0.0; // <k,v|e^{-i*G_m}j_0|k+q,c> * t_{G_m+Q}
//                             std::complex<double> ointL_kq_k = 0.0; // <k+q,c|e^{i*G_n}|k,v>
//                             std::complex<double> ointT_kq_k = 0.0; // <k+q,c|e^{i*G_n}j_0|k,v> * t_{G_n+Q}

//                             for (int i=0; i<NPW; i++){
//                                 if (dat.lat.loci[m][i]!=-1){
//                                     ointL_k_kq += conj(eigvector_k[k][v][i]) * eigvector_kq[k][c][dat.lat.loci[m][i]];
//                                     ointT_k_kq += conj(eigvector_k[k][v][i]) * eigvector_kq[k][c][dat.lat.loci[m][i]] 
//                                                                 * (dat.lat.G[i]+K +dat.lat.G[m]/2+Q/2).dot(vecqGm_t); // 
//                                 }
//                                 if (dat.lat.loci[n][i]!=-1){
//                                     ointL_kq_k += conj(eigvector_kq[k][c][dat.lat.loci[n][i]]) * eigvector_k[k][v][i];
//                                     ointT_kq_k += conj(eigvector_kq[k][c][dat.lat.loci[n][i]]) * eigvector_k[k][v][i] 
//                                                                 * (dat.lat.G[i]+K +dat.lat.G[n]/2+Q/2).dot(vecqGn_t); //
//                                 }
//                             }

//                             ointL[k][m][n][v][c] = ointL_k_kq * ointL_kq_k;
//                             ointT[k][m][n][v][c] = ointT_k_kq * ointT_kq_k;
//                         }
//                     }
//                 }
//             }
//         }
//         #pragma omp barrier


//         dat.realepsL[q] = new double **[NEPS];
//         dat.imagepsL[q] = new double **[NEPS];
//         dat.realepsT[q] = new double **[NEPS];
//         dat.imagepsT[q] = new double **[NEPS];
//         for (int m=0; m<NEPS; m++){
//             dat.realepsL[q][m] = new double *[NEPS];
//             dat.imagepsL[q][m] = new double *[NEPS];
//             dat.realepsT[q][m] = new double *[NEPS];
//             dat.imagepsT[q][m] = new double *[NEPS];
//             for (int n=0; n<=m; n++){
//                 dat.realepsL[q][m][n] = new double [NFREQ];
//                 dat.imagepsL[q][m][n] = new double [NFREQ];
//                 dat.realepsT[q][m][n] = new double [NFREQ];
//                 dat.imagepsT[q][m][n] = new double [NFREQ];
//                 #pragma omp parallel for //collapse(3)
//                 for (int f=0; f<NFREQ; f++){
//                     // temporary variables
//                     double realL1 = 0.0e0;
//                     double imagL1 = 0.0e0;
//                     double realT1 = 0.0e0;
//                     double imagT1 = 0.0e0;
//                     double depsilon = 0;
//                     double realdL, imagdL, realdT, imagdT;

//                     // sum over k-grids & bands
//                     for (int k=0; k<NKPT; k++){
//                         for (int v=0; v<NBAND_V[k]; v++){ for (int c=0; c<NBAND_C[k]; c++){
//                             depsilon = eigvalue_kq[k][c]-eigvalue_k[k][v];
//                             imagL1 += dat.lat.KW[k] * (ointL[k][m][n][v][c]).real() 
//                                                     * (*diracdelta)(depsilon-dat.freq[f]);
//                             imagT1 += dat.lat.KW[k] * (ointT[k][m][n][v][c]).real() 
//                                                     * (*diracdelta)(depsilon-dat.freq[f]);
//                             if (dat.kk==0){
//                                 realL1 += dat.lat.KW[k] * (ointL[k][m][n][v][c]).real() 
//                                                         * (depsilon) / (depsilon*depsilon-dat.freq[f]*dat.freq[f]);
//                                 realT1 += dat.lat.KW[k] * (ointT[k][m][n][v][c]).real() 
//                                                         / (depsilon) / (depsilon*depsilon-dat.freq[f]*dat.freq[f]);
//                             }
//                         }}
//                     }

//                     // longitudinal permittivity
//                     Eigen::Vector3d vecqGm = Q+dat.lat.G[m];
//                     Eigen::Vector3d vecqGn = Q+dat.lat.G[n];
//                     imagdL = 0.0e0;
//                     if (dat.freq[f]>0) imagdL = scaleL * pi * imagL1 / (vecqGm.norm()*vecqGn.norm());
//                     dat.imagepsL[q][m][n][f] = imagdL;
                    
//                     // transverse permittivity
//                     imagdT = 0.0e0;
//                     if (dat.freq[f]>0) imagdT = scaleT * pi * imagT1 / (dat.freq[f]*dat.freq[f]);
//                     dat.imagepsT[q][m][n][f] = imagdT;

//                     if (dat.kk==0){
//                         realdL = scaleL * 2 * realL1 / (vecqGm.norm()*vecqGn.norm());
//                         if (m==n) realdL += 1.0;
//                         dat.realepsL[q][m][n][f] = realdL;
//                         realdT = scaleT * 2 * realT1;
//                         if (m==n) realdT += 1.0;
//                         dat.realepsT[q][m][n][f] = realdT;
//                     }
//                 }
//             }
//         }
//         #pragma omp barrier

//         // kramer-kronig transform
//         if (dat.kk){
//             for (int m=0; m<NEPS; m++){
//                 for (int n=0; n<=m; n++){
//                     kramerskronigtransform(dat.realepsL[q][m][n],dat.imagepsL[q][m][n],dat.freq,dat.dfreq);
//                     kramerskronigtransform(dat.realepsT[q][m][n],dat.imagepsT[q][m][n],dat.freq,dat.dfreq);
//                     if (m==n){
//                         for (int f=0; f<NFREQ; f++){
//                             dat.realepsL[q][m][n][f] += 1;
//                             dat.realepsT[q][m][n][f] += 1;
//                         }
//                     }
//                 }
//             }
//         }

//         // free dynamically allocated memory
//         for (int k=0; k<NKPT; k++){
//             for (int m=0; m<NEPS; m++){
//                 for (int n=0; n<=m; n++){
//                     for (int v=0; v<NBAND_V[k]; v++){
//                         delete [] ointL[k][m][n][v];
//                         delete [] ointT[k][m][n][v];
//                     }
//                 }
//             }
//             for (int m=0; m<NBAND_C[k]; m++){
//                 delete [] eigvector_kq[k][m];
//             }
//             delete [] eigvector_kq[k];
//         }
//     }// loop q

//     // free dynamically allocated memory
//     for (int k=0; k<NKPT; k++){
//         for (int m=0; m<NEPS; m++){
//             for (int n=0; n<=m; n++){
//                 delete [] ointL[k][m][n];
//                 delete [] ointT[k][m][n];
//             }
//             delete [] ointL[k][m];
//             delete [] ointT[k][m];
//         }
//         delete [] ointL[k];
//         delete [] ointT[k];

//         for (int m=0; m<NBAND_V[k]; m++){
//             delete [] eigvector_k[k][m];
//         }
//         delete [] eigvector_k[k];
//         delete [] eigvalue_k[k];
//         delete [] eigvalue_kq[k];
//     }
//     delete [] ointL;
//     delete [] ointT;
//     delete [] eigvalue_k;
//     delete [] eigvalue_kq;
//     delete [] eigvector_k;
//     delete [] eigvector_kq;



//     NTSR = 1;
//     // initialize susceptibility tensor matrix elements
//     dat.ReXij = new double *****[NQ];
//     dat.ImXij = new double *****[NQ];
//     for (int q=0; q<NQ; q++){
//         dat.ReXij[q] = new double ****[NTSR];
//         dat.ImXij[q] = new double ****[NTSR];
//         for (int i=0; i<NTSR; i++){
//             dat.ReXij[q][i] = new double ***[NTSR];
//             dat.ImXij[q][i] = new double ***[NTSR];
//             for (int j=0; j<NTSR; j++){
//                 dat.ReXij[q][i][j] = new double **[NEPS];
//                 dat.ImXij[q][i][j] = new double **[NEPS];
//                 for (int m=0; m<NEPS; m++){
//                     dat.ReXij[q][i][j][m] = new double *[NEPS];
//                     dat.ImXij[q][i][j][m] = new double *[NEPS];
//                     for (int n=0; n<=m; n++){
//                         dat.ReXij[q][i][j][m][n] = new double [NFREQ];
//                         dat.ImXij[q][i][j][m][n] = new double [NFREQ];
//                         for (int f=0; f<NFREQ; f++){
//                             dat.ReXij[q][i][j][m][n][f] = dat.realepsL[q][m][n][f];
//                             dat.ImXij[q][i][j][m][n][f] = dat.imagepsL[q][m][n][f];
//                         }
//                     }
//                 }
//             }
//         }
//     }





//     return;
// }

// // longitudinal-longitudinal part of susceptibility tensor matrix
// void chi_LL_old(env &dat){
//     std::chrono::time_point<std::chrono::system_clock> time_0, time_1;
//     std::chrono::duration<double> elapsed_time;

//     NTSR = 3;

//     // use indirect function calls for diracdelta function
//     switch (DELTA) {
//         case 0:
//             diracdelta = &diracdelta_gaussian;
//             break;
//         case 1:
//             diracdelta = &diracdelta_lorentzian;
//             break;
//         case 2:
//             diracdelta = &diracdelta_triangle;
//             break;
//         case 3:
//             diracdelta = &diracdelta_rect;
//             break;
//     }

//     // scalefactor
//     // double SCALEFACTOR = 2.0*dat.lat.bz_volume*pow(e,2)/(pi*eV*a*angstrom);
//     double SCALEFACTOR = 8.0*pow(e,2)/(pi*eV*a*angstrom);

//     // initialize susceptibility tensor matrix elements
//     dat.ReXij = new double *****[NQ];
//     dat.ImXij = new double *****[NQ];
//     for (int q=0; q<NQ; q++){
//         dat.ReXij[q] = new double ****[NTSR];
//         dat.ImXij[q] = new double ****[NTSR];
//         for (int i=0; i<NTSR; i++){
//             dat.ReXij[q][i] = new double ***[NTSR];
//             dat.ImXij[q][i] = new double ***[NTSR];
//             for (int j=0; j<NTSR; j++){
//                 dat.ReXij[q][i][j] = new double **[NEPS];
//                 dat.ImXij[q][i][j] = new double **[NEPS];
//                 for (int m=0; m<NEPS; m++){
//                     dat.ReXij[q][i][j][m] = new double *[NEPS];
//                     dat.ImXij[q][i][j][m] = new double *[NEPS];
//                     for (int n=0; n<=m; n++){
//                         dat.ReXij[q][i][j][m][n] = new double [NFREQ];
//                         dat.ImXij[q][i][j][m][n] = new double [NFREQ];
//                     }
//                 }
//             }
//         }
//     }

//     //======================================================================
//     std::cout << "  Solving for |k,v>..." << std::endl;
//     time_0 = std::chrono::system_clock::now();
//     int NBAND_C [NKPT]; // number of conduction bands
//     double **E_k; // eigenvalue at k [NK x NC]
//     E_k = new double *[NKPT];
//     std::complex<double> ***C_k; // eigenvector at k [NK x NC]
//     C_k = new std::complex<double> **[NKPT];
//     #pragma omp parallel for
//     for (int k=0; k<NKPT; k++){
//         // vector{k}
//         Eigen::Vector3d K = dat.lat.K[k];

//         // Hamiltonian at k
//         Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K,dat.lat.atomic,dat.mat);   
        
//         // compute eigenvalue problem at k
//         Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);

//         E_k[k] = new double [NBAND];
//         for (int c=0; c<NBAND; c++){
//             E_k[k][c] = eigsolver_k.eigenvalues()(NBAND-1-c)-dat.energyoffset;
//             if (E_k[k][c]<=0) {
//                 NBAND_C[k] = c;
//                 break;
//             }
//         }
//         C_k[k] = new std::complex<double> *[NBAND_C[k]];
//         for (int c=0; c<NBAND_C[k]; c++){
//             C_k[k][c] = new std::complex<double> [NPW];
//             for (int p=0; p<NPW; p++){
//                 C_k[k][c][p] = eigsolver_k.eigenvectors().col(NBAND-1-c)(p);
//             }// loop over p
//         }// loop over c
//     }
//     #pragma omp barrier
//     time_1 = std::chrono::system_clock::now();
//     elapsed_time = time_1-time_0;
//     std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;

//     //======================================================================
//     // initialize large temporary variables
//     // overlap integral, t_{Q+G_m} * <k,v|e^{-i*G_m}j_0|k+q,c> <k+q,c|e^{i*G_n}j_0|k,v> * t_{Q+G_n}
//     std::complex<double> *****oint;
//     oint = new std::complex<double> ****[NKPT]; // overlap integral [NK x NTSR x NEPS x NC x NV]
//     for (int k=0; k<NKPT; k++){
//         oint[k] = new std::complex<double> ***[3];
//         for (int i=0; i<NTSR; i++){
//             oint[k][i] = new std::complex<double> **[NEPS];
//             for (int m=0; m<NEPS; m++){
//                 oint[k][i][m] = new std::complex<double> *[NBAND_C[k]];
//             }
//         }
//     }
//     // std::complex<double> im = std::complex<double>(0.0,1.0);


//     std::cout << "  Solving for susceptibility tensor matrix elements (LL only)..." << std::endl;
//     int NBAND_V [NKPT]; // number of valence bands
//     std::complex<double> ***C_kq; // eigenvector at k+q
//     double **E_kq; // eigenvalue at k
//     E_kq = new double *[NKPT];
//     C_kq = new std::complex<double> **[NKPT];
//     for (int k=0; k<NKPT; k++){
//         E_kq[k] = new double [NBAND];
//     }
//     // #pragma omp parallel
//     for (int q=0; q<NQ; q++){
//         std::cout << "  (" << (q+1) << "/" << NQ << ")..." << std::endl; // (q/NQ)
//         Eigen::Vector3d Q = dat.lat.Q[q];

//         std::cout << "    Solving for overlap integrals..." << std::endl;
//         time_0 = std::chrono::system_clock::now();
//         #pragma omp parallel for
//         for (int k=0; k<NKPT; k++){
            
//             // Hamiltonian at k+q
//             Eigen::Vector3d K = dat.lat.K[k];
//             Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K+Q,dat.lat.atomic,dat.mat);   

//             // compute eigenvalue problem at k+q
//             Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_kq(H);
//             for (int v=0; v<NBAND; v++){
//                 E_kq[k][v] = double(eigsolver_kq.eigenvalues()(v))-dat.energyoffset;
//                 if (E_kq[k][v]>0) {
//                     NBAND_V[k] = v;
//                     break;
//                 }
//             }
//             C_kq[k] = new std::complex<double> *[NBAND_V[k]];
//             for (int v=0; v<NBAND_V[k]; v++){
//                 C_kq[k][v] = new std::complex<double> [NPW];
//                 for (int p=0; p<NPW; p++){
//                     C_kq[k][v][p] = eigsolver_kq.eigenvectors().col(v)(p);
//                 }
//             }
            
//             // overlap integrals
//             for (int m=0; m<NEPS; m++){
//                 // std::vector<Eigen::Vector3cd> uvec_m = unitvecq(Q+dat.lat.G[m],Q); //Q+
//                 // std::vector<Eigen::Vector3cd> uvec_m = unitvec(Q+dat.lat.G[m]);
//                 for (int i=0; i<NTSR; i++){
//                     for (int c=0; c<NBAND_C[k]; c++){
//                         oint[k][i][m][c] = new std::complex<double> [NBAND_V[k]];
//                         for (int v=0; v<NBAND_V[k]; v++){
//                             // u^i_{q+g_m} * <k,c|e^{-i*(q+g_m)*r}|k+q,v>
//                             oint[k][i][m][c][v] = 0;
//                             for (int p=0; p<NPW; p++){if (dat.lat.loci[m][p]!=-1){
//                                 oint[k][i][m][c][v] += (conj(C_k[k][c][p]) * C_kq[k][v][dat.lat.loci[m][p]]);
//                             }}
//                         }//loop over v
//                     }//loop over c
//                 }//loop over i
//             }//loop over m
//         }//loop over k
//         #pragma omp barrier
//         time_1 = std::chrono::system_clock::now();
//         elapsed_time = time_1-time_0;
//         std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;

//         std::cout << "    Solving for Xijmn..." << std::endl;
//         time_0 = std::chrono::system_clock::now();
//         #pragma omp parallel for collapse(4)
//         for (int i=0; i<NTSR; i++){for (int j=0; j<NTSR; j++){for (int m=0; m<NEPS; m++){for (int n=0; n<=m; n++){
//             Eigen::Vector3d vecqgm = Q+dat.lat.G[m];
//             Eigen::Vector3d vecqgn = Q+dat.lat.G[n];
//             for (int f=0; f<NFREQ; f++){
//                 double tmp_real_1 = 0;
//                 double tmp_imag_1 = 0;
//                 double tmp_real_2 = 0;
//                 double tmp_imag_2 = 0;
//                 // sum over k,c,v
//                 for (int k=0; k<NKPT; k++){for (int c=0; c<NBAND_C[k]; c++){for (int v=0; v<NBAND_V[k]; v++){
//                     double dE = E_k[k][c]-E_kq[k][v];
//                     std::complex<double> Oij = oint[k][i][m][c][v] * conj(oint[k][j][n][c][v]);
                    
//                     // 
//                     tmp_imag_1 += dat.lat.KW[k] * Oij.real()
//                                           * (*diracdelta)(dE-dat.freq[f]);
//                     tmp_real_1 -= dat.lat.KW[k] * Oij.imag()
//                                           * (*diracdelta)(dE-dat.freq[f]);

//                     if (dat.kk==0){
//                         tmp_real_2 += dat.lat.KW[k] * Oij.real()
//                                          * (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
//                         tmp_imag_2 += dat.lat.KW[k] * Oij.imag()
//                                          * (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
//                     }
//                 }}}
//                 if (f==0){
//                     tmp_imag_1 = 0;
//                     tmp_real_1 = 0;
//                 }else{
//                     tmp_imag_1 /= (vecqgm.norm()*vecqgn.norm());
//                     // tmp_real_1 /= (vecqgm.norm()*vecqgn.norm());
//                 }
//                 if (dat.kk==0){
//                     // tmp_imag_2 /= (vecqgm.norm()*vecqgn.norm());
//                     tmp_real_2 /= (vecqgm.norm()*vecqgn.norm());
//                 }
//                 dat.ImXij[q][i][j][m][n][f] = SCALEFACTOR * (pi*tmp_imag_1+2.0*tmp_imag_2);
//                 dat.ReXij[q][i][j][m][n][f] = SCALEFACTOR * (pi*tmp_real_1+2.0*tmp_real_2);
//             }
//         }}}}
//         #pragma omp barrier
//         time_1 = std::chrono::system_clock::now();
//         elapsed_time = time_1-time_0;
//         std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;


        
//         // kramer-kronig transform
//         if(dat.kk){
//             std::cout << "    Kramers-Kronig transform..." << std::endl;
//             time_0 = std::chrono::system_clock::now();
//             #pragma omp parallel for collapse(4)
//             for (int i=0; i<NTSR; i++){ for (int j=0; j<NTSR; j++){ for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
//                 kramerskronigtransform(dat.ReXij[q][i][j][m][n],dat.ImXij[q][i][j][m][n],dat.freq,dat.dfreq);
//             }}}}
//             #pragma omp barrier
//             time_1 = std::chrono::system_clock::now();
//             elapsed_time = time_1-time_0;
//             std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;
//         }
        

//         // free dynamically allocated memory
//         for (int k=0; k<NKPT; k++){for (int i=0; i<NTSR; i++){for (int m=0; m<NEPS; m++){for (int c=0; c<NBAND_C[k]; c++){
//             delete [] oint[k][i][m][c];
//         }}}}
//         for (int k=0; k<NKPT; k++){
//             for (int v=0; v<NBAND_V[k]; v++){
//                 delete [] C_kq[k][v];
//             }
//             delete [] C_kq[k];
//         }
//     }// loop q

//     // free dynamically allocated memory
//     for (int k=0; k<NKPT; k++){
//         for (int i=0; i<NTSR; i++){
//             for (int m=0; m<NEPS; m++){
//                 delete [] oint[k][i][m];
//             }
//             delete [] oint[k][i];
//         }
//         delete [] oint[k];
//     }
//     delete [] oint;

//     for (int k=0; k<NKPT; k++){
//         for (int c=0; c<NBAND_C[k]; c++){
//             delete [] C_k[k][c];
//         }
//         delete [] C_k[k];
//         delete [] E_k[k];
//         delete [] E_kq[k];
//     }
//     delete [] E_k;
//     delete [] E_kq;
//     delete [] C_k;
//     delete [] C_kq;

//     return;
// }

