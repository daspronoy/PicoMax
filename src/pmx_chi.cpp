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
