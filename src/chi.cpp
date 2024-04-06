/*


*/
#include "pmx.hpp"

namespace pmx{

inline double (*diracdelta) (double value);

inline double diracdelta_lorentzian(double x){
    return EPSILON/pi/(pow(EPSILON,2)+pow(x,2));
}

inline double diracdelta_gaussian(double x){
    return exp(-x*x/(EPSILON*EPSILON))/(EPSILON*sqrt(pi));
}

// kramers-kronig transformation
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


// Calculate lattice longitudinal and transverse permittivity matrix elements at G1, G2
void permittivity(env &dat){

    // use indirect function calls for diracdelta function
    if (DELTA) {
        diracdelta = &diracdelta_lorentzian;
    }
    else {
        diracdelta = &diracdelta_gaussian;
    }



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

    dat.realepsL = new double ***[NQ];
    dat.imagepsL = new double ***[NQ];
    dat.realepsT = new double ***[NQ];
    dat.imagepsT = new double ***[NQ];



    //======================================================================
    eigvalue_k = new double *[NKPT];
    eigvector_k = new std::complex<double> **[NKPT];
    #pragma omp parallel for
    for (int k=0; k<NKPT; k++){
        // vector{k}
        Eigen::Vector3d K = {dat.K[k][0], dat.K[k][1], dat.K[k][2]};

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

    for (int q=0; q<NQ; q++){
        Eigen::Vector3d Q = dat.Q[q];
        #pragma omp parallel for
        for (int k=0; k<NKPT; k++){
            
            // vector{k}
            Eigen::Vector3d K = {dat.K[k][0], dat.K[k][1], dat.K[k][2]};

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
                    // Eigen::Vector3d vecqGm_t = normvec(Q+dat.G[m]);
                    // Eigen::Vector3d vecqGn_t = normvec(Q+dat.G[n]);

                    // vecqGm_t = vecqGm_t/vecqGm_t.norm();
                    // vecqGn_t = vecqGn_t/vecqGn_t.norm();

                    // TT
                    // Eigen::Vector3d vecqGm_t = perpvec(Q+dat.G[m]);
                    // Eigen::Vector3d vecqGn_t = perpvec(Q+dat.G[n]);


                    // Eigen::Vector3d vecqGm_t = tanvec(Q+dat.G[m],Q);
                    // Eigen::Vector3d vecqGn_t = tanvec(Q+dat.G[n],Q);
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
                    Eigen::Vector3d vecqGm_t = perpvec(Q+dat.G[m]);
                    Eigen::Vector3d vecqGn_t = perpvec(Q+dat.G[n]);
                    Eigen::Vector3d normQG = normvec(Q+dat.G[m]);
                    vecqGm_t = vecqGm_t.cross(normQG);
                    normQG = normvec(Q+dat.G[n]);
                    vecqGn_t = vecqGn_t.cross(normQG);

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
                                    ointT_k_kq += conj(eigvector_k[k][v][i]) * eigvector_kq[k][c][dat.loci[m][i]] 
                                                                * (dat.G[i]+K +dat.G[m]/2+Q/2).dot(vecqGm_t); // 
                                }
                                if (dat.loci[n][i]!=-1){
                                    ointL_kq_k += conj(eigvector_kq[k][c][dat.loci[n][i]]) * eigvector_k[k][v][i];
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
                    double realL1 = 0.0e0;
                    double imagL1 = 0.0e0;
                    double realT1 = 0.0e0;
                    double imagT1 = 0.0e0;
                    double depsilon = 0;
                    double realdL, imagdL, realdT, imagdT;

                    // sum over k-grids & bands
                    for (int k=0; k<NKPT; k++){
                        for (int v=0; v<NBAND_V[k]; v++){ for (int c=0; c<NBAND_C[k]; c++){
                            depsilon = eigvalue_kq[k][c]-eigvalue_k[k][v];
                            imagL1 += dat.KW[k] * (ointL[k][m][n][v][c]).real() 
                                                    * (*diracdelta)(depsilon-dat.freq[f]);
                            imagT1 += dat.KW[k] * (ointT[k][m][n][v][c]).real() 
                                                    * (*diracdelta)(depsilon-dat.freq[f]);
                            if (dat.kk==0){
                                realL1 += dat.KW[k] * (ointL[k][m][n][v][c]).real() 
                                                        * (depsilon) / (depsilon*depsilon-dat.freq[f]*dat.freq[f]);
                                realT1 += dat.KW[k] * (ointT[k][m][n][v][c]).real() 
                                                        / (depsilon) / (depsilon*depsilon-dat.freq[f]*dat.freq[f]);
                            }
                        }}
                    }

                    // longitudinal permittivity
                    Eigen::Vector3d vecqGm = Q+dat.G[m];
                    Eigen::Vector3d vecqGn = Q+dat.G[n];
                    imagdL = 0.0e0;
                    if (dat.freq[f]>0) imagdL = scaleL * pi * imagL1 / (vecqGm.norm()*vecqGn.norm());
                    dat.imagepsL[q][m][n][f] = imagdL;
                    
                    // transverse permittivity
                    imagdT = 0.0e0;
                    if (dat.freq[f]>0) imagdT = scaleT * pi * imagT1 / (dat.freq[f]*dat.freq[f]);
                    dat.imagepsT[q][m][n][f] = imagdT;

                    if (dat.kk==0){
                        realdL = scaleL * 2 * realL1 / (vecqGm.norm()*vecqGn.norm());
                        if (m==n) realdL += 1.0;
                        dat.realepsL[q][m][n][f] = realdL;
                        realdT = scaleT * 2 * realT1;
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





// full susceptibility tensor

void chi_tensor(env &dat){

    // use indirect function calls for diracdelta function
    if (DELTA) {
        diracdelta = &diracdelta_lorentzian;
    }
    else {
        diracdelta = &diracdelta_gaussian;
    }

    // scalefactor
    double SCALEFACTOR = 128.0*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));

    // initialize dielectric tensor matrix elements
    dat.realepsij = new double *****[3];
    dat.imagepsij = new double *****[3];
    for (int i=0; i<3; i++){
        dat.realepsij[i] = new double ****[3];
        dat.imagepsij[i] = new double ****[3];
        for (int j=0; j<3; j++){
            dat.realepsij[i][j] = new double ***[NQ];
            dat.imagepsij[i][j] = new double ***[NQ];
            for (int q=0; q<NQ; q++){
                dat.realepsij[i][j][q] = new double **[NEPS];
                dat.imagepsij[i][j][q] = new double **[NEPS];
                for (int m=0; m<NEPS; m++){
                    dat.realepsij[i][j][q][m] = new double *[NEPS];
                    dat.imagepsij[i][j][q][m] = new double *[NEPS];
                    for (int n=0; n<=m; n++){
                        dat.realepsij[i][j][q][m][n] = new double [NFREQ];
                        dat.imagepsij[i][j][q][m][n] = new double [NFREQ];
                    }
                }
            }
        }
    }

    //======================================================================
    int NBAND_C [NKPT]; // number of conduction bands
    double **E_k; // eigenvalue at k [NK x NC]
    E_k = new double *[NKPT];
    std::complex<double> ***C_k; // eigenvector at k [NK x NC]
    C_k = new std::complex<double> **[NKPT];
    #pragma omp parallel for
    for (int k=0; k<NKPT; k++){
        // vector{k}
        Eigen::Vector3d K = {dat.K[k][0], dat.K[k][1], dat.K[k][2]};

        // Hamiltonian at k
        Eigen::MatrixXcd H = HamiltonianEPM(dat.G,K,dat.epm);   
        
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
    

    //======================================================================
    // initialize large temporary variables
    // overlap integral, t_{Q+G_m} * <k,v|e^{-i*G_m}j_0|k+q,c> <k+q,c|e^{i*G_n}j_0|k,v> * t_{Q+G_n}
    Eigen::Vector3cd ****oint;
    oint = new Eigen::Vector3cd ***[NKPT]; // overlap integral [NK x NEPS x NC x NV]
    for (int k=0; k<NKPT; k++){
        oint[k] = new Eigen::Vector3cd **[NEPS];
        for (int m=0; m<NEPS; m++){
            oint[k][m] = new Eigen::Vector3cd *[NBAND_C[k]];
        }
    }
    std::complex<double> im = std::complex<double>(0.0,1.0);

    int NBAND_V [NKPT]; // number of valence bands
    std::complex<double> ***C_kq; // eigenvector at k+q
    double **E_kq; // eigenvalue at k
    E_kq = new double *[NKPT];
    C_kq = new std::complex<double> **[NKPT];
    for (int k=0; k<NKPT; k++){
        E_kq[k] = new double [NBAND];
    }
    for (int q=0; q<NQ; q++){
        Eigen::Vector3d Q = dat.Q[q];
        #pragma omp parallel for
        for (int k=0; k<NKPT; k++){
            
            // vector{k}
            Eigen::Vector3d K = {dat.K[k][0], dat.K[k][1], dat.K[k][2]};

            // Hamiltonian at k+q
            Eigen::MatrixXcd H = HamiltonianEPM(dat.G,K+Q,dat.epm);   

            // compute eigenvalue problem at k+q
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_kq(H);
            if(eigsolver_kq.info() != Eigen::Success){ printf("ERROR: Eigensolver failed to solve the matrix."); abort(); }
            
            
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
            
            // polarization unit vector
            // Eigen::Vector3cd **uvec_m = uvecqm(Q,dat.G);
            // Eigen::Vector3d vecd = {1/8,1/8,1/8};
            // Eigen::Vector3cd oint_k_kq, oint_kq_k;
            // overlap integrals
            for (int m=0; m<NEPS; m++){
                // std::vector<Eigen::Vector3cd> uvec_m = unitvecq(Q+dat.G[m],Q);
                // std::vector<Eigen::Vector3cd> uvec_m = unitvec(Q+dat.G[m]);
                // std::vector<Eigen::Vector3cd> uvec_n = unitvec(Q+dat.G[n]);
                for (int c=0; c<NBAND_C[k]; c++){
                    oint[k][m][c] = new Eigen::Vector3cd [NBAND_V[k]];
                    for (int v=0; v<NBAND_V[k]; v++){
                        
                        // <k+q,v|e^{i*(q+g_m)*r}j_0|k,c>
                        oint[k][m][c][v] = {0,0,0};
                        for (int p=0; p<NPW; p++){
                            if (dat.loci[m][p]!=-1){
                                oint[k][m][c][v] += (conj(C_kq[k][v][dat.loci[m][p]]) * C_k[k][c][p]) 
                                            * (dat.G[p]+K+Q/2+dat.G[m]/2);
                            }
                        }
                        

                        // //---
                        // // u^i_{q+g_m} * <k+q,v|e^{i*(q+g_m)*r}\op{j}|k,c> 
                        // // <k,c|e^{-i*(q+g_n)*r}\op{j}|k+q,v> * u^j_{q+g_n} 

                        // std::complex<double> oint_kq_k = 0.0; // <k+q,c|e^{i*G_n}j_0|k,v> * t_{G_n+Q}
                        // // Eigen::Vector3cd oint_kq_k = {0,0,0};
                        // for (int p=0; p<NPW; p++){
                        //     if (dat.loci[m][p]!=-1){
                        //         oint_kq_k += (conj(C_kq[k][v][dat.loci[m][p]]) * C_k[k][c][p]) 
                        //                     * uvec_m[i].dot(dat.G[p]+K+Q/2+dat.G[m]/2); // .dot(uvec_m[n][j]) * exp(2.0*pi*im*(Q+dat.G[n]).dot(vecd))
                        //     }
                        // }
                        // std::complex<double> oint_k_kq = 0.0; // <k,v|e^{-i*G_m}j_0|k+q,c> * t_{G_m+Q}
                        // // Eigen::Vector3cd oint_k_kq = {0,0,0};
                        // for (int p=0; p<NPW; p++){
                        //     if (dat.loci[n][p]!=-1){
                        //         oint_k_kq += (conj(C_k[k][c][p]) * C_kq[k][v][dat.loci[n][p]]) 
                        //                     * (dat.G[p]+K+Q/2+dat.G[n]/2).dot(uvec_n[j]); // uvec_m[m][i] * exp(-2.0*pi*im*(Q+dat.G[m]).dot(vecd))
                        //     }
                        // }
                        
                        // // Eigen::Matrix3cd Ointkk = oint_kq_k * oint_k_kq.transpose();
                        // // oint[i][j][k][m][n][v][c] = uvec_m[m][i].adjoint() * Ointkk.real() * uvec_m[n][j];
                        // // oint[i][j][k][m][n][c][v] = uvec_m[i].adjoint() * Ointkk * uvec_n[j];
                        // oint[i][j][k][m][n][c][v] = oint_kq_k * oint_k_kq;
                    }
                }
            }
        }
        #pragma omp barrier


        
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                std::vector<Eigen::Vector3cd> uvec_m = unitvecq(Q+dat.G[m],Q);
                std::vector<Eigen::Vector3cd> uvec_n = unitvecq(Q+dat.G[n],Q);

                #pragma omp parallel for //collapse(3)
                for (int f=0; f<NFREQ; f++){
                    // temporary variables
                    double tmpval_real[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
                    double tmpval_imag[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
                    std::complex<double> tmpval_complex[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
                    std::complex<double> Oij[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
                    double depsilon = 0;
                    double realdL, imagdL, realdT, imagdT;

                    // sum over k-grids & bands
                    for (int k=0; k<NKPT; k++){
                        for (int c=0; c<NBAND_C[k]; c++){ for (int v=0; v<NBAND_V[k]; v++){
                            Eigen::Matrix3cd Omn = oint[k][m][c][v] * oint[k][n][c][v].adjoint();
                            for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                                Oij[i][j] = uvec_m[i].adjoint() * Omn * uvec_n[j];
                            }}

                            depsilon = E_k[k][c]-E_kq[k][v];
                            // imaginary part
                            for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                                tmpval_imag[i][j] += dat.KW[k] * (Oij[i][j].real())
                                                    * (*diracdelta)(depsilon-dat.freq[f]);
                            }}
                            if (dat.kk==0){
                                for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                                    tmpval_real[i][j] += dat.KW[k] * (Oij[i][j].real())
                                                    / (depsilon) / (depsilon*depsilon-dat.freq[f]*dat.freq[f]);
                                }}
                            }

                            //
                            // for (int i=0; i<3; i++){ for (int j=0; j<3; j++){
                            //     std::complex<double> Omnij = uvec_m[i].adjoint() * Omn * uvec_n[j];
                            //     tmpval_complex[i][j] -= dat.KW[k] * Omnij 
                            //                         *( 1.0/(dat.freq[f]+im*EPSILON-depsilon) - 1.0/(depsilon+dat.freq[f]+im*EPSILON)) 
                            //                         / pow(depsilon,2);
                            // }}
                        }}
                    }


                    // permittivity
                    for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                        if (dat.freq[f]>0){
                            dat.imagepsij[i][j][q][m][n][f] = SCALEFACTOR * pi * tmpval_imag[i][j] / (dat.freq[f]*dat.freq[f]);
                        }else{
                            dat.imagepsij[i][j][q][m][n][f] = 0;
                        }
                        // dat.imagepsij[i][j][q][m][n][f] = SCALEFACTOR * tmpval_complex[i][j].imag();
                    }}
                    if (dat.kk==0){
                        for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                            dat.realepsij[i][j][q][m][n][f] = SCALEFACTOR * 2 * tmpval_real[i][j];
                            // dat.realepsij[i][j][q][m][n][f] = SCALEFACTOR * tmpval_complex[i][j].real();
                        }}
                    }
                }
            }
        }
        #pragma omp barrier

        // kramer-kronig transform
        if (dat.kk){
            for (int i=0; i<3; i++){
                for (int j=0; j<3; j++){
                    for (int m=0; m<NEPS; m++){
                        for (int n=0; n<=m; n++){
                            kramerskronigtransform(dat.realepsij[i][j][q][m][n],dat.imagepsij[i][j][q][m][n],dat.freq,dat.dfreq);
                            // if (m==n && i==j){
                            //     for (int f=0; f<NFREQ; f++){
                            //         dat.realepsij[i][j][q][m][n][f] += 1;
                            //     }
                            // }
                        }
                    }
                }
            }
        }

        // free dynamically allocated memory
        for (int k=0; k<NKPT; k++){
            for (int m=0; m<NEPS; m++){
                for (int c=0; c<NBAND_C[k]; c++){
                    delete [] oint[k][m][c];
                }
            }
        }
        for (int k=0; k<NKPT; k++){
            for (int v=0; v<NBAND_V[k]; v++){
                delete [] C_kq[k][v];
            }
            delete [] C_kq[k];
        }
    }// loop q

    // free dynamically allocated memory
    for (int k=0; k<NKPT; k++){
        for (int m=0; m<NEPS; m++){
            delete [] oint[k][m];
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


void chi_tensor_LF(env &dat){

    // use indirect function calls for diracdelta function
    if (DELTA) {
        diracdelta = &diracdelta_lorentzian;
    }
    else {
        diracdelta = &diracdelta_gaussian;
    }

    // scalefactor
    double SCALEFACTOR = 128.0*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));
    
    // initialize susceptibility tensor matrix
    dat.realepsij = new double *****[3];
    dat.imagepsij = new double *****[3];
    for (int i=0; i<3; i++){
        dat.realepsij[i] = new double ****[3];
        dat.imagepsij[i] = new double ****[3];
        for (int j=0; j<3; j++){
            dat.realepsij[i][j] = new double ***[NQ];
            dat.imagepsij[i][j] = new double ***[NQ];
            for (int q=0; q<NQ; q++){
                dat.realepsij[i][j][q] = new double **[NEPS];
                dat.imagepsij[i][j][q] = new double **[NEPS];
                for (int m=0; m<NEPS; m++){
                    dat.realepsij[i][j][q][m] = new double *[NEPS];
                    dat.imagepsij[i][j][q][m] = new double *[NEPS];
                    for (int n=0; n<=m; n++){
                        dat.realepsij[i][j][q][m][n] = new double [NFREQ];
                        dat.imagepsij[i][j][q][m][n] = new double [NFREQ];
                    }
                }
            }
        }
    }

    //======================================================================
    int NBAND_V [NKPT]; // number of valence bands
    double **E_k; // Energy at k
    std::complex<double> ***C_k; // eigenvector at k
    E_k = new double *[NKPT];
    C_k = new std::complex<double> **[NKPT];
    #pragma omp parallel for
    for (int k=0; k<NKPT; k++){
        // vector{k}
        Eigen::Vector3d K = {dat.K[k][0], dat.K[k][1], dat.K[k][2]};

        // Hamiltonian at k
        Eigen::MatrixXcd H = HamiltonianEPM(dat.G,K,dat.epm);   
        
        // compute eigenvalue problem at k
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);
        if (eigsolver_k.info() != Eigen::Success){ printf("ERROR: Eigensolver failed to solve the matrix."); abort(); }

        E_k[k] = new double [NBAND];
        for (int n=0; n<NBAND; n++){
            E_k[k][n] = eigsolver_k.eigenvalues()(n)-dat.energyoffset;
            if (E_k[k][n] > 0) {
                NBAND_V[k] = n;
                break;
            }
        }
        C_k[k] = new std::complex<double> *[NBAND_V[k]];
        for (int n=0; n<NBAND_V[k]; n++){
            C_k[k][n] = new std::complex<double> [NPW];
            for (int j=0; j<NPW; j++){
                C_k[k][n][j] = eigsolver_k.eigenvectors().col(n)(j);
            }
        }
    }
    #pragma omp barrier
    

    //======================================================================
    // initialize temporary variables for k+q
    // overlap integral, O(k,q,m,v,c) = \vec{u}_{Q+G_m} * <k,v|e^{-i*g_m}j_0|k+q+g_m,c>
    std::complex<double> *****oint;
    oint = new std::complex<double> ****[NKPT]; // overlap integral [NK x 3 x NEPS x NV x NC]
    for (int k=0; k<NKPT; k++){
        oint[k] = new std::complex<double> ***[3];
        for (int i=0; i<3; i++){
            oint[k][i] = new std::complex<double> **[NEPS];
            for (int m=0; m<NEPS; m++){
                oint[k][i][m] = new std::complex<double> *[NBAND_V[k]];
            }
        }
    }

    // initialize E_kqm [NKPT x NEPS x NPW]
    int NBAND_CB [NKPT]; // number of conduction bands
    double E_kq [NKPT][NBAND]; // energy at (k+q+g_m,c)
    // double **E_kqm; // energy at (k+q+g_m,c)
    // E_kqm = new double *[NKPT];
    // for (int k=0; k<NKPT; k++){
    //     E_kqm[k] = new double [NBAND];
    // }
    

    for (int q=0; q<NQ; q++){
        Eigen::Vector3d Q = dat.Q[q];

        // construct overlap integrals
        #pragma omp parallel for
        for (int k=0; k<NKPT; k++){
            
            // get eigenvectors
            Eigen::Vector3d K = {dat.K[k][0], dat.K[k][1], dat.K[k][2]}; // vector{k}
            std::complex<double> C_kqgm [NEPS][NBAND][NPW];
            for (int m=0; m<NEPS; m++){
                Eigen::MatrixXcd H = HamiltonianEPM(dat.G,K+Q+dat.G[m],dat.epm);
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_kq(H);
                if (m==0){
                    for (int c=0; c<NBAND; c++){
                        E_kq[k][c] = double(eigsolver_kq.eigenvalues()(NBAND-1-c))-dat.energyoffset; // "NBAND-1-c" from behind
                        if (E_kq[k][c]<=0){
                            NBAND_CB[k] = c; // number of conduction bands
                            break;
                        }
                    }// loop over c
                }
                for (int c=0; c<NBAND_CB[k]; c++){
                    for (int p=0; p<NPW; p++){
                        C_kqgm[m][c][p] = eigsolver_kq.eigenvectors().col(NBAND-1-c)(p); // "NBAND-1-c" from behind
                    }
                }// loop over c
            }// loop over m

            // phase correction (if possible?)



            // get overlap integrals
            for (int m=0; m<NEPS; m++){
                std::vector<Eigen::Vector3cd> uvec_m = unitvec(dat.G[m]+Q);
                for (int i=0; i<3; i++){
                    for (int v=0; v<NBAND_V[k]; v++){ // valence bands
                        oint[k][i][m][v] = new std::complex<double> [NBAND_CB[k]];
                        for (int c=0; c<NBAND_CB[k]; c++){ // conduction bands
                            std::complex<double> oint_k_kq = 0.0; // <k,v|e^{-i*G_m}j_0|k+q,c> * t_{G_m+Q}
                            for (int p=0; p<NPW; p++){
                                oint_k_kq += conj(C_k[k][v][p]) * C_kqgm[m][c][p] * uvec_m[i].dot(dat.G[p]+K+Q/2+dat.G[m]/2); //+dat.G[m]/2
                            }// loop over p
                            oint[k][i][m][v][c] = oint_k_kq;
                        }// loop over c
                    }// loop over v
                }// loop over i
            }// loop over m
        }// loop over k (parallel)
        #pragma omp barrier



        // construct the susceptibility tensor matrix
        for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
            #pragma omp parallel for //collapse(3)
            for (int f=0; f<NFREQ; f++){
                // temporary variables
                double tmpval_real[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
                double tmpval_imag[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
                double depsilon = 0; // E_{k+q,c}-E_{k,v}
                
                for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                    
                    // sum over k-grids & bands
                    // \sum_{k,c,v} <k,v|e^{-i*(q+g_m)*r}|k+q+g_m,c> <k+q+g_n,c|e^{i*(q+g_n)*r)}|k,v>
                    for (int k=0; k<NKPT; k++){
                        for (int v=0; v<NBAND_V[k]; v++){ for (int c=0; c<NBAND_CB[k]; c++){
                            depsilon = E_kq[k][c]-E_k[k][v];
                            // imaginary part
                            tmpval_imag[i][j] += dat.KW[k] * (oint[k][i][m][v][c] * conj(oint[k][j][n][v][c])).real()
                                                * (*diracdelta)(depsilon-dat.freq[f]);
                            // real part
                            if (dat.kk==0){
                                tmpval_real[i][j] += dat.KW[k] * (oint[k][i][m][v][c] * conj(oint[k][j][n][v][c])).real()
                                                / (depsilon) / (depsilon*depsilon-dat.freq[f]*dat.freq[f]);
                            }
                        }}
                    }// loop over k
                }}// loop over i,j

                // permittivity
                for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                    if (dat.freq[f]>0){
                        dat.imagepsij[i][j][q][m][n][f] = SCALEFACTOR * pi * tmpval_imag[i][j] / (dat.freq[f]*dat.freq[f]);
                    }else{
                        dat.imagepsij[i][j][q][m][n][f] = 0;
                    }
                }}
                if (dat.kk==0){
                    for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                        dat.realepsij[i][j][q][m][n][f] = SCALEFACTOR * 2 * tmpval_real[i][j];
                        // if (m==n && i==j) dat.realepsij[i][j][q][m][n][f] += 1;
                    }}
                }
            }// loop over f
        }}// loop over m,n
        #pragma omp barrier

        // kramer-kronig transform
        if (dat.kk){
            for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                for (int m=0; m<NEPS; m++){for (int n=0; n<=m; n++){
                    kramerskronigtransform(dat.realepsij[i][j][q][m][n],dat.imagepsij[i][j][q][m][n],dat.freq,dat.dfreq);
                    // if (m==n && i==j){
                    //     for (int f=0; f<NFREQ; f++){
                    //         dat.realepsij[i][j][q][m][n][f] += 1;
                    //     }
                    // }
                }}// loop over m,n
            }}// loop over i,j
        }

        // free dynamically allocated memory
        for (int k=0; k<NKPT; k++){
            for (int i=0; i<3; i++){
                for (int m=0; m<NEPS; m++){
                    for (int v=0; v<NBAND_V[k]; v++){
                        delete [] oint[k][i][m][v];
                    }
                }
            }
        }
    }// loop q

    // free dynamically allocated memory
    for (int k=0; k<NKPT; k++){
        for (int i=0; i<3; i++){
            for (int m=0; m<NEPS; m++){
                delete [] oint[k][i][m];
            }
            delete [] oint[k][i];
        }
        delete [] oint[k];
    }
    delete [] oint;

    for (int k=0; k<NKPT; k++){
        for (int m=0; m<NBAND_V[k]; m++){
            delete [] C_k[k][m];
        }
        delete [] C_k[k];
        delete [] E_k[k];
    }
    delete [] E_k;
    delete [] C_k;

    return;
}


void chi_tensor_LF2(env &dat){

    // use indirect function calls for diracdelta function
    if (DELTA) {
        diracdelta = &diracdelta_lorentzian;
    }
    else {
        diracdelta = &diracdelta_gaussian;
    }

    // scalefactor
    double SCALEFACTOR = 128.0*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));
    
    // initialize susceptibility tensor matrix
    dat.realepsij = new double *****[3];
    dat.imagepsij = new double *****[3];
    for (int i=0; i<3; i++){
        dat.realepsij[i] = new double ****[3];
        dat.imagepsij[i] = new double ****[3];
        for (int j=0; j<3; j++){
            dat.realepsij[i][j] = new double ***[NQ];
            dat.imagepsij[i][j] = new double ***[NQ];
            for (int q=0; q<NQ; q++){
                dat.realepsij[i][j][q] = new double **[NEPS];
                dat.imagepsij[i][j][q] = new double **[NEPS];
                for (int m=0; m<NEPS; m++){
                    dat.realepsij[i][j][q][m] = new double *[NEPS];
                    dat.imagepsij[i][j][q][m] = new double *[NEPS];
                    for (int n=0; n<=m; n++){
                        dat.realepsij[i][j][q][m][n] = new double [NFREQ];
                        dat.imagepsij[i][j][q][m][n] = new double [NFREQ];
                    }
                }
            }
        }
    }

    //======================================================================
    int NBAND_V [NKPT]; // number of valence bands
    double **E_k; // Energy at k [NKPT x NVB]
    E_k = new double *[NKPT];
    std::complex<double> ****C_k; // eigenvector at k [NKPT x NEPS x NVB x NPW]
    C_k = new std::complex<double> ***[NKPT];
    for (int k=0; k<NKPT; k++){
        E_k[k] = new double [NBAND];
        C_k[k] = new std::complex<double> **[NEPS];
    }// loop over k

    std::cout << "  Solving for |k,v>..." << std::endl;
    #pragma omp parallel for
    for (int k=0; k<NKPT; k++){
        // vec{k}
        Eigen::Vector3d K = {dat.K[k][0], dat.K[k][1], dat.K[k][2]};
        for (int m=0; m<NEPS; m++){
            // Hamiltonian at k-g_m
            Eigen::MatrixXcd H = HamiltonianEPM(dat.G,K+dat.G[m],dat.epm);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);
            if (m==0){
                for (int v=0; v<NBAND; v++){
                    E_k[k][v] = eigsolver_k.eigenvalues()(v)-dat.energyoffset;
                    if (E_k[k][v] > 0) {
                        NBAND_V[k] = v;
                        break;
                    }
                }// loop over v
            }
            C_k[k][m] = new std::complex<double> *[NBAND_V[k]];
            for (int v=0; v<NBAND_V[k]; v++){
                C_k[k][m][v] = new std::complex<double> [NPW];
                for (int p=0; p<NPW; p++){
                    C_k[k][m][v][p] = eigsolver_k.eigenvectors().col(v)(p);
                }// loop over p
            }// loop over v
        }// loop over m

        // phase compensation?
        double phase;
        std::complex<double> im = std::complex<double>(0.0,1.0);
        for (int m=0; m<NEPS; m++){
            for (int v=0; v<NBAND_V[k]; v++){
                // int idx;
                // for (int p=0; p<NPW; p++){
                //     if (abs((dat.G[m]+dat.G[p]).norm())<1e-10){
                //         idx = p;
                //         break;
                //     }
                // }
                phase = std::arg(C_k[k][m][v][m]);
                for (int p=0; p<NPW; p++){
                    C_k[k][m][v][p] *= exp(-im*phase);
                }
            }
        }

    }// loop over k
    #pragma omp barrier

    //======================================================================
    // initialize temporary variables for k+q
    // overlap integral, O(k,q,m,v,c) = \vec{u}_{Q+G_m} * <k,v|e^{-i*g_m}j_0|k+q+g_m,c>
    std::complex<double> *****oint;
    oint = new std::complex<double> ****[NKPT]; // overlap integral [NK x 3 x NEPS x NV x NC]
    for (int k=0; k<NKPT; k++){
        oint[k] = new std::complex<double> ***[3];
        for (int i=0; i<3; i++){
            oint[k][i] = new std::complex<double> **[NEPS];
            for (int m=0; m<NEPS; m++){
                oint[k][i][m] = new std::complex<double> *[NBAND_V[k]];
            }
        }
    }

    std::cout << "  Solving for susceptibility tensor matrix elements..." << std::endl;
    int NBAND_CB [NKPT]; // number of conduction bands
    double E_kq [NKPT][NBAND]; // energy at (k+q+g_m,c)
    for (int q=0; q<NQ; q++){
        std::cout << "  (" << (q+1) << "/" << NQ << ")..." << std::endl; // (q/NQ)
        Eigen::Vector3d Q = dat.Q[q];

        // construct overlap integrals
        #pragma omp parallel for
        for (int k=0; k<NKPT; k++){
            
            // Hamiltonian at {k+q} and eigenvalue decomposition
            Eigen::Vector3d K = {dat.K[k][0], dat.K[k][1], dat.K[k][2]}; // vector{k}
            Eigen::MatrixXcd H = HamiltonianEPM(dat.G,K+Q,dat.epm);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_kq(H);
            // get the energy at {k+q} and the number of conduction bands
            for (int c=0; c<NBAND; c++){
                E_kq[k][c] = double(eigsolver_kq.eigenvalues()(NBAND-1-c))-dat.energyoffset; // "NBAND-1-c" from behind
                if (E_kq[k][c]<=0){
                    NBAND_CB[k] = c; // number of conduction bands
                    break;
                }
            }// loop over c
            // get the eigenvectors at {k+q}
            std::complex<double> C_kq [NBAND_CB[k]][NPW];
            for (int c=0; c<NBAND_CB[k]; c++){
                for (int p=0; p<NPW; p++){
                    C_kq[c][p] = eigsolver_kq.eigenvectors().col(NBAND-1-c)(p); // "NBAND-1-c" from behind
                }
            }// loop over c

            // tangential vectors
            Eigen::Vector3cd **uvec_m = uvecqm(Q,dat.G);

            // get the overlap integrals
            for (int m=0; m<NEPS; m++){
                // std::vector<Eigen::Vector3cd> uvec_m = unitvec(Q-dat.G[m]);
                // std::vector<Eigen::Vector3cd> uvec_m = unitvecq(Q-dat.G[m],Q);
                for (int i=0; i<3; i++){
                    for (int v=0; v<NBAND_V[k]; v++){ // valence bands
                        oint[k][i][m][v] = new std::complex<double> [NBAND_CB[k]];
                        for (int c=0; c<NBAND_CB[k]; c++){ // conduction bands
                            std::complex<double> oint_k_kq = 0.0; // <k+g_m,v|e^{-i*(q-g_m)}j_0|k+q,c> * u_{q+g_m}
                            for (int p=0; p<NPW; p++){
                                // oint_k_kq += conj(C_k[k][m][v][p]) * C_kq[c][p] * (dat.G[p]+K+Q/2+dat.G[m]/2).dot(uvec_m[i]);
                                oint_k_kq += conj(C_k[k][m][v][p]) * C_kq[c][p] * (dat.G[p]+K+Q/2+dat.G[m]/2).dot(uvec_m[m][i]);
                                // if (dat.loci[m][p]!=-1){
                                //     oint_k_kq += conj(C_k[k][0][v][dat.loci[m][p]]) * C_kq[c][p] * (dat.G[p]+K+Q/2+dat.G[m]/2).dot(uvec_m[m][i]);
                                // }
                            }// loop over p
                            oint[k][i][m][v][c] = oint_k_kq;
                        }// loop over c
                    }// loop over v
                }// loop over i
            }// loop over m
            
        }// loop over k (parallel)
        #pragma omp barrier

        // construct the susceptibility tensor matrix
        for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
            #pragma omp parallel for //collapse(3)
            for (int f=0; f<NFREQ; f++){
                // temporary variables
                double tmpval_real[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
                double tmpval_imag[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
                double dE = 0; // E_{k+q,c}-E_{k,v}
                
                for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                    
                    // sum over k-grids & bands
                    // \sum_{k,c,v} u_{q+g_m} * <k+g_m,v|e^{-i*(q-g_m)*r}|k+q,c> <k+q,c|e^{i*(q-g_n)*r)}|k+g_n,v> * u_{q+g_n}
                    for (int k=0; k<NKPT; k++){
                        for (int v=0; v<NBAND_V[k]; v++){ for (int c=0; c<NBAND_CB[k]; c++){
                            dE = E_kq[k][c]-E_k[k][v];
                            // imaginary part
                            tmpval_imag[i][j] += dat.KW[k] * (oint[k][i][m][v][c] * conj(oint[k][j][n][v][c])).real()
                                                * (*diracdelta)(dE-dat.freq[f]);
                            // real part
                            if (dat.kk==0){
                                tmpval_real[i][j] += dat.KW[k] * (oint[k][i][m][v][c] * conj(oint[k][j][n][v][c])).real()
                                                / (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                            }
                        }}
                    }// loop over k
                }}// loop over i,j

                // permittivity
                for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                    if (dat.freq[f]>0){
                        dat.imagepsij[i][j][q][m][n][f] = SCALEFACTOR * pi * tmpval_imag[i][j] / (dat.freq[f]*dat.freq[f]);
                    }else{
                        dat.imagepsij[i][j][q][m][n][f] = 0;
                    }
                }}
                if (dat.kk==0){
                    for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                        dat.realepsij[i][j][q][m][n][f] = SCALEFACTOR * 2 * tmpval_real[i][j];
                    }}
                }
            }// loop over f
        }}// loop over m,n
        #pragma omp barrier

        // kramer-kronig transform
        if (dat.kk){
            for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                for (int m=0; m<NEPS; m++){for (int n=0; n<=m; n++){
                    kramerskronigtransform(dat.realepsij[i][j][q][m][n],dat.imagepsij[i][j][q][m][n],dat.freq,dat.dfreq);
                }}// loop over m,n
            }}// loop over i,j
        }

        // free dynamically allocated memory
        for (int k=0; k<NKPT; k++){
            for (int i=0; i<3; i++){
                for (int m=0; m<NEPS; m++){
                    for (int v=0; v<NBAND_V[k]; v++){
                        delete [] oint[k][i][m][v];
                    }
                }
            }
        }
    }// loop q

    // free dynamically allocated memory
    for (int k=0; k<NKPT; k++){
        for (int i=0; i<3; i++){
            for (int m=0; m<NEPS; m++){
                delete [] oint[k][i][m];
            }
            delete [] oint[k][i];
        }
        delete [] oint[k];
    }
    delete [] oint;

    for (int k=0; k<NKPT; k++){
        for (int m=0; m<NEPS; m++){
            for (int v=0; v<NBAND_V[k]; v++){
                delete [] C_k[k][m][v];
            }
            delete [] C_k[k][m];
        }
        delete [] C_k[k];
        delete [] E_k[k];
    }
    delete [] E_k;
    delete [] C_k;

    return;
}



void chi_tensor_LF3(env &dat){

    // use indirect function calls for diracdelta function
    if (DELTA) {
        diracdelta = &diracdelta_lorentzian;
    }
    else {
        diracdelta = &diracdelta_gaussian;
    }

    // scalefactor
    double SCALEFACTOR = 128.0*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));
    
    // initialize susceptibility tensor matrix
    dat.realepsij = new double *****[3];
    dat.imagepsij = new double *****[3];
    for (int i=0; i<3; i++){
        dat.realepsij[i] = new double ****[3];
        dat.imagepsij[i] = new double ****[3];
        for (int j=0; j<3; j++){
            dat.realepsij[i][j] = new double ***[NQ];
            dat.imagepsij[i][j] = new double ***[NQ];
            for (int q=0; q<NQ; q++){
                dat.realepsij[i][j][q] = new double **[NEPS];
                dat.imagepsij[i][j][q] = new double **[NEPS];
                for (int m=0; m<NEPS; m++){
                    dat.realepsij[i][j][q][m] = new double *[NEPS];
                    dat.imagepsij[i][j][q][m] = new double *[NEPS];
                    for (int n=0; n<=m; n++){
                        dat.realepsij[i][j][q][m][n] = new double [NFREQ];
                        dat.imagepsij[i][j][q][m][n] = new double [NFREQ];
                    }
                }
            }
        }
    }

    //======================================================================
    int NBAND_V [NKPT]; // number of valence bands
    double **E_k; // Energy at k [NKPT x NVB]
    E_k = new double *[NKPT];
    std::complex<double> ****C_k; // eigenvector at k [NKPT x NEPS x NVB x NPW]
    C_k = new std::complex<double> ***[NKPT];
    for (int k=0; k<NKPT; k++){
        E_k[k] = new double [NBAND];
        C_k[k] = new std::complex<double> **[NEPS];
    }// loop over k

    std::cout << "  Solving for |k,v>..." << std::endl;
    #pragma omp parallel for
    for (int k=0; k<NKPT; k++){
        // vec{k}
        Eigen::Vector3d K = {dat.K[k][0], dat.K[k][1], dat.K[k][2]};
        for (int m=0; m<NEPS; m++){
            // Hamiltonian at k-g_m
            Eigen::MatrixXcd H = HamiltonianEPM(dat.G,K+dat.G[m],dat.epm);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);
            if (m==0){
                for (int v=0; v<NBAND; v++){
                    E_k[k][v] = eigsolver_k.eigenvalues()(v)-dat.energyoffset;
                    if (E_k[k][v] > 0) {
                        NBAND_V[k] = v;
                        break;
                    }
                }// loop over v
            }
            C_k[k][m] = new std::complex<double> *[NBAND_V[k]];
            for (int v=0; v<NBAND_V[k]; v++){
                C_k[k][m][v] = new std::complex<double> [NPW];
                for (int p=0; p<NPW; p++){
                    C_k[k][m][v][p] = eigsolver_k.eigenvectors().col(v)(p);
                }// loop over p
            }// loop over v
        }// loop over m

        // phase compensation?
        double phase;
        std::complex<double> im = std::complex<double>(0.0,1.0);
        for (int m=0; m<NEPS; m++){
            for (int v=0; v<NBAND_V[k]; v++){
                // int idx;
                // for (int p=0; p<NPW; p++){
                //     if (abs((dat.G[m]+dat.G[p]).norm())<1e-10){
                //         idx = p;
                //         break;
                //     }
                // }
                phase = std::arg(C_k[k][m][v][m]);
                for (int p=0; p<NPW; p++){
                    C_k[k][m][v][p] *= exp(-im*phase);
                }
            }
        }

    }// loop over k
    #pragma omp barrier

    //======================================================================
    // initialize temporary variables for k+q
    // overlap integral, O(k,q,m,v,c) = \vec{u}_{Q+G_m} * <k,v|e^{-i*g_m}j_0|k+q+g_m,c>
    std::complex<double> *******oint;
    oint = new std::complex<double> ******[NKPT]; // overlap integral [NK x 3 x 3 x NEPS x NEPS x NV x NC]
    for (int k=0; k<NKPT; k++){
        oint[k] = new std::complex<double> *****[3];
        for (int i=0; i<3; i++){
            oint[k][i] = new std::complex<double> ****[3];
            for (int j=0; j<3; j++){
                oint[k][i][j] = new std::complex<double> ***[NEPS];
                for (int m=0; m<NEPS; m++){
                    oint[k][i][j][m] = new std::complex<double> **[NEPS];
                    for (int n=0; n<NEPS; n++){
                        oint[k][i][j][m][n] = new std::complex<double> *[NBAND_V[k]];
                    }
                }
            }
            
        }
    }

    std::cout << "  Solving for susceptibility tensor matrix elements..." << std::endl;
    int NBAND_CB [NKPT]; // number of conduction bands
    double E_kq [NKPT][NBAND]; // energy at (k+q+g_m,c)
    for (int q=0; q<NQ; q++){
        std::cout << "  (" << (q+1) << "/" << NQ << ")..." << std::endl; // (q/NQ)
        Eigen::Vector3d Q = dat.Q[q];

        // construct overlap integrals
        #pragma omp parallel for
        for (int k=0; k<NKPT; k++){
            
            // Hamiltonian at {k+q} and eigenvalue decomposition
            Eigen::Vector3d K = {dat.K[k][0], dat.K[k][1], dat.K[k][2]}; // vector{k}
            Eigen::MatrixXcd H = HamiltonianEPM(dat.G,K+Q,dat.epm);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_kq(H);
            // get the energy at {k+q} and the number of conduction bands
            for (int c=0; c<NBAND; c++){
                E_kq[k][c] = double(eigsolver_kq.eigenvalues()(NBAND-1-c))-dat.energyoffset; // "NBAND-1-c" from behind
                if (E_kq[k][c]<=0){
                    NBAND_CB[k] = c; // number of conduction bands
                    break;
                }
            }// loop over c
            // get the eigenvectors at {k+q}
            std::complex<double> C_kq [NBAND_CB[k]][NPW];
            for (int c=0; c<NBAND_CB[k]; c++){
                for (int p=0; p<NPW; p++){
                    C_kq[c][p] = eigsolver_kq.eigenvectors().col(NBAND-1-c)(p); // "NBAND-1-c" from behind
                }
            }// loop over c

            // tangential vectors
            Eigen::Vector3cd **uvec_m = uvecqm(Q,dat.G);

            // get the overlap integrals
            for (int m=0; m<NEPS; m++){ for (int n=0; n<NEPS; n++){
                for (int i=0; i<3; i++){ for (int j=0; j<3; j++){
                    for (int v=0; v<NBAND_V[k]; v++){ // valence bands
                        oint[k][i][j][m][n][v] = new std::complex<double> [NBAND_CB[k]];
                        for (int c=0; c<NBAND_CB[k]; c++){ // conduction bands
                            Eigen::Vector3cd oint_k_kq = {0,0,0}; // <k+g_m,v|e^{-i*(q-g_m)}j_0|k+q,c>
                            for (int p=0; p<NPW; p++){
                                oint_k_kq += (conj(C_k[k][m][v][p]) * C_kq[c][p]) * (dat.G[p]+K+Q/2+dat.G[m]/2); // .dot(uvec_m[m][i])
                            }// loop over p

                            Eigen::Vector3cd oint_kq_k = {0,0,0};
                            for (int p=0; p<NPW; p++){
                                oint_kq_k += (C_k[k][n][v][p] * conj(C_kq[c][p])) * (dat.G[p]+K+Q/2+dat.G[n]/2); // .dot(uvec_m[n][i])
                            }// loop over p

                            Eigen::Matrix3cd Ointkk = oint_k_kq * oint_kq_k.transpose();

                            oint[k][i][j][m][n][v][c] = uvec_m[m][i].adjoint() * Ointkk.real() * uvec_m[n][j];
                        }// loop over c
                    }// loop over v
                }}// loop over i,j
            }}// loop over m,n
            
        }// loop over k (parallel)
        #pragma omp barrier

        // construct the susceptibility tensor matrix
        for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
            #pragma omp parallel for //collapse(3)
            for (int f=0; f<NFREQ; f++){
                // temporary variables
                double tmpval_real[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
                double tmpval_imag[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
                double dE = 0; // E_{k+q,c}-E_{k,v}
                
                for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                    
                    // sum over k-grids & bands
                    // \sum_{k,c,v} u_{q+g_m} * <k+g_m,v|e^{-i*(q-g_m)*r}|k+q,c> <k+q,c|e^{i*(q-g_n)*r)}|k+g_n,v> * u_{q+g_n}
                    for (int k=0; k<NKPT; k++){
                        for (int v=0; v<NBAND_V[k]; v++){ for (int c=0; c<NBAND_CB[k]; c++){
                            dE = E_kq[k][c]-E_k[k][v];
                            // imaginary part
                            tmpval_imag[i][j] += dat.KW[k] * oint[k][i][j][m][n][v][c].real()
                                                * (*diracdelta)(dE-dat.freq[f]);
                            // real part
                            if (dat.kk==0){
                                tmpval_real[i][j] += dat.KW[k] * oint[k][i][j][m][n][v][c].real()
                                                / (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                            }
                        }}
                    }// loop over k
                }}// loop over i,j

                // permittivity
                for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                    if (dat.freq[f]>0){
                        dat.imagepsij[i][j][q][m][n][f] = SCALEFACTOR * pi * tmpval_imag[i][j] / (dat.freq[f]*dat.freq[f]);
                    }else{
                        dat.imagepsij[i][j][q][m][n][f] = 0;
                    }
                }}
                if (dat.kk==0){
                    for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                        dat.realepsij[i][j][q][m][n][f] = SCALEFACTOR * 2 * tmpval_real[i][j];
                    }}
                }
            }// loop over f
        }}// loop over m,n
        #pragma omp barrier

        // kramer-kronig transform
        if (dat.kk){
            for (int i=0; i<3; i++){for (int j=0; j<3; j++){
                for (int m=0; m<NEPS; m++){for (int n=0; n<=m; n++){
                    kramerskronigtransform(dat.realepsij[i][j][q][m][n],dat.imagepsij[i][j][q][m][n],dat.freq,dat.dfreq);
                }}// loop over m,n
            }}// loop over i,j
        }

        // free dynamically allocated memory
        for (int k=0; k<NKPT; k++){
            for (int i=0; i<3; i++){ for (int j=0; j<3; j++){
                for (int m=0; m<NEPS; m++){ for (int n=0; n<NEPS; n++){
                    for (int v=0; v<NBAND_V[k]; v++){
                        delete [] oint[k][i][j][m][n][v];
                    }
                }}
            }}
        }
    }// loop q

    // free dynamically allocated memory
    for (int k=0; k<NKPT; k++){
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                for (int m=0; m<NEPS; m++){
                    for (int n=0; n<NEPS; n++){
                        delete [] oint[k][i][j][m][n];
                    }
                    delete [] oint[k][i][j][m];
                }
                delete [] oint[k][i][j];
            }
            delete [] oint[k][i];
        }
        delete [] oint[k];
    }
    delete [] oint;

    for (int k=0; k<NKPT; k++){
        for (int m=0; m<NEPS; m++){
            for (int v=0; v<NBAND_V[k]; v++){
                delete [] C_k[k][m][v];
            }
            delete [] C_k[k][m];
        }
        delete [] C_k[k];
        delete [] E_k[k];
    }
    delete [] E_k;
    delete [] C_k;

    return;
}

// // Original?
// inline Eigen::Vector3d perpvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (v(1)!=0){
//         p << -v(1), v(0), 0;
//         p = p/p.norm();
//     }else if (v(1)==0 && v(2)!=0){
//         p << v(2), 0, -v(0);
//         p = p/p.norm();
//     }else if (v(1)==0 && v(2)==0){
//         p << 0, 1, 0;
//         // p = p/p.norm();
//     }
//     return p;
// }

// G-X (100) direction!!!
inline Eigen::Vector3d perpvec(Eigen::Vector3d v){
    Eigen::Vector3d p;
    if (v(2)!=0){
        p << 0, -v(2), v(1);
    }else if (v(2)==0 && v(1)!=0){
        p << -v(1), v(0), 0;
    }else if (v(1)==0 && v(2)==0){
        p << 0, 0, 1;
    }
    p = p/p.norm();
    return p;
}




// Original?
inline Eigen::Vector3d tanvec(Eigen::Vector3d v){
    Eigen::Vector3d p;
    if (v(1)!=0){
        p << -v(1), v(0), 0;
    }else if (v(1)==0 && v(2)!=0){
        p << v(2), 0, -v(0);
    }else if (v(1)==0 && v(2)==0){
        p << 0, 1, 0;
    }
    p = p/p.norm();

    return p;
}


// inline Eigen::Vector3d tanvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (v(2)!=0){
//         p << 0, -v(2), v(1);
//     }else if (v(2)==0 && v(1)!=0){
//         p << -v(1), v(0), 0;
//     }else if (v(1)==0 && v(2)==0){
//         p << 0, 1, 1;
//     }
//     p = p/p.norm();

//     return p;
// }



// // y-axis
// Eigen::Vector3d tanvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (v(2)!=0){
//         p << v(2), 0, -v(0);
//     }else if (v(2)==0 && v(1)!=0){
//         p << -v(1), v(0), 0;
//     }else if (v(2)==0 && v(1)==0){
//         p << 0, 1, 0;
//     }
//     p = p/p.norm();
//     return p;
// }


// // test
// inline Eigen::Vector3d tanvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     double tx,ty,tz;
//     if (v(2)==0 && v(2)==0){
//         p << 0,1,0;
//     }else{
//         tx = -(pow(v(1),2)+pow(v(2),2));
//         ty = v(1)*v(0);
//         tz = v(2)*v(0);
//         p << tx,ty,tz;
//         if (v(2)<0) p = -p;
//     }
//     p = p/p.norm();

//     return p;
// }



// normal vector
Eigen::Vector3d normvec(Eigen::Vector3d v){
    if (v.norm()==0){
        v << 1.0, 0.0, 0.0;
    }else{
        v = v/v.norm();
    }
    return v;
}


// // tangential vector
// Eigen::Vector3d tanvec(Eigen::Vector3d v, Eigen::Vector3d q){
//     Eigen::Vector3d p;
//     double tol = 1e-7;

//     if (abs(q(1))<tol && abs(q(2))<tol){ // [100]
//         if (v(2)!=0){
//             p << 0, -v(2), v(1); 
//         }else if (v(2)==0 && v(1)!=0){
//             p << -v(1), v(0), 0;
//         }else if (v(1)==0 && v(2)==0){
//             p << 0, 0, 1;
//         }
//     }else if (abs(q(0)-q(1))<tol && abs(q(1)-q(2))<tol && abs(q(2)-q(0))<tol){ // [111]
//         if (v(1)!=0){
//             p << -v(1), v(0), 0;
//         }else if (v(1)==0 && v(2)!=0){
//             p << v(2), 0, -v(0);
//         }else if (v(1)==0 && v(2)==0){
//             p << 0, 1, 0;
//         }
//     }else if (abs(q(0)-q(1))<tol && abs(q(2))<tol){ // [110]
//         if (v(1)!=0){
//             p << -v(1), v(0), 0;
//         }else if (v(1)==0 && v(2)!=0){
//             p << v(2), 0, -v(0);
//         }else if (v(1)==0 && v(2)==0){
//             p << 0, 1, 0;
//         }
//     }else{
//         if (v(1)!=0){
//             p << -v(1), v(0), 0;
//         }else if (v(1)==0 && v(2)!=0){
//             p << v(2), 0, -v(0);
//         }else if (v(1)==0 && v(2)==0){
//             p << 0, 1, 0;
//         }
//     }
//     p = p/p.norm();
//     return p;
// }

// tangential vector
Eigen::Vector3d tanvecq(Eigen::Vector3d nvec, Eigen::Vector3d qvec){
    Eigen::Vector3d tvec;
    Eigen::Vector3d tmp;
    double tol = 1e-7;
    double tx, ty, tz, K;

    if (abs(qvec(1))<tol && abs(qvec(2))<tol){// (100)
        if (abs(nvec(1))<tol && abs(nvec(2))<tol){
            tvec << 0,0,1;
        }else{
            tmp << 1,0,0;
            tvec = nvec.cross(tmp);
        }
    }else if (abs(qvec(0)-qvec(1))<tol && abs(qvec(0)-qvec(2))<tol){// (111)
        if (abs(nvec(1)-nvec(2))<tol && abs(nvec(0)-nvec(2))<tol){
            tmp << 1,0,0;
            tvec << nvec.cross(tmp);
        }else{
            tmp << 1,1,1;
            tvec = nvec.cross(tmp);
        }
    }else if (abs(qvec(0)-qvec(1))<tol && abs(qvec(2))<tol){// (110)
        if (abs(nvec(0)-nvec(1))<tol && abs(nvec(2))<tol){
            tvec << 0,0,1;
        }else{
            tmp << 1,1,0;
            tvec = nvec.cross(tmp);
        }
    }else{
        // if (v(1)!=0){
        //     p << -v(1), v(0), 0;
        // }else if (v(1)==0 && v(2)!=0){
        //     p << v(2), 0, -v(0);
        // }else if (v(1)==0 && v(2)==0){
        //     p << 0, 1, 0;
        // }
        tvec = tanvec(nvec);
    }

    // if (abs(v.dot(p))>tol){
    //     std::cout << "check tangential vector"; abort();
    // }


    tvec /= tvec.norm();

    return tvec;
}




// // return a perpendicular vector
// inline Eigen::Vector3d perpvec(Eigen::Vector3d v){
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




// // 2024.03.09 (100) direction
// inline Eigen::Vector3d perpvec(Eigen::Vector3d v){
//     v = v/v.norm(); // normalize
//     Eigen::Vector3d p;
//     double tx, ty, tz, K;
//     if (v(1)==0 && v(2)==0){
//         p << 0, 0, 1;
//     }else {
//         tx = 1;
//         ty = -0.5/v(2)-0.5*v(0)/v(1);
//         tz = 0.5/v(1)-0.5*v(0)/v(2);
//         // if (v(0)>0){
//         //     tx = 1;
//         //     ty = -0.5/v(2)-0.5*v(0)/v(1);
//         //     tz = 0.5/v(1)-0.5*v(0)/v(2);
//         // }else{
//         //     tx = -1;
//         //     ty = 0.5/v(2)+0.5*v(0)/v(1);
//         //     tz = -0.5/v(1)+0.5*v(0)/v(2);
//         // }
//         p << tx, ty, tz;
//     }
//     p = p/p.norm();
//     // p = p.cross(v);
//     return p;
// }


std::vector<Eigen::Vector3cd> unitvec (Eigen::Vector3d nvec){
    std::vector<Eigen::Vector3cd> uvec;
    uvec.reserve(3);

    // normal vector
    if (nvec.norm()==0){
        nvec << 1.0,0.0,0.0;
        uvec[0] = nvec;
    }else{
        nvec = nvec/nvec.norm();
        uvec[0] = nvec;
    }
    
    // tangential vectors
    Eigen::Vector3d tvec, pvec;
    tvec = tanvec(nvec);
    pvec = nvec.cross(tvec);

    std::complex<double> im = std::complex<double>(0.0,1.0);
    uvec[1] = (tvec+im*pvec)/sqrt(2);
    uvec[2] = (tvec-im*pvec)/sqrt(2);

    // double tol = 1e-10;
    // if (abs(nvec.dot(tvec))>tol || abs(nvec.dot(pvec))>tol){
    //     std::cout << "unnormal vector";
    //     abort();
    // }
    // uvec[1] = tvec;
    // uvec[2] = pvec;
    return uvec;
}


std::vector<Eigen::Vector3cd> unitvecq (Eigen::Vector3d nvec, Eigen::Vector3d qvec){
    std::vector<Eigen::Vector3cd> uvec;
    uvec.reserve(3);

    // normal vector
    if (nvec.norm()==0){
        nvec << 1.0,0.0,0.0;
        uvec[0] = nvec;
    }else{
        nvec = nvec/nvec.norm();
        uvec[0] = nvec;
    }
    
    // tangential vectors
    Eigen::Vector3d tvec, pvec;
    pvec = tanvecq(nvec,qvec); // ~ y
    tvec = pvec.cross(nvec); // x = cross(y,z)

    bool REAL = false;
    if (REAL){
        // uvec[1] = tvec;
        // uvec[2] = pvec;

        uvec[1] = (tvec-pvec)/sqrt(2);
        uvec[2] = (tvec+pvec)/sqrt(2);
    }else{
        std::complex<double> im = std::complex<double>(0.0,1.0);
        uvec[1] = (tvec+im*pvec)/sqrt(2);
        uvec[2] = (tvec-im*pvec)/sqrt(2);

        // uvec[1] = (tvec-pvec)/2+im*(tvec+pvec)/2;
        // uvec[2] = (tvec-pvec)/2-im*(tvec+pvec)/2;
    }
    return uvec;
}

Eigen::Vector3cd **uvecqm (Eigen::Vector3d Q, std::vector<Eigen::Vector3d> G){
    Eigen::Vector3cd **uvec;
    uvec = new Eigen::Vector3cd *[NEPS];
    for (int m=0; m<NEPS; m++) uvec[m] = new Eigen::Vector3cd [3];
    
    // normal vector
    for (int m=0; m<NEPS; m++){
        uvec[m][0] << Q+G[m];
        if (uvec[m][0].norm()==0){
            uvec[m][0] << 1.0,0.0,0.0;
        }else{
            uvec[m][0] /= uvec[m][0].norm(); // normalize
        }
    }

    // tangential vectors
    double tol = 1e-8;
    double tx, ty, tz, K;
    double s;
    if (abs(Q(1))<tol && abs(Q(2))<tol){ // [100]
        for (int m=0; m<NEPS; m++){
            if (abs(uvec[m][0](1))<tol && abs(uvec[m][0](2))<tol){
                uvec[m][1] << 0,1,0;
            }else{
                uvec[m][1](0) = -(pow(uvec[m][0](1),2)+pow(uvec[m][0](2),2));
                uvec[m][1](1) = uvec[m][0](1)*uvec[m][0](0);
                uvec[m][1](2) = uvec[m][0](2)*uvec[m][0](0);
                uvec[m][1] /= uvec[m][1].norm(); // normalize
            }
        }
    }else if (abs(Q(0)-Q(1))<tol && abs(Q(1)-Q(2))<tol){ // [111]
        for (int m=0; m<NEPS; m++){
            if (m==0){
                uvec[m][1] << 0.0,-1.0,1.0;
                uvec[m][1] /= uvec[m][1].norm();
            }else if (m==1){
                uvec[m][1] << 1.0,0.0,-1.0;
                uvec[m][1] /= uvec[m][1].norm();
            }else if (m==8){
                uvec[m][1] << -1.0,1.0,0.0;
                uvec[m][1] /= uvec[m][1].norm();
            }else{
                Eigen::Vector3d v = uvec[m][0].real();
                if (v(0)*v(1)*v(2)>0){
                    K = sqrt(3.0*pow(Q(0),2)-2.0*Q(0)+3);
                    if (v(0)>0){
                        tx = -(Q(0)-1)/K;
                        ty = (Q(0)+1+K)/2/K;
                        tz = (Q(0)+1-K)/2/K;
                    }else if (v(1)>0){
                        ty = -(Q(0)-1)/K;
                        tz = (Q(0)+1+K)/2/K;
                        tx = (Q(0)+1-K)/2/K;
                    }else{
                        tz = -(Q(0)-1)/K;
                        tx = (Q(0)+1+K)/2/K;
                        ty = (Q(0)+1-K)/2/K;
                    }
                    uvec[m][1] << tx,ty,tz;
                    uvec[m][1] /= uvec[m][1].norm();
                }else {
                    K = sqrt(3.0*pow(Q(0),2)+2.0*Q(0)+3);
                    if (v(0)<0){
                        tx = (Q(0)+1)/K;
                        ty = -(Q(0)-1+K)/2/K;
                        tz = (-Q(0)+1+K)/2/K;
                    }else if (v(1)<0){
                        ty = (Q(0)+1)/K;
                        tz = -(Q(0)-1+K)/2/K;
                        tx = (-Q(0)+1+K)/2/K;
                    }else{
                        tz = (Q(0)+1)/K;
                        tx = -(Q(0)-1+K)/2/K;
                        ty = (-Q(0)+1+K)/2/K;
                    }
                    uvec[m][1] << tx,ty,tz;
                    uvec[m][1] /= uvec[m][1].norm();
                }
            }
        }
    }else if (abs(Q(1)-Q(2))<tol && abs(Q(3)<tol)){// [110]
        for (int m=0; m<NEPS; m++){
            Eigen::Vector3d v = uvec[m][0].real();
            if (m==0){
                uvec[m][1] << 0.0,0.0,1.0;
            }else if (m==1 || m==2 || m==7 || m==8){
                s = G[m](0)*G[m](1)*G[m](2);
                tx = s/2+s*v(2)/2;
                ty = -s/2+s*v(2)/2;
                tz = -s*v(0);
                uvec[m][1] << tx,ty,tz;
            }else if (m==4 || m==6){
                s = G[m](0)*G[m](1)*G[m](2);
                K = sqrt(2*pow(Q(0),2)+3);
                tx = -s*(4*Q(0)*Q(0)+s*3*Q(0)+3-3*(Q(0)+s)*K);
                ty = s*(4*Q(0)*Q(0)-s*3*Q(0)+3-3*(Q(0)-s)*K);
                tz = s*2*(Q(0)*Q(0)+3);
                uvec[m][1] << tx,ty,tz;
            }else{
                s = G[m](0)*G[m](1)*G[m](2);
                K = sqrt(2*pow(Q(0),2)+3);
                tx = s*(4*Q(0)*Q(0)-s*3*Q(0)+3+3*(Q(0)+s)*K);
                ty = -s*(4*Q(0)*Q(0)+s*3*Q(0)+3+3*(Q(0)-s)*K);
                tz = s*2*(Q(0)*Q(0)+3);
                uvec[m][1] << tx,ty,tz;
            }
        }
    }else{
        for (int m=0; m<NEPS; m++){
            Eigen::Vector3d T;
            T << 0.0,0.0,1.0;
            uvec[m][1] = uvec[m][0].cross(T);
        }
    }

    //
    Eigen::Vector3cd t,p;
    std::complex<double> im = std::complex<double>(0.0,1.0);
    for (int m=0; m<NEPS; m++){
        uvec[m][1] /= uvec[m][1].norm();
        uvec[m][2] = uvec[m][0].cross(uvec[m][1]);

        t = uvec[m][1];
        p = uvec[m][2];
        uvec[m][1] = (t+im*p)/sqrt(2);
        uvec[m][2] = (t-im*p)/sqrt(2);
    }



    for (int m=0; m<NEPS; m++){
        if (abs(uvec[m][0].dot(uvec[m][1]))>tol){
            std::cout << "Warning::not normal found\n" << Q << "-" << m << std::endl;
            std::cout << uvec[m][0] << std::endl;
            std::cout << uvec[m][1] << std::endl;
            abort();
        }
        if (abs(uvec[m][0].dot(uvec[m][2]))>tol){
            std::cout << "not normal-" << Q << "-" << m << std::endl;
            std::cout << uvec[m][0] << std::endl;
            std::cout << uvec[m][2] << std::endl;
            abort();
        }
    }
    return uvec;
}

}