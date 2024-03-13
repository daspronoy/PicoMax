/*


*/
#include "pmx.hpp"

namespace pmx{

inline std::complex<double> (*diracdelta) (std::complex<double> value);

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
                    Eigen::Vector3d vecqGm_t = perpvec(Q+dat.G[m]);
                    Eigen::Vector3d vecqGn_t = perpvec(Q+dat.G[n]);
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
                            imagL1 += dat.KW[k] * (ointL[k][m][n][v][c]) 
                                                    * (*diracdelta)(depsilon-dat.freq[f]);
                            imagT1 += dat.KW[k] * (ointT[k][m][n][v][c]) 
                                                    * (*diracdelta)(depsilon-dat.freq[f]);
                            if (dat.kk==0){
                                realL1 += dat.KW[k] * (ointL[k][m][n][v][c]) 
                                                        * (depsilon) / (depsilon*depsilon-dat.freq[f]*dat.freq[f]);
                                realT1 += dat.KW[k] * (ointT[k][m][n][v][c]) 
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


// Original?
// Eigen::Vector3d perpvec(Eigen::Vector3d v){
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


// // y-axis
// Eigen::Vector3d perpvec(Eigen::Vector3d v){
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



// normal vector
Eigen::Vector3d normvec(Eigen::Vector3d v){
    if (v.norm()==0){
        v << 1, 0, 0;
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
Eigen::Vector3d tanvec(Eigen::Vector3d v, Eigen::Vector3d q){
    Eigen::Vector3d p;
    Eigen::Vector3d t;
    double tol = 1e-7;

    if (abs(q(1))<tol && abs(q(2))<tol){ // [100]
        if (abs(v(1))<tol && abs(v(2))<tol){
            p << 0,1,0;
        }else{
            t << 0,v(2),-v(1);
            p = v.cross(t);
        }
    }else if (abs(q(0)-q(1))<tol && abs(q(1)-q(2))<tol && abs(q(2)-q(0))<tol){ // [111]
        if (abs(v(1)-v(2))<tol && abs(v(2)-v(0))<tol && abs(v(0)-v(1))<tol){
            t << 1,0,0;
            p = v.cross(t);
        }else if (abs(v(0)-v(1))<tol){
            t << 1,1,0;
            p = v.cross(t);
        }else if (abs(v(1)-v(2))<tol){
            t << 0,1,1;
            p = v.cross(t);
        }else if (abs(v(2)-v(0))<tol){
            t << 1,0,1;
            p = v.cross(t);
        }
    }else if (abs(q(0)-q(1))<tol && abs(q(2))<tol){ // [110]
        if (abs(v(0)-v(1))<tol && abs(v(2))<tol){
            p << 0,0,1;
        }else{
            t << 1,1,0;
            p = v.cross(t);
        }
    }else{
        if (v(1)!=0){
            p << -v(1), v(0), 0;
        }else if (v(1)==0 && v(2)!=0){
            p << v(2), 0, -v(0);
        }else if (v(1)==0 && v(2)==0){
            p << 0, 1, 0;
        }
    }
    p = p/p.norm();
    return p;
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
// Eigen::Vector3d perpvec(Eigen::Vector3d v){
//     v = v/v.norm(); // normalize
//     Eigen::Vector3d p;
//     double tx, ty, tz, K;
//     if (v(1)==0 && v(2)==0){
//         p << 0, 1, -1;
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


// // dirac delta function
// inline std::complex<double> diracdelta(std::complex<double> x){
//     std::complex<double> delta_x;
//     if (DELTA) // lorentzian delta function
//         delta_x = EPSILON/pi/(pow(EPSILON,2)+pow(x,2));
//     else // gaussian delta function
//         delta_x = exp(-x*x/(EPSILON*EPSILON))/(EPSILON*sqrt(pi));
// 	return delta_x;
// }

inline std::complex<double> diracdelta_lorentzian(std::complex<double> x){
    return EPSILON/pi/(pow(EPSILON,2)+pow(x,2));
}

inline std::complex<double> diracdelta_gaussian(std::complex<double> x){
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

}