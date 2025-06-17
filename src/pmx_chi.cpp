/*
    pmx_chi.cpp | includes susceptibility calculation functinalities of picomax

    This file implements the calculation of the susceptibility tensor chi(q,omega) using 
    the empirical pseudopotential method (EPM). The susceptibility tensor describes the 
    linear response of a material to an external electromagnetic field.

    OVERVIEW OF CALCULATION FLOW:
    =============================
    
    1. INITIALIZATION:
       - Set up broadening functions (Gaussian, Lorentzian, etc.)
       - Define physical scaling factors
       - Allocate memory for susceptibility tensor storage
    
    2. ELECTRONIC STRUCTURE:
       - Calculate band structure at all k-points using EPM Hamiltonian H(k)
       - Separate valence and conduction bands based on energy
       - Store eigenvectors and eigenvalues for later use
    
    3. OVERLAP INTEGRALS:
       - For each q-vector, calculate H(k+q) band structure
       - Compute matrix elements ⟨k,c|e^{i(q+G)·r}|k+q,v⟩ between states
       - These represent electronic transitions between bands
    
    4. SUSCEPTIBILITY CALCULATION:
       - Use Fermi's Golden Rule with energy conservation
       - Apply broadening to discrete transitions
       - Sum over all k-points and band combinations
       - Apply appropriate scaling factors
    
    5. POST-PROCESSING:
       - Apply Kramers-Kronig relations if needed
       - Store results in output arrays
    
    KEY PHYSICS:
    ============
    - chi(q,ω) describes how electrons respond to electromagnetic fields
    - Peaks occur at inter-band transition energies  
    - Different tensor components (LL, LT, TT) represent different field orientations
    - SOC doubles the Hamiltonian size and modifies eigenvector structure

    Main functions:
    - chi_LL: Calculates longitudinal susceptibility using charge-based approach
    - chi_tensor: Full susceptibility tensor calculation
    - chi_xyz: Simplified xyz-component calculation
    - permittivity: Calculates permittivity matrix elements

    Author: Jungho Mun
    Date: June 6, 2024
*/


/*
    header files
*/
#include "pmx_chi.hpp"      // Function declarations for susceptibility calculations
#include "pmx_epm.hpp"      // Empirical pseudopotential method Hamiltonian
#include "pmx_basis.hpp"    // Basis set and G-vector utilities


/*
    function definitions
*/
namespace pmx{


/*
    SECTION 1: DIRAC DELTA FUNCTION IMPLEMENTATIONS
    
    These functions provide different approximations to the Dirac delta function δ(x)
    used for broadening discrete energy levels into continuous spectra.
    The function pointer 'diracdelta' allows runtime selection of the broadening scheme.
*/

/*
    Function pointer for dirac delta function (indirect function call)
    This allows dynamic selection of the delta function implementation at runtime
*/
inline double (*diracdelta) (double value);

/*
    Lorentzian dirac delta function
    
    Mathematical form: δ(x) = ε/π/(ε² + x²)
    
    This provides a Lorentzian broadening with full-width at half-maximum (FWHM) = 2ε
    Commonly used in spectroscopy as it has realistic line shapes with extended tails
    
    Parameters:
    - x: Energy difference (usually E_c - E_v - ħω)
    - EPSILON: Global smearing parameter controlling the width of the broadening
*/
inline double diracdelta_lorentzian(double x){
    return EPSILON/pi/(pow(EPSILON,2)+pow(x,2));
}

/*
    Gaussian dirac delta function
    
    Mathematical form: δ(x) = (1/(ε√(2π))) * exp(-x²/(2ε²))
    
    This provides Gaussian broadening, which is physically motivated by 
    thermal/phonon interactions and instrumental resolution effects
    
    Parameters:
    - x: Energy difference
    - EPSILON: Standard deviation of the Gaussian distribution
*/
inline double diracdelta_gaussian(double x){
    return exp(-0.5*x*x/(EPSILON*EPSILON))/(EPSILON*sqrt(2.0*pi));
}

/*
    Triangular dirac delta function
    
    Mathematical form: 
    δ(x) = (4/ε) * (0.5 - |x|/ε)  for |x| < ε/2
         = 0                      otherwise
    
    This provides a simple triangular broadening that goes to zero linearly.
    Less commonly used but computationally efficient.
    
    Parameters:
    - x: Energy difference
    - EPSILON: Half-width of the triangle base
*/
inline double diracdelta_triangle(double x){
    if (abs(x)<0.5*EPSILON){
        return 4.0/EPSILON*(0.5-abs(x)/EPSILON);
    }else{
        return 0;
    }
}

/*
    Rectangular (box) dirac delta function
    
    Mathematical form:
    δ(x) = 1/ε  for |x| < ε/2
         = 0    otherwise
    
    This provides uniform broadening within a fixed energy window.
    Simplest implementation but may cause artificial discontinuities.
    
    Parameters:
    - x: Energy difference  
    - EPSILON: Half-width of the rectangular window
*/
inline double diracdelta_rect(double x){
    if (abs(x)<0.5*EPSILON){
        return 1.0/EPSILON;
    }else{
        return 0;
    }
}

/*
    SECTION 2: KRAMERS-KRONIG TRANSFORMATION
    
    Hilbert transformation of susceptibility tensor matrix elements
    
    The Kramers-Kronig relations connect the real and imaginary parts of the 
    susceptibility tensor through causality constraints. This function implements
    the discrete Hilbert transform to compute:
    
    Re[χ(ω)] from Im[χ(ω')] and vice versa
    
    This is essential for obtaining the full frequency-dependent response from
    calculations that may only directly compute one component.

    References:
    [1] See the documentation Eq. XXX for explicit expressions
    [2] Phys. Rev. B 74, 035101 (2006)
    
    Parameters:
    - ReX: Array of real parts [input/output]
    - ImX: Array of imaginary parts [input/output] 
    - w: Frequency grid
    - dw: Frequency spacing
*/
void kramerskronigtransform(double *ReX, double *ImX, double *w, double dw){

    // STEP 1: Store original real part for later use in imaginary part calculation
    // This is necessary because we modify ReX[f] in the first loop, but need 
    // the original values for the second loop
    double ReX1 [NFREQ];
    for (int f=0; f<NFREQ; f++)
        ReX1[f] = ReX[f];

    // Discrete integration prefactor: dω * (2/π)
    // Factor of 2/π comes from the Kramers-Kronig relation integral
    double ds = dw*2.0/pi;

    // STEP 2: Generate real part from imaginary part using Kramers-Kronig relation
    // Re[χ(ω)] = (2/π) * P ∫ ω' Im[χ(ω')] / (ω'² - ω²) dω'
    // where P denotes the principal value integral
    for (int f=0; f<NFREQ; f++){
        double ReXf = 0.0;
        for (int s=0; s<NFREQ; s++){
            // Skip the pole at ω' = ω (principal value)
            if (s!=f)
                ReXf += w[s]*ImX[s]/(w[s]*w[s]-w[f]*w[f]);
        }
        ReX[f] += ds*ReXf;  // Add to existing real part
    }    // STEP 3: Generate imaginary part from real part using Kramers-Kronig relation  
    // Im[χ(ω)] = -(2/π) * P ∫ ω' Re[χ(ω')] / (ω'² - ω²) dω'
    // We use the original real parts (ReX1) before they were modified above
    for (int f=0; f<NFREQ; f++){
        double ImXf = 0.0;
        for (int s=0; s<NFREQ; s++){
            // Skip the pole at ω' = ω (principal value)
            if (s!=f)
                ImXf += w[s]*ReX1[s]/(w[s]*w[s]-w[f]*w[f]);
        }
        ImX[f] += ds*ImXf;  // Add to existing imaginary part
    }

    return;
}

/*    SECTION 3: LONGITUDINAL SUSCEPTIBILITY CALCULATION (CHARGE-BASED APPROACH)
    
    Calculate the susceptibility tensor matrix using an alternative expression
    which treats the longitudinal terms as charge density fluctuations.

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions
    [2] (Adler) “Quantum Theory of the Dielectric Constant in Real Solids”, Phys. Rev. 126, 413 (1962)
*/
void chi_LL(env &dat){
    // Timing variables for performance monitoring
    std::chrono::time_point<std::chrono::system_clock> time_0, time_1;
    std::chrono::duration<double> elapsed_time;

    // Set tensor dimension to 3 for full 3D calculation
    NTSR = 3;

    // STEP 1: Select broadening function based on user input
    // This determines how discrete energy levels are broadened into continuous spectra
    switch (dat.delta) {
        case 0:
            diracdelta = &diracdelta_gaussian;    // Gaussian broadening
            break;
        case 1:
            diracdelta = &diracdelta_lorentzian;  // Lorentzian broadening  
            break;
        case 2:
            diracdelta = &diracdelta_triangle;    // Triangular broadening
            break;
        case 3:
            diracdelta = &diracdelta_rect;        // Rectangular broadening
            break;
    }

    // STEP 2: Define physical scale factors for different tensor components
    // These convert from atomic units to practical units (eV, Angstrom, etc.)
    double SCALEFACTOR = 2.0*dat.lat.bz_volume*pow(e,2)/(pi*eV*a*angstrom);
    double SCALEFACTOR_LL = 2.0*dat.lat.bz_volume*pow(e,2)/(pi*eV*a*angstrom);           // Longitudinal-longitudinal
    double SCALEFACTOR_TT = 32.0*dat.lat.bz_volume*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5)); // Transverse-transverse  
    double SCALEFACTOR_LT = 8.0*dat.lat.bz_volume*pi*pow(hbar,2)*pow(e,2)/(pow(eV,2)*e_m*pow(a*angstrom,3));  // Longitudinal-transverse

    // STEP 3: Allocate memory for susceptibility tensor matrix elements
    // Structure: ReXij[q][i][j][m][n][f] and ImXij[q][i][j][m][n][f]
    // where q=q-vector, i,j=cartesian indices, m,n=G-vector indices, f=frequency
    dat.ReXij = new double *****[NQ];      // Real parts
    dat.ImXij = new double *****[NQ];      // Imaginary parts
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
                    for (int n=0; n<=m; n++){  // Only store upper triangular part (m >= n)
                        dat.ReXij[q][i][j][m][n] = new double [NFREQ];
                        dat.ImXij[q][i][j][m][n] = new double [NFREQ];
                    }
                }
            }
        }
    }

    //======================================================================
    // PHASE 1: CALCULATE CONDUCTION BAND STATES |k,c⟩
    //======================================================================
    std::cout << "  Solving for |k,v>..." << std::endl;
    time_0 = std::chrono::system_clock::now();    int NBAND_C [NKPT]; // Number of conduction bands per k-point
    double **E_k; // Conduction band eigenvalues: E_k[k-point][band]
    E_k = new double *[NKPT];
    std::complex<double> ***C_k; // Conduction band eigenvectors: C_k[k-point][band][G-vector]
    C_k = new std::complex<double> **[NKPT];
    
    // Parallelize over k-points for efficiency
    #pragma omp parallel for
    for (int k=0; k<NKPT; k++){
        // Get k-vector from the lattice structure
        Eigen::Vector3d K = dat.lat.K[k];

        // Construct Hamiltonian matrix H(k) using empirical pseudopotentials
        // This returns a matrix of size NPW×NPW (or 2*NPW×2*NPW if SOC is included)
        Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K,dat.lat.atomic,dat.mat);   
        
        // Solve eigenvalue problem: H(k)|ψ_n(k)⟩ = E_n(k)|ψ_n(k)⟩
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);

        // Allocate memory for eigenvalues at this k-point
        E_k[k] = new double [NBAND];
        
        // Extract eigenvalues and determine conduction band count
        // Note: Eigenvalues are ordered from lowest to highest
        // We start from the highest energy bands (NBAND-1-c) and work down
        for (int c=0; c<NBAND; c++){
            E_k[k][c] = eigsolver_k.eigenvalues()(NBAND-1-c)-dat.energyoffset;
            
            // Count conduction bands (those with positive energy)
            if (E_k[k][c]<=0) {
                NBAND_C[k] = c;  // c becomes the number of conduction bands
                break;
            }
        }
        
        // Store only conduction band eigenvectors to save memory
        C_k[k] = new std::complex<double> *[NBAND_C[k]];
        for (int c=0; c<NBAND_C[k]; c++){
            C_k[k][c] = new std::complex<double> [NPW];
            for (int p=0; p<NPW; p++){
                // Extract eigenvector components (plane wave coefficients)
                C_k[k][c][p] = eigsolver_k.eigenvectors().col(NBAND-1-c)(p);
            }// loop over plane wave components
        }// loop over conduction bands
    }    // Synchronization barrier to ensure all k-point calculations are complete
    #pragma omp barrier
    time_1 = std::chrono::system_clock::now();
    elapsed_time = time_1-time_0;
    std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;

    //======================================================================
    // PHASE 2: PREPARE OVERLAP INTEGRAL STORAGE 
    //======================================================================
    
    // Allocate memory for overlap integrals between valence and conduction states
    // These represent matrix elements like ⟨k,c|e^{-i(q+G_m)·r}|k+q,v⟩
    // Structure: oint[k-point][cartesian_component][G-vector][conduction_band][valence_band]
    std::complex<double> *****oint;
    oint = new std::complex<double> ****[NKPT]; 
    for (int k=0; k<NKPT; k++){
        oint[k] = new std::complex<double> ***[3];      // 3 Cartesian components (x,y,z)
        for (int i=0; i<NTSR; i++){
            oint[k][i] = new std::complex<double> **[NEPS];  // NEPS G-vectors
            for (int m=0; m<NEPS; m++){
                oint[k][i][m] = new std::complex<double> *[NBAND_C[k]];  // Conduction bands
            }
        }
    }

    //======================================================================
    // PHASE 3: MAIN SUSCEPTIBILITY CALCULATION LOOP
    //======================================================================
    std::cout << "  Solving for susceptibility tensor matrix elements (alternative form)..." << std::endl;
    
    // Variables for valence band states at k+q
    int NBAND_V [NKPT]; // Number of valence bands per k-point  
    std::complex<double> ***C_kq; // Valence band eigenvectors at k+q
    double **E_kq; // Valence band eigenvalues at k+q
    E_kq = new double *[NKPT];
    C_kq = new std::complex<double> **[NKPT];
    for (int k=0; k<NKPT; k++){
        E_kq[k] = new double [NBAND];  // Allocate for all bands initially
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

/*
    ============================================================================
    ADDITIONAL NOTES ON LOGICAL FLOW FOR REMAINING FUNCTIONS:
    ============================================================================
    
    The file contains several main calculation functions, each following similar patterns:
    
    COMMON PATTERN FOR ALL CHI CALCULATIONS:
    ---------------------------------------
    
    1. SETUP PHASE:
       - Select broadening function (Gaussian, Lorentzian, etc.)
       - Define physical scaling factors for different field orientations
       - Allocate memory for susceptibility tensor χ[q][i][j][G][G'][ω]
    
    2. ELECTRONIC STRUCTURE CALCULATION:
       - For each k-point: solve H(k)|ψ_n(k)⟩ = E_n(k)|ψ_n(k)⟩
       - Separate valence (E < 0) and conduction (E > 0) bands
       - Store eigenvectors C_k[k][n][G] and eigenvalues E_k[k][n]
    
    3. OVERLAP INTEGRAL COMPUTATION:
       - For each q-vector: calculate band structure at k+q
       - Compute matrix elements: ⟨k,c|e^{i(q+G)·r}|k+q,v⟩
       - These represent probability amplitudes for electronic transitions
    
    4. SUSCEPTIBILITY TENSOR CALCULATION:
       - Apply Fermi's Golden Rule: transitions conserve energy
       - Use delta function: δ(E_c(k) - E_v(k+q) - ħω)
       - Apply broadening to discrete energy levels
       - Sum over all k-points and band combinations
       - Include appropriate geometric factors and scaling
    
    5. TENSOR COMPONENT SPECIFIC FACTORS:
       - LL (Longitudinal-Longitudinal): charge density response
       - TT (Transverse-Transverse): current density response  
       - LT (Longitudinal-Transverse): mixed field coupling
    
    6. POST-PROCESSING:
       - Apply Kramers-Kronig relations if requested
       - Add δ(G,G') terms for macroscopic response
       - Store results in output arrays
    
    MEMORY MANAGEMENT CONSIDERATIONS:
    --------------------------------
    - Large arrays like eigenvectors require careful allocation/deallocation
    - OpenMP parallelization must account for memory usage per thread
    - Overlap integrals can be memory-intensive for large systems
    
    SPIN-ORBIT COUPLING (SOC) MODIFICATIONS:
    ---------------------------------------
    - When SOC is enabled, Hamiltonian size doubles: NPW → 2*NPW
    - Eigenvectors become spinors with both spin-up and spin-down components
    - Overlap integrals must sum over both spin components
    - Band indexing accounts for doubled Hilbert space
    
    KEY PHYSICS INSIGHTS:
    --------------------
    - χ(q,ω) peaks occur at inter-band transition energies
    - q→0 limit gives macroscopic dielectric response
    - Different tensor components probe different field geometries
    - Broadening mimics finite lifetime and instrumental resolution
    
    ============================================================================
*/

/*
    END OF FILE: pmx_chi.cpp
    
    Summary: This file implements the core susceptibility tensor calculations
    for the PicoMax empirical pseudopotential method code. The functions
    compute χ(q,ω) using different approaches (charge-based, full tensor,
    simplified xyz components) but all follow the same fundamental physics
    of electronic transitions between valence and conduction bands.
*/
